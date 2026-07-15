from __future__ import annotations

from pathlib import Path

from ophiuchus.xrd.contribution import contribution_proxy
from ophiuchus.xrd.importers import infer_xrd_range, load_xrd_file
from ophiuchus.xrd.matching import explain_multiphase, prioritize_impurities_after_main, score_candidate
from ophiuchus.xrd.peaks import find_peaks
from ophiuchus.xrd.phase_context import phase_context_cards
from ophiuchus.xrd.audit import write_analysis_audit_folder
from ophiuchus.xrd.pipeline import AnalysisResult, build_analysis_context
from ophiuchus.xrd.plotting import write_xrd_plot
from ophiuchus.xrd.report import write_reports
from ophiuchus.xrd.simulation_trust import apply_vesta_trust_check
from ophiuchus.plotting.xrd_plotter import render_xrd_from_analysis_json

from .candidate_service import CandidateStructureService, write_candidate_usage_summary
from .database import StructureLibrary
from .phase_grouping import PhaseCandidate, group_phase_candidates
from .xrd_cache import build_library_xrd_cache, load_validated_pattern, scoped_structure_ids
from ophiuchus.xrd.models import Candidate
from ophiuchus.xrd.config import XRDConfig
from ophiuchus.xrd.rietan_backend import RietanXRDBackend, upgrade_scores_with_rietan


def run_library_analysis(
    xrd_file: str | Path,
    library_db: str | Path,
    elements: list[str],
    extra_elements: list[str] | None = None,
    out_dir: str | Path = "results/ophi_library_xrd",
    project: str = "OphiProject",
    sample_id: str = "Sample",
    radiation: str = "CuKa",
    two_theta_min: float | None = None,
    two_theta_max: float | None = None,
    tolerance_deg: float = 0.20,
    max_candidates: int = 20,
    max_phases: int = 3,
    candidate_scope: str = "subsystems",
    target_candidate_id: str | None = None,
    target_formula: str | None = None,
    scientific_safe_mode: bool = True,
    force_recompute: bool = False,
    use_rietan_display: bool = False,
) -> AnalysisResult:
    allowed = sorted(set(elements) | set(extra_elements or []))
    range_info = infer_xrd_range(xrd_file)
    auto_two_theta_min = float(range_info["two_theta_min"])
    auto_two_theta_max = float(range_info["two_theta_max"])
    if two_theta_min is None or two_theta_max is None:
        two_theta_min, two_theta_max = auto_two_theta_min, auto_two_theta_max
    library = StructureLibrary(library_db)
    compatible_ids = scoped_structure_ids(library, allowed, enabled_only=True, scope=candidate_scope)
    if target_candidate_id and target_candidate_id not in compatible_ids:
        raise ValueError(
            "Selected target phase is not enabled or is outside the confirmed element scope: "
            f"{target_candidate_id}"
        )
    build_summary = build_library_xrd_cache(
        library,
        structure_ids=compatible_ids,
        radiation=radiation,
        two_theta_range=(two_theta_min, two_theta_max),
        force=force_recompute,
        scientific_safe_mode=scientific_safe_mode,
    )
    analysis_context = build_analysis_context(
        xrd_file,
        radiation=radiation,
        two_theta_range=(two_theta_min, two_theta_max),
        tolerance_deg=tolerance_deg,
        wavelength_angstrom=float(build_summary.get("wavelength_angstrom") or 1.54056),
    )
    candidate_service = CandidateStructureService(library, radiation=radiation, two_theta_range=(two_theta_min, two_theta_max), candidate_scope=candidate_scope)
    pattern = load_xrd_file(xrd_file, normalize=True)
    experimental_peaks = find_peaks(
        pattern.points,
        two_theta_min=two_theta_min,
        two_theta_max=two_theta_max,
        min_height=1.0,
        min_distance_deg=0.15,
        max_peaks=80,
        smooth_window=3,
        min_prominence=3.0,
    )
    compatible_entries = [library.get_structure(internal_id) for internal_id in compatible_ids]
    patterns_by_id = {}
    for entry in compatible_entries:
        simulated = load_validated_pattern(
            library,
            entry,
            radiation=radiation,
            two_theta_range=(two_theta_min, two_theta_max),
        )
        if simulated is not None:
            patterns_by_id[entry.internal_id] = simulated
    phase_candidates = group_phase_candidates(
        compatible_entries,
        patterns_by_id,
        target_structure_id=target_candidate_id,
        simulation_state_by_id=build_summary.get("structure_states", {}),
    )
    candidates = [_phase_to_candidate(phase) for phase in phase_candidates]
    for candidate in candidates:
        apply_vesta_trust_check(candidate, two_theta_min, two_theta_max)
        if candidate.simulated_pattern is not None:
            candidate.simulation_validation.setdefault(
                "pattern_fingerprint", candidate.simulated_pattern.pattern_fingerprint
            )
            candidate.simulation_validation.setdefault("cif_sha256", candidate.simulated_pattern.cif_sha256)
    used_structure_ids = {entry.internal_id for phase in phase_candidates for entry in phase.entries}
    source_summary = _candidate_source_summary(library, allowed, used_structure_ids, radiation, (two_theta_min, two_theta_max), candidate_scope=candidate_scope)
    scores = [score_candidate(experimental_peaks, candidate, tolerance_deg=tolerance_deg) for candidate in candidates if candidate.theory_peaks]
    candidate_usage_rows = candidate_service.list_candidate_rows(allowed, used_ids=used_structure_ids)
    candidate_usage_summary = candidate_service.usage_summary(candidate_usage_rows, used_ids=used_structure_ids)
    target_phase = next(
        (phase for phase in phase_candidates if target_candidate_id and phase.contains_structure(target_candidate_id)),
        None,
    )
    top_scores = prioritize_impurities_after_main(
        scores,
        experimental_peaks,
        main_elements=elements,
        target_candidate_id=None if target_phase is None else target_phase.phase_id,
        target_formula=target_formula,
    )[:max_candidates]
    warnings = []
    if use_rietan_display and top_scores:
        display_backend = RietanXRDBackend()
        if display_backend.available:
            top_scores, display_warnings = upgrade_scores_with_rietan(
                top_scores,
                display_backend,
                XRDConfig(
                    radiation_source=radiation,
                    two_theta_min=float(two_theta_min),
                    two_theta_max=float(two_theta_max),
                    line_model="cu_kalpha12",
                ),
                experimental_peaks,
                tolerance_deg=tolerance_deg,
                limit=4,
            )
            warnings.extend(display_warnings)
        else:
            warnings.append("VESTA/RIETAN display simulation is enabled but the local executables are unavailable.")
    target_score = next(
        (score for score in top_scores if target_phase is not None and score.candidate.candidate_id == target_phase.phase_id),
        top_scores[0] if top_scores else None,
    )
    impurity_scores = [score for score in top_scores if score is not target_score][:3]
    explanation = explain_multiphase(experimental_peaks, top_scores, max_phases=max_phases)
    contribution = contribution_proxy(explanation)
    context = phase_context_cards([score.candidate for score in top_scores[:8]])
    library_entry_ids = sorted(patterns_by_id)
    input_summary = {
        "project": project,
        "sample_id": sample_id,
        "xrd_file": str(xrd_file),
        "library_db": str(library_db),
        "library_entry_ids": library_entry_ids,
        "compatible_structure_ids": compatible_ids,
        "used_candidate_structure_ids": sorted(used_structure_ids),
        "candidate_scope": candidate_scope,
        "candidate_scope_note": "Default subsystems scope analyzes structures whose elements are all within the selected main element range; impurity elements are ignored unless explicitly enabled in a later scan.",
        "candidate_usage": candidate_usage_summary,
        "candidate_usage_rows": [row.to_dict() for row in candidate_usage_rows],
        "candidate_source_summary": source_summary,
        "library_cache_summary": build_summary,
        "number_of_experimental_points": len(pattern.points),
        "number_of_extracted_peaks": len(experimental_peaks),
        "allowed_elements": allowed,
        "radiation": radiation,
        "tolerance_deg": tolerance_deg,
        "two_theta_range": [two_theta_min, two_theta_max],
        "auto_range": range_info,
        "target_candidate_id": target_candidate_id or "",
        "target_formula": target_formula or "",
        "phase_candidate_count": len(phase_candidates),
        "scientific_safe_mode": scientific_safe_mode,
        "force_recompute": force_recompute,
        "validated_backend": build_summary.get("backend_name", ""),
        "validated_backend_version": build_summary.get("backend_version", ""),
        "display_backend": "VESTA / RIETAN-FP" if use_rietan_display else "validated screening backend",
    }
    out_path = Path(out_dir)
    plot_path = None
    presentation_plot = None
    try:
        plot_path = write_xrd_plot(out_path / "xrd_plot.png", pattern, experimental_peaks, top_scores)
        presentation_plot = write_xrd_plot(out_path / "xrd_presentation.png", pattern, experimental_peaks, top_scores)
    except Exception as exc:
        warnings.append(f"plot export failed: {exc}")
    outputs = write_reports(
        out_path,
        input_summary,
        experimental_peaks,
        top_scores,
        explanation,
        contribution=contribution,
        phase_context=context,
        xrd_plot=plot_path,
    )
    outputs["candidate_usage_summary"] = write_candidate_usage_summary(out_path / "candidate_usage_summary.csv", candidate_usage_rows)
    if presentation_plot:
        outputs["xrd_presentation_plot"] = presentation_plot
    try:
        clean_plot = render_xrd_from_analysis_json(outputs["json"], out_path / "xrd_plot_clean.png")
        outputs["xrd_plot_clean"] = clean_plot["png"]
        diagnostic_plot = render_xrd_from_analysis_json(outputs["json"], out_path / "xrd_diagnostic.png")
        outputs["xrd_diagnostic_plot"] = diagnostic_plot["png"]
    except Exception as exc:
        warnings.append(f"clean plot export failed: {exc}")
    try:
        outputs["audit_folder"] = write_analysis_audit_folder(
            out_path,
            input_summary,
            experimental_peaks,
            top_scores,
            candidates,
            outputs,
            library_entries=library.list_structures(enabled_only=False),
        )
    except Exception as exc:
        warnings.append(f"audit folder export failed: {exc}")
    for warning in build_summary.get("warnings", []):
        warnings.append(str(warning))
    return AnalysisResult(
        experimental_peaks=experimental_peaks,
        candidates=candidates,
        top_scores=top_scores,
        explanation=explanation,
        outputs=outputs,
        warnings=warnings,
        contribution=contribution,
        phase_context=context,
        phase_candidates=phase_candidates,
        target_score=target_score,
        impurity_scores=impurity_scores,
        scientific_runtime={
            "scientific_safe_mode": scientific_safe_mode,
            "force_recompute": force_recompute,
            "backend_name": build_summary.get("backend_name", ""),
            "backend_version": build_summary.get("backend_version", ""),
            "display_backend": "VESTA / RIETAN-FP" if use_rietan_display else build_summary.get("backend_name", ""),
            "freshly_simulated": build_summary.get("freshly_simulated", 0),
            "validated_cache_hits": build_summary.get("validated_cache_hits", 0),
            "radiation": build_summary.get("radiation", radiation),
            "wavelength_angstrom": build_summary.get("wavelength_angstrom", ""),
            "structure_records_checked": build_summary.get("checked", 0),
            "simulation_failures": build_summary.get("failed", 0),
            "grouped_phase_candidates": len(phase_candidates),
        },
        context=analysis_context,
    )


def _phase_to_candidate(phase: PhaseCandidate) -> Candidate:
    entry = phase.representative
    candidate = Candidate(
        candidate_id=phase.phase_id,
        formula_pretty=entry.reduced_formula or entry.formula,
        source=f"library:{entry.source}",
        source_path=phase.pattern.cif_path,
        elements=entry.elements,
        structure_hash=entry.structure_hash,
        parse_status="validated_backend_pattern",
        simulated_pattern=phase.pattern,
        phase_entry_ids=[item.internal_id for item in phase.entries],
        space_group_symbol=entry.space_group_symbol or "",
        space_group_number=entry.space_group_number,
        simulation_state=phase.simulation_state,
    )
    candidate.theory_peaks = phase.pattern.to_peaks()
    candidate.simulation_validation["pattern_fingerprint"] = phase.pattern.pattern_fingerprint
    candidate.simulation_validation["cif_sha256"] = phase.pattern.cif_sha256
    return candidate


def _candidate_source_summary(
    library: StructureLibrary,
    allowed_elements: list[str],
    used_candidate_ids: set[str],
    radiation: str,
    two_theta_range: tuple[float, float],
    candidate_scope: str = "exact",
) -> dict[str, dict[str, int]]:
    from .xrd_cache import simulation_settings_hash

    allowed = set(allowed_elements)
    settings_hash = simulation_settings_hash(radiation, two_theta_range)
    summary: dict[str, dict[str, int]] = {}
    for entry in library.list_structures(enabled_only=False):
        bucket = summary.setdefault(
            entry.source,
            {
                "loaded": 0,
                "enabled": 0,
                "simulated": 0,
                "used": 0,
                "skipped_disabled": 0,
                "skipped_element_filter": 0,
                "skipped_scope_filter": 0,
                "skipped_no_valid_xrd_cache": 0,
            },
        )
        bucket["loaded"] += 1
        if entry.enabled_for_matching:
            bucket["enabled"] += 1
        else:
            bucket["skipped_disabled"] += 1
        has_cache = bool(library.load_xrd_peaks(entry.internal_id, radiation, settings_hash))
        if has_cache:
            bucket["simulated"] += 1
        entry_elements = set(entry.elements)
        in_element_range = (not allowed) or entry_elements.issubset(allowed)
        in_scope = (not allowed) or (entry_elements.issubset(allowed) if candidate_scope == "subsystems" else entry_elements == allowed)
        if not in_element_range:
            bucket["skipped_element_filter"] += 1
        elif not in_scope:
            bucket["skipped_scope_filter"] += 1
        if entry.enabled_for_matching and not has_cache:
            bucket["skipped_no_valid_xrd_cache"] += 1
        if entry.internal_id in used_candidate_ids:
            bucket["used"] += 1
    return summary
