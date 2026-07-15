from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
from pathlib import Path

import numpy as np

from .cache import CandidateCache
from .candidates import LocalCandidateProvider, WAVELENGTH_FOR_RADIATION, simulate_or_load_peaks
from .contribution import ContributionProxy, contribution_proxy
from .importers import infer_xrd_range, load_xrd_file
from .matching import explain_multiphase, prioritize_impurities_after_main, score_candidate
from .models import Candidate, CandidateScore, MultiPhaseExplanation, Peak
from .peaks import find_peaks
from .phase_context import phase_context_cards
from .plotting import write_xrd_plot
from .report import write_reports
from .audit import write_analysis_audit_folder
from .simulation_trust import apply_vesta_trust_check
from .config import XRDConfig
from .rietan_backend import RietanXRDBackend, upgrade_scores_with_rietan
from ophiuchus.plotting.xrd_plotter import render_xrd_from_analysis_json
from ophiuchus.phase_stripping.models import AnalysisContext


@dataclass
class AnalysisResult:
    experimental_peaks: list[Peak]
    candidates: list[Candidate]
    top_scores: list[CandidateScore]
    explanation: MultiPhaseExplanation
    outputs: dict[str, str]
    warnings: list[str]
    contribution: ContributionProxy | None = None
    phase_context: list[dict[str, object]] | None = None
    phase_candidates: list[object] = field(default_factory=list)
    target_score: CandidateScore | None = None
    impurity_scores: list[CandidateScore] = field(default_factory=list)
    scientific_runtime: dict[str, object] = field(default_factory=dict)
    context: AnalysisContext | None = None


def build_analysis_context(
    xrd_file: str | Path,
    *,
    radiation: str,
    two_theta_range: tuple[float, float],
    tolerance_deg: float,
    wavelength_angstrom: float | None = None,
) -> AnalysisContext:
    raw_pattern = load_xrd_file(xrd_file, normalize=False)
    lower, upper = map(float, two_theta_range)
    selected = [point for point in raw_pattern.points if lower <= point.two_theta <= upper]
    if not selected:
        raise ValueError("The experimental XRD pattern has no points inside the active 2theta range.")
    x = np.asarray([point.two_theta for point in selected], dtype=np.float64)
    intensity = np.asarray([point.intensity for point in selected], dtype=np.float64)
    source = Path(xrd_file)
    source_fingerprint = hashlib.sha256(source.read_bytes()).hexdigest()
    data_hash = hashlib.sha256()
    data_hash.update(x.dtype.str.encode("ascii"))
    data_hash.update(np.ascontiguousarray(x).tobytes())
    data_hash.update(intensity.dtype.str.encode("ascii"))
    data_hash.update(np.ascontiguousarray(intensity).tobytes())
    return AnalysisContext(
        x=x,
        intensity=intensity,
        radiation=radiation,
        wavelength_angstrom=float(
            wavelength_angstrom
            if wavelength_angstrom is not None
            else WAVELENGTH_FOR_RADIATION.get(radiation, 1.54056)
        ),
        two_theta_range=(lower, upper),
        tolerance_deg=float(tolerance_deg),
        source_path=str(source),
        source_fingerprint=source_fingerprint,
        data_fingerprint=data_hash.hexdigest(),
    )


def run_analysis(
    xrd_file: str | Path,
    candidate_dirs: list[str | Path],
    elements: list[str],
    extra_elements: list[str] | None = None,
    out_dir: str | Path = "results/ophi_xrd",
    radiation: str = "CuKa",
    two_theta_min: float | None = None,
    two_theta_max: float | None = None,
    tolerance_deg: float = 0.20,
    max_candidates: int = 20,
    max_phases: int = 3,
    cache_path: str | Path | None = None,
    use_cache: bool = True,
    synthesis_metadata: dict[str, str] | None = None,
    target_candidate_id: str | None = None,
    target_formula: str | None = None,
    use_rietan_display: bool = False,
) -> AnalysisResult:
    allowed = set(elements) | set(extra_elements or [])
    range_info = infer_xrd_range(xrd_file)
    auto_two_theta_min = float(range_info["two_theta_min"])
    auto_two_theta_max = float(range_info["two_theta_max"])
    if two_theta_min is None or two_theta_max is None:
        two_theta_min, two_theta_max = auto_two_theta_min, auto_two_theta_max
    analysis_context = build_analysis_context(
        xrd_file,
        radiation=radiation,
        two_theta_range=(two_theta_min, two_theta_max),
        tolerance_deg=tolerance_deg,
    )
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
    provider = LocalCandidateProvider(candidate_dirs, allowed_elements=allowed or None)
    candidates = provider.iter_candidates()
    cache = CandidateCache(cache_path) if cache_path else None
    warnings: list[str] = []
    scores: list[CandidateScore] = []
    for candidate in candidates:
        try:
            is_cif = Path(candidate.source_path).suffix.lower() == ".cif"
            # CIF analysis always goes through the canonical backend. The legacy
            # peak-only cache remains useful for imported peak lists, but cannot
            # silently replace a structure-derived pattern.
            cached_peaks = (
                cache.load_pattern(candidate, radiation, (two_theta_min, two_theta_max))
                if cache and use_cache and not is_cif
                else None
            )
            if cached_peaks:
                candidate.theory_peaks = cached_peaks
            else:
                candidate.theory_peaks = simulate_or_load_peaks(
                    candidate,
                    radiation=radiation,
                    two_theta_range=(two_theta_min, two_theta_max),
                )
                if cache and candidate.theory_peaks:
                    cache.store_pattern(candidate, candidate.theory_peaks, radiation, (two_theta_min, two_theta_max))
        except RuntimeError as exc:
            warnings.append(f"{candidate.formula_pretty}: {exc}")
            continue
        if not candidate.theory_peaks:
            warnings.append(f"{candidate.formula_pretty}: no theory peaks available")
            continue
        apply_vesta_trust_check(candidate, two_theta_min, two_theta_max)
        _attach_pattern_identity(candidate)
        scores.append(score_candidate(experimental_peaks, candidate, tolerance_deg=tolerance_deg))
    top_scores = prioritize_impurities_after_main(
        scores,
        experimental_peaks,
        main_elements=elements,
        target_candidate_id=target_candidate_id,
        target_formula=target_formula,
    )[:max_candidates]
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
    explanation = explain_multiphase(experimental_peaks, top_scores, max_phases=max_phases)
    contribution = contribution_proxy(explanation)
    context = phase_context_cards([score.candidate for score in top_scores[:8]], synthesis_metadata=synthesis_metadata)
    input_summary = {
        "xrd_file": str(xrd_file),
        "number_of_experimental_points": len(pattern.points),
        "number_of_extracted_peaks": len(experimental_peaks),
        "allowed_elements": sorted(allowed),
        "radiation": radiation,
        "tolerance_deg": tolerance_deg,
        "two_theta_range": [two_theta_min, two_theta_max],
        "auto_range": range_info,
        "cache_path": str(cache_path) if cache_path else "",
        "target_candidate_id": target_candidate_id or "",
        "target_formula": target_formula or "",
        "canonical_cif_backend": "ValidatedXRDBackend",
        "legacy_cif_cache_reads_disabled": True,
        "display_backend": "VESTA / RIETAN-FP" if use_rietan_display else "canonical screening backend",
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
        out_dir,
        input_summary,
        experimental_peaks,
        top_scores,
        explanation,
        contribution=contribution,
        phase_context=context,
        xrd_plot=plot_path,
    )
    if presentation_plot:
        outputs["xrd_presentation_plot"] = presentation_plot
    try:
        clean_plot = render_xrd_from_analysis_json(outputs["json"], Path(out_dir) / "xrd_plot_clean.png")
        outputs["xrd_plot_clean"] = clean_plot["png"]
        diagnostic_plot = render_xrd_from_analysis_json(outputs["json"], Path(out_dir) / "xrd_diagnostic.png")
        outputs["xrd_diagnostic_plot"] = diagnostic_plot["png"]
    except Exception as exc:
        warnings.append(f"clean plot export failed: {exc}")
    try:
        outputs["audit_folder"] = write_analysis_audit_folder(
            out_dir,
            input_summary,
            experimental_peaks,
            top_scores,
            candidates,
            outputs,
        )
    except Exception as exc:
        warnings.append(f"audit folder export failed: {exc}")
    return AnalysisResult(
        experimental_peaks=experimental_peaks,
        candidates=candidates,
        top_scores=top_scores,
        explanation=explanation,
        outputs=outputs,
        warnings=warnings,
        contribution=contribution,
        phase_context=context,
        scientific_runtime={
            "scientific_safe_mode": True,
            "backend_name": "ValidatedXRDBackend",
            "display_backend": "VESTA / RIETAN-FP" if use_rietan_display else "ValidatedXRDBackend",
            "legacy_cif_cache_reads_disabled": True,
        },
        context=analysis_context,
    )


def _attach_pattern_identity(candidate: Candidate) -> None:
    pattern = candidate.simulated_pattern
    if pattern is not None:
        candidate.simulation_validation.setdefault("pattern_fingerprint", pattern.pattern_fingerprint)
        candidate.simulation_validation.setdefault("cif_sha256", pattern.cif_sha256)
        return
    digest = hashlib.sha256()
    for peak in sorted(candidate.theory_peaks, key=lambda item: (item.two_theta, item.intensity, item.peak_id)):
        digest.update(f"{peak.two_theta:.8f}|{peak.intensity:.8f}|{peak.hkl}\n".encode("utf-8"))
    if candidate.theory_peaks:
        candidate.simulation_validation.setdefault("pattern_fingerprint", digest.hexdigest())
