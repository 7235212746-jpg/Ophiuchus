from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from .contribution import ContributionProxy
from .models import CandidateScore, MultiPhaseExplanation, Peak


CAVEAT = (
    "This result is candidate screening only. It does not prove phase identity, "
    "does not perform Rietveld refinement, and does not estimate phase fractions. "
    "Use it to prioritize manual checks and follow-up experiments."
)


def write_reports(
    out_dir: str | Path,
    input_summary: dict[str, Any],
    experimental_peaks: list[Peak],
    top_scores: list[CandidateScore],
    explanation: MultiPhaseExplanation,
    contribution: ContributionProxy | None = None,
    phase_context: list[dict[str, object]] | None = None,
    xrd_plot: str | None = None,
) -> dict[str, str]:
    output = Path(out_dir)
    output.mkdir(parents=True, exist_ok=True)
    payload = _payload(input_summary, experimental_peaks, top_scores, explanation, contribution, phase_context, xrd_plot)
    json_path = output / "results.json"
    md_path = output / "report.md"
    peaks_path = output / "extracted_peaks.csv"
    candidates_path = output / "top_candidates.csv"
    assignments_path = output / "peak_assignments.csv"
    contribution_path = output / "contribution_proxy.csv"
    json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    md_path.write_text(_markdown(payload), encoding="utf-8")
    _write_peaks_csv(peaks_path, experimental_peaks)
    _write_candidates_csv(candidates_path, top_scores)
    _write_assignments_csv(assignments_path, experimental_peaks, top_scores)
    _write_contribution_csv(contribution_path, contribution)
    outputs = {
        "json": str(json_path),
        "markdown": str(md_path),
        "extracted_peaks": str(peaks_path),
        "top_candidates": str(candidates_path),
        "peak_assignments": str(assignments_path),
        "contribution_proxy": str(contribution_path),
    }
    if xrd_plot:
        outputs["xrd_plot"] = xrd_plot
    return outputs


def _payload(
    input_summary: dict[str, Any],
    experimental_peaks: list[Peak],
    top_scores: list[CandidateScore],
    explanation: MultiPhaseExplanation,
    contribution: ContributionProxy | None = None,
    phase_context: list[dict[str, object]] | None = None,
    xrd_plot: str | None = None,
) -> dict[str, Any]:
    return {
        "input": input_summary,
        "experimental_peaks": [_peak_dict(p) for p in experimental_peaks],
        "top_candidates": [_score_dict(score, rank) for rank, score in enumerate(top_scores, 1)],
        "possible_impurity_phases": [_score_dict(score, rank) for rank, score in enumerate(top_scores[1:4], 2)],
        "per_peak_candidate_assignments": _per_peak_assignments(experimental_peaks, top_scores),
        "multi_phase_explanation": {
            "selected_candidates": [score.candidate.formula_pretty for score in explanation.selected_candidates],
            "explained_experimental_peaks": [_peak_dict(p) for p in explanation.explained_experimental_peaks],
            "unexplained_experimental_peaks": [_peak_dict(p) for p in explanation.unexplained_experimental_peaks],
            "warnings": explanation.warnings,
        },
        "relative_contribution_proxy": {
            "values": {} if contribution is None else contribution.contributions,
            "warning": "" if contribution is None else contribution.warning,
        },
        "phase_context": phase_context or [],
        "xrd_plot": xrd_plot or "",
        "scientific_caveat": CAVEAT,
    }


def _score_dict(score: CandidateScore, rank: int) -> dict[str, Any]:
    role = "likely main phase" if rank == 1 else "possible impurity" if score.score >= 0.35 else "weak candidate"
    confidence = str(score.score_components.get("confidence_label") or ("medium" if score.score >= 0.55 else "low" if score.score >= 0.30 else "uncertain"))
    return {
        "rank": rank,
        "candidate_id": score.candidate.candidate_id,
        "formula": score.candidate.formula_pretty,
        "score": round(score.score, 4),
        "score_components": score.score_components,
        "simulation_validation": score.candidate.simulation_validation,
        "confidence_label": confidence,
        "source": score.candidate.source,
        "source_path": score.candidate.source_path,
        "space_group_symbol": score.candidate.space_group_symbol,
        "space_group_number": score.candidate.space_group_number,
        "phase_entry_ids": list(score.candidate.phase_entry_ids),
        "simulated_pattern": None
        if score.candidate.simulated_pattern is None
        else score.candidate.simulated_pattern.to_dict(),
        "role": role,
        "matched_peaks": [
            {
                "theory_two_theta": round(m.theory_peak.two_theta, 4),
                "exp_two_theta": round(m.experimental_peak.two_theta, 4),
                "delta": round(m.delta, 4),
                "theory_intensity": round(m.theory_peak.intensity, 2),
            }
            for m in score.matched_theory_peaks
        ],
        "missing_strong_peaks": [_peak_dict(p) for p in score.missing_strong_theory_peaks],
        "warnings": score.warnings,
    }


def _per_peak_assignments(experimental_peaks: list[Peak], top_scores: list[CandidateScore]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for peak in experimental_peaks:
        matches = []
        for score in top_scores:
            for match in score.matched_theory_peaks:
                if match.experimental_peak.peak_id == peak.peak_id:
                    matches.append(
                        {
                            "candidate_id": score.candidate.candidate_id,
                            "formula": score.candidate.formula_pretty,
                            "source": score.candidate.source,
                            "calculated_two_theta": round(match.theory_peak.two_theta, 4),
                            "delta_two_theta": round(match.delta, 4),
                            "hkl": match.theory_peak.hkl,
                            "calculated_relative_intensity": round(match.theory_peak.intensity, 4),
                            "confidence_level": "likely" if score.score >= 0.7 else "possible" if score.score >= 0.35 else "weak",
                        }
                    )
        rows.append(
            {
                "experimental_peak_id": peak.peak_id,
                "experimental_two_theta": round(peak.two_theta, 4),
                "experimental_intensity": round(peak.intensity, 4),
                "assignments": matches,
                "status": "assigned" if matches else "unresolved",
            }
        )
    return rows


def _peak_dict(peak: Peak) -> dict[str, Any]:
    return {
        "peak_id": peak.peak_id,
        "two_theta": round(peak.two_theta, 4),
        "intensity": round(peak.intensity, 4),
        "prominence": None if peak.prominence is None else round(peak.prominence, 4),
        "fwhm": peak.fwhm,
    }


def _markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# Ophi XRD Candidate Screener Report",
        "",
        "## Input Summary",
    ]
    for key, value in payload["input"].items():
        lines.append(f"- **{key}**: {value}")
    lines.extend(["", "## Top Candidate Phases"])
    for item in payload["top_candidates"]:
        lines.append(f"### {item['rank']}. {item['formula']} - score {item['score']}")
        lines.append(f"- Role: {item['role']}")
        lines.append(f"- Confidence: {item['confidence_label']} (screening level)")
        lines.append(f"- Source: `{item['source']}` `{item['source_path']}`")
        lines.append(f"- Structure ID: `{item['candidate_id']}`")
        lines.append(f"- Matched peaks: {len(item['matched_peaks'])}")
        if item["score_components"]:
            lines.append(
                "- Score components: "
                f"matched fraction {item['score_components'].get('matched_strong_intensity_fraction')}, "
                f"missing penalty {item['score_components'].get('missing_strong_penalty')}, "
                f"multi-peak bonus {item['score_components'].get('multi_peak_bonus')}"
            )
        validation = item.get("simulation_validation") or {}
        if validation:
            lines.append(
                "- Simulation validation: "
                f"{validation.get('status', 'unknown')}"
                + (f" against `{validation.get('reference_path')}`" if validation.get("reference_path") else "")
            )
        if item["missing_strong_peaks"]:
            missing = ", ".join(str(p["two_theta"]) for p in item["missing_strong_peaks"])
            lines.append(f"- Missing strong peaks: {missing}")
        for warning in item["warnings"]:
            lines.append(f"- Warning: {warning}")
    lines.extend(["", "## Possible Impurity Phases"])
    if payload["possible_impurity_phases"]:
        for item in payload["possible_impurity_phases"]:
            lines.append(f"- {item['formula']} ({item['source']}): score {item['score']}, matched {len(item['matched_peaks'])}, missing strong {len(item['missing_strong_peaks'])}")
    else:
        lines.append("- No possible impurity candidates above the current reporting threshold.")
    lines.extend(["", "## Per-Peak Candidate Assignments"])
    for item in payload["per_peak_candidate_assignments"][:30]:
        if item["assignments"]:
            labels = "; ".join(f"{a['formula']} Δ={a['delta_two_theta']}" for a in item["assignments"][:4])
            lines.append(f"- {item['experimental_two_theta']}: {labels}")
        else:
            lines.append(f"- {item['experimental_two_theta']}: unresolved")
    lines.extend(["", "## Multi-Phase Explanation"])
    selected = payload["multi_phase_explanation"]["selected_candidates"]
    lines.append(f"- Selected candidates: {', '.join(selected) if selected else 'none'}")
    unexplained = payload["multi_phase_explanation"]["unexplained_experimental_peaks"]
    lines.append(f"- Unexplained experimental peaks: {', '.join(str(p['two_theta']) for p in unexplained) or 'none'}")
    source_summary = payload["input"].get("candidate_source_summary")
    if source_summary:
        lines.extend(["", "## Candidate Source Summary"])
        for source, stats in source_summary.items():
            lines.append(
                f"- {source}: loaded {stats.get('loaded', 0)}, simulated {stats.get('simulated', 0)}, "
                f"enabled {stats.get('enabled', 0)}, used {stats.get('used', 0)}, "
                f"disabled {stats.get('skipped_disabled', 0)}, element-filtered {stats.get('skipped_element_filter', 0)}, "
                f"no valid XRD cache {stats.get('skipped_no_valid_xrd_cache', 0)}"
            )
    lines.extend(["", "## Relative Contribution Proxy"])
    proxy = payload["relative_contribution_proxy"]
    if proxy["values"]:
        for formula, value in proxy["values"].items():
            lines.append(f"- {formula}: {value:.1f}% screening contribution proxy")
    else:
        lines.append("- No contribution proxy available.")
    if proxy["warning"]:
        lines.append(f"- Warning: {proxy['warning']}")
    if payload.get("xrd_plot"):
        lines.extend(["", "## Plot", f"- XRD plot PNG: `{payload['xrd_plot']}`"])
    lines.extend(["", "## Phase Diagram And Conditions Context"])
    if payload["phase_context"]:
        for card in payload["phase_context"]:
            lines.append(f"### {card['formula']} ({card['chemical_system']})")
            lines.append(f"- Data status: {card['data_status']}")
            lines.append(f"- Notes: {card['stability_notes']}")
            for suggestion in card["suggestions"]:
                lines.append(f"- Suggestion: {suggestion}")
    else:
        lines.append("- No phase context cards available.")
    lines.extend(["", "## Scientific Caveats", payload["scientific_caveat"]])
    return "\n".join(lines) + "\n"


def _write_peaks_csv(path: Path, peaks: list[Peak]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["peak_id", "two_theta", "intensity", "prominence", "fwhm"])
        writer.writeheader()
        for peak in peaks:
            writer.writerow(_peak_dict(peak))


def _write_candidates_csv(path: Path, scores: list[CandidateScore]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "rank",
                "candidate_id",
                "formula",
                "score",
                "confidence_label",
                "matched_strong_intensity_fraction",
                "missing_strong_intensity_fraction",
                "missing_strong_penalty",
                "multi_peak_bonus",
                "source",
                "source_path",
                "space_group_symbol",
                "space_group_number",
                "equivalent_library_record_count",
                "equivalent_library_record_ids",
                "simulation_backend",
                "simulation_backend_version",
                "pattern_fingerprint",
                "simulation_state",
                "matched_count",
                "missing_strong_count",
                "warnings",
                "simulation_validation_status",
                "simulation_validation_reference",
            ],
        )
        writer.writeheader()
        for rank, score in enumerate(scores, 1):
            writer.writerow(
                {
                    "rank": rank,
                    "candidate_id": score.candidate.candidate_id,
                    "formula": score.candidate.formula_pretty,
                    "score": round(score.score, 4),
                    "confidence_label": score.score_components.get("confidence_label", ""),
                    "matched_strong_intensity_fraction": score.score_components.get("matched_strong_intensity_fraction", ""),
                    "missing_strong_intensity_fraction": score.score_components.get("missing_strong_intensity_fraction", ""),
                    "missing_strong_penalty": score.score_components.get("missing_strong_penalty", ""),
                    "multi_peak_bonus": score.score_components.get("multi_peak_bonus", ""),
                    "source": score.candidate.source,
                    "source_path": score.candidate.source_path,
                    "space_group_symbol": score.candidate.space_group_symbol,
                    "space_group_number": score.candidate.space_group_number,
                    "equivalent_library_record_count": len(score.candidate.phase_entry_ids),
                    "equivalent_library_record_ids": " | ".join(score.candidate.phase_entry_ids),
                    "simulation_backend": "" if score.candidate.simulated_pattern is None else score.candidate.simulated_pattern.engine_name,
                    "simulation_backend_version": "" if score.candidate.simulated_pattern is None else score.candidate.simulated_pattern.engine_version,
                    "pattern_fingerprint": "" if score.candidate.simulated_pattern is None else score.candidate.simulated_pattern.pattern_fingerprint,
                    "simulation_state": score.candidate.simulation_state,
                    "matched_count": len(score.matched_theory_peaks),
                    "missing_strong_count": len(score.missing_strong_theory_peaks),
                    "warnings": " | ".join(score.warnings),
                    "simulation_validation_status": score.candidate.simulation_validation.get("status", ""),
                    "simulation_validation_reference": score.candidate.simulation_validation.get("reference_path", ""),
                }
            )


def _write_assignments_csv(path: Path, peaks: list[Peak], scores: list[CandidateScore]) -> None:
    rows = []
    for peak in peaks:
        matches = []
        for score in scores:
            if any(m.experimental_peak.peak_id == peak.peak_id for m in score.matched_theory_peaks):
                matches.append(score.candidate.formula_pretty)
        rows.append({"peak_id": peak.peak_id, "two_theta": peak.two_theta, "intensity": peak.intensity, "candidate_matches": "; ".join(matches), "status": "explained" if matches else "unexplained"})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["peak_id", "two_theta", "intensity", "candidate_matches", "status"])
        writer.writeheader()
        writer.writerows(rows)


def _write_contribution_csv(path: Path, contribution: ContributionProxy | None) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["formula", "relative_contribution_proxy", "warning"])
        writer.writeheader()
        if contribution is None:
            return
        for formula, value in contribution.contributions.items():
            writer.writerow(
                {
                    "formula": formula,
                    "relative_contribution_proxy": round(value, 4),
                    "warning": contribution.warning,
                }
            )
