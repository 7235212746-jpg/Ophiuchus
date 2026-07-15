from __future__ import annotations

import csv
import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any

from ophiuchus.xrd.models import CandidateScore, Peak


def write_analysis_audit_folder(
    out_dir: str | Path,
    input_summary: dict[str, Any],
    experimental_peaks: list[Peak],
    top_scores: list[CandidateScore],
    all_candidates: list[Any],
    outputs: dict[str, str],
    library_entries: list[Any] | None = None,
) -> str:
    audit = Path(out_dir) / f"audit_{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}"
    audit.mkdir(parents=True, exist_ok=True)
    (audit / "input_summary.json").write_text(json.dumps(input_summary, indent=2, ensure_ascii=False), encoding="utf-8")
    _write_experimental_peaks(audit / "experimental_peaks.csv", experimental_peaks)
    _write_candidate_ranking(audit / "candidate_phase_ranking.csv", top_scores)
    _write_peak_match_table(audit / "peak_match_table.csv", top_scores)
    _write_unexplained(audit / "unexplained_peaks.csv", experimental_peaks, top_scores)
    _write_simulated_peaks(audit / "simulated_peaks_all.csv", top_scores)
    _write_simulation_validation(audit / "simulation_validation.csv", top_scores)
    if library_entries is not None:
        _write_loaded_structures(audit / "loaded_structures.csv", library_entries)
        _write_skipped_structures(audit / "skipped_structures.csv", input_summary, library_entries)
    if outputs.get("json") and Path(outputs["json"]).exists():
        shutil.copy2(outputs["json"], audit / "analysis_report.json")
    if outputs.get("candidate_usage_summary") and Path(outputs["candidate_usage_summary"]).exists():
        shutil.copy2(outputs["candidate_usage_summary"], audit / "candidate_usage_summary.csv")
    if outputs.get("xrd_plot") and Path(outputs["xrd_plot"]).exists():
        shutil.copy2(outputs["xrd_plot"], audit / "xrd_plot.png")
    if outputs.get("xrd_presentation_plot") and Path(outputs["xrd_presentation_plot"]).exists():
        shutil.copy2(outputs["xrd_presentation_plot"], audit / "xrd_presentation.png")
    if outputs.get("xrd_diagnostic_plot") and Path(outputs["xrd_diagnostic_plot"]).exists():
        shutil.copy2(outputs["xrd_diagnostic_plot"], audit / "xrd_diagnostic.png")
    return str(audit)


def _write_experimental_peaks(path: Path, peaks: list[Peak]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["peak_id", "two_theta", "intensity"])
        writer.writeheader()
        for peak in peaks:
            writer.writerow({"peak_id": peak.peak_id, "two_theta": peak.two_theta, "intensity": peak.intensity})


def _write_candidate_ranking(path: Path, scores: list[CandidateScore]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "rank",
                "candidate_id",
                "formula",
                "source",
                "source_path",
                "score",
                "confidence_label",
                "matched_strong_intensity_fraction",
                "missing_strong_penalty",
                "matched_count",
                "missing_strong_count",
                "simulation_validation_status",
                "simulation_validation_reference",
                "role",
            ],
        )
        writer.writeheader()
        for rank, score in enumerate(scores, 1):
            role = "likely_main_phase" if rank == 1 else "possible_impurity" if score.score >= 0.35 else "weak_candidate"
            writer.writerow(
                {
                    "rank": rank,
                    "candidate_id": score.candidate.candidate_id,
                    "formula": score.candidate.formula_pretty,
                    "source": score.candidate.source,
                    "source_path": score.candidate.source_path,
                    "score": round(score.score, 4),
                    "confidence_label": score.score_components.get("confidence_label", ""),
                    "matched_strong_intensity_fraction": score.score_components.get("matched_strong_intensity_fraction", ""),
                    "missing_strong_penalty": score.score_components.get("missing_strong_penalty", ""),
                    "matched_count": len(score.matched_theory_peaks),
                    "missing_strong_count": len(score.missing_strong_theory_peaks),
                    "simulation_validation_status": score.candidate.simulation_validation.get("status", ""),
                    "simulation_validation_reference": score.candidate.simulation_validation.get("reference_path", ""),
                    "role": role,
                }
            )


def _write_peak_match_table(path: Path, scores: list[CandidateScore]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["experimental_peak_id", "experimental_two_theta", "candidate_id", "formula", "source", "calculated_two_theta", "delta_two_theta", "hkl", "calculated_relative_intensity", "confidence_level"])
        writer.writeheader()
        for score in scores:
            for match in score.matched_theory_peaks:
                confidence = "likely" if score.score >= 0.7 else "possible" if score.score >= 0.35 else "weak"
                writer.writerow(
                    {
                        "experimental_peak_id": match.experimental_peak.peak_id,
                        "experimental_two_theta": round(match.experimental_peak.two_theta, 4),
                        "candidate_id": score.candidate.candidate_id,
                        "formula": score.candidate.formula_pretty,
                        "source": score.candidate.source,
                        "calculated_two_theta": round(match.theory_peak.two_theta, 4),
                        "delta_two_theta": round(match.delta, 4),
                        "hkl": match.theory_peak.hkl,
                        "calculated_relative_intensity": round(match.theory_peak.intensity, 4),
                        "confidence_level": confidence,
                    }
                )


def _write_unexplained(path: Path, experimental_peaks: list[Peak], scores: list[CandidateScore]) -> None:
    explained = {m.experimental_peak.peak_id for score in scores for m in score.matched_theory_peaks}
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["peak_id", "two_theta", "intensity"])
        writer.writeheader()
        for peak in experimental_peaks:
            if peak.peak_id not in explained:
                writer.writerow({"peak_id": peak.peak_id, "two_theta": peak.two_theta, "intensity": peak.intensity})


def _write_simulated_peaks(path: Path, scores: list[CandidateScore]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["candidate_id", "formula", "source", "peak_id", "two_theta", "relative_intensity", "hkl", "d_spacing"])
        writer.writeheader()
        seen: set[str] = set()
        for score in scores:
            if score.candidate.candidate_id in seen:
                continue
            seen.add(score.candidate.candidate_id)
            for peak in score.candidate.theory_peaks:
                writer.writerow(
                    {
                        "candidate_id": score.candidate.candidate_id,
                        "formula": score.candidate.formula_pretty,
                        "source": score.candidate.source,
                        "peak_id": peak.peak_id,
                        "two_theta": peak.two_theta,
                        "relative_intensity": peak.intensity,
                        "hkl": peak.hkl,
                        "d_spacing": peak.d_spacing,
                    }
                )


def _write_simulation_validation(path: Path, scores: list[CandidateScore]) -> None:
    fields = [
        "candidate_id",
        "formula",
        "source",
        "status",
        "reference_path",
        "matched_count",
        "reference_strong_count",
        "calculated_strong_count",
        "missing_strong_count",
        "extra_strong_count",
        "max_abs_delta_2theta",
        "max_strong_intensity_ratio_error",
        "reason",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        seen: set[str] = set()
        for score in scores:
            if score.candidate.candidate_id in seen:
                continue
            seen.add(score.candidate.candidate_id)
            validation = score.candidate.simulation_validation or {}
            writer.writerow(
                {
                    "candidate_id": score.candidate.candidate_id,
                    "formula": score.candidate.formula_pretty,
                    "source": score.candidate.source,
                    "status": validation.get("status", ""),
                    "reference_path": validation.get("reference_path", ""),
                    "matched_count": validation.get("matched_count", ""),
                    "reference_strong_count": validation.get("reference_strong_count", ""),
                    "calculated_strong_count": validation.get("calculated_strong_count", ""),
                    "missing_strong_count": validation.get("missing_strong_count", ""),
                    "extra_strong_count": validation.get("extra_strong_count", ""),
                    "max_abs_delta_2theta": validation.get("max_abs_delta_2theta", ""),
                    "max_strong_intensity_ratio_error": validation.get("max_strong_intensity_ratio_error", ""),
                    "reason": validation.get("reason", ""),
                }
            )


def _write_loaded_structures(path: Path, entries: list[Any]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["internal_id", "formula", "source", "source_id", "elements", "cached_file_path", "metadata_path", "enabled"])
        writer.writeheader()
        for entry in entries:
            writer.writerow(
                {
                    "internal_id": entry.internal_id,
                    "formula": entry.formula,
                    "source": entry.source,
                    "source_id": entry.source_id,
                    "elements": " ".join(entry.elements),
                    "cached_file_path": entry.cached_file_path,
                    "metadata_path": entry.local_metadata_path,
                    "enabled": entry.enabled_for_matching,
                }
            )


def _write_skipped_structures(path: Path, input_summary: dict[str, Any], entries: list[Any]) -> None:
    allowed = set(input_summary.get("allowed_elements") or [])
    used_ids = set(input_summary.get("library_entry_ids") or [])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["internal_id", "formula", "source", "reason"])
        writer.writeheader()
        for entry in entries:
            reason = ""
            if not entry.enabled_for_matching:
                reason = "disabled"
            elif allowed and not set(entry.elements).issubset(allowed):
                reason = "element_filter"
            elif entry.internal_id not in used_ids:
                reason = "not_used_no_valid_cache_or_filter"
            if reason:
                writer.writerow({"internal_id": entry.internal_id, "formula": entry.formula, "source": entry.source, "reason": reason})
