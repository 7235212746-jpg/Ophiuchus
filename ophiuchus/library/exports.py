from __future__ import annotations

import csv
import json
import re
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any


SCIENTIFIC_CAVEAT = (
    "This export is a screening-level XRD candidate analysis, not Rietveld refinement. "
    "Relative contribution proxy values are not Rietveld wt% or true phase fractions."
)


def safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9\u4e00-\u9fff]+", "_", value.strip()).strip("_")
    return token or "Ophi"


def next_export_path(
    folder: str | Path,
    project: str,
    sample_id: str,
    date_text: str,
    analysis_type: str,
    ext: str,
) -> Path:
    out = Path(folder)
    out.mkdir(parents=True, exist_ok=True)
    suffix = ext if ext.startswith(".") else f".{ext}"
    prefix = "_".join(safe_token(item) for item in [project, sample_id, date_text, analysis_type])
    version = 1
    while True:
        candidate = out / f"{prefix}_v{version:02d}{suffix}"
        if not candidate.exists():
            return candidate
        version += 1


def write_analysis_json(
    path: str | Path,
    project: str,
    sample_id: str,
    settings: dict[str, Any],
    selected_library_entry_ids: list[str],
    candidate_phases: list[dict[str, Any]],
    extra: dict[str, Any] | None = None,
) -> str:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "project": project,
        "sample_id": sample_id,
        "created_time": datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "settings": settings,
        "selected_library_entry_ids": selected_library_entry_ids,
        "candidate_phases": candidate_phases,
        "scientific_caveat": SCIENTIFIC_CAVEAT,
    }
    if extra:
        payload.update(extra)
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    return str(out)


def export_analysis_bundle(
    result: Any,
    export_folder: str | Path,
    project: str,
    sample_id: str,
    date_text: str | None = None,
) -> dict[str, str]:
    date_value = date_text or datetime.utcnow().strftime("%Y%m%d")
    folder = Path(export_folder)
    folder.mkdir(parents=True, exist_ok=True)

    result_outputs = getattr(result, "outputs", None) or {}
    plot_path = next_export_path(folder, project, sample_id, date_value, "xrd_candidate_screening_presentation", "png")
    source_plot = result_outputs.get("xrd_presentation_plot") or result_outputs.get("xrd_plot")
    if source_plot and Path(source_plot).exists():
        shutil.copy2(source_plot, plot_path)
    else:
        plot_path.write_bytes(b"")

    diagnostic_plot_path = next_export_path(folder, project, sample_id, date_value, "xrd_candidate_screening_diagnostic", "png")
    source_diagnostic = result_outputs.get("xrd_diagnostic_plot") or result_outputs.get("xrd_plot_clean")
    if source_diagnostic and Path(source_diagnostic).exists():
        shutil.copy2(source_diagnostic, diagnostic_plot_path)
    else:
        diagnostic_plot_path.write_bytes(b"")

    peak_table_path = next_export_path(folder, project, sample_id, date_value, "peak_table", "csv")
    _write_peak_table(peak_table_path, result)

    selected_ids = [score.candidate.candidate_id for score in result.explanation.selected_candidates]
    candidate_phases = [
        {
            "candidate_id": score.candidate.candidate_id,
            "formula": score.candidate.formula_pretty,
            "source": score.candidate.source,
            "score": round(score.score, 4),
            "matched_peak_count": len(score.matched_theory_peaks),
            "missing_strong_peak_count": len(score.missing_strong_theory_peaks),
        }
        for score in result.top_scores
    ]
    report_path = next_export_path(folder, project, sample_id, date_value, "analysis_report", "json")
    write_analysis_json(
        report_path,
        project=project,
        sample_id=sample_id,
        settings={"export_date": date_value},
        selected_library_entry_ids=selected_ids,
        candidate_phases=candidate_phases,
        extra={"warnings": result.warnings, "source_outputs": result.outputs},
    )
    return {
        "xrd_presentation": str(plot_path),
        "xrd_diagnostic": str(diagnostic_plot_path),
        "peak_table": str(peak_table_path),
        "analysis_report": str(report_path),
    }


def _write_peak_table(path: Path, result: Any) -> None:
    fieldnames = [
        "experimental_peak_id",
        "experimental_two_theta",
        "experimental_intensity",
        "candidate_phase",
        "candidate_source",
        "calculated_two_theta",
        "delta_two_theta",
        "hkl",
        "theoretical_relative_intensity",
        "assignment_status",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        assigned = set()
        for score in result.top_scores:
            for match in score.matched_theory_peaks:
                assigned.add(match.experimental_peak.peak_id)
                writer.writerow(
                    {
                        "experimental_peak_id": match.experimental_peak.peak_id,
                        "experimental_two_theta": round(match.experimental_peak.two_theta, 4),
                        "experimental_intensity": round(match.experimental_peak.intensity, 4),
                        "candidate_phase": score.candidate.formula_pretty,
                        "candidate_source": score.candidate.source,
                        "calculated_two_theta": round(match.theory_peak.two_theta, 4),
                        "delta_two_theta": round(match.delta, 4),
                        "hkl": getattr(match.theory_peak, "hkl", "") or "",
                        "theoretical_relative_intensity": round(match.theory_peak.intensity, 4),
                        "assignment_status": "matched",
                    }
                )
        for peak in result.experimental_peaks:
            if peak.peak_id in assigned:
                continue
            writer.writerow(
                {
                    "experimental_peak_id": peak.peak_id,
                    "experimental_two_theta": round(peak.two_theta, 4),
                    "experimental_intensity": round(peak.intensity, 4),
                    "candidate_phase": "",
                    "candidate_source": "",
                    "calculated_two_theta": "",
                    "delta_two_theta": "",
                    "hkl": "",
                    "theoretical_relative_intensity": "",
                    "assignment_status": "unresolved",
                }
            )
