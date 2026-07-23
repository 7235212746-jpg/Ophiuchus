from __future__ import annotations

import csv
import json
from pathlib import Path
import shutil
from typing import Any
import uuid

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from .session import PhaseStrippingSession


def export_phase_stripping_session(session: PhaseStrippingSession, folder: str | Path) -> dict[str, str]:
    """Export the accepted session state without recomputing phase profiles."""
    target = Path(folder)
    restore_empty_target = False
    if target.exists():
        if not target.is_dir() or any(target.iterdir()):
            raise FileExistsError(f"Export destination already exists and is not empty: {target}")
        target.rmdir()
        restore_empty_target = True
    payload = _session_payload(session)
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = target.parent / f".{target.name}.ophi-export-{uuid.uuid4().hex}"
    try:
        staging.mkdir()
        _write_csv(session, staging / "phase_stripping_residual.csv")
        (staging / "phase_stripping_session.json").write_text(
            json.dumps(payload, indent=2, allow_nan=False), encoding="utf-8"
        )
        _write_plot(session, staging / "phase_stripping_residual.png")
        staging.replace(target)
    except Exception:
        if staging.exists():
            shutil.rmtree(staging)
        if restore_empty_target and not target.exists():
            target.mkdir()
        raise
    csv_path = target / "phase_stripping_residual.csv"
    json_path = target / "phase_stripping_session.json"
    png_path = target / "phase_stripping_residual.png"
    return {"csv": str(csv_path), "json": str(json_path), "png": str(png_path)}


def _write_csv(session: PhaseStrippingSession, path: Path) -> None:
    contribution_columns = [f"contribution_{operation.operation_id.replace('-', '_')}" for operation in session.accepted_operations]
    headers = [
        "two_theta",
        "original_y",
        "background_y",
        "corrected_y",
        "fitted_total",
        "reconstructed_y",
        "residual_y",
        *contribution_columns,
    ]
    contributions = session.accepted_contributions
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        rows = zip(
            session.context.x,
            session.context.intensity,
            session.background_y,
            session.corrected_intensity,
            session.fitted_total,
            session.reconstructed_y,
            session.residual_y,
        )
        for index, values in enumerate(rows):
            writer.writerow([*values, *(contribution[index] for contribution in contributions)])


def _session_payload(session: PhaseStrippingSession) -> dict[str, Any]:
    serialized = session.to_dict()
    accepted_operations = []
    for item in serialized["accepted"]:
        candidate = item["candidate"]
        required = ("candidate_id", "cif_path", "structure_hash", "pattern_fingerprint")
        if not candidate or any(not candidate.get(field) for field in required):
            raise ValueError(
                "Cannot export a reproducible phase-stripping session: candidate provenance is incomplete."
            )
        accepted_operations.append(
            {
                "operation_id": item["operation_id"],
                "candidate": candidate,
                "fit": item["phase_fit"],
                "residual_fingerprint": item["residual_fingerprint"],
                "contribution": item["contribution"],
                "warnings": item["warnings"],
                "negative_area_before": item["negative_area_before"],
                "negative_area_after": item["negative_area_after"],
            }
        )
    context = serialized["context"]
    return {
        "format": "ophiuchus.phase_stripping.session.v2",
        "instrument_settings": {
            "radiation": context["radiation"],
            "wavelength_angstrom": context["wavelength_angstrom"],
            "two_theta_range": context["two_theta_range"],
            "tolerance_deg": context["tolerance_deg"],
        },
        "fit_bounds": serialized["bounds"],
        "background_model": serialized["background"],
        "accepted_operations": accepted_operations,
        "session": serialized,
    }


def _write_plot(session: PhaseStrippingSession, path: Path) -> None:
    figure = Figure(figsize=(12, 7), dpi=160, facecolor="#ffffff")
    canvas = FigureCanvasAgg(figure)
    pattern_axis, residual_axis = figure.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    x = session.context.x
    pattern_axis.set_facecolor("#ffffff")
    residual_axis.set_facecolor("#ffffff")
    pattern_axis.plot(x, session.context.intensity, color="#111827", linewidth=1.15, label="Original")
    pattern_axis.plot(x, session.background_y, color="#d48a17", linewidth=1.0, label="Estimated background")
    pattern_axis.plot(x, session.reconstructed_y, color="#3974d8", linewidth=1.1, label="Background + phase model")
    colors = ("#2a9d8f", "#8f5bd7", "#d48a17", "#0086b3", "#b84a62")
    for index, (operation, contribution) in enumerate(zip(session.accepted_operations, session.accepted_contributions)):
        pattern_axis.plot(x, session.background_y + contribution, color=colors[index % len(colors)], linewidth=0.8, alpha=0.8, label=operation.phase_fit.candidate_id)
    residual_axis.plot(x, session.residual_y, color="#c0392b", linewidth=1.0, label="Signed residual after background subtraction")
    residual_axis.axhline(0.0, color="#7b8798", linewidth=0.7)
    residual_axis.set_xlabel("2theta (deg)")
    pattern_axis.set_ylabel("Intensity")
    residual_axis.set_ylabel("Residual")
    for axis in (pattern_axis, residual_axis):
        axis.grid(color="#e6edf5", linewidth=0.6, alpha=0.8)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    pattern_axis.legend(loc="best", frameon=False, ncols=2)
    residual_axis.legend(loc="best", frameon=False)
    figure.tight_layout()
    try:
        canvas.print_png(path)
    finally:
        figure.clear()
