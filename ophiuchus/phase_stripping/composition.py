from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from ophiuchus.xrd.models import Candidate, Peak

from .session import PhasePreview, PhaseStrippingSession


@dataclass(frozen=True)
class PeakCompositionRow:
    label: str
    kind: str
    intensity: float
    explained_share_percent: float | None = None
    hkl: str = ""
    reflection_two_theta: float | None = None
    reflection_delta: float | None = None


@dataclass(frozen=True)
class PeakComposition:
    index: int
    two_theta: float
    experimental: float
    background: float
    corrected: float
    explained: float
    residual: float
    rows: tuple[PeakCompositionRow, ...]


def inspect_peak_composition(
    session: PhaseStrippingSession,
    candidates_by_id: Mapping[str, Candidate],
    two_theta: float,
    *,
    preview: PhasePreview | None = None,
) -> PeakComposition:
    """Describe the reconstruction at the experimental point nearest two_theta."""
    x = session.context.x
    index = int(np.argmin(np.abs(x - float(two_theta))))
    snapped = float(x[index])
    experimental = float(session.context.intensity[index])
    background = float(session.background_y[index])
    corrected = float(session.corrected_intensity[index])
    explained = float(session.fitted_total[index])
    residual = float(session.residual_y[index])
    rows: list[PeakCompositionRow] = [
        PeakCompositionRow("实验强度", "experimental", experimental),
        PeakCompositionRow("估计背景", "background", background),
        PeakCompositionRow("扣背景后", "corrected", corrected),
        PeakCompositionRow("已解释合计", "explained", explained),
        PeakCompositionRow("有符号残差", "residual", residual),
    ]

    for operation, contribution in zip(session.accepted_operations, session.accepted_contributions):
        fit = operation.phase_fit
        candidate = candidates_by_id.get(fit.candidate_id)
        intensity = float(contribution[index])
        rows.append(
            _phase_row(
                candidate=candidate,
                candidate_id=fit.candidate_id,
                kind="accepted_phase",
                intensity=intensity,
                explained=explained,
                snapped_two_theta=snapped,
                shift_deg=fit.shift_deg,
                preview=False,
            )
        )

    selected_preview = preview or session.current_preview
    if selected_preview is not None:
        fit = selected_preview.phase_fit
        candidate = candidates_by_id.get(fit.candidate_id)
        rows.append(
            _phase_row(
                candidate=candidate,
                candidate_id=fit.candidate_id,
                kind="preview_phase",
                intensity=float(selected_preview.contribution[index]),
                explained=explained,
                snapped_two_theta=snapped,
                shift_deg=fit.shift_deg,
                preview=True,
            )
        )

    return PeakComposition(
        index=index,
        two_theta=snapped,
        experimental=experimental,
        background=background,
        corrected=corrected,
        explained=explained,
        residual=residual,
        rows=tuple(rows),
    )


def _phase_row(
    *,
    candidate: Candidate | None,
    candidate_id: str,
    kind: str,
    intensity: float,
    explained: float,
    snapped_two_theta: float,
    shift_deg: float,
    preview: bool,
) -> PeakCompositionRow:
    label = candidate.formula_pretty if candidate is not None else candidate_id
    if preview:
        label = f"{label} (预览)"
    share = None if explained <= np.finfo(float).eps else 100.0 * intensity / explained
    nearest = _nearest_reflection(candidate, snapped_two_theta, shift_deg)
    return PeakCompositionRow(
        label=label,
        kind=kind,
        intensity=intensity,
        explained_share_percent=share,
        hkl="" if nearest is None else nearest[0].hkl,
        reflection_two_theta=None if nearest is None else nearest[1],
        reflection_delta=None if nearest is None else nearest[1] - snapped_two_theta,
    )


def _nearest_reflection(
    candidate: Candidate | None,
    two_theta: float,
    shift_deg: float,
) -> tuple[Peak, float] | None:
    if candidate is None or not candidate.theory_peaks:
        return None
    peak = min(candidate.theory_peaks, key=lambda item: abs(item.two_theta + shift_deg - two_theta))
    return peak, float(peak.two_theta + shift_deg)
