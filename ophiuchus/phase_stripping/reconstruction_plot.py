from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class ReconstructionLayer:
    label: str
    color: str
    values: np.ndarray
    stacked_lower: np.ndarray
    stacked_upper: np.ndarray
    lane_values: np.ndarray


@dataclass(frozen=True)
class ReconstructionView:
    x: np.ndarray
    original: np.ndarray
    layers: tuple[ReconstructionLayer, ...]
    phase_sum: np.ndarray
    proposed_phase_sum: np.ndarray
    residual: np.ndarray
    overflow: np.ndarray
    preview_layer: ReconstructionLayer | None = None


@dataclass(frozen=True)
class CursorBreakdown:
    index: int
    two_theta: float
    experimental: float
    phase_sum: float
    residual: float
    contributions: tuple[tuple[str, float], ...]


def build_reconstruction_view(
    x: np.ndarray,
    original: np.ndarray,
    contributions: Iterable[np.ndarray],
    *,
    labels: Sequence[str],
    colors: Sequence[str],
    preview: np.ndarray | None = None,
    preview_label: str = "Preview",
    preview_color: str = "#15a34a",
) -> ReconstructionView:
    x_values = _immutable_1d(x, "x")
    original_values = _immutable_1d(original, "original")
    if x_values.shape != original_values.shape:
        raise ValueError("x and original must have the same shape.")

    contribution_values = tuple(_immutable_1d(values, "contribution") for values in contributions)
    if len(labels) != len(contribution_values) or len(colors) != len(contribution_values):
        raise ValueError("contributions, labels and colors must have equal lengths.")
    if any(values.shape != x_values.shape for values in contribution_values):
        raise ValueError("Every contribution must have the same shape as x.")

    cumulative = np.zeros_like(original_values, dtype=float)
    layers: list[ReconstructionLayer] = []
    for label, color, values in zip(labels, colors, contribution_values):
        lower = _immutable(cumulative)
        cumulative = cumulative + values
        layers.append(
            ReconstructionLayer(
                label=str(label),
                color=str(color),
                values=values,
                stacked_lower=lower,
                stacked_upper=_immutable(cumulative),
                lane_values=_lane_values(values),
            )
        )

    phase_sum = _immutable(cumulative)
    residual = _immutable(original_values - phase_sum)
    overflow = _immutable(np.maximum(phase_sum - original_values, 0.0))
    preview_layer = None
    proposed = np.array(phase_sum, copy=True)
    if preview is not None:
        preview_values = _immutable_1d(preview, "preview")
        if preview_values.shape != x_values.shape:
            raise ValueError("preview must have the same shape as x.")
        preview_layer = ReconstructionLayer(
            label=str(preview_label),
            color=str(preview_color),
            values=preview_values,
            stacked_lower=phase_sum,
            stacked_upper=_immutable(phase_sum + preview_values),
            lane_values=_lane_values(preview_values),
        )
        proposed += preview_values

    return ReconstructionView(
        x=x_values,
        original=original_values,
        layers=tuple(layers),
        phase_sum=phase_sum,
        proposed_phase_sum=_immutable(proposed),
        residual=residual,
        overflow=overflow,
        preview_layer=preview_layer,
    )


def cursor_breakdown(view: ReconstructionView, x_value: float) -> CursorBreakdown:
    if not np.isfinite(x_value):
        raise ValueError("Cursor position must be finite.")
    index = int(np.argmin(np.abs(view.x - float(x_value))))
    return CursorBreakdown(
        index=index,
        two_theta=float(view.x[index]),
        experimental=float(view.original[index]),
        phase_sum=float(view.phase_sum[index]),
        residual=float(view.residual[index]),
        contributions=tuple((layer.label, float(layer.values[index])) for layer in view.layers),
    )


def _lane_values(values: np.ndarray) -> np.ndarray:
    maximum = float(np.max(np.abs(values))) if values.size else 0.0
    return _immutable(values / maximum if maximum > 0.0 else np.zeros_like(values, dtype=float))


def _immutable_1d(values: np.ndarray, name: str) -> np.ndarray:
    result = np.array(values, dtype=float, copy=True)
    if result.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional.")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values.")
    result.setflags(write=False)
    return result


def _immutable(values: np.ndarray) -> np.ndarray:
    result = np.array(values, dtype=float, copy=True)
    result.setflags(write=False)
    return result
