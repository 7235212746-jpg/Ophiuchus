from __future__ import annotations

from dataclasses import dataclass
import hashlib

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


@dataclass(frozen=True)
class BackgroundEstimate:
    values: np.ndarray
    method: str
    smoothness: float
    asymmetry: float
    iterations: int
    step_deg: float
    fingerprint: str

    def __post_init__(self) -> None:
        values = np.array(self.values, dtype=float, copy=True)
        values.setflags(write=False)
        object.__setattr__(self, "values", values)

    @property
    def parameters(self) -> dict[str, float | int]:
        return {
            "smoothness": self.smoothness,
            "asymmetry": self.asymmetry,
            "iterations": self.iterations,
            "step_deg": self.step_deg,
        }


def estimate_xrd_background(
    x: np.ndarray,
    intensity: np.ndarray,
    *,
    smoothness: float = 1.0e6,
    asymmetry: float = 0.001,
    iterations: int = 12,
) -> BackgroundEstimate:
    """Estimate a fixed smooth XRD background with asymmetric least squares."""
    x_values = np.asarray(x, dtype=float)
    y_values = np.asarray(intensity, dtype=float)
    if x_values.ndim != 1 or y_values.ndim != 1 or x_values.shape != y_values.shape:
        raise ValueError("x and intensity must be one-dimensional arrays with the same shape.")
    if not np.all(np.isfinite(x_values)) or not np.all(np.isfinite(y_values)):
        raise ValueError("x and intensity must contain only finite values.")
    differences = np.diff(x_values)
    if np.any(differences <= 0.0):
        raise ValueError("x must be strictly increasing.")
    if x_values.size < 3:
        raise ValueError("At least three XRD points are required for background estimation.")
    if not np.isfinite(smoothness) or smoothness <= 0.0:
        raise ValueError("smoothness must be finite and positive.")
    if not np.isfinite(asymmetry) or not 0.0 < asymmetry < 1.0:
        raise ValueError("asymmetry must lie strictly between zero and one.")
    if int(iterations) < 1:
        raise ValueError("iterations must be at least one.")

    step_deg = float(np.median(differences))
    effective_smoothness = float(np.clip(smoothness * (0.02 / step_deg) ** 4, 1.0e4, 1.0e12))
    second_difference = sparse.diags(
        [1.0, -2.0, 1.0],
        [0, 1, 2],
        shape=(x_values.size - 2, x_values.size),
        format="csc",
    )
    penalty = effective_smoothness * (second_difference.T @ second_difference)
    weights = np.ones(x_values.size, dtype=float)
    background = np.array(y_values, copy=True)
    for _ in range(int(iterations)):
        system = sparse.diags(weights, format="csc") + penalty
        background = np.asarray(spsolve(system, weights * y_values), dtype=float)
        weights = np.where(y_values > background, asymmetry, 1.0 - asymmetry)

    background.setflags(write=False)
    fingerprint = hashlib.sha256(np.ascontiguousarray(background).tobytes()).hexdigest()
    return BackgroundEstimate(
        values=background,
        method="asls",
        smoothness=float(smoothness),
        asymmetry=float(asymmetry),
        iterations=int(iterations),
        step_deg=step_deg,
        fingerprint=fingerprint,
    )
