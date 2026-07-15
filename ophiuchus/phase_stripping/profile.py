from __future__ import annotations

import numpy as np

from ophiuchus.xrd.models import Candidate


def project_candidate_profile(
    x: np.ndarray,
    candidate: Candidate,
    shift_deg: float,
    sigma_deg: float,
) -> np.ndarray:
    """Project a candidate's canonical peaks onto an experimental x grid."""
    x_values = np.asarray(x, dtype=float)
    if x_values.ndim != 1 or x_values.size == 0:
        raise ValueError("Experimental x grid must be a non-empty one-dimensional array.")
    if not np.all(np.diff(x_values) > 0.0):
        raise ValueError("Experimental x grid must be strictly increasing.")
    if not np.isfinite(sigma_deg) or sigma_deg <= 0.0:
        raise ValueError("Gaussian sigma must be positive.")

    peaks = candidate.theory_peaks
    if not peaks:
        raise ValueError("Candidate must provide canonical theory peaks.")

    profile = np.zeros_like(x_values, dtype=float)
    for peak in peaks:
        center = peak.two_theta + shift_deg
        profile += peak.intensity * np.exp(-0.5 * ((x_values - center) / sigma_deg) ** 2)

    maximum = float(profile.max())
    if maximum > 0.0:
        profile /= maximum
    return profile
