from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from ophiuchus.xrd.models import Candidate

from .models import FitBounds, PhaseFit
from .profile import project_candidate_profile


class PhaseContributionFitter:
    """Fit one complete candidate profile to an experimental residual."""

    negative_penalty_weight = 4.0

    def fit(
        self,
        residual_y: np.ndarray,
        x: np.ndarray,
        candidate: Candidate,
        bounds: FitBounds,
        fit_mask: np.ndarray | None = None,
    ) -> PhaseFit:
        x_values = np.asarray(x, dtype=float)
        residual_values = np.asarray(residual_y, dtype=float)
        if x_values.ndim != 1 or residual_values.ndim != 1:
            raise ValueError("Experimental x and residual arrays must be one-dimensional.")
        if x_values.shape != residual_values.shape:
            raise ValueError("Experimental x and residual arrays must have the same length.")
        if x_values.size == 0:
            raise ValueError("Experimental x and residual arrays must not be empty.")
        if not np.all(np.isfinite(residual_values)):
            raise ValueError("Residual values must be finite.")

        selected = self._fit_mask(fit_mask, x_values.shape)
        lower, upper = self._parameter_bounds(bounds)
        initial = self._initial_parameters(x_values, residual_values, candidate, lower, upper, selected)

        def objective(parameters: np.ndarray) -> np.ndarray:
            shift_deg, sigma_deg, scale = parameters
            contribution = scale * project_candidate_profile(x_values, candidate, shift_deg, sigma_deg)
            updated_residual = residual_values - contribution
            selected_residual = updated_residual[selected]
            newly_negative = np.maximum(-selected_residual, 0.0) - np.maximum(-residual_values[selected], 0.0)
            return np.concatenate((selected_residual, self.negative_penalty_weight * newly_negative))

        result = least_squares(objective, initial, bounds=(lower, upper), method="dogbox")
        shift_deg, sigma_deg, scale = result.x
        return PhaseFit(
            candidate_id=candidate.candidate_id,
            shift_deg=float(shift_deg),
            sigma_deg=float(sigma_deg),
            scale=float(scale),
        )

    @staticmethod
    def _fit_mask(fit_mask: np.ndarray | None, shape: tuple[int, ...]) -> np.ndarray:
        if fit_mask is None:
            return np.ones(shape, dtype=bool)
        selected = np.asarray(fit_mask, dtype=bool)
        if selected.shape != shape:
            raise ValueError("Fit mask must have the same shape as the experimental grid.")
        if not np.any(selected):
            raise ValueError("Fit mask must select at least one point.")
        return selected

    @staticmethod
    def _parameter_bounds(bounds: FitBounds) -> tuple[np.ndarray, np.ndarray]:
        lower = np.array((bounds.shift_deg[0], bounds.sigma_deg[0], bounds.scale[0]), dtype=float)
        upper = np.array((bounds.shift_deg[1], bounds.sigma_deg[1], bounds.scale[1]), dtype=float)
        if not np.all(np.isfinite(lower)) or np.isnan(upper).any():
            raise ValueError("Fit bounds must not contain NaN values.")
        if lower[2] < 0.0:
            raise ValueError("Phase scale lower bound must be non-negative.")
        if np.any(lower >= upper):
            raise ValueError("Each fit lower bound must be smaller than its upper bound.")
        return lower, upper

    @staticmethod
    def _initial_parameters(
        x: np.ndarray,
        residual_y: np.ndarray,
        candidate: Candidate,
        lower: np.ndarray,
        upper: np.ndarray,
        selected: np.ndarray,
    ) -> np.ndarray:
        shift_deg = (lower[0] + upper[0]) / 2.0
        sigma_deg = (lower[1] + upper[1]) / 2.0
        profile = project_candidate_profile(x, candidate, shift_deg, sigma_deg)
        selected_profile = profile[selected]
        scale = float(np.dot(residual_y[selected], selected_profile) / np.dot(selected_profile, selected_profile))
        scale = float(np.clip(scale, lower[2], upper[2]))
        return np.array((shift_deg, sigma_deg, scale), dtype=float)
