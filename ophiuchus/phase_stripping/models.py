from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class AnalysisContext:
    x: np.ndarray
    intensity: np.ndarray
    radiation: str
    wavelength_angstrom: float
    two_theta_range: tuple[float, float]
    tolerance_deg: float
    source_path: str
    source_fingerprint: str
    data_fingerprint: str

    def __post_init__(self) -> None:
        x = np.array(self.x, copy=True)
        intensity = np.array(self.intensity, copy=True)
        if x.ndim != 1 or intensity.ndim != 1:
            raise ValueError("Experimental x and intensity arrays must be one-dimensional.")
        if x.shape != intensity.shape:
            raise ValueError("Experimental x and intensity arrays must have the same length.")
        x.setflags(write=False)
        intensity.setflags(write=False)
        object.__setattr__(self, "x", x)
        object.__setattr__(self, "intensity", intensity)


@dataclass(frozen=True)
class FitBounds:
    shift_deg: tuple[float, float]
    sigma_deg: tuple[float, float]
    scale: tuple[float, float]


@dataclass(frozen=True)
class InstrumentSettings:
    radiation: str
    wavelength_angstrom: float
    two_theta_range: tuple[float, float]
    tolerance_deg: float


@dataclass(frozen=True)
class PhaseFit:
    candidate_id: str
    shift_deg: float
    sigma_deg: float
    scale: float


@dataclass(frozen=True)
class PhaseOperation:
    operation_id: str
    phase_fit: PhaseFit
    residual_fingerprint: str
