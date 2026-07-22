from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
from pathlib import Path
from typing import Any

import numpy as np


def _readonly_array(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=float).copy()
    array.setflags(write=False)
    return array


@dataclass(frozen=True)
class PhaseRefinementInput:
    phase_id: str
    formula: str
    cif_path: Path
    role: str
    initial_scale: float = 1.0
    cif_sha256: str = field(init=False)

    def __post_init__(self) -> None:
        path = Path(self.cif_path).resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Refinement CIF does not exist: {path}")
        if path.suffix.lower() != ".cif":
            raise ValueError(f"Multiphase refinement requires a CIF file: {path}")
        if self.role not in {"target", "impurity"}:
            raise ValueError("Phase role must be 'target' or 'impurity'.")
        if not self.phase_id.strip() or not self.formula.strip():
            raise ValueError("Phase id and formula must not be blank.")
        if not np.isfinite(self.initial_scale) or self.initial_scale <= 0.0:
            raise ValueError("Initial phase scale must be positive and finite.")
        object.__setattr__(self, "cif_path", path)
        object.__setattr__(self, "cif_sha256", hashlib.sha256(path.read_bytes()).hexdigest())


@dataclass(frozen=True)
class MultiphaseRefinementSettings:
    two_theta_min: float
    two_theta_max: float
    background_terms: int = 6
    refine_zero_shift: bool = True
    refine_profile: bool = True
    refine_lattice: bool = False
    radiation: str = "CuKalpha12"
    stability_spread_limit_percent: float = 5.0

    def __post_init__(self) -> None:
        if not (
            np.isfinite(self.two_theta_min)
            and np.isfinite(self.two_theta_max)
            and self.two_theta_min < self.two_theta_max
        ):
            raise ValueError("Multiphase refinement range must contain two finite increasing values.")
        if not 1 <= int(self.background_terms) <= 12:
            raise ValueError("RIETAN background_terms must be between 1 and 12.")
        if self.radiation not in {"CuKa", "CuKalpha12", "CuKalpha", "CuKα"}:
            raise ValueError("The first multiphase release supports Cu Kalpha radiation only.")
        if not np.isfinite(self.stability_spread_limit_percent) or self.stability_spread_limit_percent <= 0:
            raise ValueError("The stability spread limit must be positive and finite.")


@dataclass(frozen=True)
class PhaseRefinementResult:
    phase_id: str
    formula: str
    scale: float
    z: float
    molar_mass: float
    volume_angstrom3: float
    weight_percent: float
    contribution_intensity: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))
    reflection_two_theta_deg: tuple[float, ...] = ()
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "contribution_intensity", _readonly_array(self.contribution_intensity))
        object.__setattr__(self, "reflection_two_theta_deg", tuple(float(value) for value in self.reflection_two_theta_deg))
        object.__setattr__(self, "warnings", tuple(self.warnings))


@dataclass(frozen=True)
class MultiphaseRefinementResult:
    two_theta_deg: np.ndarray
    observed_intensity: np.ndarray
    calculated_intensity: np.ndarray
    residual_intensity: np.ndarray
    background_intensity: np.ndarray
    phases: tuple[PhaseRefinementResult, ...]
    rwp_percent: float
    rp_percent: float
    goodness_of_fit: float
    parameters: dict[str, float] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()
    provenance: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        arrays = tuple(
            _readonly_array(values)
            for values in (
                self.two_theta_deg,
                self.observed_intensity,
                self.calculated_intensity,
                self.residual_intensity,
                self.background_intensity,
            )
        )
        if not arrays[0].size or any(values.shape != arrays[0].shape for values in arrays[1:]):
            raise ValueError("All multiphase refinement profile arrays must have the same non-zero length.")
        if np.any(np.diff(arrays[0]) <= 0.0):
            raise ValueError("Multiphase refinement 2theta values must be strictly increasing.")
        for name, values in zip(
            (
                "two_theta_deg",
                "observed_intensity",
                "calculated_intensity",
                "residual_intensity",
                "background_intensity",
            ),
            arrays,
        ):
            object.__setattr__(self, name, values)
        object.__setattr__(self, "phases", tuple(self.phases))
        object.__setattr__(self, "warnings", tuple(self.warnings))

