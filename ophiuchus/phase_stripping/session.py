from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.signal import find_peaks

from ophiuchus.xrd.models import Candidate

from .fitting import PhaseContributionFitter
from .models import AnalysisContext, FitBounds, PhaseFit, PhaseOperation
from .profile import project_candidate_profile


DEFAULT_FIT_BOUNDS = FitBounds(
    shift_deg=(-0.25, 0.25),
    sigma_deg=(0.035, 0.25),
    scale=(0.0, float("inf")),
)


@dataclass(frozen=True)
class PhasePreview:
    phase_fit: PhaseFit
    contribution: np.ndarray
    warnings: tuple[str, ...]
    candidate_identity: _CandidateIdentity | None = None

    def __post_init__(self) -> None:
        contribution = np.array(self.contribution, dtype=float, copy=True)
        contribution.setflags(write=False)
        object.__setattr__(self, "contribution", contribution)


@dataclass(frozen=True)
class _StoredOperation:
    operation: PhaseOperation
    contribution: np.ndarray
    warnings: tuple[str, ...]
    negative_area_before: float
    negative_area_after: float
    candidate_identity: _CandidateIdentity | None = None

    def __post_init__(self) -> None:
        contribution = np.array(self.contribution, dtype=float, copy=True)
        contribution.setflags(write=False)
        object.__setattr__(self, "contribution", contribution)


@dataclass(frozen=True)
class _CandidateIdentity:
    candidate_id: str
    formula_pretty: str
    cif_path: str
    structure_hash: str | None
    pattern_fingerprint: str | None

    @classmethod
    def from_candidate(cls, candidate: Candidate) -> _CandidateIdentity:
        validation = candidate.simulation_validation
        fingerprint = validation.get("pattern_fingerprint") if isinstance(validation, dict) else None
        pattern = candidate.simulated_pattern
        if fingerprint is None and pattern is not None:
            fingerprint = getattr(pattern, "pattern_fingerprint", None)
        if fingerprint is None and candidate.theory_peaks:
            digest = hashlib.sha256()
            for peak in sorted(candidate.theory_peaks, key=lambda item: (item.two_theta, item.intensity, item.peak_id)):
                digest.update(f"{peak.two_theta:.8f}|{peak.intensity:.8f}|{peak.hkl}\n".encode("utf-8"))
            fingerprint = digest.hexdigest()
        structure_hash_value = candidate.structure_hash
        source_path = Path(candidate.source_path)
        if structure_hash_value is None and source_path.is_file():
            structure_hash_value = hashlib.sha256(source_path.read_bytes()).hexdigest()
        return cls(
            candidate_id=candidate.candidate_id,
            formula_pretty=candidate.formula_pretty,
            cif_path=candidate.source_path,
            structure_hash=structure_hash_value,
            pattern_fingerprint=None if fingerprint is None else str(fingerprint),
        )

    def to_dict(self) -> dict[str, str | None]:
        return {
            "candidate_id": self.candidate_id,
            "formula_pretty": self.formula_pretty,
            "cif_path": self.cif_path,
            "structure_hash": self.structure_hash,
            "pattern_fingerprint": self.pattern_fingerprint,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> _CandidateIdentity:
        return cls(
            candidate_id=str(data["candidate_id"]),
            formula_pretty=str(data.get("formula_pretty", "")),
            cif_path=str(data.get("cif_path", "")),
            structure_hash=None if data.get("structure_hash") is None else str(data["structure_hash"]),
            pattern_fingerprint=None
            if data.get("pattern_fingerprint") is None
            else str(data["pattern_fingerprint"]),
        )


class PhaseStrippingSession:
    """Non-destructive residual history built from immutable contributions."""

    def __init__(
        self,
        context: AnalysisContext,
        bounds: FitBounds | None = None,
        fitter: PhaseContributionFitter | None = None,
        *,
        background_y: np.ndarray | None = None,
        background_method: str = "none",
        background_parameters: dict[str, float | int] | None = None,
    ) -> None:
        self.context = context
        self.bounds = bounds or DEFAULT_FIT_BOUNDS
        self._fitter = fitter or PhaseContributionFitter()
        background = np.zeros_like(context.intensity, dtype=float) if background_y is None else np.asarray(background_y, dtype=float)
        if background.shape != context.intensity.shape:
            raise ValueError("background_y must have the same shape as the experimental intensity.")
        if not np.all(np.isfinite(background)):
            raise ValueError("background_y must contain only finite values.")
        self._background_y = np.array(background, dtype=float, copy=True)
        self._background_y.setflags(write=False)
        self.background_method = str(background_method)
        self._background_parameters = dict(background_parameters or {})
        self._accepted: list[_StoredOperation] = []
        self._redo: list[_StoredOperation] = []
        self._excluded_candidate_ids: set[str] = set()
        self._preview: PhasePreview | None = None
        self._next_operation_number = 1

    @property
    def accepted_operations(self) -> tuple[PhaseOperation, ...]:
        return tuple(stored.operation for stored in self._accepted)

    @property
    def accepted_contributions(self) -> tuple[np.ndarray, ...]:
        contributions: list[np.ndarray] = []
        for stored in self._accepted:
            contribution = np.array(stored.contribution, copy=True)
            contribution.setflags(write=False)
            contributions.append(contribution)
        return tuple(contributions)

    @property
    def accepted_peak_positions(self) -> tuple[float, ...]:
        positions: list[float] = []
        for contribution in self.accepted_contributions:
            peak_indices, _ = find_peaks(contribution)
            if contribution.size >= 2 and contribution[0] > contribution[1]:
                peak_indices = np.insert(peak_indices, 0, 0)
            if contribution.size >= 2 and contribution[-1] > contribution[-2]:
                peak_indices = np.append(peak_indices, contribution.size - 1)
            positions.extend(float(self.context.x[index]) for index in peak_indices)
        return tuple(positions)

    @property
    def excluded_candidate_ids(self) -> tuple[str, ...]:
        return tuple(sorted(self._excluded_candidate_ids))

    @property
    def current_preview(self) -> PhasePreview | None:
        return self._preview

    @property
    def background_y(self) -> np.ndarray:
        values = np.array(self._background_y, copy=True)
        values.setflags(write=False)
        return values

    @property
    def background_parameters(self) -> dict[str, float | int]:
        return dict(self._background_parameters)

    @property
    def corrected_intensity(self) -> np.ndarray:
        values = np.asarray(self.context.intensity, dtype=float) - self._background_y
        values.setflags(write=False)
        return values

    @property
    def fitted_total(self) -> np.ndarray:
        total = np.zeros_like(self.context.intensity, dtype=float)
        for stored in self._accepted:
            total += stored.contribution
        return total

    @property
    def residual_y(self) -> np.ndarray:
        return np.array(self.corrected_intensity, copy=True) - self.fitted_total

    @property
    def reconstructed_y(self) -> np.ndarray:
        return np.array(self._background_y, copy=True) + self.fitted_total

    def preview(
        self,
        candidate: Candidate,
        bounds: FitBounds | None = None,
        fit_mask: np.ndarray | None = None,
    ) -> PhasePreview:
        if candidate.candidate_id in self._excluded_candidate_ids:
            raise ValueError(f"Candidate '{candidate.candidate_id}' is excluded from this session.")
        phase_fit = self._fitter.fit(self.residual_y, self.context.x, candidate, bounds or self.bounds, fit_mask)
        contribution = phase_fit.scale * project_candidate_profile(
            self.context.x,
            candidate,
            phase_fit.shift_deg,
            phase_fit.sigma_deg,
        )
        warnings = (
            *self._negative_area_warnings(self.residual_y, contribution),
            *self._fit_boundary_warnings(phase_fit, bounds or self.bounds),
        )
        preview = PhasePreview(
            phase_fit=phase_fit,
            contribution=contribution,
            warnings=warnings,
            candidate_identity=_CandidateIdentity.from_candidate(candidate),
        )
        self._preview = preview
        return preview

    def preview_with_parameters(
        self,
        candidate: Candidate,
        *,
        scale: float,
        shift_deg: float,
        sigma_deg: float,
    ) -> PhasePreview:
        if candidate.candidate_id in self._excluded_candidate_ids:
            raise ValueError(f"Candidate '{candidate.candidate_id}' is excluded from this session.")
        values = {"scale": float(scale), "shift": float(shift_deg), "sigma": float(sigma_deg)}
        if any(not np.isfinite(value) for value in values.values()):
            raise ValueError("Manual preview parameters must be finite.")
        self._require_bounded("scale", values["scale"], self.bounds.scale)
        self._require_bounded("shift", values["shift"], self.bounds.shift_deg)
        self._require_bounded("sigma", values["sigma"], self.bounds.sigma_deg)
        phase_fit = PhaseFit(
            candidate_id=candidate.candidate_id,
            scale=values["scale"],
            shift_deg=values["shift"],
            sigma_deg=values["sigma"],
        )
        contribution = phase_fit.scale * project_candidate_profile(
            self.context.x,
            candidate,
            phase_fit.shift_deg,
            phase_fit.sigma_deg,
        )
        preview = PhasePreview(
            phase_fit=phase_fit,
            contribution=contribution,
            warnings=self._negative_area_warnings(self.residual_y, contribution),
            candidate_identity=_CandidateIdentity.from_candidate(candidate),
        )
        self._preview = preview
        return preview

    def cancel_preview(self) -> bool:
        had_preview = self._preview is not None
        self._preview = None
        return had_preview

    def accept_preview(self, preview: PhasePreview | None = None) -> PhaseOperation:
        selected_preview = preview or self._preview
        if selected_preview is None:
            raise ValueError("A phase preview is required before acceptance.")
        residual_before = self.residual_y
        residual_after = residual_before - selected_preview.contribution
        operation = PhaseOperation(
            operation_id=f"operation-{self._next_operation_number}",
            phase_fit=selected_preview.phase_fit,
            residual_fingerprint=self._array_fingerprint(residual_after),
        )
        stored = _StoredOperation(
            operation=operation,
            contribution=selected_preview.contribution,
            warnings=selected_preview.warnings,
            negative_area_before=self._negative_area(residual_before),
            negative_area_after=self._negative_area(residual_after),
            candidate_identity=selected_preview.candidate_identity,
        )
        self._accepted.append(stored)
        self._redo.clear()
        self._preview = None
        self._next_operation_number += 1
        return operation

    def undo(self) -> bool:
        if not self._accepted:
            return False
        self._redo.append(self._accepted.pop())
        self._preview = None
        return True

    def redo(self) -> bool:
        if not self._redo:
            return False
        self._accepted.append(self._redo.pop())
        self._preview = None
        return True

    def reset(self) -> None:
        self._accepted.clear()
        self._redo.clear()
        self._excluded_candidate_ids.clear()
        self._preview = None

    def exclude(self, candidate: Candidate | str) -> None:
        self._excluded_candidate_ids.add(self._candidate_id(candidate))

    def restore_excluded(self, candidate: Candidate | str) -> bool:
        candidate_id = self._candidate_id(candidate)
        if candidate_id not in self._excluded_candidate_ids:
            return False
        self._excluded_candidate_ids.remove(candidate_id)
        return True

    def rank_candidates(self, candidates: Iterable[Candidate], element_scope: Iterable[str] | None = None):
        """Rank candidates against the current residual without changing session state."""
        from .ranking import rank_candidates

        return rank_candidates(self, candidates, element_scope=element_scope)

    def to_dict(self) -> dict[str, Any]:
        return {
            "context": self._context_to_dict(self.context),
            "background": {
                "values": self._background_y.tolist(),
                "method": self.background_method,
                "parameters": dict(self._background_parameters),
                "fingerprint": self._array_fingerprint(self._background_y),
            },
            "bounds": self._bounds_to_dict(self.bounds),
            "accepted": [self._stored_to_dict(stored) for stored in self._accepted],
            "redo": [self._stored_to_dict(stored) for stored in self._redo],
            "excluded_candidate_ids": list(self.excluded_candidate_ids),
            "next_operation_number": self._next_operation_number,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> PhaseStrippingSession:
        context_data = payload["context"]
        context = AnalysisContext(
            x=np.array(context_data["x"], dtype=np.dtype(context_data["x_dtype"])),
            intensity=np.array(context_data["intensity"], dtype=np.dtype(context_data["intensity_dtype"])),
            radiation=context_data["radiation"],
            wavelength_angstrom=float(context_data["wavelength_angstrom"]),
            two_theta_range=tuple(context_data["two_theta_range"]),
            tolerance_deg=float(context_data["tolerance_deg"]),
            source_path=context_data["source_path"],
            source_fingerprint=context_data["source_fingerprint"],
            data_fingerprint=context_data["data_fingerprint"],
        )
        background_data = payload.get("background")
        session = cls(
            context,
            cls._bounds_from_dict(payload["bounds"]),
            background_y=None if background_data is None else np.array(background_data["values"], dtype=float),
            background_method="none" if background_data is None else str(background_data.get("method", "none")),
            background_parameters={} if background_data is None else dict(background_data.get("parameters", {})),
        )
        session._accepted = [cls._stored_from_dict(item) for item in payload["accepted"]]
        session._redo = [cls._stored_from_dict(item) for item in payload["redo"]]
        session._excluded_candidate_ids = set(payload["excluded_candidate_ids"])
        session._next_operation_number = int(payload["next_operation_number"])
        return session

    @staticmethod
    def _candidate_id(candidate: Candidate | str) -> str:
        return candidate if isinstance(candidate, str) else candidate.candidate_id

    @staticmethod
    def _array_fingerprint(values: np.ndarray) -> str:
        return hashlib.sha256(np.ascontiguousarray(values).tobytes()).hexdigest()

    @staticmethod
    def _negative_area(values: np.ndarray) -> float:
        return float(np.maximum(-values, 0.0).sum())

    def _negative_area_warnings(self, residual_before: np.ndarray, contribution: np.ndarray) -> tuple[str, ...]:
        before = self._negative_area(residual_before)
        after = self._negative_area(residual_before - contribution)
        if after <= before:
            return ()
        return ("Accepted contribution increases negative residual area.",)

    @staticmethod
    def _fit_boundary_warnings(fit: PhaseFit, bounds: FitBounds) -> tuple[str, ...]:
        warnings: list[str] = []
        values = (
            ("global shift", fit.shift_deg, bounds.shift_deg),
            ("peak width", fit.sigma_deg, bounds.sigma_deg),
            ("scale", fit.scale, bounds.scale),
        )
        for label, value, (lower, upper) in values:
            span = upper - lower
            tolerance = max(1e-9, abs(span) * 1e-5) if np.isfinite(span) else 1e-9
            if np.isfinite(lower) and abs(value - lower) <= tolerance:
                warnings.append(f"Fitted {label} reached its lower bound; inspect the match before accepting.")
            elif np.isfinite(upper) and abs(value - upper) <= tolerance:
                warnings.append(f"Fitted {label} reached its upper bound; inspect the match before accepting.")
        return tuple(warnings)

    @staticmethod
    def _require_bounded(name: str, value: float, bounds: tuple[float, float]) -> None:
        lower, upper = bounds
        if value < lower or value > upper:
            raise ValueError(f"Manual {name} value {value} is outside bounds [{lower}, {upper}].")

    @staticmethod
    def _context_to_dict(context: AnalysisContext) -> dict[str, Any]:
        return {
            "x": context.x.tolist(),
            "intensity": context.intensity.tolist(),
            "x_dtype": context.x.dtype.str,
            "intensity_dtype": context.intensity.dtype.str,
            "radiation": context.radiation,
            "wavelength_angstrom": context.wavelength_angstrom,
            "two_theta_range": list(context.two_theta_range),
            "tolerance_deg": context.tolerance_deg,
            "source_path": context.source_path,
            "source_fingerprint": context.source_fingerprint,
            "data_fingerprint": context.data_fingerprint,
        }

    @staticmethod
    def _bounds_to_dict(bounds: FitBounds) -> dict[str, list[float | str]]:
        return {
            "shift_deg": PhaseStrippingSession._json_numbers(bounds.shift_deg),
            "sigma_deg": PhaseStrippingSession._json_numbers(bounds.sigma_deg),
            "scale": PhaseStrippingSession._json_numbers(bounds.scale),
        }

    @staticmethod
    def _bounds_from_dict(data: dict[str, list[float | str]]) -> FitBounds:
        return FitBounds(
            shift_deg=PhaseStrippingSession._float_pair(data["shift_deg"]),
            sigma_deg=PhaseStrippingSession._float_pair(data["sigma_deg"]),
            scale=PhaseStrippingSession._float_pair(data["scale"]),
        )

    @staticmethod
    def _json_numbers(values: tuple[float, float]) -> list[float | str]:
        return [PhaseStrippingSession._json_number(value) for value in values]

    @staticmethod
    def _json_number(value: float) -> float | str:
        number = float(value)
        if np.isnan(number):
            raise ValueError("Fit bounds cannot contain NaN values.")
        if np.isposinf(number):
            return "+inf"
        if np.isneginf(number):
            return "-inf"
        return number

    @staticmethod
    def _float_pair(values: list[float | str]) -> tuple[float, float]:
        if len(values) != 2:
            raise ValueError("Each serialized fit bound must contain two values.")
        first, second = values
        return (
            PhaseStrippingSession._json_float(first),
            PhaseStrippingSession._json_float(second),
        )

    @staticmethod
    def _json_float(value: float | str) -> float:
        if value == "+inf":
            return float("inf")
        if value == "-inf":
            return float("-inf")
        if isinstance(value, str):
            raise ValueError("Serialized fit bounds must use numeric values or infinity markers.")
        number = float(value)
        if np.isnan(number):
            raise ValueError("Serialized fit bounds cannot contain NaN values.")
        return number

    @staticmethod
    def _stored_to_dict(stored: _StoredOperation) -> dict[str, Any]:
        phase_fit = stored.operation.phase_fit
        return {
            "operation_id": stored.operation.operation_id,
            "residual_fingerprint": stored.operation.residual_fingerprint,
            "phase_fit": {
                "candidate_id": phase_fit.candidate_id,
                "shift_deg": phase_fit.shift_deg,
                "sigma_deg": phase_fit.sigma_deg,
                "scale": phase_fit.scale,
            },
            "contribution": stored.contribution.tolist(),
            "warnings": list(stored.warnings),
            "negative_area_before": stored.negative_area_before,
            "negative_area_after": stored.negative_area_after,
            "candidate": None if stored.candidate_identity is None else stored.candidate_identity.to_dict(),
        }

    @staticmethod
    def _stored_from_dict(data: dict[str, Any]) -> _StoredOperation:
        fit_data = data["phase_fit"]
        phase_fit = PhaseFit(
            candidate_id=fit_data["candidate_id"],
            shift_deg=float(fit_data["shift_deg"]),
            sigma_deg=float(fit_data["sigma_deg"]),
            scale=float(fit_data["scale"]),
        )
        return _StoredOperation(
            operation=PhaseOperation(
                operation_id=data["operation_id"],
                phase_fit=phase_fit,
                residual_fingerprint=data["residual_fingerprint"],
            ),
            contribution=np.array(data["contribution"], dtype=float),
            warnings=tuple(data["warnings"]),
            negative_area_before=float(data["negative_area_before"]),
            negative_area_after=float(data["negative_area_after"]),
            candidate_identity=None
            if data.get("candidate") is None
            else _CandidateIdentity.from_dict(data["candidate"]),
        )
