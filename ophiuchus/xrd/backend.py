from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from importlib.metadata import version
from pathlib import Path

from .config import XRDConfig
from .models import Peak


BACKEND_NAME = "Ophi Validated pymatgen XRD"
BACKEND_VERSION = "validated_xrd_v3"
CU_KALPHA2_ANGSTROM = 1.54439
CU_KALPHA2_RATIO = 0.5
SYMMETRY_PROCESSING_VERSION = "pymatgen_structure_from_cif_v1"
INTENSITY_ALGORITHM_VERSION = "pymatgen_xrdcalculator_scaled_false_b0_v1"
PEAK_MERGE_ALGORITHM_VERSION = "same_line_weighted_sum_v1"


@dataclass(frozen=True)
class SimulationContext:
    structure_id: str
    source: str
    formula: str = ""
    space_group_number: int | None = None


@dataclass(frozen=True)
class SimulatedPattern:
    structure_id: str
    source: str
    cif_path: str
    cif_sha256: str
    formula: str
    space_group_number: int | None
    radiation: str
    wavelength_angstrom: float
    two_theta_min_deg: float
    two_theta_max_deg: float
    engine_name: str
    engine_version: str
    dependency_version: str
    settings_fingerprint: str
    pattern_fingerprint: str
    two_theta_deg: tuple[float, ...]
    d_spacing_angstrom: tuple[float | None, ...]
    hkl: tuple[str, ...]
    multiplicity: tuple[int | None, ...]
    line_component: tuple[str, ...]
    raw_intensity: tuple[float, ...]
    normalized_intensity: tuple[float, ...]
    profile_two_theta_deg: tuple[float, ...] = ()
    profile_normalized_intensity: tuple[float, ...] = ()

    def __post_init__(self) -> None:
        lengths = {
            len(self.two_theta_deg),
            len(self.d_spacing_angstrom),
            len(self.hkl),
            len(self.multiplicity),
            len(self.line_component),
            len(self.raw_intensity),
            len(self.normalized_intensity),
        }
        if len(lengths) != 1:
            raise ValueError("All SimulatedPattern peak arrays must have equal length.")
        if tuple(sorted(self.two_theta_deg)) != self.two_theta_deg:
            raise ValueError("SimulatedPattern.two_theta_deg must be sorted.")
        if len(self.profile_two_theta_deg) != len(self.profile_normalized_intensity):
            raise ValueError("SimulatedPattern profile arrays must have equal length.")

    @staticmethod
    def calculate_cif_sha256(cif_path: str | Path) -> str:
        digest = hashlib.sha256()
        with Path(cif_path).open("rb") as handle:
            for chunk in iter(lambda: handle.read(65536), b""):
                digest.update(chunk)
        return digest.hexdigest()

    def to_peaks(self) -> list[Peak]:
        return [
            Peak(
                peak_id=f"canonical_{index}",
                two_theta=two_theta,
                intensity=self.normalized_intensity[index - 1],
                hkl=self.hkl[index - 1],
                d_spacing=self.d_spacing_angstrom[index - 1],
                multiplicity=self.multiplicity[index - 1],
            )
            for index, two_theta in enumerate(self.two_theta_deg, 1)
        ]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)

    @classmethod
    def from_dict(cls, payload: dict[str, object]) -> "SimulatedPattern":
        tuple_fields = {
            "two_theta_deg",
            "d_spacing_angstrom",
            "hkl",
            "multiplicity",
            "line_component",
            "raw_intensity",
            "normalized_intensity",
            "profile_two_theta_deg",
            "profile_normalized_intensity",
        }
        data = dict(payload)
        for field in tuple_fields:
            data[field] = tuple(data.get(field) or [])
        return cls(**data)  # type: ignore[arg-type]


class ValidatedXRDBackend:
    name = BACKEND_NAME
    engine_version = BACKEND_VERSION

    def settings_fingerprint(self, settings: XRDConfig) -> str:
        return _fingerprint(self._settings_payload(settings))

    def simulate_cif(
        self,
        cif_path: str | Path,
        settings: XRDConfig,
        context: SimulationContext,
    ) -> SimulatedPattern:
        from pymatgen.analysis.diffraction.xrd import XRDCalculator
        from pymatgen.core import Structure

        path = Path(cif_path).resolve()
        try:
            structure = Structure.from_file(str(path))
        except Exception as exc:
            raise RuntimeError(
                f"Canonical XRD backend could not parse CIF '{path.name}': {exc}"
            ) from exc
        debye_waller = None
        if settings.debye_waller_b:
            debye_waller = {element.symbol: settings.debye_waller_b for element in structure.composition.elements}
        calculator = XRDCalculator(
            wavelength=settings.wavelength_angstrom,
            debye_waller_factors=debye_waller,
        )
        diffraction = calculator.get_pattern(
            structure,
            scaled=False,
            two_theta_range=settings.two_theta_range(),
        )
        rows: list[dict[str, object]] = []
        for index, (two_theta, raw_intensity) in enumerate(zip(diffraction.x, diffraction.y), 1):
            hkl_items = diffraction.hkls[index - 1] if index - 1 < len(diffraction.hkls) else []
            hkl_text = ";".join(str(item.get("hkl", "")) for item in hkl_items)
            multiplicity = sum(int(item.get("multiplicity", 0) or 0) for item in hkl_items) or None
            d_spacing = diffraction.d_hkls[index - 1] if index - 1 < len(diffraction.d_hkls) else None
            rows.append(
                {
                    "two_theta": float(two_theta),
                    "raw_intensity": float(raw_intensity),
                    "d_spacing": None if d_spacing is None else float(d_spacing),
                    "hkl": hkl_text,
                    "multiplicity": multiplicity,
                    "line_component": "Kalpha1",
                }
            )
        rows = self._apply_line_model(rows, settings)
        rows = self._merge_nearby_rows(rows, settings.peak_merge_tolerance_deg)
        rows.sort(key=lambda item: float(item["two_theta"]))
        max_raw = max((float(item["raw_intensity"]) for item in rows), default=0.0) or 1.0
        normalized = tuple(float(item["raw_intensity"]) / max_raw * 100.0 for item in rows)
        dependency_version = version("pymatgen")
        settings_fingerprint = self.settings_fingerprint(settings)
        two_theta = tuple(float(item["two_theta"]) for item in rows)
        raw = tuple(float(item["raw_intensity"]) for item in rows)
        pattern_fingerprint = _fingerprint(
            {
                "settings": settings_fingerprint,
                "two_theta_deg": two_theta,
                "raw_intensity": raw,
                "normalized_intensity": normalized,
            }
        )
        return SimulatedPattern(
            structure_id=context.structure_id,
            source=context.source,
            cif_path=str(path),
            cif_sha256=SimulatedPattern.calculate_cif_sha256(path),
            formula=context.formula or structure.composition.reduced_formula,
            space_group_number=context.space_group_number,
            radiation=settings.radiation_source,
            wavelength_angstrom=settings.wavelength_angstrom,
            two_theta_min_deg=settings.two_theta_min,
            two_theta_max_deg=settings.two_theta_max,
            engine_name=self.name,
            engine_version=self.engine_version,
            dependency_version=dependency_version,
            settings_fingerprint=settings_fingerprint,
            pattern_fingerprint=pattern_fingerprint,
            two_theta_deg=two_theta,
            d_spacing_angstrom=tuple(item["d_spacing"] for item in rows),
            hkl=tuple(str(item["hkl"]) for item in rows),
            multiplicity=tuple(item["multiplicity"] for item in rows),
            line_component=tuple(str(item["line_component"]) for item in rows),
            raw_intensity=raw,
            normalized_intensity=normalized,
        )

    def _settings_payload(self, settings: XRDConfig) -> dict[str, object]:
        return {
            **settings.to_dict(),
            "backend_name": self.name,
            "backend_version": self.engine_version,
            "pymatgen_version": version("pymatgen"),
            "cu_kalpha2_angstrom": CU_KALPHA2_ANGSTROM,
            "cu_kalpha2_ratio": CU_KALPHA2_RATIO,
            "symmetry_processing_version": SYMMETRY_PROCESSING_VERSION,
            "intensity_algorithm_version": INTENSITY_ALGORITHM_VERSION,
            "peak_merge_algorithm_version": PEAK_MERGE_ALGORITHM_VERSION,
        }

    def _apply_line_model(self, rows: list[dict[str, object]], settings: XRDConfig) -> list[dict[str, object]]:
        model = settings.line_model.lower().replace("-", "_")
        if model in {"kalpha1", "kalpha1_only", "cu_kalpha1"}:
            return list(rows)
        if model not in {"cu_kalpha12", "kalpha12", "cu_ka12", "cu_ka1ka2"}:
            raise ValueError(f"Unsupported XRD line model: {settings.line_model}")
        combined = list(rows)
        for row in rows:
            d_spacing = row["d_spacing"]
            if d_spacing is None:
                continue
            argument = CU_KALPHA2_ANGSTROM / (2.0 * float(d_spacing))
            if not 0.0 < argument <= 1.0:
                continue
            two_theta = math.degrees(2.0 * math.asin(argument))
            if settings.two_theta_min <= two_theta <= settings.two_theta_max:
                combined.append(
                    {
                        **row,
                        "two_theta": two_theta,
                        "raw_intensity": float(row["raw_intensity"]) * CU_KALPHA2_RATIO,
                        "hkl": f"{row['hkl']} Kalpha2".strip(),
                        "line_component": "Kalpha2",
                    }
                )
        return combined

    def _merge_nearby_rows(
        self,
        rows: list[dict[str, object]],
        tolerance_deg: float,
    ) -> list[dict[str, object]]:
        if tolerance_deg <= 0.0:
            return list(rows)
        merged: list[dict[str, object]] = []
        for component in sorted({str(row["line_component"]) for row in rows}):
            component_rows = sorted(
                (row for row in rows if str(row["line_component"]) == component),
                key=lambda row: float(row["two_theta"]),
            )
            cluster: list[dict[str, object]] = []
            for row in component_rows:
                if cluster and float(row["two_theta"]) - float(cluster[0]["two_theta"]) > tolerance_deg:
                    merged.append(self._merge_cluster(cluster))
                    cluster = []
                cluster.append(row)
            if cluster:
                merged.append(self._merge_cluster(cluster))
        return merged

    @staticmethod
    def _merge_cluster(cluster: list[dict[str, object]]) -> dict[str, object]:
        if len(cluster) == 1:
            return dict(cluster[0])
        total = sum(float(row["raw_intensity"]) for row in cluster)
        weights = [float(row["raw_intensity"]) / total if total else 1.0 / len(cluster) for row in cluster]
        hkl_values = list(dict.fromkeys(str(row["hkl"]) for row in cluster if str(row["hkl"])))
        multiplicities = [row["multiplicity"] for row in cluster if row["multiplicity"] is not None]
        d_values = [row["d_spacing"] for row in cluster]
        return {
            "two_theta": sum(float(row["two_theta"]) * weight for row, weight in zip(cluster, weights)),
            "raw_intensity": total,
            "d_spacing": (
                None
                if any(value is None for value in d_values)
                else sum(float(value) * weight for value, weight in zip(d_values, weights))
            ),
            "hkl": ";".join(hkl_values),
            "multiplicity": sum(int(value) for value in multiplicities) if multiplicities else None,
            "line_component": cluster[0]["line_component"],
        }


def _fingerprint(payload: object) -> str:
    text = json.dumps(payload, ensure_ascii=True, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()
