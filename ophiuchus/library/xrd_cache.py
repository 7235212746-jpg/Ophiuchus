from __future__ import annotations

import hashlib
from pathlib import Path

from ophiuchus.xrd.backend import BACKEND_VERSION, SimulatedPattern, SimulationContext, ValidatedXRDBackend
from ophiuchus.xrd.config import XRDConfig
from ophiuchus.xrd.models import Candidate, Peak

from .database import StructureLibrary
from .models import LibraryPeak, StructureEntry


XRD_ENGINE_VERSION = f"library_{BACKEND_VERSION}_cu_kalpha12_raw_pymatgen_b0"


def simulation_settings_hash(radiation: str, two_theta_range: tuple[float, float]) -> str:
    return simulation_settings_fingerprint(radiation, two_theta_range)[:16]


def simulation_settings_fingerprint(radiation: str, two_theta_range: tuple[float, float]) -> str:
    config = _config(radiation, two_theta_range)
    return ValidatedXRDBackend().settings_fingerprint(config)


def build_library_xrd_cache(
    library: StructureLibrary,
    structure_ids: list[str] | None = None,
    radiation: str = "CuKa",
    two_theta_range: tuple[float, float] = (10.0, 90.0),
    force: bool = False,
    scientific_safe_mode: bool = True,
) -> dict[str, object]:
    entries = [library.get_structure(item) for item in structure_ids] if structure_ids else library.list_structures(enabled_only=False)
    summary: dict[str, object] = {
        "checked": 0,
        "cached": 0,
        "simulated": 0,
        "freshly_simulated": 0,
        "validated_cache_hits": 0,
        "failed": 0,
        "scientific_safe_mode": scientific_safe_mode,
        "backend_name": ValidatedXRDBackend.name,
        "backend_version": ValidatedXRDBackend.engine_version,
        "radiation": radiation,
        "wavelength_angstrom": _config(radiation, two_theta_range).wavelength_angstrom,
        "two_theta_range": list(two_theta_range),
        "structure_states": {},
        "warnings": [],
    }
    config = _config(radiation, two_theta_range)
    backend = ValidatedXRDBackend()
    for entry in entries:
        summary["checked"] += 1
        if not force and not scientific_safe_mode and load_validated_pattern(
            library,
            entry,
            radiation=radiation,
            two_theta_range=two_theta_range,
        ) is not None:
            summary["cached"] += 1
            summary["validated_cache_hits"] += 1
            summary["structure_states"][entry.internal_id] = "validated_cache"
            continue
        try:
            candidate = _entry_to_candidate(library, entry)
            pattern = backend.simulate_cif(
                candidate.path,
                config,
                SimulationContext(
                    structure_id=entry.internal_id,
                    source=entry.source,
                    formula=entry.reduced_formula or entry.formula,
                    space_group_number=entry.space_group_number,
                ),
            )
            _store_validated_pattern(library, pattern)
            _store_peaks(
                library,
                entry,
                pattern.to_peaks(),
                radiation,
                simulation_settings_hash(radiation, two_theta_range),
                two_theta_range,
            )
            summary["simulated"] += 1
            summary["freshly_simulated"] += 1
            summary["structure_states"][entry.internal_id] = "freshly_simulated"
        except Exception as exc:
            summary["failed"] += 1
            summary["structure_states"][entry.internal_id] = "failed"
            summary["warnings"].append(f"{entry.formula} {entry.internal_id}: {exc}")
    return summary


def _store_validated_pattern(library: StructureLibrary, pattern: SimulatedPattern) -> None:
    library.store_validated_pattern(_cache_key(pattern.structure_id, pattern.cif_sha256, pattern.settings_fingerprint), pattern.to_dict())


def load_validated_pattern(
    library: StructureLibrary,
    entry: StructureEntry,
    radiation: str = "CuKalpha12",
    two_theta_range: tuple[float, float] = (10.0, 90.0),
) -> SimulatedPattern | None:
    path = Path(library.path.parent) / entry.cached_file_path
    if not path.exists():
        return None
    cif_sha256 = SimulatedPattern.calculate_cif_sha256(path)
    settings_fingerprint = ValidatedXRDBackend().settings_fingerprint(_config(radiation, two_theta_range))
    payload = library.load_validated_pattern_payload(_cache_key(entry.internal_id, cif_sha256, settings_fingerprint))
    if payload is None:
        return None
    pattern = SimulatedPattern.from_dict(payload)
    if pattern.cif_sha256 != cif_sha256 or pattern.settings_fingerprint != settings_fingerprint:
        return None
    return pattern


def _store_peaks(
    library: StructureLibrary,
    entry: StructureEntry,
    peaks: list[Peak],
    radiation: str,
    settings_hash: str,
    two_theta_range: tuple[float, float],
) -> None:
    library_peaks = [
        LibraryPeak(
            structure_internal_id=entry.internal_id,
            peak_id=peak.peak_id or f"theory_{i}",
            two_theta=peak.two_theta,
            relative_intensity=peak.intensity,
            radiation=radiation,
            settings_hash=settings_hash,
            hkl=peak.hkl,
            d_spacing=peak.d_spacing,
        )
        for i, peak in enumerate(peaks, 1)
    ]
    library.store_xrd_peaks(entry.internal_id, library_peaks, radiation, settings_hash, two_theta_range[0], two_theta_range[1])


def compatible_structure_ids(library: StructureLibrary, elements: list[str], enabled_only: bool = False) -> list[str]:
    allowed = set(elements)
    ids: list[str] = []
    for entry in library.list_structures(enabled_only=enabled_only):
        if allowed and not set(entry.elements).issubset(allowed):
            continue
        ids.append(entry.internal_id)
    return ids


def scoped_structure_ids(
    library: StructureLibrary,
    elements: list[str],
    enabled_only: bool = False,
    scope: str = "exact",
) -> list[str]:
    allowed = set(elements)
    ids: list[str] = []
    for entry in library.list_structures(enabled_only=enabled_only):
        entry_elements = set(entry.elements)
        if not allowed:
            ids.append(entry.internal_id)
        elif scope == "subsystems":
            if entry_elements.issubset(allowed):
                ids.append(entry.internal_id)
        else:
            if entry_elements == allowed:
                ids.append(entry.internal_id)
    return ids


def library_entries_to_candidates(
    library: StructureLibrary,
    elements: list[str],
    radiation: str = "CuKa",
    two_theta_range: tuple[float, float] = (10.0, 90.0),
    candidate_ids: set[str] | None = None,
) -> list[Candidate]:
    allowed = set(elements)
    candidates: list[Candidate] = []
    for entry in library.list_structures(enabled_only=True):
        if candidate_ids is not None and entry.internal_id not in candidate_ids:
            continue
        if allowed and not set(entry.elements).issubset(allowed):
            continue
        pattern = load_validated_pattern(library, entry, radiation=radiation, two_theta_range=two_theta_range)
        if pattern is None:
            continue
        candidate = _entry_to_candidate(library, entry)
        candidate.simulated_pattern = pattern
        candidate.theory_peaks = pattern.to_peaks()
        candidate.parse_status = "validated_backend_cache"
        candidates.append(candidate)
    return candidates


def _entry_to_candidate(library: StructureLibrary, entry: StructureEntry) -> Candidate:
    source_path = str(Path(library.path.parent) / entry.cached_file_path)
    return Candidate(
        candidate_id=entry.internal_id,
        formula_pretty=entry.reduced_formula or entry.formula,
        source=f"library:{entry.source}",
        source_path=source_path,
        elements=entry.elements,
        structure_hash=entry.structure_hash,
    )


def _config(radiation: str, two_theta_range: tuple[float, float]) -> XRDConfig:
    wavelength = {
        "CuKa": 1.54056,
        "CuKalpha1": 1.54056,
        "CuKalpha12": 1.54056,
        "CuKα": 1.54056,
        "CoKa": 1.78897,
        "MoKa": 0.71073,
    }.get(radiation, 1.54056)
    line_model = "cu_kalpha12" if radiation in {"CuKa", "CuKalpha12", "CuKα"} else "kalpha1"
    return XRDConfig(
        radiation_source=radiation,
        wavelength_angstrom=wavelength,
        two_theta_min=float(two_theta_range[0]),
        two_theta_max=float(two_theta_range[1]),
        line_model=line_model,
    )


def _cache_key(structure_id: str, cif_sha256: str, settings_fingerprint: str) -> str:
    return hashlib.sha256(f"{structure_id}|{cif_sha256}|{settings_fingerprint}".encode("utf-8")).hexdigest()
