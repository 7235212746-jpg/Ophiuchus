from __future__ import annotations

import hashlib
import re
from pathlib import Path

from .importers import normalize_pattern
from .config import XRDConfig
from .models import Candidate, Peak
from .peaks import find_peaks
from .backend import SimulationContext, ValidatedXRDBackend


PEAK_EXTENSIONS = {".int", ".txt", ".xy", ".dat", ".csv"}
GENERIC_PEAK_LIST_NAMES = {
    "moni",
    "relax",
    "simulated_peak_positions",
    "reference_simulated_peak_positions",
    "all_peak_positions",
    "experimental_peak_positions",
}


def formula_elements(formula: str) -> list[str]:
    return re.findall(r"[A-Z][a-z]?", formula)


def formula_from_name(path: Path) -> str:
    stem = path.stem
    stem = re.sub(r"[_\-\s]*(2theta|2θ|theta|raw|asc)$", "", stem, flags=re.I)
    match = re.search(r"([A-Z][A-Za-z0-9\.]*)", stem)
    return match.group(1) if match else stem


def formula_from_path(path: Path) -> str:
    formula = formula_from_name(path)
    if path.stem.lower() in GENERIC_PEAK_LIST_NAMES:
        parent_formula = formula_from_name(path.parent)
        if parent_formula:
            return parent_formula
    return formula


def structure_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()[:16]


def read_peak_list_file(path: str | Path, intensity_threshold: float = 0.0, max_peaks: int = 120) -> list[Peak]:
    file_path = Path(path)
    rows: list[tuple[float, float]] = []
    for line in file_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        parts = line.strip().replace(",", " ").split()
        if len(parts) < 2:
            continue
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except ValueError:
            continue
    normalized = normalize_pattern(rows)
    if _looks_like_continuous_pattern(rows):
        peaks = find_peaks(
            normalized,
            two_theta_min=min(x for x, _ in rows),
            two_theta_max=max(x for x, _ in rows),
            min_height=max(1.0, intensity_threshold),
            min_distance_deg=0.25,
            max_peaks=max_peaks,
            smooth_window=5,
        )
        peaks = [Peak(f"theory_{i}", p.two_theta, p.intensity, p.prominence, p.fwhm) for i, p in enumerate(peaks, 1)]
    else:
        peaks = [
            Peak(f"theory_{i}", point.two_theta, point.intensity)
            for i, point in enumerate(normalized, 1)
            if point.intensity >= intensity_threshold
        ]
    peaks.sort(key=lambda p: p.intensity, reverse=True)
    return sorted(peaks[:max_peaks], key=lambda p: p.two_theta)


def _looks_like_continuous_pattern(rows: list[tuple[float, float]]) -> bool:
    if len(rows) > 200:
        return True
    if len(rows) < 20:
        return False
    diffs = [round(rows[i + 1][0] - rows[i][0], 5) for i in range(min(len(rows) - 1, 50))]
    positive = [d for d in diffs if d > 0]
    if len(positive) < 10:
        return False
    median = sorted(positive)[len(positive) // 2]
    return median <= 0.1 and sum(abs(d - median) <= median * 0.2 for d in positive) >= len(positive) * 0.8


class LocalCandidateProvider:
    def __init__(
        self,
        roots: list[str | Path],
        allowed_elements: set[str] | None = None,
        include_peak_lists: bool = True,
    ) -> None:
        self.roots = [Path(root) for root in roots]
        self.allowed_elements = allowed_elements
        self.include_peak_lists = include_peak_lists

    def iter_candidates(self) -> list[Candidate]:
        candidates: list[Candidate] = []
        for root in self.roots:
            if not root.exists():
                continue
            for path in sorted(root.rglob("*")):
                if not path.is_file():
                    continue
                suffix = path.suffix.lower()
                if suffix == ".cif":
                    candidate = self._candidate_from_cif(path)
                elif self.include_peak_lists and suffix in PEAK_EXTENSIONS:
                    candidate = self._candidate_from_peak_list(path)
                else:
                    continue
                if self._allowed(candidate.elements):
                    candidates.append(candidate)
        return candidates

    def _candidate_from_cif(self, path: Path) -> Candidate:
        formula = formula_from_name(path)
        elements = sorted(set(formula_elements(formula) or _elements_from_cif_text(path)))
        return Candidate(
            candidate_id=structure_hash(path),
            formula_pretty=formula,
            source="local_cif",
            source_path=str(path),
            elements=elements,
            structure_hash=structure_hash(path),
        )

    def _candidate_from_peak_list(self, path: Path) -> Candidate:
        formula = formula_from_path(path)
        elements = sorted(set(formula_elements(formula)))
        candidate = Candidate(
            candidate_id=structure_hash(path),
            formula_pretty=formula,
            source="local_peak_list",
            source_path=str(path),
            elements=elements,
            structure_hash=structure_hash(path),
        )
        try:
            candidate.theory_peaks = read_peak_list_file(path, intensity_threshold=1.0)
        except Exception as exc:
            candidate.parse_status = "failed"
            candidate.parse_error = str(exc)
        return candidate

    def _allowed(self, elements: list[str]) -> bool:
        if not self.allowed_elements:
            return True
        return bool(elements) and set(elements).issubset(self.allowed_elements)


def _elements_from_cif_text(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    values = re.findall(r"\b([A-Z][a-z]?)\d*\b", text)
    return [value for value in values if len(value) <= 2]


def simulate_or_load_peaks(candidate: Candidate, radiation: str = "CuKa", two_theta_range: tuple[float, float] = (10, 90)) -> list[Peak]:
    if candidate.theory_peaks:
        return candidate.theory_peaks
    if candidate.path.suffix.lower() in PEAK_EXTENSIONS:
        candidate.theory_peaks = read_peak_list_file(candidate.path, intensity_threshold=1.0)
        return candidate.theory_peaks
    wavelength = WAVELENGTH_FOR_RADIATION.get(radiation, 1.54056)
    config = XRDConfig(
        radiation_source=radiation,
        wavelength_angstrom=wavelength,
        two_theta_min=float(two_theta_range[0]),
        two_theta_max=float(two_theta_range[1]),
    )
    try:
        candidate.simulated_pattern = ValidatedXRDBackend().simulate_cif(
            candidate.path,
            config,
            SimulationContext(
                structure_id=candidate.candidate_id,
                source=candidate.source,
                formula=candidate.formula_pretty,
            ),
        )
    except Exception as exc:
        candidate.parse_status = "simulation_failed"
        candidate.parse_error = _trusted_simulator_error(str(exc))
        raise RuntimeError(candidate.parse_error) from exc
    candidate.theory_peaks = candidate.simulated_pattern.to_peaks()
    candidate.parse_status = "validated_backend_simulated"
    return candidate.theory_peaks


def _trusted_simulator_error(detail: str = "") -> str:
    suffix = f" Detail: {detail}" if detail else ""
    return (
        "The canonical validated XRD backend could not simulate this CIF. "
        "Check that the CIF is complete and readable by pymatgen; no alternative simulator will be substituted."
        f"{suffix}"
    )


WAVELENGTH_FOR_RADIATION = {
    "CuKa": 1.54056,
    "CuKalpha1": 1.54056,
    "CuKalpha12": 1.54056,
    "CuKα": 1.54056,
    "CoKa": 1.78897,
    "MoKa": 0.71073,
}
