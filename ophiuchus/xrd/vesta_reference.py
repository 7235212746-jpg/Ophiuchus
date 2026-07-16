from __future__ import annotations

import csv
import math
import os
import re
from pathlib import Path

from ..paths import desktop_dir
from .models import Peak
from .peaks import find_peaks
from .importers import normalize_pattern


def find_local_vesta_reference(formula: str) -> dict[str, Path] | None:
    references = find_local_vesta_references(formula)
    return references[0] if references else None


def find_local_vesta_references(formula: str) -> list[dict[str, Path | str]]:
    roots = []
    configured = os.environ.get("OPHI_VESTA_REFERENCE_DIR")
    if configured:
        roots.append(Path(configured))
    else:
        desktop = desktop_dir()
        roots.extend(
            [
                desktop / "03_实验数据与分析" / "XRD与研相" / "XRD",
                desktop / "03_实验数据与分析" / "结构与CIF" / "结构",
                desktop / "XRD",
                desktop / "结构",
                desktop,
            ]
        )
    needle = _norm_formula_text(formula)
    matches: list[Path] = []
    override = _formula_override_reference(needle)
    if override:
        matches.append(override)
    for root in roots:
        if not root.exists():
            continue
        for path in _reference_candidates(root):
            if _path_matches_exact_formula(path, needle):
                matches.append(path)
    if not matches:
        return []
    unique_matches = _dedupe_paths(matches)
    unique_matches.sort(key=lambda path: (-1, 0, 0, 0) if override and path.resolve() == override.resolve() else _reference_sort_key(path))
    return [_reference_record(path, rank_reason="formula-specific override" if override and path.resolve() == override.resolve() else None) for path in unique_matches]


def _reference_record(pattern: Path, rank_reason: str | None = None) -> dict[str, Path | str]:
    peaks = pattern.parent / "reference_simulated_peak_positions.csv"
    if not peaks.exists():
        fallback_peaks = pattern.parent / "simulated_peak_positions.csv"
        peaks = fallback_peaks if fallback_peaks.exists() else peaks
    record: dict[str, Path | str] = {"pattern": pattern, "rank_reason": rank_reason or _reference_rank_reason(pattern)}
    if peaks.exists():
        record["peaks"] = peaks
    return record


def _formula_override_reference(needle: str) -> Path | None:
    if not needle:
        return None
    key = f"OPHI_VESTA_REFERENCE_{needle.upper()}"
    value = os.environ.get(key, "").strip()
    if not value:
        return None
    path = Path(value)
    return path if path.exists() and path.is_file() else None


def _dedupe_paths(paths: list[Path]) -> list[Path]:
    seen: set[str] = set()
    unique: list[Path] = []
    for path in paths:
        key = str(path.resolve()).lower()
        if key in seen:
            continue
        seen.add(key)
        unique.append(path)
    return unique


def read_reference_pattern(path: Path, x_min: float, x_max: float) -> tuple[list[float], list[float]]:
    rows = _read_two_column_rows(path, x_min, x_max)
    if not rows:
        return [x_min, x_max], [0.0, 0.0]
    max_y = max(y for _x, y in rows) or 1.0
    return [x for x, _y in rows], [y / max_y * 100.0 for _x, y in rows]


def read_reference_peak_positions(path: Path | None, x_min: float, x_max: float) -> list[Peak]:
    if not path or not path.exists():
        return []
    peaks: list[Peak] = []
    with path.open(newline="", encoding="utf-8", errors="ignore") as handle:
        reader = csv.DictReader(handle)
        for i, row in enumerate(reader, 1):
            lower = {str(key).strip().lower(): value for key, value in row.items() if key is not None}
            x_text = lower.get("two_theta") or lower.get("2theta") or lower.get("2-theta")
            y_text = lower.get("intensity_norm") or lower.get("intensity") or lower.get("relative_intensity")
            if not x_text or not y_text:
                continue
            try:
                x = float(x_text)
                y = float(y_text)
            except ValueError:
                continue
            if x_min <= x <= x_max:
                peaks.append(Peak(f"vesta_ref_{i}", x, y))
    return peaks


def load_local_vesta_reference_peaks(formula: str, x_min: float, x_max: float) -> tuple[list[Peak], Path | None]:
    reference = find_local_vesta_reference(formula)
    if not reference:
        return [], None
    peaks = read_reference_peak_positions(reference.get("peaks"), x_min, x_max)
    if not peaks:
        rows = _read_two_column_rows(reference["pattern"], x_min, x_max)
        if rows:
            normalized = normalize_pattern(rows)
            if _looks_like_discrete_peak_rows(rows):
                peaks = [Peak(f"vesta_ref_{i}", point.two_theta, point.intensity) for i, point in enumerate(normalized, 1) if point.intensity >= 1.0]
            else:
                peaks = find_peaks(
                    normalized,
                    two_theta_min=x_min,
                    two_theta_max=x_max,
                    min_height=1.0,
                    min_distance_deg=0.25,
                    max_peaks=120,
                    smooth_window=5,
                )
    return peaks, reference["pattern"]


def load_local_vesta_reference_for_cif(cif_file: str | Path, x_min: float, x_max: float) -> tuple[list[Peak], Path | None, str | None]:
    """Load a local VESTA reference for a CIF when its formula can be inferred.

    The formula hint intentionally prefers explicit file/data-block names before
    falling back to atom counts. User-exported VESTA references are usually
    organized by target formula, while reduced formulas reconstructed from
    partial test CIFs can be too lossy.
    """
    cif_path = Path(cif_file)
    for formula in _formula_hints_from_cif(cif_path):
        peaks, source = load_local_vesta_reference_peaks(formula, x_min, x_max)
        if peaks:
            return peaks, source, formula
    return [], None, None


def _formula_hints_from_cif(cif_path: Path) -> list[str]:
    hints: list[str] = []
    for text in (cif_path.stem, _read_data_block_name(cif_path), _read_formula_sum(cif_path), _formula_from_atom_counts(cif_path)):
        if text and _looks_like_formula(text) and text not in hints:
            hints.append(text)
    return hints


def _read_two_column_rows(path: Path, x_min: float, x_max: float) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        parts = line.strip().replace(",", " ").split()
        if len(parts) < 2:
            continue
        try:
            x = float(parts[0])
            y = float(parts[1])
        except ValueError:
            continue
        if x_min <= x <= x_max:
            rows.append((x, y))
    return rows


def _looks_like_discrete_peak_rows(rows: list[tuple[float, float]]) -> bool:
    if len(rows) < 20:
        return True
    ordered = sorted(rows)
    diffs = [ordered[i + 1][0] - ordered[i][0] for i in range(min(len(ordered) - 1, 60)) if ordered[i + 1][0] > ordered[i][0]]
    if len(diffs) < 10:
        return True
    median = sorted(diffs)[len(diffs) // 2]
    return median > 0.12


def _norm_formula_text(text: str) -> str:
    return re.sub(r"[^a-z0-9]", "", text.lower())


def _leading_formula_text(text: str) -> str | None:
    match = re.match(r"([A-Z][a-z]?(?:[0-9.]+)?(?:[A-Z][a-z]?(?:[0-9.]+)?)*)", text.strip())
    return match.group(1) if match else None


def _reference_candidates(root: Path) -> list[Path]:
    patterns = ("MONI.int", "RELAX.int", "*.int", "reference_simulated_peak_positions.csv", "simulated_peak_positions.csv")
    candidates: dict[str, Path] = {}
    recursive = root.name not in {"Desktop", "桌面"}
    for pattern in patterns:
        finder = root.rglob(pattern) if recursive else root.glob(pattern)
        for path in finder:
            if path.is_file():
                candidates[str(path.resolve()).lower()] = path
    return list(candidates.values())


def _path_matches_exact_formula(path: Path, needle: str) -> bool:
    if not needle:
        return False
    file_formula = _leading_formula_text(path.stem)
    if file_formula and _norm_formula_text(file_formula) == needle:
        return True
    if path.stem.lower() not in _GENERIC_REFERENCE_STEMS:
        return False
    folder_formula = _leading_formula_text(path.parent.name)
    return bool(folder_formula and _norm_formula_text(folder_formula) == needle)


def _reference_sort_key(path: Path) -> tuple[int, int, int, float]:
    text = f"{path.parent.name} {path.name}".lower()
    exact_int = 0 if path.suffix.lower() == ".int" else 1
    preferred_name = 0 if any(token in text for token in ("vesta", "moni", "sum", "reference")) else 1
    csv_penalty = 1 if path.name.lower() == "simulated_peak_positions.csv" else 0
    return (exact_int, preferred_name, csv_penalty, -path.stat().st_mtime)


def _reference_rank_reason(path: Path) -> str:
    text = f"{path.parent.name} {path.name}".lower()
    if "vesta" in text:
        return "exact formula VESTA file"
    if path.stem.lower() in {"moni", "relax"}:
        return "generic VESTA pattern in formula folder"
    if "reference_simulated_peak_positions" in text:
        return "reference simulated peak table"
    if "simulated_peak_positions" in text:
        return "simulated peak table in formula folder"
    return "exact formula reference file"


_GENERIC_REFERENCE_STEMS = {
    "moni",
    "relax",
    "reference_simulated_peak_positions",
    "simulated_peak_positions",
}


def _read_data_block_name(cif_path: Path) -> str | None:
    for line in cif_path.read_text(encoding="utf-8", errors="ignore").splitlines()[:20]:
        clean = line.strip()
        if clean.lower().startswith("data_"):
            return clean[5:].strip()
    return None


def _read_formula_sum(cif_path: Path) -> str | None:
    text = cif_path.read_text(encoding="utf-8", errors="ignore")
    match = re.search(r"^_chemical_formula_sum\s+(.+)$", text, re.M)
    if not match:
        return None
    raw = match.group(1).strip().strip("'\"")
    tokens = re.findall(r"([A-Z][a-z]?)\s*([0-9.]+)?", raw)
    if not tokens:
        return None
    return "".join(f"{el}{_clean_count(count)}" for el, count in tokens)


def _formula_from_atom_counts(cif_path: Path) -> str | None:
    try:
        from .simulate import parse_simple_cif

        atoms = parse_simple_cif(cif_path).atoms
    except Exception:
        atoms = []
    counts: dict[str, int] = {}
    order: list[str] = []
    for atom in atoms:
        symbol = atom.element
        if symbol not in counts:
            order.append(symbol)
        counts[symbol] = counts.get(symbol, 0) + 1
    if not counts:
        return None
    divisor = 0
    for count in counts.values():
        divisor = count if divisor == 0 else math.gcd(divisor, count)
    divisor = max(divisor, 1)
    return "".join(f"{el}{_formula_count_text(count // divisor)}" for el, count in ((el, counts[el]) for el in order))


def _clean_count(count: str) -> str:
    if not count:
        return ""
    try:
        value = float(count)
    except ValueError:
        return count
    if abs(value - round(value)) <= 1e-6:
        integer = int(round(value))
        return "" if integer == 1 else str(integer)
    return count


def _looks_like_formula(text: str) -> bool:
    return bool(re.search(r"[A-Z][a-z]?", text)) and bool(re.fullmatch(r"[A-Za-z0-9_.+-]+", text))


def _formula_count_text(count: int) -> str:
    return "" if count == 1 else str(count)
