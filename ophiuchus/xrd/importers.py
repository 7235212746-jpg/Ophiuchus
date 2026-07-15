from __future__ import annotations

import csv
import re
from pathlib import Path

from .models import PatternPoint, XrdPattern


COMMENT_PREFIXES = ("#", "//", ";")
XRD_EXTENSIONS = {".asc", ".ras", ".raw", ".xy", ".txt", ".csv", ".dat"}


def normalize_pattern(rows: list[tuple[float, float]]) -> list[PatternPoint]:
    if not rows:
        return []
    max_intensity = max(float(y) for _, y in rows)
    if max_intensity <= 0:
        return [PatternPoint(float(x), 0.0) for x, _ in rows]
    return [PatternPoint(float(x), float(y) / max_intensity * 100.0) for x, y in rows]


def load_xrd_file(path: str | Path, normalize: bool = True) -> XrdPattern:
    file_path = Path(path)
    text = _read_text(file_path)
    if _looks_like_rigaku_asc(text):
        pattern = _read_rigaku_asc(file_path, text)
    else:
        pattern = _read_two_column(file_path, text)
    points = sorted(pattern.points, key=lambda p: p.two_theta)
    if normalize:
        points = normalize_pattern([(p.two_theta, p.intensity) for p in points])
    pattern.points = points
    pattern.metadata.setdefault("filename", file_path.name)
    pattern.metadata.setdefault("path", str(file_path))
    return pattern


def infer_xrd_range(path: str | Path, padding_deg: float = 0.5) -> dict[str, object]:
    root = Path(path)
    files = [root] if root.is_file() else sorted(p for p in root.rglob("*") if p.is_file() and p.suffix.lower() in XRD_EXTENSIONS)
    ranges: list[tuple[float, float, str]] = []
    warnings: list[str] = []
    for file_path in files:
        try:
            pattern = load_xrd_file(file_path, normalize=False)
            xs = [point.two_theta for point in pattern.points]
            if xs:
                ranges.append((min(xs), max(xs), str(file_path)))
        except Exception as exc:
            warnings.append(f"{file_path}: {exc}")
    if not ranges:
        raise ValueError(f"no readable XRD files found for range detection: {root}")
    raw_min = min(item[0] for item in ranges)
    raw_max = max(item[1] for item in ranges)
    return {
        "two_theta_min": max(0.0, raw_min - padding_deg),
        "two_theta_max": raw_max + padding_deg,
        "raw_two_theta_min": raw_min,
        "raw_two_theta_max": raw_max,
        "files_read": [item[2] for item in ranges],
        "warnings": warnings,
    }


def _read_text(path: Path) -> str:
    for encoding in ("utf-8", "cp932", "gbk", "latin-1"):
        try:
            return path.read_text(encoding=encoding, errors="strict")
        except UnicodeError:
            continue
    return path.read_text(encoding="utf-8", errors="ignore")


def _looks_like_rigaku_asc(text: str) -> bool:
    return "*START" in text and "*STEP" in text and "*COUNT" in text


def _number_after_equal(line: str) -> float | None:
    if "=" not in line:
        return None
    match = re.search(r"[-+]?\d+(?:\.\d+)?", line.split("=", 1)[1])
    return float(match.group(0)) if match else None


def _metadata_value(line: str) -> str:
    return line.split("=", 1)[1].strip() if "=" in line else ""


def _read_rigaku_asc(path: Path, text: str) -> XrdPattern:
    start = step = stop = None
    count: int | None = None
    intensities: list[float] = []
    metadata: dict[str, str | float | int] = {"filename": path.name, "format": "rigaku_asc"}

    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("*SAMPLE"):
            metadata["sample"] = _metadata_value(stripped)
        elif stripped.startswith("*DATE"):
            metadata["date"] = _metadata_value(stripped)
        elif stripped.startswith("*START"):
            start = _number_after_equal(stripped)
            metadata["start"] = start
        elif stripped.startswith("*STOP"):
            stop = _number_after_equal(stripped)
            metadata["stop"] = stop
        elif stripped.startswith("*STEP"):
            step = _number_after_equal(stripped)
            metadata["step"] = step
        elif stripped.startswith("*COUNT"):
            parsed_count = _number_after_equal(stripped)
            count = int(parsed_count) if parsed_count is not None else None
            metadata["count"] = count or 0
        elif count is not None and stripped and not stripped.startswith("*"):
            intensities.extend(float(x) for x in re.findall(r"[-+]?\d+(?:\.\d+)?", stripped))
            if len(intensities) >= count:
                break

    if start is None or step is None or count is None:
        raise ValueError(f"{path} is missing START, STEP, or COUNT metadata")
    values = intensities[:count]
    points = [PatternPoint(start + i * step, y) for i, y in enumerate(values)]
    if stop is not None and points:
        metadata["computed_stop"] = points[-1].two_theta
    return XrdPattern(points=points, metadata=metadata)


def _read_two_column(path: Path, text: str) -> XrdPattern:
    rows: list[tuple[float, float]] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith(COMMENT_PREFIXES):
            continue
        normalized = re.sub(r"[,;\t]+", " ", line)
        parts = next(csv.reader([normalized], delimiter=" ", skipinitialspace=True))
        numeric: list[float] = []
        for part in parts:
            try:
                numeric.append(float(part))
            except ValueError:
                numeric = []
                break
        if len(numeric) >= 2:
            rows.append((numeric[0], numeric[1]))
    if not rows:
        raise ValueError(f"{path} does not contain readable two-column XRD data")
    return XrdPattern(
        points=[PatternPoint(x, y) for x, y in rows],
        metadata={"filename": path.name, "format": "two_column"},
    )
