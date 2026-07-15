from __future__ import annotations

import csv
import json
import os
import re
import traceback
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from statistics import median
from typing import Callable

from .config import XRDConfig
from .models import Peak
from .importers import normalize_pattern
from .peaks import find_peaks


@dataclass(frozen=True)
class ReferencePeak(Peak):
    hkl: str = ""


@dataclass(frozen=True)
class ReferencePattern:
    label: str
    source_path: str
    peaks: list[ReferencePeak]


@dataclass
class PatternComparisonReport:
    matched: list[dict[str, object]]
    missing_reference: list[dict[str, object]]
    extra_calculated: list[dict[str, object]]
    summary: dict[str, object]

    def to_dict(self) -> dict[str, object]:
        return {
            "matched": self.matched,
            "missing_reference": self.missing_reference,
            "extra_calculated": self.extra_calculated,
            "summary": self.summary,
        }


class ReferencePatternImporter:
    def read(self, path: str | Path, label: str | None = None, phase: str | None = None) -> ReferencePattern:
        file_path = Path(path)
        text = file_path.read_text(encoding="utf-8", errors="ignore")
        peaks = self._read_csv_like(file_path, text, phase=phase)
        if not peaks:
            peaks = self._read_vesta_like(text)
        if _looks_continuous_reference(peaks):
            normalized = normalize_pattern([(peak.two_theta, peak.intensity) for peak in peaks])
            peaks = [
                ReferencePeak(peak.peak_id, peak.two_theta, peak.intensity)
                for peak in find_peaks(
                    normalized,
                    two_theta_min=min(point.two_theta for point in normalized),
                    two_theta_max=max(point.two_theta for point in normalized),
                    min_height=1.0,
                    min_distance_deg=0.08,
                    max_peaks=160,
                    smooth_window=5,
                )
            ]
        if not peaks:
            raise ValueError(f"reference pattern contains no readable peaks: {file_path}")
        peaks = sorted(peaks, key=lambda peak: peak.two_theta)
        return ReferencePattern(label=label or file_path.stem, source_path=str(file_path), peaks=peaks)

    def _read_csv_like(self, file_path: Path, text: str, phase: str | None = None) -> list[ReferencePeak]:
        first_line = next((line for line in text.splitlines() if line.strip()), "")
        if "," not in first_line and "\t" not in first_line:
            return []
        delimiter = "," if "," in first_line else "\t"
        rows = csv.DictReader(text.splitlines(), delimiter=delimiter)
        peaks: list[ReferencePeak] = []
        for i, row in enumerate(rows, 1):
            lower = {str(key).strip().lower(): value for key, value in row.items() if key is not None}
            if phase and str(lower.get("phase") or "").strip() != phase:
                continue
            two_theta = lower.get("two_theta") or lower.get("2theta") or lower.get("2-theta") or lower.get("2-theta ")
            intensity = (
                lower.get("intensity")
                or lower.get("intensity_norm")
                or lower.get("i")
                or lower.get("i(f)")
                or lower.get("relative_intensity")
            )
            if two_theta is None or intensity is None:
                continue
            try:
                peaks.append(ReferencePeak(f"ref_{i}", float(two_theta), float(intensity), hkl=str(lower.get("hkl") or "")))
            except ValueError:
                continue
        return peaks

    def _read_vesta_like(self, text: str) -> list[ReferencePeak]:
        peaks: list[ReferencePeak] = []
        in_table = False
        intensity_column = 1
        for line in text.splitlines():
            clean = line.strip()
            if not clean:
                continue
            if "2-Theta" in clean or "2theta" in clean.lower() or clean.lower().startswith("two_theta"):
                in_table = True
                if "I(f)" in clean or "d(?)" in clean:
                    intensity_column = 2
                continue
            parts = clean.split()
            if len(parts) < 2:
                continue
            if not in_table and not _looks_numeric(parts[:2]):
                continue
            try:
                two_theta = float(parts[0])
                intensity = float(parts[intensity_column] if len(parts) > intensity_column else parts[1])
            except ValueError:
                continue
            hkl_match = re.search(r"\(([^)]+)\)", clean)
            peaks.append(ReferencePeak(f"ref_{len(peaks)+1}", two_theta, intensity, hkl=hkl_match.group(1).strip() if hkl_match else ""))
        return peaks


class PatternComparator:
    def __init__(self, tolerance_deg: float = 0.10, strong_intensity_threshold: float = 5.0) -> None:
        self.tolerance_deg = tolerance_deg
        self.strong_intensity_threshold = strong_intensity_threshold

    def compare(self, calculated: list[Peak], reference: list[Peak]) -> PatternComparisonReport:
        calculated_norm = _normalize_peaks(calculated)
        reference_norm = _normalize_peaks(reference)
        unused = set(range(len(calculated_norm)))
        matched: list[dict[str, object]] = []
        missing: list[dict[str, object]] = []
        for ref in reference_norm:
            best_idx = None
            best_delta = None
            for idx in list(unused):
                calc = calculated_norm[idx]
                delta = calc.two_theta - ref.two_theta
                if abs(delta) <= self.tolerance_deg and (best_delta is None or abs(delta) < abs(best_delta)):
                    best_idx = idx
                    best_delta = delta
            if best_idx is None:
                missing.append(_peak_dict(ref))
                continue
            calc = calculated_norm[best_idx]
            unused.remove(best_idx)
            ratio = calc.intensity / ref.intensity if ref.intensity else None
            matched.append(
                {
                    "reference_two_theta": ref.two_theta,
                    "calculated_two_theta": calc.two_theta,
                    "delta_2theta": best_delta,
                    "reference_intensity": ref.intensity,
                    "calculated_intensity": calc.intensity,
                    "intensity_ratio_calculated_over_reference": ratio,
                    "intensity_ratio_error": _ratio_error(ratio),
                    "reference_hkl": getattr(ref, "hkl", ""),
                }
            )
        extra = [_peak_dict(calculated_norm[idx]) for idx in sorted(unused, key=lambda item: calculated_norm[item].two_theta)]
        strong_ref = [peak for peak in reference_norm if peak.intensity >= self.strong_intensity_threshold]
        matched_ref_positions = {round(float(item["reference_two_theta"]), 5) for item in matched}
        missing_strong = [peak for peak in strong_ref if round(peak.two_theta, 5) not in matched_ref_positions]
        extra_strong = [item for item in extra if float(item["intensity"]) >= self.strong_intensity_threshold]
        deltas = [abs(float(item["delta_2theta"])) for item in matched]
        strong_ratio_errors = [
            float(item["intensity_ratio_error"])
            for item in matched
            if item.get("intensity_ratio_error") is not None and float(item["reference_intensity"]) >= self.strong_intensity_threshold
        ]
        summary = {
            "matched_count": len(matched),
            "reference_count": len(reference_norm),
            "calculated_count": len(calculated_norm),
            "missing_count": len(missing),
            "extra_count": len(extra),
            "missing_strong_count": len(missing_strong),
            "extra_strong_count": len(extra_strong),
            "median_abs_delta_2theta": median(deltas) if deltas else None,
            "max_abs_delta_2theta": max(deltas) if deltas else None,
            "median_strong_intensity_ratio_error": median(strong_ratio_errors) if strong_ratio_errors else None,
            "max_strong_intensity_ratio_error": max(strong_ratio_errors) if strong_ratio_errors else None,
            "strong_peak_coverage": (len(strong_ref) - len(missing_strong)) / len(strong_ref) if strong_ref else None,
            "position_pass": bool(deltas) and (median(deltas) <= 0.03) and (max(deltas) <= 0.10) and not missing_strong,
            "intensity_pass": bool(strong_ratio_errors) and max(strong_ratio_errors) <= 2.0,
        }
        return PatternComparisonReport(matched, missing, extra, summary)


class ValidationRunner:
    def __init__(self, simulator: Callable[[Path, XRDConfig], list[Peak]]) -> None:
        self.simulator = simulator

    def run(
        self,
        cif_file: str | Path,
        reference_file: str | Path,
        out_dir: str | Path,
        config: XRDConfig | None = None,
        label: str | None = None,
        phase: str | None = None,
    ) -> dict[str, str]:
        config = config or XRDConfig.validation_default()
        out_path = Path(out_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        cif_path = Path(cif_file)
        reference = ReferencePatternImporter().read(reference_file, label=label, phase=phase)
        try:
            calculated = self.simulator(cif_path, config)
            comparison = PatternComparator(tolerance_deg=max(0.10, config.peak_merge_tolerance_deg)).compare(calculated, reference.peaks)
            error = None
        except Exception as exc:
            calculated = []
            comparison = PatternComparisonReport([], [_peak_dict(p) for p in reference.peaks], [], {"error": str(exc), "position_pass": False})
            error = {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc(limit=5)}
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        stem = f"{label or cif_path.stem}_{stamp}"
        json_path = out_path / f"{stem}_validation.json"
        md_path = out_path / f"{stem}_validation.md"
        payload = {
            "label": label or cif_path.stem,
            "cif_file": str(cif_path),
            "reference_file": str(reference_file),
            "config": config.to_dict(),
            "config_id": config.config_id(),
            "calculated_peaks": [_peak_dict(peak) for peak in calculated],
            **comparison.to_dict(),
            "error": error,
        }
        json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
        md_path.write_text(_markdown_report(payload), encoding="utf-8")
        return {"json_report": str(json_path), "markdown_report": str(md_path)}


def simulate_cif_with_config(cif_file: str | Path, config: XRDConfig) -> list[Peak]:
    from .backend import SimulationContext, ValidatedXRDBackend

    cif_path = Path(cif_file)
    pattern = ValidatedXRDBackend().simulate_cif(
        cif_path,
        config,
        SimulationContext(structure_id=cif_path.stem, source="single_cif"),
    )
    return pattern.to_peaks()


def _normalize_peaks(peaks: list[Peak]) -> list[Peak]:
    max_i = max((peak.intensity for peak in peaks), default=0.0) or 1.0
    return [
        Peak(
            peak.peak_id,
            peak.two_theta,
            peak.intensity / max_i * 100.0,
            peak.prominence,
            peak.fwhm,
            hkl=getattr(peak, "hkl", ""),
            d_spacing=getattr(peak, "d_spacing", None),
            multiplicity=getattr(peak, "multiplicity", None),
        )
        for peak in peaks
    ]


def _peak_dict(peak: Peak) -> dict[str, object]:
    return {
        "peak_id": peak.peak_id,
        "two_theta": peak.two_theta,
        "intensity": peak.intensity,
        "hkl": getattr(peak, "hkl", ""),
        "d_spacing": getattr(peak, "d_spacing", None),
        "multiplicity": getattr(peak, "multiplicity", None),
    }


def _ratio_error(ratio: float | None) -> float | None:
    if ratio is None or ratio <= 0:
        return None
    return max(ratio, 1.0 / ratio)


def _looks_numeric(parts: list[str]) -> bool:
    try:
        for part in parts:
            float(part)
        return True
    except ValueError:
        return False


def _looks_continuous_reference(peaks: list[ReferencePeak]) -> bool:
    if len(peaks) < 200:
        return False
    ordered = sorted(peaks, key=lambda peak: peak.two_theta)
    diffs = [round(ordered[i + 1].two_theta - ordered[i].two_theta, 5) for i in range(min(len(ordered) - 1, 80))]
    positive = [diff for diff in diffs if diff > 0]
    if len(positive) < 20:
        return False
    median_step = sorted(positive)[len(positive) // 2]
    return median_step <= 0.1 and sum(abs(diff - median_step) <= max(1e-6, median_step * 0.2) for diff in positive) >= len(positive) * 0.8


def _markdown_report(payload: dict[str, object]) -> str:
    summary = payload.get("summary") or {}
    lines = [
        f"# XRD Validation Report: {payload.get('label')}",
        "",
        f"- CIF: `{payload.get('cif_file')}`",
        f"- Reference: `{payload.get('reference_file')}`",
        f"- Config: `{payload.get('config_id')}`",
        "",
        "## Summary",
    ]
    for key, value in dict(summary).items():
        lines.append(f"- {key}: {value}")
    lines.extend(["", "## Scientific Note", "Peak-position agreement is the primary validation target. Intensity agreement is diagnostic and model-dependent."])
    return "\n".join(lines) + "\n"
