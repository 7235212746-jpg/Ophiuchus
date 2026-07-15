from __future__ import annotations

import csv
from pathlib import Path

from .config import XRDConfig
from .models import Peak
from .validation import PatternComparator, ReferencePatternImporter, simulate_cif_with_config


def simulate_cif_debug_table(cif_file: str | Path, out_csv: str | Path, config: XRDConfig) -> str:
    cif_path = Path(cif_file)
    out_path = Path(out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    peaks = simulate_cif_with_config(cif_path, config)
    fieldnames = [
        "peak_id",
        "two_theta",
        "d_spacing",
        "relative_intensity",
        "hkl",
        "multiplicity",
        "wavelength_angstrom",
        "radiation_mode",
        "line_model",
        "debye_waller_b",
        "structure_source",
        "cif_path",
    ]
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for peak in peaks:
            writer.writerow(
                {
                    "peak_id": peak.peak_id,
                    "two_theta": f"{peak.two_theta:.6f}",
                    "d_spacing": "" if peak.d_spacing is None else f"{peak.d_spacing:.6f}",
                    "relative_intensity": f"{peak.intensity:.6f}",
                    "hkl": peak.hkl,
                    "multiplicity": peak.multiplicity or "",
                    "wavelength_angstrom": f"{config.wavelength_angstrom:.6f}",
                    "radiation_mode": config.radiation_source,
                    "line_model": config.line_model,
                    "debye_waller_b": f"{config.debye_waller_b:.6f}",
                    "structure_source": "cif",
                    "cif_path": str(cif_path),
                }
            )
    return str(out_path)


def compare_debug_peak_tables(
    ophi_csv: str | Path,
    reference_file: str | Path,
    out_csv: str | Path,
    position_tolerance: float = 0.05,
    phase: str | None = None,
) -> str:
    ophi_peaks = _read_ophi_debug_peaks(ophi_csv)
    reference = ReferencePatternImporter().read(reference_file, phase=phase)
    report = PatternComparator(tolerance_deg=position_tolerance, strong_intensity_threshold=5.0).compare(ophi_peaks, reference.peaks)
    out_path = Path(out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "status",
        "reference_two_theta",
        "ophi_two_theta",
        "delta_2theta",
        "reference_intensity",
        "ophi_intensity",
        "intensity_ratio_ophi_over_reference",
        "reference_hkl",
        "ophi_hkl",
        "note",
    ]
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for item in report.matched:
            ref_i = float(item["reference_intensity"])
            calc_i = float(item["calculated_intensity"])
            ratio = calc_i / ref_i if ref_i else None
            note = ""
            if ratio is not None and (ratio >= 2.0 or ratio <= 0.5) and ref_i >= 5.0:
                note = "strong intensity mismatch"
            writer.writerow(
                {
                    "status": "matched",
                    "reference_two_theta": item["reference_two_theta"],
                    "ophi_two_theta": item["calculated_two_theta"],
                    "delta_2theta": item["delta_2theta"],
                    "reference_intensity": ref_i,
                    "ophi_intensity": calc_i,
                    "intensity_ratio_ophi_over_reference": "" if ratio is None else ratio,
                    "reference_hkl": item.get("reference_hkl", ""),
                    "ophi_hkl": "",
                    "note": note,
                }
            )
        for item in report.missing_reference:
            writer.writerow(
                {
                    "status": "missing_reference",
                    "reference_two_theta": item["two_theta"],
                    "reference_intensity": item["intensity"],
                    "reference_hkl": item.get("hkl", ""),
                    "note": "reference peak not matched by Ophi",
                }
            )
        for item in report.extra_calculated:
            writer.writerow(
                {
                    "status": "extra_ophi",
                    "ophi_two_theta": item["two_theta"],
                    "ophi_intensity": item["intensity"],
                    "ophi_hkl": item.get("hkl", ""),
                    "note": "Ophi peak not matched by reference",
                }
            )
    return str(out_path)


def _read_ophi_debug_peaks(path: str | Path) -> list[Peak]:
    peaks: list[Peak] = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for i, row in enumerate(reader, 1):
            two_theta = row.get("two_theta") or row.get("2theta")
            intensity = row.get("relative_intensity") or row.get("intensity")
            if not two_theta or not intensity:
                continue
            peaks.append(
                Peak(
                    row.get("peak_id") or f"ophi_{i}",
                    float(two_theta),
                    float(intensity),
                    hkl=row.get("hkl") or "",
                    d_spacing=float(row["d_spacing"]) if row.get("d_spacing") else None,
                    multiplicity=int(float(row["multiplicity"])) if row.get("multiplicity") else None,
                )
            )
    return peaks
