from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Iterable
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.core import Structure

from ophiuchus.xrd.models import Peak
from ophiuchus.xrd.validation import PatternComparator, ReferencePatternImporter


CU_KALPHA1 = 1.54056
CU_KALPHA2 = 1.54439
CU_KALPHA2_RATIO = 0.5


def run_xrd_reproduction(
    cif_file: str | Path,
    reference_file: str | Path,
    out_dir: str | Path,
    two_theta_min: float = 10.0,
    two_theta_max: float = 90.0,
    wavelength: float = 1.54056,
    tolerance: float = 0.10,
    line_model: str = "kalpha1",
    profile_sigma: float = 0.055,
    experimental_file: str | Path | None = None,
) -> dict[str, str]:
    cif_path = Path(cif_file)
    reference_path = Path(reference_file)
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    python_peaks = simulate_with_pymatgen(
        cif_path,
        two_theta_range=(two_theta_min, two_theta_max),
        wavelength=wavelength,
        line_model=line_model,
    )
    reference = ReferencePatternImporter().read(reference_path, label=reference_path.stem)
    reference_peaks = [
        peak for peak in reference.peaks if two_theta_min <= peak.two_theta <= two_theta_max
    ]
    comparison = PatternComparator(tolerance_deg=tolerance, strong_intensity_threshold=5.0).compare(
        python_peaks,
        reference_peaks,
    )

    python_csv = out_path / "python_simulated_peaks.csv"
    reference_csv = out_path / "vesta_reference_peaks.csv"
    comparison_csv = out_path / "peak_comparison.csv"
    summary_json = out_path / "summary.json"
    plot_png = out_path / "comparison_plot.png"
    profile_csv = out_path / "profile_comparison.csv"
    profile_png = out_path / "profile_overlay.png"

    _write_peak_csv(python_csv, python_peaks, source="python_pymatgen")
    _write_peak_csv(reference_csv, reference_peaks, source="vesta_reference")
    _write_comparison_csv(comparison_csv, comparison)
    profile_summary = _write_profile_comparison(
        profile_csv,
        profile_png,
        python_peaks,
        reference_path,
        two_theta_min,
        two_theta_max,
        profile_sigma,
    )
    summary = {
        "cif_file": str(cif_path),
        "reference_file": str(reference_path),
        "experimental_file": str(experimental_file) if experimental_file else "",
        "two_theta_min": two_theta_min,
        "two_theta_max": two_theta_max,
        "wavelength": wavelength,
        "line_model": line_model,
        "profile_sigma": profile_sigma,
        "tolerance": tolerance,
        "summary": comparison.summary,
        "profile_summary": profile_summary,
        "outputs": {
            "python_peaks_csv": str(python_csv),
            "reference_peaks_csv": str(reference_csv),
            "comparison_csv": str(comparison_csv),
            "plot_png": str(plot_png),
            "profile_comparison_csv": str(profile_csv),
            "profile_plot_png": str(profile_png),
        },
    }
    summary_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    _plot_comparison(plot_png, python_peaks, reference_peaks, two_theta_min, two_theta_max)

    return {
        "python_peaks_csv": str(python_csv),
        "reference_peaks_csv": str(reference_csv),
        "comparison_csv": str(comparison_csv),
        "summary_json": str(summary_json),
        "plot_png": str(plot_png),
        "profile_comparison_csv": str(profile_csv),
        "profile_plot_png": str(profile_png),
    }


def simulate_with_pymatgen(
    cif_file: str | Path,
    two_theta_range: tuple[float, float],
    wavelength: float,
    line_model: str = "kalpha1",
) -> list[Peak]:
    structure = Structure.from_file(str(cif_file))
    calculator = XRDCalculator(wavelength=wavelength)
    pattern = calculator.get_pattern(structure, two_theta_range=two_theta_range)
    peaks: list[Peak] = []
    for idx, (two_theta, intensity) in enumerate(zip(pattern.x, pattern.y), 1):
        hkls = pattern.hkls[idx - 1] if idx - 1 < len(pattern.hkls) else []
        hkl_text = ";".join(str(item.get("hkl", "")) for item in hkls)
        multiplicity = sum(int(item.get("multiplicity", 0)) for item in hkls) or None
        d_spacing = pattern.d_hkls[idx - 1] if idx - 1 < len(pattern.d_hkls) else None
        peaks.append(
            Peak(
                peak_id=f"pymatgen_{idx}",
                two_theta=float(two_theta),
                intensity=float(intensity),
                hkl=hkl_text,
                d_spacing=float(d_spacing) if d_spacing is not None else None,
                multiplicity=multiplicity,
            )
        )
    if line_model == "cu_kalpha12":
        peaks = _add_cu_kalpha2_peaks(peaks, two_theta_range)
    elif line_model != "kalpha1":
        raise ValueError(f"Unsupported line model: {line_model}")
    return sorted(_normalize_peak_intensities(peaks), key=lambda peak: peak.two_theta)


def _add_cu_kalpha2_peaks(peaks: list[Peak], two_theta_range: tuple[float, float]) -> list[Peak]:
    out = list(peaks)
    x_min, x_max = two_theta_range
    for peak in peaks:
        if not peak.d_spacing:
            continue
        arg = CU_KALPHA2 / (2.0 * peak.d_spacing)
        if arg <= 0.0 or arg > 1.0:
            continue
        two_theta = math.degrees(2.0 * math.asin(arg))
        if x_min <= two_theta <= x_max:
            out.append(
                Peak(
                    peak_id=f"{peak.peak_id}_ka2",
                    two_theta=two_theta,
                    intensity=peak.intensity * CU_KALPHA2_RATIO,
                    hkl=f"{peak.hkl} Kalpha2".strip(),
                    d_spacing=peak.d_spacing,
                    multiplicity=peak.multiplicity,
                )
            )
    return out


def _normalize_peak_intensities(peaks: list[Peak]) -> list[Peak]:
    max_i = max((peak.intensity for peak in peaks), default=0.0) or 1.0
    return [
        Peak(
            peak.peak_id,
            peak.two_theta,
            peak.intensity / max_i * 100.0,
            peak.prominence,
            peak.fwhm,
            hkl=peak.hkl,
            d_spacing=peak.d_spacing,
            multiplicity=peak.multiplicity,
        )
        for peak in peaks
    ]


def _write_peak_csv(path: Path, peaks: Iterable[Peak], source: str) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["source", "peak_id", "two_theta", "relative_intensity", "hkl", "d_spacing", "multiplicity"],
        )
        writer.writeheader()
        for peak in peaks:
            writer.writerow(
                {
                    "source": source,
                    "peak_id": peak.peak_id,
                    "two_theta": f"{peak.two_theta:.6f}",
                    "relative_intensity": f"{peak.intensity:.6f}",
                    "hkl": peak.hkl,
                    "d_spacing": "" if peak.d_spacing is None else f"{peak.d_spacing:.6f}",
                    "multiplicity": "" if peak.multiplicity is None else peak.multiplicity,
                }
            )


def _write_comparison_csv(path: Path, comparison) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "status",
                "reference_two_theta",
                "python_two_theta",
                "delta_2theta",
                "reference_intensity",
                "python_intensity",
                "intensity_ratio_python_over_reference",
                "ratio_error",
                "reference_hkl",
                "note",
            ],
        )
        writer.writeheader()
        for row in comparison.matched:
            ratio = row.get("intensity_ratio_calculated_over_reference")
            ratio_error = row.get("intensity_ratio_error")
            writer.writerow(
                {
                    "status": "matched",
                    "reference_two_theta": row.get("reference_two_theta"),
                    "python_two_theta": row.get("calculated_two_theta"),
                    "delta_2theta": row.get("delta_2theta"),
                    "reference_intensity": row.get("reference_intensity"),
                    "python_intensity": row.get("calculated_intensity"),
                    "intensity_ratio_python_over_reference": "" if ratio is None else ratio,
                    "ratio_error": "" if ratio_error is None else ratio_error,
                    "reference_hkl": row.get("reference_hkl", ""),
                    "note": _comparison_note(row),
                }
            )
        for row in comparison.missing_reference:
            writer.writerow(
                {
                    "status": "missing_reference",
                    "reference_two_theta": row.get("two_theta"),
                    "reference_intensity": row.get("intensity"),
                    "reference_hkl": row.get("hkl", ""),
                    "note": "reference peak not reproduced by Python",
                }
            )
        for row in comparison.extra_calculated:
            writer.writerow(
                {
                    "status": "extra_python",
                    "python_two_theta": row.get("two_theta"),
                    "python_intensity": row.get("intensity"),
                    "note": "Python peak not found in reference",
                }
            )


def _comparison_note(row: dict[str, object]) -> str:
    ratio_error = row.get("intensity_ratio_error")
    ref_intensity = float(row.get("reference_intensity") or 0.0)
    if ratio_error is not None and float(ratio_error) >= 2.0 and ref_intensity >= 5.0:
        return "strong intensity mismatch"
    return ""


def _plot_comparison(
    path: Path,
    python_peaks: list[Peak],
    reference_peaks: list[Peak],
    x_min: float,
    x_max: float,
) -> None:
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(12, 6.2),
        sharex=True,
        gridspec_kw={"height_ratios": [1, 1], "hspace": 0.08},
        constrained_layout=True,
    )
    for ax, peaks, title, color in [
        (axes[0], reference_peaks, "VESTA / reference peaks", "#d62728"),
        (axes[1], python_peaks, "Python pymatgen simulation", "#111111"),
    ]:
        xs = [peak.two_theta for peak in peaks]
        ys = [peak.intensity for peak in peaks]
        ax.vlines(xs, 0, ys, color=color, linewidth=1.2)
        ax.set_ylim(0, 108)
        ax.set_ylabel("Intensity / a.u.")
        ax.text(0.015, 0.86, title, transform=ax.transAxes, fontsize=12, color=color)
        ax.spines["top"].set_linewidth(1.0)
        ax.spines["right"].set_linewidth(1.0)
        ax.spines["bottom"].set_linewidth(1.0)
        ax.spines["left"].set_linewidth(1.0)
        ax.tick_params(direction="out", length=4, width=1)
    axes[1].set_xlim(x_min, x_max)
    axes[1].set_xlabel("2theta (degree)")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def _read_continuous_reference(path: Path, x_min: float, x_max: float) -> list[tuple[float, float]]:
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
    return rows if len(rows) >= 200 else []


def _broaden_profile(peaks: list[Peak], xs: list[float], sigma: float) -> list[float]:
    ys = [0.0 for _ in xs]
    if sigma <= 0:
        sigma = 0.055
    cutoff = 5.0 * sigma
    for peak in peaks:
        for idx, x in enumerate(xs):
            dx = x - peak.two_theta
            if abs(dx) > cutoff:
                continue
            ys[idx] += peak.intensity * math.exp(-0.5 * (dx / sigma) ** 2)
    return _normalize_values(ys)


def _normalize_values(values: list[float]) -> list[float]:
    max_value = max(values, default=0.0) or 1.0
    return [value / max_value * 100.0 for value in values]


def _write_profile_comparison(
    csv_path: Path,
    plot_path: Path,
    python_peaks: list[Peak],
    reference_path: Path,
    x_min: float,
    x_max: float,
    sigma: float,
) -> dict[str, object]:
    reference_rows = _read_continuous_reference(reference_path, x_min, x_max)
    if not reference_rows:
        csv_path.write_text("two_theta,reference_intensity,python_profile_intensity,residual\n", encoding="utf-8")
        _plot_profile(plot_path, [], [], [])
        return {"available": False, "reason": "reference file is not a continuous two-column profile"}
    xs = [x for x, _y in reference_rows]
    ref_ys = _normalize_values([y for _x, y in reference_rows])
    py_ys = _broaden_profile(python_peaks, xs, sigma)
    residuals = [py - ref for py, ref in zip(py_ys, ref_ys)]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["two_theta", "reference_intensity", "python_profile_intensity", "residual"])
        writer.writeheader()
        for x, ref, py, residual in zip(xs, ref_ys, py_ys, residuals):
            writer.writerow(
                {
                    "two_theta": f"{x:.6f}",
                    "reference_intensity": f"{ref:.6f}",
                    "python_profile_intensity": f"{py:.6f}",
                    "residual": f"{residual:.6f}",
                }
            )
    _plot_profile(plot_path, xs, ref_ys, py_ys)
    rmse = math.sqrt(sum(value * value for value in residuals) / len(residuals))
    high = [(x, residual) for x, residual in zip(xs, residuals) if x >= 60.0]
    high_rmse = math.sqrt(sum(value * value for _x, value in high) / len(high)) if high else None
    return {"available": True, "rmse": rmse, "high_angle_rmse_60_90": high_rmse, "points": len(xs)}


def _plot_profile(path: Path, xs: list[float], reference_ys: list[float], python_ys: list[float]) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(12, 6.2), sharex=True, constrained_layout=True)
    if xs:
        axes[0].plot(xs, reference_ys, color="#d62728", linewidth=1.0, label="VESTA/reference profile")
        axes[0].plot(xs, python_ys, color="#111111", linewidth=0.9, alpha=0.85, label="Python broadened profile")
        axes[1].plot(xs, [py - ref for py, ref in zip(python_ys, reference_ys)], color="#355c9c", linewidth=0.8)
    axes[0].set_ylabel("Intensity / a.u.")
    axes[1].set_ylabel("Residual")
    axes[1].set_xlabel("2theta (degree)")
    axes[0].legend(loc="upper right")
    for ax in axes:
        ax.tick_params(direction="out", length=4, width=1)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Reproduce XRD from CIF with pymatgen and compare against a reference pattern.")
    parser.add_argument("--cif", required=True, help="Input CIF file.")
    parser.add_argument("--reference", required=True, help="VESTA/reference peak file: csv, int, txt, or tabular text.")
    parser.add_argument("--out", required=True, help="Output folder.")
    parser.add_argument("--two-theta-min", type=float, default=10.0)
    parser.add_argument("--two-theta-max", type=float, default=90.0)
    parser.add_argument("--wavelength", type=float, default=1.54056)
    parser.add_argument("--tolerance", type=float, default=0.10)
    parser.add_argument("--line-model", choices=("kalpha1", "cu_kalpha12"), default="kalpha1")
    parser.add_argument("--profile-sigma", type=float, default=0.055)
    parser.add_argument("--experimental", default="")
    args = parser.parse_args(argv)
    outputs = run_xrd_reproduction(
        args.cif,
        args.reference,
        args.out,
        two_theta_min=args.two_theta_min,
        two_theta_max=args.two_theta_max,
        wavelength=args.wavelength,
        tolerance=args.tolerance,
        line_model=args.line_model,
        profile_sigma=args.profile_sigma,
        experimental_file=args.experimental or None,
    )
    print(json.dumps(outputs, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
