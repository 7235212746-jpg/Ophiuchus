from __future__ import annotations

from .models import Candidate, Peak
from .vesta_reference import find_local_vesta_reference, load_local_vesta_reference_peaks


def apply_vesta_trust_check(
    candidate: Candidate,
    two_theta_min: float,
    two_theta_max: float,
    position_tolerance_deg: float = 0.12,
    strong_intensity_threshold: float = 10.0,
    intensity_ratio_limit: float = 2.5,
) -> dict[str, object]:
    reference_peaks, reference_path = load_local_vesta_reference_peaks(candidate.formula_pretty, two_theta_min, two_theta_max)
    if not reference_peaks:
        status = {
            "status": "no_reference",
            "method": "vesta_reference_check",
            "reason": "no exact local VESTA reference found",
        }
        candidate.simulation_validation = status
        return status
    if candidate.parse_status == "local_vesta_reference":
        status = {
            "status": "reference_source",
            "method": "vesta_reference_check",
            "reference_path": str(reference_path),
            "reason": "candidate peaks were loaded directly from local VESTA reference",
        }
        candidate.simulation_validation = status
        return status

    calculated = [peak for peak in candidate.theory_peaks if two_theta_min <= peak.two_theta <= two_theta_max]
    reference_info = find_local_vesta_reference(candidate.formula_pretty)
    curated_major_peaks = bool(reference_info and reference_info.get("peaks"))
    comparison = compare_peaks_to_reference(
        calculated,
        reference_peaks,
        position_tolerance_deg=position_tolerance_deg,
        strong_intensity_threshold=strong_intensity_threshold,
        intensity_ratio_limit=intensity_ratio_limit,
        allow_extra_strong=curated_major_peaks,
    )
    status_name = "passed" if comparison["position_pass"] and comparison["intensity_pass"] else "failed"
    status = {
        "status": status_name,
        "method": "vesta_reference_check",
        "reference_path": str(reference_path),
        "reference_scope": "curated_major_peaks" if curated_major_peaks else "full_pattern",
        **comparison,
    }
    if status_name == "failed":
        status["reason"] = "simulated peaks disagree with local VESTA reference"
    candidate.simulation_validation = status
    return status


def compare_peaks_to_reference(
    calculated: list[Peak],
    reference: list[Peak],
    position_tolerance_deg: float = 0.12,
    strong_intensity_threshold: float = 10.0,
    intensity_ratio_limit: float = 2.5,
    allow_extra_strong: bool = False,
) -> dict[str, object]:
    strong_reference = [peak for peak in reference if peak.intensity >= strong_intensity_threshold]
    strong_calculated = [peak for peak in calculated if peak.intensity >= strong_intensity_threshold]
    matched: list[tuple[Peak, Peak, float]] = []
    used_calc_ids: set[str] = set()
    for ref in strong_reference:
        nearest = _nearest(ref, strong_calculated, used_calc_ids)
        if nearest and abs(nearest.two_theta - ref.two_theta) <= position_tolerance_deg:
            used_calc_ids.add(nearest.peak_id)
            matched.append((ref, nearest, abs(nearest.two_theta - ref.two_theta)))
    matched_ref_ids = {ref.peak_id for ref, _calc, _delta in matched}
    missing = [peak for peak in strong_reference if peak.peak_id not in matched_ref_ids]
    matched_calc_ids = {calc.peak_id for _ref, calc, _delta in matched}
    extra = [peak for peak in strong_calculated if peak.peak_id not in matched_calc_ids]
    ratio_errors = [_ratio_error(calc.intensity, ref.intensity) for ref, calc, _delta in matched if ref.intensity > 0 and calc.intensity > 0]
    max_ratio_error = max(ratio_errors, default=1.0)
    return {
        "matched_count": len(matched),
        "reference_strong_count": len(strong_reference),
        "calculated_strong_count": len(strong_calculated),
        "missing_strong_count": len(missing),
        "extra_strong_count": len(extra),
        "max_abs_delta_2theta": round(max((delta for _ref, _calc, delta in matched), default=0.0), 4),
        "max_strong_intensity_ratio_error": round(max_ratio_error, 4),
        "position_pass": len(missing) == 0 and (allow_extra_strong or len(extra) == 0),
        "intensity_pass": max_ratio_error <= intensity_ratio_limit,
        "extra_strong_used_as_failure": not allow_extra_strong,
    }


def _nearest(reference_peak: Peak, calculated: list[Peak], used_ids: set[str]) -> Peak | None:
    available = [peak for peak in calculated if peak.peak_id not in used_ids]
    if not available:
        return None
    return min(available, key=lambda peak: abs(peak.two_theta - reference_peak.two_theta))


def _ratio_error(a: float, b: float) -> float:
    high = max(abs(a), abs(b))
    low = max(min(abs(a), abs(b)), 1e-9)
    return high / low
