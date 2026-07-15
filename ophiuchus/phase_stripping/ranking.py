from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Iterable, TYPE_CHECKING

import numpy as np
from scipy.signal import find_peaks

from ophiuchus.xrd.matching import score_candidate
from ophiuchus.xrd.models import Candidate, Peak

from .profile import project_candidate_profile

if TYPE_CHECKING:
    from .session import PhaseStrippingSession


FIXED_RANKING_SIGMA_DEG = 0.05
STRONG_PEAK_THRESHOLD = 10.0


@dataclass(frozen=True)
class CandidateEvidence:
    candidate: Candidate
    final_score: float
    peak_score: float
    strong_peak_coverage: float
    full_profile_improvement: float
    independent_peak_evidence: float
    shift_consistency: float
    element_reasonableness: float
    missing_strong_penalty: float
    over_subtraction_penalty: float
    warnings: tuple[str, ...]
    explanations: tuple[str, ...]
    provenance_ids: tuple[str, ...]


@dataclass
class _CandidateGroup:
    candidate: Candidate
    members: list[Candidate]
    candidate_ids: list[str]
    phase_entry_ids: list[str]

    @property
    def provenance_ids(self) -> tuple[str, ...]:
        return tuple(_unique([*self.candidate_ids, *self.phase_entry_ids]))

    def merge(self, candidate: Candidate) -> None:
        self.members.append(candidate)
        self.candidate_ids.append(candidate.candidate_id)
        self.phase_entry_ids.extend(candidate.phase_entry_ids)
        self.candidate = replace(self.candidate, phase_entry_ids=_unique(self.phase_entry_ids))


def rank_candidates(
    session: PhaseStrippingSession,
    candidates: Iterable[Candidate],
    element_scope: Iterable[str] | None = None,
) -> list[CandidateEvidence]:
    """Rank unique candidates against the current signed session residual."""
    candidate_list = list(candidates)
    active_candidates = [
        candidate
        for candidate in candidate_list
        if candidate.candidate_id not in session.excluded_candidate_ids
    ]
    groups = _candidate_groups(active_candidates)
    residual = session.residual_y
    residual_peaks = extract_residual_peaks(session.context.x, residual)
    accepted_peak_positions = session.accepted_peak_positions
    allowed_elements = {element.strip().lower() for element in element_scope or () if element.strip()}
    evidence = [
        _rank_candidate(
            group.candidate,
            group.provenance_ids,
            residual,
            session.context.x,
            residual_peaks,
            session.context.tolerance_deg,
            accepted_peak_positions,
            allowed_elements,
        )
        for group in groups
    ]
    return sorted(evidence, key=lambda item: (-item.final_score, item.candidate.candidate_id))


def deduplicate_candidates(candidates: Iterable[Candidate]) -> list[Candidate]:
    """Collapse repeated entries while retaining their phase-entry provenance."""
    return [group.candidate for group in _candidate_groups(candidates)]


def extract_residual_peaks(x: np.ndarray, residual_y: np.ndarray) -> list[Peak]:
    """Re-extract positive local maxima from a signed residual without changing it."""
    x_values = np.asarray(x, dtype=float)
    residual = np.asarray(residual_y, dtype=float)
    if x_values.ndim != 1 or residual.ndim != 1 or x_values.shape != residual.shape:
        raise ValueError("Residual x and intensity arrays must be one-dimensional with the same shape.")
    if x_values.size == 0:
        return []
    if not np.all(np.diff(x_values) > 0.0):
        raise ValueError("Residual x grid must be strictly increasing.")

    positive_maximum = float(np.max(residual))
    if not np.isfinite(positive_maximum) or positive_maximum <= 0.0:
        return []
    indices, _ = find_peaks(residual, height=positive_maximum * 0.02)
    return [
        Peak(peak_id=f"residual-{index}", two_theta=float(x_values[index]), intensity=float(residual[index]))
        for index in indices
    ]


def normalized_profile_fingerprint(candidate: Candidate) -> tuple[tuple[float, float], ...]:
    """Return a compact, normalized representation of the candidate's major peaks."""
    peaks = [peak for peak in candidate.theory_peaks if peak.intensity > 0.0]
    if not peaks:
        return ()
    maximum = max(peak.intensity for peak in peaks)
    return tuple(
        (round(peak.two_theta, 3), round(peak.intensity / maximum, 3))
        for peak in sorted(peaks, key=lambda item: item.two_theta)
        if peak.intensity / maximum >= 0.05
    )


def _candidate_groups(candidates: Iterable[Candidate]) -> list[_CandidateGroup]:
    groups: list[_CandidateGroup] = []
    for candidate in candidates:
        existing = next(
            (group for group in groups if any(_is_duplicate(candidate, member) for member in group.members)),
            None,
        )
        if existing is None:
            groups.append(
                _CandidateGroup(
                    candidate=replace(candidate, phase_entry_ids=list(candidate.phase_entry_ids)),
                    members=[candidate],
                    candidate_ids=[candidate.candidate_id],
                    phase_entry_ids=list(candidate.phase_entry_ids),
                )
            )
        else:
            existing.merge(candidate)
    return groups


def _is_duplicate(left: Candidate, right: Candidate) -> bool:
    if left.structure_hash and right.structure_hash and left.structure_hash == right.structure_hash:
        return True
    if _normalized_text(left.formula_pretty) != _normalized_text(right.formula_pretty):
        return False
    if _normalized_text(left.space_group_symbol) != _normalized_text(right.space_group_symbol):
        return False
    return _profiles_nearly_identical(normalized_profile_fingerprint(left), normalized_profile_fingerprint(right))


def _profiles_nearly_identical(
    left: tuple[tuple[float, float], ...],
    right: tuple[tuple[float, float], ...],
) -> bool:
    if not left or len(left) != len(right):
        return False
    return all(
        abs(left_position - right_position) <= 0.03 and abs(left_intensity - right_intensity) <= 0.08
        for (left_position, left_intensity), (right_position, right_intensity) in zip(left, right)
    )


def _rank_candidate(
    candidate: Candidate,
    provenance_ids: tuple[str, ...],
    residual: np.ndarray,
    x: np.ndarray,
    residual_peaks: list[Peak],
    tolerance_deg: float,
    accepted_peak_positions: tuple[float, ...],
    allowed_elements: set[str],
) -> CandidateEvidence:
    out_of_scope_elements = _out_of_scope_elements(candidate, allowed_elements)
    if not candidate.theory_peaks:
        warnings = ["Candidate has no canonical theory peaks."]
        if out_of_scope_elements:
            warnings.append(_element_scope_warning(out_of_scope_elements))
        return CandidateEvidence(
            candidate=candidate,
            final_score=0.0,
            peak_score=0.0,
            strong_peak_coverage=0.0,
            full_profile_improvement=0.0,
            independent_peak_evidence=0.0,
            shift_consistency=0.0,
            element_reasonableness=_element_reasonableness(candidate, allowed_elements),
            missing_strong_penalty=0.0,
            over_subtraction_penalty=0.0,
            warnings=tuple(warnings),
            explanations=("No full-profile score is available without canonical theory peaks.",),
            provenance_ids=provenance_ids,
        )

    discrete = score_candidate(
        residual_peaks,
        candidate,
        tolerance_deg=tolerance_deg,
        strong_peak_threshold=STRONG_PEAK_THRESHOLD,
    )
    matched = discrete.matched_theory_peaks
    strong_peaks = _effective_strong_peaks(candidate)
    strong_total = sum(peak.intensity for peak in strong_peaks) or 1.0
    matched_theory_ids = {match.theory_peak.peak_id for match in matched}
    matched_strong = sum(peak.intensity for peak in strong_peaks if peak.peak_id in matched_theory_ids)
    missing_strong = [peak for peak in strong_peaks if peak.peak_id not in matched_theory_ids]
    strong_coverage = _unit_interval(matched_strong / strong_total)
    missing_penalty = _unit_interval(0.7 * sum(peak.intensity for peak in missing_strong) / strong_total)
    profile_improvement, over_subtraction_penalty = _profile_metrics(residual, x, candidate)
    independent_evidence = _independent_peak_evidence(matched, accepted_peak_positions, tolerance_deg)
    shift_consistency = _shift_consistency(matched, tolerance_deg)
    element_reasonableness = _element_reasonableness(candidate, allowed_elements)

    raw_score = (
        0.24 * discrete.score
        + 0.24 * strong_coverage
        + 0.26 * profile_improvement
        + 0.16 * independent_evidence
        + 0.06 * shift_consistency
        + 0.04 * element_reasonableness
        - 0.20 * missing_penalty
        - 0.20 * over_subtraction_penalty
    )
    final_score = _unit_interval(raw_score)
    warnings = list(discrete.warnings)
    if missing_strong and not any("Missing strong" in warning for warning in warnings):
        warnings.append("Missing strong theoretical peaks reduce confidence.")
    if len(matched) <= 1:
        final_score = min(final_score * 0.35, 0.39)
        warnings.append("Candidate explains only an isolated residual peak and is ranked conservatively.")
    if over_subtraction_penalty > 0.02:
        warnings.append("Scale-only profile fit introduces negative residual area.")
    if out_of_scope_elements:
        final_score = 0.0
        warnings.append(_element_scope_warning(out_of_scope_elements))
    warnings = _unique(warnings)

    explanations = (
        f"残差离散峰得分：{discrete.score:.3f}；命中 {len(matched)} 个峰。",
        f"强峰覆盖率：{strong_coverage:.3f}；缺失强峰惩罚：{missing_penalty:.3f}。",
        f"固定峰宽、仅比例拟合的全谱改善：{profile_improvement:.3f}。",
        f"独立峰证据：{independent_evidence:.3f}；峰位偏移一致性：{shift_consistency:.3f}。",
        f"元素范围合理性：{element_reasonableness:.3f}；过度扣除惩罚：{over_subtraction_penalty:.3f}。",
    )
    return CandidateEvidence(
        candidate=candidate,
        final_score=final_score,
        peak_score=discrete.score,
        strong_peak_coverage=strong_coverage,
        full_profile_improvement=profile_improvement,
        independent_peak_evidence=independent_evidence,
        shift_consistency=shift_consistency,
        element_reasonableness=element_reasonableness,
        missing_strong_penalty=missing_penalty,
        over_subtraction_penalty=over_subtraction_penalty,
        warnings=tuple(warnings),
        explanations=explanations,
        provenance_ids=provenance_ids,
    )


def _profile_metrics(residual: np.ndarray, x: np.ndarray, candidate: Candidate) -> tuple[float, float]:
    profile = project_candidate_profile(x, candidate, shift_deg=0.0, sigma_deg=FIXED_RANKING_SIGMA_DEG)
    denominator = float(np.dot(profile, profile))
    if denominator <= 0.0:
        return 0.0, 0.0
    scale = max(0.0, float(np.dot(residual, profile) / denominator))
    updated = residual - scale * profile
    baseline_error = float(np.dot(residual, residual))
    updated_error = float(np.dot(updated, updated))
    improvement = 0.0 if baseline_error <= 0.0 else _unit_interval((baseline_error - updated_error) / baseline_error)
    negative_before = float(np.maximum(-residual, 0.0).sum())
    negative_after = float(np.maximum(-updated, 0.0).sum())
    added_negative = max(0.0, negative_after - negative_before)
    normalization = float(np.abs(residual).sum()) or 1.0
    return improvement, _unit_interval(added_negative / normalization)


def _independent_peak_evidence(
    matches: list[object],
    accepted_peak_positions: tuple[float, ...],
    tolerance_deg: float,
) -> float:
    if not matches:
        return 0.0
    matched_intensity = sum(match.theory_peak.intensity for match in matches)
    if matched_intensity <= 0.0:
        return 0.0
    independent_intensity = sum(
        match.theory_peak.intensity
        for match in matches
        if not any(abs(match.theory_peak.two_theta - position) <= tolerance_deg for position in accepted_peak_positions)
    )
    return _unit_interval(independent_intensity / matched_intensity)


def _shift_consistency(matches: list[object], tolerance_deg: float) -> float:
    if not matches:
        return 0.0
    if len(matches) == 1:
        return 0.25
    offsets = np.array(
        [match.experimental_peak.two_theta - match.theory_peak.two_theta for match in matches],
        dtype=float,
    )
    allowed_spread = max(float(tolerance_deg), 1e-9)
    return _unit_interval(1.0 - float(np.std(offsets)) / allowed_spread)


def _element_reasonableness(candidate: Candidate, allowed_elements: set[str]) -> float:
    elements = {element.strip().lower() for element in candidate.elements if element.strip()}
    if not elements:
        return 0.0
    if not allowed_elements:
        return 1.0
    return 1.0 if elements <= allowed_elements else 0.0


def _out_of_scope_elements(candidate: Candidate, allowed_elements: set[str]) -> tuple[str, ...]:
    if not allowed_elements:
        return ()
    out_of_scope = {
        element.strip().lower(): element.strip()
        for element in candidate.elements
        if element.strip() and element.strip().lower() not in allowed_elements
    }
    return tuple(out_of_scope[key] for key in sorted(out_of_scope))


def _element_scope_warning(out_of_scope_elements: tuple[str, ...]) -> str:
    return "Candidate contains elements outside current element scope: " + ", ".join(out_of_scope_elements) + "."


def _effective_strong_peaks(candidate: Candidate) -> list[Peak]:
    threshold_strong = [peak for peak in candidate.theory_peaks if peak.intensity >= STRONG_PEAK_THRESHOLD]
    if threshold_strong:
        return threshold_strong
    return [peak for peak in candidate.theory_peaks if peak.intensity > 0.0]


def _normalized_text(value: str) -> str:
    return " ".join(value.split()).lower()


def _unique(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(value for value in values if value))


def _unit_interval(value: float) -> float:
    return float(np.clip(value, 0.0, 1.0))
