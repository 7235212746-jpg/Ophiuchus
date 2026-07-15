from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from ophiuchus.xrd.candidates import LocalCandidateProvider
from ophiuchus.xrd.models import Candidate, Peak
from ophiuchus.xrd.peaks import find_peaks
from ophiuchus.xrd.refinement import RietanRefinementResult


@dataclass(frozen=True)
class ImpuritySignalEstimate:
    candidate_id: str
    formula: str
    source_label: str
    source_path: str
    signal_share_percent: float
    weight_fraction_percent: float | None
    rwp_proxy_percent: float
    rwp_improvement: float
    shift_deg: float
    fwhm_deg: float
    evidence_label: str
    peak_count: int
    is_oxide: bool
    peak_width_at_boundary: bool
    matched_residual_peaks_deg: tuple[float, ...]
    missing_strong_theory_peaks_deg: tuple[float, ...]

    def to_chinese_evidence_text(self) -> str:
        matched = _format_positions(self.matched_residual_peaks_deg)
        missing = _format_positions(self.missing_strong_theory_peaks_deg)
        return "\n".join(
            [
                f"候选相：{self.formula}",
                f"来源：{self.source_label}",
                f"结构/峰表：{self.source_path}",
                f"判断：{self.evidence_label}",
                f"加相后加权残差改善：{self.rwp_improvement:.2f} 个百分点",
                f"拟合衍射信号占比：{self.signal_share_percent:.1f}%（非 wt%）",
                f"统一峰位偏移：{self.shift_deg:+.3f}°",
                f"统一 FWHM：{self.fwhm_deg:.3f}°",
                f"命中残差峰：{matched}",
                f"缺失理论强峰：{missing}",
            ]
        )


@dataclass(frozen=True)
class SampleConclusion:
    target_formula: str
    target_evidence_label: str
    refinement: RietanRefinementResult
    target_signal_share_percent: float
    impurity_estimates: tuple[ImpuritySignalEstimate, ...]
    competing_models: bool
    oxide_formulas_checked: tuple[str, ...] = ()
    oxide_screening_warnings: tuple[str, ...] = ()
    excluded_formulas: tuple[str, ...] = ()

    def to_chinese_text(self) -> str:
        s_text = f"；S {self.refinement.s_value:.4f}" if self.refinement.s_value is not None else ""
        lines = [
            "样品结论",
            "",
            f"目标相：{self.target_formula}",
            f"判断：{self.target_evidence_label}",
            (
                f"RIETAN：Rwp {self.refinement.rwp_percent:.3f}%；"
                f"Rp {self.refinement.rp_percent:.3f}%{s_text}"
            ),
            f"目标相正积分信号解释率约 {self.target_signal_share_percent:.1f}%（诊断量，非相含量）。",
            "",
            "杂质候选与衍射信号占比（非 wt%）",
        ]
        if self.impurity_estimates:
            for index, estimate in enumerate(self.impurity_estimates, 1):
                lines.append(
                    f"{index}. {estimate.formula}：{estimate.evidence_label}；"
                    f"拟合衍射信号占比约 {estimate.signal_share_percent:.1f}%；"
                    f"加相后加权残差改善 {estimate.rwp_improvement:.2f} 个百分点。"
                )
            if self.competing_models:
                lines.append("前几种候选改善接近，当前谱图不能把它们完全区分；第一名只是最可能解释。")
            if len(self.impurity_estimates) > 1:
                lines.append("各百分比来自逐个加相的替代模型，不能相加。")
        else:
            lines.append("未找到能稳定改善目标相残差的杂质模型。")
        if self.excluded_formulas:
            lines.append("本轮人工排除：" + "、".join(self.excluded_formulas) + "。")
        if self.oxide_formulas_checked:
            oxide_names = "、".join(self.oxide_formulas_checked)
            oxide_hits = [item.formula for item in self.impurity_estimates if item.is_oxide]
            lines.extend(["", "氧化物二次筛选", f"已检查：{oxide_names}。"])
            if oxide_hits:
                lines.append(f"能改善残差的氧化物候选：{'、'.join(dict.fromkeys(oxide_hits))}。")
            else:
                lines.append("未见能稳定改善残差的常见氧化物；这不等于排除无定形氧化层或库外氧化物。")
            if self.oxide_screening_warnings:
                lines.append("筛选限制：")
                lines.extend(f"- {warning}" for warning in self.oxide_screening_warnings)
        lines.extend(
            [
                "",
                "定量边界",
                "当前不能报告可信 wt%。上面的百分比是模型解释的积分衍射信号占比，不是质量分数。",
                "可信 wt% 需要每个物相的完整可靠 CIF 或 RIR，并在同一次多相 Rietveld 精修中联合尺度因子。",
                "",
                "判断流程",
                "1. 以实验原始计数和目标 CIF 运行 RIETAN 受约束精修，同时拟合背景与尺度。",
                "2. 锁定目标相结果，逐个加入结构库候选和实验谱同目录的人工参考峰表。",
                "3. 另以受控第二遍检查目标元素的常见低能氧化物，不把氧加入主候选全集。",
                "4. 对每个候选搜索统一峰位偏移与峰宽，以加权残差改善量排序。",
                "5. 对候选解释的积分信号计算占比，并检查竞争模型、缺失强峰和异常参数。",
            ]
        )
        return "\n".join(lines)


def discover_sibling_peak_references(
    xrd_path: str | Path,
    *,
    allowed_elements: set[str],
    target_formula: str,
) -> list[Candidate]:
    source = Path(xrd_path)
    if not source.is_file():
        return []
    provider = LocalCandidateProvider(
        [source.parent],
        allowed_elements=allowed_elements,
        include_peak_lists=True,
    )
    target_key = target_formula.strip().lower()
    references: list[Candidate] = []
    seen: set[str] = set()
    for candidate in provider.iter_candidates():
        if candidate.source != "local_peak_list" or candidate.parse_status == "failed":
            continue
        if candidate.formula_pretty.strip().lower() == target_key or len(candidate.theory_peaks) < 2:
            continue
        if candidate.candidate_id in seen:
            continue
        seen.add(candidate.candidate_id)
        references.append(candidate)
    return references


def estimate_impurity_signals(
    result: RietanRefinementResult,
    target: Candidate,
    candidates: Iterable[Candidate],
    *,
    max_results: int = 3,
) -> tuple[ImpuritySignalEstimate, ...]:
    x = result.two_theta_deg
    observed_signal = result.observed_intensity - result.background_intensity
    target_profile = np.clip(result.calculated_intensity - result.background_intensity, 0.0, None)
    weights = 1.0 / np.sqrt(np.clip(result.observed_intensity, 1.0, None))
    base_proxy, base_scale = _target_only_fit(observed_signal, target_profile, weights)
    residual_peak_positions = _residual_peak_positions(
        x,
        observed_signal - base_scale * target_profile,
    )
    target_key = target.formula_pretty.strip().lower()
    estimates_by_formula: dict[str, ImpuritySignalEstimate] = {}

    for candidate in candidates:
        formula_key = candidate.formula_pretty.strip().lower()
        if not formula_key or formula_key == target_key:
            continue
        if str(candidate.simulation_validation.get("status", "")).lower() == "failed":
            continue
        peaks = [peak for peak in candidate.theory_peaks if peak.intensity > 0.0]
        if len(peaks) < 2:
            continue
        estimate = _fit_candidate(
            x,
            observed_signal,
            target_profile,
            weights,
            base_proxy,
            base_scale,
            candidate,
            peaks,
            residual_peak_positions,
        )
        if estimate.rwp_improvement < 0.25 or estimate.signal_share_percent < 0.1:
            continue
        current = estimates_by_formula.get(formula_key)
        if current is None or estimate.rwp_improvement > current.rwp_improvement:
            estimates_by_formula[formula_key] = estimate

    ordered = sorted(
        estimates_by_formula.values(),
        key=lambda item: (item.rwp_improvement, item.signal_share_percent),
        reverse=True,
    )
    return tuple(ordered[: max(0, int(max_results))])


def build_sample_conclusion(
    result: RietanRefinementResult,
    target: Candidate,
    supporting_candidates: Iterable[Candidate],
    *,
    oxide_formulas_checked: Iterable[str] = (),
    oxide_screening_warnings: Iterable[str] = (),
    excluded_formulas: Iterable[str] = (),
) -> SampleConclusion:
    estimates = estimate_impurity_signals(result, target, supporting_candidates)
    target_share = _target_signal_share_percent(result)
    explicit_scale = result.parameters.get("scale")
    target_signal_missing = (
        not result.reflection_two_theta_deg
        or target_share < 1.0
        or (explicit_scale is not None and explicit_scale <= 0.0)
    )
    if target_signal_missing:
        target_label = "目标相信号不足，低残差不能单独证明该相存在"
    elif result.rwp_percent <= 10.0 and target_share >= 50.0:
        target_label = "目标相得到强支持，并且是当前模型中的主导晶相"
    elif result.rwp_percent <= 10.0:
        target_label = "目标相得到较强支持，但当前结果不足以证明它是主导晶相"
    elif result.rwp_percent <= 20.0:
        target_label = "目标相得到支持，但仍有明显未解释信号"
    else:
        target_label = "目标相仅得到初步支持，需要检查结构、仪器和其他物相"
    competing = len(estimates) >= 2 and estimates[0].rwp_improvement - estimates[1].rwp_improvement < 0.75
    return SampleConclusion(
        target_formula=target.formula_pretty,
        target_evidence_label=target_label,
        refinement=result,
        target_signal_share_percent=target_share,
        impurity_estimates=estimates,
        competing_models=competing,
        oxide_formulas_checked=tuple(dict.fromkeys(oxide_formulas_checked)),
        oxide_screening_warnings=tuple(dict.fromkeys(oxide_screening_warnings)),
        excluded_formulas=tuple(dict.fromkeys(excluded_formulas)),
    )


def _fit_candidate(
    x: np.ndarray,
    observed_signal: np.ndarray,
    target_profile: np.ndarray,
    weights: np.ndarray,
    base_proxy: float,
    base_scale: float,
    candidate: Candidate,
    peaks: list[Peak],
    residual_peak_positions: np.ndarray,
) -> ImpuritySignalEstimate:
    best: tuple[float, float, float, float, float, np.ndarray] | None = None
    positions = np.asarray([peak.two_theta for peak in peaks[:160]], dtype=float)
    intensities = np.asarray([peak.intensity for peak in peaks[:160]], dtype=float)
    for fwhm in np.arange(0.10, 0.301, 0.02):
        for shift in np.arange(-0.25, 0.251, 0.01):
            profile = _gaussian_stick_profile(x, positions, intensities, float(shift), float(fwhm))
            proxy, target_scale, impurity_scale = _two_component_fit(
                observed_signal,
                target_profile,
                profile,
                weights,
            )
            fitted = (proxy, float(shift), float(fwhm), target_scale, impurity_scale, profile)
            if best is None or fitted[0] < best[0]:
                best = fitted
    assert best is not None
    proxy, shift, fwhm, target_scale, impurity_scale, profile = best
    target_area = _integral(x, target_scale * target_profile)
    impurity_area = _integral(x, impurity_scale * profile)
    denominator = target_area + impurity_area
    signal_share = 100.0 * impurity_area / denominator if denominator > 0.0 else 0.0
    improvement = max(0.0, base_proxy - proxy)
    if improvement >= 2.0:
        label = "较强候选"
    elif improvement >= 0.75:
        label = "可能存在"
    else:
        label = "弱候选"
    at_width_boundary = abs(fwhm - 0.30) < 1e-8 or abs(fwhm - 0.10) < 1e-8
    if at_width_boundary:
        label = "不稳定候选（峰宽达到搜索边界）"
    shifted_positions = positions + shift
    match_tolerance = max(0.18, fwhm * 0.75)
    matched_residual = tuple(
        float(position)
        for position in residual_peak_positions
        if np.any(np.abs(shifted_positions - position) <= match_tolerance)
    )
    strong_cutoff = 0.20 * float(np.max(intensities))
    missing_strong = tuple(
        float(position)
        for position, intensity in zip(shifted_positions, intensities)
        if intensity >= strong_cutoff
        and not np.any(np.abs(residual_peak_positions - position) <= match_tolerance)
    )
    return ImpuritySignalEstimate(
        candidate_id=candidate.candidate_id,
        formula=candidate.formula_pretty,
        source_label=candidate.source,
        source_path=candidate.source_path,
        signal_share_percent=float(signal_share),
        weight_fraction_percent=None,
        rwp_proxy_percent=float(proxy),
        rwp_improvement=float(improvement),
        shift_deg=float(shift),
        fwhm_deg=float(fwhm),
        evidence_label=label,
        peak_count=len(peaks),
        is_oxide="O" in set(candidate.elements),
        peak_width_at_boundary=at_width_boundary,
        matched_residual_peaks_deg=matched_residual,
        missing_strong_theory_peaks_deg=missing_strong,
    )


def _gaussian_stick_profile(
    x: np.ndarray,
    positions: np.ndarray,
    intensities: np.ndarray,
    shift: float,
    fwhm: float,
) -> np.ndarray:
    sigma = fwhm / 2.354820045
    profile = np.zeros_like(x, dtype=float)
    for position, intensity in zip(positions, intensities):
        center = position + shift
        if center < x[0] - fwhm or center > x[-1] + fwhm:
            continue
        left = int(np.searchsorted(x, center - 5.0 * sigma, side="left"))
        right = int(np.searchsorted(x, center + 5.0 * sigma, side="right"))
        if right <= left:
            continue
        profile[left:right] += intensity * np.exp(-0.5 * ((x[left:right] - center) / sigma) ** 2)
    return profile


def _target_only_fit(
    observed_signal: np.ndarray,
    target_profile: np.ndarray,
    weights: np.ndarray,
) -> tuple[float, float]:
    weighted_target = target_profile * weights
    weighted_observed = observed_signal * weights
    denominator = float(np.dot(weighted_target, weighted_target))
    scale = max(0.0, float(np.dot(weighted_target, weighted_observed)) / denominator) if denominator else 0.0
    proxy = _weighted_proxy(observed_signal, scale * target_profile, weights)
    return proxy, scale


def _two_component_fit(
    observed_signal: np.ndarray,
    target_profile: np.ndarray,
    impurity_profile: np.ndarray,
    weights: np.ndarray,
) -> tuple[float, float, float]:
    target = target_profile * weights
    impurity = impurity_profile * weights
    observed = observed_signal * weights
    gram = np.array(
        [
            [np.dot(target, target), np.dot(target, impurity)],
            [np.dot(target, impurity), np.dot(impurity, impurity)],
        ],
        dtype=float,
    )
    rhs = np.array([np.dot(target, observed), np.dot(impurity, observed)], dtype=float)
    options: list[tuple[float, float]] = [(0.0, 0.0)]
    if gram[0, 0] > 0.0:
        options.append((max(0.0, rhs[0] / gram[0, 0]), 0.0))
    if gram[1, 1] > 0.0:
        options.append((0.0, max(0.0, rhs[1] / gram[1, 1])))
    if abs(float(np.linalg.det(gram))) > 1e-18:
        both = np.linalg.solve(gram, rhs)
        if np.all(both >= 0.0):
            options.append((float(both[0]), float(both[1])))
    best = min(
        options,
        key=lambda values: float(
            np.dot(observed - values[0] * target - values[1] * impurity, observed - values[0] * target - values[1] * impurity)
        ),
    )
    calculated = best[0] * target_profile + best[1] * impurity_profile
    return _weighted_proxy(observed_signal, calculated, weights), best[0], best[1]


def _weighted_proxy(observed: np.ndarray, calculated: np.ndarray, weights: np.ndarray) -> float:
    denominator = float(np.dot(observed * weights, observed * weights))
    if denominator <= 0.0:
        return 0.0
    residual = (observed - calculated) * weights
    return float(np.sqrt(np.dot(residual, residual) / denominator) * 100.0)


def _integral(x: np.ndarray, values: np.ndarray) -> float:
    if x.size < 2:
        return 0.0
    return float(np.sum((values[:-1] + values[1:]) * np.diff(x) * 0.5))


def _target_signal_share_percent(result: RietanRefinementResult) -> float:
    observed_signal = np.clip(result.observed_intensity - result.background_intensity, 0.0, None)
    target_signal = np.clip(result.calculated_intensity - result.background_intensity, 0.0, None)
    observed_area = _integral(result.two_theta_deg, observed_signal)
    if observed_area <= 0.0:
        return 0.0
    return float(np.clip(100.0 * _integral(result.two_theta_deg, target_signal) / observed_area, 0.0, 100.0))


def _residual_peak_positions(x: np.ndarray, residual: np.ndarray) -> np.ndarray:
    positive = np.clip(np.asarray(residual, dtype=float), 0.0, None)
    maximum = float(np.max(positive)) if positive.size else 0.0
    if maximum <= 0.0:
        return np.asarray([], dtype=float)
    normalized = 100.0 * positive / maximum
    peaks = find_peaks(
        list(zip(x, normalized)),
        two_theta_min=float(x[0]),
        two_theta_max=float(x[-1]),
        min_height=2.0,
        min_distance_deg=0.12,
        max_peaks=100,
        smooth_window=3,
        min_prominence=1.0,
    )
    return np.asarray([peak.two_theta for peak in peaks], dtype=float)


def _format_positions(values: tuple[float, ...]) -> str:
    if not values:
        return "无"
    shown = "、".join(f"{value:.2f}°" for value in values[:16])
    if len(values) > 16:
        shown += f" 等 {len(values)} 个"
    return shown
