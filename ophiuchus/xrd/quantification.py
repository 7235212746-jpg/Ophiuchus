from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, replace
from enum import Enum

import numpy as np
from scipy.signal import find_peaks

from .multiphase_models import MultiphaseRefinementResult, PhaseRefinementInput


class GateLevel(Enum):
    PASS = "pass"
    WARNING = "warning"
    FAIL = "fail"


@dataclass(frozen=True)
class GateFinding:
    code: str
    level: GateLevel
    message: str


@dataclass(frozen=True)
class QuantificationAssessment:
    level: GateLevel
    allow_weight_percent: bool
    label: str
    findings: tuple[GateFinding, ...]
    stability_ranges: dict[str, tuple[float, float]]


def weight_fractions_from_zmv(
    rows: Iterable[tuple[str, float, float, float, float]],
) -> dict[str, float]:
    weighted: dict[str, float] = {}
    for phase_id, scale, z, molar_mass, volume in rows:
        values = np.asarray([scale, z, molar_mass, volume], dtype=float)
        if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
            raise ValueError("Every quantitative phase requires positive finite scale, Z, M, and V values.")
        if phase_id in weighted:
            raise ValueError(f"Duplicate quantitative phase id: {phase_id}")
        weighted[phase_id] = float(np.prod(values))
    if not weighted:
        raise ValueError("At least one quantitative phase is required.")
    total = sum(weighted.values())
    return {phase_id: value / total * 100.0 for phase_id, value in weighted.items()}


def initial_scale_variants(
    phases: Iterable[PhaseRefinementInput],
) -> tuple[tuple[PhaseRefinementInput, ...], ...]:
    phase_list = tuple(phases)
    if not 2 <= len(phase_list) <= 4:
        raise ValueError("Initial-scale stability checks require two to four phases.")
    return (
        tuple(replace(phase, initial_scale=1.0) for phase in phase_list),
        tuple(replace(phase, initial_scale=10.0 if index == 0 else 1.0) for index, phase in enumerate(phase_list)),
        tuple(replace(phase, initial_scale=1.0 if index == 0 else 5.0) for index, phase in enumerate(phase_list)),
    )


def _strong_positive_residual_peaks(result: MultiphaseRefinementResult) -> tuple[float, ...]:
    residual = np.asarray(result.residual_intensity, dtype=float)
    centered = residual - float(np.median(residual))
    noise = 1.4826 * float(np.median(np.abs(centered)))
    net_signal = np.clip(result.observed_intensity - result.background_intensity, 0.0, None)
    signal_threshold = 0.08 * float(np.max(net_signal)) if net_signal.size else 0.0
    threshold = max(1.0, 5.0 * noise, signal_threshold)
    indices, _ = find_peaks(residual, height=threshold, prominence=max(0.5, threshold * 0.5))
    return tuple(float(result.two_theta_deg[index]) for index in indices)


def assess_quantification(
    primary: MultiphaseRefinementResult,
    repeated_results: Iterable[MultiphaseRefinementResult],
    validation_states: dict[str, str],
    *,
    stability_spread_limit_percent: float = 5.0,
) -> QuantificationAssessment:
    repeats = tuple(repeated_results)
    findings: list[GateFinding] = []
    phase_ids = tuple(phase.phase_id for phase in primary.phases)
    if len(set(phase_ids)) != len(phase_ids):
        findings.append(GateFinding("duplicate_phase", GateLevel.FAIL, "精修结果中出现重复物相编号。"))

    invalid_scales = [
        phase.formula
        for phase in primary.phases
        if not np.isfinite(phase.scale) or phase.scale <= 0.0
    ]
    if invalid_scales:
        findings.append(
            GateFinding(
                "invalid_scale",
                GateLevel.FAIL,
                "以下物相的尺度因子不是正有限数：" + "、".join(invalid_scales),
            )
        )
    invalid_zmv = [
        phase.formula
        for phase in primary.phases
        if not np.all(np.isfinite([phase.z, phase.molar_mass, phase.volume_angstrom3]))
        or min(phase.z, phase.molar_mass, phase.volume_angstrom3) <= 0.0
    ]
    if invalid_zmv:
        findings.append(
            GateFinding("invalid_zmv", GateLevel.FAIL, "以下物相缺少有效 Z/M/V：" + "、".join(invalid_zmv))
        )
    weight_sum = sum(float(phase.weight_percent) for phase in primary.phases)
    if not np.isfinite(weight_sum) or abs(weight_sum - 100.0) > 0.5:
        findings.append(
            GateFinding("weight_sum_invalid", GateLevel.FAIL, f"RIETAN 物相质量分数总和为 {weight_sum:.3f}%，未闭合到 100%。")
        )

    failed_validation = [
        phase.formula
        for phase in primary.phases
        if validation_states.get(phase.phase_id, "missing").lower() != "passed"
    ]
    if failed_validation:
        findings.append(
            GateFinding(
                "pattern_validation_failed",
                GateLevel.FAIL,
                "以下物相没有通过独立计算谱验证：" + "、".join(failed_validation),
            )
        )
    zmv_difference = float(primary.provenance.get("zmv_max_difference_wt_percent", float("inf")))
    if not np.isfinite(zmv_difference) or zmv_difference > 0.5:
        findings.append(
            GateFinding(
                "zmv_crosscheck_failed",
                GateLevel.FAIL,
                f"RIETAN 原生质量分数与 CIF-ZMV 交叉核对最大偏差为 {zmv_difference:.3f} wt%。",
            )
        )

    residual_peaks = _strong_positive_residual_peaks(primary)
    if len(residual_peaks) >= 2:
        positions = "、".join(f"{value:.2f}°" for value in residual_peaks[:6])
        findings.append(
            GateFinding(
                "unexplained_residual_group",
                GateLevel.FAIL,
                f"差值曲线仍有成组强正峰：{positions}。当前物相模型不完整。",
            )
        )

    stability_ranges: dict[str, tuple[float, float]] = {}
    all_runs = (primary, *repeats)
    if len(all_runs) < 3:
        findings.append(
            GateFinding(
                "stability_missing",
                GateLevel.FAIL,
                "尚未完成等比例、目标相占优和杂质相占优三组初始尺度稳定性检查。",
            )
        )
    else:
        if any(tuple(phase.phase_id for phase in result.phases) != phase_ids for result in all_runs):
            findings.append(
                GateFinding("stability_phase_mismatch", GateLevel.FAIL, "重复精修的物相顺序或组成不一致。")
            )
        else:
            for phase_index, phase_id in enumerate(phase_ids):
                values = [float(result.phases[phase_index].weight_percent) for result in all_runs]
                low, high = min(values), max(values)
                stability_ranges[phase_id] = (low, high)
                if high - low > stability_spread_limit_percent:
                    findings.append(
                        GateFinding(
                            "initial_scale_instability",
                            GateLevel.FAIL,
                            f"{primary.phases[phase_index].formula} 在不同初始尺度下变化 {high - low:.3f} wt%。",
                        )
                    )

    if primary.rwp_percent > 20.0:
        findings.append(
            GateFinding("high_rwp", GateLevel.WARNING, f"Rwp 为 {primary.rwp_percent:.3f}%，整谱拟合仍然较差。")
        )
    if primary.goodness_of_fit > 3.0:
        findings.append(
            GateFinding(
                "high_goodness_of_fit",
                GateLevel.WARNING,
                f"GOF 为 {primary.goodness_of_fit:.3f}，误差模型或物相模型可能不充分。",
            )
        )
    for phase in primary.phases[1:]:
        if phase.weight_percent < 0.5:
            findings.append(
                GateFinding(
                    "trace_phase",
                    GateLevel.WARNING,
                    f"{phase.formula} 低于 0.5 wt%，其定量值接近当前软件的谨慎报告边界。",
                )
            )

    if any(item.level is GateLevel.FAIL for item in findings):
        level = GateLevel.FAIL
    elif findings:
        level = GateLevel.WARNING
    else:
        level = GateLevel.PASS
    allow = level is not GateLevel.FAIL
    return QuantificationAssessment(
        level=level,
        allow_weight_percent=allow,
        label="实验性定量" if allow else "不可定量",
        findings=tuple(findings),
        stability_ranges=stability_ranges,
    )
