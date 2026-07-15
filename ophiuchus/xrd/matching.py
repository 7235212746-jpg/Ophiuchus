from __future__ import annotations

from .models import Candidate, CandidateScore, MultiPhaseExplanation, Peak, PeakMatch


def score_candidate(
    experimental_peaks: list[Peak],
    candidate: Candidate,
    tolerance_deg: float = 0.20,
    strong_peak_threshold: float = 10.0,
) -> CandidateScore:
    theory_peaks = candidate.theory_peaks
    matched: list[PeakMatch] = []
    used_exp_ids: set[str] = set()
    for theory in theory_peaks:
        nearest = _nearest_peak(theory, experimental_peaks, used_exp_ids)
        if nearest and abs(nearest.two_theta - theory.two_theta) <= tolerance_deg:
            used_exp_ids.add(nearest.peak_id)
            matched.append(PeakMatch(theory, nearest, abs(nearest.two_theta - theory.two_theta)))
    strong_theory = [p for p in theory_peaks if p.intensity >= strong_peak_threshold]
    matched_theory_ids = {m.theory_peak.peak_id for m in matched}
    missing = [p for p in strong_theory if p.peak_id not in matched_theory_ids]
    explained = [m.experimental_peak for m in matched]
    explained_ids = {p.peak_id for p in explained}
    unmatched_exp = [p for p in experimental_peaks if p.peak_id not in explained_ids]

    total_strong_intensity = sum(p.intensity for p in strong_theory) or sum(p.intensity for p in theory_peaks) or 1.0
    matched_intensity = sum(m.theory_peak.intensity for m in matched)
    missing_intensity = sum(p.intensity for p in missing)
    matched_fraction = matched_intensity / total_strong_intensity
    missing_fraction = missing_intensity / total_strong_intensity
    multi_peak_bonus = 0.12 if len(matched) >= 3 else 0.06 if len(matched) == 2 else -0.08
    missing_penalty = 0.7 * missing_fraction
    score = matched_fraction - missing_penalty + multi_peak_bonus
    score = max(0.0, min(1.0, score))

    warnings: list[str] = []
    if len(matched) <= 1 and len(missing) >= 1:
        warnings.append("Candidate explains only an isolated peak; missing strong peaks reduce confidence.")
    if missing:
        warnings.append("Missing strong theoretical peaks reduce confidence.")
    validation = candidate.simulation_validation or {}
    if validation.get("status") == "failed":
        score = min(score, 0.29)
        warnings.append(f"Simulated peaks failed local VESTA reference check: {validation.get('reason', 'validation failed')}")
    elif validation.get("status") == "no_reference":
        warnings.append("No exact local VESTA reference was found; simulated peaks are model-only and require manual validation.")

    return CandidateScore(
        candidate=candidate,
        score=score,
        matched_theory_peaks=matched,
        missing_strong_theory_peaks=missing,
        explained_experimental_peaks=explained,
        unmatched_experimental_peaks=unmatched_exp,
        warnings=warnings,
        score_components={
            "matched_theory_peak_count": len(matched),
            "strong_theory_peak_count": len(strong_theory),
            "matched_strong_intensity_fraction": round(matched_fraction, 4),
            "missing_strong_intensity_fraction": round(missing_fraction, 4),
            "missing_strong_penalty": round(missing_penalty, 4),
            "multi_peak_bonus": round(multi_peak_bonus, 4),
            "final_score": round(score, 4),
            "confidence_label": _confidence_label(score, len(matched), len(missing), validation.get("status")),
            "simulation_validation_status": str(validation.get("status") or ""),
            "simulation_validation_reference": str(validation.get("reference_path") or ""),
        },
    )


def explain_multiphase(
    experimental_peaks: list[Peak],
    scores: list[CandidateScore],
    max_phases: int = 3,
    improvement_threshold: float = 0.05,
) -> MultiPhaseExplanation:
    selected: list[CandidateScore] = []
    explained_ids: set[str] = set()
    ordered = sorted(scores, key=lambda item: item.score, reverse=True)
    for _ in range(max_phases):
        best: CandidateScore | None = None
        best_gain = 0.0
        for score in ordered:
            if score in selected:
                continue
            new_peaks = [p for p in score.explained_experimental_peaks if p.peak_id not in explained_ids]
            weighted_gain = sum(p.intensity for p in new_peaks) / (sum(p.intensity for p in experimental_peaks) or 1.0)
            penalty = min(0.20, 0.03 * len(score.missing_strong_theory_peaks))
            gain = weighted_gain - penalty
            if gain > best_gain:
                best = score
                best_gain = gain
        if best is None or best_gain < improvement_threshold:
            break
        selected.append(best)
        explained_ids.update(p.peak_id for p in best.explained_experimental_peaks)
    explained = [p for p in experimental_peaks if p.peak_id in explained_ids]
    unexplained = [p for p in experimental_peaks if p.peak_id not in explained_ids]
    return MultiPhaseExplanation(
        selected_candidates=selected,
        explained_experimental_peaks=explained,
        unexplained_experimental_peaks=unexplained,
        warnings=["This is a heuristic screening result, not quantitative refinement."],
    )


def prioritize_impurities_after_main(
    scores: list[CandidateScore],
    experimental_peaks: list[Peak],
    main_elements: list[str] | set[str] | None = None,
    target_candidate_id: str | None = None,
    target_formula: str | None = None,
    main_score_threshold: float = 0.75,
) -> list[CandidateScore]:
    ordered = sorted(scores, key=lambda item: item.score, reverse=True)
    if len(ordered) <= 1:
        return ordered
    target_elements = set(main_elements or [])
    main = _choose_main_candidate(ordered, target_elements, target_candidate_id=target_candidate_id, target_formula=target_formula)
    explicit_target_found = bool(
        (target_candidate_id and main.candidate.candidate_id == target_candidate_id)
        or (target_formula and main.candidate.formula_pretty.strip().lower() == target_formula.strip().lower())
    )
    if main.score < main_score_threshold and not explicit_target_found:
        return ordered

    main_explained_ids = {peak.peak_id for peak in main.explained_experimental_peaks}
    residual_peaks = [peak for peak in experimental_peaks if peak.peak_id not in main_explained_ids]
    residual_ids = {peak.peak_id for peak in residual_peaks}
    residual_total = sum(peak.intensity for peak in residual_peaks) or 1.0

    def impurity_key(score: CandidateScore) -> tuple[float, float, int]:
        new_matches = [match for match in score.matched_theory_peaks if match.experimental_peak.peak_id in residual_ids]
        residual_gain = sum(match.experimental_peak.intensity for match in new_matches) / residual_total
        multi_peak_bonus = 0.08 if len(new_matches) >= 3 else 0.04 if len(new_matches) == 2 else 0.0
        missing_penalty = min(0.18, 0.015 * len(score.missing_strong_theory_peaks))
        missing_fraction = float(score.score_components.get("missing_strong_intensity_fraction", 0.0))
        validation_status = str(score.candidate.simulation_validation.get("status", "")).lower()
        if validation_status == "failed":
            return (-2.0 + score.score, score.score, len(new_matches))
        if missing_fraction >= 0.8 and score.score < 0.2:
            return (-1.0 + score.score, score.score, len(new_matches))
        pattern_completeness = max(0.05, 1.0 - missing_fraction) ** 2
        residual_score = (
            (residual_gain + multi_peak_bonus) * pattern_completeness
            + 0.08 * score.score
            - missing_penalty
        )
        return (residual_score, score.score, len(new_matches))

    impurities = sorted([score for score in ordered if score is not main], key=impurity_key, reverse=True)
    impurities = [score for score in impurities if score is not main]
    unique_formula: list[CandidateScore] = []
    repeated_formula: list[CandidateScore] = []
    seen_formulas: set[str] = set()
    for score in impurities:
        formula_key = score.candidate.formula_pretty.strip().lower()
        if formula_key in seen_formulas:
            repeated_formula.append(score)
            continue
        seen_formulas.add(formula_key)
        unique_formula.append(score)
    return [main, *unique_formula, *repeated_formula]


def _choose_main_candidate(
    scores: list[CandidateScore],
    target_elements: set[str],
    target_candidate_id: str | None = None,
    target_formula: str | None = None,
) -> CandidateScore:
    if target_candidate_id:
        for score in scores:
            if score.candidate.candidate_id == target_candidate_id:
                return score
    if target_formula:
        target_clean = target_formula.strip().lower()
        for score in scores:
            if score.candidate.formula_pretty.strip().lower() == target_clean:
                return score
    if not target_elements:
        return scores[0]
    exact = [score for score in scores if set(score.candidate.elements) == target_elements]
    if not exact:
        return scores[0]

    def main_key(score: CandidateScore) -> tuple[float, int, int]:
        local_bonus = 1 if score.candidate.source.startswith("library:local") or score.candidate.source.startswith("local") else 0
        return (score.score, local_bonus, len(score.matched_theory_peaks))

    best_exact = max(exact, key=main_key)
    if best_exact.score >= max(0.2, scores[0].score - 0.15):
        return best_exact
    return scores[0]


def _nearest_peak(theory: Peak, experimental: list[Peak], used_ids: set[str]) -> Peak | None:
    available = [peak for peak in experimental if peak.peak_id not in used_ids]
    if not available:
        return None
    return min(available, key=lambda peak: abs(peak.two_theta - theory.two_theta))


def _confidence_label(score: float, matched_count: int, missing_strong_count: int, validation_status: object = "") -> str:
    if validation_status == "failed":
        return "untrusted"
    if validation_status == "no_reference":
        return "model_only"
    if score >= 0.8 and matched_count >= 4 and missing_strong_count <= 1:
        return "high"
    if score >= 0.55 and matched_count >= 3:
        return "medium"
    if score >= 0.30 and matched_count >= 2:
        return "low"
    if matched_count:
        return "uncertain"
    return "unlikely"
