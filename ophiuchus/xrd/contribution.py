from __future__ import annotations

from dataclasses import dataclass

from .models import MultiPhaseExplanation


CONTRIBUTION_WARNING = (
    "Relative contribution proxy from peak matching / pattern screening, not Rietveld weight fraction."
)


@dataclass(frozen=True)
class ContributionProxy:
    contributions: dict[str, float]
    warning: str = CONTRIBUTION_WARNING


def contribution_proxy(explanation: MultiPhaseExplanation) -> ContributionProxy:
    raw: dict[str, float] = {}
    for score in explanation.selected_candidates:
        key = score.candidate.formula_pretty
        raw[key] = raw.get(key, 0.0) + sum(peak.intensity for peak in score.explained_experimental_peaks)
    total = sum(raw.values())
    if total <= 0:
        return ContributionProxy({})
    return ContributionProxy({key: value / total * 100.0 for key, value in raw.items()})
