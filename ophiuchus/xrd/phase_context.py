from __future__ import annotations

from .models import Candidate


def phase_context_cards(
    candidates: list[Candidate],
    synthesis_metadata: dict[str, str] | None = None,
) -> list[dict[str, object]]:
    metadata = synthesis_metadata or {}
    cards: list[dict[str, object]] = []
    for candidate in candidates:
        system = "-".join(sorted(set(candidate.elements)))
        suggestions = _suggestions(candidate, metadata)
        cards.append(
            {
                "formula": candidate.formula_pretty,
                "chemical_system": system,
                "data_status": "insufficient data",
                "source": "local analysis only",
                "stability_notes": "No phase-diagram data has been loaded. Ophi will not invent stability windows.",
                "synthesis_metadata": metadata,
                "suggestions": suggestions,
            }
        )
    return cards


def _suggestions(candidate: Candidate, metadata: dict[str, str]) -> list[str]:
    elements = set(candidate.elements)
    suggestions = []
    if {"Fe", "Ge"}.issubset(elements):
        suggestions.append("Fe-Ge impurity hypothesis: compare cooling rate and Fe/Ge balance before changing unrelated elements.")
    if "O" in elements:
        suggestions.append("O-containing candidate: check oxidation exposure, ampoule quality, and raw-material handling.")
    if metadata.get("temperature") or metadata.get("duration"):
        suggestions.append("Use a small annealing grid around the current condition; avoid treating one run as conclusive.")
    if not suggestions:
        suggestions.append("No thermodynamic data loaded; treat this as a candidate-screening clue and verify manually.")
    return suggestions[:3]
