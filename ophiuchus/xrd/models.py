from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class PatternPoint:
    two_theta: float
    intensity: float


@dataclass(frozen=True)
class Peak:
    peak_id: str
    two_theta: float
    intensity: float
    prominence: float | None = None
    fwhm: float | None = None
    hkl: str = ""
    d_spacing: float | None = None
    multiplicity: int | None = None


@dataclass
class XrdPattern:
    points: list[PatternPoint]
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class Candidate:
    candidate_id: str
    formula_pretty: str
    source: str
    source_path: str
    elements: list[str]
    structure_hash: str | None
    parse_status: str = "ok"
    parse_error: str | None = None
    theory_peaks: list[Peak] = field(default_factory=list)
    simulation_validation: dict[str, Any] = field(default_factory=dict)
    simulated_pattern: Any | None = None
    phase_entry_ids: list[str] = field(default_factory=list)
    space_group_symbol: str = ""
    space_group_number: int | None = None
    simulation_state: str = "unknown"

    @property
    def path(self) -> Path:
        return Path(self.source_path)


@dataclass(frozen=True)
class PeakMatch:
    theory_peak: Peak
    experimental_peak: Peak
    delta: float


@dataclass
class CandidateScore:
    candidate: Candidate
    score: float
    matched_theory_peaks: list[PeakMatch]
    missing_strong_theory_peaks: list[Peak]
    explained_experimental_peaks: list[Peak]
    unmatched_experimental_peaks: list[Peak]
    warnings: list[str] = field(default_factory=list)
    score_components: dict[str, float | int | str] = field(default_factory=dict)


@dataclass
class MultiPhaseExplanation:
    selected_candidates: list[CandidateScore]
    explained_experimental_peaks: list[Peak]
    unexplained_experimental_peaks: list[Peak]
    warnings: list[str] = field(default_factory=list)
