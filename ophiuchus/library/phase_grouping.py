from __future__ import annotations

import hashlib
from dataclasses import dataclass

from ophiuchus.xrd.backend import SimulatedPattern

from .models import StructureEntry


@dataclass(frozen=True)
class PhaseCandidate:
    phase_id: str
    representative: StructureEntry
    entries: tuple[StructureEntry, ...]
    pattern: SimulatedPattern
    simulation_state: str = "unknown"

    def contains_structure(self, structure_id: str) -> bool:
        return any(entry.internal_id == structure_id for entry in self.entries)

    @property
    def display_name(self) -> str:
        formula = self.representative.reduced_formula or self.representative.formula
        symbol = self.representative.space_group_symbol or "space group unknown"
        number = self.representative.space_group_number
        return f"{formula} - {symbol}" if number is None else f"{formula} - {symbol} (No. {number})"


def group_phase_candidates(
    entries: list[StructureEntry],
    patterns_by_id: dict[str, SimulatedPattern],
    target_structure_id: str | None = None,
    simulation_state_by_id: dict[str, str] | None = None,
) -> list[PhaseCandidate]:
    simulation_state_by_id = simulation_state_by_id or {}
    grouped: list[list[StructureEntry]] = []
    for entry in entries:
        if entry.internal_id not in patterns_by_id:
            continue
        matching_group = next((group for group in grouped if _equivalent(group[0], entry, patterns_by_id)), None)
        if matching_group is None:
            grouped.append([entry])
        else:
            matching_group.append(entry)

    phases: list[PhaseCandidate] = []
    for group in grouped:
        representative = _choose_representative(group, target_structure_id)
        digest_input = "|".join(sorted(entry.internal_id for entry in group))
        phase_id = "phase_" + hashlib.sha256(digest_input.encode("utf-8")).hexdigest()[:16]
        phases.append(
            PhaseCandidate(
                phase_id=phase_id,
                representative=representative,
                entries=tuple(group),
                pattern=patterns_by_id[representative.internal_id],
                simulation_state=simulation_state_by_id.get(representative.internal_id, "unknown"),
            )
        )
    phases.sort(key=lambda phase: (not phase.contains_structure(target_structure_id or ""), phase.display_name, phase.phase_id))
    return phases


def _equivalent(
    left: StructureEntry,
    right: StructureEntry,
    patterns_by_id: dict[str, SimulatedPattern],
) -> bool:
    if left.structure_hash and left.structure_hash == right.structure_hash:
        return True
    left_formula = left.reduced_formula or left.formula
    right_formula = right.reduced_formula or right.formula
    if left_formula != right_formula:
        return False
    if left.space_group_number != right.space_group_number:
        return False
    if left.space_group_number is None:
        return False
    if not _lattice_compatible(left, right):
        return False
    return _major_peak_fingerprint(patterns_by_id[left.internal_id]) == _major_peak_fingerprint(patterns_by_id[right.internal_id])


def _major_peak_fingerprint(pattern: SimulatedPattern) -> tuple[int, ...]:
    return tuple(
        round(two_theta / 0.05)
        for two_theta, intensity in zip(pattern.two_theta_deg, pattern.normalized_intensity)
        if intensity >= 10.0
    )


def _lattice_compatible(left: StructureEntry, right: StructureEntry, relative_tolerance: float = 0.01) -> bool:
    comparable = []
    for key in ("a", "b", "c"):
        left_value = left.lattice_parameters.get(key)
        right_value = right.lattice_parameters.get(key)
        if left_value is None or right_value is None:
            continue
        scale = max(abs(float(left_value)), abs(float(right_value)), 1e-12)
        comparable.append(abs(float(left_value) - float(right_value)) / scale <= relative_tolerance)
    return all(comparable) if comparable else True


def _choose_representative(group: list[StructureEntry], target_structure_id: str | None) -> StructureEntry:
    if target_structure_id:
        for entry in group:
            if entry.internal_id == target_structure_id:
                return entry
    priority = {"target": 0, "local": 1, "materials_project": 2}
    return min(group, key=lambda entry: (priority.get(entry.source, 9), entry.source_id, entry.internal_id))
