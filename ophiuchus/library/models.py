from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any


def normalize_elements(elements: list[str]) -> list[str]:
    return sorted({element.strip() for element in elements if element and element.strip()})


def chemical_system_from_elements(elements: list[str]) -> str:
    return "-".join(normalize_elements(elements))


@dataclass
class StructureEntry:
    internal_id: str
    source: str
    source_id: str
    formula: str
    reduced_formula: str
    elements: list[str]
    cached_file_path: str
    structure_hash: str
    chemical_system: str = ""
    space_group_symbol: str | None = None
    space_group_number: int | None = None
    crystal_system: str | None = None
    lattice_parameters: dict[str, float] = field(default_factory=dict)
    is_experimental: bool | None = None
    energy_above_hull: float | None = None
    formation_energy: float | None = None
    import_time: str = ""
    original_file_path: str = ""
    license_or_access_note: str = "User-provided or public structure; verify source license before reuse."
    local_metadata_path: str = ""
    enabled_for_matching: bool = True
    user_label: str = ""
    user_note: str = ""
    quality_flags: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.elements = normalize_elements(self.elements)
        if not self.chemical_system:
            self.chemical_system = chemical_system_from_elements(self.elements)
        if not self.import_time:
            self.import_time = datetime.utcnow().isoformat(timespec="seconds") + "Z"


@dataclass
class ProviderResult:
    provider_name: str
    provider_type: str
    source_id: str
    formula: str
    elements: list[str]
    structure_path: str
    metadata: dict[str, Any] = field(default_factory=dict)
    license_or_access_note: str = ""


@dataclass(frozen=True)
class LibraryPeak:
    structure_internal_id: str
    two_theta: float
    relative_intensity: float
    radiation: str
    settings_hash: str
    peak_id: str = ""
    hkl: str | None = None
    d_spacing: float | None = None


@dataclass
class PhaseEvidenceCard:
    evidence_id: str
    evidence_type: str
    chemical_system: str
    title: str
    source: str
    notes: str
    linked_formula: str = ""
    linked_structure_id: str = ""
    temperature_range: str = ""
    composition_range: str = ""
    file_path: str = ""
    created_time: str = ""

    def __post_init__(self) -> None:
        if not self.created_time:
            self.created_time = datetime.utcnow().isoformat(timespec="seconds") + "Z"
