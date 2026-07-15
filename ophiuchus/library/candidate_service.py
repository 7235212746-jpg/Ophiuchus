from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from .database import StructureLibrary
from .xrd_cache import simulation_settings_fingerprint, simulation_settings_hash


@dataclass(frozen=True)
class CandidateStructureRow:
    internal_id: str
    formula: str
    chemical_system: str
    source: str
    source_id: str
    structure_path: str
    structure_hash: str
    space_group_symbol: str
    space_group_number: int | None
    enabled: bool
    cache_status: str
    backend_name: str
    backend_version: str
    pattern_fingerprint: str
    peak_count: int
    last_simulated_time: str
    candidate_for_analysis: bool
    used_in_analysis: bool
    skip_reasons: list[str]

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["skip_reasons"] = " | ".join(self.skip_reasons)
        return data


class CandidateStructureService:
    def __init__(
        self,
        library: StructureLibrary,
        radiation: str = "CuKa",
        two_theta_range: tuple[float, float] = (10.0, 90.0),
        candidate_scope: str = "subsystems",
    ) -> None:
        self.library = library
        self.radiation = radiation
        self.two_theta_range = two_theta_range
        self.candidate_scope = candidate_scope
        self.settings_hash = simulation_settings_hash(radiation, two_theta_range)
        self.settings_fingerprint = simulation_settings_fingerprint(radiation, two_theta_range)

    def list_candidate_rows(self, elements: list[str], used_ids: set[str] | None = None) -> list[CandidateStructureRow]:
        used_ids = used_ids or set()
        allowed = set(elements)
        rows: list[CandidateStructureRow] = []
        for entry in self.library.list_structures(enabled_only=False):
            cache = self.library.validated_pattern_info(entry.internal_id, self.settings_fingerprint)
            entry_elements = set(entry.elements)
            in_element_scope = (not allowed) or entry_elements.issubset(allowed)
            in_candidate_scope = in_element_scope if self.candidate_scope == "subsystems" else ((not allowed) or entry_elements == allowed)
            reasons: list[str] = []
            if not entry.enabled_for_matching:
                reasons.append("disabled by user")
            if not entry.cached_file_path:
                reasons.append("no valid CIF path")
            else:
                full_path = Path(self.library.path.parent) / entry.cached_file_path
                if not full_path.exists():
                    reasons.append("no valid CIF path")
            if not in_element_scope:
                reasons.append("outside selected element scope")
            elif not in_candidate_scope:
                reasons.append("outside candidate scope")
            if cache["status"] != "validated":
                reasons.append("validated simulation missing")
            candidate_ready = not reasons
            rows.append(
                CandidateStructureRow(
                    internal_id=entry.internal_id,
                    formula=entry.reduced_formula or entry.formula,
                    chemical_system=entry.chemical_system,
                    source=entry.source,
                    source_id=entry.source_id,
                    structure_path=str(Path(self.library.path.parent) / entry.cached_file_path),
                    structure_hash=entry.structure_hash,
                    space_group_symbol=entry.space_group_symbol or "",
                    space_group_number=entry.space_group_number,
                    enabled=entry.enabled_for_matching,
                    cache_status=str(cache["status"]),
                    backend_name=str(cache["backend_name"]),
                    backend_version=str(cache["backend_version"]),
                    pattern_fingerprint=str(cache["pattern_fingerprint"]),
                    peak_count=int(cache["peak_count"]),
                    last_simulated_time=str(cache["created_time"]),
                    candidate_for_analysis=candidate_ready,
                    used_in_analysis=entry.internal_id in used_ids,
                    skip_reasons=reasons,
                )
            )
        return rows

    def usage_summary(self, rows: list[CandidateStructureRow], used_ids: set[str] | None = None) -> dict[str, int]:
        used_ids = used_ids or {row.internal_id for row in rows if row.used_in_analysis}
        return {
            "total_structures": len(rows),
            "enabled_structures": sum(1 for row in rows if row.enabled),
            "cache_ready_structures": sum(1 for row in rows if row.cache_status == "validated"),
            "candidate_ready_structures": sum(1 for row in rows if row.candidate_for_analysis),
            "used_structures": len(used_ids),
            "skipped_structures": sum(1 for row in rows if row.internal_id not in used_ids),
            "disabled_structures": sum(1 for row in rows if not row.enabled),
            "missing_cache_structures": sum(1 for row in rows if row.cache_status != "validated"),
        }


def write_candidate_usage_summary(path: str | Path, rows: list[CandidateStructureRow]) -> str:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "internal_id",
        "formula",
        "chemical_system",
        "source",
        "source_id",
        "structure_path",
        "structure_hash",
        "space_group_symbol",
        "space_group_number",
        "enabled",
        "cache_status",
        "backend_name",
        "backend_version",
        "pattern_fingerprint",
        "peak_count",
        "last_simulated_time",
        "candidate_for_analysis",
        "used_in_analysis",
        "skip_reasons",
    ]
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row.to_dict())
    return str(out)
