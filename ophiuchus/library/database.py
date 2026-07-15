from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any

from .models import LibraryPeak, PhaseEvidenceCard, StructureEntry


LIBRARY_SCHEMA_VERSION = 1


class StructureLibrary:
    def __init__(self, path: str | Path, files_root: str | Path | None = None) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.files_root = Path(files_root) if files_root else self.path.parent / "library"
        self.files_root.mkdir(parents=True, exist_ok=True)
        self._init_db()

    def _connect(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.path, timeout=5.0)
        conn.row_factory = sqlite3.Row
        return conn

    def _init_db(self) -> None:
        conn = self._connect()
        try:
            conn.executescript(
                """
                create table if not exists structures (
                    internal_id text primary key,
                    source text not null,
                    source_id text not null,
                    formula text not null,
                    reduced_formula text not null,
                    elements_json text not null,
                    chemical_system text not null,
                    space_group_symbol text,
                    space_group_number integer,
                    crystal_system text,
                    lattice_json text not null,
                    is_experimental integer,
                    energy_above_hull real,
                    formation_energy real,
                    import_time text not null,
                    original_file_path text,
                    cached_file_path text not null,
                    local_metadata_path text,
                    structure_hash text not null,
                    license_or_access_note text not null,
                    enabled_for_matching integer not null,
                    user_label text,
                    user_note text,
                    quality_flags_json text not null
                );
                create index if not exists idx_structures_system on structures(chemical_system);
                create index if not exists idx_structures_source on structures(source);
                create table if not exists xrd_peaks (
                    id integer primary key autoincrement,
                    structure_internal_id text not null,
                    peak_id text not null,
                    two_theta real not null,
                    relative_intensity real not null,
                    hkl text,
                    d_spacing real,
                    radiation text not null,
                    settings_hash text not null,
                    two_theta_min real not null,
                    two_theta_max real not null,
                    created_time text default current_timestamp,
                    unique(structure_internal_id, peak_id, radiation, settings_hash)
                );
                create index if not exists idx_xrd_lookup on xrd_peaks(structure_internal_id, radiation, settings_hash);
                create table if not exists xrd_patterns_v2 (
                    cache_key text primary key,
                    structure_internal_id text not null,
                    cif_sha256 text not null,
                    settings_fingerprint text not null,
                    backend_name text not null,
                    backend_version text not null,
                    pattern_fingerprint text not null,
                    payload_json text not null,
                    created_time text default current_timestamp,
                    last_used_time text default current_timestamp
                );
                create index if not exists idx_xrd_v2_lookup
                on xrd_patterns_v2(structure_internal_id, cif_sha256, settings_fingerprint);
                create table if not exists phase_evidence (
                    evidence_id text primary key,
                    evidence_type text not null,
                    chemical_system text not null,
                    title text not null,
                    source text not null,
                    notes text not null,
                    linked_formula text,
                    linked_structure_id text,
                    temperature_range text,
                    composition_range text,
                    file_path text,
                    created_time text not null
                );
                create table if not exists metadata (
                    key text primary key,
                    value text not null
                );
                """
            )
            conn.execute("insert or replace into metadata(key, value) values (?, ?)", ("schema_version", str(LIBRARY_SCHEMA_VERSION)))
            columns = {row["name"] for row in conn.execute("pragma table_info(structures)").fetchall()}
            if "local_metadata_path" not in columns:
                conn.execute("alter table structures add column local_metadata_path text")
            conn.commit()
        finally:
            conn.close()

    def upsert_structure(self, entry: StructureEntry) -> None:
        conn = self._connect()
        try:
            conn.execute(
                """
                insert or replace into structures (
                    internal_id, source, source_id, formula, reduced_formula, elements_json,
                    chemical_system, space_group_symbol, space_group_number, crystal_system,
                    lattice_json, is_experimental, energy_above_hull, formation_energy,
                    import_time, original_file_path, cached_file_path, local_metadata_path, structure_hash,
                    license_or_access_note, enabled_for_matching, user_label, user_note,
                    quality_flags_json
                ) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    entry.internal_id,
                    entry.source,
                    entry.source_id,
                    entry.formula,
                    entry.reduced_formula,
                    json.dumps(entry.elements, ensure_ascii=False),
                    entry.chemical_system,
                    entry.space_group_symbol,
                    entry.space_group_number,
                    entry.crystal_system,
                    json.dumps(entry.lattice_parameters, ensure_ascii=False),
                    None if entry.is_experimental is None else int(entry.is_experimental),
                    entry.energy_above_hull,
                    entry.formation_energy,
                    entry.import_time,
                    entry.original_file_path,
                    entry.cached_file_path,
                    entry.local_metadata_path,
                    entry.structure_hash,
                    entry.license_or_access_note,
                    int(entry.enabled_for_matching),
                    entry.user_label,
                    entry.user_note,
                    json.dumps(entry.quality_flags, ensure_ascii=False),
                ),
            )
            conn.commit()
        finally:
            conn.close()

    def get_structure(self, internal_id: str) -> StructureEntry:
        conn = self._connect()
        try:
            row = conn.execute("select * from structures where internal_id = ?", (internal_id,)).fetchone()
        finally:
            conn.close()
        if row is None:
            raise KeyError(f"structure not found: {internal_id}")
        return self._entry_from_row(row)

    def list_structures(
        self,
        chemical_system: str | None = None,
        source: str | None = None,
        enabled_only: bool = False,
    ) -> list[StructureEntry]:
        clauses = []
        params: list[Any] = []
        if chemical_system:
            clauses.append("chemical_system = ?")
            params.append(chemical_system)
        if source:
            clauses.append("source = ?")
            params.append(source)
        if enabled_only:
            clauses.append("enabled_for_matching = 1")
        where = "" if not clauses else " where " + " and ".join(clauses)
        conn = self._connect()
        try:
            rows = conn.execute(f"select * from structures{where} order by formula, source, source_id", params).fetchall()
        finally:
            conn.close()
        return [self._entry_from_row(row) for row in rows]

    def set_enabled(self, internal_id: str, enabled: bool) -> None:
        conn = self._connect()
        try:
            conn.execute("update structures set enabled_for_matching = ? where internal_id = ?", (int(enabled), internal_id))
            conn.commit()
        finally:
            conn.close()

    def delete_structure(self, internal_id: str) -> None:
        conn = self._connect()
        try:
            conn.execute("delete from xrd_peaks where structure_internal_id = ?", (internal_id,))
            conn.execute("delete from structures where internal_id = ?", (internal_id,))
            conn.commit()
        finally:
            conn.close()

    def store_xrd_peaks(
        self,
        structure_internal_id: str,
        peaks: list[LibraryPeak],
        radiation: str,
        settings_hash: str,
        two_theta_min: float,
        two_theta_max: float,
    ) -> None:
        conn = self._connect()
        try:
            conn.execute(
                "delete from xrd_peaks where structure_internal_id = ? and radiation = ? and settings_hash = ?",
                (structure_internal_id, radiation, settings_hash),
            )
            for i, peak in enumerate(peaks, 1):
                conn.execute(
                    """
                    insert into xrd_peaks (
                        structure_internal_id, peak_id, two_theta, relative_intensity,
                        hkl, d_spacing, radiation, settings_hash, two_theta_min, two_theta_max
                    ) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        structure_internal_id,
                        peak.peak_id or f"theory_{i}",
                        peak.two_theta,
                        peak.relative_intensity,
                        peak.hkl,
                        peak.d_spacing,
                        radiation,
                        settings_hash,
                        two_theta_min,
                        two_theta_max,
                    ),
                )
            conn.commit()
        finally:
            conn.close()

    def load_xrd_peaks(self, structure_internal_id: str, radiation: str, settings_hash: str) -> list[LibraryPeak]:
        conn = self._connect()
        try:
            rows = conn.execute(
                """
                select * from xrd_peaks
                where structure_internal_id = ? and radiation = ? and settings_hash = ?
                order by two_theta
                """,
                (structure_internal_id, radiation, settings_hash),
            ).fetchall()
        finally:
            conn.close()
        return [
            LibraryPeak(
                structure_internal_id=row["structure_internal_id"],
                two_theta=row["two_theta"],
                relative_intensity=row["relative_intensity"],
                radiation=row["radiation"],
                settings_hash=row["settings_hash"],
                peak_id=row["peak_id"],
                hkl=row["hkl"],
                d_spacing=row["d_spacing"],
            )
            for row in rows
        ]

    def xrd_status(self, structure_internal_id: str, radiation: str, settings_hash: str) -> str:
        return "cached" if self.load_xrd_peaks(structure_internal_id, radiation, settings_hash) else "missing"

    def xrd_cache_info(self, structure_internal_id: str, radiation: str, settings_hash: str) -> dict[str, Any]:
        conn = self._connect()
        try:
            row = conn.execute(
                """
                select count(*) as peak_count, max(created_time) as last_simulated_time,
                       min(two_theta_min) as two_theta_min, max(two_theta_max) as two_theta_max
                from xrd_peaks
                where structure_internal_id = ? and radiation = ? and settings_hash = ?
                """,
                (structure_internal_id, radiation, settings_hash),
            ).fetchone()
        finally:
            conn.close()
        peak_count = int(row["peak_count"] or 0)
        return {
            "status": "cached" if peak_count else "missing",
            "peak_count": peak_count,
            "last_simulated_time": row["last_simulated_time"] or "",
            "two_theta_min": row["two_theta_min"],
            "two_theta_max": row["two_theta_max"],
        }

    def store_validated_pattern(self, cache_key: str, payload: dict[str, Any]) -> None:
        conn = self._connect()
        try:
            conn.execute(
                """
                insert or replace into xrd_patterns_v2 (
                    cache_key, structure_internal_id, cif_sha256, settings_fingerprint,
                    backend_name, backend_version, pattern_fingerprint, payload_json,
                    created_time, last_used_time
                ) values (?, ?, ?, ?, ?, ?, ?, ?, current_timestamp, current_timestamp)
                """,
                (
                    cache_key,
                    payload["structure_id"],
                    payload["cif_sha256"],
                    payload["settings_fingerprint"],
                    payload["engine_name"],
                    payload["engine_version"],
                    payload["pattern_fingerprint"],
                    json.dumps(payload, ensure_ascii=False, separators=(",", ":")),
                ),
            )
            conn.commit()
        finally:
            conn.close()

    def load_validated_pattern_payload(self, cache_key: str) -> dict[str, Any] | None:
        conn = self._connect()
        try:
            row = conn.execute(
                "select payload_json from xrd_patterns_v2 where cache_key = ?",
                (cache_key,),
            ).fetchone()
            if row is not None:
                conn.execute(
                    "update xrd_patterns_v2 set last_used_time = current_timestamp where cache_key = ?",
                    (cache_key,),
                )
                conn.commit()
        finally:
            conn.close()
        return None if row is None else json.loads(row["payload_json"])

    def validated_pattern_info(self, structure_internal_id: str, settings_fingerprint: str) -> dict[str, Any]:
        conn = self._connect()
        try:
            row = conn.execute(
                """
                select backend_name, backend_version, pattern_fingerprint,
                       created_time, last_used_time, payload_json
                from xrd_patterns_v2
                where structure_internal_id = ? and settings_fingerprint = ?
                order by created_time desc
                limit 1
                """,
                (structure_internal_id, settings_fingerprint),
            ).fetchone()
        finally:
            conn.close()
        if row is None:
            return {
                "status": "missing",
                "backend_name": "",
                "backend_version": "",
                "pattern_fingerprint": "",
                "peak_count": 0,
                "created_time": "",
                "last_used_time": "",
            }
        payload = json.loads(row["payload_json"])
        return {
            "status": "validated",
            "backend_name": row["backend_name"],
            "backend_version": row["backend_version"],
            "pattern_fingerprint": row["pattern_fingerprint"],
            "peak_count": len(payload.get("two_theta_deg", [])),
            "created_time": row["created_time"] or "",
            "last_used_time": row["last_used_time"] or "",
        }

    def add_phase_evidence(self, card: PhaseEvidenceCard) -> None:
        conn = self._connect()
        try:
            conn.execute(
                """
                insert or replace into phase_evidence (
                    evidence_id, evidence_type, chemical_system, title, source, notes,
                    linked_formula, linked_structure_id, temperature_range, composition_range,
                    file_path, created_time
                ) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    card.evidence_id,
                    card.evidence_type,
                    card.chemical_system,
                    card.title,
                    card.source,
                    card.notes,
                    card.linked_formula,
                    card.linked_structure_id,
                    card.temperature_range,
                    card.composition_range,
                    card.file_path,
                    card.created_time,
                ),
            )
            conn.commit()
        finally:
            conn.close()

    def list_phase_evidence(self, chemical_system: str | None = None) -> list[PhaseEvidenceCard]:
        conn = self._connect()
        try:
            if chemical_system:
                rows = conn.execute(
                    "select * from phase_evidence where chemical_system = ? order by created_time desc",
                    (chemical_system,),
                ).fetchall()
            else:
                rows = conn.execute("select * from phase_evidence order by created_time desc").fetchall()
        finally:
            conn.close()
        return [self._evidence_from_row(row) for row in rows]

    def _entry_from_row(self, row: sqlite3.Row) -> StructureEntry:
        return StructureEntry(
            internal_id=row["internal_id"],
            source=row["source"],
            source_id=row["source_id"],
            formula=row["formula"],
            reduced_formula=row["reduced_formula"],
            elements=json.loads(row["elements_json"]),
            chemical_system=row["chemical_system"],
            space_group_symbol=row["space_group_symbol"],
            space_group_number=row["space_group_number"],
            crystal_system=row["crystal_system"],
            lattice_parameters=json.loads(row["lattice_json"]),
            is_experimental=None if row["is_experimental"] is None else bool(row["is_experimental"]),
            energy_above_hull=row["energy_above_hull"],
            formation_energy=row["formation_energy"],
            import_time=row["import_time"],
            original_file_path=row["original_file_path"] or "",
            cached_file_path=row["cached_file_path"],
            local_metadata_path=row["local_metadata_path"] or "",
            structure_hash=row["structure_hash"],
            license_or_access_note=row["license_or_access_note"],
            enabled_for_matching=bool(row["enabled_for_matching"]),
            user_label=row["user_label"] or "",
            user_note=row["user_note"] or "",
            quality_flags=json.loads(row["quality_flags_json"]),
        )

    def _evidence_from_row(self, row: sqlite3.Row) -> PhaseEvidenceCard:
        return PhaseEvidenceCard(
            evidence_id=row["evidence_id"],
            evidence_type=row["evidence_type"],
            chemical_system=row["chemical_system"],
            title=row["title"],
            source=row["source"],
            notes=row["notes"],
            linked_formula=row["linked_formula"] or "",
            linked_structure_id=row["linked_structure_id"] or "",
            temperature_range=row["temperature_range"] or "",
            composition_range=row["composition_range"] or "",
            file_path=row["file_path"] or "",
            created_time=row["created_time"],
        )
