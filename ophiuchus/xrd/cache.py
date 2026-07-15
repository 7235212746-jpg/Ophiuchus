from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any

from .candidates import LocalCandidateProvider, simulate_or_load_peaks
from .models import Candidate, Peak


CACHE_SCHEMA_VERSION = 1
XRD_SIMULATOR_VERSION = "xrd_v3_cu_kalpha12_raw_pymatgen_b0"


class CandidateCache:
    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
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
                create table if not exists candidates (
                    candidate_id text primary key,
                    formula text not null,
                    source text not null,
                    source_path text not null,
                    elements_json text not null,
                    structure_hash text,
                    parse_status text,
                    parse_error text,
                    updated_at text default current_timestamp
                );
                create table if not exists patterns (
                    pattern_key text primary key,
                    candidate_id text not null,
                    radiation text not null,
                    two_theta_min real not null,
                    two_theta_max real not null,
                    simulator_version text not null,
                    peaks_json text not null,
                    updated_at text default current_timestamp
                );
                create table if not exists metadata (
                    key text primary key,
                    value text not null
                );
                """
            )
            conn.execute(
                "insert or replace into metadata(key, value) values (?, ?)",
                ("schema_version", str(CACHE_SCHEMA_VERSION)),
            )
            conn.commit()
        finally:
            conn.close()

    def upsert_candidate(self, candidate: Candidate) -> None:
        conn = self._connect()
        try:
            conn.execute(
                """
                insert or replace into candidates
                (candidate_id, formula, source, source_path, elements_json, structure_hash, parse_status, parse_error, updated_at)
                values (?, ?, ?, ?, ?, ?, ?, ?, current_timestamp)
                """,
                (
                    candidate.candidate_id,
                    candidate.formula_pretty,
                    candidate.source,
                    candidate.source_path,
                    json.dumps(candidate.elements, ensure_ascii=False),
                    candidate.structure_hash,
                    candidate.parse_status,
                    candidate.parse_error,
                ),
            )
            conn.commit()
        finally:
            conn.close()

    def list_candidates(self) -> list[Candidate]:
        conn = self._connect()
        try:
            rows = conn.execute("select * from candidates order by formula, source_path").fetchall()
        finally:
            conn.close()
        return [self._candidate_from_row(row) for row in rows]

    def store_pattern(
        self,
        candidate: Candidate,
        peaks: list[Peak],
        radiation: str,
        two_theta_range: tuple[float, float],
        simulator_version: str = XRD_SIMULATOR_VERSION,
    ) -> None:
        key = self.pattern_key(candidate, radiation, two_theta_range, simulator_version)
        payload = [
            {
                "peak_id": peak.peak_id,
                "two_theta": peak.two_theta,
                "intensity": peak.intensity,
                "prominence": peak.prominence,
                "fwhm": peak.fwhm,
            }
            for peak in peaks
        ]
        self.upsert_candidate(candidate)
        conn = self._connect()
        try:
            conn.execute(
                """
                insert or replace into patterns
                (pattern_key, candidate_id, radiation, two_theta_min, two_theta_max, simulator_version, peaks_json, updated_at)
                values (?, ?, ?, ?, ?, ?, ?, current_timestamp)
                """,
                (
                    key,
                    candidate.candidate_id,
                    radiation,
                    two_theta_range[0],
                    two_theta_range[1],
                    simulator_version,
                    json.dumps(payload, ensure_ascii=False),
                ),
            )
            conn.commit()
        finally:
            conn.close()

    def load_pattern(
        self,
        candidate: Candidate,
        radiation: str,
        two_theta_range: tuple[float, float],
        simulator_version: str = XRD_SIMULATOR_VERSION,
    ) -> list[Peak] | None:
        key = self.pattern_key(candidate, radiation, two_theta_range, simulator_version)
        conn = self._connect()
        try:
            row = conn.execute("select peaks_json from patterns where pattern_key = ?", (key,)).fetchone()
        finally:
            conn.close()
        if row is None:
            return None
        return [
            Peak(
                str(item.get("peak_id") or f"theory_{i}"),
                float(item["two_theta"]),
                float(item["intensity"]),
                item.get("prominence"),
                item.get("fwhm"),
            )
            for i, item in enumerate(json.loads(row["peaks_json"]), 1)
        ]

    def stats(self) -> dict[str, Any]:
        conn = self._connect()
        try:
            candidates = conn.execute("select count(*) as n from candidates").fetchone()["n"]
            patterns = conn.execute("select count(*) as n from patterns").fetchone()["n"]
            last = conn.execute("select max(updated_at) as t from patterns").fetchone()["t"]
        finally:
            conn.close()
        return {"candidates": candidates, "patterns": patterns, "last_updated": last}

    def pattern_key(
        self,
        candidate: Candidate,
        radiation: str,
        two_theta_range: tuple[float, float],
        simulator_version: str = XRD_SIMULATOR_VERSION,
    ) -> str:
        source_hash = candidate.structure_hash or candidate.candidate_id
        return f"{candidate.candidate_id}|{source_hash}|{radiation}|{two_theta_range[0]:.3f}|{two_theta_range[1]:.3f}|{simulator_version}"

    def _candidate_from_row(self, row: sqlite3.Row) -> Candidate:
        return Candidate(
            candidate_id=row["candidate_id"],
            formula_pretty=row["formula"],
            source=row["source"],
            source_path=row["source_path"],
            elements=json.loads(row["elements_json"]),
            structure_hash=row["structure_hash"],
            parse_status=row["parse_status"] or "ok",
            parse_error=row["parse_error"],
        )


def build_candidate_cache(
    cache: CandidateCache,
    candidate_dirs: list[str | Path],
    allowed_elements: set[str],
    radiation: str,
    two_theta_range: tuple[float, float],
) -> dict[str, Any]:
    provider = LocalCandidateProvider(candidate_dirs, allowed_elements=allowed_elements or None)
    summary = {"scanned_candidates": 0, "stored_patterns": 0, "warnings": []}
    for candidate in provider.iter_candidates():
        summary["scanned_candidates"] += 1
        cache.upsert_candidate(candidate)
        cached = cache.load_pattern(candidate, radiation, two_theta_range)
        if cached:
            continue
        try:
            peaks = simulate_or_load_peaks(candidate, radiation=radiation, two_theta_range=two_theta_range)
            if peaks:
                cache.store_pattern(candidate, peaks, radiation, two_theta_range)
                summary["stored_patterns"] += 1
        except Exception as exc:
            summary["warnings"].append(f"{candidate.formula_pretty}: {exc}")
    return summary
