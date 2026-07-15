from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Any, Callable

from pymatgen.core import Composition

from ophiuchus.xrd.candidates import structure_hash

from .database import StructureLibrary
from .models import StructureEntry
from .providers import ProviderError, StructureProvider


MP_FIELDS = [
    "material_id",
    "formula_pretty",
    "formula_anonymous",
    "elements",
    "chemsys",
    "structure",
    "symmetry",
    "energy_above_hull",
    "formation_energy_per_atom",
    "is_stable",
]


class MaterialsProjectProvider(StructureProvider):
    name = "Materials Project"
    provider_type = "public_api"

    def __init__(
        self,
        api_key: str | None = None,
        env_path: str | Path | None = None,
        client_factory: Callable[[str], Any] | None = None,
    ) -> None:
        self.api_key = api_key if api_key is not None else self._read_api_key(env_path)
        self.client_factory = client_factory

    def _read_api_key(self, env_path: str | Path | None = None) -> str:
        env_value = os.environ.get("MP_API_KEY") or os.environ.get("PMG_MAPI_KEY") or os.environ.get("MAPI_KEY")
        if env_value:
            return env_value.strip()
        paths = []
        if env_path:
            paths.append(Path(env_path))
        else:
            paths.append(Path.cwd() / ".env")
        for path in paths:
            if not path.exists():
                continue
            for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
                clean = line.strip()
                if not clean or clean.startswith("#") or "=" not in clean:
                    continue
                key, value = clean.split("=", 1)
                if key.strip() in {"MP_API_KEY", "PMG_MAPI_KEY", "MAPI_KEY"}:
                    return value.strip().strip("'\"")
        return ""

    def is_configured(self) -> bool:
        return bool(self.api_key.strip())

    def configuration_message(self) -> str:
        if self.is_configured():
            return "Materials Project API key is configured."
        return "Set MP_API_KEY in the environment, local .env, or Ophi settings before harvesting Materials Project structures."

    def _client(self):
        if not self.is_configured():
            raise ProviderError("Materials Project API key is missing. Set MP_API_KEY or configure it in Ophi settings.")
        if self.client_factory:
            return self.client_factory(self.api_key)
        try:
            from mp_api.client import MPRester
        except Exception as exc:
            raise ProviderError(
                "Materials Project provider requires the official mp-api package. "
                "Install mp-api in the active Ophi environment before online harvesting."
            ) from exc
        return MPRester(self.api_key)

    def test_connection(self) -> dict[str, object]:
        with self._client() as client:
            client.summary.search(fields=["material_id"], chunk_size=1, num_chunks=1)
        return {"ok": True, "provider": self.name, "message": "Materials Project connection succeeded."}

    def search_by_elements(
        self,
        elements: list[str],
        max_elements: int,
        include_subsystems: bool,
        filters: dict[str, Any],
    ) -> list[Any]:
        chemsys = "-".join(sorted(elements))
        return self.search_chemical_systems([chemsys], max_entries_per_system=int(filters.get("max_entries_per_system", 50)))

    def search_chemical_systems(self, chemical_systems: list[str], max_entries_per_system: int = 50) -> list[dict[str, Any]]:
        records: list[dict[str, Any]] = []
        with self._client() as client:
            for system in chemical_systems:
                docs = client.summary.search(
                    chemsys=system,
                    fields=MP_FIELDS,
                    chunk_size=max_entries_per_system,
                    num_chunks=1,
                    all_fields=False,
                )
                records.extend(self._doc_to_dict(doc) for doc in docs[:max_entries_per_system])
        return records

    def fetch_structure(self, entry_id: str) -> Any:
        with self._client() as client:
            docs = client.summary.search(material_ids=[entry_id], fields=MP_FIELDS, chunk_size=1, num_chunks=1, all_fields=False)
        if not docs:
            raise ProviderError(f"Materials Project entry not found: {entry_id}")
        return self._doc_to_dict(docs[0])

    def harvest_to_library(
        self,
        library: StructureLibrary,
        chemical_systems: list[str],
        max_entries_per_system: int = 50,
        progress_callback: Callable[[dict[str, Any]], None] | None = None,
    ) -> dict[str, object]:
        summary: dict[str, object] = {
            "provider": "materials_project",
            "searched_systems": chemical_systems,
            "retrieved": 0,
            "imported": 0,
            "skipped_duplicates": 0,
            "failed": 0,
            "warnings": [],
        }
        existing_source_ids = {entry.source_id for entry in library.list_structures(source="materials_project")}
        try:
            with self._client() as client:
                for system in chemical_systems:
                    docs = client.summary.search(
                        chemsys=system,
                        fields=MP_FIELDS,
                        chunk_size=max_entries_per_system,
                        num_chunks=1,
                        all_fields=False,
                    )
                    system_summary = {"system": system, "retrieved": 0, "imported": 0, "skipped_duplicates": 0, "failed": 0}
                    for doc in docs[:max_entries_per_system]:
                        record = self._doc_to_dict(doc)
                        summary["retrieved"] += 1
                        system_summary["retrieved"] += 1
                        provider_id = str(record.get("material_id") or "")
                        if not provider_id:
                            summary["failed"] += 1
                            system_summary["failed"] += 1
                            summary["warnings"].append(f"{system}: record without material_id skipped")
                            continue
                        if provider_id in existing_source_ids:
                            summary["skipped_duplicates"] += 1
                            system_summary["skipped_duplicates"] += 1
                            continue
                        try:
                            entry = self._save_record(library, record)
                            library.upsert_structure(entry)
                            existing_source_ids.add(provider_id)
                            summary["imported"] += 1
                            system_summary["imported"] += 1
                        except Exception as exc:
                            summary["failed"] += 1
                            system_summary["failed"] += 1
                            summary["warnings"].append(f"{provider_id}: {exc}")
                    if progress_callback:
                        progress_callback({**system_summary, "total_imported": summary["imported"], "total_retrieved": summary["retrieved"]})
        except ProviderError:
            raise
        except Exception as exc:
            raise ProviderError(f"Materials Project harvest failed: {exc}") from exc
        return summary

    def harvest_formulas_to_library(
        self,
        library: StructureLibrary,
        formula_space_groups: dict[str, tuple[int, ...]],
        *,
        maximum_energy_above_hull: float = 0.12,
        maximum_energy_above_hull_by_formula: dict[str, float] | None = None,
        progress_callback: Callable[[dict[str, Any]], None] | None = None,
    ) -> dict[str, object]:
        summary: dict[str, object] = {
            "provider": "materials_project",
            "requested_formulas": list(formula_space_groups),
            "retrieved": 0,
            "imported": 0,
            "imported_ids": [],
            "skipped_duplicates": 0,
            "skipped_prototype": 0,
            "failed": 0,
            "warnings": [],
        }
        existing_source_ids = {entry.source_id for entry in library.list_structures(source="materials_project")}
        try:
            with self._client() as client:
                for formula, allowed_space_groups in formula_space_groups.items():
                    formula_maximum_hull = float(
                        (maximum_energy_above_hull_by_formula or {}).get(formula, maximum_energy_above_hull)
                    )
                    docs = client.summary.search(
                        formula=formula,
                        energy_above_hull=(0.0, formula_maximum_hull),
                        fields=MP_FIELDS,
                        chunk_size=100,
                        num_chunks=1,
                        all_fields=False,
                    )
                    formula_imported = 0
                    for doc in docs:
                        record = self._doc_to_dict(doc)
                        summary["retrieved"] += 1
                        record_formula = str(record.get("formula_pretty") or "")
                        symmetry = record.get("symmetry") or {}
                        space_group = _int_or_none(_symmetry_value(symmetry, "number"))
                        if _reduced_formula(record_formula) != _reduced_formula(formula) or space_group not in allowed_space_groups:
                            summary["skipped_prototype"] += 1
                            continue
                        provider_id = str(record.get("material_id") or "")
                        if not provider_id:
                            summary["failed"] += 1
                            summary["warnings"].append(f"{formula}: record without material_id skipped")
                            continue
                        if provider_id in existing_source_ids:
                            summary["skipped_duplicates"] += 1
                            continue
                        try:
                            entry = self._save_record(library, record)
                            library.upsert_structure(entry)
                            existing_source_ids.add(provider_id)
                            summary["imported"] += 1
                            summary["imported_ids"].append(entry.internal_id)
                            formula_imported += 1
                        except Exception as exc:
                            summary["failed"] += 1
                            summary["warnings"].append(f"{provider_id}: {exc}")
                    if progress_callback:
                        progress_callback(
                            {
                                "formula": formula,
                                "retrieved": len(docs),
                                "imported": formula_imported,
                                "total_imported": summary["imported"],
                            }
                        )
        except ProviderError:
            raise
        except Exception as exc:
            raise ProviderError(f"Materials Project formula harvest failed: {exc}") from exc
        return summary

    def _save_record(self, library: StructureLibrary, record: dict[str, Any]) -> StructureEntry:
        material_id = str(record["material_id"])
        folder = library.path.parent / "library"
        structures_dir = folder / "structures" / "materials_project"
        metadata_dir = folder / "metadata" / "materials_project"
        structures_dir.mkdir(parents=True, exist_ok=True)
        metadata_dir.mkdir(parents=True, exist_ok=True)
        cif_path = structures_dir / f"{material_id}.cif"
        metadata_path = metadata_dir / f"{material_id}.json"
        structure = record.get("structure")
        if structure is None:
            raise ValueError("record has no structure")
        if hasattr(structure, "to_file"):
            structure.to_file(str(cif_path), fmt="cif")
        elif hasattr(structure, "to"):
            cif_path.write_text(structure.to(fmt="cif"), encoding="utf-8")
        else:
            raise ValueError("structure object cannot be exported as CIF")
        metadata = self._metadata(record)
        metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
        elements = [str(item) for item in record.get("elements") or []]
        if not elements and hasattr(structure, "composition"):
            elements = sorted(str(el) for el in structure.composition.elements)
        symmetry = record.get("symmetry") or {}
        return StructureEntry(
            internal_id=f"mp_{material_id.replace('-', '_')}",
            source="materials_project",
            source_id=material_id,
            formula=str(record.get("formula_pretty") or material_id),
            reduced_formula=str(record.get("formula_pretty") or material_id),
            elements=elements,
            chemical_system=str(record.get("chemsys") or "-".join(sorted(elements))),
            cached_file_path=str(cif_path.relative_to(library.path.parent)),
            local_metadata_path=str(metadata_path.relative_to(library.path.parent)),
            structure_hash=structure_hash(cif_path),
            original_file_path=f"https://materialsproject.org/materials/{material_id}",
            license_or_access_note="Materials Project public API; requires user-provided free API key. Verify MP terms for reuse.",
            space_group_symbol=_symmetry_value(symmetry, "symbol"),
            space_group_number=_int_or_none(_symmetry_value(symmetry, "number")),
            crystal_system=_symmetry_value(symmetry, "crystal_system"),
            is_experimental=False,
            energy_above_hull=_float_or_none(record.get("energy_above_hull")),
            formation_energy=_float_or_none(record.get("formation_energy_per_atom")),
            user_note="stable" if record.get("is_stable") is True else "",
        )

    def _metadata(self, record: dict[str, Any]) -> dict[str, Any]:
        payload = {}
        for key, value in record.items():
            if key == "structure":
                continue
            try:
                json.dumps(value)
                payload[key] = value
            except TypeError:
                payload[key] = str(value)
        payload["retrieved_at"] = datetime.utcnow().isoformat(timespec="seconds") + "Z"
        payload["provider"] = "materials_project"
        return payload

    def _doc_to_dict(self, doc: Any) -> dict[str, Any]:
        if isinstance(doc, dict):
            return doc
        if hasattr(doc, "model_dump"):
            data = doc.model_dump()
            for field in MP_FIELDS:
                if hasattr(doc, field):
                    data[field] = getattr(doc, field)
            return data
        if hasattr(doc, "dict"):
            data = doc.dict()
            for field in MP_FIELDS:
                if hasattr(doc, field):
                    data[field] = getattr(doc, field)
            return data
        data = {}
        for field in MP_FIELDS:
            if hasattr(doc, field):
                data[field] = getattr(doc, field)
        return data


def _symmetry_value(symmetry: Any, key: str) -> Any:
    if isinstance(symmetry, dict):
        value = symmetry.get(key)
    else:
        value = getattr(symmetry, key, None)
    if value is None:
        return None
    if hasattr(value, "value"):
        return str(value.value)
    return str(value)


def _float_or_none(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _reduced_formula(formula: str) -> str:
    try:
        return Composition(formula).reduced_formula
    except Exception:
        return formula.replace(" ", "")


def _int_or_none(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None
