from __future__ import annotations

import shutil
from pathlib import Path

from ophiuchus.xrd.candidates import formula_elements, formula_from_name, structure_hash

from .database import StructureLibrary
from .models import StructureEntry
from .providers import StructureProvider


class LocalFolderProvider(StructureProvider):
    name = "local"
    provider_type = "local_import"

    def search_by_elements(self, elements, max_elements, include_subsystems, filters):
        return []

    def fetch_structure(self, entry_id: str):
        raise KeyError("LocalFolderProvider imports folders directly; fetch_structure is not used.")

    def import_folder(
        self,
        folder: str | Path,
        library: StructureLibrary,
        source: str = "local",
        access_note: str = "User-provided local CIF; verify source license before reuse.",
    ) -> dict[str, object]:
        root = Path(folder)
        summary: dict[str, object] = {"scanned": 0, "imported": 0, "skipped": 0, "warnings": []}
        if not root.is_dir():
            raise FileNotFoundError(f"local CIF folder does not exist: {root}")
        for path in sorted(root.rglob("*.cif")):
            summary["scanned"] += 1
            try:
                entry = self._entry_from_cif(path, library, source, access_note)
                target = library.path.parent / entry.cached_file_path
                target.parent.mkdir(parents=True, exist_ok=True)
                if not target.exists() or structure_hash(target) != entry.structure_hash:
                    shutil.copy2(path, target)
                library.upsert_structure(entry)
                summary["imported"] += 1
            except Exception as exc:
                summary["skipped"] += 1
                summary["warnings"].append(f"{path}: {exc}")
        return summary

    def _entry_from_cif(
        self,
        path: Path,
        library: StructureLibrary,
        source: str,
        access_note: str,
    ) -> StructureEntry:
        digest = structure_hash(path)
        formula = formula_from_name(path)
        elements = sorted(set(formula_elements(formula) or _elements_from_cif_text(path)))
        crystal = _crystal_metadata(path)
        if crystal.get("elements"):
            elements = list(crystal["elements"])
        internal_id = f"{source}_{digest}"
        relative_path = Path("library") / "structures" / source / f"{internal_id}.cif"
        return StructureEntry(
            internal_id=internal_id,
            source=source,
            source_id=path.name,
            formula=formula,
            reduced_formula=formula,
            elements=elements,
            cached_file_path=str(relative_path),
            structure_hash=digest,
            original_file_path=str(path),
            license_or_access_note=access_note,
            space_group_symbol=crystal.get("space_group_symbol"),
            space_group_number=crystal.get("space_group_number"),
            crystal_system=crystal.get("crystal_system"),
            lattice_parameters=crystal.get("lattice_parameters", {}),
        )


def _elements_from_cif_text(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    found: list[str] = []
    for token in text.replace("\n", " ").split():
        clean = "".join(ch for ch in token if ch.isalpha())
        if clean and clean[0].isupper() and len(clean) <= 2:
            found.append(clean)
    return found


def _crystal_metadata(path: Path) -> dict[str, object]:
    try:
        from pymatgen.core import Structure
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        structure = Structure.from_file(str(path))
        analyzer = SpacegroupAnalyzer(structure, symprec=0.01)
        lattice = structure.lattice
        return {
            "formula": structure.composition.reduced_formula,
            "elements": sorted(element.symbol for element in structure.composition.elements),
            "space_group_symbol": analyzer.get_space_group_symbol(),
            "space_group_number": analyzer.get_space_group_number(),
            "crystal_system": analyzer.get_crystal_system(),
            "lattice_parameters": {
                "a": float(lattice.a),
                "b": float(lattice.b),
                "c": float(lattice.c),
                "alpha": float(lattice.alpha),
                "beta": float(lattice.beta),
                "gamma": float(lattice.gamma),
                "volume": float(lattice.volume),
            },
        }
    except Exception:
        return {}
