import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from pymatgen.core import Lattice, Structure

from ophiuchus.library.database import StructureLibrary
from ophiuchus.library.mp_provider import MaterialsProjectProvider
from ophiuchus.library.providers import ProviderError


class MaterialsProjectProviderScaffoldTests(unittest.TestCase):
    def test_missing_api_key_is_not_configured(self):
        with tempfile.TemporaryDirectory() as tmp, mock.patch.dict(os.environ, {}, clear=True):
            provider = MaterialsProjectProvider(api_key=None, env_path=Path(tmp) / ".env")
        self.assertFalse(provider.is_configured())
        self.assertIn("MP_API_KEY", provider.configuration_message())

    def test_search_without_api_key_raises_actionable_error(self):
        provider = MaterialsProjectProvider(api_key="")
        with self.assertRaises(ProviderError) as ctx:
            provider.search_by_elements(["Fe", "Ge"], max_elements=2, include_subsystems=True, filters={})
        self.assertIn("Materials Project API key", str(ctx.exception))

    def test_api_key_can_come_from_environment(self):
        with mock.patch.dict(os.environ, {"MP_API_KEY": "abc123"}, clear=True):
            provider = MaterialsProjectProvider()
        self.assertTrue(provider.is_configured())

    def test_api_key_can_come_from_dotenv(self):
        with tempfile.TemporaryDirectory() as tmp:
            env_path = Path(tmp) / ".env"
            env_path.write_text("MP_API_KEY=from_dotenv\n", encoding="utf-8")
            with mock.patch.dict(os.environ, {}, clear=True):
                provider = MaterialsProjectProvider(env_path=env_path)
        self.assertEqual(provider.api_key, "from_dotenv")

    def test_harvest_mocked_mp_records_writes_cif_metadata_and_skips_duplicate(self):
        structure = Structure(Lattice.cubic(4.0), ["Fe"], [[0, 0, 0]])

        class FakeSummary:
            def search(self, **kwargs):
                return [
                    {
                        "material_id": "mp-1",
                        "formula_pretty": "Fe",
                        "formula_anonymous": "A",
                        "elements": ["Fe"],
                        "chemsys": "Fe",
                        "structure": structure,
                        "energy_above_hull": 0.0,
                        "formation_energy_per_atom": -1.0,
                        "is_stable": True,
                        "symmetry": {"symbol": "Pm-3m", "number": 221, "crystal_system": "cubic"},
                    }
                ]

        class FakeClient:
            def __init__(self):
                self.summary = FakeSummary()

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            library = StructureLibrary(root / "library.sqlite")
            provider = MaterialsProjectProvider(api_key="key", client_factory=lambda key: FakeClient())
            first = provider.harvest_to_library(library, ["Fe"], max_entries_per_system=5)
            second = provider.harvest_to_library(library, ["Fe"], max_entries_per_system=5)
            entries = library.list_structures()
            cif_path = root / entries[0].cached_file_path
            metadata_path = root / entries[0].local_metadata_path
            self.assertEqual(first["imported"], 1)
            self.assertEqual(second["skipped_duplicates"], 1)
            self.assertEqual(len(entries), 1)
            self.assertTrue(cif_path.exists())
            self.assertTrue(metadata_path.exists())
            self.assertEqual(entries[0].source, "materials_project")

    def test_harvest_saves_after_each_chemical_system_and_reports_progress(self):
        structure_fe = Structure(Lattice.cubic(4.0), ["Fe"], [[0, 0, 0]])
        structure_ge = Structure(Lattice.cubic(5.0), ["Ge"], [[0, 0, 0]])

        class FakeSummary:
            def search(self, **kwargs):
                chemsys = kwargs["chemsys"]
                if chemsys == "Fe":
                    return [{"material_id": "mp-fe", "formula_pretty": "Fe", "elements": ["Fe"], "chemsys": "Fe", "structure": structure_fe}]
                if chemsys == "Ge":
                    return [{"material_id": "mp-ge", "formula_pretty": "Ge", "elements": ["Ge"], "chemsys": "Ge", "structure": structure_ge}]
                return []

        class FakeClient:
            def __init__(self):
                self.summary = FakeSummary()

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

        progress = []
        with tempfile.TemporaryDirectory() as tmp:
            library = StructureLibrary(Path(tmp) / "library.sqlite")
            provider = MaterialsProjectProvider(api_key="key", client_factory=lambda key: FakeClient())

            def on_progress(event):
                progress.append(event)
                if event["system"] == "Fe":
                    self.assertEqual(len(library.list_structures(source="materials_project")), 1)

            summary = provider.harvest_to_library(library, ["Fe", "Ge"], max_entries_per_system=5, progress_callback=on_progress)
            entries = library.list_structures(source="materials_project")

        self.assertEqual(summary["imported"], 2)
        self.assertEqual([item["system"] for item in progress], ["Fe", "Ge"])
        self.assertEqual({entry.source_id for entry in entries}, {"mp-fe", "mp-ge"})

    def test_harvest_preserves_summarydoc_structure_object_when_model_dump_serializes_dict(self):
        structure = Structure(Lattice.cubic(4.0), ["Fe"], [[0, 0, 0]])

        class FakeDoc:
            def __init__(self, structure):
                self.material_id = "mp-summarydoc"
                self.formula_pretty = "Fe"
                self.formula_anonymous = "A"
                self.elements = ["Fe"]
                self.chemsys = "Fe"
                self.symmetry = {"symbol": "Pm-3m", "number": 221, "crystal_system": "cubic"}
                self.energy_above_hull = 0.0
                self.formation_energy_per_atom = -1.0
                self.is_stable = True
                self.structure = structure

            def model_dump(self):
                return {
                    "material_id": self.material_id,
                    "formula_pretty": self.formula_pretty,
                    "formula_anonymous": self.formula_anonymous,
                    "elements": self.elements,
                    "chemsys": self.chemsys,
                    "symmetry": self.symmetry,
                    "energy_above_hull": self.energy_above_hull,
                    "formation_energy_per_atom": self.formation_energy_per_atom,
                    "is_stable": self.is_stable,
                    "structure": self.structure.as_dict(),
                }

        class FakeSummary:
            def search(self, **kwargs):
                return [FakeDoc(structure)]

        class FakeClient:
            def __init__(self):
                self.summary = FakeSummary()

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            library = StructureLibrary(root / "library.sqlite")
            provider = MaterialsProjectProvider(api_key="key", client_factory=lambda key: FakeClient())
            summary = provider.harvest_to_library(library, ["Fe"], max_entries_per_system=5)
            entries = library.list_structures()
            cif_path = root / entries[0].cached_file_path
            self.assertEqual(summary["imported"], 1)
            self.assertEqual(summary["failed"], 0)
            self.assertTrue(cif_path.exists())

    def test_formula_harvest_filters_common_prototype_and_skips_duplicates(self):
        rocksalt = Structure(Lattice.cubic(4.3), ["Fe", "O"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        distorted = Structure(Lattice.cubic(4.4), ["Fe", "O"], [[0, 0, 0], [0.25, 0.25, 0.25]])
        magnetite = Structure(Lattice.cubic(8.4), ["Fe", "O"], [[0, 0, 0], [0.375, 0.375, 0.375]])
        calls = []

        class FakeSummary:
            def search(self, **kwargs):
                calls.append(kwargs)
                formula = kwargs["formula"]
                if formula == "FeO":
                    return [
                        {"material_id": "mp-feo", "formula_pretty": "FeO", "elements": ["Fe", "O"], "chemsys": "Fe-O", "structure": rocksalt, "energy_above_hull": 0.0, "symmetry": {"symbol": "Fm-3m", "number": 225}},
                        {"material_id": "mp-feo-wrong", "formula_pretty": "FeO", "elements": ["Fe", "O"], "chemsys": "Fe-O", "structure": distorted, "energy_above_hull": 0.0, "symmetry": {"symbol": "P-3m1", "number": 164}},
                    ]
                if formula == "Fe3O4":
                    return [{"material_id": "mp-magnetite", "formula_pretty": "Fe3O4", "elements": ["Fe", "O"], "chemsys": "Fe-O", "structure": magnetite, "energy_above_hull": 0.0, "symmetry": {"symbol": "Fd-3m", "number": 227}}]
                return []

        class FakeClient:
            def __init__(self):
                self.summary = FakeSummary()

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

        with tempfile.TemporaryDirectory() as tmp:
            library = StructureLibrary(Path(tmp) / "library.sqlite")
            provider = MaterialsProjectProvider(api_key="key", client_factory=lambda key: FakeClient())
            first = provider.harvest_formulas_to_library(library, {"FeO": (225,), "Fe3O4": (227,)})
            second = provider.harvest_formulas_to_library(library, {"FeO": (225,), "Fe3O4": (227,)})
            entries = library.list_structures(source="materials_project")

        self.assertEqual(first["imported"], 2)
        self.assertEqual(first["skipped_prototype"], 1)
        self.assertEqual(set(first["imported_ids"]), {"mp_mp_feo", "mp_mp_magnetite"})
        self.assertEqual(second["skipped_duplicates"], 2)
        self.assertEqual({entry.source_id for entry in entries}, {"mp-feo", "mp-magnetite"})
        self.assertTrue(all(call["energy_above_hull"] == (0.0, 0.12) for call in calls))


if __name__ == "__main__":
    unittest.main()
