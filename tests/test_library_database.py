import tempfile
import unittest
from pathlib import Path

from ophiuchus.library.database import StructureLibrary
from ophiuchus.library.models import LibraryPeak, PhaseEvidenceCard, StructureEntry


class StructureLibraryDatabaseTests(unittest.TestCase):
    def _entry(self) -> StructureEntry:
        return StructureEntry(
            internal_id="local_abc",
            source="local",
            source_id="Fe.cif",
            formula="Fe",
            reduced_formula="Fe",
            elements=["Fe"],
            cached_file_path="structures/local/local_abc.cif",
            structure_hash="abc",
            original_file_path="C:/tmp/Fe.cif",
            license_or_access_note="User-provided local CIF",
        )

    def test_structure_roundtrip_and_enabled_toggle(self):
        with tempfile.TemporaryDirectory() as tmp:
            library = StructureLibrary(Path(tmp) / "library.sqlite")
            library.upsert_structure(self._entry())
            library.set_enabled("local_abc", False)
            loaded = library.get_structure("local_abc")
        self.assertEqual(loaded.formula, "Fe")
        self.assertEqual(loaded.chemical_system, "Fe")
        self.assertFalse(loaded.enabled_for_matching)

    def test_list_structures_filters_by_system_and_enabled(self):
        with tempfile.TemporaryDirectory() as tmp:
            library = StructureLibrary(Path(tmp) / "library.sqlite")
            library.upsert_structure(self._entry())
            disabled = self._entry()
            disabled.internal_id = "local_disabled"
            disabled.chemical_system = "Fe"
            disabled.structure_hash = "def"
            library.upsert_structure(disabled)
            library.set_enabled("local_disabled", False)
            enabled = library.list_structures(chemical_system="Fe", enabled_only=True)
        self.assertEqual([entry.internal_id for entry in enabled], ["local_abc"])

    def test_xrd_peak_cache_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            library = StructureLibrary(Path(tmp) / "library.sqlite")
            library.upsert_structure(self._entry())
            peaks = [
                LibraryPeak("local_abc", 22.2, 100.0, "CuKa", "settings1", peak_id="p1", hkl="100"),
                LibraryPeak("local_abc", 45.0, 50.0, "CuKa", "settings1", peak_id="p2"),
            ]
            library.store_xrd_peaks("local_abc", peaks, "CuKa", "settings1", 10.0, 90.0)
            loaded = library.load_xrd_peaks("local_abc", "CuKa", "settings1")
        self.assertEqual([(p.two_theta, p.relative_intensity, p.hkl) for p in loaded], [(22.2, 100.0, "100"), (45.0, 50.0, None)])

    def test_phase_evidence_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            library = StructureLibrary(Path(tmp) / "library.sqlite")
            card = PhaseEvidenceCard(
                evidence_id="ev1",
                evidence_type="computed_stability",
                chemical_system="Fe-Ge",
                title="MP computed stability",
                source="Materials Project",
                notes="Computed evidence only; needs experimental confirmation.",
                linked_formula="Fe2Ge",
            )
            library.add_phase_evidence(card)
            cards = library.list_phase_evidence("Fe-Ge")
        self.assertEqual(len(cards), 1)
        self.assertEqual(cards[0].linked_formula, "Fe2Ge")
        self.assertIn("needs experimental confirmation", cards[0].notes)


if __name__ == "__main__":
    unittest.main()
