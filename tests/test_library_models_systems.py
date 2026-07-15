import unittest

from ophiuchus.library.models import PhaseEvidenceCard, StructureEntry
from ophiuchus.library.systems import generate_chemical_systems


class LibraryModelsSystemsTests(unittest.TestCase):
    def test_structure_entry_normalizes_elements_and_chemical_system(self):
        entry = StructureEntry(
            internal_id="s1",
            source="local",
            source_id="FeGe.cif",
            formula="GeFe2",
            reduced_formula="Fe2Ge",
            elements=["Ge", "Fe", "Fe"],
            cached_file_path="structures/local/s1.cif",
            structure_hash="abc",
        )
        self.assertEqual(entry.elements, ["Fe", "Ge"])
        self.assertEqual(entry.chemical_system, "Fe-Ge")
        self.assertTrue(entry.enabled_for_matching)

    def test_phase_evidence_card_keeps_conservative_language(self):
        card = PhaseEvidenceCard(
            evidence_id="e1",
            evidence_type="literature_note",
            chemical_system="Fe-Ge",
            title="Fe-rich binary note",
            source="manual note",
            notes="Evidence suggests Fe-rich germanides should be considered.",
        )
        self.assertEqual(card.evidence_type, "literature_note")
        self.assertIn("suggests", card.notes)

    def test_generate_chemical_systems_normal_mode(self):
        systems = generate_chemical_systems(["Zr", "Fe", "Ge"], ["O", "C"], mode="normal")
        self.assertEqual(systems[:7], ["Fe", "Ge", "Zr", "Fe-Ge", "Fe-Zr", "Ge-Zr", "Fe-Ge-Zr"])
        self.assertIn("Fe-O", systems)
        self.assertIn("Zr-O", systems)
        self.assertNotIn("C-Fe-Ge-Zr", systems)

    def test_generate_chemical_systems_broad_is_capped(self):
        systems = generate_chemical_systems(["Zr", "Fe", "Ge"], ["O", "C", "Si"], mode="broad", max_systems=10)
        self.assertEqual(len(systems), 10)
        self.assertEqual(len(set(systems)), 10)


if __name__ == "__main__":
    unittest.main()
