import tempfile
import unittest
from pathlib import Path

from ophiuchus.library.database import StructureLibrary
from ophiuchus.library.models import StructureEntry
from ophiuchus.refinement.oxide_candidates import (
    missing_common_oxide_requirements,
    select_controlled_oxide_entries,
    supplement_common_oxide_library,
)


def entry(
    internal_id: str,
    formula: str,
    elements: list[str],
    *,
    energy: float | None,
    space_group: str,
    space_group_number: int | None = None,
    enabled: bool = True,
) -> StructureEntry:
    return StructureEntry(
        internal_id=internal_id,
        source="materials_project",
        source_id=internal_id,
        formula=formula,
        reduced_formula=formula,
        elements=elements,
        cached_file_path=f"structures/{internal_id}.cif",
        structure_hash=internal_id,
        space_group_symbol=space_group,
        space_group_number=space_group_number,
        energy_above_hull=energy,
        enabled_for_matching=enabled,
    )


class ControlledOxideCandidateTests(unittest.TestCase):
    def test_selects_only_common_oxides_of_the_target_elements(self):
        entries = [
            entry("hematite", "Fe2O3", ["Fe", "O"], energy=0.0, space_group="R-3c", space_group_number=167),
            entry("zirconia", "ZrO2", ["Zr", "O"], energy=0.0, space_group="P21/c", space_group_number=14),
            entry("silica", "SiO2", ["Si", "O"], energy=0.0, space_group="P3121"),
            entry("ternary", "ZrFeO3", ["Zr", "Fe", "O"], energy=0.0, space_group="Pnma"),
            entry("metal", "Fe", ["Fe"], energy=0.0, space_group="Im-3m"),
        ]

        selected = select_controlled_oxide_entries(entries, {"Zr", "Fe", "Ge"})

        self.assertEqual([item.internal_id for item in selected], ["hematite", "zirconia"])

    def test_prefers_low_energy_distinct_polymorphs_and_caps_each_formula(self):
        entries = [
            entry("m1", "ZrO2", ["Zr", "O"], energy=0.0, space_group="P21/c", space_group_number=14),
            entry("m2", "ZrO2", ["Zr", "O"], energy=0.01, space_group="P21/c", space_group_number=14),
            entry("t", "ZrO2", ["Zr", "O"], energy=0.03, space_group="P42/nmc", space_group_number=137),
            entry("c", "ZrO2", ["Zr", "O"], energy=0.06, space_group="Fm-3m", space_group_number=225),
            entry("high", "ZrO2", ["Zr", "O"], energy=0.20, space_group="Pbca", space_group_number=61),
            entry("disabled", "ZrO2", ["Zr", "O"], energy=0.0, space_group="Pca21", space_group_number=29, enabled=False),
        ]

        selected = select_controlled_oxide_entries(
            entries,
            {"Zr"},
            max_polymorphs_per_formula=3,
        )

        self.assertEqual([item.internal_id for item in selected], ["m1", "t", "c"])

    def test_rejects_formula_match_with_nonstandard_crystal_prototype(self):
        entries = [
            entry("magnetite-wrong", "Fe3O4", ["Fe", "O"], energy=0.0, space_group="Cmcm", space_group_number=63),
            entry("magnetite", "Fe3O4", ["Fe", "O"], energy=0.01, space_group="Fd-3m", space_group_number=227),
            entry("zirconia-exotic", "ZrO2", ["Zr", "O"], energy=0.0, space_group="Pbca", space_group_number=61),
        ]

        selected = select_controlled_oxide_entries(entries, {"Fe", "Zr"})

        self.assertEqual([item.internal_id for item in selected], ["magnetite"])

    def test_keeps_standard_rocksalt_feo_despite_mp_hull_penalty(self):
        rocksalt = entry(
            "feo-rocksalt",
            "FeO",
            ["Fe", "O"],
            energy=0.4102,
            space_group="Fm-3m",
            space_group_number=225,
        )

        selected = select_controlled_oxide_entries([rocksalt], {"Fe"})

        self.assertEqual([item.internal_id for item in selected], ["feo-rocksalt"])

    def test_keeps_user_or_experimental_common_oxide_without_hull_energy(self):
        local = entry("local-ge-o2", "GeO2", ["Ge", "O"], energy=None, space_group="P42/mnm")
        local.source = "local"
        local.is_experimental = True

        selected = select_controlled_oxide_entries([local], {"Ge"})

        self.assertEqual([item.internal_id for item in selected], ["local-ge-o2"])

    def test_supplement_requests_only_missing_common_prototypes(self):
        class FakeProvider:
            def __init__(self):
                self.requested = None
                self.maximums = None

            def harvest_formulas_to_library(
                self,
                library,
                requirements,
                maximum_energy_above_hull_by_formula=None,
                progress_callback=None,
            ):
                self.requested = requirements
                self.maximums = maximum_energy_above_hull_by_formula
                return {"imported": 0, "imported_ids": [], "failed": 0, "warnings": []}

        with tempfile.TemporaryDirectory() as tmp:
            library_path = Path(tmp) / "library.sqlite"
            library = StructureLibrary(library_path)
            library.upsert_structure(
                entry("hematite", "Fe2O3", ["Fe", "O"], energy=0.0, space_group="R-3c", space_group_number=167)
            )
            provider = FakeProvider()

            before = missing_common_oxide_requirements(library, {"Fe"})
            summary = supplement_common_oxide_library(
                library_path,
                {"Fe"},
                provider=provider,
                radiation="CuKa",
                two_theta_range=(10.0, 90.0),
            )

        self.assertNotIn("Fe2O3", before)
        self.assertEqual(provider.requested, {"FeO": (225,), "Fe3O4": (227,)})
        self.assertEqual(provider.maximums["FeO"], 0.5)
        self.assertEqual(summary["requested_formulas"], ["FeO", "Fe3O4"])


if __name__ == "__main__":
    unittest.main()
