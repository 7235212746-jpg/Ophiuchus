import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from ophiuchus.library.database import StructureLibrary
from ophiuchus.library.local_provider import LocalFolderProvider
from ophiuchus.library.xrd_cache import XRD_ENGINE_VERSION, build_library_xrd_cache, library_entries_to_candidates, simulation_settings_hash


SIMPLE_CIF = """
data_Fe
_cell_length_a 4.000
_cell_length_b 4.000
_cell_length_c 4.000
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_occupancy
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_type_symbol
Fe1 1.0 0.0 0.0 0.0 Fe
"""


class LibraryXrdCacheTests(unittest.TestCase):
    def _library_with_fe(self, root: Path) -> StructureLibrary:
        source = root / "cifs"
        source.mkdir()
        (source / "Fe.cif").write_text(SIMPLE_CIF, encoding="utf-8")
        library = StructureLibrary(root / "library.sqlite")
        LocalFolderProvider().import_folder(source, library)
        return library

    def test_settings_hash_is_stable(self):
        self.assertEqual(
            simulation_settings_hash("CuKa", (10.0, 90.0)),
            simulation_settings_hash("CuKa", (10, 90)),
        )

    def test_engine_version_marks_raw_pymatgen_b0_fix(self):
        self.assertIn("raw_pymatgen", XRD_ENGINE_VERSION)
        self.assertIn("kalpha12", XRD_ENGINE_VERSION)
        self.assertIn("b0", XRD_ENGINE_VERSION)

    def test_build_library_xrd_cache_and_convert_to_candidates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            library = self._library_with_fe(root)
            summary = build_library_xrd_cache(library, radiation="CuKa", two_theta_range=(10.0, 90.0))
            candidates = library_entries_to_candidates(library, ["Fe"], radiation="CuKa", two_theta_range=(10.0, 90.0))
        self.assertEqual(summary["simulated"], 1)
        self.assertEqual(len(candidates), 1)
        self.assertEqual(candidates[0].formula_pretty, "Fe")
        self.assertTrue(candidates[0].theory_peaks)
        self.assertEqual(candidates[0].source, "library:local")

    def test_build_library_xrd_cache_does_not_replace_cif_with_vesta_reference(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            library = self._library_with_fe(root)
            xrd_root = root / "XRD"
            target = xrd_root / "Fe sum"
            target.mkdir(parents=True)
            (target / "MONI.int").write_text("44.68 100\n", encoding="utf-8")
            (target / "reference_simulated_peak_positions.csv").write_text(
                "two_theta,intensity_norm\n44.68,100\n65.02,45\n",
                encoding="utf-8",
            )

            with patch.dict("os.environ", {"OPHI_VESTA_REFERENCE_DIR": str(xrd_root), "OPHI_REQUIRE_VESTA_REFERENCE": "1"}):
                summary = build_library_xrd_cache(library, radiation="CuKa", two_theta_range=(10.0, 90.0), force=True)
                candidates = library_entries_to_candidates(library, ["Fe"], radiation="CuKa", two_theta_range=(10.0, 90.0))

        self.assertEqual(summary["simulated"], 1)
        positions = [(round(p.two_theta, 2), round(p.intensity, 2)) for p in candidates[0].theory_peaks]
        self.assertNotEqual(positions, [(44.68, 100.0), (65.02, 45.0)])
        self.assertGreater(len(positions), 2)


if __name__ == "__main__":
    unittest.main()
