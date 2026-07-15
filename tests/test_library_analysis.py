import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from ophiuchus.library.analysis import run_library_analysis
from ophiuchus.library.database import StructureLibrary
from ophiuchus.library.local_provider import LocalFolderProvider
from ophiuchus.library.xrd_cache import build_library_xrd_cache


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

CU_CIF = """
data_Cu
_cell_length_a 3.600
_cell_length_b 3.600
_cell_length_c 3.600
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
Cu1 1.0 0.0 0.0 0.0 Cu
"""


class LibraryAnalysisTests(unittest.TestCase):
    def test_run_library_analysis_uses_enabled_cached_structures(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cifs = root / "cifs"
            cifs.mkdir()
            (cifs / "Fe.cif").write_text(SIMPLE_CIF, encoding="utf-8")
            xrd = root / "sample.xy"
            rows = []
            for i in range(2000):
                x = 10.0 + i * 0.03
                y = 1.0 + 100.0 * pow(2.718281828, -((x - 22.2) ** 2) / (2 * 0.06**2))
                rows.append(f"{x:.3f} {y:.5f}")
            xrd.write_text("\n".join(rows), encoding="utf-8")
            library = StructureLibrary(root / "library.sqlite")
            LocalFolderProvider().import_folder(cifs, library)
            build_library_xrd_cache(library)
            out = root / "out"
            result = run_library_analysis(
                xrd_file=xrd,
                library_db=library.path,
                elements=["Fe"],
                out_dir=out,
                project="FeProject",
                sample_id="FeSample",
            )
            payload = json.loads(Path(result.outputs["json"]).read_text(encoding="utf-8"))
            self.assertTrue(Path(result.outputs["xrd_presentation_plot"]).exists())
            self.assertTrue(Path(result.outputs["xrd_diagnostic_plot"]).exists())
        self.assertTrue(result.top_scores)
        self.assertEqual(result.top_scores[0].candidate.source, "library:local")
        self.assertIn("library_entry_ids", payload["input"])
        self.assertTrue(payload["input"]["library_entry_ids"])
        self.assertIsNotNone(result.context)
        self.assertEqual(result.context.radiation, "CuKa")
        self.assertAlmostEqual(result.context.wavelength_angstrom, 1.54056)
        self.assertAlmostEqual(float(max(result.context.intensity)), 99.62071, places=4)

    def test_run_library_analysis_builds_cache_only_for_confirmed_element_scope(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cifs = root / "cifs"
            cifs.mkdir()
            (cifs / "Fe.cif").write_text(SIMPLE_CIF, encoding="utf-8")
            (cifs / "Cu.cif").write_text(CU_CIF, encoding="utf-8")
            xrd = root / "sample.xy"
            xrd.write_text("\n".join(f"{10 + i * 0.03:.3f} {1.0:.3f}" for i in range(1000)), encoding="utf-8")
            library = StructureLibrary(root / "library.sqlite")
            LocalFolderProvider().import_folder(cifs, library)
            result = run_library_analysis(
                xrd_file=xrd,
                library_db=library.path,
                elements=["Fe"],
                extra_elements=[],
                out_dir=root / "out",
            )
            payload = json.loads(Path(result.outputs["json"]).read_text(encoding="utf-8"))
        self.assertEqual(payload["input"]["library_cache_summary"]["checked"], 1)
        self.assertEqual(payload["input"]["candidate_source_summary"]["local"]["skipped_element_filter"], 1)

    def test_library_cache_uses_canonical_backend_despite_legacy_disable_flags(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cifs = root / "cifs"
            cifs.mkdir()
            (cifs / "Fe.cif").write_text(SIMPLE_CIF, encoding="utf-8")
            library = StructureLibrary(root / "library.sqlite")
            LocalFolderProvider().import_folder(cifs, library)
            with patch.dict("os.environ", {"OPHI_DISABLE_EXTERNAL_PYMATGEN": "1", "OPHI_DISABLE_INPROCESS_PYMATGEN": "1"}, clear=False):
                summary = build_library_xrd_cache(library)
            entry = library.list_structures()[0]
            cached = library.load_xrd_peaks(entry.internal_id, "CuKa", "unused")
        self.assertEqual(summary["simulated"], 1)
        self.assertEqual(summary["failed"], 0)
        self.assertFalse(cached)
        self.assertEqual(summary["backend_name"], "Ophi Validated pymatgen XRD")


if __name__ == "__main__":
    unittest.main()
