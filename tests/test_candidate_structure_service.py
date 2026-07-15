import csv
import json
import tempfile
import unittest
from pathlib import Path

from ophiuchus.library.analysis import run_library_analysis
from ophiuchus.library.candidate_service import CandidateStructureService
from ophiuchus.library.database import StructureLibrary
from ophiuchus.library.inspector import export_peak_inspection, inspect_peak
from ophiuchus.library.local_provider import LocalFolderProvider
from ophiuchus.library.xrd_cache import build_library_xrd_cache, simulation_settings_hash


FE_CIF = """
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
Cu1 1.0 0.0 0.0 0.0 Cu
"""


class CandidateStructureServiceTests(unittest.TestCase):
    def _library(self, root: Path) -> StructureLibrary:
        cifs = root / "cifs"
        cifs.mkdir()
        (cifs / "Fe.cif").write_text(FE_CIF, encoding="utf-8")
        (cifs / "Cu.cif").write_text(CU_CIF, encoding="utf-8")
        library = StructureLibrary(root / "library.sqlite")
        LocalFolderProvider().import_folder(cifs, library)
        build_library_xrd_cache(library, radiation="CuKa", two_theta_range=(10.0, 90.0))
        return library

    def test_candidate_service_lists_cache_status_and_skip_reasons(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            library = self._library(root)
            cu_id = next(entry.internal_id for entry in library.list_structures() if entry.reduced_formula == "Cu")
            library.set_enabled(cu_id, False)
            service = CandidateStructureService(library, radiation="CuKa", two_theta_range=(10.0, 90.0), candidate_scope="subsystems")
            rows = service.list_candidate_rows(["Fe"])
            summary = service.usage_summary(rows, used_ids={row.internal_id for row in rows if row.formula == "Fe"})

        by_formula = {row.formula: row for row in rows}
        self.assertEqual(by_formula["Fe"].cache_status, "validated")
        self.assertEqual(by_formula["Fe"].backend_name, "Ophi Validated pymatgen XRD")
        self.assertTrue(by_formula["Fe"].candidate_for_analysis)
        self.assertIn("disabled by user", by_formula["Cu"].skip_reasons)
        self.assertIn("outside selected element scope", by_formula["Cu"].skip_reasons)
        self.assertEqual(summary["total_structures"], 2)
        self.assertEqual(summary["used_structures"], 1)
        self.assertEqual(summary["skipped_structures"], 1)

    def test_library_analysis_exports_candidate_usage_summary_and_honors_disabled(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            library = self._library(root)
            cu_id = next(entry.internal_id for entry in library.list_structures() if entry.reduced_formula == "Cu")
            library.set_enabled(cu_id, False)
            xrd = root / "sample.xy"
            xrd.write_text("\n".join(f"{10 + i * 0.03:.3f} {1.0:.3f}" for i in range(1000)), encoding="utf-8")
            result = run_library_analysis(
                xrd_file=xrd,
                library_db=library.path,
                elements=["Fe", "Cu"],
                out_dir=root / "out",
            )
            usage_path = Path(result.outputs["candidate_usage_summary"])
            self.assertTrue(usage_path.exists())
            with usage_path.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            payload = json.loads(Path(result.outputs["json"]).read_text(encoding="utf-8"))
            audit_usage = Path(result.outputs["audit_folder"]) / "candidate_usage_summary.csv"
            self.assertTrue(audit_usage.exists())
            self.assertTrue((Path(result.outputs["audit_folder"]) / "xrd_presentation.png").exists())
            self.assertTrue((Path(result.outputs["audit_folder"]) / "xrd_diagnostic.png").exists())

        self.assertIn("candidate_usage", payload["input"])
        by_formula = {row["formula"]: row for row in rows}
        self.assertEqual(by_formula["Fe"]["used_in_analysis"], "True")
        self.assertEqual(by_formula["Cu"]["enabled"], "False")
        self.assertEqual(by_formula["Cu"]["used_in_analysis"], "False")

    def test_peak_inspector_returns_multiple_nearby_cached_candidates_with_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            library = self._library(root)
            settings_hash = simulation_settings_hash("CuKa", (10.0, 90.0))
            evidence = inspect_peak(library, 22.2, tolerance_deg=0.4, radiation="CuKa", settings_hash=settings_hash)

        self.assertEqual(evidence["status"], "matched")
        self.assertGreaterEqual(len(evidence["nearby_peaks"]), 2)
        first = evidence["nearby_peaks"][0]
        self.assertIn("structure_hash", first)
        self.assertIn("enabled", first)
        self.assertIn("cif_path", first)

    def test_peak_inspector_exports_json_and_csv(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            library = self._library(root)
            settings_hash = simulation_settings_hash("CuKa", (10.0, 90.0))
            evidence = inspect_peak(library, 22.2, tolerance_deg=0.4, radiation="CuKa", settings_hash=settings_hash)
            outputs = export_peak_inspection(evidence, root / "peak_inspection")

            json_payload = json.loads(Path(outputs["json"]).read_text(encoding="utf-8"))
            with Path(outputs["csv"]).open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))

        self.assertEqual(json_payload["status"], "matched")
        self.assertGreaterEqual(len(rows), 2)
        self.assertIn("structure_internal_id", rows[0])


if __name__ == "__main__":
    unittest.main()
