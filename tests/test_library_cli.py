import subprocess
import sys
import tempfile
import unittest
import os
import json
from pathlib import Path


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


class LibraryCliTests(unittest.TestCase):
    def test_library_import_list_and_cache_xrd(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cifs = root / "cifs"
            cifs.mkdir()
            (cifs / "Fe.cif").write_text(SIMPLE_CIF, encoding="utf-8")
            db = root / "library.sqlite"
            cwd = str(Path(__file__).parents[1])
            import_proc = subprocess.run(
                [sys.executable, "-m", "ophiuchus", "library-import", "--folder", str(cifs), "--library-db", str(db)],
                cwd=cwd,
                text=True,
                capture_output=True,
            )
            list_proc = subprocess.run(
                [sys.executable, "-m", "ophiuchus", "library-list", "--library-db", str(db)],
                cwd=cwd,
                text=True,
                capture_output=True,
            )
            cache_proc = subprocess.run(
                [sys.executable, "-m", "ophiuchus", "library-cache-xrd", "--library-db", str(db)],
                cwd=cwd,
                text=True,
                capture_output=True,
            )
        self.assertEqual(import_proc.returncode, 0, import_proc.stderr)
        self.assertIn('"imported": 1', import_proc.stdout)
        self.assertEqual(list_proc.returncode, 0, list_proc.stderr)
        self.assertIn("Fe", list_proc.stdout)
        self.assertEqual(cache_proc.returncode, 0, cache_proc.stderr)
        self.assertIn('"simulated": 1', cache_proc.stdout)

    def test_mp_test_without_key_returns_actionable_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            env = os.environ.copy()
            for key in ("MP_API_KEY", "PMG_MAPI_KEY", "MAPI_KEY"):
                env.pop(key, None)
            proc = subprocess.run(
                [sys.executable, "-m", "ophiuchus", "mp-test", "--env-path", str(Path(tmp) / ".env")],
                cwd=str(Path(__file__).parents[1]),
                env=env,
                text=True,
                capture_output=True,
            )
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("MP_API_KEY", proc.stderr)

    def test_xrd_validate_command_writes_reference_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = root / "reference.csv"
            reference.write_text("two_theta,intensity\n28.44,100\n47.30,60\n", encoding="utf-8")
            cif = root / "dummy.cif"
            cif.write_text("data_dummy\n", encoding="utf-8")
            out_dir = root / "validation"
            proc = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "ophiuchus",
                    "xrd-validate",
                    "--cif",
                    str(cif),
                    "--reference",
                    str(reference),
                    "--out-dir",
                    str(out_dir),
                ],
                cwd=str(Path(__file__).parents[1]),
                text=True,
                capture_output=True,
            )
            reports = list(out_dir.glob("*_validation.json"))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(len(reports), 1)

    def test_xrd_debug_simulate_and_compare_write_traceable_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cifs = root / "cifs"
            cifs.mkdir()
            cif = cifs / "Fe.cif"
            cif.write_text(SIMPLE_CIF, encoding="utf-8")
            peaks_csv = root / "debug_peaks.csv"
            ref_csv = root / "reference.csv"
            ref_csv.write_text("two_theta,intensity\n22.20,100\n", encoding="utf-8")
            compare_csv = root / "comparison.csv"
            cwd = str(Path(__file__).parents[1])
            sim_proc = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "ophiuchus",
                    "xrd-debug-simulate",
                    "--cif",
                    str(cif),
                    "--out",
                    str(peaks_csv),
                    "--two-theta-min",
                    "10",
                    "--two-theta-max",
                    "90",
                ],
                cwd=cwd,
                text=True,
                capture_output=True,
            )
            cmp_proc = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "ophiuchus",
                    "xrd-debug-compare",
                    "--ophi",
                    str(peaks_csv),
                    "--reference",
                    str(ref_csv),
                    "--out",
                    str(compare_csv),
                    "--position-tolerance",
                    "1.0",
                ],
                cwd=cwd,
                text=True,
                capture_output=True,
            )
            text = peaks_csv.read_text(encoding="utf-8")
            comparison = compare_csv.read_text(encoding="utf-8")
        self.assertEqual(sim_proc.returncode, 0, sim_proc.stderr)
        self.assertEqual(cmp_proc.returncode, 0, cmp_proc.stderr)
        self.assertIn("wavelength_angstrom", text)
        self.assertIn("d_spacing", text)
        self.assertIn("intensity_ratio_ophi_over_reference", comparison)

    def test_library_analyze_command_runs_from_cached_library(self):
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
            db = root / "library.sqlite"
            out = root / "out"
            exports = root / "exports"
            cwd = str(Path(__file__).parents[1])
            subprocess.run(
                [sys.executable, "-m", "ophiuchus", "library-import", "--folder", str(cifs), "--library-db", str(db)],
                cwd=cwd,
                text=True,
                capture_output=True,
                check=True,
            )
            proc = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "ophiuchus",
                    "library-analyze",
                    "--xrd",
                    str(xrd),
                    "--library-db",
                    str(db),
                    "--elements",
                    "Fe",
                    "--out-dir",
                    str(out),
                    "--project",
                    "FeProject",
                    "--sample-id",
                    "FeSample",
                    "--export-folder",
                    str(exports),
                    "--export-date",
                    "20260629",
                ],
                cwd=cwd,
                text=True,
                capture_output=True,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue((out / "results.json").exists())
            payload = json.loads((out / "results.json").read_text(encoding="utf-8"))
            self.assertTrue(any(path.name.endswith("_analysis_report_v01.json") for path in exports.iterdir()))
            self.assertTrue(any("xrd_candidate_screening_presentation" in path.name for path in exports.iterdir()))
            self.assertTrue(any("xrd_candidate_screening_diagnostic" in path.name for path in exports.iterdir()))
            self.assertIn("candidate_source_summary", payload["input"])
            self.assertEqual(payload["input"]["candidate_source_summary"]["local"]["used"], 1)
            self.assertIn("library analysis finished", proc.stdout)


if __name__ == "__main__":
    unittest.main()
