import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from ophiuchus.xrd.config import XRDConfig
from ophiuchus.xrd.models import Peak
from ophiuchus.xrd.validation import (
    PatternComparator,
    ReferencePatternImporter,
    ValidationRunner,
    simulate_cif_with_config,
)


class XrdValidationTests(unittest.TestCase):
    def test_xrd_config_defaults_to_cu_kalpha12_and_degrees(self):
        config = XRDConfig.validation_default()
        self.assertEqual(config.radiation_source, "CuKalpha12")
        self.assertAlmostEqual(config.wavelength_angstrom, 1.54056, places=5)
        self.assertEqual(config.line_model, "cu_kalpha12")
        self.assertEqual(config.angle_unit, "degree_2theta")
        self.assertEqual(config.debye_waller_b, 0.0)
        self.assertIn("CuKalpha12", config.config_id())

    def test_reference_pattern_importer_reads_vesta_like_peak_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "vesta_peaks.txt"
            path.write_text(
                "2-Theta    d(?)    I(f)    ( h k l )\n"
                "28.44      3.136   100.0   ( 1 1 1 )\n"
                "47.30      1.920   60.0    ( 2 2 0 )\n",
                encoding="utf-8",
            )
            pattern = ReferencePatternImporter().read(path, label="VESTA Si")
        self.assertEqual(pattern.label, "VESTA Si")
        self.assertEqual(len(pattern.peaks), 2)
        self.assertEqual(pattern.peaks[0].hkl, "1 1 1")
        self.assertAlmostEqual(pattern.peaks[0].two_theta, 28.44)

    def test_pattern_comparator_reports_position_missing_and_extra(self):
        reference = [Peak("r1", 28.44, 100.0), Peak("r2", 47.30, 60.0), Peak("r3", 56.12, 10.0)]
        calculated = [Peak("c1", 28.45, 95.0), Peak("c2", 47.38, 50.0), Peak("c3", 69.0, 30.0)]
        report = PatternComparator(tolerance_deg=0.10, strong_intensity_threshold=20.0).compare(calculated, reference)
        self.assertEqual(report.summary["matched_count"], 2)
        self.assertEqual(report.summary["missing_strong_count"], 0)
        self.assertEqual(report.summary["extra_strong_count"], 1)
        self.assertLessEqual(report.summary["median_abs_delta_2theta"], 0.08)

    def test_pattern_comparator_reports_intensity_pass_fail(self):
        reference = [Peak("r1", 35.36, 70.0), Peak("r2", 44.78, 100.0), Peak("r3", 74.80, 22.0)]
        calculated = [Peak("c1", 35.36, 70.0), Peak("c2", 44.78, 100.0), Peak("c3", 74.80, 88.0)]
        report = PatternComparator(tolerance_deg=0.03, strong_intensity_threshold=5.0).compare(calculated, reference)
        self.assertFalse(report.summary["intensity_pass"])
        self.assertGreater(report.summary["max_strong_intensity_ratio_error"], 2.0)
        near_75 = [item for item in report.matched if round(float(item["reference_two_theta"]), 2) == 74.8][0]
        self.assertAlmostEqual(near_75["intensity_ratio_calculated_over_reference"], 4.0)

    def test_validation_runner_writes_json_and_markdown_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = root / "reference.csv"
            reference.write_text("two_theta,intensity\n28.44,100\n47.30,60\n", encoding="utf-8")

            def simulator(_cif, _config):
                return [Peak("c1", 28.45, 100.0), Peak("c2", 47.31, 50.0)]

            result = ValidationRunner(simulator=simulator).run(
                cif_file=root / "dummy.cif",
                reference_file=reference,
                out_dir=root / "validation_reports",
                config=XRDConfig.validation_default(),
                label="synthetic_vesta",
            )
            json_path = Path(result["json_report"])
            md_path = Path(result["markdown_report"])
            self.assertTrue(json_path.exists())
            self.assertTrue(md_path.exists())
            payload = json.loads(json_path.read_text(encoding="utf-8"))
        self.assertEqual(payload["summary"]["matched_count"], 2)
        self.assertIn("median_abs_delta_2theta", payload["summary"])

    def test_validation_simulation_reports_malformed_cif_from_canonical_backend(self):
        with tempfile.TemporaryDirectory() as tmp:
            cif = Path(tmp) / "dummy.cif"
            cif.write_text("data_dummy\n", encoding="utf-8")
            config = XRDConfig.validation_default()
            with patch.dict(
                "os.environ",
                {
                    "OPHI_DISABLE_LOCAL_VESTA_REFERENCE": "1",
                    "OPHI_DISABLE_EXTERNAL_PYMATGEN": "1",
                    "OPHI_DISABLE_INPROCESS_PYMATGEN": "1",
                },
                clear=False,
            ):
                with self.assertRaisesRegex(RuntimeError, "Canonical XRD backend could not parse CIF"):
                    simulate_cif_with_config(cif, config)

    def test_reference_importer_can_filter_existing_vesta_phase_csv(self):
        configured = os.environ.get("OPHI_TEST_VESTA_PEAK_CSV", "")
        if not configured:
            self.skipTest("OPHI_TEST_VESTA_PEAK_CSV is not configured")
        reference = Path(configured)
        if not reference.exists():
            self.skipTest("configured VESTA reference CSV is not available")
        pattern = ReferencePatternImporter().read(reference, label="ZrFe6Ge4 VESTA", phase="VESTA simulated")
        self.assertGreaterEqual(len(pattern.peaks), 8)
        self.assertTrue(all(30.0 < peak.two_theta < 80.0 for peak in pattern.peaks))

    def test_reference_importer_extracts_peaks_from_existing_vesta_int_pattern(self):
        configured = os.environ.get("OPHI_TEST_VESTA_MONI_INT", "")
        if not configured:
            self.skipTest("OPHI_TEST_VESTA_MONI_INT is not configured")
        reference = Path(configured)
        if not reference.exists():
            self.skipTest("configured VESTA MONI.int reference is not available")
        pattern = ReferencePatternImporter().read(reference, label="ZrFe6Ge4 MONI")
        near_70 = [peak for peak in pattern.peaks if abs(peak.two_theta - 70.07) <= 0.05]
        self.assertTrue(near_70)
        self.assertGreater(near_70[0].intensity, 5.0)


if __name__ == "__main__":
    unittest.main()
