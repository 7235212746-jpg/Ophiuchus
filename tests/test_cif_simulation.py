import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from ophiuchus.xrd.cache import XRD_SIMULATOR_VERSION
from ophiuchus.xrd.candidates import LocalCandidateProvider, simulate_or_load_peaks
from ophiuchus.xrd.config import XRDConfig
from ophiuchus.xrd.pipeline import run_analysis


SIMPLE_CUBIC_CIF = """
data_test
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


class CifSimulationTests(unittest.TestCase):
    def test_local_candidate_cache_version_marks_raw_pymatgen_b0_fix(self):
        self.assertIn("raw_pymatgen", XRD_SIMULATOR_VERSION)
        self.assertIn("kalpha12", XRD_SIMULATOR_VERSION)
        self.assertIn("b0", XRD_SIMULATOR_VERSION)

    def test_default_line_model_adds_cu_kalpha2_companion_peaks(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "Fe.cif"
            path.write_text(SIMPLE_CUBIC_CIF, encoding="utf-8")
            candidate = LocalCandidateProvider([Path(tmp)], allowed_elements={"Fe"}).iter_candidates()[0]
            simulate_or_load_peaks(candidate, two_theta_range=(10, 90))
        self.assertIn("Kalpha1", candidate.simulated_pattern.line_component)
        self.assertIn("Kalpha2", candidate.simulated_pattern.line_component)

    def test_canonical_cif_simulation_produces_expected_cubic_peak(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "Fe.cif"
            path.write_text(SIMPLE_CUBIC_CIF, encoding="utf-8")
            candidate = LocalCandidateProvider([Path(tmp)], allowed_elements={"Fe"}).iter_candidates()[0]
            with patch.dict(
                "os.environ",
                {
                    "OPHI_DISABLE_EXTERNAL_PYMATGEN": "1",
                    "OPHI_DISABLE_INPROCESS_PYMATGEN": "1",
                    "OPHI_ALLOW_UNVALIDATED_FALLBACK_XRD": "1",
                },
                clear=False,
            ):
                peaks = simulate_or_load_peaks(candidate, two_theta_range=(10, 90))
        self.assertGreater(len(peaks), 3)
        self.assertTrue(any(abs(p.two_theta - 22.2) < 0.5 for p in peaks), [p.two_theta for p in peaks[:8]])
        self.assertAlmostEqual(max(p.intensity for p in peaks), 100.0)

    def test_cif_simulation_ignores_legacy_disable_flags_and_uses_canonical_backend(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "Fe.cif"
            path.write_text(SIMPLE_CUBIC_CIF, encoding="utf-8")
            candidate = LocalCandidateProvider([Path(tmp)], allowed_elements={"Fe"}).iter_candidates()[0]
            with patch.dict("os.environ", {"OPHI_DISABLE_EXTERNAL_PYMATGEN": "1", "OPHI_DISABLE_INPROCESS_PYMATGEN": "1"}, clear=False):
                peaks = simulate_or_load_peaks(candidate, two_theta_range=(10, 90))
        self.assertTrue(peaks)
        self.assertEqual(candidate.parse_status, "validated_backend_simulated")

    def test_inprocess_pymatgen_preserves_raw_pattern_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "Fe.cif"
            path.write_text(SIMPLE_CUBIC_CIF, encoding="utf-8")
            candidate = LocalCandidateProvider([Path(tmp)], allowed_elements={"Fe"}).iter_candidates()[0]
            with patch.dict("os.environ", {"OPHI_DISABLE_EXTERNAL_PYMATGEN": "1"}, clear=False):
                peaks = simulate_or_load_peaks(candidate, two_theta_range=(10, 90))
        self.assertEqual(candidate.parse_status, "validated_backend_simulated")
        self.assertTrue(all(peak.hkl for peak in peaks[:3]))
        self.assertTrue(all(peak.d_spacing for peak in peaks[:3]))

    def test_pipeline_handles_cif_only_candidate_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            xrd = root / "sample.xy"
            rows = []
            for i in range(2000):
                x = 10.0 + i * 0.03
                y = 1.0 + 100.0 * pow(2.718281828, -((x - 22.2) ** 2) / (2 * 0.06**2))
                rows.append(f"{x:.3f} {y:.5f}")
            xrd.write_text("\n".join(rows), encoding="utf-8")
            cifs = root / "cifs"
            cifs.mkdir()
            (cifs / "Fe.cif").write_text(SIMPLE_CUBIC_CIF, encoding="utf-8")
            out = root / "out"
            result = run_analysis(xrd, [cifs], ["Fe"], out_dir=out)
        self.assertTrue(result.top_scores)
        self.assertEqual(result.top_scores[0].candidate.formula_pretty, "Fe")


if __name__ == "__main__":
    unittest.main()
