import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np

from ophiuchus.phase_stripping.export import export_phase_stripping_session
from ophiuchus.phase_stripping.models import AnalysisContext, FitBounds
from ophiuchus.phase_stripping.session import PhaseStrippingSession
from ophiuchus.xrd.models import Candidate, Peak


def make_session() -> PhaseStrippingSession:
    x = np.arange(20.0, 22.001, 0.002)
    candidate = Candidate(
        candidate_id="phase-fege",
        formula_pretty="FeGe",
        source="test",
        source_path="library/FeGe.cif",
        elements=["Fe", "Ge"],
        structure_hash="structure-sha256",
        theory_peaks=[Peak("p1", 21.0, 100.0)],
        simulation_validation={"pattern_fingerprint": "pattern-sha256"},
    )
    context = AnalysisContext(
        x=x,
        intensity=1.5 * np.exp(-0.5 * ((x - 21.0) / 0.05) ** 2),
        radiation="CuKalpha12",
        wavelength_angstrom=1.54056,
        two_theta_range=(20.0, 22.0),
        tolerance_deg=0.15,
        source_path="sample.xy",
        source_fingerprint="source-sha256",
        data_fingerprint="data-sha256",
    )
    session = PhaseStrippingSession(
        context,
        bounds=FitBounds(
            shift_deg=(-0.001, 0.001),
            sigma_deg=(0.049, 0.051),
            scale=(0.0, 2.0),
        ),
    )
    session.accept_preview(session.preview(candidate))
    return session


class PhaseStrippingExportTests(unittest.TestCase):
    def test_export_writes_signed_residual_csv_reproducible_json_and_plot_without_simulation(self):
        session = make_session()

        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp) / "phase-export"
            outputs = export_phase_stripping_session(session, folder)

            with Path(outputs["csv"]).open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            payload = json.loads(Path(outputs["json"]).read_text(encoding="utf-8"))

            self.assertEqual(
                list(rows[0]),
                ["two_theta", "original_y", "fitted_total", "residual_y", "contribution_operation_1"],
            )
            self.assertEqual(len(rows), session.context.x.size)
            self.assertAlmostEqual(float(rows[0]["original_y"]), float(session.context.intensity[0]))
            self.assertAlmostEqual(float(rows[0]["fitted_total"]), float(session.fitted_total[0]))
            self.assertAlmostEqual(float(rows[0]["residual_y"]), float(session.residual_y[0]))
            self.assertEqual(payload["accepted_operations"][0]["operation_id"], "operation-1")
            self.assertEqual(payload["accepted_operations"][0]["candidate"]["candidate_id"], "phase-fege")
            self.assertEqual(payload["accepted_operations"][0]["candidate"]["cif_path"], "library/FeGe.cif")
            self.assertEqual(payload["accepted_operations"][0]["candidate"]["structure_hash"], "structure-sha256")
            self.assertEqual(payload["accepted_operations"][0]["candidate"]["pattern_fingerprint"], "pattern-sha256")
            self.assertEqual(payload["instrument_settings"]["radiation"], "CuKalpha12")
            self.assertEqual(payload["fit_bounds"]["sigma_deg"], [0.049, 0.051])
            self.assertEqual(payload["accepted_operations"][0]["fit"]["scale"], session.accepted_operations[0].phase_fit.scale)
            restored = PhaseStrippingSession.from_dict(payload["session"])
            self.assertEqual(restored.accepted_operations, session.accepted_operations)
            np.testing.assert_array_equal(restored.context.x, session.context.x)
            np.testing.assert_array_equal(restored.residual_y, session.residual_y)
            self.assertTrue(Path(outputs["png"]).is_file())
            self.assertGreater(Path(outputs["png"]).stat().st_size, 0)

    def test_export_does_not_invoke_a_simulator(self):
        session = make_session()

        with tempfile.TemporaryDirectory() as tmp:
            import ophiuchus.xrd.pymatgen_simulator as simulator

            original = simulator.simulate_cif_with_pymatgen

            def fail_if_called(*args, **kwargs):
                raise AssertionError("export must use accepted session contributions, not simulate")

            simulator.simulate_cif_with_pymatgen = fail_if_called
            try:
                export_phase_stripping_session(session, Path(tmp))
            finally:
                simulator.simulate_cif_with_pymatgen = original

    def test_export_failure_leaves_no_partial_bundle(self):
        session = make_session()

        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp) / "phase-export"
            with mock.patch("ophiuchus.phase_stripping.export._write_plot", side_effect=OSError("plot failed")):
                with self.assertRaises(OSError):
                    export_phase_stripping_session(session, folder)

            self.assertFalse(folder.exists())
            self.assertTrue(Path(export_phase_stripping_session(session, folder)["png"]).is_file())

    def test_export_rejects_legacy_operation_without_reproducible_candidate_identity(self):
        session = make_session()
        payload = session.to_dict()
        payload["accepted"][0]["candidate"] = None
        legacy = PhaseStrippingSession.from_dict(payload)

        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "candidate provenance"):
                export_phase_stripping_session(legacy, Path(tmp) / "phase-export")


if __name__ == "__main__":
    unittest.main()
