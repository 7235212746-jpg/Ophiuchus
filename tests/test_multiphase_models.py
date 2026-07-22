from pathlib import Path
import tempfile
import unittest

import numpy as np

from ophiuchus.xrd.multiphase_models import (
    MultiphaseRefinementResult,
    MultiphaseRefinementSettings,
    PhaseRefinementInput,
    PhaseRefinementResult,
)
from ophiuchus.xrd.quantification import weight_fractions_from_zmv


class MultiphaseModelTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.cif = self.root / "phase.cif"
        self.cif.write_text("data_phase\n", encoding="ascii")

    def tearDown(self) -> None:
        self.temp.cleanup()

    def test_weight_fractions_use_scale_times_zmv(self):
        values = weight_fractions_from_zmv(
            [
                ("main", 2.0, 4.0, 100.0, 50.0),
                ("imp", 1.0, 2.0, 50.0, 25.0),
            ]
        )

        self.assertAlmostEqual(values["main"], 94.117647, places=6)
        self.assertAlmostEqual(values["imp"], 5.882353, places=6)
        self.assertAlmostEqual(sum(values.values()), 100.0, places=10)

    def test_weight_fractions_reject_non_positive_scientific_inputs(self):
        with self.assertRaisesRegex(ValueError, "positive finite"):
            weight_fractions_from_zmv([("main", -1.0, 4.0, 100.0, 50.0)])
        with self.assertRaisesRegex(ValueError, "positive finite"):
            weight_fractions_from_zmv([("main", 1.0, 0.0, 100.0, 50.0)])

    def test_phase_input_requires_a_real_cif_and_valid_role(self):
        with self.assertRaises(FileNotFoundError):
            PhaseRefinementInput("x", "X", self.root / "missing.cif", "target")

        text = self.root / "phase.txt"
        text.write_text("data", encoding="ascii")
        with self.assertRaisesRegex(ValueError, "CIF"):
            PhaseRefinementInput("x", "X", text, "target")
        with self.assertRaisesRegex(ValueError, "role"):
            PhaseRefinementInput("x", "X", self.cif, "unknown")

    def test_refinement_settings_enforce_first_release_limits(self):
        with self.assertRaisesRegex(ValueError, "increasing"):
            MultiphaseRefinementSettings(90.0, 10.0)
        with self.assertRaisesRegex(ValueError, "Cu Kalpha"):
            MultiphaseRefinementSettings(10.0, 90.0, radiation="MoKa")
        with self.assertRaisesRegex(ValueError, "background"):
            MultiphaseRefinementSettings(10.0, 90.0, background_terms=13)

    def test_refinement_result_arrays_are_immutable(self):
        phase = PhaseRefinementResult(
            phase_id="main",
            formula="Main",
            scale=1.0,
            z=1.0,
            molar_mass=100.0,
            volume_angstrom3=50.0,
            weight_percent=100.0,
            contribution_intensity=np.array([2.0, 3.0, 4.0]),
        )
        result = MultiphaseRefinementResult(
            two_theta_deg=np.array([10.0, 10.1, 10.2]),
            observed_intensity=np.array([12.0, 20.0, 13.0]),
            calculated_intensity=np.array([11.0, 18.0, 13.5]),
            residual_intensity=np.array([1.0, 2.0, -0.5]),
            background_intensity=np.array([3.0, 3.1, 3.2]),
            phases=(phase,),
            rwp_percent=8.5,
            rp_percent=6.0,
            goodness_of_fit=1.1,
        )

        self.assertFalse(result.observed_intensity.flags.writeable)
        self.assertFalse(result.phases[0].contribution_intensity.flags.writeable)
        with self.assertRaises(ValueError):
            result.observed_intensity[0] = 0.0


if __name__ == "__main__":
    unittest.main()
