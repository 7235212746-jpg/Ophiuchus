import unittest

import numpy as np

from ophiuchus.phase_stripping.composition import inspect_peak_composition
from ophiuchus.phase_stripping.models import AnalysisContext
from ophiuchus.phase_stripping.session import PhaseStrippingSession
from ophiuchus.xrd.models import Candidate, Peak


def make_session():
    x = np.arange(20.0, 24.001, 0.005)
    first = Candidate(
        "phase-a",
        "FeGe",
        "library:local",
        "FeGe.cif",
        ["Fe", "Ge"],
        "hash-a",
        theory_peaks=[Peak("a1", 21.0, 100.0, hkl="100")],
    )
    second = Candidate(
        "phase-b",
        "Fe3Ge",
        "library:local",
        "Fe3Ge.cif",
        ["Fe", "Ge"],
        "hash-b",
        theory_peaks=[Peak("b1", 21.04, 100.0, hkl="110")],
    )
    intensity = 4.0 * np.exp(-0.5 * ((x - 21.0) / 0.06) ** 2) + 0.5
    context = AnalysisContext(
        x=x,
        intensity=intensity,
        radiation="CuKalpha12",
        wavelength_angstrom=1.54056,
        two_theta_range=(20.0, 24.0),
        tolerance_deg=0.15,
        source_path="sample.xy",
        source_fingerprint="source",
        data_fingerprint="data",
    )
    session = PhaseStrippingSession(
        context,
        background_y=np.full_like(x, 0.5),
        background_method="asls",
        background_parameters={"smoothness": 1.0e6, "asymmetry": 0.001, "iterations": 12},
    )
    accepted = session.preview_with_parameters(first, scale=3.0, shift_deg=0.0, sigma_deg=0.06)
    session.accept_preview(accepted)
    preview = session.preview_with_parameters(second, scale=0.5, shift_deg=-0.04, sigma_deg=0.06)
    return session, {first.candidate_id: first, second.candidate_id: second}, preview


class PeakCompositionTests(unittest.TestCase):
    def test_snaps_to_grid_and_preserves_reconstruction_identity(self):
        session, candidates, _preview = make_session()

        result = inspect_peak_composition(session, candidates, 21.002)

        self.assertAlmostEqual(result.two_theta, 21.0, places=9)
        self.assertAlmostEqual(
            result.experimental,
            result.background + result.explained + result.residual,
            places=12,
        )
        self.assertAlmostEqual(result.corrected, result.explained + result.residual, places=12)
        self.assertEqual(
            [row.kind for row in result.rows[:5]],
            ["experimental", "background", "corrected", "explained", "residual"],
        )

    def test_phase_row_uses_formula_share_and_nearest_shifted_reflection(self):
        session, candidates, _preview = make_session()

        result = inspect_peak_composition(session, candidates, 21.0)
        phase = next(row for row in result.rows if row.kind == "accepted_phase")

        self.assertEqual(phase.label, "FeGe")
        self.assertAlmostEqual(phase.intensity, 3.0, places=6)
        self.assertAlmostEqual(phase.explained_share_percent, 100.0, places=6)
        self.assertEqual(phase.hkl, "100")
        self.assertAlmostEqual(phase.reflection_two_theta, 21.0, places=9)
        self.assertAlmostEqual(phase.reflection_delta, 0.0, places=9)

    def test_preview_is_reported_separately_and_not_added_to_accepted_identity(self):
        session, candidates, preview = make_session()

        result = inspect_peak_composition(session, candidates, 21.0, preview=preview)
        preview_row = next(row for row in result.rows if row.kind == "preview_phase")

        self.assertEqual(preview_row.label, "Fe3Ge (预览)")
        self.assertAlmostEqual(preview_row.intensity, 0.5, places=6)
        self.assertAlmostEqual(
            result.experimental,
            result.background + result.explained + result.residual,
            places=12,
        )
        self.assertAlmostEqual(result.explained, 3.0, places=6)
        self.assertEqual(preview_row.hkl, "110")


if __name__ == "__main__":
    unittest.main()
