from dataclasses import FrozenInstanceError
import json
import unittest

import numpy as np

from ophiuchus.phase_stripping.models import (
    AnalysisContext,
    FitBounds,
    InstrumentSettings,
    PhaseFit,
    PhaseOperation,
)
from ophiuchus.phase_stripping.profile import project_candidate_profile
from ophiuchus.phase_stripping.fitting import PhaseContributionFitter
from ophiuchus.phase_stripping.session import PhaseStrippingSession
from ophiuchus.xrd.models import Candidate, Peak


def make_candidate(peaks: list[Peak] | None = None) -> Candidate:
    return Candidate(
        candidate_id="candidate-1",
        formula_pretty="FeGe",
        source="test",
        source_path="candidate.xy",
        elements=["Fe", "Ge"],
        structure_hash="abc123",
        theory_peaks=peaks or [],
    )


def make_context(x: np.ndarray, intensity: np.ndarray) -> AnalysisContext:
    return AnalysisContext(
        x=x,
        intensity=intensity,
        radiation="CuKalpha12",
        wavelength_angstrom=1.54056,
        two_theta_range=(float(x[0]), float(x[-1])),
        tolerance_deg=0.15,
        source_path="sample.xy",
        source_fingerprint="source-sha256",
        data_fingerprint="data-sha256",
    )


class PhaseStrippingSessionTests(unittest.TestCase):
    def test_analysis_context_preserves_unnormalized_immutable_experimental_data(self):
        x = np.array([10.0, 10.1, 10.2])
        intensity = np.array([12.0, 40.0, 18.0])

        context = AnalysisContext(
            x=x,
            intensity=intensity,
            radiation="CuKalpha12",
            wavelength_angstrom=1.54056,
            two_theta_range=(10.0, 10.2),
            tolerance_deg=0.15,
            source_path="sample.xy",
            source_fingerprint="source-sha256",
            data_fingerprint="data-sha256",
        )

        np.testing.assert_array_equal(context.x, x)
        np.testing.assert_array_equal(context.intensity, intensity)
        self.assertEqual(context.intensity.max(), 40.0)
        self.assertFalse(context.x.flags.writeable)
        self.assertFalse(context.intensity.flags.writeable)
        with self.assertRaises(FrozenInstanceError):
            context.source_path = "other.xy"
        with self.assertRaises(ValueError):
            context.intensity[0] = 0.0

    def test_analysis_context_preserves_input_dtypes_endianness_and_bytes(self):
        x = np.array([10.0, 10.1, 10.2], dtype=">f4")
        intensity = np.array([12, 40, 18], dtype=">i2")

        context = make_context(x, intensity)

        self.assertEqual(context.x.dtype.str, x.dtype.str)
        self.assertEqual(context.intensity.dtype.str, intensity.dtype.str)
        self.assertEqual(context.x.tobytes(), x.tobytes())
        self.assertEqual(context.intensity.tobytes(), intensity.tobytes())

    def test_session_data_objects_are_frozen(self):
        bounds = FitBounds(
            shift_deg=(-0.2, 0.2),
            sigma_deg=(0.02, 0.20),
            scale=(0.0, 1.0),
        )
        settings = InstrumentSettings(
            radiation="CuKalpha12",
            wavelength_angstrom=1.54056,
            two_theta_range=(10.0, 90.0),
            tolerance_deg=0.15,
        )
        fit = PhaseFit(candidate_id="candidate-1", shift_deg=0.02, sigma_deg=0.05, scale=0.8)
        operation = PhaseOperation(operation_id="operation-1", phase_fit=fit, residual_fingerprint="residual-sha256")

        for instance, attribute, value in (
            (bounds, "scale", (0.0, 2.0)),
            (settings, "radiation", "CoKa"),
            (fit, "scale", 0.5),
            (operation, "operation_id", "operation-2"),
        ):
            with self.assertRaises(FrozenInstanceError):
                setattr(instance, attribute, value)

    def test_profile_uses_all_peaks_with_uniform_shift_and_normalizes_after_summing(self):
        x = np.arange(0.0, 8.01, 0.01)
        candidate = make_candidate([Peak("p1", 2.0, 100.0), Peak("p2", 5.0, 50.0)])
        original_positions = [peak.two_theta for peak in candidate.theory_peaks]
        original_intensities = [peak.intensity for peak in candidate.theory_peaks]

        profile = project_candidate_profile(x, candidate, shift_deg=0.20, sigma_deg=0.03)

        self.assertEqual(profile.shape, x.shape)
        self.assertAlmostEqual(float(profile.max()), 1.0)
        self.assertAlmostEqual(profile[np.argmin(abs(x - 5.20))], 0.5, places=12)
        self.assertAlmostEqual(x[int(profile.argmax())], 2.20, places=2)
        self.assertEqual([peak.two_theta for peak in candidate.theory_peaks], original_positions)
        self.assertEqual([peak.intensity for peak in candidate.theory_peaks], original_intensities)

    def test_profile_normalizes_overlapping_peak_sum_once(self):
        x = np.arange(1.70, 2.31, 0.001)
        sigma_deg = 0.05
        candidate = make_candidate([Peak("p1", 2.0, 100.0), Peak("p2", 2.04, 40.0)])

        profile = project_candidate_profile(x, candidate, shift_deg=0.0, sigma_deg=sigma_deg)

        expected_sum = 100.0 * np.exp(-0.5 * ((x - 2.0) / sigma_deg) ** 2)
        expected_sum += 40.0 * np.exp(-0.5 * ((x - 2.04) / sigma_deg) ** 2)
        expected_profile = expected_sum / expected_sum.max()
        np.testing.assert_allclose(profile, expected_profile, rtol=0.0, atol=1e-12)

    def test_profile_rejects_invalid_inputs(self):
        candidate = make_candidate([Peak("p1", 2.0, 100.0)])

        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            project_candidate_profile(np.array([1.0, 1.0, 2.0]), candidate, 0.0, 0.05)
        with self.assertRaisesRegex(ValueError, "positive"):
            project_candidate_profile(np.array([1.0, 2.0]), candidate, 0.0, 0.0)
        with self.assertRaisesRegex(ValueError, "peaks"):
            project_candidate_profile(np.array([1.0, 2.0]), make_candidate(), 0.0, 0.05)

    def test_fitter_recovers_single_phase_and_session_leaves_near_zero_residual(self):
        x = np.arange(20.0, 25.001, 0.002)
        candidate = make_candidate([Peak("p1", 21.0, 100.0), Peak("p2", 23.2, 55.0)])
        expected = 1.7 * project_candidate_profile(x, candidate, shift_deg=0.08, sigma_deg=0.06)
        bounds = FitBounds(shift_deg=(-0.25, 0.25), sigma_deg=(0.035, 0.25), scale=(0.0, 3.0))

        session = PhaseStrippingSession(make_context(x, expected), bounds=bounds)
        preview = session.preview(candidate)
        session.accept_preview(preview)

        self.assertAlmostEqual(preview.phase_fit.shift_deg, 0.08, places=4)
        self.assertAlmostEqual(preview.phase_fit.sigma_deg, 0.06, places=4)
        self.assertAlmostEqual(preview.phase_fit.scale, 1.7, places=4)
        np.testing.assert_allclose(session.residual_y, np.zeros_like(expected), rtol=0.0, atol=1e-5)

    def test_session_subtracts_fixed_background_before_phase_fitting(self):
        x = np.arange(20.0, 24.001, 0.002)
        background = 30.0 + 0.5 * (x - 20.0)
        candidate = make_candidate([Peak("p1", 21.0, 100.0), Peak("p2", 23.0, 50.0)])
        phase = 2.0 * project_candidate_profile(x, candidate, shift_deg=0.0, sigma_deg=0.05)
        session = PhaseStrippingSession(
            make_context(x, background + phase),
            background_y=background,
            background_method="asls",
            background_parameters={"smoothness": 1.0e6, "asymmetry": 0.001, "iterations": 12},
        )

        np.testing.assert_allclose(session.corrected_intensity, phase, rtol=0.0, atol=1e-12)
        np.testing.assert_allclose(session.residual_y, phase, rtol=0.0, atol=1e-12)
        session.accept_preview(session.preview(candidate))

        np.testing.assert_allclose(session.residual_y, np.zeros_like(phase), rtol=0.0, atol=1e-5)
        np.testing.assert_allclose(session.reconstructed_y, background + session.fitted_total)
        self.assertFalse(session.background_y.flags.writeable)

    def test_session_serialization_preserves_background_model_exactly(self):
        x = np.array([20.0, 20.1, 20.2])
        background = np.array([10.0, 11.0, 12.0])
        session = PhaseStrippingSession(
            make_context(x, background + np.array([0.0, 3.0, 0.0])),
            background_y=background,
            background_method="asls",
            background_parameters={"smoothness": 123.0, "asymmetry": 0.01, "iterations": 8},
        )

        restored = PhaseStrippingSession.from_dict(json.loads(json.dumps(session.to_dict(), allow_nan=False)))

        np.testing.assert_array_equal(restored.background_y, background)
        np.testing.assert_array_equal(restored.corrected_intensity, session.corrected_intensity)
        self.assertEqual(restored.background_method, "asls")
        self.assertEqual(restored.background_parameters, session.background_parameters)

    def test_fitter_recovers_one_bounded_shift_for_every_candidate_peak(self):
        x = np.arange(20.0, 28.001, 0.002)
        candidate = make_candidate([Peak("p1", 21.0, 100.0), Peak("p2", 24.0, 60.0), Peak("p3", 26.5, 45.0)])
        expected_shift = 0.17
        y = 1.2 * project_candidate_profile(x, candidate, shift_deg=expected_shift, sigma_deg=0.07)
        bounds = FitBounds(shift_deg=(-0.25, 0.25), sigma_deg=(0.035, 0.25), scale=(0.0, 2.0))

        fit = PhaseContributionFitter().fit(y, x, candidate, bounds)

        self.assertAlmostEqual(fit.shift_deg, expected_shift, places=4)
        self.assertGreaterEqual(fit.shift_deg, bounds.shift_deg[0])
        self.assertLessEqual(fit.shift_deg, bounds.shift_deg[1])

    def test_fit_mask_excludes_outlier_points_from_the_fit(self):
        x = np.arange(20.0, 22.001, 0.002)
        candidate = make_candidate([Peak("p1", 21.0, 100.0)])
        profile = project_candidate_profile(x, candidate, shift_deg=0.0, sigma_deg=0.05)
        residual = profile.copy()
        residual[np.argmin(abs(x - 21.18))] += 100.0
        bounds = FitBounds(shift_deg=(-0.2, 0.2), sigma_deg=(0.049, 0.051), scale=(0.0, 2.0))
        fit_mask = x < 21.08

        fit = PhaseContributionFitter().fit(residual, x, candidate, bounds, fit_mask=fit_mask)

        self.assertAlmostEqual(fit.shift_deg, 0.0, places=3)
        self.assertAlmostEqual(fit.scale, 1.0, places=3)

    def test_fitter_rejects_shape_mismatched_or_empty_fit_masks(self):
        x = np.arange(20.0, 22.001, 0.002)
        candidate = make_candidate([Peak("p1", 21.0, 100.0)])
        bounds = FitBounds(shift_deg=(-0.2, 0.2), sigma_deg=(0.049, 0.051), scale=(0.0, 2.0))

        with self.assertRaisesRegex(ValueError, "same shape"):
            PhaseContributionFitter().fit(np.zeros_like(x), x, candidate, bounds, fit_mask=np.ones(x.size - 1))
        with self.assertRaisesRegex(ValueError, "at least one point"):
            PhaseContributionFitter().fit(np.zeros_like(x), x, candidate, bounds, fit_mask=np.zeros_like(x, dtype=bool))

    def test_removing_phase_a_from_an_overlapping_a_plus_b_pattern_leaves_phase_b(self):
        x = np.arange(20.0, 24.001, 0.002)
        candidate_a = make_candidate([Peak("a", 21.0, 100.0)])
        candidate_a.candidate_id = "phase-a"
        candidate_b = make_candidate([Peak("b", 21.18, 100.0)])
        candidate_b.candidate_id = "phase-b"
        profile_a = project_candidate_profile(x, candidate_a, shift_deg=0.0, sigma_deg=0.05)
        profile_b = project_candidate_profile(x, candidate_b, shift_deg=0.0, sigma_deg=0.05)
        original = 1.0 * profile_a + 0.8 * profile_b
        bounds = FitBounds(shift_deg=(-0.01, 0.01), sigma_deg=(0.049, 0.051), scale=(0.0, 2.0))
        session = PhaseStrippingSession(make_context(x, original), bounds=bounds)

        session.accept_preview(session.preview(candidate_a))

        b_center = int(np.argmin(abs(x - 21.18)))
        self.assertGreater(session.residual_y[b_center], 0.45)
        self.assertGreater(session.residual_y[b_center], session.residual_y[int(np.argmin(abs(x - 21.0)))])

    def test_over_subtraction_keeps_negative_values_and_warns(self):
        x = np.arange(20.0, 22.001, 0.002)
        candidate = make_candidate([Peak("p1", 21.0, 100.0)])
        profile = project_candidate_profile(x, candidate, shift_deg=0.0, sigma_deg=0.05)
        bounds = FitBounds(shift_deg=(-0.001, 0.001), sigma_deg=(0.049, 0.051), scale=(0.8, 1.0))
        session = PhaseStrippingSession(make_context(x, 0.1 * profile), bounds=bounds)

        preview = session.preview(candidate)
        session.accept_preview(preview)

        self.assertTrue(preview.warnings)
        self.assertLess(float(session.residual_y.min()), 0.0)
        self.assertTrue(np.array_equal(session.residual_y < 0.0, session.residual_y < 0.0))

    def test_history_undo_redo_and_reset_recompute_residuals_deterministically(self):
        x = np.arange(20.0, 29.001, 0.002)
        candidates = [
            make_candidate([Peak("a", 21.0, 100.0)]),
            make_candidate([Peak("b", 24.0, 100.0)]),
            make_candidate([Peak("c", 27.0, 100.0)]),
        ]
        for index, candidate in enumerate(candidates, start=1):
            candidate.candidate_id = f"phase-{index}"
        original = sum(
            project_candidate_profile(x, candidate, shift_deg=0.0, sigma_deg=0.05)
            for candidate in candidates
        )
        bounds = FitBounds(shift_deg=(-0.001, 0.001), sigma_deg=(0.049, 0.051), scale=(0.0, 1.1))
        session = PhaseStrippingSession(make_context(x, original), bounds=bounds)
        original_x_bytes = session.context.x.tobytes()
        original_y_bytes = session.context.intensity.tobytes()

        session.accept_preview(session.preview(candidates[0]))
        residual_after_a = session.residual_y
        session.accept_preview(session.preview(candidates[1]))
        residual_after_b = session.residual_y
        session.accept_preview(session.preview(candidates[2]))
        residual_after_c = session.residual_y

        session.undo()
        np.testing.assert_array_equal(session.residual_y, residual_after_b)
        session.redo()
        np.testing.assert_array_equal(session.residual_y, residual_after_c)
        session.undo()
        session.undo()
        np.testing.assert_array_equal(session.residual_y, residual_after_a)
        session.reset()
        np.testing.assert_array_equal(session.residual_y, original)
        self.assertEqual(session.context.x.tobytes(), original_x_bytes)
        self.assertEqual(session.context.intensity.tobytes(), original_y_bytes)

    def test_session_serialization_round_trip_preserves_arrays_and_state_exactly(self):
        x = np.arange(20.0, 25.001, 0.002)
        candidate = make_candidate([Peak("p1", 21.0, 100.0), Peak("p2", 23.0, 60.0)])
        candidate.candidate_id = "phase-a"
        original = 1.3 * project_candidate_profile(x, candidate, shift_deg=0.05, sigma_deg=0.06)
        session = PhaseStrippingSession(make_context(x, original))
        session.accept_preview(session.preview(candidate))
        session.exclude("phase-b")

        payload = json.loads(json.dumps(session.to_dict(), allow_nan=False))
        restored = PhaseStrippingSession.from_dict(payload)

        self.assertEqual(restored.excluded_candidate_ids, ("phase-b",))
        self.assertEqual(restored.accepted_operations, session.accepted_operations)
        np.testing.assert_array_equal(restored.context.x, session.context.x)
        np.testing.assert_array_equal(restored.context.intensity, session.context.intensity)
        np.testing.assert_array_equal(restored.fitted_total, session.fitted_total)
        np.testing.assert_array_equal(restored.residual_y, session.residual_y)

    def test_session_serialization_restores_original_dtype_endianness_and_bytes(self):
        x = np.array([20.0, 20.1, 20.2], dtype=">f4")
        intensity = np.array([12, 40, 18], dtype=">i2")
        session = PhaseStrippingSession(make_context(x, intensity))

        payload = json.loads(json.dumps(session.to_dict(), allow_nan=False))
        restored = PhaseStrippingSession.from_dict(payload)

        self.assertEqual(payload["context"]["x_dtype"], x.dtype.str)
        self.assertEqual(payload["context"]["intensity_dtype"], intensity.dtype.str)
        self.assertEqual(restored.context.x.dtype.str, x.dtype.str)
        self.assertEqual(restored.context.intensity.dtype.str, intensity.dtype.str)
        self.assertEqual(restored.context.x.tobytes(), x.tobytes())
        self.assertEqual(restored.context.intensity.tobytes(), intensity.tobytes())

    def test_bounds_serialization_distinguishes_finite_positive_and_negative_infinity(self):
        bounds = FitBounds(
            shift_deg=(float("-inf"), float("inf")),
            sigma_deg=(0.035, 0.25),
            scale=(float("-inf"), float("inf")),
        )
        x = np.array([20.0, 20.1, 20.2])
        session = PhaseStrippingSession(make_context(x, np.array([1.0, 2.0, 1.0])), bounds=bounds)

        payload = json.loads(json.dumps(session.to_dict(), allow_nan=False))
        restored = PhaseStrippingSession.from_dict(payload)

        self.assertEqual(payload["bounds"]["shift_deg"], ["-inf", "+inf"])
        self.assertEqual(payload["bounds"]["sigma_deg"], [0.035, 0.25])
        self.assertEqual(payload["bounds"]["scale"], ["-inf", "+inf"])
        self.assertTrue(np.isneginf(restored.bounds.shift_deg[0]))
        self.assertTrue(np.isposinf(restored.bounds.shift_deg[1]))
        self.assertTrue(np.isneginf(restored.bounds.scale[0]))
        self.assertTrue(np.isposinf(restored.bounds.scale[1]))

    def test_bounds_serialization_rejects_nan(self):
        x = np.array([20.0, 20.1, 20.2])
        bounds = FitBounds(shift_deg=(float("nan"), 0.25), sigma_deg=(0.035, 0.25), scale=(0.0, 1.0))
        session = PhaseStrippingSession(make_context(x, np.array([1.0, 2.0, 1.0])), bounds=bounds)

        with self.assertRaisesRegex(ValueError, "NaN"):
            session.to_dict()

    def test_manual_preview_is_bounded_and_can_be_cancelled_without_changing_residual(self):
        x = np.arange(20.0, 22.001, 0.002)
        candidate = make_candidate([Peak("p1", 21.0, 100.0)])
        session = PhaseStrippingSession(make_context(x, np.ones_like(x)))
        original = session.residual_y.copy()

        preview = session.preview_with_parameters(candidate, scale=0.5, shift_deg=0.1, sigma_deg=0.06)

        self.assertIs(session.current_preview, preview)
        np.testing.assert_array_equal(session.residual_y, original)
        self.assertTrue(session.cancel_preview())
        self.assertIsNone(session.current_preview)
        np.testing.assert_array_equal(session.residual_y, original)
        with self.assertRaisesRegex(ValueError, "scale"):
            session.preview_with_parameters(candidate, scale=-0.1, shift_deg=0.0, sigma_deg=0.06)
        with self.assertRaisesRegex(ValueError, "shift"):
            session.preview_with_parameters(candidate, scale=0.5, shift_deg=1.0, sigma_deg=0.06)

    def test_reset_restores_excluded_candidates_as_well_as_contributions(self):
        x = np.arange(20.0, 22.001, 0.002)
        candidate = make_candidate([Peak("p1", 21.0, 100.0)])
        session = PhaseStrippingSession(make_context(x, np.ones_like(x)))
        session.exclude(candidate)

        session.reset()

        self.assertEqual(session.excluded_candidate_ids, ())

    def test_preview_warns_when_fitted_shift_or_width_reaches_a_bound(self):
        x = np.arange(20.0, 22.001, 0.002)
        candidate = make_candidate([Peak("p1", 21.0, 100.0)])

        class BoundaryFitter:
            def fit(self, *_args, **_kwargs):
                return PhaseFit(candidate.candidate_id, shift_deg=-0.25, sigma_deg=0.25, scale=0.5)

        session = PhaseStrippingSession(make_context(x, np.ones_like(x)), fitter=BoundaryFitter())
        preview = session.preview(candidate)

        self.assertTrue(any("global shift" in warning and "lower bound" in warning for warning in preview.warnings))
        self.assertTrue(any("peak width" in warning and "upper bound" in warning for warning in preview.warnings))


if __name__ == "__main__":
    unittest.main()
