import unittest

import numpy as np

from ophiuchus.phase_stripping.reconstruction_plot import (
    build_reconstruction_view,
    cursor_breakdown,
)


class ReconstructionPlotTests(unittest.TestCase):
    def test_builds_additive_stacks_signed_residual_and_overflow(self):
        x = np.array([10.0, 20.0, 30.0])
        original = np.array([3.0, 5.0, 4.0])
        first = np.array([1.0, 2.0, 1.0])
        second = np.array([0.5, 4.0, 0.0])

        view = build_reconstruction_view(
            x,
            original,
            [first, second],
            labels=["FeGe", "Sn"],
            colors=["#008f7a", "#d88700"],
        )

        np.testing.assert_allclose(view.phase_sum, [1.5, 6.0, 1.0])
        np.testing.assert_allclose(view.residual, [1.5, -1.0, 3.0])
        np.testing.assert_allclose(view.overflow, [0.0, 1.0, 0.0])
        np.testing.assert_allclose(view.layers[0].stacked_lower, [0.0, 0.0, 0.0])
        np.testing.assert_allclose(view.layers[0].stacked_upper, first)
        np.testing.assert_allclose(view.layers[1].stacked_lower, first)
        np.testing.assert_allclose(view.layers[1].stacked_upper, first + second)
        self.assertEqual(view.layers[0].label, "FeGe")
        self.assertEqual(view.layers[1].color, "#d88700")

    def test_background_is_part_of_reconstruction_but_not_phase_sum(self):
        x = np.array([10.0, 20.0, 30.0])
        original = np.array([11.0, 17.0, 13.0])
        background = np.array([10.0, 11.0, 12.0])
        phase = np.array([1.0, 4.0, 1.0])

        view = build_reconstruction_view(
            x,
            original,
            [phase],
            labels=["FeGe"],
            colors=["#111111"],
            background=background,
        )

        np.testing.assert_array_equal(view.background, background)
        np.testing.assert_array_equal(view.corrected, original - background)
        np.testing.assert_array_equal(view.phase_sum, phase)
        np.testing.assert_array_equal(view.reconstructed, background + phase)
        np.testing.assert_array_equal(view.residual, original - background - phase)
        np.testing.assert_array_equal(view.overflow, [0.0, 0.0, 0.0])

    def test_phase_lanes_use_local_normalization_without_changing_raw_values(self):
        x = np.array([10.0, 20.0, 30.0])
        original = np.array([100.0, 100.0, 100.0])
        weak = np.array([0.0, 2.0, 1.0])
        strong = np.array([0.0, 100.0, 50.0])

        view = build_reconstruction_view(
            x,
            original,
            [weak, strong],
            labels=["weak", "strong"],
            colors=["#111111", "#222222"],
        )

        np.testing.assert_allclose(view.layers[0].lane_values, [0.0, 1.0, 0.5])
        np.testing.assert_allclose(view.layers[1].lane_values, [0.0, 1.0, 0.5])
        np.testing.assert_array_equal(view.layers[0].values, weak)
        np.testing.assert_array_equal(view.layers[1].values, strong)
        self.assertFalse(view.layers[0].values.flags.writeable)

    def test_preview_is_not_added_to_accepted_sum_but_has_a_proposed_sum(self):
        x = np.array([10.0, 20.0])
        original = np.array([5.0, 8.0])
        accepted = np.array([1.0, 3.0])
        preview = np.array([2.0, 4.0])

        view = build_reconstruction_view(
            x,
            original,
            [accepted],
            labels=["accepted"],
            colors=["#111111"],
            preview=preview,
            preview_label="preview",
            preview_color="#00aa55",
        )

        np.testing.assert_array_equal(view.phase_sum, accepted)
        np.testing.assert_array_equal(view.proposed_phase_sum, accepted + preview)
        np.testing.assert_array_equal(view.residual, original - accepted)
        self.assertEqual(view.preview_layer.label, "preview")

    def test_cursor_breakdown_uses_nearest_grid_point_and_reports_every_component(self):
        view = build_reconstruction_view(
            np.array([10.0, 20.0, 30.0]),
            np.array([3.0, 5.0, 4.0]),
            [np.array([1.0, 2.0, 1.0]), np.array([0.5, 1.0, 0.0])],
            labels=["FeGe", "Sn"],
            colors=["#111111", "#222222"],
        )

        values = cursor_breakdown(view, 19.6)

        self.assertEqual(values.index, 1)
        self.assertEqual(values.two_theta, 20.0)
        self.assertEqual(values.experimental, 5.0)
        self.assertEqual(values.phase_sum, 3.0)
        self.assertEqual(values.residual, 2.0)
        self.assertEqual(values.contributions, (("FeGe", 2.0), ("Sn", 1.0)))

    def test_rejects_mismatched_shapes_and_metadata(self):
        with self.assertRaisesRegex(ValueError, "same shape"):
            build_reconstruction_view(
                np.array([1.0, 2.0]),
                np.array([1.0]),
                [],
                labels=[],
                colors=[],
            )
        with self.assertRaisesRegex(ValueError, "labels and colors"):
            build_reconstruction_view(
                np.array([1.0]),
                np.array([1.0]),
                [np.array([1.0])],
                labels=[],
                colors=[],
            )


if __name__ == "__main__":
    unittest.main()
