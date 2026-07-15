import tempfile
import tkinter as tk
import unittest
from unittest import mock

import numpy as np
from matplotlib import rcParams

from ophiuchus.phase_stripping.models import AnalysisContext
from ophiuchus.phase_stripping.window import PhaseStrippingWindow, open_or_raise_phase_stripping_window
from ophiuchus.xrd.models import Candidate, Peak


def make_window_inputs():
    x = np.arange(20.0, 24.001, 0.005)
    candidate = Candidate(
        candidate_id="phase-fege",
        formula_pretty="FeGe",
        source="library:local",
        source_path="FeGe.cif",
        elements=["Fe", "Ge"],
        structure_hash="structure-sha256",
        theory_peaks=[
            Peak("p1", 21.0, 100.0, hkl="100", d_spacing=4.2),
            Peak("p2", 23.0, 60.0, hkl="110", d_spacing=3.1),
        ],
        simulation_validation={"pattern_fingerprint": "pattern-sha256"},
        space_group_symbol="P4/mmm",
        simulation_state="ready",
    )
    intensity = 2.0 * np.exp(-0.5 * ((x - 21.0) / 0.06) ** 2)
    intensity += 1.2 * np.exp(-0.5 * ((x - 23.0) / 0.06) ** 2)
    context = AnalysisContext(
        x=x,
        intensity=intensity,
        radiation="CuKalpha12",
        wavelength_angstrom=1.54056,
        two_theta_range=(20.0, 24.0),
        tolerance_deg=0.15,
        source_path="sample.xy",
        source_fingerprint="source-sha256",
        data_fingerprint="data-sha256",
    )
    return context, [candidate]


class PhaseStrippingWindowTests(unittest.TestCase):
    def setUp(self):
        self.root = tk.Tk()
        self.root.withdraw()

    def tearDown(self):
        self.root.destroy()

    def test_constructs_three_panes_with_inherited_context_and_disabled_accept(self):
        context, candidates = make_window_inputs()
        window = PhaseStrippingWindow(self.root, context, candidates, element_scope=("Fe", "Ge"))
        self.root.update_idletasks()
        try:
            self.assertEqual(len(window.panes.panes()), 3)
            self.assertIs(window.session.context, context)
            np.testing.assert_array_equal(window.session.context.intensity, context.intensity)
            self.assertEqual(window.accept_button.instate(["disabled"]), True)
            self.assertEqual(len(window.figure.axes), 3)
            self.assertIs(window.axis, window.pattern_axis)
            self.assertTrue(window.pattern_axis.get_shared_x_axes().joined(window.pattern_axis, window.contribution_axis))
            self.assertTrue(window.pattern_axis.get_shared_x_axes().joined(window.pattern_axis, window.residual_axis))
            self.assertGreaterEqual(len(window.pattern_axis.lines), 2)
            self.assertGreaterEqual(len(window.residual_axis.collections), 2)
            self.assertIn("Microsoft YaHei UI", rcParams["font.sans-serif"])
            self.assertEqual(len(window.candidate_tree.get_children()), 1)
        finally:
            window.destroy()

    def test_default_layout_keeps_composition_panel_readable_and_readout_on_its_own_row(self):
        context, candidates = make_window_inputs()
        window = PhaseStrippingWindow(self.root, context, candidates, element_scope=("Fe", "Ge"))
        window.deiconify()
        window.update()
        try:
            self.assertGreaterEqual(window.right_pane.winfo_width(), 400)
            self.assertLessEqual(
                sum(int(window.composition_tree.column(column, "width")) for column in window.composition_tree["columns"]),
                350,
            )
            self.assertEqual(int(window.readout_label.grid_info()["row"]), 3)
            self.assertLessEqual(window.winfo_width(), window.winfo_screenwidth() - 60)
            self.assertLessEqual(window.winfo_height(), window.winfo_screenheight() - 90)
        finally:
            window.destroy()

    def test_preview_accept_undo_redo_refreshes_controls_and_signed_residual(self):
        context, candidates = make_window_inputs()
        window = PhaseStrippingWindow(self.root, context, candidates, element_scope=("Fe", "Ge"))
        self.root.update_idletasks()
        try:
            item = window.candidate_tree.get_children()[0]
            window.candidate_tree.selection_set(item)
            window._on_candidate_selected()
            window._auto_fit_preview()
            self.assertFalse(window.accept_button.instate(["disabled"]))
            window._accept_preview()
            self.assertEqual(len(window.session.accepted_operations), 1)
            self.assertFalse(window.undo_button.instate(["disabled"]))
            window._undo()
            self.assertEqual(len(window.session.accepted_operations), 0)
            window._redo()
            self.assertEqual(len(window.session.accepted_operations), 1)
            self.assertEqual(window.session.residual_y.shape, context.intensity.shape)
        finally:
            window._closing = True
            window.destroy()

    def test_click_populates_peak_composition_with_formula_and_reconstruction_rows(self):
        context, candidates = make_window_inputs()
        window = PhaseStrippingWindow(self.root, context, candidates, element_scope=("Fe", "Ge"))
        item = window.candidate_tree.get_children()[0]
        window.candidate_tree.selection_set(item)
        window._on_candidate_selected()
        window._auto_fit_preview()
        window._accept_preview()
        event = type("ChartEvent", (), {"inaxes": window.pattern_axis, "xdata": 21.002, "ydata": 2.0})()

        window._on_chart_click(event)
        rows = [window.composition_tree.item(row, "values") for row in window.composition_tree.get_children()]

        try:
            self.assertTrue(any(values[0] == "实验强度" for values in rows))
            self.assertTrue(any(values[0] == "已解释合计" for values in rows))
            self.assertTrue(any(values[0] == "FeGe" and values[3] == "100" for values in rows))
            self.assertIn("21.0000", window.readout_var.get())
        finally:
            window._closing = True
            window.destroy()

    def test_negligible_point_contribution_is_displayed_as_zero(self):
        self.assertEqual(PhaseStrippingWindow._format_point_intensity(4.075e-26), "0")
        self.assertEqual(PhaseStrippingWindow._format_point_intensity(-3.2e-12), "0")

    def test_chart_readout_uses_zero_for_negligible_reconstruction(self):
        context, candidates = make_window_inputs()
        window = PhaseStrippingWindow(self.root, context, candidates, element_scope=("Fe", "Ge"))
        item = window.candidate_tree.get_children()[0]
        window.candidate_tree.selection_set(item)
        window._on_candidate_selected()
        window._auto_fit_preview()
        window._accept_preview()
        event = type("ChartEvent", (), {"inaxes": window.pattern_axis, "xdata": 22.0, "ydata": 0.0})()

        window._on_chart_click(event)

        try:
            self.assertIn("| 解释 0 |", window.readout_var.get())
        finally:
            window._closing = True
            window.destroy()

    def test_open_or_raise_reuses_existing_window(self):
        context, candidates = make_window_inputs()
        first = open_or_raise_phase_stripping_window(
            self.root, None, context=context, candidates=candidates, element_scope=("Fe", "Ge")
        )
        second = open_or_raise_phase_stripping_window(
            self.root, first, context=context, candidates=candidates, element_scope=("Fe", "Ge")
        )
        try:
            self.assertIs(first, second)
        finally:
            first._closing = True
            first.destroy()

    def test_close_with_accepted_operation_can_discard_without_export(self):
        context, candidates = make_window_inputs()
        window = PhaseStrippingWindow(self.root, context, candidates, element_scope=("Fe", "Ge"))
        item = window.candidate_tree.get_children()[0]
        window.candidate_tree.selection_set(item)
        window._on_candidate_selected()
        window._auto_fit_preview()
        window._accept_preview()

        with mock.patch("ophiuchus.phase_stripping.window.messagebox.askyesnocancel", return_value=False):
            window._request_close()

        self.assertFalse(window.winfo_exists())
        self.assertIsNone(window._last_export_path)


if __name__ == "__main__":
    unittest.main()
