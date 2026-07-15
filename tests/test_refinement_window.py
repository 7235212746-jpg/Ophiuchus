import tkinter as tk
import time
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from ophiuchus.phase_stripping.models import AnalysisContext
from ophiuchus.refinement.window import RefinementWindow, SampleConclusionDialog
from ophiuchus.refinement.conclusion import (
    ImpuritySignalEstimate,
    SampleConclusion,
    build_sample_conclusion,
)
from ophiuchus.refinement.oxide_candidates import ControlledOxideLoadResult
from ophiuchus.xrd.models import Candidate, Peak
from ophiuchus.xrd.refinement import RietanRefinementResult


def make_inputs():
    x = np.linspace(20.0, 24.0, 401)
    intensity = 100.0 * np.exp(-0.5 * ((x - 22.0) / 0.08) ** 2) + 10.0
    context = AnalysisContext(
        x=x,
        intensity=intensity,
        radiation="CuKa",
        wavelength_angstrom=1.540593,
        two_theta_range=(20.0, 24.0),
        tolerance_deg=0.15,
        source_path="sample.xy",
        source_fingerprint="source",
        data_fingerprint="data",
    )
    candidate = Candidate(
        "target",
        "FeGe",
        "library:local",
        "target.cif",
        ["Fe", "Ge"],
        "structure",
        theory_peaks=[Peak("p1", 22.0, 100.0)],
    )
    return context, [candidate]


class FakeBackend:
    available = True

    def refine(self, cif_path, x, intensity, settings):
        observed = np.asarray(intensity, dtype=float)
        background = np.full_like(observed, 10.0)
        calculated = observed * 0.95
        return RietanRefinementResult(
            two_theta_deg=x,
            observed_intensity=observed,
            calculated_intensity=calculated,
            residual_intensity=observed - calculated,
            background_intensity=background,
            reflection_two_theta_deg=(22.0,),
            rwp_percent=8.5,
            rp_percent=6.1,
            goodness_of_fit=1.2,
            s_value=1.0954,
            parameters={"zero_shift": 0.01, "scale": 1e-5},
            provenance={"engine": "RIETAN-FP constrained refinement"},
        )


class RefinementWindowTests(unittest.TestCase):
    def setUp(self):
        self.root = tk.Tk()
        self.root.withdraw()

    def tearDown(self):
        self.root.destroy()

    def test_constructs_constrained_workbench_with_locked_advanced_parameters(self):
        context, candidates = make_inputs()
        window = RefinementWindow(self.root, context, candidates, backend=FakeBackend())
        self.root.update_idletasks()
        try:
            self.assertEqual(window.title(), "RIETAN-FP 受约束精修")
            self.assertEqual(len(window.figure.axes), 2)
            self.assertIn("FeGe", window.target_var.get())
            self.assertIn("原子坐标", window.locked_parameters_var.get())
            self.assertEqual(window.mode_combo["values"], ("背景 + 尺度", "加零点偏移", "加峰宽 U/V/W"))
            self.assertTrue(hasattr(window, "control_canvas"))
            self.assertTrue(hasattr(window, "control_scrollbar"))
            self.assertTrue(window.export_button.instate(["disabled"]))
            self.assertTrue(window.conclusion_button.instate(["disabled"]))
            self.assertFalse(window.start_button.instate(["disabled"]))
        finally:
            window.destroy()

    def test_destroy_cancels_pending_queue_poll_callback(self):
        context, candidates = make_inputs()
        window = RefinementWindow(self.root, context, candidates, backend=FakeBackend())
        callback_id = window._poll_after_id
        window.canvas.draw_idle()
        draw_id = window.canvas._idle_draw_id

        window.destroy()

        pending = self.root.tk.call("after", "info")
        self.assertNotIn(callback_id, pending)
        self.assertNotIn(draw_id, pending)

    def test_result_plot_contains_observed_calculated_background_difference_and_ticks(self):
        context, candidates = make_inputs()
        window = RefinementWindow(self.root, context, candidates, backend=FakeBackend())
        try:
            result = FakeBackend().refine("target.cif", context.x, context.intensity, SimpleNamespace())
            window._show_result(result)

            labels = {line.get_label() for line in window.pattern_axis.lines}
            self.assertIn("实验谱", labels)
            self.assertIn("计算谱", labels)
            self.assertIn("背景", labels)
            self.assertIn("差值", {line.get_label() for line in window.residual_axis.lines})
            self.assertGreaterEqual(len(window.pattern_axis.collections), 1)
            self.assertIn("Rwp 8.500%", window.metrics_var.get())
            self.assertIn("S 1.0954", window.metrics_var.get())
            self.assertIn("GofF 1.2000", window.metrics_var.get())
            self.assertFalse(window.export_button.instate(["disabled"]))
            self.root.update_idletasks()
            self.assertGreater(window.control_canvas.yview()[0], 0.0)
        finally:
            window.destroy()

    def test_conclusion_button_opens_scrollable_sample_report_dialog(self):
        context, candidates = make_inputs()
        window = RefinementWindow(self.root, context, candidates, backend=FakeBackend())
        try:
            result = FakeBackend().refine("target.cif", context.x, context.intensity, SimpleNamespace())
            conclusion = build_sample_conclusion(result, candidates[0], [])
            window._show_result(result, conclusion)
            window.conclusion_button.invoke()
            self.root.update_idletasks()

            self.assertFalse(window.conclusion_button.instate(["disabled"]))
            self.assertIsNotNone(window.conclusion_dialog)
            self.assertEqual(window.conclusion_dialog.title(), "Ophi 样品结论")
            report = window.conclusion_dialog.report_text.get("1.0", "end")
            self.assertIn("不能报告可信 wt%", report)
            self.assertIn("判断流程", report)
        finally:
            if getattr(window, "conclusion_dialog", None) is not None:
                window.conclusion_dialog.destroy()
            window.destroy()

    def test_worker_runs_controlled_oxide_second_pass_from_library(self):
        context, candidates = make_inputs()
        oxide = Candidate(
            "hematite",
            "Fe2O3",
            "library:materials_project",
            "hematite.cif",
            ["Fe", "O"],
            "oxide",
            theory_peaks=[Peak("o1", 21.0, 100.0), Peak("o2", 23.0, 50.0)],
        )
        oxide_result = ControlledOxideLoadResult((oxide,), ("Fe2O3",), ())
        extra = Candidate("tin", "Sn", "library", "tin.cif", ["Sn"], "tin")
        window = RefinementWindow(
            self.root,
            context,
            candidates,
            backend=FakeBackend(),
            supporting_candidates=[extra],
            oxide_library_path="library.sqlite",
        )
        try:
            with patch("ophiuchus.refinement.window.discover_sibling_peak_references", return_value=[]) as discover, patch(
                    "ophiuchus.refinement.window.load_controlled_oxide_candidates",
                    return_value=oxide_result,
                ) as loader:
                window._worker(candidates[0], SimpleNamespace())

            kind, payload = window._queue.get_nowait()
            self.assertEqual(kind, "ok")
            self.assertEqual(payload[1].oxide_formulas_checked, ("Fe2O3",))
            self.assertEqual(discover.call_args.kwargs["allowed_elements"], {"Fe", "Ge"})
            loader.assert_called_once_with(
                "library.sqlite",
                {"Fe", "Ge"},
                radiation="CuKa",
                two_theta_range=(20.0, 24.0),
            )
        finally:
            window.destroy()

    def test_conclusion_dialog_copies_report_and_opens_excludable_candidate_evidence(self):
        context, candidates = make_inputs()
        result = FakeBackend().refine("target.cif", context.x, context.intensity, SimpleNamespace())
        estimate = ImpuritySignalEstimate(
            candidate_id="imp",
            formula="Fe3Ge",
            source_label="library:materials_project",
            source_path="imp.cif",
            signal_share_percent=8.3,
            weight_fraction_percent=None,
            rwp_proxy_percent=5.2,
            rwp_improvement=3.1,
            shift_deg=-0.13,
            fwhm_deg=0.30,
            evidence_label="不稳定候选（峰宽达到搜索边界）",
            peak_count=5,
            is_oxide=False,
            peak_width_at_boundary=True,
            matched_residual_peaks_deg=(35.7, 44.4),
            missing_strong_theory_peaks_deg=(65.4,),
        )
        conclusion = SampleConclusion(
            target_formula="ZrFe6Ge4",
            target_evidence_label="目标相得到强支持",
            refinement=result,
            target_signal_share_percent=62.9,
            impurity_estimates=(estimate,),
            competing_models=False,
        )
        excluded = []
        self.root.deiconify()
        dialog = SampleConclusionDialog(self.root, conclusion, on_exclude_formula=excluded.append)
        try:
            dialog.copy_button.invoke()
            self.root.update()
            self.assertIn("ZrFe6Ge4", self.root.clipboard_get())

            dialog.evidence_button.invoke()
            self.root.update()
            evidence = dialog.evidence_dialog
            self.assertIsNotNone(evidence)
            self.assertEqual(evidence.title(), "Ophi 候选证据")
            self.assertTrue(evidence.exclude_button.winfo_ismapped())
            gain_bbox = evidence.tree.bbox("imp", "#3")
            self.assertTrue(gain_bbox)
            self.assertLessEqual(gain_bbox[0] + gain_bbox[2], evidence.tree.winfo_width())
            detail = evidence.detail_text.get("1.0", "end")
            self.assertIn("命中残差峰", detail)
            self.assertIn("35.70", detail)
            evidence.exclude_button.invoke()
            self.assertEqual(excluded, ["Fe3Ge"])
        finally:
            if getattr(dialog, "evidence_dialog", None) is not None:
                dialog.evidence_dialog.destroy()
            dialog.destroy()

    def test_excluding_candidate_recalculates_conclusion_without_rerunning_rietan(self):
        context, candidates = make_inputs()
        target = candidates[0]
        excluded_candidate = Candidate("a", "Fe3Ge", "library", "a.cif", ["Fe", "Ge"], "a")
        excluded_candidate.theory_peaks = [Peak("a1", 21.0, 100.0), Peak("a2", 23.0, 50.0)]
        retained_candidate = Candidate("b", "Fe13Ge3", "library", "b.cif", ["Fe", "Ge"], "b")
        retained_candidate.theory_peaks = [Peak("b1", 21.5, 100.0), Peak("b2", 23.5, 50.0)]
        window = RefinementWindow(self.root, context, candidates, backend=FakeBackend())
        try:
            result = FakeBackend().refine("target.cif", context.x, context.intensity, SimpleNamespace())
            window.result = result
            window._active_candidate = target
            window._conclusion_candidates = [excluded_candidate, retained_candidate]
            window._exclude_candidate_formula("Fe3Ge")
            deadline = time.time() + 5.0
            while window._queue.empty() and time.time() < deadline:
                time.sleep(0.02)
            kind, conclusion = window._queue.get_nowait()

            self.assertEqual(kind, "conclusion_ok")
            self.assertIn("Fe3Ge", conclusion.excluded_formulas)
            self.assertNotIn("Fe3Ge", [item.formula for item in conclusion.impurity_estimates])
            self.assertEqual(window.backend.__class__, FakeBackend)
        finally:
            window.destroy()


if __name__ == "__main__":
    unittest.main()
