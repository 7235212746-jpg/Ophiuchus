import tempfile
import unittest
from pathlib import Path

import numpy as np

from ophiuchus.refinement.conclusion import (
    build_sample_conclusion,
    discover_sibling_peak_references,
    estimate_impurity_signals,
)
from ophiuchus.xrd.models import Candidate, Peak
from ophiuchus.xrd.refinement import RietanRefinementResult


def gaussian(x, center, amplitude, fwhm=0.18):
    sigma = fwhm / 2.354820045
    return amplitude * np.exp(-0.5 * ((x - center) / sigma) ** 2)


class RefinementConclusionTests(unittest.TestCase):
    def test_discovers_peak_reference_beside_experimental_pattern_and_excludes_target(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            xrd = root / "980.asc"
            xrd.write_text("3 1\n4 2\n", encoding="utf-8")
            (root / "Fe1.67Ge.txt").write_text(
                "2-Theta d I\n35.743 2.51 8\n44.415 2.038 100\n45.020 2.012 75\n",
                encoding="utf-8",
            )
            (root / "MONI.int").write_text("35.36 100\n44.78 80\n", encoding="utf-8")

            references = discover_sibling_peak_references(
                xrd,
                allowed_elements={"Zr", "Fe", "Ge"},
                target_formula="ZrFe6Ge4",
            )

        self.assertEqual([candidate.formula_pretty for candidate in references], ["Fe1.67Ge"])
        self.assertEqual(len(references[0].theory_peaks), 3)

    def test_incremental_fit_recovers_impurity_signal_share_without_calling_it_weight_fraction(self):
        x = np.linspace(20.0, 70.0, 2501)
        background = np.full_like(x, 30.0)
        target_profile = gaussian(x, 30.0, 1000.0)
        impurity_profile = gaussian(x, 45.0, 100.0) + gaussian(x, 55.0, 50.0)
        observed = background + target_profile + impurity_profile
        result = RietanRefinementResult(
            two_theta_deg=x,
            observed_intensity=observed,
            calculated_intensity=background + target_profile,
            residual_intensity=impurity_profile,
            background_intensity=background,
            reflection_two_theta_deg=(30.0,),
            rwp_percent=8.0,
            rp_percent=6.0,
            goodness_of_fit=1.4,
            s_value=1.1832,
        )
        target = Candidate("target", "Target", "local", "target.cif", ["A"], None)
        impurity = Candidate("imp", "Impurity", "local_peak_list", "impurity.txt", ["B"], None)
        impurity.theory_peaks = [Peak("i1", 45.0, 100.0), Peak("i2", 55.0, 50.0)]

        estimates = estimate_impurity_signals(result, target, [target, impurity])

        self.assertEqual(len(estimates), 1)
        self.assertEqual(estimates[0].formula, "Impurity")
        self.assertGreater(estimates[0].rwp_improvement, 1.0)
        self.assertGreater(estimates[0].signal_share_percent, 5.0)
        self.assertLess(estimates[0].signal_share_percent, 25.0)
        self.assertIsNone(estimates[0].weight_fraction_percent)

    def test_impurity_evidence_lists_matched_residual_and_missing_strong_theory_peaks(self):
        x = np.linspace(20.0, 70.0, 2501)
        background = np.full_like(x, 30.0)
        target_profile = gaussian(x, 30.0, 1000.0)
        impurity_profile = gaussian(x, 45.0, 100.0) + gaussian(x, 55.0, 50.0)
        result = RietanRefinementResult(
            two_theta_deg=x,
            observed_intensity=background + target_profile + impurity_profile,
            calculated_intensity=background + target_profile,
            residual_intensity=impurity_profile,
            background_intensity=background,
            reflection_two_theta_deg=(30.0,),
            rwp_percent=8.0,
            rp_percent=6.0,
            goodness_of_fit=1.4,
        )
        target = Candidate("target", "Target", "local", "target.cif", ["A"], None)
        impurity = Candidate("imp", "Impurity", "library:materials_project", "impurity.cif", ["B"], None)
        impurity.theory_peaks = [
            Peak("i1", 45.0, 100.0),
            Peak("i2", 55.0, 50.0),
            Peak("i3", 62.0, 80.0),
        ]

        estimate = estimate_impurity_signals(result, target, [impurity])[0]
        evidence = estimate.to_chinese_evidence_text()

        self.assertEqual(estimate.candidate_id, "imp")
        self.assertEqual(estimate.source_label, "library:materials_project")
        self.assertTrue(any(abs(value - 45.0) < 0.1 for value in estimate.matched_residual_peaks_deg))
        self.assertTrue(any(abs(value - 55.0) < 0.1 for value in estimate.matched_residual_peaks_deg))
        self.assertTrue(any(abs(value - 62.0) < 0.1 for value in estimate.missing_strong_theory_peaks_deg))
        self.assertIn("命中残差峰", evidence)
        self.assertIn("缺失理论强峰", evidence)
        self.assertIn("统一峰位偏移", evidence)
        self.assertIn("FWHM", evidence)

    def test_conclusion_text_reports_target_evidence_competing_models_and_quantitation_limit(self):
        x = np.array([20.0, 30.0, 40.0])
        result = RietanRefinementResult(
            two_theta_deg=x,
            observed_intensity=np.array([10.0, 20.0, 10.0]),
            calculated_intensity=np.array([10.0, 18.0, 10.0]),
            residual_intensity=np.array([0.0, 2.0, 0.0]),
            background_intensity=np.array([5.0, 5.0, 5.0]),
            reflection_two_theta_deg=(30.0,),
            rwp_percent=7.13,
            rp_percent=4.69,
            goodness_of_fit=15.08,
            s_value=3.884,
        )
        target = Candidate("target", "ZrFe6Ge4", "local", "target.cif", ["Zr", "Fe", "Ge"], None)

        conclusion = build_sample_conclusion(
            result,
            target,
            [],
            oxide_formulas_checked=("Fe2O3", "ZrO2", "GeO2"),
        )
        text = conclusion.to_chinese_text()

        self.assertIn("ZrFe6Ge4", text)
        self.assertIn("Rwp 7.130%", text)
        self.assertIn("不能报告可信 wt%", text)
        self.assertIn("氧化物二次筛选", text)
        self.assertIn("Fe2O3、ZrO2、GeO2", text)
        self.assertIn("未见能稳定改善残差的常见氧化物", text)
        self.assertIn("判断流程", text)

    def test_low_rwp_without_target_signal_does_not_claim_strong_target_support(self):
        x = np.linspace(20.0, 50.0, 601)
        background = np.full_like(x, 20.0)
        observed = background + gaussian(x, 35.0, 100.0)
        result = RietanRefinementResult(
            two_theta_deg=x,
            observed_intensity=observed,
            calculated_intensity=background,
            residual_intensity=observed - background,
            background_intensity=background,
            reflection_two_theta_deg=(30.0,),
            rwp_percent=6.0,
            rp_percent=5.0,
            goodness_of_fit=1.2,
            parameters={"scale": 0.0},
        )
        target = Candidate("target", "Target", "local", "target.cif", ["A"], None)

        conclusion = build_sample_conclusion(result, target, [])

        self.assertEqual(conclusion.target_signal_share_percent, 0.0)
        self.assertNotIn("强支持", conclusion.target_evidence_label)
        self.assertIn("目标相信号不足", conclusion.target_evidence_label)

    def test_candidate_at_peak_width_search_boundary_is_not_labeled_strong(self):
        x = np.linspace(20.0, 60.0, 2001)
        background = np.full_like(x, 20.0)
        target_profile = gaussian(x, 30.0, 500.0, fwhm=0.18)
        broad_impurity = gaussian(x, 45.0, 160.0, fwhm=0.55) + gaussian(x, 52.0, 80.0, fwhm=0.55)
        observed = background + target_profile + broad_impurity
        result = RietanRefinementResult(
            two_theta_deg=x,
            observed_intensity=observed,
            calculated_intensity=background + target_profile,
            residual_intensity=broad_impurity,
            background_intensity=background,
            reflection_two_theta_deg=(30.0,),
            rwp_percent=8.0,
            rp_percent=6.0,
            goodness_of_fit=1.4,
        )
        target = Candidate("target", "Target", "local", "target.cif", ["A"], None)
        impurity = Candidate("imp", "Impurity", "local", "imp.cif", ["B"], None)
        impurity.theory_peaks = [Peak("i1", 45.0, 100.0), Peak("i2", 52.0, 50.0)]

        estimates = estimate_impurity_signals(result, target, [impurity])

        self.assertEqual(len(estimates), 1)
        self.assertTrue(estimates[0].peak_width_at_boundary)
        self.assertIn("不稳定候选", estimates[0].evidence_label)
        self.assertNotIn("较强候选", estimates[0].evidence_label)

    def test_osmium_candidate_is_not_misreported_as_an_oxide(self):
        x = np.linspace(20.0, 60.0, 2001)
        background = np.full_like(x, 20.0)
        target_profile = gaussian(x, 30.0, 500.0)
        osmium_profile = gaussian(x, 45.0, 120.0) + gaussian(x, 52.0, 60.0)
        result = RietanRefinementResult(
            two_theta_deg=x,
            observed_intensity=background + target_profile + osmium_profile,
            calculated_intensity=background + target_profile,
            residual_intensity=osmium_profile,
            background_intensity=background,
            reflection_two_theta_deg=(30.0,),
            rwp_percent=8.0,
            rp_percent=6.0,
            goodness_of_fit=1.4,
        )
        target = Candidate("target", "Target", "local", "target.cif", ["A"], None)
        osmium = Candidate("os", "Os", "local", "os.cif", ["Os"], None)
        osmium.theory_peaks = [Peak("o1", 45.0, 100.0), Peak("o2", 52.0, 50.0)]

        conclusion = build_sample_conclusion(
            result,
            target,
            [osmium],
            oxide_formulas_checked=("Fe2O3",),
        )

        self.assertIn("Os", conclusion.to_chinese_text())
        self.assertIn("未见能稳定改善残差的常见氧化物", conclusion.to_chinese_text())


if __name__ == "__main__":
    unittest.main()
