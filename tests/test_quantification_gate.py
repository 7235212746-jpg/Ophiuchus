from dataclasses import replace
from pathlib import Path
import tempfile
import unittest

import numpy as np

from ophiuchus.xrd.multiphase_models import (
    MultiphaseRefinementResult,
    PhaseRefinementInput,
    PhaseRefinementResult,
)
from ophiuchus.xrd.quantification import (
    GateLevel,
    assess_quantification,
    initial_scale_variants,
)


def make_result(
    weights=(70.0, 30.0),
    *,
    scales=(2.0, 1.0),
    residual=None,
    rwp=8.0,
    zmv_difference=0.1,
):
    x = np.linspace(10.0, 20.0, 101)
    residual_values = np.zeros_like(x) if residual is None else np.asarray(residual, dtype=float)
    observed = np.full_like(x, 100.0)
    background = np.full_like(x, 5.0)
    calculated = observed - residual_values
    phases = (
        PhaseRefinementResult("main", "Main", scales[0], 1.0, 100.0, 50.0, weights[0], reflection_two_theta_deg=(12.0, 15.0)),
        PhaseRefinementResult("imp", "Impurity", scales[1], 1.0, 100.0, 50.0, weights[1], reflection_two_theta_deg=(13.0, 17.0)),
    )
    return MultiphaseRefinementResult(
        x,
        observed,
        calculated,
        residual_values,
        background,
        phases,
        rwp,
        6.0,
        1.2,
        provenance={"zmv_max_difference_wt_percent": zmv_difference},
    )


class QuantificationGateTests(unittest.TestCase):
    def test_negative_scale_hides_weight_fraction(self):
        result = make_result(scales=(-1.0, 1.0))

        assessment = assess_quantification(result, [make_result(), make_result()], {"main": "passed", "imp": "passed"})

        self.assertIs(assessment.level, GateLevel.FAIL)
        self.assertFalse(assessment.allow_weight_percent)
        self.assertTrue(any(item.code == "invalid_scale" for item in assessment.findings))

    def test_missing_stability_runs_is_a_hard_failure(self):
        assessment = assess_quantification(make_result(), [], {"main": "passed", "imp": "passed"})

        self.assertIs(assessment.level, GateLevel.FAIL)
        self.assertTrue(any(item.code == "stability_missing" for item in assessment.findings))

    def test_unstable_repeat_runs_fail_quantification(self):
        assessment = assess_quantification(
            make_result((70.0, 30.0)),
            [make_result((55.0, 45.0)), make_result((82.0, 18.0))],
            {"main": "passed", "imp": "passed"},
            stability_spread_limit_percent=5.0,
        )

        self.assertIs(assessment.level, GateLevel.FAIL)
        self.assertTrue(any(item.code == "initial_scale_instability" for item in assessment.findings))
        self.assertEqual(assessment.stability_ranges["main"], (55.0, 82.0))

    def test_stable_validated_repeats_allow_experimental_weight_fraction(self):
        assessment = assess_quantification(
            make_result((70.0, 30.0)),
            [make_result((70.5, 29.5)), make_result((69.8, 30.2))],
            {"main": "passed", "imp": "passed"},
            stability_spread_limit_percent=5.0,
        )

        self.assertIs(assessment.level, GateLevel.PASS)
        self.assertTrue(assessment.allow_weight_percent)
        self.assertEqual(assessment.label, "实验性定量")

    def test_unexplained_residual_peak_group_is_a_hard_failure(self):
        residual = np.zeros(101)
        residual[20] = 15.0
        residual[21] = 5.0
        residual[70] = 12.0
        assessment = assess_quantification(
            make_result(residual=residual),
            [make_result(residual=residual), make_result(residual=residual)],
            {"main": "passed", "imp": "passed"},
        )

        self.assertIs(assessment.level, GateLevel.FAIL)
        self.assertTrue(any(item.code == "unexplained_residual_group" for item in assessment.findings))

    def test_validation_and_zmv_disagreement_are_hard_failures(self):
        assessment = assess_quantification(
            make_result(zmv_difference=2.0),
            [make_result(zmv_difference=2.0), make_result(zmv_difference=2.0)],
            {"main": "passed", "imp": "failed"},
        )

        self.assertIs(assessment.level, GateLevel.FAIL)
        codes = {item.code for item in assessment.findings}
        self.assertIn("pattern_validation_failed", codes)
        self.assertIn("zmv_crosscheck_failed", codes)

    def test_initial_scale_variants_are_deterministic_and_preserve_roles(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = root / "first.cif"
            second = root / "second.cif"
            first.write_text("data_first\n", encoding="ascii")
            second.write_text("data_second\n", encoding="ascii")
            phases = (
                PhaseRefinementInput("main", "Main", first, "target"),
                PhaseRefinementInput("imp", "Impurity", second, "impurity"),
            )

            variants = initial_scale_variants(phases)

        self.assertEqual(len(variants), 3)
        self.assertEqual([[item.initial_scale for item in variant] for variant in variants], [[1.0, 1.0], [10.0, 1.0], [1.0, 5.0]])
        self.assertTrue(all(variant[0].role == "target" and variant[1].role == "impurity" for variant in variants))


if __name__ == "__main__":
    unittest.main()
