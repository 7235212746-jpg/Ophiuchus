import json
import random
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from ophiuchus.xrd.candidates import LocalCandidateProvider, read_peak_list_file
from ophiuchus.xrd.matching import explain_multiphase, prioritize_impurities_after_main, score_candidate
from ophiuchus.xrd.models import Candidate, Peak
from ophiuchus.xrd.report import write_reports
from ophiuchus.xrd.simulation_trust import apply_vesta_trust_check


class CandidateMatchingReportTests(unittest.TestCase):
    def test_local_provider_filters_by_allowed_elements_and_int_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "ZrFe6Ge4.cif").write_text(
                "_cell_length_a 5\n_atom_site_type_symbol\nZr\nFe\nGe\n",
                encoding="utf-8",
            )
            (root / "Fe2O3.int").write_text("20 100\n31 55\n44 15\n", encoding="utf-8")
            provider = LocalCandidateProvider([root], allowed_elements={"Zr", "Fe", "Ge"})
            candidates = list(provider.iter_candidates())
        formulas = {c.formula_pretty for c in candidates}
        self.assertIn("ZrFe6Ge4", formulas)
        self.assertNotIn("Fe2O3", formulas)

    def test_read_peak_list_file_normalizes_intensity(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "phase.int"
            path.write_text("GENERAL$\n3\n20 10 0\n25 20 0\n", encoding="utf-8")
            peaks = read_peak_list_file(path)
        self.assertEqual([round(p.two_theta, 2) for p in peaks], [20.0, 25.0])
        self.assertAlmostEqual(max(p.intensity for p in peaks), 100.0)

    def test_dense_int_pattern_is_peak_picked_before_matching(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "MONI.int"
            rows = []
            for i in range(3000):
                x = 10.0 + i * 0.02
                y = 0.1
                for center, amp in [(20.0, 100.0), (30.0, 60.0), (44.4, 30.0)]:
                    y += amp * pow(2.718281828, -((x - center) ** 2) / (2 * 0.05**2))
                rows.append(f"{x:.3f} {y:.6f} 0")
            path.write_text("GENERAL$\n" + "\n".join(rows), encoding="utf-8")
            peaks = read_peak_list_file(path, intensity_threshold=5.0)
        self.assertLessEqual(len(peaks), 10)
        self.assertTrue(any(abs(p.two_theta - 20.0) < 0.08 for p in peaks))
        self.assertTrue(any(abs(p.two_theta - 30.0) < 0.08 for p in peaks))
        self.assertTrue(any(abs(p.two_theta - 44.4) < 0.08 for p in peaks))

    def test_generic_peak_list_name_uses_parent_sample_formula(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "ZrFe6Ge4 in different Temperature"
            root.mkdir()
            (root / "MONI.int").write_text("20 100\n30 50\n", encoding="utf-8")
            provider = LocalCandidateProvider([root], allowed_elements={"Zr", "Fe", "Ge"})
            candidates = list(provider.iter_candidates())
        self.assertEqual(len(candidates), 1)
        self.assertEqual(candidates[0].formula_pretty, "ZrFe6Ge4")

    def test_candidate_with_multiple_real_matches_beats_isolated_match(self):
        exp = [Peak("e1", 20.0, 100.0), Peak("e2", 30.0, 70.0), Peak("e3", 44.4, 35.0)]
        good = Candidate("good", "GoodPhase", "local", "good.int", ["Zr", "Fe"], None)
        good.theory_peaks = [Peak("t1", 20.02, 100.0), Peak("t2", 30.01, 60.0), Peak("t3", 52.0, 5.0)]
        trap = Candidate("trap", "TrapPhase", "local", "trap.int", ["Fe", "Ge"], None)
        trap.theory_peaks = [Peak("t1", 44.38, 100.0), Peak("t2", 28.0, 90.0), Peak("t3", 35.0, 75.0)]
        good_score = score_candidate(exp, good, tolerance_deg=0.2)
        trap_score = score_candidate(exp, trap, tolerance_deg=0.2)
        self.assertGreater(good_score.score, trap_score.score)
        self.assertGreaterEqual(len(trap_score.missing_strong_theory_peaks), 2)

    def test_vesta_mismatch_caps_score_and_reports_warning(self):
        exp = [Peak("e1", 20.0, 100.0), Peak("e2", 30.0, 80.0), Peak("e3", 40.0, 60.0)]
        candidate = Candidate("bad", "BadPhase", "library:local", "bad.cif", ["A"], None)
        candidate.theory_peaks = [Peak("t1", 20.0, 100.0), Peak("t2", 30.0, 80.0), Peak("t3", 40.0, 60.0)]
        candidate.simulation_validation = {
            "status": "failed",
            "reference_path": "BadPhase VESTA.int",
            "reason": "simulated peaks disagree with local VESTA reference",
        }

        score = score_candidate(exp, candidate, tolerance_deg=0.2)

        self.assertLessEqual(score.score, 0.29)
        self.assertEqual(score.score_components["confidence_label"], "untrusted")
        self.assertTrue(any("VESTA" in warning for warning in score.warnings))

    def test_no_vesta_reference_is_reported_as_model_only_not_high_confidence(self):
        exp = [Peak("e1", 20.0, 100.0), Peak("e2", 30.0, 80.0), Peak("e3", 40.0, 60.0), Peak("e4", 50.0, 40.0)]
        candidate = Candidate("model", "ModelOnly", "library:materials_project", "model.cif", ["A"], None)
        candidate.theory_peaks = [
            Peak("t1", 20.0, 100.0),
            Peak("t2", 30.0, 80.0),
            Peak("t3", 40.0, 60.0),
            Peak("t4", 50.0, 40.0),
        ]
        candidate.simulation_validation = {"status": "no_reference"}

        score = score_candidate(exp, candidate, tolerance_deg=0.2)

        self.assertEqual(score.score_components["confidence_label"], "model_only")
        self.assertTrue(any("VESTA" in warning for warning in score.warnings))

    def test_vesta_trust_check_detects_extra_strong_simulated_peak(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference_dir = root / "refs"
            reference_dir.mkdir()
            (reference_dir / "ZrFe6Ge4 VESTA.int").write_text("20 100\n30 80\n", encoding="utf-8")
            candidate = Candidate("bad", "ZrFe6Ge4", "library:local", "bad.cif", ["Zr", "Fe", "Ge"], None)
            candidate.theory_peaks = [Peak("t1", 20.0, 100.0), Peak("t2", 30.0, 80.0), Peak("t3", 70.0, 95.0)]
            with patch.dict("os.environ", {"OPHI_VESTA_REFERENCE_DIR": str(reference_dir)}):
                status = apply_vesta_trust_check(candidate, 10, 90)

        self.assertEqual(status["status"], "failed")
        self.assertGreaterEqual(status["extra_strong_count"], 1)
        self.assertEqual(candidate.simulation_validation["status"], "failed")

    def test_impurity_priority_targets_residual_peaks_after_strong_main_phase(self):
        exp = [
            Peak("e_main_1", 20.0, 100.0),
            Peak("e_main_2", 30.0, 90.0),
            Peak("e_extra", 50.0, 70.0),
        ]
        main = Candidate("main", "Main", "local", "main.int", ["A"], None)
        main.theory_peaks = [Peak("m1", 20.0, 100.0), Peak("m2", 30.0, 90.0)]
        high_total_overlap = Candidate("overlap", "Overlap", "local", "overlap.int", ["A"], None)
        high_total_overlap.theory_peaks = [Peak("o1", 20.0, 100.0), Peak("o2", 30.0, 80.0), Peak("o3", 70.0, 20.0)]
        residual_impurity = Candidate("residual", "ResidualImpurity", "local", "residual.int", ["B"], None)
        residual_impurity.theory_peaks = [Peak("r1", 50.0, 100.0)]

        scores = [score_candidate(exp, item, tolerance_deg=0.2) for item in [main, high_total_overlap, residual_impurity]]
        ordered = prioritize_impurities_after_main(scores, exp)

        self.assertEqual(ordered[0].candidate.formula_pretty, "Main")
        self.assertEqual(ordered[1].candidate.formula_pretty, "ResidualImpurity")

    def test_impurity_priority_rejects_dense_candidate_missing_most_of_its_strong_pattern(self):
        exp = [
            Peak("e_main_1", 20.0, 100.0),
            Peak("e_main_2", 30.0, 90.0),
            Peak("e_extra_1", 40.0, 70.0),
            Peak("e_extra_2", 50.0, 70.0),
            Peak("e_extra_3", 60.0, 70.0),
        ]
        main = Candidate("main", "Main", "local", "main.int", ["A"], None)
        main.theory_peaks = [Peak("m1", 20.0, 100.0), Peak("m2", 30.0, 90.0)]
        dense = Candidate("dense", "DenseFalsePositive", "local", "dense.int", ["B"], None)
        dense.theory_peaks = [
            Peak("d0", 40.0, 100.0),
            Peak("d1", 50.0, 100.0),
            Peak("d2", 60.0, 100.0),
        ] + [Peak(f"d{i}", 70.0 + i, 100.0) for i in range(3, 10)]
        supported = Candidate("supported", "SupportedImpurity", "local", "supported.int", ["B"], None)
        supported.theory_peaks = [Peak("s0", 50.0, 100.0)]

        scores = [score_candidate(exp, item, tolerance_deg=0.2) for item in [main, dense, supported]]
        ordered = prioritize_impurities_after_main(scores, exp, target_candidate_id="main")

        self.assertEqual(ordered[0].candidate.formula_pretty, "Main")
        self.assertEqual(ordered[1].candidate.formula_pretty, "SupportedImpurity")

    def test_main_phase_prefers_complete_main_element_system_over_single_element_tie(self):
        exp = [Peak("e1", 20.0, 100.0), Peak("e2", 30.0, 90.0)]
        single = Candidate("zr", "Zr", "library:materials_project", "zr.cif", ["Zr"], None)
        single.theory_peaks = [Peak("z1", 20.0, 100.0), Peak("z2", 30.0, 90.0)]
        ternary = Candidate("main", "ZrFe6Ge4", "library:local", "main.cif", ["Fe", "Ge", "Zr"], None)
        ternary.theory_peaks = [Peak("m1", 20.0, 100.0), Peak("m2", 30.0, 90.0)]

        scores = [score_candidate(exp, item, tolerance_deg=0.2) for item in [single, ternary]]
        ordered = prioritize_impurities_after_main(scores, exp, main_elements=["Zr", "Fe", "Ge"])

        self.assertEqual(ordered[0].candidate.formula_pretty, "ZrFe6Ge4")

    def test_randomized_scoring_penalizes_missing_strong_peaks(self):
        rng = random.Random(271828)
        for _ in range(20):
            base = sorted(rng.uniform(15, 80) for _ in range(5))
            exp = [Peak(f"e{i}", x + rng.uniform(-0.03, 0.03), rng.uniform(20, 100)) for i, x in enumerate(base)]
            strong = Candidate("strong", "Strong", "local", "", ["A"], None)
            strong.theory_peaks = [Peak(f"s{i}", x, rng.uniform(30, 100)) for i, x in enumerate(base[:4])]
            isolated = Candidate("isolated", "Isolated", "local", "", ["A"], None)
            isolated.theory_peaks = [Peak("i0", base[0], 100.0)] + [
                Peak(f"i{j}", rng.uniform(82, 89), rng.uniform(40, 100)) for j in range(1, 5)
            ]
            self.assertGreater(
                score_candidate(exp, strong).score,
                score_candidate(exp, isolated).score,
            )

    def test_multiphase_and_reports_have_expected_shape(self):
        exp = [Peak("e1", 20.0, 100.0), Peak("e2", 30.0, 80.0), Peak("e3", 50.0, 40.0)]
        c1 = Candidate("c1", "A", "local", "a.int", ["A"], None)
        c1.theory_peaks = [Peak("t1", 20.0, 100.0), Peak("t2", 30.0, 50.0)]
        c2 = Candidate("c2", "B", "local", "b.int", ["B"], None)
        c2.theory_peaks = [Peak("t1", 50.0, 100.0)]
        explanation = explain_multiphase(exp, [score_candidate(exp, c1), score_candidate(exp, c2)], max_phases=2)
        with tempfile.TemporaryDirectory() as tmp:
            outputs = write_reports(Path(tmp), {"xrd_file": "synthetic.xy"}, exp, [score_candidate(exp, c1)], explanation)
            payload = json.loads(Path(outputs["json"]).read_text(encoding="utf-8"))
            report = Path(outputs["markdown"]).read_text(encoding="utf-8")
        self.assertIn("top_candidates", payload)
        self.assertIn("candidate screening only", report)
        self.assertTrue(explanation.selected_candidates)


if __name__ == "__main__":
    unittest.main()
