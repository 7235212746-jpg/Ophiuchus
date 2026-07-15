import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from ophiuchus.xrd.cache import CandidateCache, build_candidate_cache
from ophiuchus.xrd.contribution import contribution_proxy
from ophiuchus.xrd.models import Candidate, CandidateScore, MultiPhaseExplanation, Peak, PeakMatch, PatternPoint, XrdPattern
from ophiuchus.xrd.phase_context import phase_context_cards
from ophiuchus.xrd.plotting import select_dashboard_plot_scores, write_xrd_plot
from ophiuchus.xrd.audit import write_analysis_audit_folder


class CachePlotContributionTests(unittest.TestCase):
    def test_cache_roundtrip_candidate_and_peaks(self):
        with tempfile.TemporaryDirectory() as tmp:
            cache = CandidateCache(Path(tmp) / "cache.sqlite")
            candidate = Candidate("c1", "ZrFe6Ge4", "local_peak_list", "phase.int", ["Zr", "Fe", "Ge"], "hash")
            peaks = [Peak("p1", 20.0, 100.0), Peak("p2", 30.0, 50.0)]
            cache.upsert_candidate(candidate)
            cache.store_pattern(candidate, peaks, radiation="CuKa", two_theta_range=(10, 90))
            loaded = cache.load_pattern(candidate, radiation="CuKa", two_theta_range=(10, 90))
        self.assertEqual([(p.two_theta, p.intensity) for p in loaded], [(20.0, 100.0), (30.0, 50.0)])

    def test_build_candidate_cache_from_local_peak_lists(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "candidates"
            source.mkdir()
            (source / "ZrFe6Ge4.int").write_text("20 100\n30 60\n", encoding="utf-8")
            cache = CandidateCache(root / "cache.sqlite")
            summary = build_candidate_cache(cache, [source], {"Zr", "Fe", "Ge"}, radiation="CuKa", two_theta_range=(10, 90))
            candidates = cache.list_candidates()
        self.assertEqual(summary["stored_patterns"], 1)
        self.assertEqual(candidates[0].formula_pretty, "ZrFe6Ge4")

    def test_plot_png_is_written(self):
        with tempfile.TemporaryDirectory() as tmp:
            pattern = XrdPattern([PatternPoint(10 + i, float(i % 10 + 1)) for i in range(50)])
            peaks = [Peak("e1", 20.0, 100.0), Peak("e2", 30.0, 80.0)]
            candidate = Candidate("c1", "PhaseA", "local", "", ["A"], None)
            candidate.theory_peaks = [Peak("t1", 20.1, 100.0), Peak("t2", 32.0, 40.0)]
            score = CandidateScore(candidate, 0.8, [PeakMatch(candidate.theory_peaks[0], peaks[0], 0.1)], [], [peaks[0]], [peaks[1]])
            out = write_xrd_plot(Path(tmp) / "plot.png", pattern, peaks, [score])
            self.assertTrue(Path(out).exists())
            self.assertGreater(Path(out).stat().st_size, 1000)

    def test_plotting_cli_renders_from_analysis_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = {
                "input": {"two_theta_range": [10, 50]},
                "experimental_peaks": [
                    {"peak_id": "e1", "two_theta": 20.0, "intensity": 100.0},
                    {"peak_id": "e2", "two_theta": 30.0, "intensity": 70.0},
                ],
                "top_candidates": [
                    {
                        "rank": 1,
                        "candidate_id": "c1",
                        "formula": "PhaseA",
                        "score": 0.8,
                        "source": "test",
                        "matched_peaks": [{"theory_two_theta": 20.1, "theory_intensity": 100.0}],
                        "missing_strong_peaks": [],
                    }
                ],
            }
            analysis = root / "analysis.json"
            out = root / "plot.png"
            analysis.write_text(json.dumps(payload), encoding="utf-8")
            proc = subprocess.run(
                [sys.executable, "-m", "ophiuchus.plotting.render_xrd", "--analysis-json", str(analysis), "--out", str(out)],
                cwd=str(Path(__file__).parents[1]),
                text=True,
                capture_output=True,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue(out.exists())
            self.assertGreater(out.stat().st_size, 1000)

    def test_dashboard_plot_selects_target_plus_three_impurities(self):
        scores = []
        for idx, formula in enumerate(["Target", "ImpurityA", "ImpurityB", "ImpurityC", "ImpurityD"]):
            candidate = Candidate(f"c{idx}", formula, "local", "", ["A"], None)
            scores.append(CandidateScore(candidate, 1.0 - idx * 0.1, [], [], [], []))
        selected = select_dashboard_plot_scores(scores)
        self.assertEqual([item[0] for item in selected], ["目标模拟峰", "可能杂质 1", "可能杂质 2", "可能杂质 3"])
        self.assertEqual([item[1].candidate.formula_pretty for item in selected], ["Target", "ImpurityA", "ImpurityB", "ImpurityC"])

    def test_contribution_proxy_is_normalized_and_warned(self):
        c1 = Candidate("c1", "Main", "local", "", ["A"], None)
        c2 = Candidate("c2", "Impurity", "local", "", ["A"], None)
        e1 = Peak("e1", 20, 100)
        e2 = Peak("e2", 30, 50)
        s1 = CandidateScore(c1, 0.8, [], [], [e1], [])
        s2 = CandidateScore(c2, 0.4, [], [], [e2], [])
        explanation = MultiPhaseExplanation([s1, s2], [e1, e2], [])
        proxy = contribution_proxy(explanation)
        self.assertAlmostEqual(sum(proxy.contributions.values()), 100.0)
        self.assertIn("not Rietveld", proxy.warning)

    def test_phase_context_cards_are_conservative(self):
        candidate = Candidate("c1", "Fe1.67Ge", "local", "", ["Fe", "Ge"], None)
        cards = phase_context_cards([candidate], synthesis_metadata={"temperature": "900 C", "duration": "72 h"})
        self.assertEqual(cards[0]["chemical_system"], "Fe-Ge")
        self.assertEqual(cards[0]["data_status"], "insufficient data")
        self.assertTrue(cards[0]["suggestions"])

    def test_audit_exports_simulation_validation_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate = Candidate("c1", "ZrFe6Ge4", "library:local", "phase.cif", ["Zr", "Fe", "Ge"], None)
            candidate.simulation_validation = {
                "status": "failed",
                "reference_path": "ZrFe6Ge4 VESTA.int",
                "missing_strong_count": 1,
                "extra_strong_count": 2,
                "max_strong_intensity_ratio_error": 3.2,
                "reason": "simulated peaks disagree with local VESTA reference",
            }
            score = CandidateScore(candidate, 0.2, [], [], [], [], score_components={"confidence_label": "untrusted"})
            audit = Path(write_analysis_audit_folder(root, {"xrd_file": "sample.xy"}, [], [score], [candidate], {}))
            validation_csv = audit / "simulation_validation.csv"
            ranking_csv = audit / "candidate_phase_ranking.csv"
            self.assertTrue(validation_csv.exists())
            validation_text = validation_csv.read_text(encoding="utf-8")
            ranking_text = ranking_csv.read_text(encoding="utf-8")
        self.assertIn("failed", validation_text)
        self.assertIn("ZrFe6Ge4 VESTA.int", validation_text)
        self.assertIn("simulation_validation_status", ranking_text)


if __name__ == "__main__":
    unittest.main()
