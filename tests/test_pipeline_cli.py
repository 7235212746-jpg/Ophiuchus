import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from ophiuchus.xrd.cache import CandidateCache
from ophiuchus.xrd.pipeline import run_analysis


class PipelineCliTests(unittest.TestCase):
    def _make_case(self, root: Path) -> tuple[Path, Path, Path]:
        xrd = root / "sample.xy"
        rows = []
        for i in range(2000):
            x = 10.0 + i * 0.03
            y = 1.0
            for center, amp in [(20.0, 100), (30.0, 70), (44.4, 35)]:
                y += amp * pow(2.718281828, -((x - center) ** 2) / (2 * 0.06**2))
            rows.append(f"{x:.3f} {y:.5f}")
        xrd.write_text("\n".join(rows), encoding="utf-8")
        cands = root / "candidates"
        cands.mkdir()
        (cands / "ZrFe6Ge4.int").write_text("20.01 100\n30.02 70\n44.39 20\n", encoding="utf-8")
        (cands / "Fe2O3.int").write_text("44.39 100\n65.0 80\n70.0 60\n", encoding="utf-8")
        out = root / "results"
        return xrd, cands, out

    def test_run_analysis_writes_reports_from_generated_case(self):
        with tempfile.TemporaryDirectory() as tmp:
            xrd, cands, out = self._make_case(Path(tmp))
            result = run_analysis(
                xrd_file=xrd,
                candidate_dirs=[cands],
                elements=["Zr", "Fe", "Ge"],
                out_dir=out,
                max_candidates=5,
            )
            payload = json.loads(Path(result.outputs["json"]).read_text(encoding="utf-8"))
            audit = Path(result.outputs["audit_folder"])
            self.assertGreaterEqual(len(result.experimental_peaks), 3)
            self.assertEqual(payload["top_candidates"][0]["formula"], "ZrFe6Ge4")
            self.assertIn("possible_impurity_phases", payload)
            self.assertIn("per_peak_candidate_assignments", payload)
            self.assertTrue((audit / "candidate_phase_ranking.csv").exists())
            self.assertTrue((audit / "peak_match_table.csv").exists())
            self.assertIsNotNone(result.context)
            self.assertAlmostEqual(float(np.max(result.context.intensity)), 99.62071, places=4)
            self.assertEqual(result.context.source_path, str(xrd))
            self.assertTrue(result.context.source_fingerprint)
            self.assertTrue(result.context.data_fingerprint)

    def test_cli_analyze_command_writes_json_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            xrd, cands, out = self._make_case(Path(tmp))
            proc = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "ophiuchus",
                    "analyze",
                    "--xrd",
                    str(xrd),
                    "--cif-dir",
                    str(cands),
                    "--elements",
                    "Zr",
                    "Fe",
                    "Ge",
                    "--out-dir",
                    str(out),
                ],
                cwd=str(Path(__file__).parents[1]),
                text=True,
                capture_output=True,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertTrue((out / "results.json").exists())
            self.assertIn("report.md", proc.stdout)

    def test_run_analysis_uses_cache_and_writes_plot_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            xrd, cands, out = self._make_case(root)
            cache_path = root / "ophi_cache.sqlite"
            result = run_analysis(
                xrd_file=xrd,
                candidate_dirs=[cands],
                elements=["Zr", "Fe", "Ge"],
                out_dir=out,
                max_candidates=5,
                cache_path=cache_path,
            )
            cache_stats = CandidateCache(cache_path).stats()
            payload = json.loads(Path(result.outputs["json"]).read_text(encoding="utf-8"))
            self.assertGreaterEqual(cache_stats["patterns"], 1)
            self.assertTrue(Path(result.outputs["xrd_plot"]).exists())
            self.assertTrue(Path(result.outputs["xrd_plot_clean"]).exists())
            self.assertTrue(Path(result.outputs["contribution_proxy"]).exists())
            self.assertIn("relative_contribution_proxy", payload)
            self.assertIn("auto_range", payload["input"])

    def test_cli_cache_build_command_creates_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, cands, _ = self._make_case(root)
            cache_path = root / "cache.sqlite"
            proc = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "ophiuchus",
                    "cache-build",
                    "--cif-dir",
                    str(cands),
                    "--elements",
                    "Zr",
                    "Fe",
                    "Ge",
                    "--cache-path",
                    str(cache_path),
                ],
                cwd=str(Path(__file__).parents[1]),
                text=True,
                capture_output=True,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertGreaterEqual(CandidateCache(cache_path).stats()["patterns"], 1)
            self.assertIn("stored_patterns", proc.stdout)


if __name__ == "__main__":
    unittest.main()
