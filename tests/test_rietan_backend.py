import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock
from types import SimpleNamespace

from ophiuchus.xrd.backend import SimulationContext
from ophiuchus.xrd.config import XRDConfig
from ophiuchus.xrd.rietan_backend import (
    RietanXRDBackend,
    build_rietan_command,
    parse_rietan_profile,
    parse_rietan_reflections,
    patch_rietan_input,
    profile_peaks,
    upgrade_scores_with_rietan,
)
from ophiuchus.xrd.models import Candidate, CandidateScore, Peak, PatternPoint, XrdPattern


INS_TEXT = """! analytical mode\r
NMODE = 0: Rietveld analysis.\r
NMODE = 1! Simulation.\r
NPRINT = 0! Minimal.\r
NPRINT = 1! Standard.\r
NPRINT = 1\r
   DEG1 = 3.0: Minimum.\r
   DEG2 = 120.0: Maximum.\r
NPAT = 1: Gnuplot.\r
"""

LST_TEXT = """
            No. Phase     h    k    l  Code     2-theta          d     Ical |F(crys)| |F(magn)|     POF      FWHM    m       Dd/d
             1      1     1    0    0    +1      20.000    4.00000       25   10.0000      -      1.000   0.07500    2   0.010000
             2      1     2    0    0    +1      40.000    2.00000      100   20.0000      -      1.000   0.08000    4   0.005000

       *** End of job ***
           --- RIETAN-FP v3.12  Copyright 2009-2021 by Fujio Izumi ---
"""

GPD_TEXT = """# 2-theta and calculated intensity
10.0 0.0
20.0 1.0
30.0 0.0
40.0 4.0
50.0 0.0

# reflection section
1 0 0 20.0 20.0 4.0
"""


class RietanBackendTests(unittest.TestCase):
    def test_patch_switches_only_active_mode_and_requested_range(self):
        patched = patch_rietan_input(INS_TEXT, 10.0, 90.0)

        self.assertIn("NMODE = 0! Rietveld analysis.", patched)
        self.assertIn("NMODE = 1: Simulation.", patched)
        self.assertIn("DEG1 = 10.000000:", patched)
        self.assertIn("DEG2 = 90.000000:", patched)
        self.assertIn("NPAT = 1:", patched)
        self.assertIn("NPRINT = 1\r\n", patched)
        self.assertNotIn("NPRINT = 1:", patched)
        self.assertNotIn("NMODE = 0:", patched)

    def test_parses_profile_and_uses_rietan_profile_maxima_with_hkl(self):
        profile = parse_rietan_profile(GPD_TEXT)
        reflections = parse_rietan_reflections(LST_TEXT)
        peaks = profile_peaks(profile, reflections, intensity_threshold=1.0)

        self.assertEqual(profile.two_theta_deg, (10.0, 20.0, 30.0, 40.0, 50.0))
        self.assertEqual(profile.normalized_intensity, (0.0, 25.0, 0.0, 100.0, 0.0))
        self.assertEqual([peak.two_theta for peak in peaks], [20.0, 40.0])
        self.assertEqual([peak.intensity for peak in peaks], [25.0, 100.0])
        self.assertEqual([peak.hkl for peak in peaks], ["(1 0 0)", "(2 0 0)"])
        self.assertEqual([peak.multiplicity for peak in peaks], [2, 4])

    def test_command_uses_relative_sample_files_required_by_rietan(self):
        command = build_rietan_command(Path("C:/Tools/RIETAN.exe"), "ophi")

        self.assertEqual(command[0], str(Path("C:/Tools/RIETAN.exe")))
        self.assertEqual(command[1], "ophi.ins")
        self.assertEqual(command[-1], "ophi.exp")
        self.assertEqual(len(command), 19)

    def test_simulate_runs_vesta_then_rietan_and_preserves_profile(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cif = root / "phase.cif"
            cif.write_text("data_phase\n", encoding="utf-8")
            vesta = root / "VESTA.exe"
            rietan = root / "RIETAN.exe"
            vesta.write_bytes(b"")
            rietan.write_bytes(b"")

            def fake_run(command, **kwargs):
                work = Path(kwargs["cwd"])
                if Path(command[0]).name.lower() == "vesta.exe":
                    (work / "ophi.ins").write_text(INS_TEXT, encoding="utf-8", newline="")
                    # VESTA CLI may report a non-zero status even after a successful save.
                    return subprocess.CompletedProcess(command, 1, stdout="Saved data to: ophi.ins", stderr="")
                (work / "ophi.gpd").write_text(GPD_TEXT, encoding="ascii", newline="")
                return subprocess.CompletedProcess(command, 0, stdout=LST_TEXT, stderr="")

            backend = RietanXRDBackend(vesta_exe=vesta, rietan_exe=rietan)
            with mock.patch("ophiuchus.xrd.rietan_backend.subprocess.run", side_effect=fake_run) as run:
                pattern = backend.simulate_cif(
                    cif,
                    XRDConfig(two_theta_min=10.0, two_theta_max=50.0),
                    SimulationContext("phase", "local", "FeGe"),
                )

        self.assertEqual(run.call_count, 2)
        self.assertEqual(run.call_args_list[0].args[0][1:4], ["-nogui", "-i", str(cif.resolve())])
        self.assertEqual(run.call_args_list[1].args[0][1], "ophi.ins")
        self.assertTrue(run.call_args_list[1].kwargs["env"]["RIETAN"].endswith(os.sep))
        self.assertEqual(pattern.engine_name, "VESTA / RIETAN-FP")
        self.assertEqual(pattern.engine_version, "RIETAN-FP 3.12")
        self.assertEqual(pattern.two_theta_deg, (20.0, 40.0))
        self.assertEqual(pattern.normalized_intensity, (25.0, 100.0))
        self.assertEqual(pattern.profile_two_theta_deg, (10.0, 20.0, 30.0, 40.0, 50.0))
        self.assertEqual(pattern.profile_normalized_intensity, (0.0, 25.0, 0.0, 100.0, 0.0))

    def test_upgrade_replaces_display_candidate_peaks_and_rescores(self):
        with tempfile.TemporaryDirectory() as tmp:
            cif = Path(tmp) / "phase.cif"
            cif.write_text("data_phase\n", encoding="utf-8")
            candidate = Candidate("phase", "FeGe", "local", str(cif), ["Fe", "Ge"], "hash")
            candidate.theory_peaks = [Peak("old", 25.0, 100.0)]
            original = CandidateScore(candidate, 0.0, [], [], [], [])
            exact_peak = Peak("exact", 40.0, 100.0, hkl="(2 0 0)")
            pattern = SimpleNamespace(
                to_peaks=lambda: [exact_peak],
                engine_name="VESTA / RIETAN-FP",
                engine_version="RIETAN-FP 3.12",
                pattern_fingerprint="pattern",
                cif_sha256="cif",
            )
            backend = SimpleNamespace(simulate_cif=mock.Mock(return_value=pattern))

            upgraded, warnings = upgrade_scores_with_rietan(
                [original],
                backend,
                XRDConfig(two_theta_min=10.0, two_theta_max=50.0),
                [Peak("exp", 40.01, 100.0)],
                tolerance_deg=0.2,
                limit=1,
            )

        self.assertEqual(warnings, [])
        self.assertEqual(upgraded[0].candidate.theory_peaks, [exact_peak])
        self.assertEqual(len(upgraded[0].matched_theory_peaks), 1)
        self.assertEqual(upgraded[0].candidate.parse_status, "vesta_rietan_simulated")
        self.assertEqual(upgraded[0].candidate.simulation_validation["engine_name"], "VESTA / RIETAN-FP")

    def test_dashboard_uses_exact_rietan_profile_instead_of_rebroadening_sticks(self):
        from ophiuchus.xrd.plotting import build_xrd_figure

        profile_x = (10.0, 20.0, 30.0, 40.0, 50.0)
        profile_y = (0.0, 25.0, 0.0, 100.0, 0.0)
        candidate = Candidate("phase", "FeGe", "local", "phase.cif", ["Fe", "Ge"], "hash")
        candidate.theory_peaks = [Peak("p1", 20.0, 25.0), Peak("p2", 40.0, 100.0)]
        candidate.simulated_pattern = SimpleNamespace(
            profile_two_theta_deg=profile_x,
            profile_normalized_intensity=profile_y,
        )
        score = CandidateScore(candidate, 1.0, [], [], [], [])
        experimental = XrdPattern([PatternPoint(10.0, 0.0), PatternPoint(50.0, 1.0)])

        figure = build_xrd_figure(experimental, [], [score])
        try:
            target_line = next(
                line for line in figure.axes[1].lines if str(line.get_label()).startswith("Target simulated")
            )
            self.assertEqual(tuple(float(value) for value in target_line.get_xdata()), profile_x)
            self.assertEqual(tuple(float(value) for value in target_line.get_ydata()), profile_y)
        finally:
            import matplotlib.pyplot as plt

            plt.close(figure)


if __name__ == "__main__":
    unittest.main()
