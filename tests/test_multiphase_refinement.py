from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock

import numpy as np

from ophiuchus.xrd.multiphase_models import MultiphaseRefinementSettings, PhaseRefinementInput
from ophiuchus.xrd.multiphase_refinement import (
    RietanMultiphaseBackend,
    parse_multiphase_output,
    patch_multiphase_refinement_input,
)


COMBINED_INS = """NMODE = 0: Rietveld analysis.
NMODE = 1! Simulation.
NPRINT = 1
NINT = 0! RIETAN format.
NINT = 1: General format.
NRANGE = 0: Refine background.
NRANGE = 2! Fixed background.
NPRFN = 0
NASYM = 1
SHIFT0 0.0 0.0 0.0 0.0 1110
BACKGROUND 0 0 0 0 0 0 0 0 0 0 0 0 111111111111
! Parameters @1
SCALE@1 1.0 1
GAUSS01@1 0.01 -0.006 0.003 0.0 1110
LORENTZ01@1 0.04 0.0 -0.04 0.0 1010
ASYM01@1 0.03 0 0 0 0 0 100000
PREF@1 1 1 0 0 0 0 000000
CELLQ@1 5 5 8 90 90 120 0 1111110
Zr1@1/Zr 1 0 0 0 0.5 01111
# End Parameters @1
! Parameters @2
SCALE@2 1.0 1
GAUSS01@2 0.01 -0.006 0.003 0.0 1110
LORENTZ01@2 0.04 0.0 -0.04 0.0 1010
ASYM01@2 0.03 0 0 0 0 0 100000
PREF@2 1 1 0 0 0 0 000000
CELLQ@2 4 4 6 90 90 90 0 1111110
Sn1@2/Sn 1 0 0 0 0.5 01111
# End Parameters @2
NAUTO = 0! simultaneous
NAUTO = 2: automatic incremental
NUPDT = 0
NPAT = 1
"""

GPD = """# 2-theta, observed intensity, calculated intensity, residual, background
10.0 12.0 11.0 1.0 3.0
10.1 20.0 18.0 2.0 3.1
10.2 13.0 13.5 -0.5 3.2

# h, k, l, 2-theta (shifted), 2-theta (not shifted), d, sin(theta), FWHM*cos(theta), and beta*cos(theta) for phase No. 1
1 0 0 10.05 10.06 8.0 0.1 0.002 0.003

# h, k, l, 2-theta (shifted), 2-theta (not shifted), d, sin(theta), FWHM*cos(theta), and beta*cos(theta) for phase No. 2
0 1 0 10.15 10.16 7.0 0.1 0.002 0.003
"""

LST = """Phase #1: Main
Phase #2: Impurity
Cycle # 9
Rwp = 9.0 Rp = 7.0 RR = 8.0 Re = 5.0 S = 1.8 GofF = 3.2
Cycle # 10
Rwp = 8.5 Rp = 6.0 RR = 7.7 Re = 5.0 S = 1.2 GofF = 1.44

*** Final parameters and their estimated standard uncertainties ***
 No. A SIGMA DELTA.A/SIGMA
 21 2.000000E+00 1.0E-3 0.1 Scale factor, s
 80 1.000000E+00 1.0E-3 0.1 Scale factor, s

Effective radii (R), particle absorption factors (tau), and mass/mole fractions uncorrected and corrected for microabsorption
 Phase R mu/rho mu [mu-mu(mean)]R tau w X w(cor) X(cor)
 Main 5.0 10.0 20.0 0.0 1.0 0.94117647 0.8 0.94 0.8
 Impurity 5.0 10.0 20.0 0.0 1.0 0.05882353 0.2 0.06 0.2
"""


class MultiphaseRefinementTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.main_cif = self.root / "main.cif"
        self.impurity_cif = self.root / "impurity.cif"
        self.main_cif.write_text("data_main\n", encoding="ascii")
        self.impurity_cif.write_text("data_impurity\n", encoding="ascii")
        self.phases = (
            PhaseRefinementInput("main", "Main", self.main_cif, "target", initial_scale=2.0),
            PhaseRefinementInput("imp", "Impurity", self.impurity_cif, "impurity", initial_scale=1.0),
        )

    def tearDown(self) -> None:
        self.temp.cleanup()

    def test_patch_first_stage_refines_only_background_and_all_scales(self):
        patched = patch_multiphase_refinement_input(
            COMBINED_INS,
            self.phases,
            MultiphaseRefinementSettings(10.0, 90.0, background_terms=4, refine_zero_shift=False, refine_profile=False),
        )

        self.assertIn("NMODE = 0:", patched)
        self.assertIn("NINT = 1:", patched)
        self.assertIn("BACKGROUND 0 0 0 0 0 0 0 0 0 0 0 0 111100000000", patched)
        self.assertIn("SCALE@1 2 1", patched)
        self.assertIn("SCALE@2 1 1", patched)
        self.assertIn("SHIFT0 0.0 0.0 0.0 0.0 0000", patched)
        self.assertIn("GAUSS01@1 0.01 -0.006 0.003 0.0 0000", patched)
        self.assertIn("CELLQ@1 5 5 8 90 90 120 0 0000000", patched)
        self.assertIn("Zr1@1/Zr 1 0 0 0 0.5 00000", patched)
        self.assertIn("Sn1@2/Sn 1 0 0 0 0.5 00000", patched)

    def test_patch_removes_stale_atomic_constraints_but_keeps_profile_constraints(self):
        combined = COMBINED_INS + """
A(O2@1,B)=A(O1@1,B); A(O3@1,B)=A(O1@1,B)
A(O2@1,X)=A(O1@1,X)
A(GAUSS01@2,1)=A(GAUSS01@1,1)
"""

        patched = patch_multiphase_refinement_input(
            combined,
            self.phases,
            MultiphaseRefinementSettings(10.0, 90.0),
        )

        self.assertNotIn("A(O2@1,B)", patched)
        self.assertNotIn("A(O2@1,X)", patched)
        self.assertIn("A(GAUSS01@2,1)=A(GAUSS01@1,1)", patched)

    def test_parser_assigns_scales_native_weights_and_reflections_by_phase_order(self):
        result = parse_multiphase_output(
            GPD,
            LST,
            self.phases,
            [(4.0, 100.0, 50.0), (2.0, 50.0, 25.0)],
        )

        self.assertAlmostEqual(result.rwp_percent, 8.5)
        self.assertAlmostEqual(result.phases[0].scale, 2.0)
        self.assertAlmostEqual(result.phases[1].scale, 1.0)
        self.assertAlmostEqual(result.phases[0].weight_percent, 94.117647)
        self.assertAlmostEqual(result.phases[1].weight_percent, 5.882353)
        self.assertEqual(result.phases[0].reflection_two_theta_deg, (10.05,))
        self.assertEqual(result.phases[1].reflection_two_theta_deg, (10.15,))

    def test_backend_runs_one_rietan_job_and_records_engine_provenance(self):
        rietan = self.root / "RIETAN.exe"
        cif2ins = self.root / "cif2ins.exe"
        template = self.root / "template.ins"
        rietan.write_bytes(b"")
        cif2ins.write_bytes(b"")
        template.write_text("official template", encoding="ascii")
        job_dirs: list[Path] = []

        def fake_export(phase, phase_number, work, *_args, **_kwargs):
            path = Path(work) / f"phase@{phase_number}.ins"
            path.write_text(f"phase {phase_number}", encoding="ascii")
            return path

        def fake_combine(_template, _phase_texts):
            return COMBINED_INS

        def fake_run(command, **kwargs):
            work = Path(kwargs["cwd"])
            job_dirs.append(work)
            self.assertEqual(Path(command[0]).name.lower(), "rietan.exe")
            self.assertTrue((work / "multi_phase.int").is_file())
            (work / "multi_phase.gpd").write_text(GPD, encoding="ascii")
            return subprocess.CompletedProcess(command, 0, stdout=LST, stderr="")

        backend = RietanMultiphaseBackend(
            rietan_exe=rietan,
            cif2ins_exe=cif2ins,
            template_path=template,
        )
        with mock.patch("ophiuchus.xrd.multiphase_refinement.export_phase_input", side_effect=fake_export), mock.patch(
            "ophiuchus.xrd.multiphase_refinement.combine_phase_inputs", side_effect=fake_combine
        ), mock.patch("ophiuchus.xrd.multiphase_refinement.subprocess.run", side_effect=fake_run), mock.patch(
            "ophiuchus.xrd.multiphase_refinement.structure_zmv",
            side_effect=[(4.0, 100.0, 50.0), (2.0, 50.0, 25.0)],
        ):
            result = backend.refine(
                self.phases,
                np.array([10.0, 10.1, 10.2]),
                np.array([12.0, 20.0, 13.0]),
                MultiphaseRefinementSettings(10.0, 10.2),
            )

        self.assertEqual(result.provenance["profile_engine"], "RIETAN-FP multiphase refinement")
        self.assertEqual(result.provenance["phase_count"], 2)
        self.assertTrue(job_dirs)
        self.assertFalse(job_dirs[0].exists())


if __name__ == "__main__":
    unittest.main()
