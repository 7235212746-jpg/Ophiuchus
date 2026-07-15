from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock

import numpy as np

from ophiuchus.xrd.refinement import (
    RefinementSettings,
    RietanRefinementBackend,
    parse_refinement_output,
    patch_refinement_input,
    refinement_trust_warnings,
    write_xy_intensity,
)


REFINEMENT_INS = """NMODE = 0: Rietveld analysis.\r
NMODE = 1! Simulation.\r
NPRINT = 1\r
NINT = 0! RIETAN format.\r
NINT = 1: General format.\r
NRANGE = 0: Refine background.\r
NRANGE = 2! Fixed background.\r
NPRFN = 1\r
SHIFTN  0.0 0.0 0.0 0.0  1000\r
BKGD  0 0 0 0 0 0 0 0 0 0 0 0  111111111100\r
SCALE  3.5E-5  1\r
FWHM12  0.005 -0.001 0.005 0.0  1110\r
ASYM12  1.0 0.1 -0.04 0.0  1110\r
ETA12  0.6 0.1 0.5 0.1  1111\r
PREF  1.0 1.0 0.0 0.0 0.0 0.0  000000\r
CELLQ  9.0 9.0 5.4 90 90 120 0.0  1111110\r
Sn1/Sn  1.0 0.0 0.0 0.0 0.5  01111\r
} End of lines for label/species, A(I), and ID(I)\r
NAUTO = 0! simultaneous\r
NAUTO = 2: automatic incremental\r
NPAT = 1: Gnuplot\r
"""

REFINEMENT_GPD = """# 2-theta, observed intensity, calculated intensity, residual, background
10.0 12.0 11.0 1.0 3.0
10.1 20.0 18.0 2.0 3.1
10.2 13.0 13.5 -0.5 3.2

# reflections
1 0 0 10.1 10.1 2.0
"""

REFINEMENT_LST = """
Cycle # 9
Rwp =  9.100    Rp =  7.200    RR = 8.0    Re = 5.0    S = 1.8200    GofF = 3.3124
Cycle # 10
Rwp =  8.568    Rp =  6.038    RR = 7.7    Re = 10.268    S = 0.8345    GofF = 0.6964

*** Final parameters and their estimated standard uncertainties ***
1  1.127494E-02  1.0E-05  1.0 Peak-shift parameter, t0
9  21.0957  0.10  0.07 Background parameter, b0
21  9.986060E-06  2.2E-08  0.01 Scale factor, s
22  1.102781E-02  8.4E-06  5.9 FWHM parameter, U
23  1.723084E-03  8.1E-06  3.2 FWHM parameter, V
24  4.850867E-03  1.7E-06  1.7 FWHM parameter, W
58  8.87529  1.0E-04  2.0 Lattice parameter, a
59  5.45213  1.0E-04  2.0 Lattice parameter, b
60  7.15271  1.0E-04  2.0 Lattice parameter, c
61  90.0000  1.0E-04  2.0 Lattice parameter, alpha
62  90.0000  1.0E-04  2.0 Lattice parameter, beta
63  120.0000  1.0E-04  2.0 Lattice parameter, gamma
"""


class RietanRefinementTests(unittest.TestCase):
    def test_writes_general_xy_data_without_normalizing_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "ophi.int"
            write_xy_intensity(path, np.array([10.0, 10.1]), np.array([120.0, 85.5]))

            self.assertEqual(
                path.read_text(encoding="ascii"),
                "GENERAL\n2\n10.000000 120.000000\n10.100000 85.500000\n",
            )

    def test_patches_a_conservative_single_phase_refinement(self):
        patched = patch_refinement_input(
            REFINEMENT_INS,
            RefinementSettings(two_theta_min=10.0, two_theta_max=90.0, background_terms=6),
        )

        self.assertIn("NMODE = 0:", patched)
        self.assertIn("NMODE = 1!", patched)
        self.assertIn("NINT = 1:", patched)
        self.assertIn("NRANGE = 0:", patched)
        self.assertIn("SHIFTN  0.0 0.0 0.0 0.0  1000", patched)
        self.assertIn("BKGD  0 0 0 0 0 0 0 0 0 0 0 0  111111000000", patched)
        self.assertIn("FWHM12  0.005 -0.001 0.005 0.0  1110", patched)
        self.assertIn("ASYM12  1.0 0.1 -0.04 0.0  0000", patched)
        self.assertIn("ETA12  0.6 0.1 0.5 0.1  0000", patched)
        self.assertIn("CELLQ  9.0 9.0 5.4 90 90 120 0.0  0000000", patched)
        self.assertIn("Sn1/Sn  1.0 0.0 0.0 0.0 0.5  00000", patched)
        self.assertIn("NAUTO = 2:", patched)

    def test_optional_lattice_refinement_never_unlocks_structure_parameters(self):
        patched = patch_refinement_input(
            REFINEMENT_INS,
            RefinementSettings(
                two_theta_min=10.0,
                two_theta_max=90.0,
                background_terms=6,
                refine_lattice=True,
            ),
        )

        self.assertIn("CELLQ  9.0 9.0 5.4 90 90 120 0.0  1111110", patched)
        self.assertIn("Sn1/Sn  1.0 0.0 0.0 0.0 0.5  00000", patched)

    def test_first_stage_can_fit_only_background_and_scale(self):
        patched = patch_refinement_input(
            REFINEMENT_INS,
            RefinementSettings(
                two_theta_min=10.0,
                two_theta_max=90.0,
                background_terms=4,
                refine_zero_shift=False,
                refine_profile=False,
            ),
        )

        self.assertIn("SHIFTN  0.0 0.0 0.0 0.0  0000", patched)
        self.assertIn("BKGD  0 0 0 0 0 0 0 0 0 0 0 0  111100000000", patched)
        self.assertIn("SCALE  3.5E-5  1", patched)
        self.assertIn("FWHM12  0.005 -0.001 0.005 0.0  0000", patched)

    def test_parses_profile_metrics_parameters_and_keeps_arrays_immutable(self):
        result = parse_refinement_output(REFINEMENT_GPD, REFINEMENT_LST)

        np.testing.assert_allclose(result.two_theta_deg, [10.0, 10.1, 10.2])
        np.testing.assert_allclose(result.observed_intensity, [12.0, 20.0, 13.0])
        np.testing.assert_allclose(result.calculated_intensity, [11.0, 18.0, 13.5])
        np.testing.assert_allclose(result.residual_intensity, [1.0, 2.0, -0.5])
        np.testing.assert_allclose(result.background_intensity, [3.0, 3.1, 3.2])
        self.assertEqual(result.reflection_two_theta_deg, (10.1,))
        self.assertAlmostEqual(result.rwp_percent, 8.568)
        self.assertAlmostEqual(result.rp_percent, 6.038)
        self.assertAlmostEqual(result.s_value, 0.8345)
        self.assertAlmostEqual(result.goodness_of_fit, 0.6964)
        self.assertAlmostEqual(result.parameters["zero_shift"], 0.01127494)
        self.assertAlmostEqual(result.parameters["scale"], 9.986060e-6)
        self.assertAlmostEqual(result.parameters["fwhm_u"], 1.102781e-2)
        self.assertAlmostEqual(result.parameters["cell_a"], 8.87529)
        self.assertAlmostEqual(result.parameters["cell_alpha"], 90.0)
        self.assertFalse(result.observed_intensity.flags.writeable)
        with self.assertRaises(ValueError):
            result.observed_intensity[0] = 0.0

    def test_backend_copies_unicode_cif_to_ascii_job_and_runs_vesta_then_rietan(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            unicode_dir = root / "结构"
            unicode_dir.mkdir()
            cif = unicode_dir / "目标相.cif"
            cif.write_text("data_phase\n", encoding="utf-8")
            vesta = root / "VESTA.exe"
            rietan = root / "RIETAN.exe"
            vesta.write_bytes(b"")
            rietan.write_bytes(b"")
            job_dirs: list[Path] = []

            def fake_run(command, **kwargs):
                work = Path(kwargs["cwd"])
                job_dirs.append(work)
                if Path(command[0]).name.lower() == "vesta.exe":
                    self.assertEqual(Path(command[3]).name, "phase.cif")
                    self.assertNotIn("结构", str(command[3]))
                    (work / "ophi.ins").write_text(REFINEMENT_INS, encoding="utf-8", newline="")
                    return subprocess.CompletedProcess(command, 0, stdout="saved", stderr="")
                self.assertTrue((work / "ophi.int").is_file())
                patched = (work / "ophi.ins").read_text(encoding="utf-8")
                self.assertIn("NMODE = 0:", patched)
                (work / "ophi.gpd").write_text(REFINEMENT_GPD, encoding="ascii", newline="")
                return subprocess.CompletedProcess(command, 0, stdout=REFINEMENT_LST, stderr="")

            backend = RietanRefinementBackend(vesta_exe=vesta, rietan_exe=rietan)
            with mock.patch("ophiuchus.xrd.refinement.subprocess.run", side_effect=fake_run) as run:
                result = backend.refine(
                    cif,
                    np.array([10.0, 10.1, 10.2]),
                    np.array([12.0, 20.0, 13.0]),
                    RefinementSettings(10.0, 10.2),
                )

        self.assertEqual(run.call_count, 2)
        self.assertEqual(result.rwp_percent, 8.568)
        self.assertEqual(result.provenance["engine"], "RIETAN-FP constrained refinement")
        self.assertTrue(job_dirs)
        self.assertFalse(job_dirs[0].exists())

    def test_trust_warnings_flag_large_shift_and_implausible_high_angle_width(self):
        result = parse_refinement_output(REFINEMENT_GPD, REFINEMENT_LST)
        parameters = dict(result.parameters)
        parameters.update({"zero_shift": 0.21, "fwhm_u": -0.0145, "fwhm_v": 0.445, "fwhm_w": -0.035})
        result = type(result)(
            two_theta_deg=result.two_theta_deg,
            observed_intensity=result.observed_intensity,
            calculated_intensity=result.calculated_intensity,
            residual_intensity=result.residual_intensity,
            background_intensity=result.background_intensity,
            reflection_two_theta_deg=result.reflection_two_theta_deg,
            rwp_percent=result.rwp_percent,
            rp_percent=result.rp_percent,
            goodness_of_fit=result.goodness_of_fit,
            parameters=parameters,
        )

        warnings = refinement_trust_warnings(result, RefinementSettings(10.0, 90.0))

        self.assertTrue(any("zero shift" in warning.lower() for warning in warnings))
        self.assertTrue(any("fwhm" in warning.lower() for warning in warnings))


if __name__ == "__main__":
    unittest.main()
