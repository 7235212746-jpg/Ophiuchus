import os
import tempfile
import unittest
import zipfile
from pathlib import Path

import numpy as np

from ophiuchus.xrd.multiphase_models import MultiphaseRefinementSettings, PhaseRefinementInput
from ophiuchus.xrd.multiphase_refinement import RietanMultiphaseBackend
from ophiuchus.xrd.rietan_assets import discover_rietan_archive
from ophiuchus.xrd.backend import SimulationContext
from ophiuchus.xrd.config import XRDConfig
from ophiuchus.xrd.rietan_backend import RietanXRDBackend, discover_rietan_executable, discover_vesta_executable


@unittest.skipUnless(
    os.environ.get("OPHI_RUN_RIETAN_INTEGRATION") == "1",
    "real RIETAN integration disabled",
)
class RietanMultiphaseIntegrationTests(unittest.TestCase):
    def _assert_user_cif_high_angle_reference(
        self,
        formula: str,
        expected_sha256: str,
        expected_peaks: tuple[tuple[float, float], ...],
    ) -> None:
        vesta = discover_vesta_executable()
        rietan = discover_rietan_executable()
        matches = [
            path
            for path in (Path.home() / "Desktop").rglob(f"{formula}.cif")
            if ".worktrees" not in str(path)
        ]
        if vesta is None or rietan is None or not matches:
            self.skipTest(f"VESTA, RIETAN, or {formula}.cif is unavailable")
        pattern = RietanXRDBackend(vesta_exe=vesta, rietan_exe=rietan).simulate_cif(
            matches[0],
            XRDConfig(two_theta_min=10.0, two_theta_max=90.0, intensity_threshold=0.5),
            SimulationContext(formula, "RIETAN 3.12 golden reference", formula),
        )
        self.assertEqual(pattern.cif_sha256, expected_sha256)
        for expected_angle, expected_intensity in expected_peaks:
            index = min(
                range(len(pattern.two_theta_deg)),
                key=lambda item: abs(pattern.two_theta_deg[item] - expected_angle),
            )
            self.assertAlmostEqual(pattern.two_theta_deg[index], expected_angle, delta=0.011)
            self.assertAlmostEqual(pattern.normalized_intensity[index], expected_intensity, delta=0.5)

    def test_official_three_phase_example_runs_and_preserves_high_angle_reflections(self):
        archive = discover_rietan_archive()
        rietan = discover_rietan_executable()
        if archive is None or rietan is None:
            self.skipTest("official RIETAN archive or executable is unavailable")

        suffixes = {
            "Cu3Fe4_PO4_6@1.cif",
            "Cu3_PO4_2@2.cif",
            "Cu2P2O7@3.cif",
            "multi_phase.int",
            "template.ins",
        }
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with zipfile.ZipFile(archive) as bundle:
                for member in bundle.namelist():
                    name = Path(member).name
                    if "Cu3Fe4P6_combins" in member and name in suffixes:
                        (root / name).write_bytes(bundle.read(member))

            lines = (root / "multi_phase.int").read_text(encoding="utf-8").splitlines()
            count_text, start_text, step_text = lines[1].split()[:3]
            count = int(count_text)
            start = float(start_text)
            step = float(step_text)
            y = np.asarray([float(value) for line in lines[2:] for value in line.split()][:count])
            x = start + np.arange(count, dtype=float) * step
            phases = (
                PhaseRefinementInput("phase-1", "Cu3Fe4(PO4)6", root / "Cu3Fe4_PO4_6@1.cif", "target"),
                PhaseRefinementInput("phase-2", "Cu3(PO4)2", root / "Cu3_PO4_2@2.cif", "impurity"),
                PhaseRefinementInput("phase-3", "Cu2P2O7", root / "Cu2P2O7@3.cif", "impurity"),
            )
            result = RietanMultiphaseBackend(
                rietan_exe=rietan,
                template_path=root / "template.ins",
            ).refine(
                phases,
                x,
                y,
                MultiphaseRefinementSettings(float(x.min()), float(x.max())),
            )

        self.assertEqual(len(result.phases), 3)
        self.assertTrue(np.isfinite(result.rwp_percent))
        self.assertLess(result.rwp_percent, 5.0)
        self.assertAlmostEqual(sum(phase.weight_percent for phase in result.phases), 100.0, places=4)
        self.assertTrue(all(any(angle > 60.0 for angle in phase.reflection_two_theta_deg) for phase in result.phases))
        self.assertEqual(result.provenance["profile_engine"], "RIETAN-FP multiphase refinement")

    def test_zr3v3gesn4_rietan_312_high_angle_reference(self):
        self._assert_user_cif_high_angle_reference(
            "Zr3V3GeSn4",
            "156fbce64724a5d26c0f1a970fea747718a6be33dd791a06ec3815b0763c6117",
            ((67.21, 10.827), (70.93, 16.148), (80.54, 13.112), (85.44, 5.644)),
        )

    def test_zrfe6ge4_rietan_312_high_angle_reference(self):
        self._assert_user_cif_high_angle_reference(
            "ZrFe6Ge4",
            "362b46a32755d09f0a9076e5104279151c55ca06ddc4c1b74164d49606862a9b",
            ((62.59, 10.207), (70.07, 22.166), (74.80, 33.474), (81.11, 24.881), (88.61, 25.820)),
        )


if __name__ == "__main__":
    unittest.main()
