from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock

from ophiuchus.xrd.multiphase_input import combine_phase_inputs, export_phase_input
from ophiuchus.xrd.multiphase_models import PhaseRefinementInput
from ophiuchus.xrd.rietan_backend import discover_cif2ins_executable


TEMPLATE = """Title
case 1, 2
! Elements @N
'Zr' /
# End Elements @N
Data concerning crystalline phases contained in the sample {
! Phase(s)
! Phase @N
PHNAME@N = 'template'
# End Phase @N
# When two or more phases are included in the sample, repeat their data below.
} End of information about phases.
! Profile
! Parameters @N
SCALE@N 1.0 1
X@N/X 1 0 0 0 0 00000
# End Parameters @N
} End of lines for label/species, A(I), and ID(I)
! Linear constraints
A(X@1,B)=A(X@1,B)
! Constraints @N
A(FWHM@N,1)=A(FWHM@1,1)
# End Constraints @N
} End of linear constraints.
NPHASE@ = 1: Number of phases contained in this sample.
Tail
"""

PHASE1 = """! Elements @1
'Zr' 'Fe' /
# End Elements @1
! Phase @1
PHNAME@1 = 'main'
# End Phase @1
! Parameters @1
SCALE@1 1.0 1
Zr1@1/Zr 1 0 0 0 0 00000
# End Parameters @1
"""

PHASE2 = """! Elements @2
'Sn' 'Fe' /
# End Elements @2
! Phase @2
PHNAME@2 = 'impurity'
# End Phase @2
! Parameters @2
SCALE@2 0.2 1
Sn1@2/Sn 1 0 0 0 0 00000
# End Parameters @2
"""


class MultiphaseInputTests(unittest.TestCase):
    def test_combines_numbered_official_phase_sections(self):
        combined = combine_phase_inputs(TEMPLATE, [PHASE1, PHASE2])

        self.assertIn("NPHASE@ = 2:", combined)
        self.assertLess(combined.index("! Phase @1"), combined.index("! Phase @2"))
        self.assertLess(combined.index("! Parameters @1"), combined.index("! Parameters @2"))
        self.assertEqual(combined.count("Zr1@1/Zr"), 1)
        self.assertEqual(combined.count("Sn1@2/Sn"), 1)
        self.assertIn("'Fe'  'Sn'  'Zr' /", combined)
        self.assertNotIn("! Constraints @N", combined)
        self.assertIn("A(FWHM@2,1)=A(FWHM@1,1)", combined)

    def test_rejects_a_vesta_template_without_combination_markers(self):
        with self.assertRaisesRegex(ValueError, "multiphase template"):
            combine_phase_inputs("NMODE = 0", [PHASE1, PHASE2])

    def test_rejects_non_contiguous_phase_numbers(self):
        with self.assertRaisesRegex(ValueError, "phase 2"):
            combine_phase_inputs(TEMPLATE, [PHASE1, PHASE2.replace("@2", "@3")])

    def test_discovers_cif2ins_next_to_configured_rietan(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rietan = root / "RIETAN.exe"
            cif2ins = root / "cif2ins.exe"
            rietan.write_bytes(b"")
            cif2ins.write_bytes(b"")

            self.assertEqual(discover_cif2ins_executable(rietan_exe=rietan), cif2ins.resolve())

    def test_export_uses_ascii_job_name_for_unicode_cif(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source_dir = root / "结构"
            source_dir.mkdir()
            source = source_dir / "目标相.cif"
            source.write_text("data_phase\n", encoding="utf-8")
            phase = PhaseRefinementInput("main", "Main", source, "target")
            tool = root / "cif2ins.exe"
            tool.write_bytes(b"")
            template = root / "template_source.ins"
            template.write_text(TEMPLATE, encoding="utf-8")
            work = root / "job"
            work.mkdir()

            def fake_run(command, **kwargs):
                self.assertEqual(command[1:4], ["0", "phase@1.cif", "template.ins"])
                self.assertEqual(kwargs["cwd"], work)
                self.assertTrue((work / "phase@1.cif").is_file())
                self.assertNotIn("结构", str(work / "phase@1.cif"))
                (work / "phase@1.ins").write_text(PHASE1.replace("@1", "@N"), encoding="utf-8")
                return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

            with mock.patch("ophiuchus.xrd.multiphase_input.subprocess.run", side_effect=fake_run):
                result = export_phase_input(phase, 1, work, tool, template)

            self.assertEqual(result, work / "phase@1.ins")
            self.assertIn("! Phase @1", result.read_text(encoding="utf-8"))
            self.assertNotIn("@N", result.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
