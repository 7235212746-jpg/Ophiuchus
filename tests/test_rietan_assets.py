import tempfile
import unittest
import zipfile
from pathlib import Path


from ophiuchus.xrd.rietan_assets import install_multiphase_template_from_archive


VALID_TEMPLATE = """\
! Elements @N
  'Cu' /
# End Elements @N
! Phase @N
PHNAME@ = 'phase':
# End Phase @N
! Parameters @N
SCALE@ = 1.0:
# End Parameters @N
! Constraints @N
SCALE@ = SCALE@1:
# End Constraints @N
NPHASE@ = 1:
"""


class RietanAssetTests(unittest.TestCase):
    def test_installs_the_official_combins_template_from_nested_archive(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            archive = root / "Windows_versions.zip"
            member = (
                "Windows_versions/Windows_versions/RIETAN_VENUS_examples/"
                "Cu3Fe4P6_combins/template.ins"
            )
            with zipfile.ZipFile(archive, "w") as bundle:
                bundle.writestr(member, VALID_TEMPLATE)

            installed = install_multiphase_template_from_archive(archive, root / "runtime")

            self.assertEqual(installed, root / "runtime" / "template_multiphase.ins")
            self.assertEqual(installed.read_text(encoding="utf-8"), VALID_TEMPLATE)

    def test_rejects_a_template_without_official_multiphase_markers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            archive = root / "Windows_versions.zip"
            with zipfile.ZipFile(archive, "w") as bundle:
                bundle.writestr(
                    "RIETAN_VENUS_examples/Cu3Fe4P6_combins/template.ins",
                    "NMODE = 0:\n",
                )

            with self.assertRaisesRegex(ValueError, "multiphase template"):
                install_multiphase_template_from_archive(archive, root / "runtime")

            self.assertFalse((root / "runtime" / "template_multiphase.ins").exists())


if __name__ == "__main__":
    unittest.main()
