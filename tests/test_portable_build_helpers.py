import tempfile
import unittest
from pathlib import Path

from tools.portable_build_helpers import discover_python_modules


class PortableBuildHelperTests(unittest.TestCase):
    def test_discovers_modules_below_pep420_namespace_directories(self):
        with tempfile.TemporaryDirectory() as temporary:
            site_packages = Path(temporary)
            package = site_packages / "pymatgen"
            elasticity = package / "analysis" / "elasticity"
            elasticity.mkdir(parents=True)
            (elasticity / "__init__.py").write_text("", encoding="utf-8")
            (elasticity / "elastic.py").write_text("VALUE = 1\n", encoding="utf-8")

            modules = discover_python_modules(site_packages, "pymatgen")

        self.assertIn("pymatgen.analysis.elasticity", modules)
        self.assertIn("pymatgen.analysis.elasticity.elastic", modules)


if __name__ == "__main__":
    unittest.main()
