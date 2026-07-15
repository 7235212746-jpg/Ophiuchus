import re
import unittest
from pathlib import Path


class ShareReadinessTests(unittest.TestCase):
    def test_start_script_discovers_python_without_author_specific_paths(self):
        root = Path(__file__).resolve().parents[1]
        script = (root / "start_ophiuchus.bat").read_text(encoding="utf-8").lower()

        self.assertIsNone(re.search(r"c:\\users\\[^\\\r\n]+\\", script))
        self.assertIn("ophi_python", script)
        self.assertIn(r"%userprofile%\anaconda3\envs\ophi\python.exe", script)
        self.assertIn(r"%userprofile%\miniconda3\envs\ophi\python.exe", script)
        self.assertIn("conda run -n ophi python", script)

    def test_readme_has_no_author_specific_absolute_paths(self):
        root = Path(__file__).resolve().parents[1]
        readme = (root / "README.md").read_text(encoding="utf-8").lower()

        self.assertIsNone(re.search(r"c:\\users\\[^\\\r\n]+\\", readme))

    def test_test_runner_discovers_python_without_author_specific_paths(self):
        root = Path(__file__).resolve().parents[1]
        script = (root / "run_tests.bat").read_text(encoding="utf-8").lower()

        self.assertIsNone(re.search(r"c:\\users\\[^\\\r\n]+\\", script))
        self.assertIn("ophi_python", script)
        self.assertIn("conda run -n ophi python", script)

    def test_gitignore_excludes_local_research_data_and_internal_work_notes(self):
        root = Path(__file__).resolve().parents[1]
        ignored = (root / ".gitignore").read_text(encoding="utf-8").splitlines()

        for required in [
            "data/",
            "results/",
            "exports/",
            "tmp_vesta_probe/",
            "xrd_sandbox/inputs/",
            "xrd_sandbox/runs/",
            ".superpowers/",
            "docs/superpowers/",
        ]:
            self.assertIn(required, ignored)


if __name__ == "__main__":
    unittest.main()
