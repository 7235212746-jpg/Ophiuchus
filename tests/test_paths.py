import tempfile
import unittest
from pathlib import Path

from ophiuchus.paths import first_existing_directory


class PathResolutionTests(unittest.TestCase):
    def test_first_existing_directory_uses_parent_for_missing_saved_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp) / "existing"
            parent.mkdir()

            resolved = first_existing_directory(parent / "missing" / "pattern.asc", fallback=Path(tmp))

            self.assertEqual(resolved, parent)

    def test_first_existing_directory_uses_fallback_when_saved_path_has_no_existing_parent(self):
        with tempfile.TemporaryDirectory() as tmp:
            fallback = Path(tmp) / "fallback"
            fallback.mkdir()

            resolved = first_existing_directory(Path("Z:/definitely/missing/path"), fallback=fallback)

            self.assertEqual(resolved, fallback)


if __name__ == "__main__":
    unittest.main()
