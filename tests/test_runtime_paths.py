import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock


from ophiuchus.app import default_app_state_path, default_cache_path, default_env_path, default_library_path
from ophiuchus.help_support import manual_path
from ophiuchus.runtime import resource_root, user_data_root


class RuntimePathTests(unittest.TestCase):
    def test_source_mode_uses_repository_for_resources_and_mutable_data(self):
        expected = Path(__file__).resolve().parents[1]

        self.assertEqual(resource_root(), expected)
        self.assertEqual(user_data_root(), expected)

    def test_frozen_mode_separates_bundled_resources_from_local_user_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            resources = root / "runtime" / "_internal"
            local = root / "LocalAppData"
            resources.mkdir(parents=True)
            with (
                mock.patch("ophiuchus.runtime.sys.frozen", True, create=True),
                mock.patch("ophiuchus.runtime.sys._MEIPASS", str(resources), create=True),
                mock.patch.dict(os.environ, {"LOCALAPPDATA": str(local)}, clear=False),
            ):
                self.assertEqual(resource_root(), resources)
                self.assertEqual(user_data_root(), local / "Ophiuchus")
                self.assertEqual(default_cache_path(), local / "Ophiuchus" / "data" / "ophi_xrd_cache.sqlite")
                self.assertEqual(default_library_path(), local / "Ophiuchus" / "data" / "ophi_library.sqlite")
                self.assertEqual(default_env_path(), local / "Ophiuchus" / ".env")
                self.assertEqual(default_app_state_path(), local / "Ophiuchus" / "data" / "ophi_app_state.json")
                self.assertEqual(manual_path(), resources / "docs" / "Ophiuchus_操作手册.md")


if __name__ == "__main__":
    unittest.main()
