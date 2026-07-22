import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock


from ophiuchus.portable_entry import main, portable_health_report


class PortableEntryTests(unittest.TestCase):
    def test_health_report_checks_core_runtime_without_exposing_secrets(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            resources = root / "resources"
            user_data = root / "user-data"
            manual = resources / "docs" / "Ophiuchus_操作手册.md"
            manual.parent.mkdir(parents=True)
            manual.write_text("manual", encoding="utf-8")
            with (
                mock.patch("ophiuchus.portable_entry.resource_root", return_value=resources),
                mock.patch("ophiuchus.portable_entry.user_data_root", return_value=user_data),
                mock.patch("ophiuchus.portable_entry.discover_vesta_executable", return_value=None),
                mock.patch("ophiuchus.portable_entry.discover_rietan_executable", return_value=None),
                mock.patch("ophiuchus.portable_entry.discover_cif2ins_executable", return_value=None),
                mock.patch("ophiuchus.portable_entry.discover_multiphase_template", return_value=None),
            ):
                report = portable_health_report()

        self.assertTrue(report["core_ready"])
        self.assertTrue(report["manual_available"])
        self.assertTrue(report["user_data_writable"])
        self.assertEqual(set(report["dependencies"]), {"numpy", "scipy", "matplotlib", "pymatgen", "mp-api"})
        self.assertEqual(report["optional_engines"]["rietan"], False)
        serialized = json.dumps(report).lower()
        self.assertNotIn("api_key", serialized)
        self.assertNotIn("mp_api_key", serialized)

    def test_health_check_command_writes_json_and_returns_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "health.json"
            with mock.patch(
                "ophiuchus.portable_entry.portable_health_report",
                return_value={"core_ready": True, "manual_available": True, "user_data_writable": True},
            ):
                exit_code = main(["--health-check", str(target)])

            self.assertEqual(exit_code, 0)
            self.assertTrue(json.loads(target.read_text(encoding="utf-8"))["core_ready"])


if __name__ == "__main__":
    unittest.main()
