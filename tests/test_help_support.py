import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from ophiuchus.help_support import CONTACT_EMAIL, manual_path, open_contact_email, open_manual


class HelpSupportTests(unittest.TestCase):
    def test_contact_email_is_the_maintainer_address(self):
        self.assertEqual(CONTACT_EMAIL, "wanyc@issp.u-tokyo.ac.jp")

    def test_manual_path_lives_in_project_docs(self):
        root = Path(r"C:\Ophiuchus")

        self.assertEqual(manual_path(root), root / "docs" / "Ophiuchus_操作手册.md")

    def test_open_manual_uses_windows_file_association(self):
        with tempfile.TemporaryDirectory() as tmp:
            guide = Path(tmp) / "guide.md"
            guide.write_text("guide", encoding="utf-8")
            with mock.patch("ophiuchus.help_support.os.name", "nt"), mock.patch(
                "ophiuchus.help_support.os.startfile", create=True
            ) as startfile:
                open_manual(guide)

        startfile.assert_called_once_with(str(guide))

    def test_open_manual_rejects_missing_file(self):
        missing = Path(tempfile.gettempdir()) / "missing-ophiuchus-guide.md"
        if missing.exists():
            os.remove(missing)

        with self.assertRaises(FileNotFoundError):
            open_manual(missing)

    def test_open_contact_email_uses_encoded_mailto_url(self):
        with mock.patch("ophiuchus.help_support.webbrowser.open", return_value=True) as open_url:
            open_contact_email()

        url = open_url.call_args.args[0]
        self.assertTrue(url.startswith("mailto:wanyc@issp.u-tokyo.ac.jp?subject="))
        self.assertIn("Ophiuchus%20", url)


if __name__ == "__main__":
    unittest.main()
