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

    def test_shortcut_installer_resolves_project_and_desktop_at_runtime(self):
        root = Path(__file__).resolve().parents[1]
        installer = (root / "install_desktop_shortcut.ps1").read_text(encoding="utf-8").lower()

        self.assertIn("$psscriptroot", installer)
        self.assertIn('specialfolders.item("desktop")', installer)
        self.assertIn("ophiuchus.exe", installer)
        self.assertIn("start_ophiuchus.bat", installer)
        self.assertIsNone(re.search(r"c:\\users\\[^\\\r\n]+\\", installer))

    def test_windows_batch_launchers_use_crlf_and_repository_enforces_it(self):
        root = Path(__file__).resolve().parents[1]
        attributes = (root / ".gitattributes").read_text(encoding="utf-8")

        self.assertIn("*.bat text eol=crlf", attributes)
        for name in [
            "start_ophiuchus.bat",
            "run_tests.bat",
            "install_desktop_shortcut.bat",
            "build_launcher_exe.bat",
            "build_portable.bat",
        ]:
            raw = (root / name).read_bytes()
            self.assertNotIn(b"\n", raw.replace(b"\r\n", b""), name)

    def test_windows_exe_launcher_is_portable_and_rebuildable(self):
        root = Path(__file__).resolve().parents[1]
        source = (root / "tools" / "OphiuchusLauncher.cs").read_text(encoding="utf-8")
        builder = (root / "build_launcher_exe.ps1").read_text(encoding="utf-8").lower()

        self.assertIn("AppDomain.CurrentDomain.BaseDirectory", source)
        self.assertIn("start_ophiuchus.bat", source)
        self.assertIn("/target:winexe", builder)
        self.assertIn("$psscriptroot", builder)
        self.assertIsNone(re.search(r"c:\\users\\[^\\\r\n]+\\", source + builder))

        executable = root / "Ophiuchus.exe"
        self.assertTrue(executable.is_file())
        self.assertEqual(executable.read_bytes()[:2], b"MZ")

    def test_windows_launcher_prefers_self_contained_runtime_before_source_fallback(self):
        root = Path(__file__).resolve().parents[1]
        source = (root / "tools" / "OphiuchusLauncher.cs").read_text(encoding="utf-8")

        portable = source.index(r'runtime\OphiuchusApp.exe')
        source_fallback = source.index("start_ophiuchus.bat")
        self.assertLess(portable, source_fallback)
        self.assertIn("File.Exists(portableApp)", source)

    def test_portable_builder_uses_pyinstaller_health_gate_and_explicit_release_files(self):
        root = Path(__file__).resolve().parents[1]
        builder = (root / "build_portable.ps1").read_text(encoding="utf-8")
        spec = (root / "packaging" / "OphiuchusPortable.spec").read_text(encoding="utf-8")
        entry = (root / "tools" / "OphiuchusPortableEntry.py").read_text(encoding="utf-8")
        lock = (root / "requirements-portable-lock.txt").read_text(encoding="utf-8")

        self.assertIn("PyInstaller", builder)
        self.assertIn("--health-check", builder)
        self.assertIn("SHA256SUMS.txt", builder)
        self.assertIn("$releaseFiles", builder)
        self.assertNotIn('Copy-Item -Path (Join-Path $PSScriptRoot "*")', builder)
        self.assertIn("collect_all", spec)
        self.assertIn("copy_metadata", spec)
        self.assertIn("Ophiuchus_操作手册.md", spec)
        self.assertIn("root = Path(SPECPATH).resolve().parent", spec)
        self.assertNotIn("root = Path(SPECPATH).resolve().parent.parent", spec)
        self.assertIn('root / "tools" / "OphiuchusPortableEntry.py"', spec)
        self.assertIn("from ophiuchus.portable_entry import main", entry)
        self.assertIn("WaitForExit", builder)
        self.assertIn('Join-Path $pythonRoot "Library\\bin"', builder)
        self.assertIn('"PySide6"', spec)
        self.assertIn("requirements-portable-lock.txt", builder)
        self.assertIn("--target", builder)
        self.assertIn("OPHI_PORTABLE_SITE_PACKAGES", builder)
        self.assertIn("OPHI_PORTABLE_SITE_PACKAGES", spec)
        self.assertIn('discover_python_modules(portable_site_packages, "pymatgen")', spec)
        self.assertTrue(all("==" in line for line in lock.splitlines() if line and not line.startswith("#")))

    def test_gitignore_excludes_portable_build_outputs(self):
        root = Path(__file__).resolve().parents[1]
        ignored = (root / ".gitignore").read_text(encoding="utf-8").splitlines()

        self.assertIn("build/portable/", ignored)
        self.assertIn("build/portable-site-packages/", ignored)
        self.assertIn("dist/", ignored)

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
