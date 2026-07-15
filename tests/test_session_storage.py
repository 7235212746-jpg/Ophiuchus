from pathlib import Path
from unittest import mock
import tempfile
import unittest

from ophiuchus.session_storage import TransientAnalysisStore


class TransientAnalysisStoreTests(unittest.TestCase):
    def test_commit_keeps_one_current_session_and_remaps_output_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            store = TransientAnalysisStore(Path(tmp) / "analysis-session")
            first_pending = store.begin()
            first_report = first_pending / "reports" / "analysis.json"
            first_report.parent.mkdir(parents=True)
            first_report.write_text("first", encoding="utf-8")
            current, outputs = store.commit({"json": str(first_report)})

            self.assertEqual(current, store.current_path)
            self.assertEqual(Path(outputs["json"]), current / "reports" / "analysis.json")
            self.assertEqual(Path(outputs["json"]).read_text(encoding="utf-8"), "first")
            self.assertFalse(store.pending_path.exists())
            self.assertFalse(store.previous_path.exists())

            second_pending = store.begin()
            second_report = second_pending / "reports" / "analysis.json"
            second_report.parent.mkdir(parents=True)
            second_report.write_text("second", encoding="utf-8")
            current, outputs = store.commit({"json": str(second_report)})

            self.assertEqual(Path(outputs["json"]).read_text(encoding="utf-8"), "second")
            self.assertFalse(store.previous_path.exists())
            self.assertEqual(sorted(path.name for path in current.iterdir()), ["reports"])

    def test_rollback_preserves_the_previous_successful_session(self):
        with tempfile.TemporaryDirectory() as tmp:
            store = TransientAnalysisStore(Path(tmp) / "analysis-session")
            initial = store.begin() / "analysis.json"
            initial.write_text("kept", encoding="utf-8")
            current, _ = store.commit({"json": str(initial)})

            failed = store.begin() / "analysis.json"
            failed.write_text("discarded", encoding="utf-8")
            restored = store.rollback()

            self.assertEqual(restored, current)
            self.assertEqual((store.current_path / "analysis.json").read_text(encoding="utf-8"), "kept")
            self.assertFalse(store.pending_path.exists())

    def test_save_current_copies_every_file_and_writes_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            store = TransientAnalysisStore(root / "analysis-session")
            pending = store.begin()
            (pending / "nested").mkdir()
            (pending / "analysis.json").write_text("analysis", encoding="utf-8")
            (pending / "nested" / "plot.png").write_bytes(b"plot")
            store.commit({"json": str(pending / "analysis.json")})

            destination = root / "saved-analysis"
            saved = store.save_current(destination)

            self.assertEqual(saved, destination)
            self.assertEqual((destination / "analysis.json").read_text(encoding="utf-8"), "analysis")
            self.assertEqual((destination / "nested" / "plot.png").read_bytes(), b"plot")
            self.assertTrue((destination / "manifest.json").is_file())

    def test_rollback_recovers_current_from_an_interrupted_commit(self):
        with tempfile.TemporaryDirectory() as tmp:
            store = TransientAnalysisStore(Path(tmp) / "analysis-session")
            store.root.mkdir(parents=True)
            store.previous_path.mkdir()
            (store.previous_path / "analysis.json").write_text("recover me", encoding="utf-8")
            store.pending_path.mkdir()

            restored = store.rollback()

            self.assertEqual(restored, store.current_path)
            self.assertEqual((store.current_path / "analysis.json").read_text(encoding="utf-8"), "recover me")
            self.assertFalse(store.previous_path.exists())

    def test_post_commit_backup_cleanup_failure_does_not_report_commit_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            store = TransientAnalysisStore(Path(tmp) / "analysis-session")
            first = store.begin() / "analysis.json"
            first.write_text("first", encoding="utf-8")
            store.commit({"json": first})
            second = store.begin() / "analysis.json"
            second.write_text("second", encoding="utf-8")
            original_remove = store._remove_owned_path

            def fail_only_for_committed_backup(path):
                if path == store.previous_path and path.exists() and store.current_path.exists():
                    raise PermissionError("locked backup")
                return original_remove(path)

            with mock.patch.object(store, "_remove_owned_path", side_effect=fail_only_for_committed_backup):
                current, outputs = store.commit({"json": second})

            self.assertEqual((current / "analysis.json").read_text(encoding="utf-8"), "second")
            self.assertEqual(Path(outputs["json"]), current / "analysis.json")

    def test_commit_rejects_missing_output_before_replacing_current(self):
        with tempfile.TemporaryDirectory() as tmp:
            store = TransientAnalysisStore(Path(tmp) / "analysis-session")
            first = store.begin() / "analysis.json"
            first.write_text("first", encoding="utf-8")
            store.commit({"json": first})
            store.begin()

            with self.assertRaises(FileNotFoundError):
                store.commit({"json": store.pending_path / "missing.json"})

            self.assertEqual((store.current_path / "analysis.json").read_text(encoding="utf-8"), "first")

    def test_save_failure_leaves_no_partial_destination_and_can_retry(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            store = TransientAnalysisStore(root / "analysis-session")
            pending = store.begin()
            (pending / "analysis.json").write_text("analysis", encoding="utf-8")
            store.commit({"json": pending / "analysis.json"})
            destination = root / "saved-analysis"

            with mock.patch.object(store, "_manifest_entry", side_effect=OSError("hash failed")):
                with self.assertRaises(OSError):
                    store.save_current(destination)

            self.assertFalse(destination.exists())
            self.assertEqual(store.save_current(destination), destination)

    def test_save_rejects_destination_inside_current_and_source_links(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            store = TransientAnalysisStore(root / "analysis-session")
            pending = store.begin()
            (pending / "analysis.json").write_text("analysis", encoding="utf-8")
            store.commit({"json": pending / "analysis.json"})

            with self.assertRaises(ValueError):
                store.save_current(store.current_path / "nested-save")

            link = store.current_path / "linked"
            try:
                link.symlink_to(root, target_is_directory=True)
            except OSError:
                self.skipTest("Creating symlinks is not permitted on this Windows installation.")
            with self.assertRaises(ValueError):
                store.save_current(root / "saved-analysis")


if __name__ == "__main__":
    unittest.main()
