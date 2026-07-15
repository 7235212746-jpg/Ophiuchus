from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
from typing import Mapping
import uuid


class TransientAnalysisStore:
    """Own a single replaceable analysis bundle outside permanent results."""

    def __init__(self, root: str | Path | None = None) -> None:
        default_root = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local")) / "Ophiuchus" / "analysis-session"
        self.root = Path(root) if root is not None else default_root
        self.pending_path = self.root / "pending"
        self.current_path = self.root / "current"
        self.previous_path = self.root / "previous"

    def begin(self) -> Path:
        self.root.mkdir(parents=True, exist_ok=True)
        self._recover_interrupted_commit()
        self._remove_owned_path(self.pending_path)
        self.pending_path.mkdir()
        return self.pending_path

    def commit(self, output_paths: Mapping[str, str | Path]) -> tuple[Path, dict[str, str]]:
        if not self.pending_path.is_dir():
            raise RuntimeError("begin() must create a pending analysis session before commit().")
        remapped_outputs = self._remap_output_paths(output_paths)
        for name, value in output_paths.items():
            if not Path(value).exists():
                raise FileNotFoundError(f"Analysis output '{name}' does not exist: {value}")
        self._recover_interrupted_commit()
        self._remove_owned_path(self.previous_path)
        moved_current = False
        try:
            if self.current_path.exists():
                self.current_path.replace(self.previous_path)
                moved_current = True
            self.pending_path.replace(self.current_path)
        except Exception:
            if moved_current and self.previous_path.exists() and not self.current_path.exists():
                self.previous_path.replace(self.current_path)
            raise
        try:
            self._remove_owned_path(self.previous_path)
        except OSError:
            # The new current bundle is already committed. A locked backup is
            # harmless and will be cleaned or recovered on the next operation.
            pass
        return self.current_path, remapped_outputs

    def rollback(self) -> Path | None:
        self._remove_owned_path(self.pending_path)
        self._recover_interrupted_commit()
        return self.current_path if self.current_path.exists() else None

    def save_current(self, destination: str | Path) -> Path:
        if not self.current_path.is_dir():
            raise RuntimeError("No current analysis session is available to save.")
        target = Path(destination)
        if target.exists():
            raise FileExistsError(f"Save destination already exists: {target}")
        current_resolved = self.current_path.resolve()
        target_resolved = target.resolve()
        try:
            target_resolved.relative_to(current_resolved)
        except ValueError:
            pass
        else:
            raise ValueError("Save destination cannot be inside the transient current session.")
        self._assert_plain_tree(self.current_path)
        target.parent.mkdir(parents=True, exist_ok=True)
        staging = target.parent / f".{target.name}.ophi-saving-{uuid.uuid4().hex}"
        try:
            shutil.copytree(self.current_path, staging, symlinks=True)
            source_manifest = staging / "manifest.json"
            if source_manifest.exists():
                preserved = staging / "source_manifest.json"
                if preserved.exists():
                    raise FileExistsError("The session already contains both manifest.json and source_manifest.json.")
                source_manifest.replace(preserved)
            manifest = {
                "format": "ophiuchus.analysis-session-manifest.v1",
                "files": [self._manifest_entry(path, staging) for path in sorted(staging.rglob("*")) if path.is_file()],
            }
            (staging / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
            staging.replace(target)
        except Exception:
            self._remove_owned_path(staging)
            raise
        return target

    def _remap_output_paths(self, output_paths: Mapping[str, str | Path]) -> dict[str, str]:
        pending = self.pending_path.resolve()
        current = self.current_path.resolve()
        remapped: dict[str, str] = {}
        for name, value in output_paths.items():
            source = Path(value).resolve()
            try:
                relative = source.relative_to(pending)
            except ValueError as exc:
                raise ValueError(f"Output path is outside the pending analysis session: {source}") from exc
            remapped[str(name)] = str(current / relative)
        return remapped

    def _recover_interrupted_commit(self) -> None:
        if self.previous_path.exists() and not self.current_path.exists():
            self.previous_path.replace(self.current_path)
        elif self.previous_path.exists() and self.current_path.exists():
            try:
                self._remove_owned_path(self.previous_path)
            except OSError:
                pass

    @staticmethod
    def _assert_plain_tree(root: Path) -> None:
        for directory, names, filenames in os.walk(root, followlinks=False):
            base = Path(directory)
            for name in [*names, *filenames]:
                path = base / name
                is_junction = bool(getattr(os.path, "isjunction", lambda _path: False)(path))
                if path.is_symlink() or is_junction:
                    raise ValueError(f"Transient session contains a link and cannot be saved safely: {path}")

    @staticmethod
    def _manifest_entry(path: Path, root: Path) -> dict[str, str]:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        return {"path": path.relative_to(root).as_posix(), "sha256": digest}

    @staticmethod
    def _remove_owned_path(path: Path) -> None:
        if path.is_dir():
            shutil.rmtree(path)
        elif path.exists():
            path.unlink()
