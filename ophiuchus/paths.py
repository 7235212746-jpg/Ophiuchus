from __future__ import annotations

import os
from pathlib import Path


def desktop_dir() -> Path:
    """Return the current user's real desktop without assuming OneDrive."""
    candidates: list[Path] = []
    if os.name == "nt":
        try:
            import winreg

            with winreg.OpenKey(
                winreg.HKEY_CURRENT_USER,
                r"Software\Microsoft\Windows\CurrentVersion\Explorer\User Shell Folders",
            ) as key:
                value, _kind = winreg.QueryValueEx(key, "Desktop")
            candidates.append(Path(os.path.expandvars(str(value))))
        except (OSError, TypeError, ValueError):
            pass
    candidates.extend((Path.home() / "Desktop", Path.home() / "OneDrive" / "Desktop"))
    for candidate in candidates:
        if candidate.is_dir():
            return candidate.resolve()
    return Path.home().resolve()


def first_existing_directory(path: str | Path | None, *, fallback: str | Path) -> Path:
    """Resolve a dialog start directory, walking up stale saved paths as needed."""
    fallback_path = Path(fallback).expanduser()
    if path:
        candidate = Path(path).expanduser()
        if candidate.is_file():
            return candidate.parent.resolve()
        if candidate.is_dir():
            return candidate.resolve()
        for parent in candidate.parents:
            if parent.is_dir() and parent != Path(parent.anchor):
                return parent.resolve()
    if fallback_path.is_file():
        return fallback_path.parent.resolve()
    if fallback_path.is_dir():
        return fallback_path.resolve()
    return desktop_dir()
