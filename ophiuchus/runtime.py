from __future__ import annotations

import os
from pathlib import Path
import sys


def is_frozen() -> bool:
    return bool(getattr(sys, "frozen", False))


def resource_root() -> Path:
    if is_frozen() and getattr(sys, "_MEIPASS", None):
        return Path(str(sys._MEIPASS)).resolve()  # type: ignore[attr-defined]
    return Path(__file__).resolve().parents[1]


def user_data_root() -> Path:
    if not is_frozen():
        return resource_root()
    local = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    return (local / "Ophiuchus").resolve()
