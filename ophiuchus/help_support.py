from __future__ import annotations

import os
import webbrowser
from pathlib import Path
from urllib.parse import quote


CONTACT_EMAIL = "wanyc@issp.u-tokyo.ac.jp"
MANUAL_FILENAME = "Ophiuchus_操作手册.md"


def manual_path(root: Path | None = None) -> Path:
    project = root or Path(__file__).resolve().parents[1]
    return project / "docs" / MANUAL_FILENAME


def open_manual(path: Path | None = None) -> None:
    target = Path(path) if path is not None else manual_path()
    if not target.is_file():
        raise FileNotFoundError(f"Ophiuchus user manual not found: {target}")
    if os.name == "nt":
        os.startfile(str(target))
    else:
        webbrowser.open(target.resolve().as_uri())


def open_contact_email(email: str = CONTACT_EMAIL) -> None:
    subject = quote("Ophiuchus 使用反馈")
    webbrowser.open(f"mailto:{email}?subject={subject}")
