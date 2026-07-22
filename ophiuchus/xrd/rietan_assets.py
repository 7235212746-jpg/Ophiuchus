from __future__ import annotations

import os
from pathlib import Path
import re
import zipfile


OFFICIAL_RIETAN_PAGE_URL = "https://jp-minerals.org/rietan/"
OFFICIAL_RIETAN_ARCHIVE_NAME = "Windows_versions.zip"
_OFFICIAL_TEMPLATE_SUFFIX = "rietan_venus_examples/cu3fe4p6_combins/template.ins"
_REQUIRED_TEMPLATE_MARKERS = (
    "Elements",
    "Phase",
    "Parameters",
    "Constraints",
)


def validate_multiphase_template_text(text: str) -> None:
    for marker in _REQUIRED_TEMPLATE_MARKERS:
        pattern = rf"(?mis)^\s*!\s*{marker}\s+@N\b.*?^\s*#\s*End\s+{marker}\s+@N\b"
        if re.search(pattern, text) is None:
            raise ValueError(f"Official multiphase template is missing the {marker} @N marker block.")
    if re.search(r"(?mi)^\s*NPHASE@\s*=\s*1\s*:", text) is None:
        raise ValueError("Official multiphase template is missing the NPHASE@ = 1 declaration.")


def install_multiphase_template_from_archive(
    archive_path: str | Path,
    destination_dir: str | Path,
) -> Path:
    archive = Path(archive_path).expanduser().resolve()
    destination = Path(destination_dir).expanduser().resolve()
    if not archive.is_file():
        raise FileNotFoundError(f"RIETAN official archive does not exist: {archive}")

    with zipfile.ZipFile(archive) as bundle:
        matches = [
            name
            for name in bundle.namelist()
            if name.replace("\\", "/").lower().endswith(_OFFICIAL_TEMPLATE_SUFFIX)
        ]
        if len(matches) != 1:
            raise ValueError(
                "RIETAN official archive does not contain exactly one Cu3Fe4P6_combins multiphase template."
            )
        text = bundle.read(matches[0]).decode("utf-8-sig", errors="replace")

    validate_multiphase_template_text(text)
    destination.mkdir(parents=True, exist_ok=True)
    installed = destination / "template_multiphase.ins"
    temporary = destination / ".template_multiphase.ins.tmp"
    temporary.write_bytes(text.encode("utf-8"))
    os.replace(temporary, installed)
    return installed


def default_rietan_archive_candidates() -> tuple[Path, ...]:
    local = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    return (
        local / "Ophiuchus" / "downloads" / OFFICIAL_RIETAN_ARCHIVE_NAME,
        Path.home() / "Downloads" / OFFICIAL_RIETAN_ARCHIVE_NAME,
        Path.home() / "Desktop" / OFFICIAL_RIETAN_ARCHIVE_NAME,
        Path.home() / "OneDrive" / "Desktop" / OFFICIAL_RIETAN_ARCHIVE_NAME,
    )


def discover_rietan_archive() -> Path | None:
    return next((path.resolve() for path in default_rietan_archive_candidates() if path.is_file()), None)
