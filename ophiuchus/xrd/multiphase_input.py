from __future__ import annotations

import os
from pathlib import Path
import re
import shutil
import subprocess

from .multiphase_models import PhaseRefinementInput


def _block_pattern(label: str, marker: str) -> re.Pattern[str]:
    return re.compile(
        rf"(?ms)^[ \t]*![ \t]*{re.escape(label)}[ \t]+@{re.escape(marker)}[^\r\n]*\r?\n"
        rf".*?^[ \t]*#[ \t]*End[ \t]+{re.escape(label)}[ \t]+@{re.escape(marker)}[^\r\n]*(?:\r?\n|$)"
    )


def _extract_block(text: str, label: str, marker: str, *, context: str) -> str:
    match = _block_pattern(label, marker).search(text)
    if match is None:
        raise ValueError(f"{context} is missing the required {label} @{marker} multiphase marker block.")
    return match.group(0)


def _replace_generic_block(template: str, label: str, replacement: str) -> str:
    pattern = _block_pattern(label, "N")
    match = pattern.search(template)
    if match is None:
        raise ValueError(f"Official multiphase template is missing the required {label} @N marker block.")
    return template[: match.start()] + replacement + template[match.end() :]


def _element_symbols(text: str, phase_number: int) -> set[str]:
    block = _extract_block(text, "Elements", str(phase_number), context=f"RIETAN phase {phase_number} input")
    lines = block.splitlines()[1:-1]
    return set(re.findall(r"'[^']+'", "\n".join(lines)))


def combine_phase_inputs(template_text: str, phase_texts: list[str]) -> str:
    if not 2 <= len(phase_texts) <= 4:
        raise ValueError("Multiphase RIETAN input requires two to four phases.")
    for label in ("Elements", "Phase", "Parameters", "Constraints"):
        _extract_block(template_text, label, "N", context="Official multiphase template")
    if not re.search(r"(?m)^\s*NPHASE@\s*=\s*1\s*:", template_text):
        raise ValueError("Official multiphase template is missing the NPHASE@ = 1 declaration.")

    phase_blocks: list[str] = []
    parameter_blocks: list[str] = []
    elements: set[str] = set()
    for phase_number, phase_text in enumerate(phase_texts, 1):
        context = f"RIETAN phase {phase_number} input"
        phase_blocks.append(_extract_block(phase_text, "Phase", str(phase_number), context=context).rstrip())
        parameter_blocks.append(
            _extract_block(phase_text, "Parameters", str(phase_number), context=context).rstrip()
        )
        elements.update(_element_symbols(phase_text, phase_number))

    combined = _replace_generic_block(template_text, "Phase", "\n\n".join(phase_blocks) + "\n")
    combined = _replace_generic_block(combined, "Parameters", "\n\n".join(parameter_blocks) + "\n")

    constraint_template = _extract_block(combined, "Constraints", "N", context="Official multiphase template")
    constraints = "\n\n".join(
        constraint_template.replace("@N", f"@{phase_number}").rstrip()
        for phase_number in range(2, len(phase_texts) + 1)
    )
    combined = _replace_generic_block(combined, "Constraints", constraints + ("\n" if constraints else ""))

    element_line = "  " + "  ".join(sorted(elements, key=str.casefold)) + " /\n"
    element_replacement = "! Elements @N\n" + element_line + "# End Elements @N\n"
    combined = _replace_generic_block(combined, "Elements", element_replacement)
    combined, count = re.subn(
        r"(?m)^(\s*NPHASE@\s*=\s*)1(\s*:)",
        rf"\g<1>{len(phase_texts)}\g<2>",
        combined,
        count=1,
    )
    if count != 1:
        raise ValueError("Official multiphase template NPHASE declaration could not be updated.")
    return combined.replace("\r\n", "\n").replace("\r", "\n").replace("\n", "\r\n")


def export_phase_input(
    phase: PhaseRefinementInput,
    phase_number: int,
    work_dir: str | Path,
    cif2ins_exe: str | Path,
    template_path: str | Path,
    *,
    timeout_seconds: float = 60.0,
) -> Path:
    if not 1 <= int(phase_number) <= 4:
        raise ValueError("RIETAN phase number must be between 1 and 4.")
    work = Path(work_dir).resolve()
    tool = Path(cif2ins_exe).resolve()
    template = Path(template_path).resolve()
    if not work.is_dir():
        raise FileNotFoundError(f"RIETAN work directory does not exist: {work}")
    if not tool.is_file():
        raise FileNotFoundError(f"cif2ins executable does not exist: {tool}")
    if not template.is_file():
        raise FileNotFoundError(f"RIETAN multiphase template does not exist: {template}")

    stem = f"phase@{phase_number}"
    local_cif = work / f"{stem}.cif"
    local_template = work / "template.ins"
    output = work / f"{stem}.ins"
    shutil.copy2(phase.cif_path, local_cif)
    if template != local_template:
        shutil.copy2(template, local_template)

    command = [
        str(tool),
        "0",
        local_cif.name,
        local_template.name,
        output.name,
        "report.tex",
        "report.pdf",
        "structure.pdf",
        "result.lst",
        "mscs.pdf",
        "density.pdf",
    ]
    environment = os.environ.copy()
    environment["CIF2INS"] = str(tool.parent)
    completed = subprocess.run(
        command,
        cwd=work,
        env=environment,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=float(timeout_seconds),
    )
    if completed.returncode != 0 or not output.is_file() or output.stat().st_size == 0:
        detail = (completed.stderr or completed.stdout or "no .ins file was produced").strip()
        raise RuntimeError(f"cif2ins could not export RIETAN phase {phase_number}: {detail}")
    numbered = output.read_text(encoding="utf-8", errors="replace").replace("@N", f"@{phase_number}")
    output.write_bytes(numbered.encode("utf-8"))
    return output
