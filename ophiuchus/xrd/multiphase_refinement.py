from __future__ import annotations

from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
from typing import Iterable

import numpy as np
from pymatgen.core import Structure

from .multiphase_input import combine_phase_inputs, export_phase_input
from .multiphase_models import (
    MultiphaseRefinementResult,
    MultiphaseRefinementSettings,
    PhaseRefinementInput,
    PhaseRefinementResult,
)
from .quantification import weight_fractions_from_zmv
from .refinement import write_xy_intensity
from .rietan_backend import (
    build_rietan_command,
    discover_cif2ins_executable,
    discover_multiphase_template,
    discover_rietan_executable,
)


_FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


def _activate_choice(text: str, name: str, selected: str) -> str:
    pattern = re.compile(rf"(?m)^(\s*{re.escape(name)}\s*=\s*)([^\s:!]+)([:!])(.*)$")
    plain_selected = bool(re.search(rf"(?m)^\s*{re.escape(name)}\s*=\s*{re.escape(selected)}\s*$", text))

    def replace_choice(match: re.Match[str]) -> str:
        marker = ":" if not plain_selected and match.group(2) == selected else "!"
        return f"{match.group(1)}{match.group(2)}{marker}{match.group(4)}"

    return pattern.sub(replace_choice, text)


def _replace_ids(text: str, label_pattern: str, ids: str) -> str:
    pattern = re.compile(rf"(?m)^(\s*(?:{label_pattern})\b.*?\s)([0-3]{{{len(ids)}}})(\s*(?:#.*)?)$")
    return pattern.sub(lambda match: f"{match.group(1)}{ids}{match.group(3)}", text)


def _replace_background_ids(text: str, ids: str) -> str:
    pattern = re.compile(r"(?ms)^(\s*(?:BKGD|BACKGROUND)\b.*?\s)([0-3]{12})(\s*(?:#.*)?)$")
    patched, count = pattern.subn(lambda match: f"{match.group(1)}{ids}{match.group(3)}", text, count=1)
    if count != 1:
        raise RuntimeError("RIETAN multiphase template is missing a 12-term background parameter block.")
    return patched


def _patch_scale(text: str, phase_number: int, value: float) -> str:
    pattern = re.compile(
        rf"(?m)^(\s*SCALE@{phase_number}\s+)({_FLOAT})(\s+)([0-3])(\s*(?:#.*)?)$"
    )
    patched, count = pattern.subn(
        lambda match: f"{match.group(1)}{value:.8g}{match.group(3)}1{match.group(5)}",
        text,
        count=1,
    )
    if count != 1:
        raise RuntimeError(f"RIETAN multiphase input is missing SCALE@{phase_number}.")
    return patched


def _patch_phase_profile_ids(text: str, phase_number: int, refine_profile: bool) -> str:
    labels = (
        "GAUSS01",
        "LORENTZ01",
        "ASYM01",
        "GAUSS00",
        "LORENTZ00",
        "ASYM00",
        "FWHM12",
        "DECAY12",
        "ETA12",
        "FWHM3",
        "M3",
    )
    pattern = re.compile(
        rf"(?m)^(\s*(?:{'|'.join(labels)})@{phase_number}\b.*?\s)([0-3]+)(\s*(?:#.*)?)$"
    )

    def replace_profile(match: re.Match[str]) -> str:
        current = match.group(2)
        if not refine_profile:
            ids = "0" * len(current)
        elif phase_number == 1:
            ids = "".join("1" if value != "0" else "0" for value in current)
        else:
            ids = "".join("2" if value != "0" else "0" for value in current)
        return f"{match.group(1)}{ids}{match.group(3)}"

    return pattern.sub(replace_profile, text)


def _lock_structure_parameters(text: str) -> str:
    pattern = re.compile(
        r"(?m)^(\s*[^#!\r\n\s]+@\d+/[^\s]+\s+.*?\s)([0-3]{5,7})(\s*(?:#.*)?)$"
    )
    return pattern.sub(
        lambda match: f"{match.group(1)}{'0' * len(match.group(2))}{match.group(3)}",
        text,
    )


def patch_multiphase_refinement_input(
    text: str,
    phases: Iterable[PhaseRefinementInput],
    settings: MultiphaseRefinementSettings,
) -> str:
    phase_list = tuple(phases)
    if not 2 <= len(phase_list) <= 4:
        raise ValueError("Multiphase RIETAN refinement requires two to four phases.")
    patched = text.replace("\r\n", "\n").replace("\r", "\n")
    for name, selected in (
        ("NMODE", "0"),
        ("NPRINT", "1"),
        ("NINT", "1"),
        ("NRANGE", "0"),
        ("NLESQ", "0"),
        ("NAUTO", "2"),
        ("NPRFN", "0"),
        ("NASYM", "1"),
        ("NPAT", "1"),
    ):
        patched = _activate_choice(patched, name, selected)
    patched = _replace_background_ids(
        patched,
        "1" * settings.background_terms + "0" * (12 - settings.background_terms),
    )
    patched = _replace_ids(
        patched,
        r"SHIFT0|SHIFTN",
        "1000" if settings.refine_zero_shift else "0000",
    )
    for phase_number, phase in enumerate(phase_list, 1):
        patched = _patch_scale(patched, phase_number, phase.initial_scale)
        patched = _patch_phase_profile_ids(patched, phase_number, settings.refine_profile)
        patched = _replace_ids(patched, rf"PREF@{phase_number}", "000000")
        patched = _replace_ids(
            patched,
            rf"CELLQ@{phase_number}",
            "1111110" if settings.refine_lattice else "0000000",
        )
        patched = _replace_ids(patched, rf"ANISTR0@{phase_number}|ANISOBR12@{phase_number}", "00")
    patched = _lock_structure_parameters(patched)
    required = (
        r"(?m)^\s*NMODE\s*=\s*0(?::|\s*$)",
        r"(?m)^\s*NINT\s*=\s*1(?::|\s*$)",
        r"(?m)^\s*NRANGE\s*=\s*0(?::|\s*$)",
        r"(?m)^\s*NPRFN\s*=\s*0(?::|\s*$)",
    )
    if any(not re.search(pattern, patched) for pattern in required):
        raise RuntimeError("RIETAN multiphase template could not activate the conservative refinement modes.")
    return patched.replace("\n", "\r\n")


def _parse_profile(gpd_text: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rows: list[tuple[float, float, float, float, float]] = []
    for raw_line in gpd_text.splitlines():
        line = raw_line.strip()
        if not line:
            if rows:
                break
            continue
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            row = tuple(float(value) for value in parts[:5])
        except ValueError:
            continue
        rows.append(row)  # type: ignore[arg-type]
    if not rows:
        raise RuntimeError("RIETAN multiphase .gpd did not contain an observed/calculated profile.")
    columns = tuple(np.asarray(values, dtype=float) for values in zip(*rows))
    return columns  # type: ignore[return-value]


def _parse_phase_reflections(gpd_text: str, phase_count: int) -> list[tuple[float, ...]]:
    positions: list[list[float]] = [[] for _ in range(phase_count)]
    active: int | None = None
    header = re.compile(r"phase\s+No\.\s*(\d+)", re.IGNORECASE)
    for raw_line in gpd_text.splitlines():
        line = raw_line.strip()
        if line.startswith("#"):
            match = header.search(line)
            active = int(match.group(1)) - 1 if match else None
            continue
        if active is None or not 0 <= active < phase_count or not line:
            continue
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            positions[active].append(float(parts[3]))
        except ValueError:
            continue
    return [tuple(values) for values in positions]


def _parse_final_metrics(lst_text: str) -> tuple[float, float, float]:
    matches = re.findall(
        rf"Rwp\s*=\s*({_FLOAT}).*?Rp\s*=\s*({_FLOAT}).*?GofF\s*=\s*({_FLOAT})",
        lst_text,
        flags=re.IGNORECASE,
    )
    if not matches:
        raise RuntimeError("RIETAN multiphase output did not contain final Rwp/Rp/GofF statistics.")
    return tuple(float(value) for value in matches[-1])  # type: ignore[return-value]


def _parse_final_scales(lst_text: str, phase_count: int) -> list[float]:
    marker = "*** Final parameters and their estimated standard uncertainties ***"
    if marker not in lst_text:
        raise RuntimeError("RIETAN multiphase output did not contain a final parameter table.")
    final = lst_text.rsplit(marker, 1)[-1]
    scales = [
        float(match.group(1))
        for match in re.finditer(
            rf"(?m)^\s*\d+\s+({_FLOAT})(?:\s+{_FLOAT}){{0,2}}\s+Scale factor,\s*s\s*$",
            final,
            flags=re.IGNORECASE,
        )
    ]
    if len(scales) < phase_count:
        raise RuntimeError("RIETAN multiphase final table did not contain one scale factor per phase.")
    return scales[:phase_count]


def _parse_native_weight_percent(lst_text: str, phase_count: int) -> list[float]:
    marker = "Effective radii (R), particle absorption factors (tau), and mass/mole fractions"
    index = lst_text.find(marker)
    if index < 0:
        raise RuntimeError("RIETAN multiphase output did not contain its native mass-fraction table.")
    rows: list[float] = []
    header_seen = False
    for raw_line in lst_text[index:].splitlines()[1:]:
        line = raw_line.strip()
        if not line:
            if rows:
                break
            continue
        if not header_seen:
            if re.search(r"\bPhase\b.*\bw\b", line):
                header_seen = True
            continue
        numbers = [float(value) for value in re.findall(_FLOAT, line)]
        if len(numbers) < 9:
            if rows:
                break
            continue
        rows.append(numbers[-4] * 100.0)
        if len(rows) == phase_count:
            break
    if len(rows) != phase_count:
        raise RuntimeError("RIETAN native mass-fraction table did not contain one row per phase.")
    return rows


def parse_multiphase_output(
    gpd_text: str,
    lst_text: str,
    phases: Iterable[PhaseRefinementInput],
    zmv_values: Iterable[tuple[float, float, float]],
) -> MultiphaseRefinementResult:
    phase_list = tuple(phases)
    zmv_list = tuple(zmv_values)
    if len(phase_list) != len(zmv_list):
        raise ValueError("Each RIETAN phase requires one Z/M/V tuple.")
    x, observed, calculated, residual, background = _parse_profile(gpd_text)
    reflections = _parse_phase_reflections(gpd_text, len(phase_list))
    scales = _parse_final_scales(lst_text, len(phase_list))
    native_weights = _parse_native_weight_percent(lst_text, len(phase_list))
    rwp, rp, gof = _parse_final_metrics(lst_text)
    hill_howard = weight_fractions_from_zmv(
        (
            phase.phase_id,
            scale,
            z,
            molar_mass,
            volume,
        )
        for phase, scale, (z, molar_mass, volume) in zip(phase_list, scales, zmv_list)
    )
    warnings: list[str] = []
    max_difference = max(
        abs(native - hill_howard[phase.phase_id])
        for phase, native in zip(phase_list, native_weights)
    )
    if max_difference > 0.5:
        warnings.append(
            f"RIETAN native mass fractions differ from the CIF ZMV cross-check by up to {max_difference:.3f} wt%."
        )
    phase_results = tuple(
        PhaseRefinementResult(
            phase_id=phase.phase_id,
            formula=phase.formula,
            scale=scale,
            z=z,
            molar_mass=molar_mass,
            volume_angstrom3=volume,
            weight_percent=native_weight,
            reflection_two_theta_deg=reflection_positions,
        )
        for phase, scale, (z, molar_mass, volume), native_weight, reflection_positions in zip(
            phase_list,
            scales,
            zmv_list,
            native_weights,
            reflections,
        )
    )
    return MultiphaseRefinementResult(
        two_theta_deg=x,
        observed_intensity=observed,
        calculated_intensity=calculated,
        residual_intensity=residual,
        background_intensity=background,
        phases=phase_results,
        rwp_percent=rwp,
        rp_percent=rp,
        goodness_of_fit=gof,
        warnings=tuple(warnings),
        provenance={"zmv_max_difference_wt_percent": max_difference},
    )


def structure_zmv(cif_path: str | Path) -> tuple[float, float, float]:
    structure = Structure.from_file(str(Path(cif_path)))
    reduced, factor = structure.composition.get_reduced_composition_and_factor()
    z = float(factor)
    molar_mass = float(reduced.weight)
    volume = float(structure.volume)
    values = np.asarray([z, molar_mass, volume], dtype=float)
    if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
        raise RuntimeError(f"CIF did not provide positive finite Z/M/V values: {cif_path}")
    return z, molar_mass, volume


class RietanMultiphaseBackend:
    def __init__(
        self,
        *,
        rietan_exe: str | Path | None = None,
        cif2ins_exe: str | Path | None = None,
        template_path: str | Path | None = None,
        timeout_seconds: float = 240.0,
    ) -> None:
        self.rietan_exe = discover_rietan_executable(rietan_exe)
        self.cif2ins_exe = discover_cif2ins_executable(cif2ins_exe, rietan_exe=self.rietan_exe)
        self.template_path = discover_multiphase_template(template_path, rietan_exe=self.rietan_exe)
        self.timeout_seconds = float(timeout_seconds)

    @property
    def available(self) -> bool:
        return self.rietan_exe is not None and self.cif2ins_exe is not None and self.template_path is not None

    def refine(
        self,
        phases: Iterable[PhaseRefinementInput],
        x: np.ndarray,
        intensity: np.ndarray,
        settings: MultiphaseRefinementSettings,
    ) -> MultiphaseRefinementResult:
        phase_list = tuple(phases)
        if not self.available:
            raise RuntimeError(
                "Multiphase RIETAN is unavailable; configure RIETAN.exe, cif2ins.exe, and the official multiphase template."
            )
        if not 2 <= len(phase_list) <= 4:
            raise ValueError("Multiphase RIETAN refinement requires two to four phases.")
        if phase_list[0].role != "target" or any(phase.role != "impurity" for phase in phase_list[1:]):
            raise ValueError("The first phase must be the target and all remaining phases must be impurities.")
        hashes = [phase.cif_sha256 for phase in phase_list]
        if len(set(hashes)) != len(hashes):
            raise ValueError("Duplicate CIF structures cannot enter the same multiphase refinement.")
        x_values = np.asarray(x, dtype=float)
        y_values = np.asarray(intensity, dtype=float)
        if x_values.ndim != 1 or y_values.ndim != 1 or x_values.shape != y_values.shape:
            raise ValueError("Experimental refinement arrays must be equal one-dimensional arrays.")
        selected = (x_values >= settings.two_theta_min) & (x_values <= settings.two_theta_max)
        if int(np.count_nonzero(selected)) < 3:
            raise ValueError("The active multiphase refinement range contains fewer than three points.")
        assert self.rietan_exe is not None
        assert self.cif2ins_exe is not None
        assert self.template_path is not None

        with tempfile.TemporaryDirectory(prefix="ophi_multi_refine_") as temp:
            work = Path(temp)
            phase_inputs = [
                export_phase_input(
                    phase,
                    phase_number,
                    work,
                    self.cif2ins_exe,
                    self.template_path,
                    timeout_seconds=min(self.timeout_seconds, 60.0),
                ).read_text(encoding="utf-8", errors="replace")
                for phase_number, phase in enumerate(phase_list, 1)
            ]
            template_text = self.template_path.read_text(encoding="utf-8", errors="replace")
            combined = combine_phase_inputs(template_text, phase_inputs)
            patched = patch_multiphase_refinement_input(combined, phase_list, settings)
            (work / "multi_phase.ins").write_bytes(patched.encode("utf-8"))
            write_xy_intensity(work / "multi_phase.int", x_values[selected], y_values[selected])
            environment = os.environ.copy()
            environment["RIETAN"] = str(self.rietan_exe.parent) + os.sep
            command = build_rietan_command(self.rietan_exe, "multi_phase")
            completed = subprocess.run(
                command,
                cwd=work,
                env=environment,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=self.timeout_seconds,
            )
            output = f"{completed.stdout or ''}\n{completed.stderr or ''}"
            gpd = work / "multi_phase.gpd"
            if (
                completed.returncode != 0
                or not gpd.is_file()
                or re.search(r"error occurred|forrtl: severe|singular matrix", output, re.IGNORECASE)
            ):
                detail = "\n".join(line for line in output.splitlines()[-18:] if line.strip())
                raise RuntimeError(f"RIETAN-FP multiphase refinement failed: {detail or 'no fitted profile was produced'}")
            zmv = [structure_zmv(phase.cif_path) for phase in phase_list]
            result = parse_multiphase_output(
                gpd.read_text(encoding="utf-8", errors="replace"),
                output,
                phase_list,
                zmv,
            )

        provenance = {
            **result.provenance,
            "profile_engine": "RIETAN-FP multiphase refinement",
            "phase_count": len(phase_list),
            "phase_cif_sha256": {phase.phase_id: phase.cif_sha256 for phase in phase_list},
            "rietan_executable": str(self.rietan_exe),
            "cif2ins_executable": str(self.cif2ins_exe),
            "template_path": str(self.template_path),
            "input_sha256": hashlib.sha256(patched.encode("utf-8")).hexdigest(),
            "settings_sha256": hashlib.sha256(
                json.dumps(settings.__dict__, sort_keys=True, separators=(",", ":")).encode("ascii")
            ).hexdigest(),
        }
        return replace(result, provenance=provenance)
