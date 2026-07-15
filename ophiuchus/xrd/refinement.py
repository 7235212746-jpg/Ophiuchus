from __future__ import annotations

from dataclasses import dataclass, field, replace
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import numpy as np

from .rietan_backend import build_rietan_command, discover_rietan_executable, discover_vesta_executable


_FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


@dataclass(frozen=True)
class RefinementSettings:
    two_theta_min: float
    two_theta_max: float
    background_terms: int = 6
    refine_zero_shift: bool = True
    refine_profile: bool = True
    refine_lattice: bool = False
    radiation: str = "CuKa"

    def __post_init__(self) -> None:
        if not (
            np.isfinite(self.two_theta_min)
            and np.isfinite(self.two_theta_max)
            and self.two_theta_min < self.two_theta_max
        ):
            raise ValueError("Refinement range must contain two finite increasing values.")
        if not 1 <= int(self.background_terms) <= 12:
            raise ValueError("RIETAN background_terms must be between 1 and 12.")
        if self.radiation not in {"CuKa", "CuKalpha12", "CuKα"}:
            raise ValueError("The constrained RIETAN refinement currently supports Cu Kalpha only.")


@dataclass(frozen=True)
class RietanRefinementResult:
    two_theta_deg: np.ndarray
    observed_intensity: np.ndarray
    calculated_intensity: np.ndarray
    residual_intensity: np.ndarray
    background_intensity: np.ndarray
    reflection_two_theta_deg: tuple[float, ...]
    rwp_percent: float
    rp_percent: float
    goodness_of_fit: float
    s_value: float | None = None
    parameters: dict[str, float] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()
    provenance: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        arrays: list[np.ndarray] = []
        for name in (
            "two_theta_deg",
            "observed_intensity",
            "calculated_intensity",
            "residual_intensity",
            "background_intensity",
        ):
            values = np.array(getattr(self, name), dtype=float, copy=True)
            if values.ndim != 1 or not np.all(np.isfinite(values)):
                raise ValueError(f"{name} must be a finite one-dimensional array.")
            values.setflags(write=False)
            object.__setattr__(self, name, values)
            arrays.append(values)
        if len({array.size for array in arrays}) != 1 or not arrays[0].size:
            raise ValueError("All RIETAN refinement profile arrays must have the same non-zero length.")
        if np.any(np.diff(arrays[0]) <= 0.0):
            raise ValueError("RIETAN refinement 2theta values must be strictly increasing.")


def write_xy_intensity(path: str | Path, x: np.ndarray, intensity: np.ndarray) -> str:
    x_values = np.asarray(x, dtype=float)
    y_values = np.asarray(intensity, dtype=float)
    if x_values.ndim != 1 or y_values.ndim != 1 or x_values.shape != y_values.shape:
        raise ValueError("Experimental 2theta and intensity arrays must be equal one-dimensional arrays.")
    if x_values.size < 2 or not np.all(np.isfinite(x_values)) or not np.all(np.isfinite(y_values)):
        raise ValueError("Experimental refinement data must contain at least two finite points.")
    if np.any(np.diff(x_values) <= 0.0):
        raise ValueError("Experimental 2theta values must be strictly increasing.")
    target = Path(path)
    target.write_text(
        "GENERAL\n"
        f"{x_values.size}\n"
        + "".join(f"{two_theta:.6f} {value:.6f}\n" for two_theta, value in zip(x_values, y_values)),
        encoding="ascii",
        newline="",
    )
    return str(target)


def patch_refinement_input(text: str, settings: RefinementSettings) -> str:
    patched = text.replace("\r\n", "\n").replace("\r", "\n")
    for name, selected in (
        ("NMODE", "0"),
        ("NPRINT", "1"),
        ("NINT", "1"),
        ("NRANGE", "0"),
        ("NLESQ", "0"),
        ("NAUTO", "2"),
        ("NPAT", "1"),
    ):
        patched = _activate_choice(patched, name, selected)

    patched = _replace_parameter_ids(patched, "SHIFTN", "1000" if settings.refine_zero_shift else "0000")
    patched = _replace_parameter_ids(
        patched,
        "BKGD",
        "1" * settings.background_terms + "0" * (12 - settings.background_terms),
    )
    patched = _replace_parameter_ids(patched, "SCALE", "1")
    patched = _replace_parameter_ids(patched, "FWHM12", "1110" if settings.refine_profile else "0000")
    patched = _replace_parameter_ids(patched, "ASYM12", "0000")
    patched = _replace_parameter_ids(patched, "ETA12", "0000")
    patched = _replace_parameter_ids(patched, "ANISOBR12", "00")
    patched = _replace_parameter_ids(patched, "PREF", "000000")
    patched = _replace_parameter_ids(patched, "CELLQ", "1111110" if settings.refine_lattice else "0000000")
    patched = _lock_structure_parameters(patched)

    required = {
        "NMODE": r"(?m)^\s*NMODE\s*=\s*0:",
        "NINT": r"(?m)^\s*NINT\s*=\s*1:",
        "NRANGE": r"(?m)^\s*NRANGE\s*=\s*0:",
        "SCALE": r"(?m)^\s*SCALE\b",
        "FWHM12": r"(?m)^\s*FWHM12\b",
    }
    missing = [name for name, pattern in required.items() if not re.search(pattern, patched)]
    if missing:
        raise RuntimeError("VESTA RIETAN template is missing required refinement entries: " + ", ".join(missing))
    return patched.replace("\n", "\r\n")


def _activate_choice(text: str, name: str, selected: str) -> str:
    pattern = re.compile(rf"(?m)^(\s*{re.escape(name)}\s*=\s*)([^\s:!]+)([:!])(.*)$")
    plain_selected = bool(re.search(rf"(?m)^\s*{re.escape(name)}\s*=\s*{re.escape(selected)}\s*$", text))

    def replace(match: re.Match[str]) -> str:
        marker = ":" if not plain_selected and match.group(2) == selected else "!"
        return f"{match.group(1)}{match.group(2)}{marker}{match.group(4)}"

    return pattern.sub(replace, text)


def _replace_parameter_ids(text: str, label: str, ids: str) -> str:
    pattern = re.compile(rf"(?m)^(\s*{re.escape(label)}\b.*?\s)([0-3]{{{len(ids)}}})(\s*(?:#.*)?)$")
    return pattern.sub(lambda match: f"{match.group(1)}{ids}{match.group(3)}", text, count=1)


def _lock_structure_parameters(text: str) -> str:
    start = text.find("# Label/(chemical species")
    end = text.find("} End of lines for label/species", start)
    if start < 0 or end < 0:
        return re.sub(
            r"(?m)^(\s*[A-Za-z0-9_.+-]+/[A-Za-z0-9_.+-]+\s+.*?\s)([0-3]{5,7})(\s*)$",
            lambda match: f"{match.group(1)}{'0' * len(match.group(2))}{match.group(3)}",
            text,
        )
    section = text[start:end]
    section = re.sub(
        r"(?m)^(\s*[^#!\r\n]+/[^\s]+.*?\s)([0-3]{5,7})(\s*)$",
        lambda match: f"{match.group(1)}{'0' * len(match.group(2))}{match.group(3)}",
        section,
    )
    return text[:start] + section + text[end:]


def parse_refinement_output(gpd_text: str, lst_text: str) -> RietanRefinementResult:
    profile_rows: list[tuple[float, float, float, float, float]] = []
    reflection_positions: list[float] = []
    in_profile = True
    for raw_line in gpd_text.splitlines():
        line = raw_line.strip()
        if not line:
            if profile_rows:
                in_profile = False
            continue
        if line.startswith("#"):
            continue
        parts = line.split()
        try:
            values = [float(item) for item in parts]
        except ValueError:
            continue
        if in_profile and len(values) >= 5:
            profile_rows.append(tuple(values[:5]))
        elif not in_profile and len(values) >= 4:
            reflection_positions.append(float(values[3]))
    if len(profile_rows) < 3:
        raise RuntimeError("RIETAN refinement .gpd did not contain an observed/calculated profile.")

    metric_matches = list(
        re.finditer(
            rf"Rwp\s*=\s*({_FLOAT}).*?Rp\s*=\s*({_FLOAT}).*?\bS\s*=\s*({_FLOAT}).*?GofF\s*=\s*({_FLOAT})",
            lst_text,
            re.IGNORECASE,
        )
    )
    if not metric_matches:
        raise RuntimeError("RIETAN refinement output did not contain final Rwp/Rp/GofF statistics.")
    metric = metric_matches[-1]
    parameters = _parse_final_parameters(lst_text)
    data = np.asarray(profile_rows, dtype=float)
    return RietanRefinementResult(
        two_theta_deg=data[:, 0],
        observed_intensity=data[:, 1],
        calculated_intensity=data[:, 2],
        residual_intensity=data[:, 3],
        background_intensity=data[:, 4],
        reflection_two_theta_deg=tuple(reflection_positions),
        rwp_percent=float(metric.group(1)),
        rp_percent=float(metric.group(2)),
        goodness_of_fit=float(metric.group(4)),
        s_value=float(metric.group(3)),
        parameters=parameters,
    )


def _parse_final_parameters(text: str) -> dict[str, float]:
    final = text.rsplit("Final parameters and their estimated standard uncertainties", 1)[-1]
    descriptions = {
        "Peak-shift parameter, t0": "zero_shift",
        "Scale factor, s": "scale",
        "FWHM parameter, U": "fwhm_u",
        "FWHM parameter, V": "fwhm_v",
        "FWHM parameter, W": "fwhm_w",
        "Lattice parameter, a": "cell_a",
        "Lattice parameter, b": "cell_b",
        "Lattice parameter, c": "cell_c",
        "Lattice parameter, alpha": "cell_alpha",
        "Lattice parameter, beta": "cell_beta",
        "Lattice parameter, gamma": "cell_gamma",
    }
    parameters: dict[str, float] = {}
    for line in final.splitlines():
        match = re.match(rf"\s*\d+\s+({_FLOAT})\b", line)
        if not match:
            continue
        for description, key in sorted(descriptions.items(), key=lambda item: len(item[0]), reverse=True):
            if line.rstrip().endswith(description):
                parameters[key] = float(match.group(1))
                break
    return parameters


def refinement_trust_warnings(
    result: RietanRefinementResult,
    settings: RefinementSettings,
) -> tuple[str, ...]:
    warnings: list[str] = []
    shift = result.parameters.get("zero_shift")
    if shift is not None and abs(shift) > 0.15:
        warnings.append(
            f"Refined zero shift is large ({shift:.4f} deg); verify sample displacement and target-phase identity."
        )
    profile_values = tuple(result.parameters.get(name) for name in ("fwhm_u", "fwhm_v", "fwhm_w"))
    if all(value is not None for value in profile_values):
        u, v, w = (float(value) for value in profile_values)
        theta = np.radians(np.linspace(settings.two_theta_min, settings.two_theta_max, 240) / 2.0)
        width_squared = u * np.tan(theta) ** 2 + v * np.tan(theta) + w
        if np.any(width_squared <= 0.0):
            warnings.append("Refined Caglioti FWHM becomes non-positive inside the active 2theta range.")
        else:
            maximum_width = float(np.sqrt(width_squared).max())
            if maximum_width > 0.5:
                warnings.append(
                    f"Refined FWHM reaches {maximum_width:.3f} deg; an instrument-standard calibration is required."
                )
    if result.rwp_percent > 30.0:
        warnings.append(
            f"Rwp remains high ({result.rwp_percent:.2f}%); this single phase does not explain the complete pattern."
        )
    return tuple(warnings)


class RietanRefinementBackend:
    def __init__(
        self,
        *,
        vesta_exe: str | Path | None = None,
        rietan_exe: str | Path | None = None,
        timeout_seconds: float = 180.0,
    ) -> None:
        self.vesta_exe = discover_vesta_executable(vesta_exe)
        self.rietan_exe = discover_rietan_executable(rietan_exe)
        self.timeout_seconds = float(timeout_seconds)

    @property
    def available(self) -> bool:
        return self.vesta_exe is not None and self.rietan_exe is not None

    def refine(
        self,
        cif_path: str | Path,
        x: np.ndarray,
        intensity: np.ndarray,
        settings: RefinementSettings,
    ) -> RietanRefinementResult:
        if not self.available:
            raise RuntimeError("RIETAN refinement is unavailable; configure both VESTA.exe and RIETAN.exe.")
        source_cif = Path(cif_path).resolve()
        if not source_cif.is_file():
            raise FileNotFoundError(f"Refinement CIF does not exist: {source_cif}")
        x_values = np.asarray(x, dtype=float)
        y_values = np.asarray(intensity, dtype=float)
        if x_values.shape != y_values.shape:
            raise ValueError("Experimental refinement arrays must have equal shapes.")
        selected = (x_values >= settings.two_theta_min) & (x_values <= settings.two_theta_max)
        if int(np.count_nonzero(selected)) < 3:
            raise ValueError("The active refinement range contains fewer than three experimental points.")
        assert self.vesta_exe is not None
        assert self.rietan_exe is not None

        with tempfile.TemporaryDirectory(prefix="ophi_refine_") as temp:
            work = Path(temp)
            local_cif = work / "phase.cif"
            shutil.copy2(source_cif, local_cif)
            ins_path = work / "ophi.ins"
            vesta_command = [
                str(self.vesta_exe),
                "-nogui",
                "-i",
                str(local_cif),
                "-save",
                "format=rietan",
                str(ins_path),
            ]
            exported = subprocess.run(
                vesta_command,
                cwd=work,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=min(self.timeout_seconds, 60.0),
            )
            if not ins_path.is_file() or ins_path.stat().st_size == 0:
                detail = (exported.stderr or exported.stdout or "no .ins file was produced").strip()
                raise RuntimeError(f"VESTA could not export the refinement input: {detail}")
            original_ins = ins_path.read_text(encoding="utf-8", errors="replace")
            patched_ins = patch_refinement_input(original_ins, settings)
            ins_path.write_bytes(patched_ins.encode("utf-8"))
            write_xy_intensity(work / "ophi.int", x_values[selected], y_values[selected])

            environment = os.environ.copy()
            environment["RIETAN"] = str(self.rietan_exe.parent) + os.sep
            command = build_rietan_command(self.rietan_exe)
            refined = subprocess.run(
                command,
                cwd=work,
                env=environment,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=self.timeout_seconds,
            )
            output = f"{refined.stdout or ''}\n{refined.stderr or ''}"
            gpd_path = work / "ophi.gpd"
            if (
                refined.returncode != 0
                or not gpd_path.is_file()
                or re.search(r"error occurred|forrtl: severe|singular matrix", output, re.IGNORECASE)
            ):
                detail = "\n".join(line for line in output.splitlines()[-16:] if line.strip())
                raise RuntimeError(f"RIETAN-FP constrained refinement failed: {detail or 'no fitted profile was produced'}")
            result = parse_refinement_output(
                gpd_path.read_text(encoding="utf-8", errors="replace"),
                output,
            )

        provenance = {
            "engine": "RIETAN-FP constrained refinement",
            "cif_path": str(source_cif),
            "cif_sha256": _sha256_file(source_cif),
            "vesta_executable": str(self.vesta_exe),
            "rietan_executable": str(self.rietan_exe),
            "input_sha256": hashlib.sha256(patched_ins.encode("utf-8")).hexdigest(),
            "settings_sha256": hashlib.sha256(
                json.dumps(
                    {
                        "two_theta_min": settings.two_theta_min,
                        "two_theta_max": settings.two_theta_max,
                        "background_terms": settings.background_terms,
                        "refine_zero_shift": settings.refine_zero_shift,
                        "refine_profile": settings.refine_profile,
                        "refine_lattice": settings.refine_lattice,
                        "radiation": settings.radiation,
                    },
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode("ascii")
            ).hexdigest(),
        }
        return replace(
            result,
            warnings=refinement_trust_warnings(result, settings),
            provenance=provenance,
        )


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()
