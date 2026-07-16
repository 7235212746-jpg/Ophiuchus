from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
from typing import Iterable

import numpy as np
from scipy.signal import find_peaks

from ..paths import desktop_dir
from .backend import SimulatedPattern, SimulationContext
from .config import XRDConfig
from .models import Peak


ENGINE_NAME = "VESTA / RIETAN-FP"
BACKEND_VERSION = "vesta_rietan_profile_v1"
RIETAN_CU_KALPHA1_ANGSTROM = 1.540593
_OUTPUT_SUFFIXES = (
    "ins",
    "int",
    "bkg",
    "itx",
    "hkl",
    "xyz",
    "fos",
    "ffe",
    "fba",
    "ffi",
    "ffo",
    "vesta",
    "plt",
    "gpd",
    "alb",
    "prf",
    "inflip",
    "exp",
)


@dataclass(frozen=True)
class RietanReflection:
    two_theta: float
    d_spacing: float
    intensity: float
    hkl: str
    multiplicity: int
    line_component: str


@dataclass(frozen=True)
class RietanProfile:
    two_theta_deg: tuple[float, ...]
    raw_intensity: tuple[float, ...]
    normalized_intensity: tuple[float, ...]


def discover_vesta_executable(configured: str | Path | None = None) -> Path | None:
    desktop = desktop_dir()
    candidates = [
        configured,
        os.environ.get("OPHI_VESTA_EXE"),
        desktop / "02_科研软件" / "便携软件" / "VESTA-win64" / "VESTA.exe",
        desktop / "VESTA-win64" / "VESTA.exe",
        Path(os.environ.get("PROGRAMFILES", r"C:\Program Files")) / "VESTA-win64" / "VESTA.exe",
    ]
    return _first_existing_file(candidates)


def discover_rietan_executable(configured: str | Path | None = None) -> Path | None:
    local = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    candidates = [
        configured,
        os.environ.get("OPHI_RIETAN_EXE"),
        local / "Ophiuchus" / "tools" / "rietan" / "RIETAN.exe",
        Path(os.environ.get("PROGRAMFILES", r"C:\Program Files")) / "RIETAN_VENUS" / "RIETAN.exe",
    ]
    return _first_existing_file(candidates)


def _first_existing_file(candidates: Iterable[str | Path | None]) -> Path | None:
    for candidate in candidates:
        if not candidate:
            continue
        path = Path(candidate).expanduser()
        if path.is_file():
            return path.resolve()
    return None


def patch_rietan_input(text: str, two_theta_min: float, two_theta_max: float) -> str:
    if not (np.isfinite(two_theta_min) and np.isfinite(two_theta_max) and two_theta_min < two_theta_max):
        raise ValueError("RIETAN simulation range must contain two finite increasing values.")
    patched = _activate_choice(text, "NMODE", "1")
    patched = _activate_choice(patched, "NPRINT", "1")
    patched = _activate_choice(patched, "NPAT", "1")
    patched = _replace_assignment(patched, "DEG1", float(two_theta_min))
    patched = _replace_assignment(patched, "DEG2", float(two_theta_max))
    if not re.search(r"(?m)^\s*NMODE\s*=\s*1:", patched):
        raise RuntimeError("VESTA RIETAN template does not contain a selectable NMODE = 1 entry.")
    if not re.search(r"(?m)^\s*DEG1\s*=", patched) or not re.search(r"(?m)^\s*DEG2\s*=", patched):
        raise RuntimeError("VESTA RIETAN template does not contain DEG1/DEG2 simulation limits.")
    return patched.replace("\r\n", "\n").replace("\r", "\n").replace("\n", "\r\n")


def _activate_choice(text: str, name: str, selected: str) -> str:
    pattern = re.compile(rf"(?m)^(\s*{re.escape(name)}\s*=\s*)([^\s:!]+)([:!])(.*)$")
    plain_selected = bool(
        re.search(rf"(?m)^\s*{re.escape(name)}\s*=\s*{re.escape(selected)}\s*$", text)
    )

    def replace(match: re.Match[str]) -> str:
        marker = ":" if not plain_selected and match.group(2) == selected else "!"
        return f"{match.group(1)}{match.group(2)}{marker}{match.group(4)}"

    return pattern.sub(replace, text)


def _replace_assignment(text: str, name: str, value: float) -> str:
    pattern = re.compile(rf"(?m)^(\s*{re.escape(name)}\s*=\s*)[^:!\r\n]+([:!])(.*)$")
    return pattern.sub(lambda match: f"{match.group(1)}{value:.6f}:{match.group(3)}", text, count=1)


def parse_rietan_reflections(text: str) -> tuple[RietanReflection, ...]:
    reflections: list[RietanReflection] = []
    in_table = False
    for line in text.splitlines():
        if "No. Phase" in line and "2-theta" in line and "Ical" in line:
            in_table = True
            continue
        if not in_table:
            continue
        parts = line.split()
        if len(parts) < 15 or not parts[0].isdigit():
            if reflections and not line.strip():
                break
            continue
        try:
            h, k, l = (int(parts[2]), int(parts[3]), int(parts[4]))
            code = parts[5]
            reflections.append(
                RietanReflection(
                    two_theta=float(parts[6]),
                    d_spacing=float(parts[7]),
                    intensity=float(parts[8]),
                    hkl=f"({h} {k} {l})",
                    multiplicity=int(parts[13]),
                    line_component="Kalpha2" if code.endswith("2") else "Kalpha1",
                )
            )
        except (ValueError, IndexError):
            continue
    if not reflections:
        raise RuntimeError("RIETAN output did not contain a readable reflection table.")
    return tuple(reflections)


def parse_rietan_profile(text: str) -> RietanProfile:
    x: list[float] = []
    y: list[float] = []
    for line in text.splitlines():
        clean = line.strip()
        if not clean or clean.startswith("#"):
            if x and not clean:
                break
            continue
        parts = clean.split()
        if len(parts) != 2:
            if x:
                break
            continue
        try:
            point = (float(parts[0]), float(parts[1]))
        except ValueError:
            if x:
                break
            continue
        x.append(point[0])
        y.append(point[1])
    if len(x) < 3 or any(second <= first for first, second in zip(x, x[1:])):
        raise RuntimeError("RIETAN .gpd did not contain a valid increasing simulated profile.")
    maximum = max(y)
    if maximum <= 0.0:
        raise RuntimeError("RIETAN simulated profile has no positive intensity.")
    normalized = tuple(float(value) / maximum * 100.0 for value in y)
    return RietanProfile(tuple(x), tuple(y), normalized)


def profile_peaks(
    profile: RietanProfile,
    reflections: tuple[RietanReflection, ...],
    *,
    intensity_threshold: float = 1.0,
) -> list[Peak]:
    values = np.asarray(profile.normalized_intensity, dtype=float)
    indices, _ = find_peaks(values, prominence=1e-6)
    selected = list(int(index) for index in indices if values[index] >= intensity_threshold)
    if values.size >= 2 and values[0] > values[1] and values[0] >= intensity_threshold:
        selected.insert(0, 0)
    if values.size >= 2 and values[-1] > values[-2] and values[-1] >= intensity_threshold:
        selected.append(values.size - 1)
    peaks: list[Peak] = []
    for number, index in enumerate(selected, 1):
        position = profile.two_theta_deg[index]
        reflection = min(reflections, key=lambda item: abs(item.two_theta - position))
        peaks.append(
            Peak(
                peak_id=f"rietan_profile_{number}",
                two_theta=position,
                intensity=float(values[index]),
                hkl=reflection.hkl,
                d_spacing=reflection.d_spacing,
                multiplicity=reflection.multiplicity,
            )
        )
    return peaks


def build_rietan_command(executable: Path, sample: str = "ophi") -> list[str]:
    return [str(executable), *(f"{sample}.{suffix}" for suffix in _OUTPUT_SUFFIXES)]


def upgrade_scores_with_rietan(
    scores,
    backend,
    settings: XRDConfig,
    experimental_peaks: list[Peak],
    *,
    tolerance_deg: float,
    limit: int = 4,
):
    """Re-simulate only final display candidates and recompute their evidence scores."""
    from .matching import score_candidate

    upgraded = []
    warnings: list[str] = []
    for index, score in enumerate(scores):
        candidate = score.candidate
        if index >= limit or Path(candidate.source_path).suffix.lower() != ".cif":
            upgraded.append(score)
            continue
        try:
            pattern = backend.simulate_cif(
                candidate.source_path,
                settings,
                SimulationContext(
                    structure_id=candidate.candidate_id,
                    source=candidate.source,
                    formula=candidate.formula_pretty,
                    space_group_number=candidate.space_group_number,
                ),
            )
        except Exception as exc:
            warnings.append(f"{candidate.formula_pretty}: VESTA/RIETAN display simulation failed: {exc}")
            upgraded.append(score)
            continue
        candidate.simulated_pattern = pattern
        candidate.theory_peaks = pattern.to_peaks()
        candidate.parse_status = "vesta_rietan_simulated"
        candidate.simulation_state = "vesta_rietan_exact"
        candidate.simulation_validation.update(
            {
                "status": "passed",
                "reason": "Simulated directly through VESTA and RIETAN-FP.",
                "engine_name": pattern.engine_name,
                "engine_version": pattern.engine_version,
                "pattern_fingerprint": pattern.pattern_fingerprint,
                "cif_sha256": pattern.cif_sha256,
            }
        )
        upgraded.append(score_candidate(experimental_peaks, candidate, tolerance_deg=tolerance_deg))
    return upgraded, warnings


class RietanXRDBackend:
    name = ENGINE_NAME
    engine_version = BACKEND_VERSION

    def __init__(
        self,
        *,
        vesta_exe: str | Path | None = None,
        rietan_exe: str | Path | None = None,
        timeout_seconds: float = 120.0,
    ) -> None:
        self.vesta_exe = discover_vesta_executable(vesta_exe)
        self.rietan_exe = discover_rietan_executable(rietan_exe)
        self.timeout_seconds = float(timeout_seconds)

    @property
    def available(self) -> bool:
        return self.vesta_exe is not None and self.rietan_exe is not None

    def settings_fingerprint(self, settings: XRDConfig) -> str:
        return _fingerprint(
            {
                **settings.to_dict(),
                "backend_name": self.name,
                "backend_version": self.engine_version,
                "vesta_exe": "" if self.vesta_exe is None else str(self.vesta_exe),
                "rietan_exe": "" if self.rietan_exe is None else str(self.rietan_exe),
                "rietan_cu_kalpha1_angstrom": RIETAN_CU_KALPHA1_ANGSTROM,
            }
        )

    def simulate_cif(
        self,
        cif_path: str | Path,
        settings: XRDConfig,
        context: SimulationContext,
    ) -> SimulatedPattern:
        if not self.available:
            raise RuntimeError("VESTA/RIETAN backend is unavailable; configure both VESTA.exe and RIETAN.exe.")
        if settings.radiation_source not in {"CuKa", "CuKalpha12", "CuKα"}:
            raise RuntimeError("The VESTA/RIETAN backend currently supports the Cu Kalpha1+2 template only.")
        cif = Path(cif_path).resolve()
        if not cif.is_file():
            raise FileNotFoundError(f"CIF file does not exist: {cif}")
        assert self.vesta_exe is not None
        assert self.rietan_exe is not None
        with tempfile.TemporaryDirectory(prefix="ophi_rietan_") as temp:
            work = Path(temp)
            ins_path = work / "ophi.ins"
            vesta_command = [
                str(self.vesta_exe),
                "-nogui",
                "-i",
                str(cif),
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
                raise RuntimeError(f"VESTA could not export a RIETAN input file: {detail}")
            ins_text = ins_path.read_text(encoding="utf-8", errors="replace")
            patched = patch_rietan_input(ins_text, settings.two_theta_min, settings.two_theta_max)
            ins_path.write_bytes(patched.encode("utf-8"))

            environment = os.environ.copy()
            environment["RIETAN"] = str(self.rietan_exe.parent) + os.sep
            command = build_rietan_command(self.rietan_exe)
            simulated = subprocess.run(
                command,
                cwd=work,
                env=environment,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=self.timeout_seconds,
            )
            output = f"{simulated.stdout or ''}\n{simulated.stderr or ''}"
            gpd_path = work / "ophi.gpd"
            if simulated.returncode != 0 or not gpd_path.is_file() or re.search(r"error occurred|forrtl: severe", output, re.I):
                detail = "\n".join(line for line in output.splitlines()[-12:] if line.strip())
                raise RuntimeError(f"RIETAN-FP simulation failed: {detail or 'no profile was produced'}")
            reflections = parse_rietan_reflections(output)
            profile = parse_rietan_profile(gpd_path.read_text(encoding="utf-8", errors="replace"))
            peaks = profile_peaks(profile, reflections, intensity_threshold=settings.intensity_threshold)
            if not peaks:
                raise RuntimeError("RIETAN-FP produced a profile but no peaks above the configured threshold.")

        reflection_for_peak = [min(reflections, key=lambda item: abs(item.two_theta - peak.two_theta)) for peak in peaks]
        profile_lookup = {value: index for index, value in enumerate(profile.two_theta_deg)}
        raw = tuple(profile.raw_intensity[profile_lookup[peak.two_theta]] for peak in peaks)
        normalized = tuple(peak.intensity for peak in peaks)
        settings_fingerprint = self.settings_fingerprint(settings)
        version_match = re.search(r"RIETAN-FP\s+v([0-9.]+)", output)
        dependency_version = f"RIETAN-FP {version_match.group(1)}" if version_match else "RIETAN-FP unknown"
        pattern_fingerprint = _fingerprint(
            {
                "settings": settings_fingerprint,
                "profile_two_theta_deg": profile.two_theta_deg,
                "profile_normalized_intensity": profile.normalized_intensity,
            }
        )
        return SimulatedPattern(
            structure_id=context.structure_id,
            source=context.source,
            cif_path=str(cif),
            cif_sha256=SimulatedPattern.calculate_cif_sha256(cif),
            formula=context.formula,
            space_group_number=context.space_group_number,
            radiation="CuKalpha12",
            wavelength_angstrom=RIETAN_CU_KALPHA1_ANGSTROM,
            two_theta_min_deg=settings.two_theta_min,
            two_theta_max_deg=settings.two_theta_max,
            engine_name=self.name,
            engine_version=dependency_version,
            dependency_version=dependency_version,
            settings_fingerprint=settings_fingerprint,
            pattern_fingerprint=pattern_fingerprint,
            two_theta_deg=tuple(peak.two_theta for peak in peaks),
            d_spacing_angstrom=tuple(peak.d_spacing for peak in peaks),
            hkl=tuple(peak.hkl for peak in peaks),
            multiplicity=tuple(peak.multiplicity for peak in peaks),
            line_component=tuple(item.line_component for item in reflection_for_peak),
            raw_intensity=raw,
            normalized_intensity=normalized,
            profile_two_theta_deg=profile.two_theta_deg,
            profile_normalized_intensity=profile.normalized_intensity,
        )


def _fingerprint(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()
