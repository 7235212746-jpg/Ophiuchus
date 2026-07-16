from __future__ import annotations

import json
import os
import queue
import re
import shutil
import subprocess
import threading
import time
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

from .help_support import CONTACT_EMAIL, open_contact_email, open_manual
from .library.analysis import run_library_analysis
from .library.candidate_service import CandidateStructureService
from .library.database import StructureLibrary
from .library.inspector import inspect_peak
from .library.local_provider import LocalFolderProvider
from .library.mp_provider import MaterialsProjectProvider
from .library.providers import ProviderError
from .library.systems import generate_chemical_systems
from .library.xrd_cache import build_library_xrd_cache, scoped_structure_ids, simulation_settings_hash
from .periodic_table import (
    element_scope_mismatch,
    infer_elements_from_xrd_path,
    parse_element_symbols,
    show_periodic_table,
)
from .paths import desktop_dir, first_existing_directory
from .phase_stripping.window import PhaseStrippingWindow, open_or_raise_phase_stripping_window
from .refinement.window import RefinementWindow, open_or_raise_refinement_window
from .refinement.oxide_candidates import supplement_common_oxide_library
from .session_storage import TransientAnalysisStore
from .theme import COLORS, FONTS, SPACING
from .xrd.cache import CandidateCache, build_candidate_cache
from .xrd.candidates import structure_hash
from .xrd.pipeline import run_analysis
from .xrd.vesta_reference import find_local_vesta_references
from .xrd.rietan_backend import discover_rietan_executable
from .xrd.refinement import RietanRefinementBackend


CANDIDATE_EXTENSIONS = {".cif", ".int", ".txt", ".xy", ".dat", ".csv"}


def workbench_sections() -> list[dict[str, str]]:
    return [
        {"key": "projects", "label": "Projects"},
        {"key": "samples", "label": "Samples"},
        {"key": "library", "label": "Library"},
        {"key": "xrd", "label": "XRD Analysis"},
        {"key": "phase", "label": "Phase Evidence"},
        {"key": "settings", "label": "Settings"},
    ]


def project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def default_candidate_dir(root: Path | None = None) -> Path:
    desktop = desktop_dir()
    candidates = (
        desktop / "03_实验数据与分析" / "结构与CIF" / "结构",
        desktop / "结构",
        (root or project_root()) / "data" / "library",
        (root or project_root()) / "data",
    )
    return next((path for path in candidates if path.is_dir()), desktop)


def default_xrd_dir() -> Path:
    desktop = desktop_dir()
    candidates = (
        desktop / "03_实验数据与分析" / "XRD与研相" / "XRD",
        desktop / "XRD",
        desktop,
    )
    return next(path for path in candidates if path.is_dir())


def initial_window_geometry(screen_width: int, screen_height: int) -> str:
    width = min(1240, max(980, screen_width - 80))
    height = min(700, max(620, int(screen_height * 0.78)))
    left = max(0, (screen_width - width) // 2)
    top = max(0, (screen_height - height) // 2)
    return f"{width}x{height}+{left}+{top}"


def default_cache_path(root: Path | None = None) -> Path:
    base = root or project_root()
    return base / "data" / "ophi_xrd_cache.sqlite"


def default_library_path(root: Path | None = None) -> Path:
    base = root or project_root()
    return base / "data" / "ophi_library.sqlite"


def default_env_path(root: Path | None = None) -> Path:
    base = root or project_root()
    return base / ".env"


def default_app_state_path(root: Path | None = None) -> Path:
    base = root or project_root()
    return base / "data" / "ophi_app_state.json"


APP_STATE_KEYS = {
    "xrd_file",
    "candidate_dir",
    "elements",
    "extra_elements",
    "out_dir",
    "cache_path",
    "library_path",
    "target_phase_label",
    "harvest_mode",
    "database_provider",
    "database_endpoint",
    "vesta_exe",
    "rietan_exe",
    "vesta_reference_dir",
    "scientific_safe_mode",
}


def load_app_state(state_path: str | Path) -> dict[str, str]:
    path = Path(state_path)
    if not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
    except Exception:
        return {}
    if not isinstance(payload, dict):
        return {}
    return {key: str(value) for key, value in payload.items() if key in APP_STATE_KEYS and value is not None}


def save_app_state(state_path: str | Path, state: dict[str, object]) -> None:
    path = Path(state_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {key: str(value) for key, value in state.items() if key in APP_STATE_KEYS and value is not None}
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def repair_app_state_paths(state: dict[str, str], root: Path | None = None) -> dict[str, str]:
    """Repair paths saved before the portable project folder was moved."""
    base = root or project_root()
    repaired = dict(state)
    xrd_file = repaired.get("xrd_file", "").strip()
    if xrd_file and not Path(xrd_file).is_file():
        repaired["xrd_file"] = ""

    candidate_dir = repaired.get("candidate_dir", "").strip()
    if candidate_dir and not Path(candidate_dir).is_dir():
        repaired["candidate_dir"] = str(default_candidate_dir(base))

    out_dir = repaired.get("out_dir", "").strip()
    if not out_dir or not Path(out_dir).is_dir():
        suffix = Path(out_dir).name if out_dir else ""
        repaired["out_dir"] = (
            str(base / "results" / suffix)
            if suffix and suffix.lower() != "results"
            else str(base / "results")
        )

    cache_path = repaired.get("cache_path", "").strip()
    if not cache_path or not Path(cache_path).is_file():
        repaired["cache_path"] = str(default_cache_path(base))

    library_path = repaired.get("library_path", "").strip()
    if not library_path or not Path(library_path).is_file():
        repaired["library_path"] = str(default_library_path(base))

    reference_dir = repaired.get("vesta_reference_dir", "").strip()
    if reference_dir and not Path(reference_dir).is_dir():
        repaired["vesta_reference_dir"] = str(default_xrd_dir())
    for key in ("vesta_exe", "rietan_exe"):
        configured = repaired.get(key, "").strip()
        if configured and not Path(configured).is_file():
            repaired.pop(key, None)
    return repaired


def read_local_env(env_path: str | Path) -> dict[str, str]:
    path = Path(env_path)
    if not path.exists():
        return {}
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        clean = line.strip()
        if not clean or clean.startswith("#") or "=" not in clean:
            continue
        key, value = clean.split("=", 1)
        values[key.strip()] = value.strip().strip("'\"")
    return values


def update_local_env(env_path: str | Path, updates: dict[str, str]) -> None:
    path = Path(env_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    existing_lines = path.read_text(encoding="utf-8", errors="ignore").splitlines() if path.exists() else []
    remaining = dict(updates)
    output: list[str] = []
    for line in existing_lines:
        clean = line.strip()
        if not clean or clean.startswith("#") or "=" not in clean:
            output.append(line)
            continue
        key = clean.split("=", 1)[0].strip()
        if key in remaining:
            output.append(f"{key}={remaining.pop(key).strip()}")
        else:
            output.append(line)
    for key, value in remaining.items():
        output.append(f"{key}={value.strip()}")
    path.write_text("\n".join(output).rstrip() + "\n", encoding="utf-8")


def save_mp_api_key_to_env(env_path: str | Path, api_key: str) -> None:
    if api_key.strip():
        update_local_env(env_path, {"MP_API_KEY": api_key.strip()})


def normalize_database_provider(provider: str) -> str:
    value = provider.strip().lower()
    aliases = {
        "materials project": "materials_project",
        "mp": "materials_project",
        "materials_project": "materials_project",
        "cod": "cod",
        "aflow": "aflow",
        "oqmd": "oqmd",
        "nomad/optimade": "nomad_optimade",
        "nomad_optimade": "nomad_optimade",
    }
    return aliases.get(value, value.replace(" ", "_"))


def save_database_api_config(env_path: str | Path, provider: str, api_key: str, endpoint: str = "") -> dict[str, object]:
    provider_id = normalize_database_provider(provider)
    current = read_local_env(env_path)
    updates = {"OPHI_STRUCTURE_DATABASE": provider_id}
    effective_key = api_key.strip() or current.get("MP_API_KEY", "")
    if provider_id == "materials_project":
        if effective_key:
            updates["MP_API_KEY"] = effective_key
    if endpoint.strip():
        updates["OPHI_STRUCTURE_DATABASE_ENDPOINT"] = endpoint.strip()
    elif current.get("OPHI_STRUCTURE_DATABASE_ENDPOINT"):
        updates["OPHI_STRUCTURE_DATABASE_ENDPOINT"] = current["OPHI_STRUCTURE_DATABASE_ENDPOINT"]
    update_local_env(env_path, updates)
    return {"configured": bool(effective_key) if provider_id == "materials_project" else False, "provider": provider_id}


def load_database_api_config(env_path: str | Path) -> dict[str, str]:
    values = read_local_env(env_path)
    return {
        "provider": values.get("OPHI_STRUCTURE_DATABASE", "materials_project"),
        "api_key": values.get("MP_API_KEY") or values.get("PMG_MAPI_KEY") or values.get("MAPI_KEY") or "",
        "endpoint": values.get("OPHI_STRUCTURE_DATABASE_ENDPOINT", "https://api.materialsproject.org"),
    }


def find_vesta_executable(env_path: str | Path | None = None) -> str:
    values = read_local_env(env_path or default_env_path())
    configured = values.get("OPHI_VESTA_EXE", "")
    if configured and Path(configured).exists():
        return configured
    desktop = desktop_dir()
    candidates = [
        Path(os.environ.get("PROGRAMFILES", r"C:\Program Files")) / "VESTA-win64" / "VESTA.exe",
        Path(os.environ.get("PROGRAMFILES", r"C:\Program Files")) / "VESTA" / "VESTA.exe",
        Path(os.environ.get("PROGRAMFILES(X86)", r"C:\Program Files (x86)")) / "VESTA-win64" / "VESTA.exe",
        desktop / "02_科研软件" / "便携软件" / "VESTA-win64" / "VESTA.exe",
        desktop / "VESTA.exe",
        desktop / "VESTA-win64" / "VESTA.exe",
        Path.home() / "OneDrive" / "Desktop" / "VESTA.exe",
        Path.home() / "OneDrive" / "Desktop" / "VESTA-win64" / "VESTA.exe",
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    return ""


def load_vesta_config(env_path: str | Path | None = None) -> dict[str, str]:
    path = env_path or default_env_path()
    apply_vesta_env_config(path)
    values = read_local_env(path)
    configured_reference = values.get("OPHI_VESTA_REFERENCE_DIR", "").strip()
    reference_dir = Path(configured_reference) if configured_reference else None
    return {
        "vesta_exe": find_vesta_executable(path),
        "rietan_exe": str(discover_rietan_executable(values.get("OPHI_RIETAN_EXE", "")) or ""),
        "reference_dir": str(reference_dir if reference_dir and reference_dir.is_dir() else default_xrd_dir()),
    }


def apply_vesta_env_config(env_path: str | Path | None = None) -> None:
    values = read_local_env(env_path or default_env_path())
    for key, value in values.items():
        if key in {"OPHI_VESTA_EXE", "OPHI_RIETAN_EXE", "OPHI_VESTA_REFERENCE_DIR"} or key.startswith("OPHI_VESTA_REFERENCE_"):
            if value and Path(value).exists():
                os.environ[key] = value
            else:
                os.environ.pop(key, None)


def save_vesta_config(
    env_path: str | Path,
    vesta_exe: str,
    reference_dir: str,
    rietan_exe: str = "",
) -> dict[str, object]:
    updates: dict[str, str] = {}
    if vesta_exe.strip():
        updates["OPHI_VESTA_EXE"] = vesta_exe.strip()
    if reference_dir.strip():
        updates["OPHI_VESTA_REFERENCE_DIR"] = reference_dir.strip()
    if rietan_exe.strip():
        updates["OPHI_RIETAN_EXE"] = rietan_exe.strip()
    if updates:
        update_local_env(env_path, updates)
    return {
        "vesta_exe_exists": bool(vesta_exe.strip()) and Path(vesta_exe.strip()).exists(),
        "rietan_exe_exists": bool(rietan_exe.strip()) and Path(rietan_exe.strip()).exists(),
        "reference_dir_exists": bool(reference_dir.strip()) and Path(reference_dir.strip()).exists(),
    }


def save_formula_vesta_reference(env_path: str | Path, formula: str, reference_file: str) -> dict[str, object]:
    clean_formula = re.sub(r"[^A-Za-z0-9]", "", formula.strip()).upper()
    path = Path(reference_file.strip())
    if not clean_formula:
        raise ValueError("formula is required")
    if not reference_file.strip():
        raise ValueError("reference file is required")
    key = f"OPHI_VESTA_REFERENCE_{clean_formula}"
    update_local_env(env_path, {key: str(path)})
    os.environ[key] = str(path)
    return {"key": key, "reference_exists": path.exists() and path.is_file()}


def launch_vesta(vesta_exe: str, target_file: str = "") -> None:
    exe = Path(vesta_exe)
    if not exe.exists():
        raise FileNotFoundError(f"VESTA executable not found: {vesta_exe}")
    args = [str(exe)]
    if target_file.strip():
        target = Path(target_file)
        if not target.exists():
            raise FileNotFoundError(f"VESTA target file not found: {target_file}")
        args.append(str(target))
    subprocess.Popen(args)


def assess_vesta_preflight(vesta_exe: str, reference_dir: str, formulas: list[str]) -> dict[str, object]:
    clean_formulas = []
    for formula in formulas:
        clean = str(formula).strip()
        if clean and clean not in clean_formulas:
            clean_formulas.append(clean)
    previous = os.environ.get("OPHI_VESTA_REFERENCE_DIR")
    if reference_dir.strip():
        os.environ["OPHI_VESTA_REFERENCE_DIR"] = reference_dir.strip()
    try:
        reference_paths: dict[str, str] = {}
        reference_candidate_counts: dict[str, int] = {}
        reference_candidates: dict[str, list[str]] = {}
        missing: list[str] = []
        for formula in clean_formulas:
            references = find_local_vesta_references(formula)
            reference_candidate_counts[formula] = len(references)
            reference_candidates[formula] = [str(item["pattern"]) for item in references[:5]]
            if references:
                reference_paths[formula] = str(references[0]["pattern"])
            else:
                missing.append(formula)
    finally:
        if previous is None:
            os.environ.pop("OPHI_VESTA_REFERENCE_DIR", None)
        else:
            os.environ["OPHI_VESTA_REFERENCE_DIR"] = previous
    return {
        "vesta_exe_exists": bool(vesta_exe.strip()) and Path(vesta_exe.strip()).exists(),
        "reference_dir_exists": bool(reference_dir.strip()) and Path(reference_dir.strip()).exists(),
        "formulas_checked": len(clean_formulas),
        "references_found": len(reference_paths),
        "missing_references": missing,
        "reference_paths": reference_paths,
        "reference_candidate_counts": reference_candidate_counts,
        "reference_candidates": reference_candidates,
    }


def vesta_preflight_lines(status: dict[str, object]) -> list[str]:
    formulas_checked = int(status.get("formulas_checked") or 0)
    references_found = int(status.get("references_found") or 0)
    missing = [str(item) for item in status.get("missing_references", [])]
    lines = [
        "VESTA 预检：",
        f"- VESTA 程序：{'已找到' if status.get('vesta_exe_exists') else '未找到'}",
        f"- VESTA 参考目录：{'已找到' if status.get('reference_dir_exists') else '未找到'}",
        f"- 参考峰覆盖：{references_found}/{formulas_checked}",
    ]
    if missing:
        lines.append("- 缺少 VESTA 参考：" + ", ".join(missing[:8]))
    counts = status.get("reference_candidate_counts") or {}
    multi = [f"{formula}({count})" for formula, count in counts.items() if int(count) > 1]
    if multi:
        lines.append("- 多个 VESTA 备选：" + ", ".join(multi[:8]))
        candidates = status.get("reference_candidates") or {}
        for formula in [item.split("(", 1)[0] for item in multi[:3]]:
            paths = [str(path) for path in candidates.get(formula, [])[:3]]
            if paths:
                lines.append(f"  {formula} 备选：")
                lines.extend(f"  - {path}" for path in paths)
    return lines


def materials_project_harvest_from_app(
    library_db: str | Path,
    elements: str,
    impurities: str,
    mode: str,
    api_key: str = "",
    env_path: str | Path | None = None,
    dry_run: bool = False,
    max_entries_per_system: int = 50,
    progress_callback=None,
) -> dict[str, object]:
    systems = generate_chemical_systems(elements.split(), impurities.split(), mode=mode, max_systems=80)
    if dry_run:
        return {"chemical_systems": systems, "count": len(systems), "mode": mode}
    provider = MaterialsProjectProvider(api_key=api_key or None, env_path=env_path)
    library = StructureLibrary(library_db)
    summary = provider.harvest_to_library(
        library,
        systems,
        max_entries_per_system=max_entries_per_system,
        progress_callback=progress_callback,
    )
    cache_summary = build_library_xrd_cache(library)
    return {**summary, "chemical_systems": systems, "xrd_cache": cache_summary}


def import_local_cifs_to_library(folder: str | Path, library_db: str | Path) -> dict[str, object]:
    library = StructureLibrary(library_db)
    return LocalFolderProvider().import_folder(folder, library)


def import_target_cif_to_library(cif_file: str | Path, library_db: str | Path) -> dict[str, object]:
    path = Path(cif_file)
    if not path.is_file() or path.suffix.lower() != ".cif":
        raise FileNotFoundError(f"target phase CIF does not exist: {path}")
    library = StructureLibrary(library_db)
    provider = LocalFolderProvider()
    entry = provider._entry_from_cif(
        path,
        library,
        source="target",
        access_note="User-selected target phase CIF; verify source license before reuse.",
    )
    target = library.path.parent / entry.cached_file_path
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists() or structure_hash(target) != entry.structure_hash:
        shutil.copy2(path, target)
    library.upsert_structure(entry)
    formula = entry.reduced_formula or entry.formula
    return {
        "scanned": 1,
        "imported": 1,
        "skipped": 0,
        "internal_id": entry.internal_id,
        "formula": formula,
        "elements": entry.elements,
        "label": f"{formula} | {entry.source} | {entry.internal_id}",
        "cached_file_path": str(target),
    }


def build_library_cache_from_folder(library_db: str | Path) -> dict[str, object]:
    library = StructureLibrary(library_db)
    return build_library_xrd_cache(library)


def library_manager_rows(
    library_db: str | Path,
    radiation: str = "CuKa",
    two_theta_range: tuple[float, float] = (10.0, 90.0),
) -> list[dict[str, object]]:
    library = StructureLibrary(library_db)
    service = CandidateStructureService(library, radiation=radiation, two_theta_range=two_theta_range, candidate_scope="subsystems")
    rows = []
    for row in service.list_candidate_rows([]):
        rows.append(
            {
                "internal_id": row.internal_id,
                "formula": row.formula,
                "chemical_system": row.chemical_system,
                "source": row.source,
                "source_id": row.source_id,
                "enabled": row.enabled,
                "xrd_cache": row.cache_status,
                "space_group_symbol": row.space_group_symbol,
                "space_group_number": row.space_group_number,
                "backend_name": row.backend_name,
                "backend_version": row.backend_version,
                "pattern_fingerprint": row.pattern_fingerprint,
                "candidate_ready": row.candidate_for_analysis,
                "skip_reasons": " | ".join(row.skip_reasons),
                "last_simulated_time": row.last_simulated_time,
                "structure_hash": row.structure_hash,
                "access_note": "",
                "cached_file_path": row.structure_path,
                "metadata_path": "",
                "original_file_path": row.structure_path,
            }
        )
    return rows


def set_library_entry_enabled(library_db: str | Path, internal_id: str, enabled: bool) -> None:
    StructureLibrary(library_db).set_enabled(internal_id, enabled)


def inspect_peak_from_app(
    library_db: str | Path,
    experimental_two_theta: float,
    tolerance_deg: float = 0.2,
    radiation: str = "CuKa",
    two_theta_range: tuple[float, float] = (10.0, 90.0),
) -> dict[str, object]:
    settings_hash = simulation_settings_hash(radiation, two_theta_range)
    return inspect_peak(
        StructureLibrary(library_db),
        experimental_two_theta=experimental_two_theta,
        tolerance_deg=tolerance_deg,
        radiation=radiation,
        settings_hash=settings_hash,
    )


def run_library_analysis_from_app(
    xrd_file: str | Path,
    library_db: str | Path,
    elements: str,
    extra_elements: str,
    out_dir: str | Path,
    target_candidate_id: str = "",
    scientific_safe_mode: bool = True,
    force_recompute: bool = False,
    use_rietan_display: bool = False,
):
    return run_library_analysis(
        xrd_file=xrd_file,
        library_db=library_db,
        elements=elements.split(),
        extra_elements=[],
        out_dir=out_dir,
        project="Ophiuchus",
        sample_id=Path(xrd_file).stem,
        target_candidate_id=target_candidate_id or None,
        scientific_safe_mode=scientific_safe_mode,
        force_recompute=force_recompute,
        use_rietan_display=use_rietan_display,
    )


def library_target_phase_options(library_db: str | Path, elements: str) -> list[dict[str, str]]:
    library = StructureLibrary(library_db)
    target_elements = set(elements.split())
    rows: list[dict[str, str]] = []
    for entry in library.list_structures(enabled_only=True):
        entry_elements = set(entry.elements)
        if target_elements and entry_elements != target_elements:
            continue
        formula = entry.reduced_formula or entry.formula
        label = f"{formula} | {entry.source} | {entry.internal_id}"
        rows.append(
            {
                "label": label,
                "id": entry.internal_id,
                "formula": formula,
                "source": entry.source,
                "element_count": str(len(entry_elements)),
                "is_exact": "1" if target_elements and entry_elements == target_elements else "0",
                "is_local": "1" if entry.source == "local" else "0",
            }
        )
    rows.sort(key=lambda row: (row["is_exact"] != "1", row["is_local"] != "1", row["formula"], row["id"]))
    return rows


def resolve_target_phase_selection(rows: list[dict[str, str]], current_label: str) -> tuple[str, str, dict[str, str]]:
    options = {row["label"]: row["id"] for row in rows}
    labels = list(options)
    if not labels:
        return "", "", {}
    if current_label not in options:
        return "", "", options
    return current_label, options[current_label], options


def build_output_dir(xrd_file: str | Path, root: Path | None = None) -> Path:
    base = root or project_root()
    stem = Path(xrd_file).stem if str(xrd_file).strip() else "ophi_xrd"
    safe = re.sub(r"[^A-Za-z0-9\u4e00-\u9fff]+", "_", stem).strip("_") or "ophi_xrd"
    return base / "results" / safe


def validate_analysis_inputs(xrd_file: str, candidate_dir: str, elements: str, out_dir: str) -> list[str]:
    errors: list[str] = []
    if not xrd_file.strip():
        errors.append("请选择实验 XRD 文件。")
    elif not Path(xrd_file).is_file():
        errors.append("实验 XRD 文件不存在，请重新选择。")

    if not candidate_dir.strip():
        errors.append("请选择候选结构/模拟峰文件夹。")
    elif not Path(candidate_dir).is_dir():
        errors.append("候选文件夹不存在，请重新选择。")
    elif not any(path.is_file() and path.suffix.lower() in CANDIDATE_EXTENSIONS for path in Path(candidate_dir).rglob("*")):
        errors.append("候选文件夹里没有找到 CIF 或峰表文件。")

    if not elements.strip():
        errors.append("请填写主元素，例如 Zr Fe Ge。")
    else:
        try:
            selected_elements = parse_element_symbols(elements)
        except ValueError as exc:
            errors.append(f"主元素包含无效符号：{exc}")
        else:
            inferred_elements = infer_elements_from_xrd_path(xrd_file) if xrd_file.strip() else ()
            if inferred_elements:
                missing, _extra = element_scope_mismatch(selected_elements, inferred_elements)
                if missing:
                    errors.append(
                        "实验谱名称推断出的元素没有全部包含在主元素范围中："
                        + " ".join(missing)
                        + "。请检查元素周期表或文件名。"
                    )

    return errors


class OphiuchusApp(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("Ophiuchus")
        self.geometry(initial_window_geometry(self.winfo_screenwidth(), self.winfo_screenheight()))
        self.minsize(1000, 620)
        self.configure(bg=COLORS["background"])
        self.result_queue: queue.Queue[tuple[str, object]] = queue.Queue()
        self.running = False
        self.analysis_store = TransientAnalysisStore()
        self.current_analysis_result = None
        self.phase_stripping_window: PhaseStrippingWindow | None = None
        self.refinement_window: RefinementWindow | None = None
        self._vesta_dialog: tk.Toplevel | None = None
        self._help_dialog: tk.Toplevel | None = None
        self.last_saved_analysis_path: Path | None = None
        self._analysis_run_id = 0
        self._closed = False
        self._poll_after_id: str | None = None
        self._build_style()
        self._build_layout()
        self._load_saved_database_config()
        self._restore_app_state()
        self._refresh_target_phase_options()
        self.protocol("WM_DELETE_WINDOW", self._on_close)
        self._schedule_result_poll()

    def _build_style(self) -> None:
        style = ttk.Style(self)
        style.theme_use("clam")
        style.configure("Shell.TFrame", background=COLORS["background"])
        style.configure("Panel.TFrame", background=COLORS["glass"], relief="flat", borderwidth=1)
        style.configure("Subtle.TFrame", background=COLORS["panel_alt"], relief="flat", borderwidth=1)
        style.configure("TLabel", background=COLORS["background"], foreground=COLORS["text"], font=FONTS["body_zh"])
        style.configure("Panel.TLabel", background=COLORS["glass"], foreground=COLORS["text"], font=FONTS["body_zh"])
        style.configure("Subtle.TLabel", background=COLORS["panel_alt"], foreground=COLORS["text"], font=FONTS["body_zh"])
        style.configure("Muted.TLabel", background=COLORS["background"], foreground=COLORS["muted"], font=FONTS["small"])
        style.configure("PanelMuted.TLabel", background=COLORS["glass"], foreground=COLORS["muted"], font=FONTS["small"])
        style.configure("Title.TLabel", background=COLORS["background"], foreground=COLORS["text"], font=FONTS["title"])
        style.configure("Section.TLabel", background=COLORS["glass"], foreground=COLORS["text"], font=FONTS["section"])
        style.configure("Step.TLabel", background=COLORS["glass"], foreground=COLORS["accent"], font=FONTS["button"])
        style.configure("TButton", font=FONTS["button"], padding=(12, 8), background=COLORS["panel_alt"], foreground=COLORS["text"], borderwidth=1, relief="flat")
        style.map("TButton", background=[("active", COLORS["accent_soft"]), ("disabled", COLORS["border_soft"])])
        style.configure("Primary.TButton", font=FONTS["button"], padding=(14, 10), background=COLORS["accent"], foreground="#ffffff", borderwidth=0, relief="flat")
        style.map("Primary.TButton", background=[("active", COLORS["accent_hover"]), ("disabled", COLORS["border"])])
        style.configure("TEntry", fieldbackground=COLORS["input"], foreground=COLORS["text"], insertcolor=COLORS["text"], bordercolor=COLORS["border"], lightcolor=COLORS["border_soft"], darkcolor=COLORS["border"])
        style.configure("TNotebook", background=COLORS["glass"], borderwidth=0)
        style.configure("TNotebook.Tab", background=COLORS["panel_alt"], foreground=COLORS["muted"], padding=(12, 7), font=FONTS["button"])
        style.map("TNotebook.Tab", background=[("selected", COLORS["accent_soft"])], foreground=[("selected", COLORS["text"])])
        style.configure("Treeview", background=COLORS["panel"], fieldbackground=COLORS["panel"], foreground=COLORS["text"], rowheight=28, font=FONTS["body"])
        style.configure("Treeview.Heading", background=COLORS["panel_alt"], foreground=COLORS["muted"], font=FONTS["button"])
        style.configure("Panel.TCheckbutton", background=COLORS["glass"], foreground=COLORS["text"], font=FONTS["body_zh"])
        style.map("Panel.TCheckbutton", background=[("active", COLORS["glass"])])

    def _build_layout(self) -> None:
        root = ttk.Frame(self, style="Shell.TFrame", padding=24)
        root.pack(fill="both", expand=True)

        self._build_header(root)

        body = ttk.Frame(root, style="Shell.TFrame")
        body.pack(fill="both", expand=True, pady=(18, 0))
        body.columnconfigure(0, weight=0, minsize=170)
        body.columnconfigure(1, weight=1, minsize=480)
        body.columnconfigure(2, weight=1, minsize=520)
        body.rowconfigure(0, weight=1)

        self._build_sidebar(body)
        self._build_workflow_panel(body)
        self._build_result_panel(body)

    def _build_sidebar(self, parent: ttk.Frame) -> None:
        self.sidebar = ttk.Frame(parent, style="Panel.TFrame", padding=14)
        self.sidebar.grid(row=0, column=0, sticky="nsew", padx=(0, 16))
        ttk.Label(self.sidebar, text="Workspace", style="Section.TLabel").pack(anchor="w")
        ttk.Label(self.sidebar, text="Local evidence workbench", style="PanelMuted.TLabel").pack(anchor="w", pady=(2, 14))
        self.nav_buttons: dict[str, ttk.Button] = {}
        for section in workbench_sections():
            button = ttk.Button(
                self.sidebar,
                text=section["label"],
                command=lambda key=section["key"]: self._select_section(key),
            )
            button.pack(fill="x", pady=3)
            self.nav_buttons[section["key"]] = button
        footer = ttk.Frame(self.sidebar, style="Panel.TFrame")
        footer.pack(fill="x", side="bottom", pady=(16, 0))
        self.help_button = ttk.Button(footer, text="帮助与关于", command=self._open_help_dialog)
        self.help_button.pack(fill="x", pady=(0, 10))
        ttk.Label(
            footer,
            text="Screening evidence only. No fake wt%.",
            style="PanelMuted.TLabel",
            wraplength=130,
        ).pack(anchor="w")

    def _open_help_dialog(self) -> None:
        if self._help_dialog is not None and self._help_dialog.winfo_exists():
            self._help_dialog.deiconify()
            self._help_dialog.lift()
            self._help_dialog.focus_force()
            return

        dialog = tk.Toplevel(self)
        self._help_dialog = dialog
        dialog.title("帮助与关于")
        dialog.configure(bg=COLORS["background"])
        dialog.transient(self)
        dialog.resizable(False, False)
        dialog.protocol("WM_DELETE_WINDOW", self._close_help_dialog)

        content = ttk.Frame(dialog, style="Shell.TFrame", padding=24)
        content.pack(fill="both", expand=True)
        ttk.Label(content, text="Ophiuchus", style="Title.TLabel").pack(anchor="w")
        ttk.Label(
            content,
            text="本地优先的材料科研工作台",
            style="Muted.TLabel",
        ).pack(anchor="w", pady=(2, 18))
        ttk.Label(
            content,
            text=(
                "用于 XRD 候选相筛选、结构库证据整理、残差谱分析与受约束精修。\n"
                "筛选结果是科研证据，不替代人工晶体学判断或完整定量精修。"
            ),
            style="TLabel",
            justify="left",
        ).pack(anchor="w")
        ttk.Label(content, text="维护与反馈", style="Section.TLabel").pack(anchor="w", pady=(20, 4))
        ttk.Label(content, text=CONTACT_EMAIL, style="TLabel").pack(anchor="w")

        self.help_status_var = tk.StringVar(value="完整说明见操作手册。")
        ttk.Label(content, textvariable=self.help_status_var, style="Muted.TLabel").pack(
            anchor="w", pady=(12, 0)
        )

        actions = ttk.Frame(content, style="Shell.TFrame")
        actions.pack(fill="x", pady=(20, 0))
        ttk.Button(actions, text="打开操作手册", command=self._open_user_manual).pack(side="left")
        ttk.Button(actions, text="复制邮箱", command=self._copy_contact_email).pack(side="left", padx=(8, 0))
        ttk.Button(actions, text="发送邮件", command=self._send_contact_email).pack(side="left", padx=(8, 0))
        ttk.Button(actions, text="关闭", command=self._close_help_dialog).pack(side="right")

        dialog.update_idletasks()
        width = max(600, dialog.winfo_reqwidth())
        height = max(330, dialog.winfo_reqheight())
        x = self.winfo_rootx() + max(30, (self.winfo_width() - width) // 2)
        y = self.winfo_rooty() + max(30, (self.winfo_height() - height) // 2)
        dialog.geometry(f"{width}x{height}+{x}+{y}")
        dialog.lift()
        dialog.focus_force()

    def _open_user_manual(self) -> None:
        try:
            open_manual()
        except (FileNotFoundError, OSError) as exc:
            messagebox.showerror("无法打开操作手册", str(exc), parent=self._help_dialog or self)
            return
        self.help_status_var.set("已使用系统默认程序打开操作手册。")

    def _copy_contact_email(self) -> None:
        self.clipboard_clear()
        self.clipboard_append(CONTACT_EMAIL)
        self.update_idletasks()
        self.help_status_var.set("邮箱已复制到剪贴板。")

    def _send_contact_email(self) -> None:
        try:
            open_contact_email()
        except OSError as exc:
            messagebox.showerror("无法打开邮件程序", str(exc), parent=self._help_dialog or self)
            return
        self.help_status_var.set("已调用系统默认邮件程序。")

    def _close_help_dialog(self) -> None:
        if self._help_dialog is not None and self._help_dialog.winfo_exists():
            self._help_dialog.destroy()
        self._help_dialog = None

    def _build_header(self, root: ttk.Frame) -> None:
        ttk.Label(root, text="Ophiuchus", style="Title.TLabel").pack(anchor="w")
        ttk.Label(
            root,
            text="你的本地材料科研流程助手。先从 XRD 候选相筛选开始，帮你判断哪些相值得查、哪些强峰缺失、哪些峰仍未解释。",
            style="Muted.TLabel",
        ).pack(anchor="w", pady=(3, 0))

    def _build_workflow_panel(self, parent: ttk.Frame) -> None:
        self.workspace_panel = ttk.Frame(parent, style="Panel.TFrame", padding=18)
        panel = self.workspace_panel
        panel.grid(row=0, column=1, sticky="nsew", padx=(0, 16))

        ttk.Label(panel, text="本次任务", style="Section.TLabel").pack(anchor="w")
        ttk.Label(
            panel,
            text="按顺序选择实验谱、候选来源和元素范围。Ophi 会生成候选排序和透明报告。",
            style="PanelMuted.TLabel",
        ).pack(anchor="w", pady=(4, 16))

        self.xrd_var = tk.StringVar()
        self.dir_var = tk.StringVar(value=str(default_candidate_dir()))
        self.elements_var = tk.StringVar(value="Zr Fe Ge")
        self.extra_var = tk.StringVar(value="O C Si Al Cu Sn Hf Sc I")
        self.out_var = tk.StringVar(value=str(project_root() / "results"))
        self.cache_var = tk.StringVar(value=str(default_cache_path()))
        self.library_var = tk.StringVar(value=str(default_library_path()))
        self.target_phase_var = tk.StringVar()
        self.scientific_safe_mode_var = tk.BooleanVar(value=True)
        self.force_recompute_var = tk.BooleanVar(value=False)
        self.target_phase_options: dict[str, str] = {}
        self.mp_api_key_var = tk.StringVar()
        self.harvest_mode_var = tk.StringVar(value="normal")
        self.database_provider_var = tk.StringVar(value="Materials Project")
        self.database_endpoint_var = tk.StringVar(value="https://api.materialsproject.org")
        vesta_config = load_vesta_config(default_env_path())
        self.vesta_exe_var = tk.StringVar(value=vesta_config["vesta_exe"])
        self.rietan_exe_var = tk.StringVar(value=vesta_config["rietan_exe"])
        self.vesta_reference_dir_var = tk.StringVar(value=vesta_config["reference_dir"])
        self.vesta_formula_var = tk.StringVar()
        self.vesta_formula_reference_var = tk.StringVar()
        self.current_plot_path: str | None = None
        self.plot_image = None

        action_row = ttk.Frame(panel, style="Panel.TFrame")
        action_row.pack(fill="x", pady=(0, 14))
        for column in range(3):
            action_row.columnconfigure(column, weight=1, uniform="workflow_actions")
        self.run_button = ttk.Button(action_row, text="开始候选筛选", style="Primary.TButton", command=self._run)
        self.run_button.grid(row=0, column=0, sticky="ew", padx=(0, 5), pady=(0, 6))
        self.database_button = ttk.Button(action_row, text="连接 API / 数据库", command=self._open_database_dialog)
        self.database_button.grid(row=0, column=1, sticky="ew", padx=5, pady=(0, 6))
        self.vesta_button = ttk.Button(action_row, text="VESTA / RIETAN 设置", command=self._open_vesta_dialog)
        self.vesta_button.grid(row=0, column=2, sticky="ew", padx=5, pady=(0, 6))
        self.cache_button = ttk.Button(action_row, text="构建/更新缓存", command=self._build_cache)
        self.cache_button.grid(row=1, column=0, sticky="ew", padx=(0, 5), pady=(0, 6))
        self.library_button = ttk.Button(action_row, text="导入结构库", command=self._import_library)
        self.library_button.grid(row=1, column=1, sticky="ew", padx=5, pady=(0, 6))
        self.library_run_button = ttk.Button(action_row, text="结构库分析", command=self._run_library)
        self.library_run_button.grid(row=1, column=2, sticky="ew", padx=(5, 0), pady=(0, 6))
        self.save_analysis_button = ttk.Button(action_row, text="保存本次分析", command=self._save_current_analysis)
        self.save_analysis_button.grid(row=2, column=0, sticky="ew", padx=(0, 5), pady=(0, 6))
        self.phase_stripping_button = ttk.Button(
            action_row,
            text="手动物相剥离 / 残差谱分析",
            command=self._open_phase_stripping,
        )
        self.phase_stripping_button.grid(row=2, column=1, sticky="ew", padx=5, pady=(0, 6))
        self.refinement_button = ttk.Button(
            action_row,
            text="RIETAN 受约束精修",
            command=self._open_refinement,
        )
        self.refinement_button.grid(row=2, column=2, sticky="ew", padx=(5, 0), pady=(0, 6))
        self.save_analysis_button.state(["disabled"])
        self.phase_stripping_button.state(["disabled"])
        self.refinement_button.state(["disabled"])
        ttk.Button(action_row, text="打开保存位置", command=self._open_output).grid(
            row=3, column=0, columnspan=3, sticky="ew"
        )

        trust_row = ttk.Frame(panel, style="Panel.TFrame")
        trust_row.pack(fill="x", pady=(0, 12))
        ttk.Checkbutton(
            trust_row,
            text="Scientific Safe Mode（CIF 每次重新计算）",
            variable=self.scientific_safe_mode_var,
            style="Panel.TCheckbutton",
        ).pack(side="left")
        ttk.Checkbutton(
            trust_row,
            text="本次强制重算",
            variable=self.force_recompute_var,
            style="Panel.TCheckbutton",
        ).pack(side="left", padx=(16, 0))

        scroll_shell = ttk.Frame(panel, style="Panel.TFrame")
        scroll_shell.pack(fill="both", expand=True)
        self.workflow_canvas = tk.Canvas(scroll_shell, bg=COLORS["glass"], bd=0, highlightthickness=0)
        self.workflow_scrollbar = ttk.Scrollbar(scroll_shell, orient="vertical", command=self.workflow_canvas.yview)
        self.workflow_canvas.configure(yscrollcommand=self.workflow_scrollbar.set)
        self.workflow_canvas.pack(side="left", fill="both", expand=True)
        self.workflow_scrollbar.pack(side="right", fill="y", padx=(8, 0))
        self.workflow_content = ttk.Frame(self.workflow_canvas, style="Panel.TFrame")
        self.workflow_window_id = self.workflow_canvas.create_window((0, 0), window=self.workflow_content, anchor="nw")
        self.workflow_content.bind("<Configure>", self._sync_workflow_scrollregion)
        self.workflow_canvas.bind("<Configure>", self._sync_workflow_width)
        self.workflow_canvas.bind("<Enter>", self._bind_workflow_mousewheel)
        self.workflow_canvas.bind("<Leave>", self._unbind_workflow_mousewheel)

        steps = self.workflow_content
        self._file_step(steps, "1", "选择实验 XRD", "支持 Rigaku .asc 和常见两列文本。", self.xrd_var, self._pick_xrd)
        self._file_step(steps, "2", "选择候选结构/模拟峰文件夹", "建议先用桌面的 结构 文件夹，也可以选含 .int/.csv 的实验目录。", self.dir_var, self._pick_dir)
        self._element_step(steps)
        self._entry_step(steps, "", "可能杂质/替代元素", "额外元素", self.extra_var)
        self._file_step(steps, "4", "本地候选缓存", "可先构建缓存；正式筛选会自动复用已模拟的峰表。", self.cache_var, self._pick_cache)
        self._file_step(steps, "5", "本地结构库", "Phase 2 结构库数据库；导入本地 CIF 后可构建库级 XRD 缓存。", self.library_var, self._pick_library)
        self._target_phase_step(steps)
        self._harvest_step(steps)
        self._file_step(
            steps,
            "8",
            "选择默认保存位置",
            "分析默认只保留一个临时会话；点击“保存本次分析”后才会复制到这里。",
            self.out_var,
            self._pick_out,
        )

    def _sync_workflow_scrollregion(self, _event=None) -> None:
        self.workflow_canvas.configure(scrollregion=self.workflow_canvas.bbox("all"))

    def _sync_workflow_width(self, event) -> None:
        self.workflow_canvas.itemconfigure(self.workflow_window_id, width=event.width)

    def _bind_workflow_mousewheel(self, _event=None) -> None:
        self.workflow_canvas.bind_all("<MouseWheel>", self._on_workflow_mousewheel)

    def _unbind_workflow_mousewheel(self, _event=None) -> None:
        self.workflow_canvas.unbind_all("<MouseWheel>")

    def _on_workflow_mousewheel(self, event) -> None:
        self.workflow_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

    def _load_saved_database_config(self) -> None:
        config = load_database_api_config(default_env_path())
        provider_labels = {
            "materials_project": "Materials Project",
            "cod": "COD",
            "aflow": "AFLOW",
            "oqmd": "OQMD",
            "nomad_optimade": "NOMAD/OPTIMADE",
        }
        self.database_provider_var.set(provider_labels.get(config["provider"], "Materials Project"))
        self.database_endpoint_var.set(config["endpoint"] or "https://api.materialsproject.org")
        if config["api_key"]:
            self.mp_api_key_var.set(config["api_key"])
            self.status_var.set("已加载本地 Materials Project API 配置。")

    def _restore_app_state(self) -> None:
        state = repair_app_state_paths(load_app_state(default_app_state_path()))
        if not state:
            return
        mapping = {
            "xrd_file": self.xrd_var,
            "candidate_dir": self.dir_var,
            "elements": self.elements_var,
            "extra_elements": self.extra_var,
            "out_dir": self.out_var,
            "cache_path": self.cache_var,
            "library_path": self.library_var,
            "target_phase_label": self.target_phase_var,
            "harvest_mode": self.harvest_mode_var,
            "database_provider": self.database_provider_var,
            "database_endpoint": self.database_endpoint_var,
            "vesta_exe": self.vesta_exe_var,
            "rietan_exe": self.rietan_exe_var,
            "vesta_reference_dir": self.vesta_reference_dir_var,
        }
        for key, var in mapping.items():
            if key in state:
                var.set(state[key])
        if "scientific_safe_mode" in state:
            self.scientific_safe_mode_var.set(state["scientific_safe_mode"].strip().lower() in {"1", "true", "yes", "on"})
        self.status_var.set("已恢复上次关闭时的工作界面。")

    def _collect_app_state(self) -> dict[str, str]:
        return {
            "xrd_file": self.xrd_var.get(),
            "candidate_dir": self.dir_var.get(),
            "elements": self.elements_var.get(),
            "extra_elements": self.extra_var.get(),
            "out_dir": self.out_var.get(),
            "cache_path": self.cache_var.get(),
            "library_path": self.library_var.get(),
            "target_phase_label": self.target_phase_var.get(),
            "harvest_mode": self.harvest_mode_var.get(),
            "database_provider": self.database_provider_var.get(),
            "database_endpoint": self.database_endpoint_var.get(),
            "vesta_exe": self.vesta_exe_var.get(),
            "rietan_exe": self.rietan_exe_var.get(),
            "vesta_reference_dir": self.vesta_reference_dir_var.get(),
            "scientific_safe_mode": self.scientific_safe_mode_var.get(),
        }

    def _save_app_state(self) -> None:
        save_app_state(default_app_state_path(), self._collect_app_state())

    def _on_close(self) -> None:
        if self._phase_stripping_is_open():
            self.phase_stripping_window._request_close()
            if self._phase_stripping_is_open():
                return
        if self._refinement_is_open():
            self.refinement_window._request_close()
        try:
            self._save_app_state()
        finally:
            self.destroy()

    def _file_step(self, parent: ttk.Frame, number: str, title: str, hint: str, var: tk.StringVar, command) -> None:
        box = ttk.Frame(parent, style="Subtle.TFrame", padding=10)
        box.pack(fill="x", pady=5)
        label = f"{number}. {title}" if number else title
        ttk.Label(box, text=label, style="Subtle.TLabel").pack(anchor="w")
        ttk.Label(box, text=hint, style="Subtle.TLabel").pack(anchor="w", pady=(2, 5))
        row = ttk.Frame(box, style="Subtle.TFrame")
        row.pack(fill="x")
        entry = ttk.Entry(row, textvariable=var)
        entry.pack(side="left", fill="x", expand=True)
        ttk.Button(row, text="选择", command=command).pack(side="left", padx=(8, 0))

    def _entry_step(self, parent: ttk.Frame, number: str, title: str, label: str, var: tk.StringVar) -> None:
        box = ttk.Frame(parent, style="Subtle.TFrame", padding=10)
        box.pack(fill="x", pady=5)
        heading = f"{number}. {title}" if number else title
        ttk.Label(box, text=heading, style="Subtle.TLabel").pack(anchor="w")
        ttk.Label(box, text=label, style="Subtle.TLabel").pack(anchor="w", pady=(2, 5))
        ttk.Entry(box, textvariable=var).pack(fill="x")

    def _element_step(self, parent: ttk.Frame) -> None:
        box = ttk.Frame(parent, style="Subtle.TFrame", padding=10)
        box.pack(fill="x", pady=5)
        ttk.Label(box, text="3. 确认元素范围", style="Subtle.TLabel").pack(anchor="w")
        ttk.Label(box, text="主元素（可手动输入，也可从周期表多选）", style="Subtle.TLabel").pack(anchor="w", pady=(2, 5))
        row = ttk.Frame(box, style="Subtle.TFrame")
        row.pack(fill="x")
        ttk.Entry(row, textvariable=self.elements_var).pack(side="left", fill="x", expand=True)
        self.periodic_table_button = ttk.Button(row, text="元素周期表", command=self._open_periodic_table)
        self.periodic_table_button.pack(side="left", padx=(8, 0))

    def _open_periodic_table(self) -> None:
        selected = show_periodic_table(self, self.elements_var.get())
        if selected is None:
            return
        self.elements_var.set(" ".join(selected))
        self._refresh_target_phase_options()
        self.status_var.set("已更新主元素范围；请重新明确选择目标相。")

    def _target_phase_step(self, parent: ttk.Frame) -> None:
        box = ttk.Frame(parent, style="Subtle.TFrame", padding=10)
        box.pack(fill="x", pady=5)
        ttk.Label(box, text="6. 指定目标相", style="Subtle.TLabel").pack(anchor="w")
        ttk.Label(box, text="目标相由你指定；剩余候选只作为非目标杂质解释未覆盖峰。", style="Subtle.TLabel").pack(anchor="w", pady=(2, 5))
        self.target_phase_combo = ttk.Combobox(
            box,
            textvariable=self.target_phase_var,
            values=(),
            state="readonly",
            height=10,
        )
        self.target_phase_combo.pack(fill="x")
        row = ttk.Frame(box, style="Subtle.TFrame")
        row.pack(fill="x", pady=(8, 0))
        ttk.Button(row, text="刷新目标相", command=self._refresh_target_phase_options).pack(side="left", fill="x", expand=True)
        self.import_target_button = ttk.Button(row, text="导入目标相 CIF", command=self._import_target_phase_cif)
        self.import_target_button.pack(side="left", fill="x", expand=True, padx=(8, 0))

    def _harvest_step(self, parent: ttk.Frame) -> None:
        box = ttk.Frame(parent, style="Subtle.TFrame", padding=10)
        box.pack(fill="x", pady=5)
        ttk.Label(box, text="7. 连接结构数据库并采集", style="Subtle.TLabel").pack(anchor="w")
        ttk.Label(box, text="当前启用 Materials Project；COD/AFLOW/OQMD/NOMAD 入口先保留，后续逐步接入。", style="Subtle.TLabel").pack(anchor="w", pady=(2, 5))
        ttk.Combobox(
            box,
            textvariable=self.database_provider_var,
            values=("Materials Project", "COD", "AFLOW", "OQMD", "NOMAD/OPTIMADE"),
            state="readonly",
        ).pack(fill="x", pady=(0, 6))
        ttk.Entry(box, textvariable=self.mp_api_key_var, show="*").pack(fill="x")
        row = ttk.Frame(box, style="Subtle.TFrame")
        row.pack(fill="x", pady=(8, 0))
        ttk.Combobox(row, textvariable=self.harvest_mode_var, values=("conservative", "normal", "broad"), width=14, state="readonly").pack(side="left")
        ttk.Button(row, text="保存 API", command=self._save_database_api).pack(side="left", padx=(8, 0))
        ttk.Button(row, text="连接测试", command=self._test_database_api).pack(side="left", padx=(8, 0))
        ttk.Button(row, text="预览采集", command=self._preview_mp_harvest).pack(side="left", padx=(8, 0))
        ttk.Button(row, text="按当前元素采集", command=self._harvest_from_selected_database).pack(side="left", padx=(8, 0))

    def _build_result_panel(self, parent: ttk.Frame) -> None:
        panel = ttk.Frame(parent, style="Panel.TFrame", padding=18)
        panel.grid(row=0, column=2, sticky="nsew")
        panel.rowconfigure(3, weight=1)
        panel.columnconfigure(0, weight=1)

        ttk.Label(panel, text="分析状态与结果", style="Section.TLabel").grid(row=0, column=0, sticky="w")
        ttk.Label(
            panel,
            text="这里会显示 Ophi 对候选相的判断。注意：这是候选筛选，不是确定性定相。",
            style="PanelMuted.TLabel",
        ).grid(row=1, column=0, sticky="w", pady=(4, 14))

        self.status_var = tk.StringVar(value="等待输入：请选择实验谱和候选来源。")
        ttk.Label(panel, textvariable=self.status_var, style="Step.TLabel").grid(row=2, column=0, sticky="w", pady=(0, 10))

        self.tabs = ttk.Notebook(panel)
        self.tabs.grid(row=3, column=0, sticky="nsew")
        self.summary = self._make_text_tab(self.tabs, "结果")
        self._build_xrd_plot_tab(self.tabs)
        self._build_library_tab(self.tabs)
        self._build_candidate_phase_tab(self.tabs)
        self._build_peak_inspector_tab(self.tabs)
        self.used_structures_text = self._make_text_tab(self.tabs, "已用结构")
        self.report_paths = self._make_text_tab(self.tabs, "导出")
        self.log_text = self._make_text_tab(self.tabs, "日志")
        self._write_summary(
            "Ophi 已准备好。\n\n"
            "建议流程：\n"
            "1. 选择一个实验 XRD 文件。\n"
            "2. 选择候选结构或模拟峰文件夹。\n"
            "3. 确认主元素和可能杂质元素。\n"
            "4. 可先点击“构建/更新缓存”。\n"
            "5. 点击“开始候选筛选”。\n\n"
            "输出会包括候选相排序、相对贡献 proxy、XRD 图片、缺失强峰、未解释峰和报告文件路径。"
        )
        self._write_text_widget(self.candidate_phase_text, "完成一次结构库分析后，这里会显示候选相排名、置信标签和评分分解。")
        self._write_text_widget(self.used_structures_text, "完成一次结构库分析后，这里会显示实际使用的结构和跳过原因摘要。")
        self._write_report_paths("报告还未生成。完成一次筛选后，这里会列出 JSON、Markdown、CSV 和图片路径。")
        self._write_log("等待任务。")

    def _build_xrd_plot_tab(self, tabs: ttk.Notebook) -> None:
        frame = ttk.Frame(tabs, style="Panel.TFrame", padding=(0, 8, 0, 0))
        tabs.add(frame, text="XRD 图")
        frame.rowconfigure(1, weight=1)
        frame.columnconfigure(0, weight=1)
        actions = ttk.Frame(frame, style="Panel.TFrame")
        actions.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        ttk.Label(actions, text="实验谱 + 目标模拟峰 + 三个可能杂质", style="Panel.TLabel").pack(side="left")
        self.save_plot_button = ttk.Button(actions, text="保存图片", command=self._save_current_plot)
        self.save_plot_button.pack(side="right")
        self.save_plot_button.state(["disabled"])
        self.plot_canvas = tk.Label(
            frame,
            bg=COLORS["panel"],
            fg=COLORS["muted"],
            text="完成一次分析后，这里会显示 XRD 对比图。",
            anchor="center",
            padx=12,
            pady=12,
            font=FONTS["body_zh"],
        )
        self.plot_canvas.grid(row=1, column=0, sticky="nsew")

    def _build_candidate_phase_tab(self, tabs: ttk.Notebook) -> None:
        frame = ttk.Frame(tabs, style="Panel.TFrame", padding=(0, 8, 0, 0))
        tabs.add(frame, text="候选相")
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        columns = ("phase", "space_group", "score", "simulation", "source")
        self.candidate_phase_tree = ttk.Treeview(frame, columns=columns, show="tree headings", height=12)
        headings = {
            "#0": "记录",
            "phase": "PhaseCandidate",
            "space_group": "Space Group",
            "score": "Score",
            "simulation": "Simulation",
            "source": "Source",
        }
        widths = {"#0": 85, "phase": 145, "space_group": 105, "score": 70, "simulation": 120, "source": 150}
        self.candidate_phase_tree.heading("#0", text=headings["#0"])
        self.candidate_phase_tree.column("#0", width=widths["#0"], stretch=False)
        for column in columns:
            self.candidate_phase_tree.heading(column, text=headings[column])
            self.candidate_phase_tree.column(column, width=widths[column], anchor="w", stretch=True)
        self.candidate_phase_tree.grid(row=0, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(frame, orient="vertical", command=self.candidate_phase_tree.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        self.candidate_phase_tree.configure(yscrollcommand=scroll.set)
        self.candidate_phase_details: dict[str, str] = {}
        self.candidate_phase_tree.bind("<<TreeviewSelect>>", self._show_candidate_phase_detail)
        self.candidate_phase_text = tk.Text(
            frame,
            bg=COLORS["panel"],
            fg=COLORS["text"],
            insertbackground=COLORS["text"],
            relief="flat",
            wrap="word",
            height=7,
            padx=10,
            pady=10,
            font=FONTS["small"],
        )
        self.candidate_phase_text.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(8, 0))

    def _show_candidate_phase_detail(self, _event=None) -> None:
        selected = self.candidate_phase_tree.selection()
        if not selected:
            return
        self._write_text_widget(
            self.candidate_phase_text,
            self.candidate_phase_details.get(str(selected[0]), "没有可显示的结构溯源。"),
        )

    def _select_section(self, key: str) -> None:
        tab_by_key = {
            "projects": "结果",
            "samples": "结果",
            "library": "结构库",
            "xrd": "结果",
            "phase": "查峰",
            "settings": "导出",
        }
        target = tab_by_key.get(key, "结果摘要")
        for tab_id in self.tabs.tabs():
            if self.tabs.tab(tab_id, "text") == target:
                self.tabs.select(tab_id)
                break
        labels = {section["key"]: section["label"] for section in workbench_sections()}
        self.status_var.set(f"当前工作区：{labels.get(key, key)}")

    def _build_library_tab(self, tabs: ttk.Notebook) -> None:
        frame = ttk.Frame(tabs, style="Panel.TFrame", padding=(0, 8, 0, 0))
        tabs.add(frame, text="结构库")
        frame.rowconfigure(1, weight=1)
        frame.columnconfigure(0, weight=1)
        actions = ttk.Frame(frame, style="Panel.TFrame")
        actions.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        ttk.Button(actions, text="刷新结构库", command=self._refresh_library_table).pack(side="left")
        ttk.Button(actions, text="启用/禁用选中", command=self._toggle_selected_library_entry).pack(side="left", padx=(8, 0))
        ttk.Button(actions, text="打开 CIF 位置", command=self._open_selected_cif_location).pack(side="left", padx=(8, 0))
        ttk.Button(actions, text="查看 CIF", command=self._view_selected_cif_text).pack(side="left", padx=(8, 0))
        ttk.Button(actions, text="查看模拟峰", command=self._view_selected_peak_table).pack(side="left", padx=(8, 0))
        ttk.Button(actions, text="重算 XRD", command=self._recompute_selected_xrd).pack(side="left", padx=(8, 0))
        columns = ("formula", "space_group", "system", "source", "source_id", "enabled", "cache", "backend", "ready")
        self.library_tree = ttk.Treeview(frame, columns=columns, show="headings", height=12)
        headings = {
            "formula": "Formula",
            "space_group": "Space Group",
            "system": "System",
            "source": "Source",
            "source_id": "Source ID",
            "enabled": "Enabled",
            "cache": "XRD Cache",
            "backend": "Backend",
            "ready": "Candidate",
        }
        widths = {
            "formula": 110,
            "space_group": 105,
            "system": 100,
            "source": 105,
            "source_id": 110,
            "enabled": 70,
            "cache": 82,
            "backend": 150,
            "ready": 82,
        }
        for column in columns:
            self.library_tree.heading(column, text=headings[column])
            self.library_tree.column(column, width=widths[column], anchor="w", stretch=True)
        self.library_tree.grid(row=1, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(frame, orient="vertical", command=self.library_tree.yview)
        scroll.grid(row=1, column=1, sticky="ns")
        self.library_tree.configure(yscrollcommand=scroll.set)
        self.library_tree.bind("<<TreeviewSelect>>", self._show_selected_library_entry)
        self.library_detail = tk.Text(
            frame,
            bg=COLORS["panel"],
            fg=COLORS["text"],
            insertbackground=COLORS["text"],
            relief="flat",
            wrap="word",
            height=7,
            padx=10,
            pady=10,
            font=FONTS["small"],
        )
        self.library_detail.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(8, 0))

    def _build_peak_inspector_tab(self, tabs: ttk.Notebook) -> None:
        frame = ttk.Frame(tabs, style="Panel.TFrame", padding=(0, 8, 0, 0))
        tabs.add(frame, text="查峰")
        frame.rowconfigure(1, weight=1)
        frame.columnconfigure(0, weight=1)
        controls = ttk.Frame(frame, style="Panel.TFrame")
        controls.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        ttk.Label(controls, text="Experimental 2theta", style="Panel.TLabel").pack(side="left")
        self.inspect_theta_var = tk.StringVar()
        ttk.Entry(controls, textvariable=self.inspect_theta_var, width=12).pack(side="left", padx=(8, 8))
        ttk.Button(controls, text="检查峰", command=self._inspect_peak).pack(side="left")
        self.inspector_text = tk.Text(
            frame,
            bg=COLORS["panel"],
            fg=COLORS["text"],
            insertbackground=COLORS["text"],
            relief="flat",
            wrap="word",
            height=22,
            padx=14,
            pady=14,
            font=FONTS["body_zh"],
        )
        self.inspector_text.grid(row=1, column=0, sticky="nsew")
        self._write_inspector("选择或输入一个实验峰 2theta，Ophi 会从本地结构库缓存中列出附近理论峰。")

    def _make_text_tab(self, tabs: ttk.Notebook, title: str) -> tk.Text:
        frame = ttk.Frame(tabs, style="Panel.TFrame", padding=(0, 8, 0, 0))
        tabs.add(frame, text=title)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        text = tk.Text(
            frame,
            bg=COLORS["panel"],
            fg=COLORS["text"],
            insertbackground=COLORS["text"],
            relief="flat",
            wrap="word",
            height=28,
            padx=14,
            pady=14,
            font=FONTS["body_zh"],
        )
        text.grid(row=0, column=0, sticky="nsew")
        return text

    def _pick_xrd(self) -> None:
        path = filedialog.askopenfilename(
            initialdir=str(first_existing_directory(self.xrd_var.get(), fallback=default_xrd_dir())),
            filetypes=[("XRD files", "*.asc *.ras *.raw *.xy *.txt *.csv *.dat"), ("All files", "*.*")],
        )
        if path:
            self.xrd_var.set(path)
            self.status_var.set("已选择实验谱。下一步选择候选结构/模拟峰文件夹。")

    def _pick_dir(self) -> None:
        path = filedialog.askdirectory(
            initialdir=str(first_existing_directory(self.dir_var.get(), fallback=default_candidate_dir()))
        )
        if path:
            self.dir_var.set(path)
            self.status_var.set("已选择候选来源。确认元素范围后即可开始筛选。")

    def _pick_out(self) -> None:
        path = filedialog.askdirectory(
            initialdir=str(first_existing_directory(self.out_var.get(), fallback=project_root() / "results"))
        )
        if path:
            self.out_var.set(path)

    def _pick_cache(self) -> None:
        path = filedialog.asksaveasfilename(
            initialdir=str(first_existing_directory(self.cache_var.get(), fallback=default_cache_path().parent)),
            initialfile=default_cache_path().name,
            defaultextension=".sqlite",
            filetypes=[("SQLite cache", "*.sqlite *.db"), ("All files", "*.*")],
        )
        if path:
            self.cache_var.set(path)

    def _pick_library(self) -> None:
        path = filedialog.asksaveasfilename(
            initialdir=str(first_existing_directory(self.library_var.get(), fallback=default_library_path().parent)),
            initialfile=default_library_path().name,
            defaultextension=".sqlite",
            filetypes=[("SQLite library", "*.sqlite *.db"), ("All files", "*.*")],
        )
        if path:
            self.library_var.set(path)
            self._refresh_target_phase_options()

    def _refresh_target_phase_options(self) -> None:
        try:
            rows = library_target_phase_options(self.library_var.get(), self.elements_var.get())
        except Exception as exc:
            self.target_phase_options = {}
            if hasattr(self, "target_phase_combo"):
                self.target_phase_combo.configure(values=())
            self.target_phase_var.set("")
            self._write_log(f"目标相列表无法读取：{exc}")
            return
        selected_label, _selected_id, options = resolve_target_phase_selection(rows, self.target_phase_var.get())
        self.target_phase_options = options
        labels = list(options)
        if hasattr(self, "target_phase_combo"):
            self.target_phase_combo.configure(values=labels)
        self.target_phase_var.set(selected_label)
        if labels:
            if hasattr(self, "log_text"):
                self._write_log(f"目标相候选已刷新：{len(labels)} 个。当前目标相：{self.target_phase_var.get()}")

    def _import_target_phase_cif(self) -> None:
        if not self.library_var.get().strip():
            self.status_var.set("还不能导入目标相：请先选择结构库数据库。")
            messagebox.showwarning("还缺少必要输入", "请先选择本地结构库数据库。")
            return
        path = filedialog.askopenfilename(
            initialdir=str(first_existing_directory(self.dir_var.get(), fallback=default_candidate_dir())),
            title="选择目标相 CIF",
            filetypes=[("CIF files", "*.cif"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            imported = import_target_cif_to_library(path, self.library_var.get())
            if not self.elements_var.get().strip() and imported.get("elements"):
                self.elements_var.set(" ".join(imported["elements"]))
            build_library_cache_from_folder(self.library_var.get())
            self._refresh_library_table()
            self._refresh_target_phase_options()
            label = str(imported.get("label") or "")
            if label in self.target_phase_options:
                self.target_phase_var.set(label)
            self.status_var.set("目标相已导入并选中。")
            self._write_log(f"目标相 CIF 已导入：\n{path}\n当前目标相：{self.target_phase_var.get()}")
        except Exception as exc:
            guidance = self._friendly_error(exc)
            self.status_var.set("目标相导入失败。")
            self._write_log(guidance)
            messagebox.showerror("目标相导入失败", guidance)

    def _run(self) -> None:
        if self.running:
            return
        if self._phase_stripping_is_open():
            messagebox.showwarning("残差谱会话仍在打开", "请先关闭当前手动物相剥离窗口，再开始新的主分析。")
            return
        if self._refinement_is_open():
            messagebox.showwarning("精修会话仍在打开", "请先关闭当前 RIETAN 精修窗口，再开始新的主分析。")
            return
        errors = validate_analysis_inputs(
            self.xrd_var.get(),
            self.dir_var.get(),
            self.elements_var.get(),
            self.out_var.get(),
        )
        if errors:
            self.status_var.set("还不能开始：请先补齐必要输入。")
            messagebox.showwarning("还缺少必要输入", "\n".join(errors))
            return
        try:
            selected_elements = parse_element_symbols(self.elements_var.get())
            selected_extras = parse_element_symbols(self.extra_var.get())
        except ValueError as exc:
            self.status_var.set("还不能开始：元素范围含有无效符号。")
            messagebox.showwarning("元素范围无效", str(exc))
            return
        preflight = self._vesta_preflight_for_candidate_dir()
        self.running = True
        self.run_button.state(["disabled"])
        self.cache_button.state(["disabled"])
        self.library_button.state(["disabled"])
        self.library_run_button.state(["disabled"])
        self.status_var.set("正在读取实验谱并扫描候选来源...")
        summary_lines = ["开始候选筛选。", "", "步骤：", "- 读取实验谱", "- 提取实验峰", "- 扫描候选结构/模拟峰", "- 匹配候选相", "- 生成报告", ""]
        summary_lines.extend(vesta_preflight_lines(preflight))
        self._write_summary("\n".join(summary_lines))
        self._write_log("开始分析。缓存路径：\n" + self.cache_var.get() + "\n\n" + "\n".join(vesta_preflight_lines(preflight)))
        self._analysis_run_id += 1
        run_id = self._analysis_run_id
        snapshot = {
            "xrd_file": self.xrd_var.get(),
            "candidate_dir": self.dir_var.get(),
            "elements": selected_elements,
            "extra_elements": selected_extras,
            "cache_path": self.cache_var.get().strip(),
        }
        self.save_analysis_button.state(["disabled"])
        self.phase_stripping_button.state(["disabled"])
        self.refinement_button.state(["disabled"])
        thread = threading.Thread(target=self._worker, args=(run_id, snapshot), daemon=True)
        thread.start()

    def _build_cache(self) -> None:
        errors: list[str] = []
        candidate_dir = self.dir_var.get().strip()
        if not candidate_dir:
            errors.append("请选择候选结构/模拟峰文件夹。")
        elif not Path(candidate_dir).is_dir():
            errors.append("候选文件夹不存在，请重新选择。")
        elif not any(path.is_file() and path.suffix.lower() in CANDIDATE_EXTENSIONS for path in Path(candidate_dir).rglob("*")):
            errors.append("候选文件夹里没有找到 CIF 或峰表文件。")
        if not self.elements_var.get().strip():
            errors.append("请填写主元素，例如 Zr Fe Ge。")
        if not self.cache_var.get().strip():
            errors.append("请选择本地缓存文件。")
        if errors:
            self.status_var.set("还不能构建缓存：请先确认候选来源和元素范围。")
            messagebox.showwarning("还缺少必要输入", "\n".join(errors))
            return
        self.running = True
        self.run_button.state(["disabled"])
        self.cache_button.state(["disabled"])
        self.library_button.state(["disabled"])
        self.library_run_button.state(["disabled"])
        self.status_var.set("正在构建候选缓存...")
        self._write_summary("开始构建/更新候选缓存。\n\nOphi 会扫描候选文件夹，把可读取的 CIF 或峰表保存到本地缓存，后续筛选会更快。")
        self._write_log("开始构建缓存。缓存路径：\n" + self.cache_var.get())
        thread = threading.Thread(target=self._cache_worker, daemon=True)
        thread.start()

    def _import_library(self) -> None:
        candidate_dir = self.dir_var.get().strip()
        if not candidate_dir or not Path(candidate_dir).is_dir():
            self.status_var.set("还不能导入结构库：请先选择包含 CIF 的文件夹。")
            messagebox.showwarning("还缺少必要输入", "请选择候选结构/模拟峰文件夹，里面应包含本地 CIF。")
            return
        if not self.library_var.get().strip():
            self.status_var.set("还不能导入结构库：请先选择结构库数据库。")
            messagebox.showwarning("还缺少必要输入", "请选择本地结构库数据库。")
            return
        self.running = True
        self.run_button.state(["disabled"])
        self.cache_button.state(["disabled"])
        self.library_button.state(["disabled"])
        self.library_run_button.state(["disabled"])
        self.status_var.set("正在导入本地 CIF 到结构库...")
        self._write_summary("开始导入结构库。\n\nOphi 会复制本地 CIF 到结构库目录，并保存来源、结构哈希和访问备注。")
        self._write_log("开始导入结构库。数据库：\n" + self.library_var.get())
        thread = threading.Thread(target=self._library_worker, daemon=True)
        thread.start()

    def _run_library(self) -> None:
        if self.running:
            return
        if self._phase_stripping_is_open():
            messagebox.showwarning("残差谱会话仍在打开", "请先关闭当前手动物相剥离窗口，再开始新的结构库分析。")
            return
        if self._refinement_is_open():
            messagebox.showwarning("精修会话仍在打开", "请先关闭当前 RIETAN 精修窗口，再开始新的结构库分析。")
            return
        errors: list[str] = []
        if not self.xrd_var.get().strip() or not Path(self.xrd_var.get()).is_file():
            errors.append("请选择实验 XRD 文件。")
        if not self.library_var.get().strip():
            errors.append("请选择本地结构库数据库。")
        if not self.elements_var.get().strip():
            errors.append("请填写主元素，例如 Zr Fe Ge。")
        try:
            selected_elements = parse_element_symbols(self.elements_var.get())
        except ValueError as exc:
            selected_elements = ()
            errors.append(f"主元素包含无效符号：{exc}")
        inferred = infer_elements_from_xrd_path(self.xrd_var.get()) if self.xrd_var.get().strip() else ()
        if inferred and selected_elements:
            missing, _extra = element_scope_mismatch(selected_elements, inferred)
            if missing:
                errors.append("实验谱名称推断出的元素没有全部包含在主元素范围中：" + " ".join(missing))
        self._refresh_target_phase_options()
        if not self.target_phase_var.get().strip():
            errors.append("请选择目标相。可以先点“刷新目标相”。")
        if errors:
            self.status_var.set("还不能进行结构库分析：请先补齐必要输入。")
            messagebox.showwarning("还缺少必要输入", "\n".join(errors))
            return
        preview_lines = [
            "开始结构库分析。",
            "",
            "本轮只分析主元素相关结构，不扫描额外杂质元素。",
            f"Scientific Safe Mode：{'开启（直接重算 CIF）' if self.scientific_safe_mode_var.get() else '关闭（允许使用指纹验证缓存）'}",
            f"本次强制重算：{'是' if self.force_recompute_var.get() else '否'}",
        ]
        log_lines = [
            "开始结构库分析。仅使用主元素范围：",
            self.elements_var.get(),
            "结构库：",
            self.library_var.get(),
            "目标相：",
            self.target_phase_var.get(),
        ]
        try:
            library = StructureLibrary(self.library_var.get())
            elements = self.elements_var.get().split()
            scoped_ids = scoped_structure_ids(library, elements, enabled_only=True, scope="subsystems")
            preflight = self._vesta_preflight_for_target_phase()
            preview_lines.append(f"预计分析结构数：{len(scoped_ids)}")
            preview_lines.append("范围规则：元素集合必须完全落在主元素内，例如 Zr Fe Ge 会分析 Zr、Fe-Ge、Fe-Ge-Zr 等子体系。")
            preview_lines.append(f"目标相：{self.target_phase_var.get()}")
            preview_lines.extend(["", *vesta_preflight_lines(preflight)])
            log_lines.append(f"预计分析结构数：{len(scoped_ids)}")
            log_lines.extend(["", *vesta_preflight_lines(preflight)])
        except Exception as exc:
            preview_lines.append("预计分析结构数：暂时无法读取结构库，开始后会在日志里报告具体错误。")
            log_lines.append(f"结构库预检查失败：{exc}")
        self.running = True
        self.run_button.state(["disabled"])
        self.cache_button.state(["disabled"])
        self.library_button.state(["disabled"])
        self.library_run_button.state(["disabled"])
        self.status_var.set("正在使用本地结构库分析实验谱...")
        self._write_summary("\n".join(preview_lines))
        self._write_log("\n".join(log_lines))
        self._analysis_run_id += 1
        run_id = self._analysis_run_id
        snapshot = {
            "xrd_file": self.xrd_var.get(),
            "library_db": self.library_var.get(),
            "elements": tuple(selected_elements),
            "target_candidate_id": self.target_phase_options.get(self.target_phase_var.get(), ""),
            "scientific_safe_mode": self.scientific_safe_mode_var.get(),
            "force_recompute": self.force_recompute_var.get(),
        }
        self.save_analysis_button.state(["disabled"])
        self.phase_stripping_button.state(["disabled"])
        self.refinement_button.state(["disabled"])
        thread = threading.Thread(target=self._library_analysis_worker, args=(run_id, snapshot), daemon=True)
        thread.start()

    def _save_mp_key(self) -> None:
        key = self.mp_api_key_var.get().strip()
        if not key:
            self._write_log("MP API key 为空，未保存。")
            return
        save_mp_api_key_to_env(default_env_path(), key)
        self._write_log(f"MP API key 已保存到本地 .env：{default_env_path()}\n.env 已加入 .gitignore。")

    def _open_database_dialog(self) -> None:
        self._load_saved_database_config()
        dialog = tk.Toplevel(self)
        dialog.title("连接 API / 数据库")
        dialog.configure(bg=COLORS["background"])
        dialog.transient(self)
        dialog.grab_set()
        dialog.geometry("680x350")
        box = ttk.Frame(dialog, style="Panel.TFrame", padding=18)
        box.pack(fill="both", expand=True, padx=18, pady=18)
        ttk.Label(box, text="结构数据库连接", style="Section.TLabel").pack(anchor="w")
        ttk.Label(box, text="当前先启用 Materials Project；其他常用数据库保留入口。", style="PanelMuted.TLabel").pack(anchor="w", pady=(4, 14))
        ttk.Label(box, text="Database", style="Panel.TLabel").pack(anchor="w")
        ttk.Combobox(
            box,
            textvariable=self.database_provider_var,
            values=("Materials Project", "COD", "AFLOW", "OQMD", "NOMAD/OPTIMADE"),
            state="readonly",
        ).pack(fill="x", pady=(2, 10))
        ttk.Label(box, text="API Key / Token", style="Panel.TLabel").pack(anchor="w")
        ttk.Entry(box, textvariable=self.mp_api_key_var, show="*").pack(fill="x", pady=(2, 10))
        ttk.Label(box, text="Endpoint", style="Panel.TLabel").pack(anchor="w")
        ttk.Entry(box, textvariable=self.database_endpoint_var).pack(fill="x", pady=(2, 14))
        row = ttk.Frame(box, style="Panel.TFrame")
        row.pack(fill="x")
        ttk.Button(row, text="保存 API", command=self._save_database_api).pack(side="left")
        ttk.Button(row, text="连接测试", command=self._test_database_api).pack(side="left", padx=(8, 0))
        actions = ttk.Frame(box, style="Panel.TFrame")
        actions.pack(fill="x", pady=(10, 0))
        ttk.Button(
            actions,
            text="按当前元素采集",
            style="Primary.TButton",
            command=self._harvest_from_selected_database,
        ).pack(side="left")
        ttk.Button(actions, text="补齐常见氧化物", command=self._supplement_common_oxides).pack(
            side="left", padx=(8, 0)
        )

    def _open_vesta_dialog(self) -> None:
        if self._vesta_dialog is not None and self._vesta_dialog.winfo_exists():
            self._present_vesta_dialog(self._vesta_dialog)
            return
        config = load_vesta_config(default_env_path())
        if not self.vesta_exe_var.get().strip():
            self.vesta_exe_var.set(config["vesta_exe"])
        if not self.vesta_reference_dir_var.get().strip():
            self.vesta_reference_dir_var.set(config["reference_dir"])
        if not self.rietan_exe_var.get().strip():
            self.rietan_exe_var.set(config["rietan_exe"])
        dialog = tk.Toplevel(self)
        self._vesta_dialog = dialog
        dialog.title("VESTA / RIETAN 设置")
        dialog.configure(bg=COLORS["background"])
        dialog.transient(self)
        width, height = 640, 530
        self.update_idletasks()
        x = max(0, self.winfo_rootx() + (self.winfo_width() - width) // 2)
        y = max(0, self.winfo_rooty() + (self.winfo_height() - height) // 2)
        x = min(x, max(0, self.winfo_screenwidth() - width))
        y = min(y, max(0, self.winfo_screenheight() - height - 48))
        dialog.geometry(f"{width}x{height}+{x}+{y}")

        def close_dialog() -> None:
            self._vesta_dialog = None
            dialog.destroy()

        dialog.protocol("WM_DELETE_WINDOW", close_dialog)
        box = ttk.Frame(dialog, style="Panel.TFrame", padding=18)
        box.pack(fill="both", expand=True, padx=18, pady=18)
        ttk.Label(box, text="VESTA / RIETAN-FP 设置", style="Section.TLabel").pack(anchor="w")
        ttk.Label(
            box,
            text="此处只管理程序路径和参考文件。开始候选筛选或结构库分析时，Ophi 会自动调用 VESTA / RIETAN-FP 生成最终展示谱。",
            style="PanelMuted.TLabel",
            wraplength=500,
        ).pack(anchor="w", pady=(4, 14))
        self._vesta_path_row(box, "VESTA.exe", self.vesta_exe_var, self._pick_vesta_exe)
        self._vesta_path_row(box, "RIETAN.exe", self.rietan_exe_var, self._pick_rietan_exe)
        self._vesta_path_row(box, "VESTA 参考目录", self.vesta_reference_dir_var, self._pick_vesta_reference_dir)
        ttk.Label(box, text="公式级标准参考", style="Panel.TLabel").pack(anchor="w", pady=(8, 0))
        formula_row = ttk.Frame(box, style="Panel.TFrame")
        formula_row.pack(fill="x", pady=(3, 8))
        ttk.Entry(formula_row, textvariable=self.vesta_formula_var, width=18).pack(side="left")
        ttk.Entry(formula_row, textvariable=self.vesta_formula_reference_var).pack(side="left", fill="x", expand=True, padx=(8, 0))
        ttk.Button(formula_row, text="选择参考", command=self._pick_formula_vesta_reference).pack(side="left", padx=(8, 0))
        status = ttk.Label(box, text=self._vesta_status_text(), style="PanelMuted.TLabel", wraplength=500)
        status.pack(anchor="w", pady=(10, 0))

        def save_and_refresh() -> None:
            saved = save_vesta_config(
                default_env_path(),
                self.vesta_exe_var.get(),
                self.vesta_reference_dir_var.get(),
                self.rietan_exe_var.get(),
            )
            if self.vesta_exe_var.get().strip():
                os.environ["OPHI_VESTA_EXE"] = self.vesta_exe_var.get().strip()
            if self.rietan_exe_var.get().strip():
                os.environ["OPHI_RIETAN_EXE"] = self.rietan_exe_var.get().strip()
            if self.vesta_reference_dir_var.get().strip():
                os.environ["OPHI_VESTA_REFERENCE_DIR"] = self.vesta_reference_dir_var.get().strip()
            status.configure(text=self._vesta_status_text())
            self._write_log(
                "VESTA 配置已保存。\n"
                f"VESTA.exe: {self.vesta_exe_var.get() or '未设置'}\n"
                f"RIETAN.exe: {self.rietan_exe_var.get() or '未设置'}\n"
                f"参考目录: {self.vesta_reference_dir_var.get() or '未设置'}\n"
                f"VESTA: {saved['vesta_exe_exists']} | RIETAN: {saved['rietan_exe_exists']} | 参考目录: {saved['reference_dir_exists']}"
            )

        actions = ttk.Frame(box, style="Panel.TFrame")
        actions.pack(fill="x", pady=(16, 0))
        ttk.Button(actions, text="保存配置", style="Primary.TButton", command=save_and_refresh).pack(side="left")
        ttk.Button(actions, text="保存公式参考", command=self._save_formula_vesta_reference_from_dialog).pack(side="left", padx=(8, 0))
        ttk.Button(actions, text="打开 VESTA", command=self._launch_vesta_empty).pack(side="left", padx=(8, 0))
        ttk.Button(actions, text="打开目标相文件", command=self._launch_vesta_target_phase).pack(side="left", padx=(8, 0))
        ttk.Button(actions, text="关闭", command=close_dialog).pack(side="right")
        self._present_vesta_dialog(dialog)

    @staticmethod
    def _present_vesta_dialog(dialog: tk.Toplevel) -> None:
        dialog.deiconify()
        dialog.lift()
        dialog.attributes("-topmost", True)
        dialog.focus_force()

        def release_topmost() -> None:
            if dialog.winfo_exists():
                dialog.attributes("-topmost", False)

        dialog.after(800, release_topmost)

    def _vesta_path_row(self, parent: ttk.Frame, label: str, variable: tk.StringVar, picker) -> None:
        ttk.Label(parent, text=label, style="Panel.TLabel").pack(anchor="w")
        row = ttk.Frame(parent, style="Panel.TFrame")
        row.pack(fill="x", pady=(3, 8))
        ttk.Entry(row, textvariable=variable).pack(side="left", fill="x", expand=True)
        ttk.Button(row, text="选择", command=picker).pack(side="left", padx=(8, 0))

    def _pick_vesta_exe(self) -> None:
        path = filedialog.askopenfilename(
            title="选择 VESTA.exe",
            initialdir=str(first_existing_directory(self.vesta_exe_var.get(), fallback=desktop_dir())),
            filetypes=[("VESTA executable", "VESTA.exe"), ("Executable", "*.exe"), ("All files", "*.*")],
        )
        if path:
            self.vesta_exe_var.set(path)

    def _pick_rietan_exe(self) -> None:
        path = filedialog.askopenfilename(
            title="选择 RIETAN.exe",
            initialdir=str(first_existing_directory(self.rietan_exe_var.get(), fallback=desktop_dir())),
            filetypes=[("RIETAN executable", "RIETAN.exe"), ("Executable", "*.exe"), ("All files", "*.*")],
        )
        if path:
            self.rietan_exe_var.set(path)

    def _pick_vesta_reference_dir(self) -> None:
        path = filedialog.askdirectory(
            title="选择 VESTA 参考峰目录",
            initialdir=str(
                first_existing_directory(self.vesta_reference_dir_var.get(), fallback=default_xrd_dir())
            ),
        )
        if path:
            self.vesta_reference_dir_var.set(path)

    def _pick_formula_vesta_reference(self) -> None:
        path = filedialog.askopenfilename(
            title="选择该公式的标准 VESTA 参考文件",
            initialdir=str(
                first_existing_directory(self.vesta_reference_dir_var.get(), fallback=default_xrd_dir())
            ),
            filetypes=[("VESTA / peak files", "*.int *.csv *.txt *.dat *.xy"), ("All files", "*.*")],
        )
        if path:
            self.vesta_formula_reference_var.set(path)
            if not self.vesta_formula_var.get().strip():
                stem = Path(path).stem.split()[0]
                self.vesta_formula_var.set(stem)

    def _save_formula_vesta_reference_from_dialog(self) -> None:
        try:
            saved = save_formula_vesta_reference(
                default_env_path(),
                self.vesta_formula_var.get(),
                self.vesta_formula_reference_var.get(),
            )
            formula = self.vesta_formula_var.get().strip()
            status = assess_vesta_preflight(self.vesta_exe_var.get(), self.vesta_reference_dir_var.get(), [formula])
            self._write_log(
                f"公式级 VESTA 参考已保存：{formula}\n"
                f"{saved['key']}={self.vesta_formula_reference_var.get()}\n"
                + "\n".join(vesta_preflight_lines(status))
            )
            self.status_var.set(f"{formula} 的 VESTA 标准参考已保存。")
        except Exception as exc:
            messagebox.showerror("无法保存公式参考", str(exc))

    def _vesta_status_text(self) -> str:
        exe = self.vesta_exe_var.get().strip()
        rietan = self.rietan_exe_var.get().strip()
        reference_dir = self.vesta_reference_dir_var.get().strip()
        return (
            f"VESTA 程序：{'已找到' if exe and Path(exe).exists() else '未找到'}\n"
            f"RIETAN-FP：{'已找到，可直接模拟' if rietan and Path(rietan).exists() else '未找到'}\n"
            f"参考目录：{'已找到' if reference_dir and Path(reference_dir).exists() else '未找到'}"
        )

    def _launch_vesta_empty(self) -> None:
        try:
            launch_vesta(self.vesta_exe_var.get().strip())
            self._write_log("已请求打开 VESTA。")
        except Exception as exc:
            messagebox.showerror("无法打开 VESTA", str(exc))

    def _launch_vesta_target_phase(self) -> None:
        try:
            target = self._selected_target_phase_file()
            launch_vesta(self.vesta_exe_var.get().strip(), target)
            self._write_log(f"已请求用 VESTA 打开目标相文件：\n{target}")
        except Exception as exc:
            messagebox.showerror("无法打开目标相文件", str(exc))

    def _selected_target_phase_file(self) -> str:
        internal_id = self.target_phase_options.get(self.target_phase_var.get(), "")
        if not internal_id:
            raise RuntimeError("请先刷新并选择目标相。")
        library = StructureLibrary(self.library_var.get())
        entry = library.get_structure(internal_id)
        path = Path(library.path.parent) / entry.cached_file_path
        if not path.exists():
            raise FileNotFoundError(f"目标相文件不存在：{path}")
        return str(path)

    def _vesta_preflight_for_target_phase(self) -> dict[str, object]:
        formula = self.target_phase_var.get().split("|", 1)[0].strip()
        return assess_vesta_preflight(
            self.vesta_exe_var.get(),
            self.vesta_reference_dir_var.get(),
            [formula] if formula else [],
        )

    def _vesta_preflight_for_candidate_dir(self) -> dict[str, object]:
        formulas: list[str] = []
        try:
            allowed = set(self.elements_var.get().split()) | set(self.extra_var.get().split())
            provider = LocalCandidateProvider([self.dir_var.get()], allowed_elements=allowed or None)
            for candidate in provider.iter_candidates()[:40]:
                if candidate.formula_pretty not in formulas:
                    formulas.append(candidate.formula_pretty)
        except Exception as exc:
            self._write_log(f"VESTA 预检无法扫描候选目录：{exc}")
        return assess_vesta_preflight(
            self.vesta_exe_var.get(),
            self.vesta_reference_dir_var.get(),
            formulas[:20],
        )

    def _save_database_api(self) -> None:
        provider_id = normalize_database_provider(self.database_provider_var.get())
        if provider_id != "materials_project":
            self._write_log(f"{self.database_provider_var.get()} 入口已保留，但当前版本还没有接入采集 API。")
            self.status_var.set("该数据库入口已保留，当前先使用 Materials Project。")
            return
        saved = save_database_api_config(
            default_env_path(),
            self.database_provider_var.get(),
            self.mp_api_key_var.get(),
            endpoint=self.database_endpoint_var.get(),
        )
        if saved["configured"]:
            self.status_var.set("Materials Project API 已保存。")
            self._write_log(f"数据库 API 已保存到本地 .env：{default_env_path()}\n采集时会按当前元素自动去重导入。")
        else:
            self.status_var.set("API key 为空，未完成连接配置。")
            self._write_log("API key 为空。请填入 Materials Project API key 后再保存。")

    def _test_database_api(self) -> None:
        provider_id = normalize_database_provider(self.database_provider_var.get())
        if provider_id != "materials_project":
            self.status_var.set("该数据库暂未接入在线测试。")
            self._write_log(f"{self.database_provider_var.get()} 暂未接入；当前可用数据库是 Materials Project。")
            return
        self._test_mp()

    def _harvest_from_selected_database(self) -> None:
        provider_id = normalize_database_provider(self.database_provider_var.get())
        if provider_id != "materials_project":
            self.status_var.set("该数据库暂未接入自动采集。")
            self._write_log(f"{self.database_provider_var.get()} 入口已保留；当前版本先支持 Materials Project 自动采集。")
            return
        if self.mp_api_key_var.get().strip():
            save_database_api_config(default_env_path(), self.database_provider_var.get(), self.mp_api_key_var.get(), endpoint=self.database_endpoint_var.get())
        self._harvest_mp()

    def _test_mp(self) -> None:
        provider = MaterialsProjectProvider(api_key=self.mp_api_key_var.get().strip() or None, env_path=default_env_path())
        try:
            result = provider.test_connection()
            self._write_log(result["message"])
            self.status_var.set("Materials Project 连接测试成功。")
        except ProviderError as exc:
            self.status_var.set("Materials Project 连接测试失败。")
            self._write_log(str(exc))

    def _preview_mp_harvest(self) -> None:
        preview = materials_project_harvest_from_app(
            library_db=self.library_var.get(),
            elements=self.elements_var.get(),
            impurities=self.extra_var.get(),
            mode=self.harvest_mode_var.get(),
            api_key=self.mp_api_key_var.get().strip(),
            env_path=default_env_path(),
            dry_run=True,
        )
        systems = preview["chemical_systems"]
        lines = [
            f"MP 采集预览：{preview['count']} 个 chemical systems",
            f"模式：{preview['mode']}",
            "",
            ", ".join(systems[:80]),
        ]
        self._write_log("\n".join(lines))
        self.status_var.set("已生成 Materials Project 采集预览。")

    def _harvest_mp(self) -> None:
        self.running = True
        self.run_button.state(["disabled"])
        self.cache_button.state(["disabled"])
        self.library_button.state(["disabled"])
        self.library_run_button.state(["disabled"])
        self.status_var.set("正在从 Materials Project 采集结构...")
        self._write_log("开始 MP 采集。若缺少 mp-api、网络或 API key，Ophi 会记录清楚的失败原因。")
        thread = threading.Thread(target=self._mp_harvest_worker, daemon=True)
        thread.start()

    def _supplement_common_oxides(self) -> None:
        if normalize_database_provider(self.database_provider_var.get()) != "materials_project":
            self.status_var.set("常见氧化物补库当前使用 Materials Project。")
            return
        elements = set(parse_element_symbols(self.elements_var.get()))
        if not elements:
            messagebox.showwarning("没有主元素", "请先在主界面选择主元素，再补齐对应氧化物。")
            return
        if not self.library_var.get().strip():
            messagebox.showwarning("没有结构库", "请先选择本地结构库数据库。")
            return
        if self.mp_api_key_var.get().strip():
            save_database_api_config(
                default_env_path(),
                self.database_provider_var.get(),
                self.mp_api_key_var.get(),
                endpoint=self.database_endpoint_var.get(),
            )
        snapshot = {
            "library_path": self.library_var.get().strip(),
            "elements": elements,
            "api_key": self.mp_api_key_var.get().strip(),
        }
        self.running = True
        for button in (self.run_button, self.cache_button, self.library_button, self.library_run_button):
            button.state(["disabled"])
        self.status_var.set("正在精确查询并补齐常见氧化物结构...")
        self._write_log("开始补齐常见氧化物：只查询当前主元素的允许公式和常见空间群。")
        threading.Thread(target=self._oxide_supplement_worker, args=(snapshot,), daemon=True).start()

    def _refresh_library_table(self) -> None:
        for item in self.library_tree.get_children():
            self.library_tree.delete(item)
        try:
            rows = library_manager_rows(self.library_var.get())
        except Exception as exc:
            self._write_library_detail(f"结构库无法读取：{exc}")
            return
        for row in rows:
            self.library_tree.insert(
                "",
                "end",
                iid=row["internal_id"],
                values=(
                    row["formula"],
                    f"{row['space_group_symbol'] or '?'} ({row['space_group_number'] or '?'})",
                    row["chemical_system"],
                    row["source"],
                    row["source_id"],
                    "yes" if row["enabled"] else "no",
                    row["xrd_cache"],
                    row["backend_version"] or "n/a",
                    "ready" if row["candidate_ready"] else "skipped",
                ),
            )
        self._write_library_detail(f"结构库条目：{len(rows)}\n数据库：{self.library_var.get()}")

    def _selected_library_id(self) -> str | None:
        selected = self.library_tree.selection()
        return str(selected[0]) if selected else None

    def _toggle_selected_library_entry(self) -> None:
        internal_id = self._selected_library_id()
        if not internal_id:
            self._write_library_detail("请先在结构库表格里选择一个条目。")
            return
        rows = {row["internal_id"]: row for row in library_manager_rows(self.library_var.get())}
        current = bool(rows[internal_id]["enabled"])
        set_library_entry_enabled(self.library_var.get(), internal_id, not current)
        self._refresh_library_table()
        if internal_id in self.library_tree.get_children():
            self.library_tree.selection_set(internal_id)
            self._show_selected_library_entry()

    def _show_selected_library_entry(self, _event=None) -> None:
        internal_id = self._selected_library_id()
        if not internal_id:
            return
        rows = {row["internal_id"]: row for row in library_manager_rows(self.library_var.get())}
        row = rows.get(internal_id)
        if not row:
            return
        detail = [
            f"ID: {row['internal_id']}",
            f"Formula: {row['formula']}",
            f"Space group: {row['space_group_symbol'] or '?'} (No. {row['space_group_number'] or '?'})",
            f"System: {row['chemical_system']}",
            f"Source: {row['source']} / {row['source_id']}",
            f"Enabled: {row['enabled']}",
            f"XRD cache: {row['xrd_cache']}",
            f"Backend: {row['backend_name'] or 'n/a'} {row['backend_version'] or ''}",
            f"Pattern fingerprint: {row['pattern_fingerprint'] or 'n/a'}",
            f"Candidate ready: {row['candidate_ready']}",
            f"Skip reasons: {row['skip_reasons'] or 'none'}",
            f"Last simulated: {row['last_simulated_time'] or 'n/a'}",
            f"Structure hash: {row['structure_hash']}",
            f"CIF path: {row['cached_file_path']}",
            f"Metadata path: {row['metadata_path']}",
            f"Original path: {row['original_file_path']}",
            f"Access note: {row['access_note']}",
        ]
        self._write_library_detail("\n".join(detail))

    def _selected_library_row(self) -> dict[str, object] | None:
        internal_id = self._selected_library_id()
        if not internal_id:
            self._write_library_detail("请先在结构库表格里选择一个条目。")
            return None
        rows = {row["internal_id"]: row for row in library_manager_rows(self.library_var.get())}
        row = rows.get(internal_id)
        if row is None:
            self._write_library_detail("选中的结构已经不在当前结构库中，请刷新。")
        return row

    def _open_selected_cif_location(self) -> None:
        row = self._selected_library_row()
        if not row:
            return
        path = Path(str(row["cached_file_path"]))
        if not path.exists():
            self._write_library_detail(f"CIF 文件不存在：{path}")
            return
        os.startfile(str(path.parent))

    def _view_selected_cif_text(self) -> None:
        row = self._selected_library_row()
        if not row:
            return
        path = Path(str(row["cached_file_path"]))
        if not path.exists():
            self._write_library_detail(f"CIF 文件不存在：{path}")
            return
        text = path.read_text(encoding="utf-8", errors="ignore")
        self._write_library_detail(f"CIF: {path}\n\n{text[:12000]}")

    def _view_selected_peak_table(self) -> None:
        row = self._selected_library_row()
        if not row:
            return
        internal_id = str(row["internal_id"])
        radiation = "CuKa"
        settings_hash = simulation_settings_hash(radiation, (10.0, 90.0))
        peaks = StructureLibrary(self.library_var.get()).load_xrd_peaks(internal_id, radiation, settings_hash)
        if not peaks:
            self._write_library_detail("这个结构还没有有效的模拟 XRD 缓存。请点“重算 XRD”。")
            return
        lines = [f"模拟峰表：{row['formula']} / {internal_id}", "", "2theta    intensity    hkl    d"]
        for peak in peaks[:160]:
            lines.append(f"{peak.two_theta:7.3f}  {peak.relative_intensity:9.3f}  {peak.hkl or ''}  {'' if peak.d_spacing is None else f'{peak.d_spacing:.4f}'}")
        self._write_library_detail("\n".join(lines))

    def _recompute_selected_xrd(self) -> None:
        row = self._selected_library_row()
        if not row:
            return
        library = StructureLibrary(self.library_var.get())
        summary = build_library_xrd_cache(library, structure_ids=[str(row["internal_id"])], force=True)
        self._refresh_library_table()
        self._write_library_detail(
            "选中结构 XRD 重算完成。\n"
            f"checked: {summary['checked']}\n"
            f"simulated: {summary['simulated']}\n"
            f"cached: {summary['cached']}\n"
            f"failed: {summary['failed']}\n"
            + "\n".join(str(item) for item in summary.get("warnings", []))
        )

    def _inspect_peak(self) -> None:
        try:
            theta = float(self.inspect_theta_var.get().strip())
        except ValueError:
            self._write_inspector("请输入有效的实验峰 2theta 数值。")
            return
        try:
            evidence = inspect_peak_from_app(self.library_var.get(), theta, tolerance_deg=0.2)
        except Exception as exc:
            self._write_inspector(f"Peak Inspector 无法读取结构库：{exc}")
            return
        lines = [
            f"Experimental 2theta: {evidence['experimental_two_theta']:.4f}",
            f"Tolerance: ±{evidence['tolerance_deg']:.3f} deg",
            f"Status: {evidence['status']}",
            "",
            "Nearby simulated peaks:",
        ]
        nearby = evidence["nearby_peaks"]
        if not nearby:
            lines.append("- none")
        for item in nearby[:20]:
            hkl = item["hkl"] or "n/a"
            lines.append(
                f"- {item['formula']} | {item['source']} | {item['structure_internal_id']} | "
                f"calc {item['calculated_two_theta']:.4f} | delta {item['delta']:.4f} | "
                f"I {item['theoretical_relative_intensity']:.1f} | hkl {hkl} | "
                f"enabled {item['enabled']} | strong peaks {item['strong_theoretical_peak_count']}"
            )
            lines.append(f"  CIF: {item['cif_path']}")
            if item["top_strong_theoretical_peaks"]:
                preview = ", ".join(f"{p['two_theta']:.2f}/{p['intensity']:.0f}" for p in item["top_strong_theoretical_peaks"][:5])
                lines.append(f"  Strong peak preview: {preview}")
        if evidence["warnings"]:
            lines.extend(["", "Warnings:"])
            lines.extend(f"- {warning}" for warning in evidence["warnings"])
        lines.extend(["", str(evidence["scientific_note"])])
        self._write_inspector("\n".join(lines))

    def _worker(self, run_id: int, snapshot: dict[str, object]) -> None:
        pending = None
        try:
            pending = self.analysis_store.begin()
            result = run_analysis(
                xrd_file=str(snapshot["xrd_file"]),
                candidate_dirs=[str(snapshot["candidate_dir"])],
                elements=list(snapshot["elements"]),
                extra_elements=list(snapshot["extra_elements"]),
                out_dir=pending,
                cache_path=str(snapshot["cache_path"]) or None,
                use_rietan_display=True,
            )
            if run_id != self._analysis_run_id:
                self.analysis_store.rollback()
                self.result_queue.put(("analysis_stale", {"run_id": run_id}))
                return
            _current, remapped = self.analysis_store.commit(result.outputs)
            result.outputs = remapped
            self.result_queue.put(("analysis_ok", {"run_id": run_id, "result": result}))
        except Exception as exc:
            if pending is not None:
                self.analysis_store.rollback()
            self.result_queue.put(("analysis_error", {"run_id": run_id, "error": exc}))

    def _library_worker(self) -> None:
        try:
            imported = import_local_cifs_to_library(self.dir_var.get(), self.library_var.get())
            cached = build_library_cache_from_folder(self.library_var.get())
            self.result_queue.put(("library_ok", {"imported": imported, "cached": cached}))
        except Exception as exc:
            self.result_queue.put(("library_error", exc))

    def _library_analysis_worker(self, run_id: int, snapshot: dict[str, object]) -> None:
        pending = None
        try:
            pending = self.analysis_store.begin()
            result = run_library_analysis_from_app(
                str(snapshot["xrd_file"]),
                str(snapshot["library_db"]),
                " ".join(snapshot["elements"]),
                "",
                str(pending),
                str(snapshot["target_candidate_id"]),
                bool(snapshot["scientific_safe_mode"]),
                bool(snapshot["force_recompute"]),
                True,
            )
            if run_id != self._analysis_run_id:
                self.analysis_store.rollback()
                self.result_queue.put(("analysis_stale", {"run_id": run_id}))
                return
            _current, remapped = self.analysis_store.commit(result.outputs)
            result.outputs = remapped
            self.result_queue.put(("analysis_ok", {"run_id": run_id, "result": result}))
        except Exception as exc:
            if pending is not None:
                self.analysis_store.rollback()
            self.result_queue.put(("analysis_error", {"run_id": run_id, "error": exc}))

    def _mp_harvest_worker(self) -> None:
        try:
            def on_progress(event):
                self.result_queue.put(("mp_progress", event))

            result = materials_project_harvest_from_app(
                library_db=self.library_var.get(),
                elements=self.elements_var.get(),
                impurities=self.extra_var.get(),
                mode=self.harvest_mode_var.get(),
                api_key=self.mp_api_key_var.get().strip(),
                env_path=default_env_path(),
                dry_run=False,
                progress_callback=on_progress,
            )
            self.result_queue.put(("mp_ok", result))
        except Exception as exc:
            self.result_queue.put(("mp_error", exc))

    def _oxide_supplement_worker(self, snapshot: dict[str, object] | None = None) -> None:
        try:
            data = snapshot or {
                "library_path": self.library_var.get(),
                "elements": set(parse_element_symbols(self.elements_var.get())),
                "api_key": self.mp_api_key_var.get().strip(),
            }

            def on_progress(event):
                self.result_queue.put(("oxide_supplement_progress", event))

            provider = MaterialsProjectProvider(
                api_key=str(data.get("api_key") or "") or None,
                env_path=default_env_path(),
            )
            result = supplement_common_oxide_library(
                str(data["library_path"]),
                set(data["elements"]),
                provider=provider,
                radiation="CuKa",
                two_theta_range=(10.0, 90.0),
                progress_callback=on_progress,
            )
            self.result_queue.put(("oxide_supplement_ok", result))
        except Exception as exc:
            self.result_queue.put(("oxide_supplement_error", exc))

    def _cache_worker(self) -> None:
        try:
            cache = CandidateCache(self.cache_var.get())
            allowed = set(self.elements_var.get().split()) | set(self.extra_var.get().split())
            summary = build_candidate_cache(
                cache,
                [self.dir_var.get()],
                allowed,
                radiation="CuKa",
                two_theta_range=(10.0, 90.0),
            )
            self.result_queue.put(("cache_ok", {"summary": summary, "stats": cache.stats()}))
        except Exception as exc:
            self.result_queue.put(("cache_error", exc))

    def _poll_result(self) -> None:
        self._poll_after_id = None
        if self._closed:
            return
        try:
            kind, payload = self.result_queue.get_nowait()
        except queue.Empty:
            self._schedule_result_poll()
            return
        if kind == "mp_progress":
            self._show_mp_progress(payload)
            self._schedule_result_poll()
            return
        if kind == "oxide_supplement_progress":
            self.status_var.set(f"正在补齐氧化物：{payload.get('formula', '')}")
            self._write_log(
                f"氧化物补库：{payload.get('formula', '')}，检索 {payload.get('retrieved', 0)}，"
                f"导入 {payload.get('imported', 0)}。"
            )
            self._schedule_result_poll()
            return
        if kind in {"analysis_ok", "analysis_error", "analysis_stale"}:
            payload_run_id = int(payload.get("run_id", -1))
            if payload_run_id != self._analysis_run_id:
                self._schedule_result_poll()
                return
            if kind == "analysis_stale":
                self._schedule_result_poll()
                return
        self.running = False
        self.run_button.state(["!disabled"])
        self.cache_button.state(["!disabled"])
        self.library_button.state(["!disabled"])
        self.library_run_button.state(["!disabled"])
        if kind == "analysis_ok":
            self.force_recompute_var.set(False)
            self._show_result(payload["result"])
        elif kind == "cache_ok":
            self._show_cache_result(payload)
        elif kind == "library_ok":
            self._show_library_result(payload)
        elif kind == "library_error":
            self.status_var.set("结构库导入失败。请检查 CIF 文件夹。")
            guidance = self._friendly_error(payload)
            self._write_summary(guidance)
            self._write_log(guidance)
            messagebox.showerror("结构库导入失败", guidance)
        elif kind == "mp_ok":
            self._show_mp_result(payload)
        elif kind == "mp_error":
            self.status_var.set("Materials Project 采集失败。")
            guidance = self._friendly_error(payload)
            self._write_log(guidance)
            messagebox.showerror("MP 采集失败", guidance)
        elif kind == "oxide_supplement_ok":
            imported = int(payload.get("imported") or 0)
            requested = [str(item) for item in payload.get("requested_formulas", [])]
            missing = [str(item) for item in payload.get("missing_after", [])]
            self.status_var.set(f"常见氧化物补库完成：新增 {imported} 个结构。")
            lines = [
                "常见氧化物补库完成。",
                f"请求公式：{', '.join(requested) if requested else '无需补充'}",
                f"新增结构：{imported}",
                f"仍缺少：{', '.join(missing) if missing else '无'}",
            ]
            self._write_log("\n".join(lines))
            self._refresh_target_phase_options()
            self._refresh_library_table()
        elif kind == "oxide_supplement_error":
            self.status_var.set("常见氧化物补库失败。")
            guidance = self._friendly_error(payload)
            self._write_log(guidance)
            messagebox.showerror("氧化物补库失败", guidance)
        elif kind == "cache_error":
            self.status_var.set("缓存构建失败。请按提示检查候选来源。")
            guidance = self._friendly_error(payload)
            self._write_summary(guidance)
            self._write_log(guidance)
            messagebox.showerror("缓存构建失败", guidance)
        elif kind == "analysis_error":
            self.status_var.set("分析失败。请按提示检查输入。")
            guidance = self._friendly_error(payload["error"])
            self._write_summary(guidance)
            self._write_log(guidance)
            messagebox.showerror("分析失败", guidance)
        else:
            self.status_var.set("任务失败。请查看日志。")
            guidance = self._friendly_error(payload)
            self._write_log(guidance)
        if self.current_analysis_result is not None:
            self.save_analysis_button.state(["!disabled"])
            if getattr(self.current_analysis_result, "context", None) is not None:
                self.phase_stripping_button.state(["!disabled"])
            if self._refinement_candidates(self.current_analysis_result):
                self.refinement_button.state(["!disabled"])
        self._schedule_result_poll()

    def _schedule_result_poll(self) -> None:
        if not self._closed and self.winfo_exists():
            self._poll_after_id = self.after(120, self._poll_result)

    def destroy(self) -> None:
        self._closed = True
        if self._poll_after_id is not None:
            try:
                self.after_cancel(self._poll_after_id)
            except tk.TclError:
                pass
            self._poll_after_id = None
        try:
            super().destroy()
        except tk.TclError:
            pass

    def _show_result(self, result) -> None:
        self.current_analysis_result = result
        self.save_analysis_button.state(["!disabled"])
        if getattr(result, "context", None) is not None:
            self.phase_stripping_button.state(["!disabled"])
        if self._refinement_candidates(result):
            self.refinement_button.state(["!disabled"])
        self.status_var.set("分析完成：当前结果为临时会话，点击“保存本次分析”才会永久保存。")
        lines = [
            "分析完成。\n",
            f"实验峰数量：{len(result.experimental_peaks)}",
            f"候选数量：{len(result.candidates)}",
            "",
            "优先检查的候选相：",
        ]
        runtime = getattr(result, "scientific_runtime", {}) or {}
        if runtime:
            lines[3:3] = [
                f"计算后端：{runtime.get('backend_name', 'unknown')} {runtime.get('backend_version', '')}",
                f"辐射源：{runtime.get('radiation', 'n/a')}；波长：{runtime.get('wavelength_angstrom', 'n/a')} Å",
                f"运行模式：{'Scientific Safe Mode' if runtime.get('scientific_safe_mode') else 'Validated Cache Mode'}",
                f"本次直接计算：{runtime.get('freshly_simulated', 0)}；验证缓存命中：{runtime.get('validated_cache_hits', 0)}",
                f"结构记录：{runtime.get('structure_records_checked', 0)}；归并后相候选：{runtime.get('grouped_phase_candidates', len(result.candidates))}；失败：{runtime.get('simulation_failures', 0)}",
            ]
        if not result.top_scores:
            lines.append("- 没有得到可用候选。请检查候选文件夹、元素过滤或 CIF/峰表文件。")
        for rank, score in enumerate(result.top_scores[:6], 1):
            missing = ", ".join(f"{p.two_theta:.2f}" for p in score.missing_strong_theory_peaks[:6])
            if len(score.missing_strong_theory_peaks) > 6:
                missing += " ..."
            lines.append(
                f"{rank}. {score.candidate.formula_pretty} | score {score.score:.3f} | "
                f"matched {len(score.matched_theory_peaks)} | missing strong {len(score.missing_strong_theory_peaks)}"
            )
            if missing:
                lines.append(f"   缺失强峰：{missing}")
        lines.extend(["", f"仍未解释峰数量：{len(result.explanation.unexplained_experimental_peaks)}"])
        if result.explanation.unexplained_experimental_peaks:
            preview = ", ".join(f"{p.two_theta:.2f}" for p in result.explanation.unexplained_experimental_peaks[:14])
            lines.append(f"未解释峰预览：{preview}")
        if result.contribution and result.contribution.contributions:
            lines.extend(["", "相对贡献 proxy（不是 Rietveld 定量）："])
            for formula, value in result.contribution.contributions.items():
                lines.append(f"- {formula}: {value:.1f}%")
        trust_lines = self._simulation_validation_summary_lines(result)
        if trust_lines:
            lines.extend(["", "模拟峰可信度 / VESTA 对照：", *trust_lines])
        lines.extend(["", "临时会话文件（下次成功分析会替换）："])
        for name, path in result.outputs.items():
            lines.append(f"- {name}: {path}")
        if result.warnings:
            lines.extend(["", "注意："])
            lines.extend(f"- {warning}" for warning in result.warnings)
        lines.extend(
            [
                "",
                "科学边界：这是候选筛选，不是确定性定相；高分表示值得优先检查，不代表该相已被确认。",
            ]
        )
        self._write_summary("\n".join(lines))
        self._write_candidate_phase_panel(result)
        self._write_used_structures_panel(result)
        report_lines = ["当前临时会话文件（尚未永久保存）："]
        for name, path in result.outputs.items():
            report_lines.append(f"- {name}: {path}")
        if "xrd_presentation_plot" in result.outputs:
            report_lines.extend(
                [
                    "",
                    "Presentation 图已生成：适合汇报展示实验谱、目标模拟峰和可能杂质候选。",
                    "Diagnostic 图已生成：适合检查实验峰、候选峰和缺失峰证据。",
                ]
            )
            self._display_xrd_plot(result.outputs["xrd_presentation_plot"])
        elif "xrd_plot" in result.outputs:
            report_lines.extend(["", "图片报告已生成：xrd_plot PNG 会展示实验谱和候选峰。"])
            self._display_xrd_plot(result.outputs["xrd_plot"])
        elif "xrd_plot_clean" in result.outputs:
            report_lines.extend(["", "图片报告已生成：xrd_plot_clean PNG 会展示实验峰和候选峰位置。"])
            self._display_xrd_plot(result.outputs["xrd_plot_clean"])
        self._write_report_paths("\n".join(report_lines))
        log_lines = [
            "分析完成。",
            f"缓存路径：{self.cache_var.get()}",
            f"候选数：{len(result.candidates)}",
            f"实验峰数：{len(result.experimental_peaks)}",
        ]
        if result.warnings:
            log_lines.extend(["", "警告：", *[f"- {warning}" for warning in result.warnings]])
        self._write_log("\n".join(log_lines))

    def _simulation_validation_summary_lines(self, result) -> list[str]:
        counts: dict[str, int] = {}
        failed: list[str] = []
        for score in result.top_scores:
            status = str(score.candidate.simulation_validation.get("status") or "not_checked")
            counts[status] = counts.get(status, 0) + 1
            if status == "failed":
                failed.append(score.candidate.formula_pretty)
        if not counts:
            return []
        order = ["reference_source", "passed", "failed", "no_reference", "not_checked"]
        lines = [", ".join(f"{key}: {counts[key]}" for key in order if key in counts)]
        if failed:
            lines.append("VESTA 对照失败候选：" + ", ".join(failed[:8]))
        return lines

    def _write_candidate_phase_panel(self, result) -> None:
        for item in self.candidate_phase_tree.get_children():
            self.candidate_phase_tree.delete(item)
        self.candidate_phase_details = {}
        phases = {phase.phase_id: phase for phase in getattr(result, "phase_candidates", [])}
        if not result.top_scores:
            self._write_text_widget(self.candidate_phase_text, "没有可显示的候选相。")
            return
        for rank, score in enumerate(result.top_scores[:20], 1):
            components = score.score_components or {}
            phase = phases.get(score.candidate.candidate_id)
            parent_id = f"phase-row-{rank}"
            space_group = (
                f"{score.candidate.space_group_symbol or '?'} "
                f"(No. {score.candidate.space_group_number or '?'})"
            )
            self.candidate_phase_tree.insert(
                "",
                "end",
                iid=parent_id,
                text=f"#{rank}",
                open=rank == 1,
                values=(
                    score.candidate.formula_pretty,
                    space_group,
                    f"{score.score:.4f}",
                    score.candidate.simulation_state,
                    score.candidate.source,
                ),
            )
            pattern = score.candidate.simulated_pattern
            detail_lines = [
                f"PhaseCandidate: {score.candidate.candidate_id}",
                f"Formula: {score.candidate.formula_pretty}",
                f"Space group: {space_group}",
                f"Score: {score.score:.4f} ({components.get('confidence_label', 'n/a')})",
                f"Simulation: {score.candidate.simulation_state}",
                f"Source: {score.candidate.source}",
                "   components: "
                f"matched fraction {components.get('matched_strong_intensity_fraction', 'n/a')}, "
                f"missing penalty {components.get('missing_strong_penalty', 'n/a')}, "
                f"multi-peak bonus {components.get('multi_peak_bonus', 'n/a')}",
                f"   matched peaks {len(score.matched_theory_peaks)}, "
                f"missing strong peaks {len(score.missing_strong_theory_peaks)}",
            ]
            if pattern is not None:
                detail_lines.extend(
                    [
                        f"Backend: {pattern.engine_name} {pattern.engine_version}",
                        f"CIF SHA256: {pattern.cif_sha256}",
                        f"Pattern fingerprint: {pattern.pattern_fingerprint}",
                        f"Exact CIF: {pattern.cif_path}",
                    ]
                )
            self.candidate_phase_details[parent_id] = "\n".join(detail_lines)
            if phase is None:
                continue
            for entry_index, entry in enumerate(phase.entries, 1):
                child_id = f"{parent_id}-entry-{entry_index}"
                cif_path = Path(self.library_var.get()).parent / entry.cached_file_path
                child_space_group = f"{entry.space_group_symbol or '?'} (No. {entry.space_group_number or '?'})"
                self.candidate_phase_tree.insert(
                    parent_id,
                    "end",
                    iid=child_id,
                    text="source",
                    values=(
                        entry.internal_id,
                        child_space_group,
                        "",
                        phase.simulation_state if entry.internal_id == phase.representative.internal_id else "equivalent_record",
                        f"{entry.source} / {entry.source_id}",
                    ),
                )
                self.candidate_phase_details[child_id] = "\n".join(
                    [
                        f"StructureEntry: {entry.internal_id}",
                        f"Formula: {entry.reduced_formula or entry.formula}",
                        f"Source: {entry.source} / {entry.source_id}",
                        f"Space group: {child_space_group}",
                        f"Structure hash: {entry.structure_hash}",
                        f"Enabled: {entry.enabled_for_matching}",
                        f"Exact CIF: {cif_path}",
                        f"Original source: {entry.original_file_path or 'n/a'}",
                    ]
                )
        first = self.candidate_phase_tree.get_children()
        if first:
            self.candidate_phase_tree.selection_set(first[0])
            self._show_candidate_phase_detail()

    def _write_used_structures_panel(self, result) -> None:
        rows = []
        try:
            import csv

            usage_path = result.outputs.get("candidate_usage_summary")
            if usage_path and Path(usage_path).exists():
                with Path(usage_path).open(encoding="utf-8") as handle:
                    rows = list(csv.DictReader(handle))
        except Exception:
            rows = []
        if not rows:
            self._write_text_widget(self.used_structures_text, "没有找到 candidate_usage_summary.csv。")
            return
        used = [row for row in rows if row.get("used_in_analysis") == "True"]
        skipped = [row for row in rows if row.get("used_in_analysis") != "True"]
        lines = [
            "结构使用摘要",
            "",
            f"总结构：{len(rows)}",
            f"实际用于匹配：{len(used)}",
            f"跳过：{len(skipped)}",
            "",
            "Used structures:",
        ]
        for row in used[:40]:
            lines.append(
                f"- {row.get('formula')} | {row.get('internal_id')} | {row.get('source')} | "
                f"cache {row.get('cache_status')} | path {row.get('structure_path')}"
            )
        if len(used) > 40:
            lines.append(f"... 还有 {len(used) - 40} 个 used 结构，详见 candidate_usage_summary.csv")
        lines.extend(["", "Skipped preview:"])
        for row in skipped[:40]:
            lines.append(f"- {row.get('formula')} | {row.get('internal_id')} | {row.get('skip_reasons')}")
        if len(skipped) > 40:
            lines.append(f"... 还有 {len(skipped) - 40} 个 skipped 结构，详见 candidate_usage_summary.csv")
        self._write_text_widget(self.used_structures_text, "\n".join(lines))

    def _display_xrd_plot(self, path: str) -> None:
        plot_path = Path(path)
        if not plot_path.exists():
            self.current_plot_path = None
            self.plot_canvas.configure(image="", text="XRD 图没有生成。请查看日志。")
            self.save_plot_button.state(["disabled"])
            return
        try:
            from PIL import Image, ImageTk

            with Image.open(plot_path) as source_image:
                image = source_image.copy()
            image.thumbnail((760, 520))
            self.plot_image = ImageTk.PhotoImage(image)
            image.close()
            self.plot_canvas.configure(image=self.plot_image, text="")
            self.current_plot_path = str(plot_path)
            self.save_plot_button.state(["!disabled"])
            for tab_id in self.tabs.tabs():
                if self.tabs.tab(tab_id, "text") == "XRD 图":
                    self.tabs.select(tab_id)
                    break
        except Exception as exc:
            self.current_plot_path = str(plot_path)
            self.plot_canvas.configure(image="", text=f"图片预览失败，但文件已生成：\n{plot_path}\n{exc}")
            self.save_plot_button.state(["!disabled"])

    def _save_current_analysis(self) -> None:
        if self.current_analysis_result is None or not self.analysis_store.current_path.is_dir():
            messagebox.showwarning("没有可保存的分析", "请先完成一次候选筛选或结构库分析。")
            return
        initial = Path(self.out_var.get()) if self.out_var.get().strip() else project_root() / "results"
        if not initial.is_dir():
            initial = initial.parent if initial.parent.is_dir() else project_root() / "results"
        parent = filedialog.askdirectory(title="选择本次分析的保存位置", initialdir=str(initial), parent=self)
        if not parent:
            return
        stem = re.sub(r"[^A-Za-z0-9\u4e00-\u9fff]+", "_", Path(self.xrd_var.get()).stem).strip("_") or "ophi_xrd"
        base = Path(parent) / f"{stem}_analysis_{time.strftime('%Y%m%d_%H%M%S')}"
        target = base
        index = 2
        while target.exists():
            target = base.with_name(f"{base.name}_{index}")
            index += 1
        try:
            saved = self.analysis_store.save_current(target)
        except Exception as exc:
            messagebox.showerror("保存分析失败", str(exc))
            self.status_var.set("保存失败；临时分析仍然保留。")
            return
        self.last_saved_analysis_path = saved
        self.out_var.set(str(Path(parent)))
        self.status_var.set(f"本次分析已显式保存：{saved}")
        self._write_log(f"本次分析已显式保存，未重新运行模拟：\n{saved}")
        messagebox.showinfo("保存完成", f"完整分析与 manifest 已保存到：\n{saved}")

    def _phase_stripping_is_open(self) -> bool:
        if self.phase_stripping_window is None:
            return False
        try:
            return bool(self.phase_stripping_window.winfo_exists())
        except tk.TclError:
            return False

    def _open_phase_stripping(self) -> None:
        result = self.current_analysis_result
        context = getattr(result, "context", None) if result is not None else None
        if context is None:
            messagebox.showwarning("还没有可用分析", "请先完成一次候选筛选或结构库分析。")
            return
        candidates = [candidate for candidate in result.candidates if candidate.theory_peaks]
        if not candidates:
            messagebox.showwarning("没有可剥离候选", "当前分析没有带规范模拟峰的候选相。")
            return
        self.phase_stripping_window = open_or_raise_phase_stripping_window(
            self,
            self.phase_stripping_window,
            context=context,
            candidates=candidates,
            element_scope=parse_element_symbols(self.elements_var.get()),
            on_closed=self._phase_stripping_closed,
        )
        self.status_var.set("手动物相剥离窗口已打开；主分析结果不会在窗口打开时被替换。")

    def _phase_stripping_closed(self) -> None:
        self.phase_stripping_window = None

    def _refinement_is_open(self) -> bool:
        if self.refinement_window is None:
            return False
        try:
            return bool(self.refinement_window.winfo_exists())
        except tk.TclError:
            return False

    @staticmethod
    def _refinement_candidates(result) -> list:
        if result is None or getattr(result, "context", None) is None:
            return []
        ordered = [score.candidate for score in getattr(result, "top_scores", [])]
        ordered.extend(getattr(result, "candidates", []))
        candidates = []
        seen: set[str] = set()
        for candidate in ordered:
            if candidate.candidate_id in seen or Path(candidate.source_path).suffix.lower() != ".cif":
                continue
            seen.add(candidate.candidate_id)
            candidates.append(candidate)
        return candidates

    def _open_refinement(self) -> None:
        result = self.current_analysis_result
        context = getattr(result, "context", None) if result is not None else None
        candidates = self._refinement_candidates(result)
        if context is None or not candidates:
            messagebox.showwarning("还不能精修", "请先完成分析，并确保目标候选来自 CIF。")
            return
        backend = RietanRefinementBackend(
            vesta_exe=self.vesta_exe_var.get().strip() or None,
            rietan_exe=self.rietan_exe_var.get().strip() or None,
        )
        supporting_candidates = [
            score.candidate for score in getattr(result, "top_scores", [])[:12]
        ]
        self.refinement_window = open_or_raise_refinement_window(
            self,
            self.refinement_window,
            context=context,
            candidates=candidates,
            backend=backend,
            supporting_candidates=supporting_candidates,
            oxide_library_path=self.library_var.get().strip() or None,
            on_closed=self._refinement_closed,
        )
        self.status_var.set("RIETAN-FP 受约束精修窗口已打开；主分析结果保持不变。")

    def _refinement_closed(self) -> None:
        self.refinement_window = None

    def _save_current_plot(self) -> None:
        if not self.current_plot_path or not Path(self.current_plot_path).exists():
            messagebox.showwarning("没有可保存图片", "请先完成一次分析，生成 XRD 图预览。")
            return
        default_name = f"{Path(self.xrd_var.get()).stem or 'ophi_xrd'}_dashboard_xrd.png"
        target = filedialog.asksaveasfilename(
            initialdir=str(
                first_existing_directory(self.out_var.get(), fallback=project_root() / "results")
            ),
            initialfile=default_name,
            defaultextension=".png",
            filetypes=[("PNG image", "*.png"), ("All files", "*.*")],
        )
        if not target:
            return
        shutil.copy2(self.current_plot_path, target)
        self.status_var.set("XRD 图已保存。")
        self._write_log(f"XRD 图已保存：{target}")

    def _show_cache_result(self, payload: dict[str, object]) -> None:
        summary = payload["summary"]
        stats = payload["stats"]
        self.status_var.set("缓存构建完成：后续筛选会优先复用缓存。")
        lines = [
            "缓存构建/更新完成。",
            "",
            f"扫描候选：{summary['scanned_candidates']}",
            f"本次新增/更新峰表：{summary['stored_patterns']}",
            f"缓存候选总数：{stats['candidates']}",
            f"缓存峰表总数：{stats['patterns']}",
            f"缓存路径：{self.cache_var.get()}",
        ]
        warnings = summary.get("warnings") or []
        if warnings:
            lines.extend(["", "注意："])
            lines.extend(f"- {warning}" for warning in warnings)
        self._write_summary("\n".join(lines))
        self._write_log("\n".join(lines))

    def _show_library_result(self, payload: dict[str, object]) -> None:
        imported = payload["imported"]
        cached = payload["cached"]
        self.status_var.set("结构库已更新：本地 CIF 已导入并尝试构建库级 XRD 缓存。")
        lines = [
            "结构库导入完成。",
            "",
            f"扫描 CIF：{imported['scanned']}",
            f"导入结构：{imported['imported']}",
            f"跳过结构：{imported['skipped']}",
            f"库级 XRD 检查：{cached['checked']}",
            f"新模拟缓存：{cached['simulated']}",
            f"已有缓存：{cached['cached']}",
            f"失败：{cached['failed']}",
            f"结构库：{self.library_var.get()}",
        ]
        warnings = list(imported.get("warnings") or []) + list(cached.get("warnings") or [])
        if warnings:
            lines.extend(["", "注意："])
            lines.extend(f"- {warning}" for warning in warnings)
        self._write_summary("\n".join(lines))
        self._write_log("\n".join(lines))
        self._refresh_library_table()
        self._refresh_target_phase_options()

    def _show_mp_result(self, payload: dict[str, object]) -> None:
        self.status_var.set("Materials Project 采集完成。")
        lines = [
            "Materials Project 采集完成。",
            "",
            f"搜索系统数：{len(payload.get('chemical_systems', []))}",
            f"检索记录：{payload.get('retrieved', 0)}",
            f"新增结构：{payload.get('imported', 0)}",
            f"跳过重复：{payload.get('skipped_duplicates', 0)}",
            f"失败：{payload.get('failed', 0)}",
        ]
        cache = payload.get("xrd_cache") or {}
        if cache:
            lines.extend(["", f"XRD 缓存新模拟：{cache.get('simulated', 0)}", f"XRD 缓存已有：{cache.get('cached', 0)}"])
        warnings = payload.get("warnings") or []
        if warnings:
            lines.extend(["", "注意："])
            lines.extend(f"- {warning}" for warning in warnings)
        self._write_summary("\n".join(lines))
        self._write_log("\n".join(lines))
        self._refresh_library_table()
        self._refresh_target_phase_options()

    def _show_mp_progress(self, payload: dict[str, object]) -> None:
        line = (
            f"MP 采集中：{payload.get('system')} | "
            f"检索 {payload.get('retrieved', 0)} | "
            f"新增 {payload.get('imported', 0)} | "
            f"跳过重复 {payload.get('skipped_duplicates', 0)} | "
            f"失败 {payload.get('failed', 0)} | "
            f"累计新增 {payload.get('total_imported', 0)}"
        )
        self.status_var.set(line)
        self._write_log(line)
        self._refresh_library_table()

    def _friendly_error(self, error: object) -> str:
        text = str(error)
        if "Materials Project API key is missing" in text or "MP_API_KEY" in text:
            return (
                "Materials Project 还没有可用 API key。\n"
                "请在“免费公开结构采集”输入框填入你的 MP key 后点“保存 Key”，或在本地 .env 里写入 MP_API_KEY。"
            )
        if "mp-api" in text:
            return "当前 Python 环境缺少 Materials Project 官方 mp-api 包。请在 Ophi 使用的 Anaconda 环境中安装 mp-api 后重试。"
        if "Materials Project harvest failed" in text:
            return f"Materials Project 在线采集失败：{text}\n请检查网络、API key 权限，以及本次元素系统是否过宽。"
        if "does not contain readable" in text:
            return "实验谱文件无法读取。请确认它是 Rigaku .asc 或两列 2theta/intensity 文本。"
        if "No such file" in text or "not found" in text:
            return "有文件路径不存在。请重新选择实验谱或候选文件夹。"
        if "CIF" in text or "cif" in text:
            return f"CIF 解析或模拟失败：{text}\n可以先换用已有 .int/.csv 模拟峰文件夹。"
        return f"分析过程中出现问题：{text}\n请检查输入文件、候选目录和元素范围。"

    def _open_output(self) -> None:
        path = self.last_saved_analysis_path or (Path(self.out_var.get()) if self.out_var.get().strip() else None)
        if path is None or not path.exists():
            messagebox.showwarning("保存位置不存在", "还没有已保存的分析目录。请先点击“保存本次分析”。")
            return
        os.startfile(path)

    def _write_summary(self, text: str) -> None:
        self._write_text_widget(self.summary, text)

    def _write_report_paths(self, text: str) -> None:
        self._write_text_widget(self.report_paths, text)

    def _write_log(self, text: str) -> None:
        self._write_text_widget(self.log_text, text)

    def _write_library_detail(self, text: str) -> None:
        self._write_text_widget(self.library_detail, text)

    def _write_inspector(self, text: str) -> None:
        self._write_text_widget(self.inspector_text, text)

    def _write_text_widget(self, widget: tk.Text, text: str) -> None:
        widget.delete("1.0", "end")
        widget.insert("end", text)
        widget.see("1.0")


def launch_app() -> None:
    OphiuchusApp().mainloop()
