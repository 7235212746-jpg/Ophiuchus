from __future__ import annotations

import argparse
import importlib
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
import platform
import sys

from .help_support import MANUAL_FILENAME
from .runtime import is_frozen, resource_root, user_data_root
from .xrd.rietan_backend import (
    discover_cif2ins_executable,
    discover_multiphase_template,
    discover_rietan_executable,
    discover_vesta_executable,
)


_DEPENDENCIES = {
    "numpy": "numpy",
    "scipy": "scipy",
    "matplotlib": "matplotlib",
    "pymatgen": "pymatgen",
    "mp-api": "mp_api",
}


def _dependency_report() -> dict[str, dict[str, object]]:
    report: dict[str, dict[str, object]] = {}
    for distribution, module_name in _DEPENDENCIES.items():
        try:
            importlib.import_module(module_name)
            package_version = version(distribution)
            report[distribution] = {"available": True, "version": package_version}
        except (ImportError, PackageNotFoundError) as exc:
            report[distribution] = {"available": False, "error": str(exc)}
    return report


def _check_user_data_writable(path: Path) -> bool:
    try:
        path.mkdir(parents=True, exist_ok=True)
        probe = path / ".portable_write_probe"
        probe.write_text("ok", encoding="ascii")
        probe.unlink()
        return True
    except OSError:
        return False


def portable_health_report() -> dict[str, object]:
    dependencies = _dependency_report()
    try:
        import tkinter

        tkinter_version: float | None = float(tkinter.TkVersion)
    except ImportError:
        tkinter_version = None
    rietan = discover_rietan_executable()
    manual = resource_root() / "docs" / MANUAL_FILENAME
    data_root = user_data_root()
    core_ready = tkinter_version is not None and all(
        bool(item.get("available")) for item in dependencies.values()
    )
    return {
        "core_ready": core_ready,
        "frozen": is_frozen(),
        "python": platform.python_version(),
        "tkinter": tkinter_version,
        "dependencies": dependencies,
        "resource_root": str(resource_root()),
        "manual_available": manual.is_file(),
        "user_data_root": str(data_root),
        "user_data_writable": _check_user_data_writable(data_root),
        "optional_engines": {
            "vesta": discover_vesta_executable() is not None,
            "rietan": rietan is not None,
            "cif2ins": discover_cif2ins_executable(rietan_exe=rietan) is not None,
            "multiphase_template": discover_multiphase_template(rietan_exe=rietan) is not None,
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--health-check")
    args, unknown = parser.parse_known_args(argv)
    if args.health_check:
        target = Path(args.health_check).expanduser().resolve()
        target.parent.mkdir(parents=True, exist_ok=True)
        report = portable_health_report()
        target.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
        return 0 if all(
            bool(report.get(key)) for key in ("core_ready", "manual_available", "user_data_writable")
        ) else 1
    if unknown:
        raise SystemExit(f"Unsupported portable launcher arguments: {' '.join(unknown)}")
    from .app import launch_app

    launch_app()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
