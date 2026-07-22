import os
import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_all, collect_data_files, collect_submodules, copy_metadata


root = Path(SPECPATH).resolve().parent
portable_site_packages = Path(os.environ["OPHI_PORTABLE_SITE_PACKAGES"]).resolve()
sys.path.insert(0, str(root))
from tools.portable_build_helpers import discover_python_modules

datas = [
    (str(root / "docs" / "Ophiuchus_操作手册.md"), "docs"),
    (str(root / "README.md"), "."),
    (str(root / "VERSION"), "."),
]
binaries = []
hiddenimports = []

for package in ("matplotlib", "pymatgen", "mp_api", "emmet"):
    try:
        datas += collect_data_files(package)
    except Exception:
        pass

monty_datas, monty_binaries, monty_hidden = collect_all("monty")
datas += monty_datas
binaries += monty_binaries
hiddenimports += monty_hidden
hiddenimports += collect_submodules(
    "mp_api.client",
    filter=lambda name: ".contribs" not in name and not name.endswith("._server_utils"),
)
hiddenimports += discover_python_modules(portable_site_packages, "pymatgen")
hiddenimports += [
    "matplotlib.backends.backend_tkagg",
]

for distribution in ("numpy", "scipy", "pandas", "matplotlib", "pymatgen", "mp-api"):
    try:
        datas += copy_metadata(distribution)
    except Exception:
        pass

a = Analysis(
    [str(root / "tools" / "OphiuchusPortableEntry.py")],
    pathex=[str(portable_site_packages), str(root)],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["pytest", "IPython", "jupyter", "notebook", "PySide6", "PyQt5", "PyQt6", "wx"],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="OphiuchusApp",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name="OphiuchusApp",
)
