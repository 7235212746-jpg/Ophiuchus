from pathlib import Path

from PyInstaller.utils.hooks import collect_all, copy_metadata


root = Path(SPECPATH).resolve().parent.parent
datas = [
    (str(root / "docs" / "Ophiuchus_操作手册.md"), "docs"),
    (str(root / "README.md"), "."),
    (str(root / "VERSION"), "."),
]
binaries = []
hiddenimports = []

for package in ("matplotlib", "pymatgen", "mp_api", "emmet", "monty"):
    try:
        package_datas, package_binaries, package_hidden = collect_all(package)
        datas += package_datas
        binaries += package_binaries
        hiddenimports += package_hidden
    except Exception:
        pass

for distribution in ("numpy", "scipy", "pandas", "matplotlib", "pymatgen", "mp-api"):
    try:
        datas += copy_metadata(distribution)
    except Exception:
        pass

a = Analysis(
    [str(root / "ophiuchus" / "portable_entry.py")],
    pathex=[str(root)],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["pytest", "IPython", "jupyter", "notebook"],
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
