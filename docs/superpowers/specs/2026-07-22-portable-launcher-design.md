# Ophiuchus Portable Launcher Design

## Goal

Build a Windows 10/11 x64 portable distribution that starts on a computer without Python, Conda, NumPy, SciPy, matplotlib, pymatgen, or mp-api installed. Missing separately licensed VESTA/RIETAN programs must disable only their dependent workflows and must never crash the base application.

## Chosen Approach

Use a PyInstaller `onedir` scientific runtime behind the existing small .NET launcher. `onedir` starts faster, is easier to inspect, and is less likely to trigger antivirus heuristics than a self-extracting `onefile` executable. A network bootstrap installer was rejected because it cannot provide offline/reproducible startup and would execute package installation on the recipient's computer.

The release layout is:

```text
Ophiuchus_Portable/
  Ophiuchus.exe
  runtime/OphiuchusApp.exe
  runtime/_internal/...
  docs/Ophiuchus_操作手册.md
  README.md
  install_desktop_shortcut.*
  SHA256SUMS.txt
```

## Runtime Boundaries

- The portable runtime contains the Python interpreter and all open Python dependencies required by Ophiuchus.
- It does not contain `.env`, API keys, app state, SQLite libraries, experiment files, exports, VESTA, RIETAN, or RIETAN example/template files.
- Frozen user state, cache, local structure library, and local configuration live under `%LOCALAPPDATA%\Ophiuchus`.
- Bundled documentation is read from the frozen resource directory.
- Source-tree execution retains the existing project-relative development paths.

## Launcher Behavior

The .NET `Ophiuchus.exe` first looks for `runtime\OphiuchusApp.exe`. If found, it launches that executable directly with the portable package root as its working directory. If absent, it falls back to `start_ophiuchus.bat` for source checkouts. If neither route exists, it shows a Windows error dialog with the missing path instead of failing silently.

## Health Check And Safety

`OphiuchusApp.exe --health-check <json-path>` imports tkinter and the required scientific libraries without opening the GUI. It writes package versions, resource-path readiness, writable user-data readiness, and optional VESTA/RIETAN capability states. The build script runs this probe before packaging and stops on any failed core dependency.

The build script uses fixed repository-owned `build/portable` and `dist/Ophiuchus_Portable` paths, deletes only those verified paths, creates `SHA256SUMS.txt`, and excludes ignored/local research files by constructing the distribution from an explicit allowlist.

## External Scientific Engines

VESTA and RIETAN remain external because their redistribution terms are separate from Ophiuchus. Their absence is reported in settings and leaves candidate/library workflows available. Multiphase quantification remains disabled until RIETAN, `cif2ins.exe`, and the validated official multiphase template are present. The existing in-app installer extracts only the official template from a user-obtained RIETAN archive.

## Verification

1. Unit tests cover frozen resource/user-data path selection and launcher preference order.
2. The health check is run from the built runtime with the Conda environment removed from `PATH`.
3. The built executable is launched once on this machine and its process/window is observed.
4. The release allowlist is inspected for secrets, personal data, local databases, and RIETAN/VESTA binaries.
5. SHA-256 hashes are generated for all top-level release artifacts and the final zip.
