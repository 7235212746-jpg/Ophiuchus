# Ophiuchus Portable Launcher Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and verify a Windows x64 Ophiuchus distribution that requires no preinstalled Python or scientific Python libraries.

**Architecture:** Freeze a dedicated GUI entry point with PyInstaller `onedir`, keep mutable state under LocalAppData, and make the existing .NET launcher prefer the frozen runtime while retaining source-checkout fallback. Build the release from an explicit allowlist and verify it through a machine-readable health probe.

**Tech Stack:** Python 3.11, PyInstaller, PowerShell, C#/.NET Framework launcher, unittest.

## Global Constraints

- Do not redistribute VESTA, RIETAN, user CIFs, experiment data, API keys, `.env`, app state, or SQLite databases.
- Preserve source-tree startup and tests.
- Frozen cache, config, app state, and library files must live under `%LOCALAPPDATA%\Ophiuchus`.
- Build output must be `onedir`, not self-extracting `onefile`.
- A failed core health check must prevent release packaging.

---

### Task 1: Frozen Runtime Paths

**Files:**
- Create: `ophiuchus/runtime.py`
- Modify: `ophiuchus/app.py`
- Modify: `ophiuchus/help_support.py`
- Test: `tests/test_runtime_paths.py`

**Interfaces:**
- Produces: `is_frozen()`, `resource_root()`, `user_data_root()`.

- [ ] Write tests that patch `sys.frozen`, `sys._MEIPASS`, and `LOCALAPPDATA` and assert frozen resources and mutable data use different roots.
- [ ] Run `python -m unittest tests.test_runtime_paths -v` and confirm the missing-module failure.
- [ ] Implement the path helpers and route cache/library/env/app-state defaults through `user_data_root()`.
- [ ] Route the bundled manual through `resource_root()` while preserving explicit-root tests.
- [ ] Run runtime, app-experience, help-support, and path tests.
- [ ] Commit `Prepare frozen runtime paths`.

### Task 2: Portable Entry And Health Probe

**Files:**
- Create: `ophiuchus/portable_entry.py`
- Create: `tests/test_portable_entry.py`

**Interfaces:**
- Produces: `portable_health_report()` and `main(argv=None)`.

- [ ] Write a failing test asserting the report contains core dependency versions, writable data status, and separate optional engine states without secrets.
- [ ] Run the focused test and confirm the missing-module failure.
- [ ] Implement the health report and `--health-check <path>` JSON writer; default execution calls `launch_app()`.
- [ ] Run the focused tests and a source-mode health probe.
- [ ] Commit `Add portable runtime health probe`.

### Task 3: Launcher Preference And Packaging Definition

**Files:**
- Modify: `tools/OphiuchusLauncher.cs`
- Create: `packaging/OphiuchusPortable.spec`
- Create: `build_portable.ps1`
- Create: `build_portable.bat`
- Modify: `.gitignore`
- Modify: `tests/test_share_readiness.py`

**Interfaces:**
- Consumes: `ophiuchus.portable_entry`.
- Produces: `dist/Ophiuchus_Portable/Ophiuchus.exe`, `runtime/OphiuchusApp.exe`, and `SHA256SUMS.txt`.

- [ ] Write failing share-readiness tests for frozen-runtime preference, explicit release allowlist, and excluded local files.
- [ ] Run the focused tests and confirm the assertions fail.
- [ ] Update the C# launcher to start `runtime\OphiuchusApp.exe` first and retain batch fallback.
- [ ] Add a PyInstaller spec collecting tkinter, matplotlib, scipy, pymatgen, mp-api, package metadata, and Ophiuchus source/resources.
- [ ] Add guarded PowerShell/batch builders that create only approved output paths and generate SHA-256 manifests.
- [ ] Run the focused tests and rebuild the .NET launcher.
- [ ] Commit `Build self-contained portable launcher`.

### Task 4: Real Portable Build Gate

**Files:**
- Modify: `README.md`
- Modify: `docs/Ophiuchus_操作手册.md`

**Interfaces:**
- Consumes: PyInstaller output and health probe.
- Produces: verified portable folder and zip.

- [ ] Install PyInstaller only in the local build environment and record its version in build output.
- [ ] Run `build_portable.ps1`; stop if the frozen health probe fails.
- [ ] Run the frozen probe again with a minimal system `PATH` and inspect its JSON.
- [ ] Launch the built top-level `Ophiuchus.exe`, verify a visible responsive window, then close it normally.
- [ ] Scan the package for `.env`, SQLite, API-key patterns, user paths, VESTA/RIETAN executables, and research-data extensions outside documentation.
- [ ] Document portable startup, optional-engine setup, size expectations, and data locations.
- [ ] Run the full test suite and `git diff --check`.
- [ ] Commit `Verify portable Ophiuchus distribution`.
