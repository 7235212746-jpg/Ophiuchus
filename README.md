# Ophiuchus

Current version: v0.2.0-dev

中文操作说明见 [`docs/Ophiuchus_操作手册.md`](docs/Ophiuchus_操作手册.md)。软件内也可通过左侧“帮助与关于”打开手册。

Maintainer contact: `wanyc@issp.u-tokyo.ac.jp`

For a computer without Python or scientific libraries, extract `Ophiuchus_Portable.zip` completely and double-click `Ophiuchus.exe`. Keep the adjacent `runtime` folder in place; the portable build contains its own pinned Python, NumPy, SciPy, Matplotlib, pymatgen, and mp-api runtime. User settings, API credentials, the local structure library, caches, and deliberately saved results are kept under `%LOCALAPPDATA%\Ophiuchus`, not inside the release folder.

VESTA, RIETAN-FP, and `cif2ins` are optional external crystallographic engines with separate distribution terms and are not included in the portable archive. Their dependent controls remain unavailable or explicitly fall back until the user configures legal local copies. Core import, candidate screening, pymatgen XRD simulation, plotting, and Materials Project access do not require a local Python installation.

For source development on Windows, double-click `Ophiuchus.exe` to launch without a console window. The project folder can live anywhere. Run `install_desktop_shortcut.bat` after cloning or moving the folder to create a desktop shortcut that resolves the current path at runtime. Rebuild the launcher with `build_launcher_exe.bat` after changing its wrapper source. Build the self-contained archive with `build_portable.bat`; its dependency versions are pinned in `requirements-portable-lock.txt`, and the build fails unless the frozen runtime passes its health check.

Ophiuchus is a local-first materials research workflow project. Phase 1 is the Ophi XRD Candidate Screener. Phase 2 is building the local structure library and evidence workbench foundation.

## What v0.1.0 Does

- Imports experimental XRD files, including Rigaku-style `.asc` and common two-column text files.
- Extracts experimental peaks without requiring scipy.
- Scans local candidate folders for `.cif` files and simple peak-list files such as `.int`, `.txt`, `.xy`, `.dat`, and `.csv`.
- Simulates CIF-derived candidate patterns through the canonical `ValidatedXRDBackend`; a simulation failure is reported instead of silently switching to an unvalidated lightweight fallback.
- Scores candidate phases by matched theoretical peaks and missing strong theoretical peaks.
- Builds a small greedy multi-phase explanation.
- Writes `results.json`, `report.md`, `extracted_peaks.csv`, `top_candidates.csv`, and `peak_assignments.csv`.
- Provides both a CLI and a simple desktop window.

## Phase 2 Library Foundation

The current development build adds the first pieces of the Phase 2 local structure library:

- SQLite local structure library at `data\ophi_library.sqlite`.
- Legal local CIF import through `library-import` and the desktop "导入结构库" button.
- Provider architecture for local, Materials Project, COD, AFLOW, OQMD, NOMAD/OPTIMADE, and ICSD manual/API positions.
- Materials Project public API harvesting with user-provided API key, local CIF/metadata export, duplicate skipping, and XRD cache building.
- Library-level simulated XRD cache for imported structures.
- Deterministic export naming helpers and local Peak Inspector evidence lookup.
- Phase Evidence Card storage foundation.
- Desktop Library Manager tab with structure table, enabled/disabled toggle, source/access notes, and XRD cache status.
- Desktop Peak Inspector tab that checks a typed experimental 2theta against local cached simulated peaks.
- Desktop workbench shell with left navigation, center workflow controls, and right evidence/results panels.
- Centralized light UI theme with Chinese-capable fonts, soft gray-blue background, white panels, subtle borders, and wider right-side analysis tabs.
- Main-screen API/database connection entry. Materials Project can be configured now; COD, AFLOW, OQMD, and NOMAD/OPTIMADE are visible placeholders for later provider work.
- Right-side XRD plot dashboard showing experimental XRD, the target simulated pattern, and up to three likely impurity/candidate patterns. Use the adjacent save button when you want to export the preview image.
- Modal 118-element periodic table for building the analysis scope, with path/formula mismatch protection and no silent target-phase replacement.
- Independent three-pane manual phase-stripping and signed-residual workspace with bounded preview, accept, exclude, undo, redo, reset, and explicit export.
- Transactional GUI output: one temporary current analysis is retained locally and becomes permanent only after `保存本次分析` is clicked.

Local CIF import works offline. Materials Project harvesting works when `mp-api` is installed, network access is available, and you provide your own free MP API key. Restricted databases such as ICSD remain manual/legal-import only; Ophi does not scrape restricted sources.

## Scientific Boundary

Ophi XRD Candidate Screener performs peak-based candidate screening. It is not a replacement for manual crystallographic judgment, Rietveld refinement, or authorized PDF card databases. A high score means that a candidate phase is worth checking; it does not prove the phase is present.

## Run The Desktop App

```powershell
cd C:\path\to\Ophiuchus
python -m ophiuchus app
```

## Run The CLI

```powershell
cd C:\path\to\Ophiuchus
python -m ophiuchus analyze `
  --xrd "C:\path\to\sample.asc" `
  --cif-dir "C:\path\to\candidate_folder" `
  --elements Zr Fe Ge `
  --extra-elements O C Si Al Cu Sn Hf Sc I `
  --out-dir ".\results\sample"
```

## Structure Library CLI

Import legally obtained local CIF files:

```powershell
python -m ophiuchus library-import `
  --folder "C:\path\to\cifs" `
  --library-db ".\data\ophi_library.sqlite"
```

List library entries:

```powershell
python -m ophiuchus library-list --library-db ".\data\ophi_library.sqlite"
```

Build simulated XRD cache for library entries:

```powershell
python -m ophiuchus library-cache-xrd --library-db ".\data\ophi_library.sqlite"
```

In the desktop app, the left navigation exposes Projects, Samples, Library, XRD Analysis, Phase Evidence, and Settings. The center workspace keeps the active XRD/library controls. The right-side `结构库` tab shows imported structures, enabled state, and XRD cache status. The `Peak Inspector` tab accepts an experimental 2theta value and lists nearby local simulated peaks from the structure library cache.

Select elements with the periodic-table button and choose the target phase explicitly before analysis. After a successful run, `手动物相剥离 / 残差谱分析` opens the independent residual workspace. Preview a candidate before accepting it; Undo, Redo, and Reset always recompute from the immutable original data. See `MANUAL_PHASE_STRIPPING.md` for the equations, warnings, scientific limits, and full workflow.

Desktop analysis output is temporary by default and is limited to one transactional session under `%LOCALAPPDATA%\Ophiuchus\analysis-session`. Click `保存本次分析` to choose a permanent folder. The residual workspace has its own explicit CSV/JSON/PNG export. Existing `results` content is not deleted automatically, and CLI commands continue to write directly to their requested `--out-dir`.

The main workspace also has `连接 API / 数据库`. Choose `Materials Project`, paste your API key, save it locally, test the connection, then use `按当前元素采集` to harvest structures for the current element fields. The importer skips duplicate MP structures by provider id.

Analyze an experimental XRD file using enabled local library entries and export fixed-name outputs:

```powershell
python -m ophiuchus library-analyze `
  --xrd "C:\path\to\sample.xy" `
  --library-db ".\data\ophi_library.sqlite" `
  --elements Zr Fe Ge `
  --out-dir ".\results\sample_library" `
  --project "ZrFe6Ge4" `
  --sample-id "ZFG_980C72h" `
  --export-folder ".\exports"
```

## Phase 4 Multi-Phase Candidate Screening

The library analysis workflow is now the preferred Phase 4 path for experimental XRD screening. It uses the unified local `StructureLibrary`, respects each structure's enabled/disabled state, rebuilds or reads the library XRD cache, and exports the exact candidate usage evidence.

Key outputs include:

- `candidate_usage_summary.csv`: every library structure, source, path, hash, enabled state, cache status, used/skipped state, and skip reason.
- `top_candidates.csv`: ranked candidate phases with conservative confidence labels and score components.
- `peak_assignments.csv`: experimental peaks and candidate phase matches.
- `results.json`: full machine-readable evidence, including used structure IDs and candidate usage summary.
- `audit_*/peak_match_table.csv`: detailed per-peak assignment evidence.
- `xrd_presentation.png`: clean report-style plot for the experimental pattern, target reference/simulated peaks, and likely impurity candidates.
- `xrd_diagnostic.png`: evidence-oriented plot for checking detected peaks, candidate matches, missing peaks, and unresolved regions.

The desktop `结构库` tab functions as the Candidate Phase Panel. It shows imported/cached structures, source IDs, enabled state, XRD cache status, candidate-ready state, skip reasons, and CIF paths. Disable structures here before running analysis when they should not participate.

The `Peak Inspector` tab checks a selected experimental 2theta against enabled cached structures and reports all nearby simulated peaks, structure IDs, CIF paths, hkl, intensity, enabled state, and strong-peak context.

See `CANDIDATE_MATCHING.md` and `LIBRARY_USAGE.md` for the scoring model, traceability rules, and scientific limits.

## Online Structure Harvesting

Materials Project harvesting reads the key from the desktop app, `MP_API_KEY`, or a local `.env` file. The `.env` file is ignored by git.

```powershell
$env:MP_API_KEY="your_key_here"
```

Test the key and official `mp-api` connection:

```powershell
python -m ophiuchus mp-test
```

Preview the chemical systems Ophi will query before spending API calls:

```powershell
python -m ophiuchus mp-harvest `
  --library-db ".\data\ophi_library.sqlite" `
  --elements Zr Fe Ge `
  --impurities O C Si Al Cu `
  --mode normal `
  --dry-run
```

Harvest structures into the local library, write CIF/metadata files, skip duplicates, and build the library XRD cache:

```powershell
python -m ophiuchus mp-harvest `
  --library-db ".\data\ophi_library.sqlite" `
  --elements Zr Fe Ge `
  --impurities O C Si Al Cu `
  --mode conservative `
  --max-entries-per-system 25
```

Saved MP structures live under `data\library\structures\materials_project`; metadata lives under `data\library\metadata\materials_project`. Use conservative mode first for real samples, then broaden only when the first pass misses plausible impurity or substitution phases.

## XRD Calculator Validation

### Direct VESTA / RIETAN-FP display simulation

The desktop app now uses a two-stage calculation so library analysis remains responsive:

1. The full structure library is screened with the validated cached backend.
2. The selected target phase and up to three final displayed candidates are re-simulated through the local `VESTA.exe -> RIETAN-FP 3.12` route.

VESTA standardizes the CIF and writes the RIETAN input. Ophi switches that input to simulation mode, runs RIETAN in an OS temporary folder, and reads the generated continuous `.gpd` profile. The dashboard plots that profile directly; it does not broaden the RIETAN peak list a second time.

Open `VESTA / RIETAN` on the main window to inspect or change both executable paths. The compact RIETAN runtime is installed outside OneDrive at:

```text
%LOCALAPPDATA%\Ophiuchus\tools\rietan\RIETAN.exe
```

Temporary `.ins`, `.lst`, and `.gpd` files are removed after each calculation. The scientific report records the actual display engine. If direct simulation fails, Ophi reports the failure and keeps the screening result explicitly labeled as the fallback; it does not relabel pymatgen output as VESTA/RIETAN.

### Guarded multiphase Rietveld refinement

After an analysis has produced CIF-backed candidates, open `RIETAN 受约束精修` and switch from `单相确认` to `多相实验定量`. Select one target and one to three impurity phases. Ophi then runs three real RIETAN-FP jobs with equal, target-dominant, and impurity-dominant starting scales. It uses RIETAN's final scale factors and CIF-derived `Z * M * V` values for the Hill-Howard mass-fraction cross-check; candidate score or peak area is never presented as wt%.

Exact wt% is hidden when any hard gate fails, including invalid scale/ZMV values, failed independent pattern validation, an incomplete phase model indicated by grouped positive residual peaks, disagreement with RIETAN's native mass fractions, or unstable results across the three starting conditions. A non-failing result is still labeled `实验性定量` until the complete instrument/sample workflow is validated against measured standards of known composition. Amorphous content and phases absent from the selected model are not quantified.

The multiphase input is built from the official RIETAN `Cu3Fe4P6_combins` template and `cif2ins.exe`. Ophi does not redistribute that template. In `VESTA / RIETAN 设置`, click `安装多相支持` and select the official `Windows_versions.zip` when it is not found automatically. The installer validates the official marker structure and extracts only `template.ins` beside the configured RIETAN executable. Obtain RIETAN and its current citation instructions from the [official RIETAN-FP distribution page](https://jp-minerals.org/rietan/).

All RIETAN work files remain in an operating-system temporary directory. PNG, profile/reflection CSV files, and provenance JSON are written only after the user clicks the refinement window's export button. The JSON records engine paths, input/settings hashes, CIF SHA-256 values, gate findings, and repeat-run stability ranges.

Ophi can compare a CIF-derived simulated pattern against a VESTA/reference peak table exported by the user. This is the preferred trust check before relying on CIF-derived screening.

For your existing VESTA pattern, prefer the complete `.int` file rather than a selected top-peak CSV:

```powershell
python -m ophiuchus xrd-validate `
  --cif ".\data\library\structures\local\local_362b46a32755d09f.cif" `
  --reference "C:\path\to\MONI.int" `
  --out-dir ".\results\validation_reports" `
  --label "ZrFe6Ge4_vs_MONI_full"
```

Selected peak CSV files such as `all_peak_positions.csv` are still supported with `--phase`, but they are better for figure labels than full scientific validation because they may omit lower-ranked peaks.

By default, validation and CIF analysis try to use the base Anaconda `pymatgen` installation through a separate Python process. This avoids the current `ophi` environment hard-exit issue and prevents the simplified fallback simulator from over-amplifying high-angle peaks. If your Anaconda is installed somewhere else, set `OPHI_EXTERNAL_PYMATGEN_PYTHON` to that `python.exe`. Use `OPHI_DISABLE_EXTERNAL_PYMATGEN=1` only when intentionally testing fallback behavior.

Validation reports are written as JSON and Markdown. See `SCIENTIFIC_AUDIT.md` for current validation status and limitations. The current ZrFe6Ge4 vs VESTA `MONI.int` check passes peak-position and strong-peak coverage, but exact relative intensities remain model-dependent unless VESTA and pymatgen are compared with the same line-broadening/profile settings.

For transparent single-CIF debugging, export the Ophi simulated peak table:

```powershell
python -m ophiuchus xrd-debug-simulate `
  --cif ".\data\library\structures\local\local_362b46a32755d09f.cif" `
  --out ".\results\validation_reports\ZrFe6Ge4_debug_B4_peaks.csv" `
  --two-theta-min 10 `
  --two-theta-max 90 `
  --debye-waller-b 4
```

Then compare it with a VESTA/reference export:

```powershell
python -m ophiuchus xrd-debug-compare `
  --ophi ".\results\validation_reports\ZrFe6Ge4_debug_B4_peaks.csv" `
  --reference "C:\path\to\MONI.int" `
  --out ".\results\validation_reports\ZrFe6Ge4_debug_B4_vs_MONI.csv" `
  --position-tolerance 0.05
```

See `VESTA_VALIDATION.md` and `LIBRARY_TRACEABILITY.md` for the current trust-repair evidence.

## Dependency Notes

The core can parse experimental data, peak-list candidates, and many explicit CIF files with the Python standard library. Installing `pymatgen` in base Anaconda is strongly recommended for fuller crystallographic correctness, especially CIFs that depend on symmetry operations or complex occupancy conventions. Optional scientific plotting in later versions may use `matplotlib`, `scipy`, `sumo`, and `pyprocar`.
