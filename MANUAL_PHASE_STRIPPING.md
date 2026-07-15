# Manual Phase Stripping And Residual Spectrum

Ophiuchus provides a manual, non-destructive workspace for testing whether a
candidate phase can explain part of an experimental XRD pattern. It is a
screening and evidence tool. It is not Rietveld refinement, and the fitted
scale is not a quantitative phase fraction or weight percent.

## Opening The Workspace

1. Import the experimental XRD file and select the complete element system.
2. Select the target phase explicitly, then run local-candidate or structure-library analysis.
3. When analysis succeeds, click `手动物相剥离 / 残差谱分析` in the main window.

The button remains disabled until the main analysis has supplied the original,
unnormalized XRD grid and canonical candidate patterns. Starting another main
analysis is blocked while the stripping window is open so that the two windows
cannot silently refer to different runs.

## Three-Pane Workflow

- **Candidate pane:** search and sort phases, inspect evidence, preview a fit,
  accept it, exclude it from ranking, or open its CIF.
- **Plot pane:** compare the immutable original pattern, fitted total, signed
  residual, accepted phase contributions, and the current preview. Negative
  residual intensity remains visible.
- **Control pane:** choose full-spectrum or selected-range fitting, inspect and
  adjust bounded scale/shift/width parameters, then use Preview, Accept,
  Cancel, Undo, Redo, Reset, and Export.

Preview never changes the accepted session. Accept records one immutable
operation. Undo, Redo, and Reset recompute from the original intensity array;
they do not repeatedly subtract and add rounded curves.

## Scientific Model

For accepted candidate phases `k`, Ophiuchus uses

```text
c_k(x) = s_k P_k(x - delta_k, sigma_k)
y_residual(x) = y_original(x) - sum_k c_k(x)
```

`P_k` is projected from the complete canonical peak list produced by the same
validated XRD backend used by the main analysis. One global 2theta shift is
applied to the whole phase; individual peaks cannot move independently.
Scale is constrained to be non-negative, and shift and peak width are bounded.
The original intensity array is loaded without normalization, frozen, and never
clipped or overwritten.

The candidate list is re-ranked against the current residual after every
accepted or undone operation. Ranking rewards multi-peak residual coverage and
profile improvement, and penalizes missing strong theoretical peaks,
one-peak coincidences, chemically implausible candidates, and excessive
negative residual area. Duplicate structures/patterns are collapsed while
their provenance IDs remain available.

## Warnings And Interpretation

Inspect the fit before accepting when Ophiuchus reports any of these conditions:

- fitted shift, width, or scale reaches a configured bound;
- subtraction increases negative residual area;
- a candidate is supported by too few independent peaks;
- strong theoretical peaks are absent from the residual;
- the candidate lies outside the selected element scope.

A useful subtraction means that the candidate can explain coherent features
under the chosen profile assumptions. It does not prove uniqueness, determine
a crystal-structure refinement, or produce a defensible phase fraction.

## Saving And Exporting

Normal GUI analysis writes one transactional current session under
`%LOCALAPPDATA%\Ophiuchus\analysis-session`. A successful new analysis replaces
the previous temporary session. Existing project `results` files are never
deleted or expanded automatically.

- Click `保存本次分析` in the main window to copy the complete current analysis
  to a permanent folder that you choose.
- Click `Export` in the stripping window to create an explicit bundle containing
  `phase_stripping_residual.csv`, `phase_stripping_session.json`, and
  `phase_stripping_residual.png`.

The CSV contains 2theta, original intensity, fitted total, signed residual, and
one column per accepted contribution. The JSON records operation order,
candidate/CIF provenance and hashes, instrument settings, fit bounds, and fit
parameters. Export reuses stored contributions and does not rerun simulation.

The command-line analysis interface keeps its existing behavior and writes
directly to the specified `--out-dir`.

## VESTA Boundary

Ophiuchus can validate canonical peak tables against user-exported VESTA or
other reference patterns. VESTA's Powder Diffraction command prepares RIETAN-FP
input and launches an external RIETAN-FP executable; it is not a standalone
documented VESTA command-line XRD engine. The local VESTA folder currently does
not contain `rietan.exe` or `rietan64.exe`, so Ophiuchus does not silently claim
to use VESTA for simulation. See `VESTA_VALIDATION.md` for the existing
reference comparison and its limits.
