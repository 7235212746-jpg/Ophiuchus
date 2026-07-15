# Ophi Direct XRD Backend Replacement

## Scope

Ophi now uses one canonical CIF-to-XRD implementation for local candidates,
structure-library analysis, cache generation, plotting, and report export.
Reference patterns exported by VESTA remain validation evidence only. They are
never substituted for a CIF-derived result.

## Canonical Data Contract

`ophiuchus.xrd.backend.SimulatedPattern` is the single scientific pattern model.
It records:

- source structure ID, CIF path, and full SHA-256;
- formula and space-group number;
- radiation, wavelength, scan range, and line model fingerprint;
- backend, backend version, and pymatgen dependency version;
- full sorted 2theta, d-spacing, hkl, multiplicity, line-component,
  raw-intensity, and normalized-intensity arrays;
- a deterministic pattern fingerprint.

Near-degenerate reflections from the same line component are merged once using
the configured angular tolerance and intensity-weighted position; their raw
intensities and multiplicities are summed before normalization. No relative-
intensity threshold is then applied. Weak high-angle reflections remain in the
canonical array and can be inspected or replotted.

## Simulation Rules

- Structure factors are calculated by pymatgen's `XRDCalculator` with
  `scaled=False`.
- Cu Kalpha1 is calculated at 1.54056 angstrom.
- In the Cu Kalpha1+2 model, Kalpha2 is generated from each d-spacing at
  1.54439 angstrom with a 0.5 intensity ratio.
- Intensities are normalized only after all requested line components are
  present.
- A malformed or unreadable CIF stops that candidate with an actionable error.
  Ophi does not fall back to an unvalidated hand-written simulator.

## Cache Safety

The validated cache lives in `xrd_patterns_v2`. Its key includes the structure
ID, full CIF SHA-256, and the complete simulation-settings fingerprint. The
payload stores the full canonical pattern without rounding or peak filtering.

Scientific Safe Mode is enabled by default in the desktop app. It recomputes
each in-scope CIF during the analysis. When Safe Mode is disabled, only an exact
validated v2 cache match can be reused. Legacy peak-only cache rows cannot
replace CIF-derived patterns in an analysis.

## Phase Candidates

Library records are converted to phase-level candidates before scoring.
Identical structure hashes are grouped. Records with the same reduced formula,
space group, and major-peak fingerprint may also be grouped. Polymorphs with
different space groups remain separate candidates.

An explicitly selected target phase is always kept as the target, even when its
score is low. Impurity ranking then favors distinct non-target phases that
explain experimental residual peaks. Each grouped candidate retains all source
record IDs for provenance.

## Plot and Report Contract

The Matplotlib dashboard uses the same canonical arrays used for matching.
Target peak positions and intensities are not recalculated and are not loaded
from a VESTA reference file. The target is shown with its simulated profile and
peak labels. Up to three phase-level impurities are shown as separate peak rows;
their strongest peaks receive dashed guides without numeric peak labels.

`results.json` stores the full canonical pattern for every scored candidate.
`candidates.csv` includes space group, equivalent library record IDs, backend
version, and pattern fingerprint so a result can be audited and replotted.

## Validation Performed

- Real `Zr3V3GeSn4.cif`, 10-90 degrees, Cu Kalpha1+2: 146 merged canonical reflections,
  including reflections below 1 percent relative intensity.
- Real `Ti3V3GeSn4.cif` from the previously failing P1 workflow: near-degenerate
  reflections are merged from 919 raw line records to 151 canonical Kalpha1+2
  records over the experimental range. All eight curated VESTA major peaks are
  recovered within 0.0062 degrees; the maximum major-peak intensity ratio error
  is 1.2269.
- The real Ti structure-library run analyzed 82 records with zero failures in
  18.8 seconds in Safe Mode; the validated-cache rerun hit all 82 records in
  about 9 seconds.
- Direct backend, local candidate route, structure-library route, JSON report,
  and Matplotlib target stems are asserted to use identical 2theta arrays.
- A fake matching VESTA reference is asserted not to alter simulation, caching,
  or plotting.
- Curated top-N VESTA peak tables are treated as major-peak coverage evidence,
  not as exhaustive proof that every unlisted calculated reflection is wrong.
- Duplicate library records are asserted to group while different-space-group
  polymorphs remain separate.

## Scientific Boundary

This replacement removes inconsistent software paths; it does not turn candidate
screening into definitive phase identification or quantitative phase analysis.
Preferred orientation, absorption, crystallite size, strain, occupancy errors,
instrument broadening, and sample displacement can still change experimental
intensities. VESTA comparisons should therefore be treated as independent
validation evidence, not as a hidden answer key.
