# Library Usage

Date: 2026-07-01

Ophi Phase 4 uses one unified local structure library path for XRD candidate screening.

```text
Local CIF import / Materials Project harvest
        -> StructureLibrary SQLite
        -> cached CIF and metadata files
        -> library XRD cache
        -> enabled candidate structures
        -> multi-phase XRD screening
        -> usage summary, ranking, peak assignments, and audit export
```

## Candidate Phase Panel

In the desktop app, open the `结构库` tab to see candidate structures. Each row shows:

- formula and chemical system
- source and source ID
- enabled state
- XRD cache status
- candidate-ready state

Selecting a row shows:

- structure ID
- structure hash
- CIF path
- cache status
- last simulated time
- skip reasons, if any

Use `启用/禁用选中` before analysis to decide which structures may participate.

## How A Structure Becomes A Candidate

A structure is candidate-ready when all of these are true:

- it is enabled by the user
- its elements are inside the selected analysis element scope
- its CIF/cache path exists
- it has a valid simulated XRD cache for the current radiation and 2theta range

Structures are skipped with explicit reasons such as:

- disabled by user
- outside selected element scope
- outside candidate scope
- no valid CIF path
- simulation cache missing

## Analysis Outputs

Every library analysis writes:

- `candidate_usage_summary.csv`: all structures, used/skipped state, and skip reasons
- `top_candidates.csv`: ranked candidates and score components
- `peak_assignments.csv`: experimental peak assignments
- `results.json`: complete machine-readable evidence
- `audit_*`: detailed tables for loaded structures, simulated peaks, peak matches, and unresolved peaks

The most important file for traceability is `candidate_usage_summary.csv`. It lets the user verify:

1. which CIF/structures were in the local library
2. which were candidate-ready
3. which were actually used in matching
4. which were skipped and why

## Peak Inspector

The `Peak Inspector` tab reads the same local simulated XRD cache. Enter an experimental 2theta value to list nearby simulated peaks from enabled structures. The result includes structure ID, formula, source, CIF path, calculated 2theta, delta, hkl, intensity, and strong-peak context.
