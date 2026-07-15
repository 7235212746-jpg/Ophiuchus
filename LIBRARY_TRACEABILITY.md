# Library Traceability

Date: 2026-06-30

## Unified Flow

Ophi now treats the local structure library as the traceable structure source for library analysis:

```text
Manual CIF import
Materials Project harvest
Existing local CIF import
        -> StructureLibrary SQLite index
        -> copied/cached CIF files
        -> simulated XRD cache
        -> enabled candidate screening
        -> report source summary and per-candidate provenance
```

## Required Stored Evidence

Each structure entry stores:

- internal structure id
- source and source_id
- formula/reduced formula
- element list and chemical system
- original CIF/source path
- copied local CIF path
- metadata JSON path when available
- enabled/disabled state
- access/license note

Each simulated XRD cache row stores:

- structure id
- 2theta
- relative intensity
- radiation and settings hash
- hkl when available
- d-spacing when available

## UI Visibility

The Structure Library tab now shows source id in the table and displays full traceability details in the selected-entry panel:

- source/source_id
- enabled state
- XRD cache status
- CIF path
- metadata path
- original path
- access note

The tab also includes actions:

- refresh library
- enable/disable selected structure
- open CIF location
- view raw CIF text
- view simulated peak table
- recompute selected XRD

## Analysis Traceability

Library analysis reports now include `candidate_source_summary` in `results.json` and the Markdown report. The summary records, by source:

- loaded
- enabled
- simulated
- used
- skipped because disabled
- skipped by element filter
- skipped because no valid XRD cache

Each top candidate already includes candidate id, source, and source path. The candidate CSV now includes `candidate_id`.

## Verified State

The current default library contains:

- local structures: 3
- Materials Project structures: 21
- total structures: 24

A real Materials Project import test retrieved 21 structures and imported 21 with 0 failures. XRD cache generation simulated 21 new patterns with 0 failures.
