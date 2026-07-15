# Ophiuchus Implementation Status

This document tracks what Ophiuchus can currently do and how the complete workflow is implemented. Each implementation item ends with the version where it was added or last materially changed.

## Current Capability

Ophiuchus is a local-first materials research workflow project. In v0.1.0, its working vertical slice is the Ophi XRD Candidate Screener: it imports experimental XRD files, extracts peaks, scans local candidate structures or peak-list files, simulates CIF peak positions with pymatgen when available or a lightweight fallback simulator when not, ranks candidate phases, produces a small multi-phase heuristic explanation, and writes transparent reports. (v0.1.0)

## Complete Workflow

1. The user starts the desktop shell with `start_ophiuchus.bat` or `python -m ophiuchus app`. (v0.1.0)
2. The user selects an experimental XRD file and a local candidate directory containing CIF files or peak-list files. (v0.1.0)
3. Ophiuchus imports Rigaku `.asc` files or common two-column text files, normalizes intensity to max 100, and extracts experimental peaks. (v0.1.0)
4. Ophiuchus scans local candidates, filters by allowed elements, reads sparse peak lists directly, and peak-picks dense simulated patterns such as `MONI.int`. (v0.1.0)
5. For CIF candidates, Ophiuchus uses pymatgen XRD simulation when installed; if pymatgen is unavailable, it uses a lightweight fallback simulator for explicit atom-site CIFs. (v0.1.0)
6. Candidate scoring rewards matched theory peaks and penalizes missing strong theory peaks, so isolated accidental matches do not receive unrealistic confidence. (v0.1.0)
7. The report includes top candidates, missing strong peaks, a heuristic multi-phase explanation, unexplained experimental peaks, JSON output, Markdown output, and CSV tables. (v0.1.0)

## Versioning Rule

- Larger feature phases update the middle version number, for example `v0.1.0` to `v0.2.0`. (v0.1.0)
- Smaller fixes update the patch version, for example `v0.1.0` to `v0.1.1`. (v0.1.0)
- Only `0.x` milestone versions are archived in `Ophiuchus_Versions`; small patch updates may overwrite the active project folder. (v0.1.0)
- Every release must update `docs/version_history/CHANGELOG.md`. (v0.1.0)
- Every feature update must update this implementation status document, and each implementation item must end with the version where it was added or materially changed. (v0.1.0)
