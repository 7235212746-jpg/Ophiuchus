# Ophiuchus Version History

## v0.5.0-dev - 2026-07-22

- Added guarded two-to-four-phase RIETAN-FP refinement with explicit target/impurity roles. (v0.5.0-dev)
- Added official `cif2ins + combins` multiphase input construction and validated template installation from RIETAN's official archive. (v0.5.0-dev)
- Added Hill-Howard/ZMV mass-fraction cross-checking, three-start stability tests, and a hard gate that hides untrustworthy wt%. (v0.5.0-dev)
- Added direct RIETAN total-profile, background, residual, and per-phase reflection rendering with explicit-only export. (v0.5.0-dev)
- Added opt-in real-engine tests for the official three-phase example and RIETAN 3.12 high-angle references for Zr3V3GeSn4 and ZrFe6Ge4. (v0.5.0-dev)

## v0.4.0-dev - 2026-07-01

- Added Phase 4 multi-phase candidate screening traceability. (v0.4.0-dev)
- Added `CandidateStructureService` for unified library-to-analysis candidate visibility, cache status, enabled state, used/skipped state, and skip reasons. (v0.4.0-dev)
- Added `candidate_usage_summary.csv` export and `candidate_usage` JSON evidence for library XRD analysis. (v0.4.0-dev)
- Added conservative score components and confidence labels to candidate ranking outputs. (v0.4.0-dev)
- Enhanced Peak Inspector evidence with structure hash, CIF path, enabled state, and strong-peak context. (v0.4.0-dev)
- Updated the desktop structure library table to show candidate-ready status and skip reasons. (v0.4.0-dev)
- Added explicit presentation and diagnostic XRD plot outputs, including audit-folder copies and deterministic export bundle names. (v0.4.0-dev)
- Added `CANDIDATE_MATCHING.md` and `LIBRARY_USAGE.md` for screening limits and library traceability. (v0.4.0-dev)

## v0.1.0 - 2026-06-27

- Created the first Ophiuchus desktop research workflow project. (v0.1.0)
- Added the Ophi XRD Candidate Screener vertical slice for local XRD import, peak extraction, candidate scoring, multi-phase heuristic explanation, and reports. (v0.1.0)
- Added local CIF peak simulation through pymatgen when available and a lightweight fallback simulator for explicit P1/VESTA-style CIFs when pymatgen is unavailable. (v0.1.0)
- Added a simple local desktop window and CLI that share the same analysis pipeline. (v0.1.0)
