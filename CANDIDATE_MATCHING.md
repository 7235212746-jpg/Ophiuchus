# Candidate Matching

Date: 2026-07-01

Ophi Phase 4 performs XRD candidate screening. It does not perform Rietveld refinement, phase-fraction estimation, or definitive phase confirmation.

## What The Score Means

Each enabled candidate phase with a valid simulated XRD cache is compared against detected experimental peaks within the configured 2theta tolerance.

The score is decomposed into:

- matched theoretical peak count
- strong theoretical peak count
- matched strong-intensity fraction
- missing strong-intensity fraction
- missing strong-peak penalty
- multi-peak bonus
- final score
- conservative confidence label

The confidence labels are screening labels only:

- `high`: many strong peaks match and very few strong theoretical peaks are missing
- `medium`: several peaks match, but evidence is not definitive
- `low`: some support exists, but the phase needs manual checking
- `uncertain`: only weak or isolated support exists
- `unlikely`: no useful peak support was found

Ophi intentionally avoids words such as confirmed unless a user manually records that conclusion outside the automated score.

## Multi-Phase Explanation

The multi-phase mode is a greedy screening heuristic:

1. Rank all enabled cached candidate phases.
2. Choose the explicit target phase when provided and sufficiently supported.
3. Treat other candidates as possible impurities if they explain residual experimental peaks.
4. Penalize candidates with many strong theoretical peaks absent from the experiment.
5. Report remaining unresolved experimental peaks.

This is useful for prioritizing which phases to inspect next. It is not quantitative refinement.

## Peak Assignments

Peak assignments are deterministic. For each experimental peak, Ophi records nearby simulated peaks from ranked candidates, including:

- candidate structure ID
- formula
- source and source path
- calculated 2theta
- delta 2theta
- hkl when available
- simulated relative intensity
- screening confidence

Overlapping assignments are expected in multi-phase materials. Ophi reports overlaps instead of forcing a single answer.

## Scientific Limits

XRD peak screening can miss preferred orientation, texture, strain, profile broadening, background errors, amorphous phases, and instrumental effects. Candidate scores should be used as a shortlist for manual crystallographic review, VESTA/reference comparison, and later refinement where appropriate.
