# MMM Eigenvalue Study

Spectral analysis of the MMM (Magnetic Moment Method) interaction matrix `N`
and the deflated BiCGSTAB convergence behavior on the C-type electromagnet
geometry (mu_r = 100000, 4x4x2 and 6x6x1 element grids).

## Scripts

| Script | Purpose |
|---|---|
| `beautiful_ugly_viz.py` | Per-element spectrum visualization on coarse vs fine meshes (the "beautiful" well-conditioned cases vs "ugly" near-zero eigenvalues that break iterative solvers). |
| `block_spectrum.py` | Per-element 6x6 block diagonal eigenvalue extraction (MSC hex DOF). |
| `deflated_bicgstab_real.py` | Real-arithmetic deflated BiCGSTAB implementation using the dominant left/right eigenvectors as the deflation basis (pure-numpy study; the C++ runtime deflation API was removed 2026-06-09 -- see note below). |
| `hmatrix_truncation_spectrum.py` | How ACA tolerance truncation shifts the H-matrix's effective spectrum. |

## Output

Generated PNGs sit next to their producing script:
`beautiful_ugly_4x4x2_mu100000.png`, `beautiful_ugly_6x6x1_mu100000.png`,
`bicgstab_dynamics_4x4x2_mu100000.png`, etc.

## Context

This study informs the MSC nullspace-deflation strategy documented in
`docs/solver/MSC_NULLSPACE_DEFLATION.md` and the broader Phase B
equivalence_source C++ work.  Not a production sample -- research
exploration that feeds into solver design decisions.

**Outcome (2026-06-09):** the design decision this study fed into is made --
the project consolidated on the **HDiv-VIM** operator (`radia.vim`), where
the loop space is `ker(B)` (field-null by construction via de Rham), so NO
runtime loop deflation / loop-star gauge / loop projection is needed. The C++
solver APIs (`SetHACApKDeflation`, `SetDeflateNullspace`, `SolveLoopStar`,
`SetLoopProjection`) were removed; the scripts here that exercised them were
deleted. The remaining scripts are pure-numpy spectral analysis (no removed
API) that reproduce the nullspace / conditioning theory motivating the choice.
