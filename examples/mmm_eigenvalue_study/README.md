# MMM Eigenvalue Study

Spectral analysis of the MMM (Magnetic Moment Method) interaction matrix `N`
and the deflated BiCGSTAB convergence behavior on the C-type electromagnet
geometry (mu_r = 100000, 4x4x2 and 6x6x1 element grids).

## Scripts

| Script | Purpose |
|---|---|
| `beautiful_ugly_viz.py` | Per-element spectrum visualization on coarse vs fine meshes (the "beautiful" well-conditioned cases vs "ugly" near-zero eigenvalues that break iterative solvers). |
| `belt_loop_validation.py` | Verify the dominant-eigenvalue belt is consistent across mesh refinement. |
| `block_spectrum.py` | Per-element 6x6 block diagonal eigenvalue extraction (MSC hex DOF). |
| `cost_deflation_benchmark.py` | Compare iteration count + wall-clock cost of plain BiCGSTAB vs deflated BiCGSTAB across mu_r ranges. |
| `deflated_bicgstab_real.py` | Real-arithmetic deflated BiCGSTAB implementation using the dominant left/right eigenvectors as the deflation basis. |
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
