# Stage 5 finding: HACApK-MNA (nested BiCGSTAB) does not cross over

Helical solenoid, 10 turns x 20 seg/turn (200 connected segments),
nwinc^2 parallel sub-filaments.  Freq 10 kHz.  Dense path is the
stock C++ `MNASolver` (LU).  HACApK path is
`PEECCircuitSolver(use_hacapk=True)` with outer rtol 1e-4, inner
rtol 1e-8, maxiter 5000 (the tightest tolerances that actually
converge on this geometry).

## Measurements (LAB, 2026-04-16)

| Mode   | N_fil | t_solve | peak mem | |Z11|    | rel-err vs dense |
|--------|------:|--------:|---------:|---------:|------------------:|
| dense  |   200 |  0.01 s |   116 MB | 0.5739   | --                |
| dense  |   800 |  0.12 s |   165 MB | 0.5682   | --                |
| dense  |  1800 |  1.07 s |   367 MB | 0.5660   | --                |
| dense  |  3200 |  5.59 s |   907 MB | 0.5652   | --                |
| dense  |  5000 | 20.36 s |  2034 MB | 0.5648   | --                |
| hacapk |   200 |  1.36 s |   109 MB | 0.5737   | 4.94e-04          |

At N=200 the HACApK driver is already **170x slower** than dense LU.
Beyond N=200 the outer BiCGSTAB stalls well above the loose 1e-4
target and/or wall time explodes (wall time ~300x dense at N=500,
previously measured at ~30 s).

## Why

The pure L matvec **does** scale (stage 3 bench shows 255x speedup
and 100x memory reduction at N=12500).  What does NOT scale is the
nodal BiCGSTAB on top:

1. Each outer matvec `Y @ v = A Z^{-1} A^T v` requires a full inner
   BiCGSTAB on the N-by-N branch system `Z I = A^T v`.
2. scipy BiCGSTAB on the inner problem breaks down around residual
   1e-7 on a strongly coupled helical coil (rho ~ 0).  This caps
   the usable outer tolerance at ~1e-4 for the bench to converge
   at all.
3. Even at outer=1e-4, the outer loop needs 700+ matvecs, each
   triggering a full inner solve with 20+ HACApK matvecs.  Product:
   ~60,000 real H-matrix matvecs per port impedance evaluation.

The pathology is visible even when the inner solve is replaced
with a dense LU factorization of Z: nested Krylov on Y takes 2700
outer matvecs at outer=1e-6, whereas direct BiCGSTAB on the
assembled Y matrix converges in ~100 iterations.  Inner noise,
however small, destroys outer Krylov convergence efficiency.

## Implication for the paper plan (memory/hacapk_peec_prima_plan.md)

The "HACApK-MNA crossover N" measurement we expected to land
somewhere above 2000 does not exist with the current architecture.
Reporting this as a negative result is valuable: it motivates the
stage-6+ redesign rather than hiding the limitation.

## Suggested next steps (stage 6+)

1. **FGMRES outer with BiCGSTAB inner.**  scipy.sparse.linalg.lgmres
   is flexible with respect to preconditioner changes and should
   tolerate the inner residual jitter that kills BiCGSTAB-on-BiCGSTAB.
2. **Saddle-point formulation.**  Solve the coupled 2-by-2 block
   system `[Z_branch, A^T; A, 0] [I; V] = [0; i_ext]` directly with
   a single Krylov method + block preconditioner (e.g. Schur
   complement approximation using `diag(A Z_diag^{-1} A^T)`).
3. **Assemble Y as its own H-matrix.**  Since A is sparse and Z has
   an H-matrix, `Y = A Z^{-1} A^T` can be built in H-matrix form
   (inverse H-matrix of Z composed with sparse A).  Solve Y V = b
   with direct H-matrix back-substitution, no nested Krylov.

The current commit keeps the nested driver in the library as-is,
so the stage-6 successor can diff the finding numerically.
