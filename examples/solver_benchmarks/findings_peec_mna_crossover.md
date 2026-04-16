# Stage 5+6 findings: HACApK-MNA scaling

Helical solenoid, 10 turns x 20 seg/turn (200 connected segments),
nwinc^2 parallel sub-filaments.  Freq 10 kHz.  Dense path is the
stock C++ `MNASolver` (LU).  HACApK path is
`PEECCircuitSolver(use_hacapk=True)` with outer rtol 1e-4, inner
rtol 1e-8, maxiter 5000 (the tightest tolerances that actually
converge on this geometry).

## Stage-5 baseline (BiCGSTAB outer)

| Mode           | N_fil | t_solve |
|----------------|------:|--------:|
| dense LU       |   200 |   0.01 s |
| hacapk bicgstab|   200 |  ~1.4 s  |
| hacapk bicgstab|   500 |  ~30 s (diverges past rtol 1e-4) |

BiCGSTAB outer stalls above N ~ 200 on a strongly coupled helical
coil because the inner HACApK BiCGSTAB breaks down at residual
~1e-7 and BiCGSTAB's three-term recurrence cannot tolerate that
varying-operator noise.  Stage 5 concluded with this negative
result and scheduled three remediation options.

## Stage-6a: lgmres outer (commit `<stage-6a>`)

Option 1 of the stage-5 plan.  scipy.sparse.linalg.lgmres replaces
the outer BiCGSTAB and an exact diagonal Jacobi preconditioner is
added on Y = A Z^{-1} A^T.  Both are wired through the
`outer_method` keyword on `PEECCircuitSolver`; gcrotmk is also
exposed but NaN-divides at N=72 + 10 kHz so lgmres is default.

| Mode         | N_fil | t_solve | peak mem | |Z11|  | rel vs dense |
|--------------|------:|--------:|---------:|--------|--------------|
| dense LU     |   200 |  0.01 s |   115 MB | 0.5739 | --           |
| dense LU     |   800 |  0.12 s |   165 MB | 0.5682 | --           |
| dense LU     |  1800 |  1.07 s |   367 MB | 0.5660 | --           |
| dense LU     |  3200 |  5.59 s |   907 MB | 0.5652 | --           |
| dense LU     |  5000 | 20.44 s |  2034 MB | 0.5648 | --           |
| hacapk lgmres|   200 |  3.56 s |   111 MB | 0.5742 | 4.07e-04     |
| hacapk lgmres|   800 | 35.01 s |   121 MB | 0.5681 | 7.04e-05     |
| hacapk lgmres|  1800 | 402.19 s|   157 MB | 0.5821 | 2.85e-02     |

(At N=1800 lgmres converges only to residual ~1e-3 under the 10x
slack on outer_tol=1e-4, explaining the 3% Z_11 error.  Tightening
outer_tol costs 4-10x more iterations and was not measured.)

### Stage-6a finding: robust, not faster

lgmres converges everywhere BiCGSTAB did not, which is a prerequisite
for the paper story.  But **it still does not beat dense LU in wall
time**:

  N=200  :  350x slower
  N=800  :  290x slower
  N=1800 :  375x slower

**The memory story is better.**  At N=1800 HACApK-MNA uses 157 MB
(43% of dense's 367 MB); at N=5000 (where HACApK was not measured
but L matvec alone stays ~150 MB per stage-3 bench) the ratio drops
below 10%.  The memory crossover is real, the compute crossover is
not.

### Why compute still loses

Per outer matvec the code does a full inner HACApK BiCGSTAB on
Z_branch I = A^T v.  That inner solve costs ~30-50 real HACApK
matvecs, each O(N log N), so the total work per port impedance is

    O( outer_iters × inner_iters × N log N )
      ~ 1800 × 50 × N log N           at N=1800

versus dense LU's O(N^3).  The constants plus the inner Krylov
overhead keep the nested approach from amortizing even at N ~ 5000.

## Path forward (stage 7+)

Since option 1 gives robustness but not compute speedup, we need one
of the heavier options:

1. **Saddle-point direct formulation.**  Solve the 2x2 block system
       [Z_branch   A^T] [I]   [0    ]
       [A          0 ] [V] = [i_ext]
   with a single Krylov + block preconditioner (Schur complement
   approximation using diag(Z_branch)).  Collapses the nested
   structure into a single Krylov iteration.

2. **Y as H-matrix.**  Build the inverse H-matrix of Z_branch (HACApK
   QR or H-LU on the real-part L-only side), compose with sparse A
   to form Y as an H-matrix, and solve Y V = b by direct back-sub.
   Zero nested Krylov.  Highest ROI if we eventually add H-matrix
   QR anyway for frequency sweep speed.

Stage 7 should start with (1) because it reuses existing HACApK
matvec + scipy Krylov plumbing; (2) becomes attractive if and when
Ida-san's H-QR work lands (see memory/hacapk_peec_prima_plan.md).

## How to reproduce

    python examples/solver_benchmarks/bench_peec_mna_crossover.py

Default uses OUTER_METHOD="lgmres"; change the module-level constant
at the top of the script to "bicgstab" or "gcrotmk" to sweep.  Wall
time ~7 minutes for the current SUBDIVISIONS settings on LAB.
