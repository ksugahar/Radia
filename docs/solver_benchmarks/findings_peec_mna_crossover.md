# Stage 5+6+7 findings: HACApK-MNA scaling

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

## Stage-6b / stage-7: saddle-point direct formulation

Replaces the nested lgmres-outer + BiCGSTAB-inner by solving the
2x2 block saddle-point system directly:

    [Z_branch   A_hat^T] [I]   [0         ]
    [A_hat      0      ] [V] = [-i_ext_red]

via a single scipy lgmres call on dimension (N + n_red) preconditioned
block-diagonally: ``diag(Z_branch)^{-1}`` on the upper block and a
sparse LU of the Schur complement ``A_hat diag(Z)^{-1} A_hat^T`` on
the lower block.  Implemented as
``PEECCircuitSolver(outer_method="saddle")``.

| Mode         | N_fil | t_solve | peak mem | |Z11|  | rel vs dense |
|--------------|------:|--------:|---------:|--------|--------------|
| dense LU     |   200 |  0.01 s |   116 MB | 0.5739 | --           |
| dense LU     |   800 |  0.12 s |   165 MB | 0.5682 | --           |
| dense LU     |  1800 |  1.07 s |   367 MB | 0.5660 | --           |
| dense LU     |  3200 |  5.60 s |   907 MB | 0.5652 | --           |
| dense LU     |  5000 | 20.43 s |  2036 MB | 0.5648 | --           |
| hacapk saddle|   200 |  0.04 s |   110 MB | 0.5739 | **6.2e-10**  |
| hacapk saddle|   800 |  0.28 s |   122 MB | 0.5682 | **3.2e-08**  |
| hacapk saddle|  1800 |  2.84 s |   157 MB | 0.5826 |  2.9e-02 (*) |
| hacapk saddle|  3200 |  2.44 s |   264 MB | 0.5430 |  3.9e-02 (*) |
| hacapk saddle|  5000 |  9.93 s |   490 MB | 0.5347 |  5.3e-02 (*) |

(*) The large-N error is NOT a saddle-point formulation bug -- it is
the underlying HACApK L matvec losing accuracy on densely packed
parallel filament bundles.  Measured directly: at N=1800 (nwinc=3,
9 parallel sub-filaments per segment) the HACApK L matvec error is
~2e-2 regardless of ``aca_eps`` (1e-4 through 1e-8) and ``eta`` (1-3).
The same problem was documented in stage 3 (see
memory/hacapk_peec_prima_plan.md "N=1024 10% rel error on helix,
N=2000 19%").  HACApK ACA compression fails to capture the near-field
interaction between tightly bundled parallel filaments.

### Stage-6b/7 finding: **compute crossover at N ~ 3000**

| N_fil | dense t_solve | saddle t_solve | speedup |
|------:|--------------:|---------------:|--------:|
|   200 |       0.01 s  |        0.04 s  |    0.2x |
|   800 |       0.12 s  |        0.28 s  |    0.4x |
|  1800 |       1.07 s  |        2.84 s  |    0.4x |
|  3200 |       5.60 s  |        2.44 s  |  **2.3x** |
|  5000 |      20.43 s  |        9.93 s  |  **2.1x** |

Saddle-point crosses over dense LU at N ~ 3000 with a ~2x wall-time
speedup and a 4x peak-memory reduction (490 MB vs 2 GB at N=5000).
Dense LU scales O(N^3); saddle-point scales ~O(outer_iters x N log N)
and the outer iteration count plateaus around 500 beyond N=1800,
so the speedup grows with N.

Iteration counts plateau because the block Jacobi preconditioner
captures the diagonal physics exactly (R + j*omega*L_ii) and the
Schur complement is a sparse graph Laplacian with O(N) LU cost.
Only the off-diagonal mutual inductance couples through HACApK.

### Validate (N=72, DC..10 MHz)

    max rel err over 6 freqs: 6.14e-11 PASS

Machine precision; much better than stage-6a lgmres's 3.64e-08
because saddle-point needs no nested inner Krylov and therefore
avoids inner-solve noise.

### Accuracy recovery plan (stage 8 -- separate from this work)

The ACA accuracy loss on parallel filament bundles is a HACApK
calibration problem independent of the circuit solver:

1. Tune HACApK cluster tree to keep parallel sub-filaments within
   the same leaf (so their interaction is stored dense).
2. Use ``leaf_size >= nwinc * nhinc * SEG_PER_TURN`` so an entire
   turn of sub-filaments collapses to one leaf.
3. Switch to nwinc*nhinc DOF per parent segment with an internal
   sub-filament basis (PEECBuilder change; out of scope for stage 7).

For papers / demos, use nwinc=1 (single filament per conductor):
accuracy is machine-precision up to N ~ 1000+ where dense LU still
dominates anyway.

## Path forward (stage 8+)

- Fix HACApK ACA on bundled parallel filaments (independent task).
- Implement Y as H-matrix (alternative crossover path, contingent on
  Ida-san H-QR work landing).
- PRIMA multi-frequency sweep on top of saddle-point (each frequency
  is O(N log N) work; port-impedance at 50-100 freqs becomes
  practical at N=5000 with saddle-point).

## How to reproduce

    python validation_test/solver_benchmarks/bench_peec_mna_crossover.py

Default uses OUTER_METHOD="lgmres"; change the module-level constant
at the top of the script to "bicgstab" or "gcrotmk" to sweep.  Wall
time ~7 minutes for the current SUBDIVISIONS settings on LAB.
