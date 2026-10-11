# Hiruma SparseSolv benchmark

This directory is the durable validation lane for the 30 kHz, three-material
Hiruma eddy-current problem used to measure Compact AMS + COCR.  It replaces
the retired `src/ext/sparsesolv/examples/hiruma` location.

The repository tracks the smallest case, `meshes/mesh1_2.5T.vol`, as the
repeatable validation fixture.  The larger 3.5T, 4.5T, 5.5T and 20.5T meshes
are optional scaling inputs because together they add about 78 MB.  Place them
in `meshes/` and use `--all` when a full performance sweep is required.

Run the numerical gate from the repository root:

```powershell
python validation_test/sparsesolv/hiruma/bench_compact_ams.py --verify-baseline
```

The gate fixes mesh hashes and finite-element sizes exactly, while accepting a
small iteration/residual range.  Setup and solve times are recorded but are not
gated because they depend on the machine and concurrent load.  Heavy scaling
runs belong on mdx; the tracked 2.5T case is suitable as a validation host smoke run.

The fixture was converted from the existing laboratory `mesh1_2.5T.msh`
dataset.  `prepare_hiruma_vol.py` preserves its points and volume elements and
adds explicit air/material interface triangles.  The resulting `.vol` passes:

```powershell
check-vol validation_test/sparsesolv/hiruma/meshes/mesh1_2.5T.vol `
  --contract validation_test/sparsesolv/hiruma/hiruma_vol_label_contract.json
```

The assembled complex system before and after adding the interface descriptors
is identical to floating-point precision.  The source and final SHA-256 values,
mesh sizes, golden range, and one observed timing are recorded in
`compact_ams_baseline.json`; the complete mesh check is stored in
`mesh1_2.5T.vol-check.json`.

Additional comparisons remain available as `bench_ams_vs_abmc.py` and
`bench_cocr_vs_gmres.py`.  They write transient result JSON below
`private-runtime-path`, never into the repository.

## Order 2: AMS as the BDDC wirebasket solver

`bench_bddc_coarse.py` solves the same problem on `HCurl(order=2,
nograds=True)` with NGSolve BDDC + COCR (tol 1e-8).  The system carries an
`eps*nu` mass with eps = 1e-6: without it the sigma=0 air leaves the
lowest-order gradients in the wirebasket kernel and the direct wirebasket
factorization breaks down (COCR stops at iteration 0).  At order 2 the
wirebasket is the lowest-order edge block, so importing the sparsesolv module
lets BDDC use `coarsetype="sparsesolv_ams"` instead of its direct inverse.
HCurl also puts the face dofs of badly shaped faces into the wirebasket
(124,040 of them on the 2.5T mesh); `--edge-wirebasket` returns them to the
interface.  The larger meshes are the optional scaling inputs described above.

worker-a, 8 threads, 2026-09-27, one binary (`bddc_coarse/`):

| mesh | dofs | coarse solver | iterations | setup s | solve s | wall s | peak GB |
|---|---|---|---|---|---|---|---|
| 2.5T | 680,041 | direct (sparse Cholesky) | 217 | 71.2 | 101.5 | 174.3 | 11.6 |
| 2.5T | 680,041 | AMS, 4 cycles, edge wirebasket | 594 | 2.1 | 70.5 | 77.3 | 1.28 |
| 3.5T | 864,781 | direct | 231 | 109.8 | 156.5 | 268.4 | 15.8 |
| 3.5T | 864,781 | AMS, 4 cycles, edge wirebasket | 783 | 2.5 | 143.1 | 150.7 | 1.58 |
| 5.5T | 1,456,057 | direct | fails (process exit 0xC0000409 during the wirebasket factorization, three runs) | | | | |
| 5.5T | 1,456,057 | AMS, 4 cycles, edge wirebasket | 584 | 4.0 | 180.2 | 192.9 | 2.58 |

The 5.5T direct failure is not the wirebasket size: with `--edge-wirebasket`
(331,595 wirebasket dofs) it exits the same way, while the 3.5T default
wirebasket (353,473 dofs) factors on the same host.  The same exit code has
been seen from sparse Cholesky on a small SPD system on worker-b, so the cause is
unresolved.

Conductor loss agrees with the direct runs to 1e-12 relative and the magnetic
energy to 6e-9 (the solver tolerance is 1e-8).  Another session's 4-thread job
shared worker-a during all runs.  The direct solver's cost is the wirebasket
factorization (setup and memory) and its triangular solves (about 470 ms of
each iteration at 2.5T); PARDISO was slower (5.5 s per iteration).

What the exploration (`summary_worker_a_20260927.json`) settled:

- One AMS cycle per wirebasket solve is too weak: none of the one-cycle runs
  converged within 1000-1500 iterations.  `cycles=k` runs k stationary steps on the wirebasket system,
  which keeps the preconditioner fixed and complex symmetric: the solve time
  is flat from k=2 to k=4 at 2.5T and k=4 holds the iteration count on the
  finer meshes.
- The AMS surrogate needs no extra diagonal shift (`eps`): the system's own
  `eps*nu` mass is in the wirebasket matrix; a 1e-6 relative shift is six
  orders stronger and costs 70 % more iterations (1481 vs 861 at k=2).
- Keeping the face dofs in the wirebasket helps the direct solver (217 vs
  463 iterations) but not AMS, which can only smooth them.
- Extra smoothing, cycle type 7 and a complex relaxation weight did not
  reduce the iteration count.  More than half of an AMS cycle (57 %) is the
  five auxiliary AMG V-cycles on the 23.7k-vertex hierarchies.

### Auxiliary AMG V-cycles (`bddc_coarse/amg_20260927/`)

Per-level timers (`AMSCoarseStats()["amg_levels"]`) showed the nodal (Pi)
AMGs coarsening past the point where a step still shrinks (Pz: 17 levels,
506 -> 504 -> 503 ... rows) and a sparse Cholesky solve of 0.28 ms for the two
right-hand sides of a 427-500-row coarsest level.  With `lean_coarse` (the
wirebasket default) the auxiliary AMGs stop once a step keeps more than 80 %
of a level of at most 5000 rows and solve a coarsest level of at most 1024
rows with a dense inverse (in-place Gauss-Jordan on the caller's
TaskManager), accepted only if it reproduces a test vector; the nearly
singular gradient coarsest level fails that check and keeps sparse Cholesky
(storing that factorization's operator densely diverged).  The dual V-cycle
also runs a loop serially by work (fewer than 16k nonzeros) instead of rows.

worker-a, 8 threads, AMS 4 cycles, edge wirebasket, one binary:

| mesh | lean_coarse | iterations | solve s | nodal AMG s | gradient AMG s |
|---|---|---|---|---|---|
| 2.5T | 0 (median of 4) | 596 | 67.5 | 14.8 | 6.6 |
| 2.5T | 1 (median of 4) | 595.5 | 65.5 | 11.7 | 6.8 |
| 3.5T | 0 / 1 | 783 / 782 | 140.7 / 134.6 | 26.4 / 22.6 | 10.7 / 10.8 |
| 5.5T | 0 / 1 | 582 / 583 | 176.3 / 170.8 | 31.5 / 28.0 | 13.7 / 13.3 |

Before these changes the same host measured 17.2 s nodal and 7.0 s gradient
AMG at 2.5T, so the auxiliary V-cycles cost 24 % less and the solve 5-7 %
less.  One 2.5T `lean_coarse=1` run took 953 iterations (loss differing at
1e-11); three repeats took 595-596 and it was not reproduced.  Also tried and
dropped: float32 level values (all AMGs: 47 -> 70 iterations on a test box;
nodal only: no gain, the levels are cache resident) and the additive 0(1+2)0
AMS cycle (825 vs 595 iterations).

## Independent-host evidence

`compact_ams_results_compute.json` records the 2026-09-11 compute host run. Its
fixture hash and element/DOF counts match the tracked baseline. The run passes
the existing acceptance range: 144 iterations and true residual
9.891989345431059e-11. The system shift is zero; epsilon=1e-6 is applied only
in the preconditioner. Recorded setup/solve times are 1.747/3.213 seconds and
are observations, not performance gates.

The operator reported Radia 4.95.90 and NGSolve 6.2.2606. Those versions are
not embedded in the collected JSON; neither a source commit nor a native
binary hash was captured there. This is independent-host numerical evidence
for the reported runtime, not acceptance of a new main build, MATLAB parity,
or release-quad completion. The original result values are preserved.
