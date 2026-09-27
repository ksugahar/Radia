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
runs belong on mdx; the tracked 2.5T case is suitable as a LAB smoke run.

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
`C:\temp\radia-validation`, never into the repository.

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

mdx1, 8 threads, 2026-09-27, one binary (`bddc_coarse/`):

| mesh | dofs | coarse solver | iterations | setup s | solve s | wall s | peak GB |
|---|---|---|---|---|---|---|---|
| 2.5T | 680,041 | direct (sparse Cholesky) | 217 | 71.2 | 101.5 | 174.3 | 11.6 |
| 2.5T | 680,041 | AMS, 4 cycles, edge wirebasket | 594 | 2.1 | 70.5 | 77.3 | 1.28 |
| 3.5T | 864,781 | direct | 231 | 109.8 | 156.5 | 268.4 | 15.8 |
| 3.5T | 864,781 | AMS, 4 cycles, edge wirebasket | 783 | 2.5 | 143.1 | 150.7 | 1.58 |
| 5.5T | 1,456,057 | direct | fails (process exit 0xC0000409 during the wirebasket factorization, three runs) | | | | |
| 5.5T | 1,456,057 | AMS, 4 cycles, edge wirebasket | 584 | 4.0 | 180.2 | 192.9 | 2.58 |

Conductor loss agrees with the direct runs to 1e-12 relative and the magnetic
energy to 6e-9 (the solver tolerance is 1e-8).  Another session's 4-thread job
shared mdx1 during all runs.  The direct solver's cost is the wirebasket
factorization (setup and memory) and its triangular solves (about 470 ms of
each iteration at 2.5T); PARDISO was slower (5.5 s per iteration).

What the exploration (`summary_mdx1_20260927.json`) settled:

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

## Independent-host evidence

`compact_ams_results_hibino.json` records the 2026-09-11 Hibino run. Its
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
