# C-type dipole, first order, finite box: reduced-A with AMS and mixed Omega

This lane asks a solver question, not an accuracy-certificate question: on
one finite-box mesh of straight tetrahedra at first order, with the coil
meshed, how fast and how accurately do Radia's two static formulations reach
the C-type gap field?  It complements
[`../c_type_three_engine`](../c_type_three_engine/README.md), which owns the
open-boundary (Kelvin) accuracy certificate on Cubit meshes and forbids a
finite air box.  Everything physical is shared with that lane: the ESRF
Example-5 C-yoke, the rounded-rectangle solid current of `CoilBuilder`
(-2000 A), the CEFC B-H table `em_sample_bh.txt` through the same monotone
PCHIP B(H) law, and the 81-point gap stencil with its median-plane projection
and 10 mm gap core.

## Mesh

`build_box_mesh.py` builds the iron of `radia.esrf_examples` (Netgen/OCC, the
explicit alternative to the Cubit authority), the swept coil solid, a refined
air block spanning the gap, and an air box of half-width 0.5 m.  Materials are
`iron`, `coil`, `air`; every iron face is `iron_air_interface`, the box faces
are `outer`, and one box corner is the `GND` vertex.  `--scale` multiplies all
mesh sizes; the `.json` beside the `.vol` carries counts, labels and the
SHA-256 that `run_p1_box.py` verifies.

| scale | tetrahedra | vertices | HCurl P1 unknowns |
|---|---|---|---|
| 3.0 | 48,297 | 8,309 | 56,690 |
| 1.5 | 67,397 | 11,603 | 79,209 |
| 0.75 | 195,309 | 33,552 | 229,732 |
| 0.5 | 530,833 | 90,791 | 623,564 |

## Engines

**reduced_a** -- lowest-order Nedelec `A_r` with `B = B_s + curl A_r`,
`nograds` gauge plus the `1e-6 nu_0` mass regularisation, `A_r x n = 0` on the
box.  The exact Radia source is projected once to a second-order discontinuous
field on the iron (the only place the reduced right-hand side integrates it)
and evaluated exactly at the observation points.  Nonlinear loop:

* `--nonlinear-method newton` (default): element-constant flux density,
  reluctivity from the exact inverse of the PCHIP law and its differential
  reluctivity from that tabulated inverse, Jacobian
  `nu I + (dH/dB - nu) B B^T / |B|^2`, Armijo backtracking on the true
  residual.  Stops when the relative residual is below `--newton-tolerance`
  and the element flux change is below `--tolerance` times `B_sat`.
* `--nonlinear-method picard`: the damped / constrained-Anderson update of
  the three-engine lane (kept for the cross-check; on this mesh family it
  needs about a hundred solves at relaxation 0.1 and oscillates at 0.3, in the
  production solver exactly as here).

Linear solves (`--reduced-a-solver`):

* `ams` (default): Radia's compiled auxiliary-space Maxwell preconditioner
  (`radia.sparsesolv_ngsolve.HypreBasedAMSPreconditioner`, built once and
  updated in place at every Newton step, constructed outside `TaskManager` as
  its contract requires) with conjugate gradients on the mass-regularised
  system (`--gauge-epsilon`, default 1e-6).
* `iccg`: Radia's compiled shifted incomplete-Cholesky CG of the same module
  (`SparseSolvSolver("ICCG")`, `--ic-shift` 1.05, ABMC parallel triangular
  solves), refactorised per system.  It wants the **ungauged** system
  (`--gauge-epsilon 0`): the reduced right-hand side lies in the range of the
  singular curl-curl operator, so CG converges on it, whereas the ε-mass term
  that AMS needs makes IC(0)-CG stall.  Diagonal scaling together with ABMC
  ordering multiplied the iteration count twentyfold in the diagnostic, so it
  stays off.
* `direct`: METIS SPD PARDISO, the cross-check.

Every iterative solve is continued until the **true** relative residual on the
free DOFs is below `--cg-tolerance` (the Krylov loop's own measure is the
preconditioned norm and is not the contract).

**mixed_omega** -- first-order total/reduced scalar potential of
`radia.kelvin_solver` on the same mesh: air and coil are the reduced region,
iron the total region, the total-Hodge split supplies the interface trace and
the harmonic remainder, natural `B.n = 0` on the box and the `GND` point gauge.
Its nonlinear loop is the production Picard iteration (damped or Anderson).

## Reference and metrics

`--reference` takes a three-engine result JSON of the Kelvin lane.  The
comparison reports, for every engine here and every engine there, the
gap-core relative RMS of the median-plane-projected B and the centre field
difference.  That deviation bundles the box truncation with the first-order
discretisation and is not a statement about either alone; the Kelvin lane's
certificate bounds the reference itself to 0.35%.

## Timing contract

One process, `--threads` fixed (8 on the development host), mesh loaded once,
source evaluated once, phases timed separately (assembly, preconditioner
update, Krylov solve, line search).  A timing run is only valid when nothing
else computes on the host; runs that overlapped other jobs are kept as
convergence evidence, not as timings.

## Run

```powershell
python validation_test/c_type_p1_ams_box/build_box_mesh.py --output C:/temp/ctype_box/meshes/ctype_box_s075.vol --scale 0.75

python validation_test/c_type_p1_ams_box/run_p1_box.py `
  --vol C:/temp/ctype_box/meshes/ctype_box_s075.vol `
  --output C:/temp/ctype_box/results/s075_newton.json `
  --mode nonlinear --threads 8 --engines reduced_a `
  --reference validation_test/c_type_three_engine/results/c_type_20260903_nonlinear_bdm2_finer_mdx.json
```

`tests/test_c_type_p1_box_lane.py` pins the lane's contracts on a small
synthetic mesh: the flux-side law equals the production inverse, AMS+CG
matches the direct solve under the true-residual contract, Newton is one exact
step on a linear law and converges with full steps on the sample table, Picard
and Newton agree on a mildly saturated cube, and the mesh contract rejects a
mesh without the lane's labels.

## Results (LAB, 2026-09-22, uncontended sequential runs, 8 threads)

`results/lab_20260922_<mesh>_<variant>.json`, meshes reproducible from
`results/meshes/*.json` (scale, sizes, SHA-256).  "tight" stops at relative
residual 1e-6 and |dB|/B_sat 2e-5; "loose" at relative residual 1e-3 with no
flux-change condition.  "no-coil" is the same magnet on a mesh whose coil
volume is air, which the reduced formulation allows.  The accuracy column is
the gap-core relative RMS against the Kelvin lane's finest order-2 reduced-A
field (open boundary; certificate envelope 0.35%).

| mesh | tets | HCurl P1 DOF | tight (s) | Newton / CG | loose (s) | Newton / CG | vs Kelvin reduced-A | centre Bz (T) |
|---|---|---|---|---|---|---|---|---|
| s15 | 67,397 | 79,209 | 11.2 | 11 / 779 | 7.8 | 7 / 510 | 0.299% | 0.249187 |
| s075 | 195,309 | 229,732 | 33.0 | 11 / 773 | 22.3 | 7 / 508 | 0.068% | 0.249321 |
| s05 | 530,833 | 623,564 | 95.9 | 11 / 827 | 62.9 | 7 / 539 | 0.045% | 0.249327 |
| s075 no-coil | 147,822 | 174,190 | 28.7 | 10 / 710 | 13.9 | 6 / 441 | 0.064% | 0.249358 |
| s05 no-coil | 432,455 | 508,616 | 79.9 | 11 / 827 | 51.9 | 7 / 542 | 0.044% | 0.249329 |

The AMS-CG iteration count per Newton step is mesh-independent (65 to 84 from
79k to 624k DOF), so the runtime grows linearly with the element count.  Of
the 95.9 s at s05: AMS update 27.6 s, AMS-CG 47.2 s (57 ms per iteration),
Jacobian assembly 10.4 s, line-search residual evaluations 9.8 s; mesh load
2.7 s and the one-off source projection 9.9 s are outside the engine runtime.

Shifted ICCG on the ungauged system (`results/lab_20260922_*_iccg_*.json`):

| mesh | loose (s) | Newton / ICCG | ms per ICCG iteration | tight (s) | Newton / ICCG |
|---|---|---|---|---|---|
| s15 | 2.9 | 7 / 740 | 2.3 | 42.7 | 11 / 8,783 |
| s075 | 11.0 | 7 / 1,081 | 6.3 | not run | |
| s05 | 91.4 | 7 / 1,456 | 50 | not run | |

At the loose rule ICCG is 2.7x faster than AMS on 67k tets and 2x on 195k, with
iteration counts growing like 1/h as expected of IC(0); on 531k tets it loses
because the cost of one IC(0) triangular solve jumps eightfold for 2.7x
unknowns (2.3 → 6.3 → 50 ms), which points at the parallel substitution in
`sparsesolv` rather than at the algorithm.  The tight rule is not ICCG's
territory: the late Newton right-hand sides of the singular system carry a
consistency floor near 1e-8 and CG grinds thousands of iterations against it.
The damped Picard update of the production solver, run on the same mesh as a
control, did not converge within 80 iterations at relaxation 0.1 and cycled
at 0.3; Newton is what makes first-order reduced-A usable on this mesh family.

The mixed Omega engine on the same box: linear order 1 runs, with the exact
Radia `H_s` evaluated at every quadrature point of the reduced region as the
dominant cost (12 s of 15.5 s at 48k tets, one-off per problem); its nonlinear
Picard loop at relaxation 0.3 did not reach the 2e-5 criterion within 80
iterations on the 67k mesh and a constrained-Anderson run on the 48k mesh was
stopped after 3 CPU-hours without converging.  A Newton loop for that engine is
the open item.
