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

The Python API `TotalAP1Box(mesh, current_cf, **settings)` also runs a
meshed-current total-A comparison with `A_total x n = 0` on `outer`.
It assembles `integral_coil J dot v`, without any analytical source field.
The caller owns current continuity and cross-section normalization checks;
the solver's linear/Newton field agreement, zero-current and polarity
controls are covered by the focused tests. The existing CLI engine list
is unchanged.

For raw element-field comparisons, identical coordinates are not sufficient
when a probe lies on a shared face: the lowest-order curl is discontinuous
and point locators can choose different neighbouring cells. Freeze the tet4
connectivity and material ID, evaluate at its centroid (barycentric weights
all 1/4), and verify that the point locator returns that cell. Keep the
original face-probe diagnostic rather than silently replacing it. A tiny
coordinate offset can remain inside a locator's tolerance and is not a
reliable side-selection contract. Report gap samples and material-stratified
samples separately; neither establishes agreement at every mesh element.

Closed meshed conductors can use `radia.meshed_current.solve_closed_coil_current`:
an RT0/P0 mixed solve enforces zero cell divergence and insulating walls, including
walls whose labels are shared with internal material interfaces. The source is
normalized using a complete straight driven leg and its physical length. Only
one face-connected conductor on straight tetrahedra is accepted. Pass the returned
`current` to `TotalAP1Box`; use a caller-owned `ngsolve.TaskManager()`.
A continuous H1 potential gradient only enforces weak current continuity and
must not be assumed to preserve normal flux between cells or at conductor walls.

The material interpolation is part of the comparison identity, not merely the
tabulated values. `SoftIronLaw` retains the monotone PCHIP B(H) inverse.
`PiecewiseLinearIronLaw` explicitly selects linear H(B) between strictly increasing
samples starting at (0,0), with a vacuum-slope continuation above the table.
Neither interpolation is silently substituted for the other.

The reduced-A option `--outer-boundary natural_total` imposes the weak
condition `n x H_total = 0`, instead of the default `source_flux` condition
`A_r x n = 0`. It frees the outer tangential DOFs and adds
`nu_0 integral_outer (n x B_s) . v` to the right-hand side (with the
opposite sign in the nonlinear residual). These are different finite-box
problems. The option does not change the mixed-Omega engine. Compare total
`B`, not gauge-dependent `A`, and record the boundary policy with results.
A curl-free uniform vacuum source must be cancelled by the natural-total
solution; the default source-flux boundary retains it. This signed control
is covered by the focused tests.

**reduced_a** -- lowest-order Nedelec `A_r` with `B = B_s + curl A_r`,
`nograds` gauge plus the `1e-6 nu_0` mass regularisation, `A_r x n = 0` on the
box.  The exact Radia source is projected once to a second-order discontinuous
field on the iron (the only place the reduced right-hand side integrates it)
and evaluated exactly at the observation points.  `--arc-rel-tol` sets
Radia's `PrcArc` relative tolerance of the arc-current source quadrature
(default 1e-9, accepted range [1e-12, 1e-3]); the value is recorded with the
coil in the result JSON.  Nonlinear loop:

* `--nonlinear-method newton` (default): element-constant flux density,
  reluctivity interpolated from a dense sampling of the production PCHIP
  inverse, and differential reluctivity from numerical differentiation of
  that table. The approximate Jacobian is
  `nu I + (dH/dB - nu) B B^T / |B|^2`, Armijo backtracking on the true
  residual.  Stops when the relative residual is below `--newton-tolerance`
  and the element flux change is below `--tolerance` times `B_sat`.
  The interpolation also enters the residual; it is not an exact evaluation
  of the production inverse. The unit test bounds its discrepancy at six
  sampled flux densities to relative tolerance 1e-4.
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
  `--gauge-epsilon 0 --ams-beta-zero` selects the pure curl-curl mode:
  gradient correction and its Galerkin matrix/solver are omitted; fine-space
  smoothing and vector nodal corrections remain. This requires a native
  build exposing `beta_zero=True`. Neither the physical operator nor the
  preconditioner receives a mass shift. The caller must supply a compatible
  right-hand side; this option does not project incompatible residuals.
  It is an independent AMS mode, not an incomplete-Cholesky shift.
  The mode follows the omission described by hypre's
  [`HYPRE_AMSSetBetaPoissonMatrix(NULL)`](https://hypre.readthedocs.io/en/latest/api-sol-parcsr.html#c.HYPRE_AMSSetBetaPoissonMatrix);
  Radia still uses its own compact implementation, not a hypre binding.
  With `--gauge-epsilon 0 --ams-preconditioner-shift 0.01`, only the AMS
  hierarchy sees `K + sigma nu_0 M`; CG still solves the original singular,
  compatible `K A = f`. The shifted matrix has independent storage and is
  updated with the tangent. This mode uses the existing compiled AMS; it
  does not disable the gradient correction or require a native rebuild.
  `--ams-project-gradients` optionally wraps this as `P0 B_shift P0`, where
  `P0` is the Euclidean complement of gradients from interior H1 nodes. This
  diagnostic uses a cached sparse nodal factorization, leaves `K` and `f`
  unchanged, and still checks their original residual. It costs additional
  time and memory. On the coarse nonlinear case, the unprojected shifted
  route passed the loose `1e-3` Newton residual gate but stalled at a `1e-8`
  gate; the projected route reached the tighter gate. This is a validation
  option, not a new production solver default.
* `iccg`: Radia's compiled shifted incomplete-Cholesky CG of the same module
  (`SparseSolvSolver("ICCG")`, `--ic-shift` 1.05, ABMC parallel triangular
  solves), refactorised per system.  It wants the **ungauged** system
  (`--gauge-epsilon 0`): the reduced right-hand side lies in the range of the
  singular curl-curl operator, so CG converges on it, whereas the ε-mass term
  used by the default AMS route makes IC(0)-CG stall in this lane. Diagonal scaling together with ABMC
  ordering multiplied the iteration count twentyfold in the diagnostic, so it
  stays off.
* `direct`: METIS SPD PARDISO, the cross-check.

Every iterative solve must meet the **true** relative residual on the free
DOFs, below `--cg-tolerance`. AMS uses NGSolve's CG recurrence with this check
at every iteration. Its preconditioned residual norm is not the acceptance
criterion: stopping on that norm and restarting with successively tighter
relative tolerances can over-solve a singular system and amplify roundoff.
An iteration limit is a failure, not convergence. Shift selection and timing
remain case-dependent; a successful shifted solve is not a general performance
claim. The IC factorization shift and the AMS mass shift are different operations.

The CLI writes `completed` and `passed` in its result. `passed` means all
requested engines converged, not that an accuracy or timing certificate was
met. Nonconvergence exits nonzero after saving the results. Runtime exceptions
also save a failed JSON with the exception and options before propagating.

### Beta-zero measurements (2026-09-23)

The [measurement record](results/beta_zero_ams_20260923.json) compares two
AMS modes on mdx2, NGSolve 6.2.2606, MKL with one thread. Both use a true
linear residual tolerance of `1e-6`, nonlinear residual tolerance of `1e-8`,
and flux-change tolerance of `2e-5`. Values below are single measurements.

| Elements | NGSolve threads | Beta-zero Newton / whole process (s) | Penalty AMS Newton / whole process (s) |
|---|---|---|---|
| 67,397 | 2 | 9.85 / 41.97 | 12.55 / 42.39 |
| 67,397 | 8 | 5.64 / 18.00 | 7.10 / 19.60 |
| 195,309 | 8 | 21.51 / 72.44 | 27.09 / 78.89 |

Newton time includes assembly, preconditioner setup, CG and line search,
but excludes source projection and element helper setup. Beta-zero's
67,397-element Newton time was 13.99 s at one thread. On the fixed field
stencil the two modes differ by relative L2 norm below `5.4e-8`; this measures
agreement between these Radia modes, not error against an exact solution.
Both large-mesh two-thread runs reached the 150 s wall-time limit, so their
speed and convergence are not established. A separate small-mesh beta-zero
run with linear tolerance `1e-9` stalled at Newton step 7; omitting gradient
correction does not eliminate finite-precision compatibility problems.

The native/Python option is opt-in. MATLAB can pass it through the existing
`radia.python.sparsesolv` batch entry point; the native MEX AMS options and
the production solver defaults are unchanged.

**mixed_omega** -- first-order total/reduced scalar potential of
`radia.kelvin_solver` on the same mesh: air and coil are the reduced region,
iron the total region, the total-Hodge split supplies the interface trace and
the harmonic remainder, natural `B.n = 0` on the box and the `GND` point gauge.
Its default nonlinear loop is the production Picard iteration (damped or
Anderson). `--omega-nonlinear-method newton` selects quadrature-based PCHIP
coenergy Newton, with residual backtracking and both equation-residual and
vector-field-change convergence checks. `--omega-order 2` enables the P2
scalar-potential comparison; reduced-A remains first order. Newton does not
use Anderson mixing. These switches do not imply a timing certificate.

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

### Optional Inexact Newton Inner Solves

`--inexact-linear --ams-update-every 2` enables an experimental static
Newton configuration. The inner true-relative-residual target is
`max(cg_tolerance, min(0.01, 0.1 * nonlinear_residual / initial_residual))`.
The final nonlinear residual and flux-change gates remain unchanged; the
fixed inner tolerance is restored after every solve, including failures.
Each history row records `linear_tolerance` alongside the achieved true
residual. This is opt-in and does not establish suitability for complex
frequency-domain systems. Compare both preparation-inclusive and engine
times, and verify fields before interpreting any speedup.

`--cg-check-interval N` optionally reduces repeated original-equation matrix
products during AMS-CG. The default is 1. The initial and iteration-limit
checks are mandatory, as is the independent true-residual check after CG
returns. Skipped checks never authorize convergence. A larger interval may
add Krylov iterations and is not necessarily faster. History records the
number of true-residual checks; residual samples are not an every-iteration
trace when N exceeds 1. Other linear backends reject nondefault intervals.

`--ams-print-level 1` enables native setup and application timing output.
Redirect stdout to retain it. Application breakdowns sample the first 25
applications after each hierarchy construction/update, not the entire solve.
Gradient/nodal correction times include residual computation, restriction,
AMG application and prolongation; they are not isolated coarsest-level solve
times. Keep the JSON phase totals as the whole-run timing authority.
Native builds with correction subphase diagnostics additionally report
residual, restriction, auxiliary-cycle and prolongation samples. These
subphases are contained within gradient/nodal correction times and must not
be added again to the overall total. The auxiliary-cycle time still includes
all AMG levels; it is not a coarsest-level-only measurement.

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
