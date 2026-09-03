# HDiv-Type VIM

HDiv-VIM is Radia's production magnet and soft-iron demagnetization route.  It keeps the
magnetization, material state, charge map, and reduced-FEM handoff in NGSolve
mesh/function-space vocabulary, then uses Radia's C++ charge-Gram H-matrix for
the open-boundary integral operator.

Design documents in this directory:
[PRODUCTIONIZATION.md](PRODUCTIONIZATION.md) (productionization plan),
[TOPOLOGY_OPTIMIZATION_ISOCHRONOUS.md](TOPOLOGY_OPTIMIZATION_ISOCHRONOUS.md)
(topology optimization of isochronous-magnet iron: verified foundations,
gL objective/adjoint formulation, staged execution plan, CUDA lane).
Result-bearing showcase notebooks:
[isochronous_topopt.ipynb](isochronous_topopt.ipynb) (density adjoint gate,
constrained trust-region SLP design loop, exact-void iron-only verification;
sidecar `isochronous_topopt_result.json`; golden-band lane
`validation_test/isochronous_topopt/`), and
[c_type_three_formulation_tosca_mixed.ipynb](c_type_three_formulation_tosca_mixed.ipynb)
(Cubit/ACIS C-yoke; HDiv-MMM, HCurl reduced-A, and H1 TOSCA mixed
total/reduced Omega fixed-mesh acceptance; sidecar
`c_type_three_formulation_tosca_mixed_result.json`; validation lane
`validation_test/c_type_three_engine/`).

## NGSolve Family Convention

Radia's established production path constructs `HDiv(mesh, order=1|2)` without
the `RT=True` flag.  In NGSolve this is the Brezzi--Douglas--Marini family, so
the production spaces are **BDM1 and BDM2**.  Raviart--Thomas is a different
family and must be requested explicitly as `HDiv(mesh, order=p, RT=True)`.
Historical Radia material that called the production spaces RT1/RT2 was using
the wrong name; the solver and stored results were BDM.

## Discrete de Rham and charge-map contract

In three-dimensional finite-element exterior calculus, magnetization is a
2-form whose vector proxy belongs to `HDiv`.  NGSolve supplies the compatible
discrete sequence

```text
H1 --grad--> HCurl --curl--> HDiv --div--> L2,
                        div curl = 0.
```

`HCurl` is not an additional state unknown in an HDiv-VIM material solve.  It
is the upstream space that characterizes charge-free HDiv modes: a
boundary-compatible `curl(a)` lies in the kernel of the complete magnetic
charge map.  Radia discretizes that map as

```text
q = B m = (-div M, normal trace M.n),
```

with interface jumps included when a deliberately broken HDiv space is used.
Sharing oriented normal-trace moments on an ordinary interior facet prevents a
numerical charge sheet there, but sharing the facet DoF is a consequence of
HDiv conformity rather than the definition of the de Rham structure.  A
physical normal jump is represented by separate body spaces or by the checked
`internal_interfaces=True` broken-space path.

The Coulomb energy is pulled back from charge space to magnetization space:

```text
E_demag(m) = 0.5 * (B m)^T G (B m),
N = B^T G B.
```

Consequently `N` is symmetric positive semidefinite and vanishes on `ker(B)`.
The production solve does not construct or constrain a loop/cotree basis; the
positive material Hodge operator `W` selects the unique physical solution of
the symmetric positive-definite system `W + B^T G B`.  On a multiply connected
domain, harmonic/cohomology modes remain a topology-owned finite-dimensional
part of `ker(B)` and must not be confused with spurious charge modes.

## Live API

```python
import radia as rad
import radia.vim as vim

iron = vim.MeshSoftIron(mesh, mu_r=1000)
model = rad.ObjCnt([iron, source])
result = rad.Solve(model, 1e-3, 100, 2, demag_backend="hdiv")
field = rad.Fld(model, "b", [0, 0, 0.02])
```

For direct VIM use, call `radia.vim.Solve(mesh, mu_r=... | bh_table=...,
H_ext=..., image=...)`.

`rad.Solve` also accepts a container with multiple independently registered
`vim.MeshSoftIron` bodies.  Eligible pure TET/HEX/WEDGE bodies dispatch through
`vim.SolveCoupled`; each body keeps its own HDiv normal trace, and `rad.Fld` on
the top container sums all persistent C++ BDM1/BDM2 fields plus ordinary Radia
source objects.  A single registered iron retains the ordinary `vim.Solve`
result shape.

## Permanent-magnet model ladder

Radia deliberately exposes four permanent-magnet levels.  Select the least
complex level that represents the physics; all four use the same HDiv charge
and field machinery.

| Level | Model | Production API | Use when |
|---|---|---|---|
| 1 | Fixed/given magnetization | `vim.MagnetizationSource(mesh, M_given)` | The manufactured magnetization distribution is prescribed and does not change. |
| 2 | Linear recoil/demagnetization | `vim.Solve(mesh, mu_r=mu_rec, B_r=B_r, H_ext=...)` | Reversible load-line motion and self-demagnetization are needed, but no hysteretic state is required. |
| 3 | Simplified Play hysteresis | `vim.PlayHysteresisMaterial(...)` with `vim.SolveHysteresis(...)` | A compact engineering history model is sufficient. |
| 4 | Full B-input EnergyStop | `vim.EnergyStopMaterial(...)` with `vim.SolveHysteresis(...)` | Irreversible demagnetization, vector/rotational histories, convex energy structure, and explicit restart state are required. |

Level 2 uses the recoil law

```text
B = mu0 * mu_rec * H + B_r
```

where `B_r` is in tesla and may be either a constant vector or a spatial
NGSolve `CoefficientFunction`.  The implementation is not a Python-side
iteration: it shifts the right-hand side of the existing symmetric C++ HDiv
system exactly, so the unknown and returned field remain the total
magnetization.  `mu_rec` must be greater than one.  The rigid `mu_rec = 1`
limit belongs to level 1 rather than a numerically singular recoil solve.
The spatial field must belong to one physically continuous magnet body.  If
two segments have a jump in normal magnetization, they require separate HDiv
spaces to retain the interface charge.  Level 1 already supports this by using
one `MagnetizationSource` per fixed segment; mutually coupled level-2 recoil
segments use `vim.SolveCoupled` and must not be approximated as one conforming
space.

```python
with ng.TaskManager():
    pm = vim.Solve(
        pm_mesh,
        mu_r=1.05,                         # recoil relative permeability
        B_r=ng.CoefficientFunction((0, 0, 1.2)),  # T, spatial CF allowed
        H_ext=coil_field,                  # optional; zero is the default
    )
```

For a recoil magnet interacting with nonlinear iron, or for segmented recoil
magnets whose normal magnetization jumps at an interface, give every body its
own mesh and HDiv space:

```python
bodies = [
    vim.CoupledBody(pm_mesh, "pm", mu_r=1.05, B_r=(0, 0, 1.2)),
    vim.CoupledBody(iron_mesh, "iron", bh_table=bh_table),
]
with ng.TaskManager():
    coupled = vim.SolveCoupled(bodies)

H = vim.FieldFromCoupledSolution(coupled, observation_points)
```

The block Gauss--Seidel solve builds each geometry-only ChargeGram once, then
reuses it while persistent C++ field CoefficientFunctions exchange the body
fields.  Separate spaces preserve interface surface charge.  Failure to reach
the global coefficient fixed point raises rather than returning a partial
coupling.

Given or manufactured magnetization distributions use a separate source-owned
HDiv space:

```python
import ngsolve as ng
import radia.vim as vim

with ng.TaskManager():
    pm = vim.MagnetizationSource(pm_mesh, M_given, order=1)
    result = vim.Solve(
        iron_mesh,
        mu_r=1000,
        magnetization_sources=[pm],
        H_ext=coil_field,  # optional when the prescribed sources are sufficient
    )

H_pm = pm.Field(points)  # A/m
H_pm_cf = pm.field_cf    # native NGSolve CoefficientFunction
```

`MagnetizationSource` performs a true mass-Riesz/L2 projection of `M_given`
into an independent HDiv space, then materializes the corresponding immutable
C++ charge source without building a ChargeGram H-matrix.  `Solve` assembles
the source field directly into the iron LinearForm and L2-projects the combined
applied field.  PM coefficients are never solve unknowns.  The PM and iron must
use separate mesh objects/spaces; this preserves the physical normal jump and
surface charge even when the bodies touch.  Multiple sources superpose, and
the result retains them in `_magnetization_sources`.

The 3D source supports BDM1/BDM2 TET/HEX/WEDGE, Curve(2), and IMA on the same
geometry/field kernels as the solve.  Planar 2D keeps its established
`magnets=[(mesh, M), ...]` path.  This API represents a fixed prescribed
magnetization; it does not advance a material history.
If `M_given` itself has an internal normal discontinuity (for example, distinct
magnet segments), represent each segment as a separate `MagnetizationSource`;
one conforming HDiv space intentionally enforces normal continuity inside its
own source mesh.

For a permanent magnet whose state changes under reverse field, use the C++
vector B-input Stop law on the PM mesh:

```python
material = vim.EnergyStopMaterial(
    eta_T,
    g_tables,                 # [(r_T, g_A_per_m), ...], monotone per branch
    alpha=5.0,                # positive reversible reluctivity floor
    gamma=0.0,                # hard Stop; positive values use the convex prox
    b_max=1.5,                # calibrated operating limit (optional)
)

with ng.TaskManager():
    solver = vim.HDivSolver(pm_mesh, order=2)
    result = solver.SolveHysteresis(
        applied_field_steps,  # 3-vectors or NGSolve CoefficientFunctions
        material=material,
        initial_b_path=manufacturing_B_history,
    )
    # Continue without replaying or hiding the constitutive history.
    continued = solver.SolveHysteresis(
        next_steps, material=material, initial_state=result["state"]
    )

H_demag = vim.FieldFromSolution(continued, observation_points)
```

Each Stop state lives in the fixed ball `|s_k| <= eta_k`.  Non-negative,
non-decreasing radial tables give convex branch energies; malformed tables are
rejected at construction.  Trial `forward(B, state)` calls are pure, and state
is committed only after the coupled HDiv step converges.  The C++ batch kernel
runs under TaskManager with the GIL released.  `initial_b_path` is an explicit
constitutive initializer for a manufactured magnet state, not a substitute for
modelling the magnetizing fixture.  Reverse-field/unload and restart are locked
by `validation_test/hysteresis/test_energy_stop_irreversible_pm.py`.
The returned final state owns the same persistent C++ field evaluator as
`vim.Solve`, so `FieldFromSolution` evaluates its external demagnetizing field
without collapsing the BDM1/BDM2 magnetization to element constants.

This history-dependent path currently solves PM self-demagnetization under an
arbitrary prescribed NGSolve applied field.  BDM1 keeps one committed
constitutive state per element.  BDM2 stores and updates state on an NGSolve
`IntegrationRuleSpace` and returns the constitutive source through the matching
weak-form transpose; the public step result remains element-averaged for a
stable reporting contract.

For a stateful EnergyStop/Play magnet interacting with nonlinear iron, use the
history-specific block contract:

```python
pm = vim.CoupledHistoryBody(
    pm_mesh, "pm", material, order=2,
    initial_b_path=manufacturing_B_history,
)
iron = vim.CoupledBody(iron_mesh, "iron", bh_table=bh_table, order=2)

with ng.TaskManager():
    coupled_history = vim.SolveCoupledHysteresis(
        [pm], [iron], applied_field_steps,
    )

H = vim.FieldFromCoupledHysteresis(coupled_history, observation_points)
```

At every physical field step, each outer PM trial restarts from that PM's same
previously committed constitutive state.  All trial states are accepted
together only after the global PM/iron coefficient fixed point converges, so
block iterations cannot advance any hysteresis history multiple times.  Each
body owns one `HDivSolver`; its ChargeGram is built once and reused across
outer iterations and history steps.  Add further independently meshed
`CoupledHistoryBody` objects to the first list for segmented stateful magnets.

For level 3, construct `vim.PlayHysteresisMaterial(K, eta, f_k_tables)` and
pass it as `material=` to the same `SolveHysteresis` stepping API.  It retains
branch history and is useful as the simplified engineering model.  It is not
the level-4 claim: EnergyStop additionally fixes every vector Stop state to a
bounded domain, derives the branch update from a convex energy/proximal law,
and exposes the explicit manufacturing/restart state used for irreversible
demagnetization studies.

The public VIM operator and material-solve contract supports BDM1 and BDM2 on
pure TET/HEX/WEDGE meshes, including mapped/non-affine and true Curve(2) HEX
BDM2.  For that
lane the C++ ChargeGram uses complete-host tensor source rules for smooth
pairs and reflection-invariant whole-host Duffy rules for self and adjacent
volume/boundary pairs.  Volume and surface charge therefore share one complete
host representation instead of independently subdivided clouds, preserving
the large BDM2 cancellation.  WEDGE BDM2 uses the corresponding tri-P2 by
z-P2 charge basis.

The mapped-HEX production gate is recorded by
`validate_mapped_hex_bdm2_production.py`.  On mdx its 756-DoF non-affine
trilinear body had generalized spectrum `[-8.53e-16, 0.999899]` with the
q9/q12 default,
linear and nonlinear solves converged, and the q9/q12 material response was
within `5.28e-4` in the mass norm of q10/q16.  A separate q10/q16 versus
q11/q20 run measured `3.94e-4`; the default remains q9/q12 because q10 and q11
cost about 3x and 7.5x more build time without changing the material response
at the 0.1% level.  Mapped BDM2 shape derivatives fail loudly until this
composite quadrature has a consistent differentiated form.

The distinct true-curvature gate regenerates a four-HEX cylinder through
Cubit 2025.12 with Netgen export order 2.  `check-vol` confirms all four cells
are non-affine with mesh curve order 2 and minimum scaled Jacobian `0.33473`.
On mdx, BDM2 has 396 DoF, the linear mass-Riesz CG solve converges in 19
iterations, the equivalent nonlinear Energy-Newton solve converges, and the
prescribed-source `rad.Fld` agrees with an independent NGSolve boundary
integral to `1.06e-10` maximum relative error.  The generated `.vol` remains
an external validation artifact; the Cubit builder and result JSON are tracked.

`radia.vim.H1HodgeDemagOperator(hdiv, h1, definedon=...)` remains an independent
finite-domain H1 Omega cross-check and assembles `C.T K^-1 C` without forming a
dense parent matrix.  Its shared-mesh BDM2 response reduction feeds HCurl
eddy-bubble solve with a snapshot residual below `3e-12`, mixed residual below
`3e-16`, and positive Joule loss.  This proves the response/coupling path and
does not replace the open-boundary ChargeGram accuracy gate.  Partial-material response solves use
the HDiv space's `FreeDofs()` in the generic NGSolve CG path.  The one-call
mixed builder accepts `hdiv_definedon=mesh.Materials("body")` together with
`demag_operator_factory=lambda hdiv: H1HodgeDemagOperator(...)`; the factory
receives the exact restricted BDM space built by `NgsolveBDMEddyBubbleVIM`.

Observable accuracy is checked separately from that mixed-mechanics gate.  On
the static circular-disk lane, the fine axisymmetric Q2 reference is
`1.0916977441`.  The 3-D finite-domain H1-Hodge sequence reduces error from
4.65% (coarse BDM1/H1-p2) to 2.24% (coarse BDM2/H1-p3), 1.44% (fine
BDM2/H1-p3), and 0.98% (fine BDM2/H1-p4).  This establishes a strict h/p
accuracy ladder for that geometry; it is not a universal open-boundary claim.

For a magnetic and conductive body, use the same material label in
`SharedMeshMaterialModel.magnetic_regions` and `conductive_regions`.
At frequencies where `EddySIBCApplicability` selects the surface route,
`mixed.solve_frequency_local_esim(model, frequency)` now places the local
nonlinear ESIM Karl iteration around the complete HDiv-MMM/HCurl-VIM solve.
Each update reconstructs the solution-shaped tangential field, solves or
interpolates the one-dimensional nonlinear B-H skin cell, assembles the local
`SurfaceImpedanceGram`, and resolves both magnetic and eddy-current blocks.
Nonlinear superposition is invalid, so this API accepts exactly one physical
excitation at a time.

The regenerated TET iron-cube validation uses the same `iron` region for BDM1
magnetization and HCurl bulk/bridge/SIBC current.  Its 50, 1000, and 5000 A/m
ladder converges in 2--5 Karl updates, keeps every local Gram passive, produces
positive Joule loss, and replays the converged Gram in the fixed mixed solver
to machine precision.  The mean surface resistance decreases with increasing
field and the departure from linear SIBC grows, as expected for the supplied
monotone saturating B-H curve.  This establishes nonlinear skin coupling; the
ordinary bulk HDiv operator still uses a fixed low-field scalar permeability.
A simultaneous update of bulk nonlinear B-H and local ESIM is a separate
double-nonlinear constitutive iteration and remains open.

Geometry order and HDiv order are independent Piola-FEM choices.  The
authoritative table is `radia.vim.hdiv_capabilities()`: in 2D, BDM1 supports
geometry orders 1/2 (Q2 recommended) and BDM2 supports 1/2/3 (Q3 recommended);
in 3D, TET/HEX/WEDGE BDM1 and BDM2 expose geometry orders 1/2 (P2/Q2
recommended).  A combination outside the table fails loudly.  In particular,
BDM1/Q3 is intentionally excluded in 2D: under the Q3 contravariant Piola map,
the BDM1 space does not reproduce a uniform physical field to roundoff, and an
ellipse torque check showed no accuracy gain over BDM1/Q2.
For a geometrically and topologically symmetric reduced/full hex pair, the
persistent field evaluator must agree at the roundoff contract (`< 10 eps`
relative error), not merely within a percent.  The mapped BDM2 prescribed-source
gate is `2.75 eps` on mdx.  Two independently converged mass-Riesz CG material
solves currently agree to `3.28e-13` in the sampled field; that is a separate
Krylov/reduction error budget and is not represented as a field-evaluator pass.

Every supported 3D solve also materializes one immutable C++ field evaluator.
TET BDM1/BDM2 sources retain analytic volume/triangle kernels; HEX/WEDGE and
curved sources retain the NGSolve quadrature cloud.  Repeated `rad.Fld` calls pass contiguous NumPy target
arrays directly to that evaluator, with no repeated source packing.  IMA terms
are accumulated inside the same TaskManager region.  Ordinary batches use the
exact direct source sum.  Only very large source-target work is considered for
the quadrupole source tree, and auto selection requires both a direct-reference
probe below the configured tolerance and a measured speed benefit.  The result
records `field_evaluator_stats` and `field_evaluator_build_wall_s`.  IMA stays
on the direct evaluator in automatic mode so the reduced/full roundoff contract
is not weakened by different source-tree truncations.

The implementation follows NGSolve's Python-front-end/C++-execution boundary
and treats NGSolve as the finite-element source of truth.
Python declares `HDiv`, `L2`, `SurfaceL2`, coefficient functions, and bilinear
forms and prepares the one-time sparse charge topology.  NGSolve assembles the
forms in C++; pybind extracts their native sparse matrices directly.  The
persistent C++ operator then owns B/BT, geometric and material mass matrices,
the NGSolve `BaseMatrix`, Krylov iterations, and the immutable field source.
Only vectors and target arrays cross the NumPy boundary; Python and SciPy are
not in the per-iteration solve or repeated-field path.
High-order material states use NGSolve mapped interpolation and
`IntegrationRuleSpace`; Python does not reconstruct HDiv/HCurl orientation,
Piola maps, or hidden local DOF transforms from `CalcShape` and `GetDofNrs`.

The production 3D solve uses symmetric C++ CG on the SPD `W + B^T G B`
system.  The public default is `preconditioner="auto"`: linear solves and
small tet nonlinear solves use the exact mass-Riesz map, while nonlinear
energy-Newton on pure hex/wedge meshes and medium-or-larger tet meshes uses the
exact diagonal of `W + B^T G B`.  The diagonal branch is still the same
symmetric energy-Newton system; it changes only the inner Krylov
preconditioner and avoids one PARDISO phase-33 mass solve per CG iteration.
The result artifact records both `preconditioner_requested` and the resolved
`preconditioner`, plus `preconditioner_policy` so benchmark JSON shows why
`auto` chose that branch.  The current tet switch point is `6000` HDiv DOF
(`RADIA_HDIV_AUTO_JACOBI_TET_NFACE` overrides it for measurement sweeps).
`preconditioner="mass-riesz"` and `preconditioner="jacobi"` remain explicit
diagnostic overrides.

For flat pure-hex `.vol` meshes, the BDM1/BDM2 charge-basis path follows the
NGSolve reference ordering directly: the Q2 geometry lattice is built from the
linear `.vol` vertices, and the Q1/Q2 shape-moment to monomial map is applied as a
cached block-diagonal sparse transform.  Curved `.vol` meshes use `GetTrafo` as
the geometry source of truth.  Mapped BDM2 materializes immutable complete-host
tensor rules once, so parallel H-matrix entry fill does not recreate them.

The charge-Gram build caches the symmetrized host-pair block
`0.5*(AB + BA^T)`.  Production evaluates both directed FAR quadrature rules;
mapped BDM2 near pairs additionally use the whole-host Duffy construction.
Storing only an upper triangle makes a matrix algebraically symmetric, but a
one-sided finite quadrature rule is not invariant when a reduced image model is
replaced by its explicit reflected mesh.  The bidirectional rule keeps the
prescribed-source multicell HEX full-vs-image `rad.Fld` contract below `10 eps`
at the normal quadrature order.  `RADIA_HDIV_HEX_FAR_ONESIDED`,
`RADIA_HDIV_WEDGE_FAR_ONESIDED`, and `RADIA_HDIV_HO_FAR_ONESIDED` are retained
only for diagnostic timing experiments and must not be used for release
results. The effective cache, quadrature, image-block, profiling, and native
thread settings are recorded in `hmat_stats`; diagnostic paths set
`release_claim_eligible=0`. The C++ linear solve reports `solve_*` timing
fields, while NumPy buffers replace Python lists at pybind boundaries. Hot
block-cache counters remain opt-in so ordinary timing runs avoid per-entry
atomic overhead.

The first mdx cube timing sweep on 2026-07-08 reached structured-hex N=40
(`1.56M` HDiv DOF) in about `355 s` wall time and `116 s` H-matrix build time.
The matching Netgen tet comparison showed that unstructured tet runs need a
tighter ACA tolerance (`gram_eps=1e-8` to `1e-10`) to keep mass-Riesz CG near
25--30 iterations; too-loose `1e-4` compression can inflate CG to hundreds of
iterations or fail at 4000 iterations.  Earlier one-sided FAR timing figures
are not production claims because that rule did not satisfy explicit-reflection
invariance.  The
remaining performance focus is split by regime: charge-Gram build for linear
or small-tet runs, and H-matrix apply count / nonlinear globalization for large
hex-wedge energy-Newton runs.

The 2026-07-09 mdx nonlinear energy-Newton preconditioner sweep showed why the
`auto` policy is element- and size-aware.  For structured hex, the tuned
diagonal branch reduced N=16 from `247.8 s` to `33.7 s` and completed N=20 in
`79.9 s`, with `M_avg_z` matching the mass-Riesz tuned result at the
sub-ppm-to-few-ppm level.  For wedge, N=10 moved from `60.2 s` to `24.5 s`.
For tet, exact mass-Riesz is still faster on small meshes (`maxh=0.20`,
`4314` DOF: `4.22 s` mass-Riesz vs `4.79 s` diagonal), but the diagonal branch
wins by `maxh=0.15` (`8193` DOF: `8.82 s` vs `13.21 s`).  The supporting JSON
and CSV are stored under the manuscript `results/` directory as
`hdiv_energy_newton_preconditioner_mdx_20260709_*`.  The cube benchmark
drivers expose these timing defaults as `--nonlinear-profile mdx-scaling`
(`nl_tol=3e-4`, `newton_continuation=2`, `newton_reuse_tangent_steps=3`,
`preconditioner=auto`); the default `strict` profile keeps solver-regression
settings.

The wedge/prism path uses the same bidirectional FAR block-serving
infrastructure but is less optimized than the structured-hex path.  Wedge also
uses translated host-block reuse by default:
`RADIA_HDIV_WEDGE_TRANS_CACHE=all` / `2` reuses cell-cell, cell-face, and
face-face template blocks; `1` falls back to the older cell-cell-only subset;
`0` disables the translation cache.  The earlier face-bearing drift was traced
to non-deterministic tie breaks in wedge near quadrature.  After stabilizing the
source-site and target-corner choices, a dense N=2 wedge charge-Gram comparison
between `0` and `all` matched to `8.6e-17` relative Frobenius error, and an
N=4 solve matched `M_avg_z` to `3.2e-16` relative while cutting local H-matrix
build from `2.14 s` to `0.315 s`.  The wedge charge-basis build now mirrors the
hex fast path on linear prism meshes: Q2 lattice nodes are interpolated from
mesh vertices, and the small L2/SurfaceL2-to-monomial transforms are cached and
applied as sparse block transforms instead of solving per element.  On mdx this
cut wedge N=12 basis construction from about `59 s` to `0.75 s`.  With both
changes, mdx N=8/N=12 translation-cache A/B at `gram_eps=1e-8` moved from
`12.0/44.2 s` to `2.52/6.59 s` wall time; H-matrix build moved from
`9.31/34.84 s` to `1.04/2.34 s`, with relative `M_avg_z` drift
`4.8e-13/2.5e-8`.  The follow-up apply/solve profile added optional
`hmatvec_*` timing fields (`RADIA_HDIV_HMATVEC_STATS=1`) and limits MKL's
thread count inside HACApK leaf dgemv calls through
`RADIA_HACAPK_MATVEC_MKL_THREADS` (default `1`) so TaskManager remains the
outer parallel layer while PARDISO mass-Riesz solves can still use the process
thread setting.  On mdx, the updated all/default wedge sweep reached N=20
(`293.6k` HDiv DOF) in `25.3 s` wall time with `6.7 s` H-matrix build and
`8.1 s` solve time; the charge-Gram H-matvec portion of that solve was
`0.64 s`, leaving PARDISO mass-Riesz factor/apply as the next solve-side
target.  Wedge also tolerated `gram_eps=1e-4` in earlier cube checks, cutting
build further while moving `M_avg_z` by only `2e-5`--`4e-5` relative to the
`1e-8` reference.  Treat that as a promising wedge timing lane, not yet the
universal default for all wedge geometries.

## Why This Is The Main Route

- `N = B^T G B` is symmetric and loop-free by construction: loops are
  `ker(B)`, so charge-free modes are field-null without hand-built loop bases.
- NGSolve `Mesh`, `GridFunction`, `CoefficientFunction`, `BilinearForm`, and
  `TaskManager` are shared with reduced FEM and motor workflows.
- Prescribed magnetization is projected once into a source-owned HDiv space;
  its native C++ field CoefficientFunction couples to a separate iron space
  without sampling through Radia objects.
- Curved/high-order geometry uses the same finite-element geometry path instead
  of translating through a separate object discretization.
- TET/HEX/WEDGE and 2D planar validation live in `validation_test/feec/`.

## Validation Expectations

Fast tests belong in `tests/`; heavier numerical checks belong in
`validation_test/feec/` and JSON artifacts labelled with the actual validation
host (`mdx` or `hibino`).  Required HDiv gates:

- analytic demag factors where available;
- nonlinear BH convergence metadata;
- image symmetry checked against an explicitly mirrored full model on truly
  symmetric meshes;
- `rad.Fld` after `rad.Solve(..., demag_backend="hdiv")`;
- persistent-field direct/tree accuracy and scaling through
  `validation_test/feec/bench_hdiv_field_evaluator_scaling.py` on an idle
  compute host;
- prescribed-source L2 residual, direct-field/native-CF equality, source
  immutability, superposition, and iron-response equality;
- energy-Stop hard projection/proximal stationarity, non-negative vector-loop
  dissipation, reverse-field remanence loss, and state-restart reproducibility;
- BDM1/BDM2 flat/curved TET/HEX/WEDGE operator accuracy and cost, including
  mapped-HEX BDM2 spectrum, quadrature convergence, linear/nonlinear solve,
  IMA, and fail-loud shape-derivative gates, plus charge-Gram H-matrix build stats and
  memory/timing on hibino first, or on mdx only when hibino is unavailable and
  the mdx CI queue is idle, for large runs;
- 2D planar motor saliency checks for the motor lane.
- reduced-motor torque agreement among Maxwell stress, magnetization-volume
  coupling, and fixed-current coenergy while reusing one charge Gram;
- native planar source/target-frame `CoefficientFunction` agreement with the
  explicit rigid-coordinate transform.

## Public Docs

Current user-facing notebooks in this directory should be result-bearing and
paired with synchronized JSON sidecars.  Old migration notes and comparison
archives are not part of the public docs surface; recover them from git history
only if needed.
