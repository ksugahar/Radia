# HDiv-MMM implementation review

This note reviews the production HDiv-MMM implementation across the public
Python API, NGSolve assembly layer, pybind11 boundary, C++/HACApK kernels,
field evaluation, permanent-magnet models, MATLAB boundary, and focused
validation. It also reconciles the ESRF example-5 C-yoke investigation that
exposed several failure modes.

## Review status

- Revised: 2026-08-29 against `origin/main` at `22fc5630e`.
- HDiv hardening anchor: `87309e662` (`fix(hdiv): contain HACApK callback failures`).
- Claude review source: branch `claude/cubit-vol-main-path`, commits
  `2c939a8ed`, `4b30479b6`, and `cf79734b9`. Their measurements and useful
  regressions are reconciled below; none of those three commits is merged on
  `main` as written.
- Scope: correctness, deterministic execution, NGSolve/TaskManager policy,
  public capability, field and material contracts, maintainability,
  performance, MATLAB parity, and reproducibility.

The current source sizes are:

| File | Lines | Role |
|---|---:|---|
| `src/core/rad_hdiv_vim.cpp` | 2,284 | Analytic charge and field kernels |
| `src/core/rad_hacapk_hdiv.cpp` | 11,857 | H-matrix construction and material solves |
| `src/core/rad_hacapk_hdiv_entry.cpp` | 251 | Element-specific entry strategies |
| `src/core/rad_hacapk_hdiv.h` | 1,277 | Charge-Gram state and public C++ contract |
| `src/core/rad_hdiv_field_evaluator.cpp` | 1,260 | Persistent direct/tree field evaluator |
| `src/core/rad_hdiv_hysteresis.cpp` | 244 | Native vector Play/EnergyStop updates |
| `src/radia/vim/_vim.py` | 2,246 | NGSolve charge-map and geometry orchestration |
| `src/radia/vim/_solve.py` | 1,341 | Linear/nonlinear production solve orchestration |
| `src/radia/vim/_hysteresis.py` | 778 | Stateful B-input material stepping |
| `src/radia/vim/_field_batch.py` | 735 | Persistent field and vector-potential surfaces |

Wall-clock performance claims are intentionally limited to measurements made
on `mdx` or `hibino`. LAB is used for build and focused correctness tests.

## Executive assessment

The reviewed core is materially safer than the version that started this
investigation. HACApK callbacks contain exceptions before returning to C;
process-wide callback state is serialized; linear and nonlinear gathers are
deterministic; the heavy build/solve/field bindings release the GIL; and the
C++ kernels self-wrap TaskManager regions while remaining safe inside a
caller-owned NGSolve region. The public solve keeps NGSolve as the owner of FE
spaces, orientation, Piola maps, mapped evaluation, and weak-form assembly.
This Python-orchestration/C++-execution split is NGSolve-native; "complete C++"
must not be misread as moving NGSolve FE plumbing into Radia.

The implementation is production-capable for its declared BDM1/BDM2 lanes,
but it does not yet justify an unconditional "complete for every HDiv case"
claim. Four P1 findings remain:

1. Three explicit full-versus-IMA field regressions fail the required
   `< 10 eps` contract in isolated pytest processes.
2. Broken RT0 was reintroduced into the public capability table and
   `DemagOperator`, contradicting the BDM1/BDM2-only production policy.
3. Fine TET C-yoke runs reportedly lose positive definiteness under refinement;
   the failing mesh and driver are not tracked, so the defect is not replayable.
4. The corrected C-yoke comparison still needs one identical B-H interpolant
   and committed result artifacts before FEM-versus-HDiv accuracy is quoted.

Known, fail-loud limitations are not hidden failures: mapped/non-affine HEX
BDM2 material solves, 3D pyramid/mixed meshes, and 2D history solves remain
unsupported. Large IMA field maps remain direct because the tree route is
disabled whenever images are present. MATLAB owns native field-evaluator and
EnergyStop MEX surfaces, but the high-level `radia.vim` solve is still an
explicit in-process Python fallback because NGSolve setup is Python-owned.

## 0. Cross-layer production review

### 0.1 Findings by priority

| ID | Priority | Finding | Evidence and required disposition |
|---|---|---|---|
| F1 | P1 | The IMA field contract is red. | Isolated measurements were `2.0248220545e-14` and `4.9318891895e-14` relative error for the two HEX `rad.Fld`/`FieldFromSolution` gates, and `2.3931079340e-15` component error for curved TET BDM2. The limit is `2.2204460493e-15` (`10 eps`). Preserve the limit; align solve reduction and source accumulation order rather than loosening it. |
| F2 | P1 | RT0 is publicly advertised again despite the BDM1/BDM2-only decision. | `_capabilities.py` exposes 3D TET/HEX order 0 and `DemagOperator` documents an order-0 broken-interface path. `HDivSolver` and field evaluation accept only orders 1 and 2. Remove the public RT0 entries/path and retain any topology-only experiment outside the production API. |
| F3 | P1 | Fine-TET loss of SPD is an unclosed correctness report. | Claude reported `p^T A p = -1.90e5` at 8.75 mm and `-7.51e8` at 7.0 mm, but the mesh/configuration is absent from `main`. Commit the reproducer before changing quadrature, ACA, or CG. |
| F4 | P1 | The C-yoke accuracy comparison is not self-replaying. | The Cubit journals exist only on the Claude branch and `src/radia/esrf_examples.py` was untracked. The B-H interpolants also differed. Promote a validation driver, one material law, configs, and JSON results together. |
| F5 | P2 | Mapped/non-affine HEX BDM2 is operator-only, not a material solve. | `Solve` rejects it before wrong physics; mapped HEX BDM1, affine HEX BDM2, TET BDM2, and WEDGE BDM2 are the current alternatives. This is correctly documented and tested, but it remains a major completeness boundary. |
| F6 | P2 | IMA disables tree acceleration for field maps. | `HDivFieldEvaluator::AlgorithmFor` returns `Direct` whenever images exist. This protects full/reduced roundoff parity, but large IMA observation maps cannot use the otherwise guarded treecode. Any image-aware acceleration needs a common full/reduced grouping and the F1 contract first. |
| F7 | P2 | Exact vector-potential evaluation is narrower than H-field evaluation. | Exact `A` uses straight TET BDM1 equivalent currents. BDM2, curved, HEX, and WEDGE use NGSolve-mapped quadrature clouds assembled in Python. This is valid as an explicit converged quadrature route, not an all-topology exact/native claim. |
| F8 | P2 | MATLAB parity is partial at the method level. | `HDivFieldEvaluator` and `EnergyStopMaterial` have native checked MEX handles. The parity manifest classifies `vim/__init__.py` as `python-fallback`; complete solve/mesh/form orchestration is not a native MATLAB API. |
| F9 | P3 | Configuration provenance and class ownership remain broad. | Fourteen `RADIA_HDIV_*` variables remain, and `RadHACApKChargeGram` still owns entry, build, cache, solve, derivative, field, and diagnostics. Keep diagnostic switches out of release claims and split only at measured ownership boundaries. |

### 0.2 Production capability matrix

| Dimension/topology | HDiv space | Geometry | Material/history status | Field status |
|---|---|---|---|---|
| 3D TET | BDM1, BDM2 | affine or P2 curved | linear, nonlinear, IMA, recoil PM, Play, EnergyStop | persistent native direct/tree; curved direct leaves retained exactly; exact vector `A` only for straight BDM1 |
| 3D HEX | BDM1 | mapped/affine, geometry order 1 or 2 | linear/nonlinear/IMA/history | affine polynomial decomposition; mapped/curved source cloud |
| 3D HEX | BDM2 | geometry order 1 or 2 | material solve only when the mapping is affine; diagnostic `ChargeGram` otherwise | persistent native evaluator from configured source representation |
| 3D WEDGE | BDM1, BDM2 | geometry order 1 or 2 | linear, nonlinear, IMA, recoil PM, Play, EnergyStop | persistent mapped source representation; native repeated evaluation |
| 2D TRI/QUAD/mixed | BDM1 with geometry 1/2; BDM2 with geometry 1/2/3 | Q2/Q3 production pairings are explicit | linear/nonlinear and IMA; no `SolveHysteresis` | persistent native planar field and `Az` evaluator |
| 3D mixed/pyramid | none | none | fail loud pending NGSolve HDiv pyramid | fixed-magnetization Radia geometry is separate and is not a material solve |

The four permanent-magnet levels are implemented on the same method surface:

| Level | API | Review result |
|---|---|---|
| 1 fixed/given magnetization | `MagnetizationSource` | L2 projection, superposition, IMA reconstruction, TET/HEX/WEDGE and curved BDM2 source tests pass; it is a 3D API. |
| 2 linear recoil | `Solve(mu_r=..., B_r=...)` | BDM1/BDM2 shifted symmetric-system tests and spatial `B_r` pass; segmented magnets require independent spaces. |
| 3 simplified Play | `PlayHysteresisMaterial` + `SolveHysteresis` | Native batched constitutive update, persistent state, BDM1 element and BDM2 quadrature layouts pass. |
| 4 B-input EnergyStop | `EnergyStopMaterial` + `SolveHysteresis` | Convex-table validation, hard projection, stationarity, non-negative loop dissipation, reverse-field irreversible loss, and restart pass. |

### 0.3 Claude branch reconciliation

The Claude branch contains useful evidence, but it must not be merged wholesale:

- Its earlier build-state, deterministic-scatter, and strategy work reached
  `main` under the reviewed `3c015cc46`/`ce162a466`/`87309e662` sequence.
- `4b30479b6` correctly added a negative-side annihilation regression, but its
  implementation accepts every folded diagonal within `1e-12 * max(diag)`.
  That band is much wider than machine roundoff for a local charge and can hide
  a real negative direction. Current `main` instead scales the allowance by
  the image-free direct self-energy, image count, and machine epsilon. The
  negative fixture is retained with that narrower contract.
- `2c939a8ed` supplies useful C-yoke Cubit journals, but its policy text changes
  "build123d or Cubit, chosen by fit" into a universal `Cubit > build123d`
  ranking. That contradicts the current CAD policy and is not accepted. The
  journals should be promoted only with `check-vol`, labels, a tracked driver,
  and result artifacts.
- `cf79734b9` is the original review record. Its useful measurements are
  retained here; stale implementation claims and the broad roundoff fix are
  superseded by this review.

## 1. C++ implementation review

### 1.1 Resolved findings

| ID | Resolution | Evidence |
|---|---|---|
| R1 | HACApK callback failures are caught at the C ABI, stored as `exception_ptr`, and rethrown after fill workers join. Build state is restored by a base-class scope guard. | `87309e662`; `test_fill_exception_restores_chargegram_and_global_hacapk_state` |
| R2 | Concurrent H-matrix builds and destruction are serialized around HACApK's process-wide callback state; each build remains TaskManager-parallel internally. | `3c015cc46`, `87309e662`; six repeated two-thread builds in `test_concurrent_chargegram_builds_keep_callback_state_isolated` |
| R3 | Linear, demag-apply, and nonlinear transpose scatters use deterministic CSR gathers. Fixed 4,096-entry blocks and compensated ordered reduction are used for dot products. | `ce162a466`; no `AtomicAdd` remains in the audited HDiv files |
| R4 | The nonlinear energy-Newton route has a same-process bitwise regression after an intervening solve on another mesh. | `test_alternating_nonlinear_solves_are_bitwise_deterministic` |
| R5 | `GetInteractionMatrixElementRaw` performs one immutable strategy dispatch instead of an eight-mode branch cascade in the hot path. | `rad_hacapk_hdiv_entry.cpp`, `ce162a466` |
| R6 | CG fails immediately on a non-positive or non-finite `p^T A p`, and reports the final true residual instead of disguising breakdown as slow convergence. | `b65f04b7e` |
| R7 | Image-folded diagonal normalization uses an image-free direct self-energy scale to accept only roundoff-sized cancellation. Invalid negative or non-finite diagonals still fail. | `87309e662`; positive- and negative-roundoff fixtures plus genuine-negative rejection |
| R8 | Degenerate triangle kernels return the zero limit instead of dividing by a zero normal. | `b65f04b7e` |
| R9 | The fixed 84-entry HEX affine polynomial scratch is tied to total degree six by `static_assert` and guarded at the shared multiplication choke point. | `b65f04b7e` |
| R10 | Per-entry environment lookups used by high-order/image paths are cached. | `b65f04b7e`, `ce162a466`, `87309e662` |
| R11 | Block-PCG zero-rank deflation uses a typed local exception; unrelated numerical failures can no longer be swallowed because their text happens to match. | `SolveConfiguredLinearMaterialAutoPrecMany` |

R1 needs precise wording. The original review said that exceptions "unwind the
pure-C frames." That behavior was the defect, not the final design. The current
callback catches every exception before returning to HACApK C code, records the
first failure, returns NaN only as an internal abort sentinel, waits for all
workers, frees the failed build, restores global and object state, and then
rethrows on the build thread.

R3 is verified for the audited solve paths, but its test boundary should not be
overstated. The nonlinear bitwise test covers the production energy-Newton path
on a small TET mesh. It does not prove bitwise identity for every element,
curvature, IMA, BLAS, or processor combination; those routes retain their
separate numerical regressions.

### 1.2 Hot-path initialization correction

The earlier review called the scratch change "scratch sizing" and presented
the values as storage reductions. The arrays remain fixed-size stack arrays;
the optimization reduces the bytes initialized on each call. The corrected
description is:

| Scratch | Bytes previously zeroed | Bytes zeroed at degree 2 |
|---|---:|---:|
| `face_moments[4][1330]` | 42,560 | 320 |
| per-moment `poly[19][19]` | 2,888 | 48 |
| `poly2_mul_linear` temporary | 2,888 | 80 |
| `TetMomentMemo::seen` | 6,859 | 27 |

`TetMomentMemo::val` is still allocated but is intentionally left
uninitialized; it is read only after the corresponding `seen` entry is set.
This is an initialization-traffic optimization, not a reduction of the stack
frame to 27 bytes.

### 1.3 Design debt that remains

The entry-strategy extraction completed the high-value part of the former
"god class" finding, but `RadHACApKChargeGram` is still a large owner of build,
cache, element, solve, field-evaluation, and diagnostic state. Future
decomposition should follow ownership boundaries rather than adding another
mode flag:

- immutable element geometry and entry kernels;
- H-matrix build/cache policy;
- linear and nonlinear solve operators;
- persistent field evaluation;
- diagnostics and artifact provenance.

The audited source contains 14 `RADIA_HDIV_*` environment variables. They are
not all numerical switches: two collect cache statistics and one injects a
test failure. The remaining path, cache, quadrature, and preconditioner knobs
must either be diagnostic-only and rejected for release claims, or be surfaced
through `SolverConfig` and copied into result provenance. Calling all 14
"silent numerical behavior" was too broad; leaving the production subset
unrecorded would still be a reproducibility defect.

## 2. ESRF example-5 C-yoke investigation

The figures in this section came from the Claude C-yoke campaign. They are
measurement records, but they are not yet self-replaying evidence on `main`
because `src/radia/esrf_examples.py` remains untracked in the Claude worktree.
The Cubit journals were committed on the Claude branch in `2c939a8ed` and
`4b30479b6`, but that branch has not been merged and the journals are also
absent from `main`.

The starting discrepancy was 1.409849% full-vector B RMS between BDM2 and
Kelvin FEM. The pole-edge region held 95.0905% of the vector-error energy,
while the central flat agreed to 0.110757%. Mesh refinement, quadrature order,
`leaf_size`, Kelvin gap, and IMA count had been varied without closing it.

### 2.1 Pole geometry defect

`_example5_iron()` built the pole chamfer with a three-section
`netgen.occ.ThruSections`. The wrapper does not expose OCC's ruled-loft switch,
so the resulting surface is smooth rather than the intended planar-faced
`ObjMltExtRtg` solid. The 34 x 24 mm to 50 x 40 mm section change has a kink at
z = 13 mm that the smooth loft overshoots.

| z (mm) | Measured area | Ruled reference | Ratio |
|---:|---:|---:|---:|
| 5.0 | 886.6 mm2 | 816.0 mm2 | +8.7% |
| 8.0 | 1,471.6 mm2 | 1,200.0 mm2 | +22.6% |
| 13.0 | 1,999.4 mm2 | 2,000.0 mm2 | approximately 1.000 |
| 18.0 | 2,146.2 mm2 | 2,000.0 mm2 | +7.3% |

The measured pole volume was 7.65% high and its bounding box was 56.93 x
46.93 mm instead of 50 x 40 mm. With the same mesher, material, coil, order,
and 25 mm target size, replacing only this geometry changed full-vector B by
7.2429% RMS. The maximum local change was 22.1482% at (30, 0, 0) mm in the
fringe; the B magnitude at the origin changed by -0.5224%.

The geometry defect alone cannot explain a difference between two routes that
consume the same CAD. The working interpretation is that it activates a
discretization mismatch: curved FEM follows the spline after `mesh.Curve`,
while the tested BDM2 path used affine facets. The measured BDM2 refinement
correction was anti-correlated with the FEM-required correction (cosine
-0.920220). That is evidence for the interpretation, not a proof; the
corrected-CAD rerun is the deciding experiment.

### 2.2 B-H interpolation mismatch

Both routes receive `get_esrf_bh_table(5)` but do not solve the same material
law between table nodes:

| Route | Interpolant |
|---|---|
| FEM | `ng.BSpline(2, ...)`, piecewise linear |
| HDiv-MMM | SciPy `PchipInterpolator`, monotone cubic |

At table nodes they agreed to 1.7e-16. Against the analytic `MatSatIsoFrm`
law over B = 1.7 to 2.1 T, the measured interpolation errors were:

| Route | Maximum relative error | RMS relative error |
|---|---:|---:|
| FEM | 1.460e-3 | 4.285e-4 |
| HDiv-MMM | 5.778e-5 | 1.603e-5 |

Only 20 of 221 table points cover B = 1.7 to 2.0 T because the table is
geometric in `mu0*H` across eleven decades. Differential permeability differed
by as much as 70% near H = 339 A/m. Unlike the CAD issue, this is a direct
cross-route inconsistency. Accuracy must be remeasured only after both routes
use the same interpolation contract.

### 2.3 Coil check

The source audit found no route difference:

| Item | Measured value | Reference |
|---|---:|---|
| Current | -2,000 A | `ex.current = -2000` |
| Closure | closed; 0.0 mm gap | closed loop |
| Racetrack radius/width | 22.5 / 35 mm | inner/outer radii 5 / 40 mm |

Both routes used the same `CoilBuilder.to_radia` solid-current field in the
iron. The coil is therefore not a current explanation for the discrepancy.

## 3. Cubit and HEX/TET cross-check

### 3.1 CAD route

The corrected Cubit/ACIS solid measured 1,446,095.333333 mm3, equal to the
analytic volume at the reported precision, with the intended bounding box and
ten planar faces. `export netgen` preserved that volume to about 1.5e-12%, and
`check-vol` reported no errors.

The review also measured a 14-solid STEP for which
`netgen.occ.OCCGeometry(step)` exposed one solid, while the STEP contained 14
`MANIFOLD_SOLID_BREP` entities and OCP returned all 14. For this benchmark the
supported route is therefore Cubit/ACIS directly to Netgen `.vol`, with
`check-vol`, rather than Cubit to STEP to `netgen.occ`.

### 3.2 Current HEX limitation

The original statement "a chamfered pole cannot be meshed into affine hexes"
was too absolute. What was measured is narrower: the current Cubit
decomposition and mesh family leave cells crossing the 45-degree chamfer as
non-affine trapezoidal prisms. Their mapping residual was 0.745, 0.666, and
0.405 at 25, 12.5, and 6.25 mm, respectively, against the 1e-10 affine gate.
Removing the chamfer reduced the residual to 2.2e-14, but that changes the
benchmark and is not a production remedy.

A united C-yoke solid did not hex-mesh with the tested `auto`, `sweep`,
`webcut_cyl_auto`, or `polyhedron` routes. Decomposition creates
same-material internal interfaces, which the exporter writes with
`DomainIn == DomainOut == 1`. At present the exporter does not omit them.

### 3.3 Measured cross-check

Corrected Cubit/ACIS quarter yoke, `image='+x-z'`, order 1:

| Case | Elements | HDiv DoF | B magnitude at origin | Vector-B RMS vs reference |
|---|---:|---:|---:|---:|
| HEX 8.75 mm, reference | 372 | 9,232 | 0.249965 T | - |
| HEX 12.5 mm | 176 | 4,432 | 0.249908 T | 0.1898% |
| TET 12.5 mm | 1,314 | 8,208 | 0.249708 T | 0.4788% |

At comparable DoF, TET and HEX differed by 0.48% RMS and 0.10% at the centre.
This shows no material corruption from the current same-material interface at
these resolutions, but it does not prove that such interfaces are harmless for
all orders, image groups, or material assignments.

Finer TET cases reported CG breakdown with `p^T A p = -1.90e5` at 8.75 mm and
`-7.51e8` at 7.0 mm. The worsening sign and magnitude under refinement identify
an operator/SPD defect rather than a request for more CG iterations. The TET
C-yoke lane is therefore certified only down to the last passing 12.5 mm case
until the failing Gram/operator is materialized and compared with an
independent dense or NGSolve reference.

## 4. Image-folded roundoff contract

An antisymmetric image can annihilate a charge on its fixed plane. In exact
arithmetic the folded self-energy is zero; finite-precision analytic
integration can leave a residue on either side of zero. Rejecting every
negative result breaks valid reduced models, while accepting an arbitrary
negative diagonal hides loss of positive semidefiniteness.

Current `main` uses the image-free direct self-energy and image count to form a
machine-epsilon-scale cancellation bound. A repeated-image fixture still
rejects a genuinely negative diagonal. Two fixed-plane fixtures cover positive
and negative roundoff.

The negative fixture was re-evaluated on the 4.95.66 native build:

| Quantity | Value |
|---|---:|
| Folded annihilated entry | -2.7755575615628914e-17 |
| Positive companion entry | 1.9252635826150971e-1 |
| Absolute ratio | 1.4416506844184082e-16 |

The earlier review's `-5e-7` wording is not retained because it contradicts
both this fixture and the accompanying `1e-16`-scale statement. Any larger
mesh-specific residue must be preserved as a result artifact before it is used
as evidence.

The diagonal-normalization fixture is green, but it does not close the public
field contract. On the 2026-08-29 LAB build, the following tests failed both in
the combined validation process and when run alone:

| Test | Measured error | Required limit |
|---|---:|---:|
| HEX `rad.Fld`, one reflected cell | 2.024822054511977e-14 relative | 2.220446049250313e-15 |
| HEX `FieldFromSolution`, multicell reflection | 4.931889189540598e-14 relative | 2.220446049250313e-15 |
| curved TET BDM2 IMA field | 2.393107934040017e-15 componentwise | 2.220446049250313e-15 |

These are deterministic arithmetic-order discrepancies, not the previously
reported same-process state contamination: each result reproduced in a fresh
pytest process. They are also not permission to relax the gate. A matching
full/reduced mesh must share a reflection-invariant operator, solve reduction,
and field accumulation order closely enough to satisfy the declared roundoff
contract.

## 5. Corrections retained from the original review

1. The non-convergence contract was already present in `src/radia/vim/_solve.py`.
   The missing behavior was immediate `p^T A p` breakdown detection.
2. A zero folded diagonal is legitimate for an antisymmetric fixed-plane
   charge. The gate must distinguish annihilation roundoff from an invalid
   negative self-energy.
3. The affine-HEX conclusion must be scoped to the measured Cubit
   decomposition and mesh family, not stated as a mathematical impossibility.
4. The chamfer taper, not only the erroneous spline surface, creates the
   non-affine cells.
5. A comparison script used a silent fallback and reported 0.000000% by
   evaluating only the coil field. That result is invalid and must not be
   reused.
6. The Cubit MCP disconnection was ultimately traced to an expired Cubit
   license, not an `export netgen` crash.

## 6. Remaining work and acceptance criteria

| Priority | Item | Acceptance criterion |
|---|---|---|
| P1 | Full-versus-IMA `rad.Fld` roundoff | Make all three isolated failures in section 4 pass below `10 eps` without weakening tolerances. Compare solved coefficient vectors before debugging source evaluation, then align directed block symmetrization and full/reduced field summation order. |
| P1 | Remove production RT0 | Delete the 3D order-0 entries from `hdiv_capabilities`, remove the order-0 `DemagOperator` production path and dedicated order-0 tests/docs, and keep public `Solve`, operator, field, and MATLAB inventory consistently BDM1/BDM2. |
| P1 | Fine-TET operator indefiniteness | Commit the failing mesh/configuration and result JSON; materialize the relevant Gram/operator block; locate a negative mode; compare it with dense analytic assembly or an independent NGSolve weak-form route; add a focused regression. |
| P1 | Reproducible ESRF comparison | Promote `src/radia/esrf_examples.py` into a tracked API or validation driver, commit configs and result JSON, use the corrected Cubit CAD, and run both routes with one B-H interpolation contract. |
| P2 | Mapped HEX BDM2 material solve | Build one composite mapped charge representation that preserves volume/surface cancellation, then require spectrum, linear/nonlinear solve, IMA, field, curved, and shape-derivative gates before removing the fail-loud guard. |
| P2 | Image-aware field acceleration | Design grouping that is invariant under explicit reflection and reduced IMA representation; prove `<10 eps` direct parity before enabling tree/H-matrix evaluation for image-bearing field maps. |
| P2 | Vector-potential topology coverage | Add exact/native BDM2 and HEX/WEDGE/curved source representations only with independent NGSolve mapped-volume convergence and A/B route checks. Keep the current quadrature construction explicit until then. |
| P2 | MATLAB method parity | Preserve the native field/EnergyStop handles, but do not claim native MATLAB HDiv solve parity while `vim-public` is classified as Python fallback. Promote stable numeric/artifact boundaries with MATLAB regression tests. |
| P2 | B-H contract | Select one interpolation rule, test value and differential-permeability parity at nodes and between nodes, then rerun the C-yoke comparison. |
| P2 | Configuration provenance | Classify all 14 `RADIA_HDIV_*` variables; keep fault injection/test telemetry private, expose supported tuning through `SolverConfig`, and serialize resolved values into result artifacts. |
| P3 | Same-material interfaces | Teach the exporter to omit internal interfaces only when material and orientation contracts permit it; retain an A/B regression before changing production meshes. |
| P3 | Class ownership | Continue decomposition only along measured ownership boundaries; do not replace the old branch cascade with another flag registry. |

## 7. Focused verification

This revision was verified on LAB with the native module loaded from this
worktree:

- `Build.ps1 -Verbose`: PASS after the C++ review hardening;
- 71 focused production tests: PASS in 33.05 s;
- loop-free, symmetry-loop, PSD, high-order TET, linear recoil, and irreversible
  EnergyStop validation: PASS;
- NGSolve HDiv pyramid tripwire: expected xfail;
- 3 IMA field roundoff tests: FAIL, including isolated-process reruns, with the
  measurements recorded in section 4;
- batched/block-PCG true-residual and constrained-H-matrix checks: 2 tests PASS.

The focused commands were:

```powershell
python -m pytest -q `
  tests/test_hdiv_vim_capabilities.py `
  tests/test_hdiv_chargegram_build_safety.py `
  tests/test_hdiv_same_process_determinism.py `
  tests/test_hdiv_vim_2d_orders.py `
  tests/test_hdiv_vim_2d_ima.py `
  tests/test_hdiv_vim_hex_wedge_rt2.py `
  tests/test_hdiv_field_evaluator.py `
  tests/test_hdiv_vim_magnetization_source.py `
  tests/test_hdiv_vim_linear_recoil.py `
  tests/test_hdiv_vim_energy_stop.py `
  tests/test_hdiv_vim_hysteresis_rt2.py `
  tests/test_hdiv_vim_coupled.py

python -m pytest -q `
  validation_test/feec/test_hdiv_radfld_contract.py `
  validation_test/feec/test_hdiv_vim_highorder_cpp.py `
  validation_test/feec/test_hdiv_vim_curved_ima_roundoff.py `
  validation_test/feec/test_hdiv_vim_loop_free.py `
  validation_test/feec/test_hdiv_vim_symmetry_loops.py `
  validation_test/feec/test_hdiv_vim_psd.py `
  validation_test/feec/test_hdiv_pyramid_gate.py `
  validation_test/hysteresis/test_linear_recoil_permanent_magnet.py `
  validation_test/hysteresis/test_energy_stop_irreversible_pm.py

python -m pytest -q `
  tests/test_isochronous_topopt.py::test_native_batched_multi_rhs_is_row_major_and_true_residual `
  tests/test_topology_optimization.py::test_configured_hmatrix_prunes_inactive_principal_submatrix_exactly
```

Any timing, fine-TET, or corrected C-yoke claim belongs in
`validation_test/` and must run on an idle `mdx` or `hibino`, with the machine,
native build identity, element/geometry order, image group, material
interpolant, ACA settings, DoF, build/apply/solve timing, and result checks
recorded in JSON.
