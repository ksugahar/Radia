# HDiv-MMM implementation review

This note reviews the production HDiv-MMM charge-Gram implementation and the
ESRF example-5 C-yoke cross-route investigation that exposed several of its
failure modes.

## Review status

- Revised: 2026-08-28 against `origin/main` at `d4ba797d1`.
- HDiv hardening anchor: `87309e662` (`fix(hdiv): contain HACApK callback failures`).
- C-yoke evidence source: Claude worktree commit `cf79734b9` and its parent
  measurements. Those measurements are retained below, but the complete
  Python reproduction harness is not yet tracked on `main`.
- Scope: correctness, deterministic execution, NGSolve/TaskManager policy,
  maintainability, performance, and reproducibility.

The current source sizes are:

| File | Lines | Role |
|---|---:|---|
| `src/core/rad_hdiv_vim.cpp` | 2,284 | Analytic charge and field kernels |
| `src/core/rad_hacapk_hdiv.cpp` | 11,857 | H-matrix construction and material solves |
| `src/core/rad_hacapk_hdiv_entry.cpp` | 251 | Element-specific entry strategies |
| `src/core/rad_hacapk_hdiv.h` | 1,277 | Charge-Gram state and public C++ contract |

Wall-clock performance claims are intentionally limited to measurements made
on `mdx` or `hibino`. LAB is used for build and focused correctness tests.

## Executive assessment

The reviewed core is materially safer than the version that started this
investigation. HACApK build callbacks no longer let C++ exceptions unwind
through C or TaskManager frames; process-wide callback state is serialized;
linear and nonlinear sparse gathers are deterministic; and element entry
selection is isolated behind immutable strategies. The image-folded diagonal
gate now distinguishes a physically invalid negative self-energy from the
roundoff residue of a charge annihilated by antisymmetry.

This is not yet an unconditional production certificate for every mesh and
material route. Two release claims remain open:

1. Fine TET C-yoke runs lose positive definiteness under refinement and need a
   reproducible operator-level diagnosis.
2. The corrected CAD comparison must be rerun with one identical B-H
   interpolant on both routes before quoting FEM-versus-HDiv accuracy.

The validated HEX lane and the existing production regression suite are not
blocked by those two C-yoke research findings. The findings do block a blanket
claim that all TET refinements and the current ESRF cross-route comparison are
production-certified.

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
| P1 | Fine-TET operator indefiniteness | Commit the failing mesh/configuration and result JSON; materialize the relevant Gram/operator block; locate a negative mode; compare it with dense analytic assembly or an independent NGSolve weak-form route; add a focused regression. |
| P1 | Reproducible ESRF comparison | Promote `src/radia/esrf_examples.py` into a tracked API or validation driver, commit configs and result JSON, use the corrected Cubit CAD, and run both routes with one B-H interpolation contract. |
| P2 | B-H contract | Select one interpolation rule, test value and differential-permeability parity at nodes and between nodes, then rerun the C-yoke comparison. |
| P2 | Configuration provenance | Classify all 14 `RADIA_HDIV_*` variables; keep fault injection/test telemetry private, expose supported tuning through `SolverConfig`, and serialize resolved values into result artifacts. |
| P3 | Same-material interfaces | Teach the exporter to omit internal interfaces only when material and orientation contracts permit it; retain an A/B regression before changing production meshes. |
| P3 | Class ownership | Continue decomposition only along measured ownership boundaries; do not replace the old branch cascade with another flag registry. |

## 7. Focused verification

This revision was verified on LAB with the native module loaded from this
worktree:

- `Build.ps1 -Verbose`: PASS;
- Charge-Gram safety and same-process determinism: 8 tests PASS;
- batched/block-PCG true-residual and constrained-H-matrix checks: 2 tests PASS.

The focused commands were:

```powershell
python -m pytest -q `
  tests/test_hdiv_chargegram_build_safety.py `
  tests/test_hdiv_same_process_determinism.py

python -m pytest -q `
  tests/test_isochronous_topopt.py::test_native_batched_multi_rhs_is_row_major_and_true_residual `
  tests/test_topology_optimization.py::test_configured_hmatrix_prunes_inactive_principal_submatrix_exactly
```

Any timing, fine-TET, or corrected C-yoke claim belongs in
`validation_test/` and must run on an idle `mdx` or `hibino`, with the machine,
native build identity, element/geometry order, image group, material
interpolant, ACA settings, DoF, build/apply/solve timing, and result checks
recorded in JSON.
