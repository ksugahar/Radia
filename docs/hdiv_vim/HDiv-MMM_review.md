# HDiv-MMM implementation review

This note reviews the production HDiv-MMM implementation across the public
Python API, NGSolve assembly layer, pybind11 boundary, C++/HACApK kernels,
field evaluation, permanent-magnet models, MATLAB boundary, and focused
validation. It also reconciles the ESRF example-5 C-yoke investigation that
exposed several failure modes.

## Review status

- Revised: 2026-09-04 after the released-binary TOSCA-style mixed
  total/reduced-Omega C-yoke mesh campaign on mdx and Hibino.
- The `v4.95.71` four-level result cited below is evidence for the former
  global reduced-Omega formulation only.  It is retained as history, but it
  is **not** an accuracy certificate for the current mixed formulation.  The
  current mixed route now has a checked four-level, independent-host
  certificate on `radia 4.95.77`. It is a numerical convergence envelope, not
  an analytic-truth claim. Any later solver implementation hash requires its
  own certificate rather than inheriting this one.
- The `22fc5630e..ec57769de` HEAD increment consists of Eqnedit64 pull
  requests `#34` and `#35`. It changes no HDiv source, test, validation,
  MATLAB, or HDiv documentation path, so the earlier HDiv measurements remain
  current at this HEAD.
- HDiv hardening anchor: `87309e662` (`fix(hdiv): contain HACApK callback failures`).
- Claude review source: branch `claude/cubit-vol-main-path`, commits
  `2c939a8ed`, `4b30479b6`, and `cf79734b9`. Their measurements and useful
  regressions are reconciled below; none of those three commits is merged on
  `main` as written.
- Latest Claude implementation source: an uncommitted shared-tree patch on
  `backup/main-pre-release-20260821` at `cbc029319`, documented in local
  handover notes dated 2026-08-28. The BDM2 TET directional-moment correction
  is isolated on this review branch as `51dce89c1` and is included in this
  revision.
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

The implementation is production-capable for the declared BDM1/BDM2 primal
solve and field lanes, but it does not justify an unconditional "complete for
every HDiv case" claim. The BDM2 TET directional defect found by this review
is corrected, and the mapped/non-affine HEX BDM2 material lane is now promoted
by the cancellation-preserving composite quadrature described below. Four P1
findings remain:

1. Three explicit full-versus-IMA field regressions fail the required
   `< 10 eps` contract in isolated pytest processes.
2. Broken RT0 was reintroduced into the public capability table and
   `DemagOperator`, contradicting the BDM1/BDM2-only production policy.
3. Fine TET C-yoke runs reportedly lose positive definiteness under refinement;
   the failing mesh and driver are not tracked, so the defect is not replayable.
4. The symmetric ACA/H-matrix ChargeGram is indefinite on the ESRF #3 BDM2
   response-iron mesh.  The raw analytic Gram and a bounded exact-dense
   materialization are positive, so this is a compression defect, not a
   material, mass-Riesz, or PARDISO defect.  It blocks H-matrix material
   release until a PSD-preserving compression or certificate is implemented.
Known, fail-loud limitations are not hidden failures: mapped HEX BDM2 shape
derivatives, 3D pyramid/mixed meshes, and 2D history solves remain unsupported.
Large IMA field maps remain direct because the tree route is
disabled whenever images are present. MATLAB owns native field-evaluator and
EnergyStop MEX surfaces, but the high-level `radia.vim` solve is still an
explicit in-process Python fallback because NGSolve setup is Python-owned.

### 0.0 2026-09-02 formulation-status correction

The C-yoke review previously used the phrase "three-formulation certificate"
for a calculation whose H1 route was the former global reduced-Omega model.
That route lacks the required physical-air/Kelvin source-potential jump and
must not be presented as the current TOSCA-style mixed total/reduced-Omega
formulation.  The current implementation uses a reduced physical-air scalar
potential, a total iron/Kelvin scalar potential, and two independently
projected source traces:

```
Gamma_iron-air:   phi_total - phi_reduced =  Phi_source
Gamma_kelvin:     phi_total - phi_reduced = -Phi_source
```

The second equation is the orientation-reversing Kelvin pullback of the
source 0-form.  `project_source_interface_potential` makes each trace a
measured gate; an interface that links current or has nontrivial cohomology
fails loudly and requires an explicit cut representative or HCurl reduced-A.
The saddle system is symmetric indefinite and therefore uses PARDISO rather
than an SPD-only CG path.

On the shared exact Cubit C-yoke mesh, all three routes pass a 1% gap-core
relative-RMS B gate on Hibino.  The current v4 artifacts report a 0.45977%
maximum in the linear order-3 run and 0.16023% in the nonlinear order-2 run.
The nonlinear HDiv-MMM, HCurl reduced-A, and H1 TOSCA mixed total/reduced Omega
runs converge and take 11.27 s, 216.93 s, and 134.71 s respectively.  The
associated source-trace residuals are below the separate 5% cut/cohomology
gate.  These values establish fixed-mesh cross-formulation agreement; they
are deliberately not an analytic-truth or mesh-convergence claim.

The released-binary campaign then ran four exact Cubit/ACIS levels
(`434/542/846/1,688` iron elements and `24,134/36,208/54,886/89,454`
Kelvin-domain elements). All three nonlinear routes converged at every level.
At the finest level, the maximum pairwise parity-projected gap-core RMS
difference was `0.27714%`; the maximum discretisation uncertainty was
`0.17601%`, giving a conservative combined numerical envelope of `0.35399%`.
The independently repeated finest solve on Hibino reproduced the mdx field
within `5.25e-14` relative RMS. This closes the mixed-formulation numerical
acceptance for the exact `radia 4.95.77` implementation contract; it remains
an agreement certificate, not analytic absolute truth.

## 0. Cross-layer production review

### 0.1 Findings by priority

| ID | Priority | Finding | Evidence and required disposition |
|---|---|---|---|
| F1 | Resolved | BDM2 TET directional ChargeGram derivatives used the wrong degree-one moment order on the reviewed baseline. | `TetPotentialMomentsDirectionalUpTo1` stores degree-one moments in `z,y,x` (`PotentialMomentIndex`) order, while two consumers in `rad_hacapk_hdiv.cpp` used `mv[k+1]` as `x,y,z`. On the same 1-cell BDM2 case, the `ec57769de` baseline differs from finite differences by `1.327698e-1` for the complete Gram and `4.145764e-1` for the volume block; fifth-degree homogeneity is wrong by `3.727428e-1`. The correction in `51dce89c1`, included in this revision, reduces these to `4.019822e-9`, `2.778490e-9`, and `4.002814e-16`. |
| F2 | P1, partially resolved | The field-evaluator IMA contract is green for mapped HEX BDM2 prescribed sources; independent solve parity remains a separate numerical lane. | On the current mdx production body, prescribed full/reduced fields differ by `2.7506 eps`, below the `10 eps` limit. Independently converged mass-Riesz CG full/reduced solves differ by `3.2835e-13` in sampled field. Three legacy focused checks were rerun and remain narrowly red: single-cell HEX `2.02e-14`, multicell HEX `4.93e-14`, and curved TET BDM2 `2.3931e-15` against a `2.2204e-15` limit. Preserve the field limit and fix those paths rather than loosening their tolerances or relabeling Krylov/reduction error as evaluator error. |
| F3 | P1 | RT0 is publicly advertised again despite the BDM1/BDM2-only decision. | `_capabilities.py` exposes 3D TET/HEX order 0 and `DemagOperator` documents an order-0 broken-interface path. `HDivSolver` and field evaluation accept only orders 1 and 2. Remove the public RT0 entries/path and retain any topology-only experiment outside the production API. |
| F4 | Resolved on `v4.95.71` | The released operator completes the finer C-yoke TET lane without loss of SPD. | The 1,688-element iron mesh solves on mdx and hibino, all three nonlinear routes converge, and the final three mesh levels pass the contraction/order gate. The older untracked `p^T A p < 0` report is not used as current evidence. |
| F5 | Resolved for the released `4.95.77` implementation contract | The C-type comparison is the TOSCA-style mixed total/reduced-Omega route with two required source-trace jumps, and its three-route nonlinear BDM2 mesh certificate passes. | `validation_test/c_type_three_engine/` owns the exact Cubit/ACIS mesh, shared CoilBuilder, PCHIP B(H) law, Kelvin contract, checkpoints, and portable JSON gates. Four levels on mdx plus the independent Hibino finest replay yield 0.27714% maximum finest pairwise gap-core RMS, a 0.35399% combined numerical envelope, and `5.25e-14` cross-host RMS. The old global-Omega four-level artifact is historical only. A later implementation hash must rerun this campaign before it can make the same claim. |
| F6 | Resolved for primal solve/field; derivative open | Mapped/non-affine HEX BDM2 is a production material lane. | Complete-host tensor source rules preserve smooth-pair charge cancellation; reflection-invariant whole-host Duffy rules handle self and adjacent pairs. On mdx the 756-DoF q9/q12 operator has spectrum `[-8.53e-16, 0.999899]`, linear/nonlinear solves converge, and its material response differs from q10/q16 by `5.28e-4` in mass norm. q10/q16 differs from q11/q20 by `3.94e-4`. An independent Cubit 2025.12 Curve(2) four-HEX gate also passes linear/nonlinear solve and field checks. Shape derivatives fail loudly until the composite rule is differentiated. |
| F7 | P2 | IMA disables tree acceleration for field maps. | `HDivFieldEvaluator::AlgorithmFor` returns `Direct` whenever images exist. This protects full/reduced roundoff parity, but large IMA observation maps cannot use the otherwise guarded treecode. Any image-aware acceleration needs a common full/reduced grouping and the F2 contract first. |
| F8 | P2 | Exact vector-potential evaluation is narrower than H-field evaluation. | Exact `A` uses straight TET BDM1 equivalent currents. BDM2, curved, HEX, and WEDGE use NGSolve-mapped quadrature clouds assembled in Python. This is valid as an explicit converged quadrature route, not an all-topology exact/native claim. |
| F9 | P2 | MATLAB parity is partial at the method level. | `HDivFieldEvaluator` and `EnergyStopMaterial` have native checked MEX handles. The parity manifest classifies `vim/__init__.py` as `python-fallback`; complete solve/mesh/form orchestration is not a native MATLAB API. |
| F10 | P3 | Configuration provenance and class ownership remain broad. | Fourteen `RADIA_HDIV_*` variables remain, and `RadHACApKChargeGram` still owns entry, build, cache, solve, derivative, field, and diagnostics. Keep diagnostic switches out of release claims and split only at measured ownership boundaries. |
| F11 | P2, integration open | Kelvin and FFAG rotational periodic identifications can coexist only as distinct identification classes. | The FFAG cyclic helper accepts an explicit `identnr`, while Kelvin identification refuses to treat an unrelated cyclic pair as Kelvin-ready. Focused tests cover separate IDs and the non-poisoning predicate. An end-to-end sector FFAG magnet solve with paired periodic faces, Kelvin exterior, and an independent FEM/HDiv B comparison remains the next application validation. |

### 0.2 Production capability matrix

| Dimension/topology | HDiv space | Geometry | Material/history status | Field status |
|---|---|---|---|---|
| 3D TET | BDM1, BDM2 | affine or P2 curved | linear, nonlinear, IMA, recoil PM, Play, EnergyStop | persistent native direct/tree; curved direct leaves retained exactly; exact vector `A` only for straight BDM1 |
| 3D HEX | BDM1 | mapped/affine, geometry order 1 or 2 | linear/nonlinear/IMA/history | affine polynomial decomposition; mapped/curved source cloud |
| 3D HEX | BDM2 | affine or mapped geometry order 1 or 2 | linear, nonlinear, IMA, recoil PM, Play, EnergyStop; mapped shape derivative fails loudly | persistent native evaluator from configured source representation |
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

### 0.4 Latest Claude BDM2 directional implementation

The latest Claude HDiv implementation originated outside a branch commit. It
was a small part of the much larger dirty shared tree at
`S:\Radia\01_GitHub`, based on `backup/main-pre-release-20260821` at
`cbc029319`. Commit `51dce89c1` isolates only the patch identified in
the local handover notes; no other co-located WIP was copied.

The correction is mathematically and structurally appropriate:

1. `MomentIndex3` is made `constexpr`.
2. `MomentIndex3Linear` derives the physical x, y, and z degree-one slots from
   that canonical indexer instead of hard-coding `{3,2,1}`.
3. Both BDM2 TET volume-moment contractions use the derived map instead of
   `mv[k+1]` and `dm[k+1]`.
4. The misleading `[1,x,y,z]` comment is replaced with the actual
   `PotentialMomentIndex` storage contract.

An independent LAB comparison used the same mesh, deformation field,
quadrature, charge basis, and finite-difference step for clean `main` and the
Claude candidate:

| Quantity | `origin/main` source | Claude candidate | Reference scale |
|---|---:|---:|---:|
| complete Gram derivative relative error | 1.3276980331e-1 | 4.0198217292e-9 | 1.1070212280e-1 norm |
| volume self-block relative error | 4.1457643916e-1 | 2.7784903499e-9 | 2.5207637398e-3 norm |
| fifth-degree homogeneity relative error | 3.7274276552e-1 | 4.0028143694e-16 | 4.0072631527e-2 norm |
| rigid-translation maximum derivative | 2.1996594761e-19 | 2.1996594761e-19 | exact zero |

The candidate test
`test_native_tet_directional_derivative_matches_fd_at_bdm2` covers the right
concepts and its tolerances reject the measured unfixed errors by a wide
margin. Cross-worktree execution needs care: `tests/conftest.py` inserts its
own repository's `src` at the front of `sys.path`, so pointing pytest at the
shared-tree test from another worktree still loads the shared-tree native
module. The comparison above therefore used one standalone driver with an
explicit module path for each build. Commit `51dce89c1` then transplanted the
patch and test together, added aggregate relative-error gates, rebuilt the
native module in this worktree, and passed the new BDM2 test plus the adjacent
BDM1 self-block and complete-Gram/Piola regressions in 4.14 s.

The handover also reported an access violation in the zero-coupling candidate
Schur test. That report is stale against current `main`: after a clean
`Build.ps1 -RadiaOnly -Rebuild`,
`test_native_candidate_schur_reports_zero_coupling_rank_and_stable_iters`
passes in 2.47 s. No crash finding is carried forward without a reproducer on
the current native build.

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

The historical figures below came from the Claude C-yoke campaign. The
replacement self-replaying route is now `validation_test/c_type_three_engine/`:
Cubit/ACIS is the CAD authority, HDiv receives an iron-only `.vol`, and both
FEM formulations receive one periodic spherical Kelvin `.vol`. No finite
outer air box is part of the comparison.

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

### 3.2 Curved HEX production path

The original statement "a chamfered pole cannot be meshed into affine hexes"
was too absolute. What was measured is narrower: the current Cubit
decomposition and mesh family leave cells crossing the 45-degree chamfer as
non-affine trapezoidal prisms. Their mapping residual was 0.745, 0.666, and
0.405 at 25, 12.5, and 6.25 mm, respectively, against the 1e-10 affine gate.
Removing the chamfer reduced the residual to 2.2e-14, but that changes the
benchmark and is not a production remedy. The remedy is the mapped BDM2
composite charge representation now used by `Solve` and `ChargeGram`.

For smooth host pairs, the C++ kernel uses a complete Q2 tensor source rule so
the volume and surface modes share one host representation. For self and
adjacent interactions it sweeps the complete host to six faces and each face
to four edges with target-anchored, reflection-invariant Duffy rules. This
removes the fixed sub-tet/sub-triangle diagonals that broke the large
volume/surface cancellation on non-affine cells. The operator remains
`B.T G B`, so loop-free nullspaces and symmetry are preserved by construction.
The q9/q12 default is materialized once and reused during parallel H-matrix
fill.

The tracked mdx quadrature-convergence evidence uses a reflection-symmetric
eight-cell non-affine trilinear body (mesh curve order 1, 756 BDM2 DoF). Its
q9/q12 build takes 46.74 s and has spectrum `[-8.53e-16, 0.999899]`;
q10/q16 takes 140.60 s and q11/q20 takes 351.50 s.
The material mass-norm differences are `5.28e-4` for q9/q12 to q10/q16 and
`3.94e-4` for q10/q16 to q11/q20. Linear, energy-Newton, and IMA solves pass.
This certifies the primal material/field path; mapped BDM2 topology derivatives
remain fail-loud pending a differentiated composite rule.

True curved geometry is covered separately by a Cubit 2025.12 generated
four-HEX cylinder. `check-vol` records curve order 2, four of four non-affine
cells, minimum scaled Jacobian `0.33473`, and volume error `-0.2197%`. On mdx,
the 396-DoF BDM2 linear solve converges in 19 mass-Riesz CG iterations, the
equivalent nonlinear Energy-Newton path converges, and prescribed-source
`rad.Fld` differs from an independent NGSolve boundary integral by at most
`1.06e-10` relative. This closes primal solve/field production support for
Curve(2) HEX; mapped BDM2 shape derivatives remain deliberately fail-loud.

A united C-yoke solid did not hex-mesh with the tested `auto`, `sweep`,
`webcut_cyl_auto`, or `polyhedron` routes. Decomposition creates
same-material internal interfaces. The Cubit 2025.12 exporter now removes
their surface elements and unused face descriptors while retaining the shared
volume nodes, edges, and mesh-support curves. The generated C-yoke `.vol`
therefore contains only physical `iron_air_interface`, `kelvin_int`, and
`kelvin_ext` boundaries.

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

### 3.4 Canonical Kelvin three-engine result

The canonical C-yoke route intentionally has two meshes derived from the same
Cubit/ACIS solid:

- HDiv-MMM: 542 TET iron elements, no air mesh, Coulomb ChargeGram open
  boundary;
- reduced-A and Omega-reduced-Omega: 36,208 TET elements over iron, a locally
  refined physical-air sphere, and a translated Kelvin sphere.

The Kelvin mesh has 462 periodic point pairs. Their fitted translation is
`(0.6600000000000019, 1.87e-17, 0) m`, with maximum pair error
`2.00e-15 m`. `Periodic(H1)` slaves 462 free DOFs, and the functional
`kelvin_ext`/`kelvin_int` trace-norm ratio is `0.9999999999999987`. Both
meshes have zero missing reflected vertices and zero
missing reflected elements. The 10 mm air gap contains 2,082 elements and a
maximum z span of 5 mm; this is a local gap-resolution gate, not a global air
box.

With order-2 HDiv and FEM spaces, `gram_eps=1e-14`, and linear `mu_r=1000`,
the maximum parity-projected gap-core pairwise relative RMS is 0.12113%.
Off-plane reflection errors are `1.95e-10` for HDiv, `4.67e-10` for
reduced-A, and `1.95e-10` for Omega. The three LAB runtimes are 14.01 s,
12.00 s, and 2.60 s respectively; these are correctness timings, not release
performance claims.

The first nonlinear order-1 smoke exposed two missing production contracts:
the reduced-A Picard path used an undefined gauge coefficient, and the Omega
Picard path constructed raw `H1` instead of the periodic Kelvin H1 factory.
Both are corrected and regression-locked. After correction all three engines
converge and retain `1.6e-10` to `2.2e-10` reflection error. HDiv and reduced-A
agree to 0.49342% in the gap core, while Omega remains 5.85068% away from
HDiv. Replacing Omega's linear table interpolation with the same monotone PCHIP
and vacuum-slope continuation used by HDiv changes that result only to 5.84969%.
The interpolation mismatch was real but was not the source of the order-1
field discrepancy.

The direct nonlinear order-2 primary comparison closes the discrepancy. Both
engines converged with the shared PCHIP material law and `gram_eps=1e-14`; the
parity-projected gap-core relative RMS is 0.18032%, with a maximum vector
difference of `8.2457e-4 T`. HDiv used 10,860 DoF, five Newton iterations,
1,638 inner linear iterations. Omega used 50,322 DoF and 17 Picard iterations.
Their off-plane reflection errors are `1.85e-10` and about `1.83e-10`,
respectively. This establishes the accuracy comparison and identifies the
order-1 result as a discretization failure.

The exact `v4.95.70` PyPI wheel was then run three times without checkpoint
reuse on each idle 38-core Xeon Platinum 8368 host. Median HDiv/Omega runtimes
were 12.09/42.59 s on mdx and 11.31/43.28 s on hibino. Thus HDiv was 3.52x and
3.83x faster in this fixed order-2 nonlinear comparison while using 4.63x fewer
DoFs. Every run produced the same 0.1803201266% gap-core discrepancy and the
same nonlinear iteration counts. The raw artifacts and their hashes are
indexed by
`validation_test/c_type_three_engine/results/mdx_hibino_20260830_nonlinear_order2_summary.json`.

The subsequent `v4.95.71` four-level campaign belongs to the historical global
reduced-Omega formulation. Its refinement arithmetic is retained for
diagnostic archaeology only: the formulation lacks the physical-air/Kelvin
source-potential jump and is not an acceptance route. The current accepted
campaign is the `4.95.77` mixed-formulation evidence reported in section 0.0
and stored as `c_type_20260903_nonlinear_bdm2_mesh_convergence_certificate.json`.

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
| HEX `rad.Fld`, one reflected cell | 2.014001326215872e-14 relative | 2.220446049250313e-15 |
| HEX `FieldFromSolution`, multicell reflection | 4.931889189540598e-14 relative | 2.220446049250313e-15 |
| curved TET BDM2 IMA field | 2.393107934040017e-15 componentwise | 2.220446049250313e-15 |

These are deterministic arithmetic-order discrepancies, not the previously
reported same-process state contamination: each result reproduced in a fresh
pytest process. They are also not permission to relax the gate. A matching
full/reduced mesh must share a reflection-invariant operator, solve reduction,
and field accumulation order closely enough to satisfy the declared roundoff
contract.

The new mapped HEX BDM2 production validation separates two error sources. A
prescribed symmetric magnetization, which bypasses iterative material solves,
gives maximum pointwise vector-relative error `8.0078e-16` on mdx; the global
componentwise maximum absolute error is `2.7506 eps` of the field scale.
The persistent field evaluator therefore satisfies the `<10 eps` contract on
that topology. Two independent mass-Riesz CG material solves give
`3.2835e-13` field difference and `2.8799e-13` average-magnetization
difference. Those values are within the production solve gate (`1e-10`) but
are not field-evaluator roundoff evidence. Future strict solve parity must
align the reduced/full linear solve or use an independently justified common
solution representation; it must not weaken the field gate.

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

### 5.1 ESRF #3 material-operator PSD boundary (2026-09-04)

The ESRF #3 hybrid-undulator response mesh is the decisive counterexample to
the assumption that a symmetrized ACA matrix is enough for CG.  On the 144-cell
curved Q2 HEX iron mesh (14,040 HDiv DoF, 8,640 charge DoF), the first
mass-Riesz PCG direction has `p^T G_H p = -110486.500638061`, while the fixed
order raw analytic oracle gives `p^T G_raw p = +78453.287592292`.  The
compressed material solve therefore fails at iteration zero; increasing the
iteration limit or changing a CG tolerance cannot repair it.

An explicit, memory-bounded `gram_backend="exact-dense"` diagnostic was added
for medium meshes.  It materializes the normalized physical Gram directly from
the analytic entry oracle, requires a caller-provided cap, and preserves the
same NGSolve charge map, mass matrix, and configured constraint semantics.  At
a 1 GiB cap, #3 used 569.53125 MiB and took 95.50 s to materialize.  It gives
`p^T G_dense p = +78453.287528763`, agrees with the raw oracle to `8.1e-10`
relative, and the linear material recovery converges in 58 iterations to
`9.59e-10` true residual with `4.48e-8` coefficient recovery error.

The public `vim.Solve` path was also profiled on the same response mesh with
`gram_backend="exact-dense"`, order 2, curve order 2, eight threads, and the
same explicit 1 GiB cap.  Its zero-load setup completed in 149.71 s (149.16 s
charge-Gram stage), used no H-matrix statistics, and required zero linear
iterations as expected for a zero source.  This is a public-path smoke and a
cost boundary, not a scalable timing claim.

This validates the FE/mass/PARDISO path and provides a safe three-engine
medium-mesh reference.  It is not a scalable fallback and must not be selected
implicitly.  The default H-matrix backend remains release-blocked for material
solves until it has a genuine PSD-preserving construction or an acceptance
certificate that rejects this failure before CG starts.

## 6. Remaining work and acceptance criteria

| Priority | Item | Acceptance criterion |
|---|---|---|
| P1 | Full-versus-IMA `rad.Fld` roundoff | Make all three isolated failures in section 4 pass below `10 eps` without weakening tolerances. Compare solved coefficient vectors before debugging source evaluation, then align directed block symmetrization and full/reduced field summation order. |
| P1 | Remove production RT0 | Delete the 3D order-0 entries from `hdiv_capabilities`, remove the order-0 `DemagOperator` production path and dedicated order-0 tests/docs, and keep public `Solve`, operator, field, and MATLAB inventory consistently BDM1/BDM2. |
| P1 | Nonlinear C-yoke memory evidence | Four-level accuracy and repeated timing are closed on mdx and hibino for `v4.95.71`. Add measured process peak memory to a future scaling campaign before making a memory-efficiency claim. reduced-A remains an independent third-formulation audit rather than the primary production route. |
| P1 | H-matrix material PSD certificate | On ESRF #3, prove the compressed Gram is PSD before a material-CG solve, or replace independent ACA blocks with a PSD-preserving compression. The H-matrix backend is release-blocked until the counterexample in section 5.1 passes without exact-dense materialization. |
| Resolved for primal path | Mapped HEX BDM2 material solve | The composite mapped charge representation passes spectrum, linear/nonlinear solve, IMA, field, and quadrature-convergence gates on mdx. |
| P2 | Mapped HEX BDM2 shape derivative | Differentiate the same complete-host tensor and whole-host Duffy representation, then lock it against finite differences before enabling topology optimization. The current API fails loudly. |
| P2 | Image-aware field acceleration | Design grouping that is invariant under explicit reflection and reduced IMA representation; prove `<10 eps` direct parity before enabling tree/H-matrix evaluation for image-bearing field maps. |
| P2 | Vector-potential topology coverage | Add exact/native BDM2 and HEX/WEDGE/curved source representations only with independent NGSolve mapped-volume convergence and A/B route checks. Keep the current quadrature construction explicit until then. |
| P2 | MATLAB method parity | Preserve the native field/EnergyStop handles, but do not claim native MATLAB HDiv solve parity while `vim-public` is classified as Python fallback. Promote stable numeric/artifact boundaries with MATLAB regression tests. |
| Resolved | reduced-A B-H contract | reduced-A now inverts the shared monotone PCHIP B(H) law by checked scalar root solves, uses the same vacuum-slope continuation, and passes the four-level three-formulation certificate. |
| P2 | Configuration provenance | Classify all 14 `RADIA_HDIV_*` variables; keep fault injection/test telemetry private, expose supported tuning through `SolverConfig`, and serialize resolved values into result artifacts. |
| Resolved | Same-material interfaces | The exporter removes only `DomainIn == DomainOut > 0` seams, remaps retained descriptors, and the regenerated C-yoke mesh passes strict labels, adjacency, exact reflection, and Kelvin identification gates. |
| Resolved | Fine-TET operator indefiniteness | The released `v4.95.71` operator solves the 1,688-element finest C-yoke iron mesh on both mdx and hibino and the final three levels satisfy the contraction/order gate. |
| P3 | Class ownership | Continue decomposition only along measured ownership boundaries; do not replace the old branch cascade with another flag registry. |

## 7. Focused verification

This revision was verified on LAB with the native module loaded from this
worktree, then repeated with the exact `v4.95.70` timing wheel and the final
`v4.95.71` four-level certificate wheel on mdx and hibino:

- `Build.ps1 -Verbose`: PASS after the C++ review hardening;
- latest `origin/main` clean native rebuild with
  `Build.ps1 -RadiaOnly -Rebuild`: PASS;
- 71 focused production tests: PASS in 33.05 s;
- latest-HEAD smoke (`ec57769de`): 15 focused build-safety, deterministic,
  capability, and field-evaluator tests PASS in 3.78 s;
- BDM2 TET directional finite-difference comparison: the `ec57769de` baseline
  is wrong by 13.28% for the complete Gram and 41.46% for the volume block;
  commit `51dce89c1`, included here, reduces them to `4.02e-9` and `2.78e-9`
  relative error;
- isolated support commit `51dce89c1`: new BDM2 regression plus adjacent BDM1
  self-block and complete-Gram/Piola regressions, 3 tests PASS in 4.14 s;
- complete topology-optimization regression after the native rebuild:
  160 tests PASS in 56.53 s;
- current-main zero-coupling candidate Schur regression: PASS in 2.47 s after
  the clean rebuild; the older access-violation report is not reproduced;
- loop-free, symmetry-loop, PSD, high-order TET, linear recoil, and irreversible
  EnergyStop validation: PASS;
- NGSolve HDiv pyramid tripwire: expected xfail;
- 3 IMA field roundoff tests: FAIL, including isolated-process reruns, with the
  measurements recorded in section 4;
- batched/block-PCG true-residual and constrained-H-matrix checks: 2 tests PASS.
- exact-dense ChargeGram entry, public `DemagOperator`/`vim.Solve`, and
  configured-principal-submatrix semantics: 9 tests PASS in 5.87 s on LAB;
- ESRF #3 exact-dense material diagnostic: PASS under its explicit 1 GiB cap;
  the parallel H-matrix run fails at PCG iteration zero as recorded in section
  5.1 and remains a release blocker rather than a passing validation result.
- canonical Cubit C-yoke mesh: PASS with exact reflected topology, 462 Kelvin
  point pairs, and `2.00e-15 m` maximum translation error;
- order-2 linear Kelvin three-engine comparison: PASS, 0.12113% maximum
  gap-core pairwise relative RMS;
- order-1 nonlinear Kelvin smoke: all engines converged; FAIL accuracy at
  5.85068% because Omega remains outside the 3% gate.
- order-1 nonlinear primary pair with the shared PCHIP law: converged; FAIL
  accuracy at 5.84969%, ruling out interpolation choice as the material cause;
- order-2 nonlinear primary pair with the shared PCHIP law: PASS at 0.18032%
  gap-core relative RMS on all six remote runs; median HDiv/Omega timing was
  12.09/42.59 s on mdx and 11.31/43.28 s on hibino, with 10,860/50,322 DoF.
- historical four-level order-2 nonlinear global-Omega certificate on
  `v4.95.71`: PASS for that retired formulation only. It is not evidence for
  the current TOSCA mixed route and must not be used in release material.
- current v4 TOSCA mixed nonlinear BDM2 mesh certificate on `radia 4.95.77`:
  PASS. Four Cubit levels converge for exactly `hdiv_mmm`, `reduced_a`, and
  `mixed_total_reduced_omega`. The finest maximum pairwise gap-core RMS is
  0.27714%, the combined numerical envelope is 0.35399%, and the mdx/Hibino
  replay is `5.25e-14` relative RMS. A global reduced-Omega calculation is
  historical evidence only.
- release-qud: PASS for `radia 4.95.71`; package versions and production file
  hashes agree across LAB, the 100-machine, mdx, and hibino.

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
`validation_test/` and must run on hibino first, or on mdx only when hibino is
unavailable and the mdx CI queue is idle, with the machine, native build
identity, element/geometry order, image group, material
interpolant, ACA settings, DoF, build/apply/solve timing, and result checks
recorded in JSON.
