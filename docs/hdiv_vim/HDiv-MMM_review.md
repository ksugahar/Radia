# HDiv-MMM implementation review

This note reviews the production HDiv-MMM implementation across the public
Python API, NGSolve assembly layer, pybind11 boundary, C++/HACApK kernels,
field evaluation, permanent-magnet models, MATLAB boundary, and focused
validation. It also reconciles the ESRF example-5 C-yoke investigation that
exposed several failure modes.

## Review status

- Revised: 2026-09-04 after the released-binary TOSCA-style mixed
  total/reduced-Omega C-yoke mesh campaign on mdx and Hibino, and the
  symmetric ChargeGram diagonal-leaf repair for ESRF Example #3.
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
4. The symmetric ACA/H-matrix ChargeGram previously made diagonal,
   zero-width clusters low rank on the ESRF #3 BDM2 response-iron mesh.
   The symmetric diagonal-leaf repair removes that construction defect and a
   focused regression locks it. The full nonlinear three-engine #3 run remains
   the acceptance gate; a single positive quadratic direction is not a global
   PSD certificate.
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
| F12 | Implementation resolved; remote application rerun open | The reduced-A Kelvin failure was a gauge-conditioning defect, not a periodic-identification or source-pullback defect. | The order-2 HCurl interface slaves all edge and face trace DOFs. On the real Example #6 mesh, `eps=1e-10` makes both Periodic and non-Periodic BDDC local inverses non-finite; strengthening only the Kelvin region does not repair them. A uniform `1e-6` gauge gives a finite 705,838-DoF solve in 28 CG iterations with `5.31e-8` true relative residual. `VectorPotentialSolver` now defaults every physical and Kelvin block to `1e-6`, permits an explicit independent `kelvin_eps`, enables BDDC for large Kelvin systems, and rejects non-finite or inaccurate linear solves. The released-wheel #6/#7 three-formulation reruns remain required. |

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

The audited source contains 16 supported `RADIA_HDIV_*` environment controls.
They are now exhaustively classified and guarded by
`tests/test_hdiv_environment_policy.py`; adding an unclassified control fails
the fast test lane. The obsolete `RADIA_HDIV_HEX_CACHE_STATS` alias was removed.

| Class | Controls | Contract |
|---|---|---|
| Diagnostic counters | `RADIA_HDIV_BLOCK_CACHE_STATS`, `RADIA_HDIV_HMATVEC_STATS` | Opt-in instrumentation; invalidates production timing claims. |
| Failure injection | `RADIA_HDIV_TEST_FAIL_FILL_AFTER` | Test-only build failure; no public API and no successful result artifact. |
| Performance/cache A/B | `RADIA_HDIV_HEX_BLOCK_CACHE_LIMIT`, `RADIA_HDIV_WEDGE_TRANS_CACHE`, `RADIA_HDIV_DISABLE_TRANS_CACHE`, `RADIA_HDIV_DISABLE_CONGRUENT_CACHE` | Diagnostic benchmark paths; effective values are copied into `hmat_stats`. |
| Numerical/path A/B | `RADIA_HDIV_CURVED_DIRECT`, `RADIA_HDIV_HEX_FAR_ONESIDED`, `RADIA_HDIV_WEDGE_FAR_ONESIDED`, `RADIA_HDIV_HEX_DISTORTED_FAR_FACTOR`, `RADIA_HDIV_HO_FAR_ONESIDED`, `RADIA_HDIV_DISABLE_HO_ANALYTIC_BLOCK`, `RADIA_HDIV_DISABLE_HO_IMAGE_BLOCK`, `RADIA_HDIV_DISABLE_HO_IMAGE_FAR`, `RADIA_HDIV_HEX_CLUSTER_RADIUS`, `RADIA_HDIV_HEX_PAIR_DUFFY`, `RADIA_HDIV_HEX_GLPAIR_W_N` | Non-production comparison paths; `hmat_stats.nonproduction_numerical_override_active` and `release_claim_eligible` make them fail-loud in provenance. |
| Preconditioner A/B | `RADIA_HDIV_AUTO_JACOBI_TET_NFACE` | Measurement-only auto-policy threshold; the resolved threshold and branch are recorded in `preconditioner_policy`. |

The native `stats()` surface records every effective C++ cache, quadrature,
image-block, profiling, and MKL-thread choice. `vim.Solve` already embeds that
dictionary as `hmat_stats`, while the resolved Python preconditioner decision
is carried separately. Thus a successful solve no longer leaves a hidden
release-relevant HDiv path unreported. A/B validation must use subprocesses
because the native hot-path settings are intentionally read once per process.

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
implicitly.

### 5.2 Symmetric diagonal-leaf repair (2026-09-04)

The failure was not caused by IMA: the #3 asset has no image reduction, and
independent fourfold/sixfold cyclic IMA checks remain within `8.43e-5` relative
field RMS. The failing H-matrix contained same-cluster diagonal ACA leaves.
For example, a size-18/rank-3 leaf contributed `1824.028974949` to the
quadratic form where its raw analytic self block contributed
`41092.606726508`. In symmetric storage an off-diagonal leaf is reflected,
but a diagonal leaf is applied once; unconstrained ACA there is neither
symmetric nor positive semidefinite.

`cHACApK_count_lntmx` and `cHACApK_generate_leafmtx` now keep each
same-cluster diagonal leaf dense only when Radia requests symmetric fill.
Off-diagonal leaves remain ACA-compressed. The native regression uses repeated
charge centres and asserts that there are no low-rank diagonal leaves and that
the materialized H-matrix quadratic matches the raw analytic quadratic.

On the fixed first PCG direction of the #3 mesh, the repaired `leaf=16`
operator gives `+79506.431938098`, replacing the former negative
`-110486.500638061`. The calibrated production default is `leaf=64`: it gives
`+78435.277824058` against `+78453.287592292` from the raw oracle, a `0.023%`
gap, with 1,084 low-rank and 1,458 dense leaves (94.65 MB, 16.62% of dense
storage). The repair fixes the identified construction defect; remote
nonlinear three-engine evidence is still required before closing Example #3.

### 5.3 ESRF coil-yoke three-engine readiness (2026-09-04)

ESRF Examples 6 and 7 are coil-driven nonlinear quadrupoles.  The validation
must use one mesh-free `CoilBuilder` solid-current object tree for every
formulation: HDiv-MMM receives the checked iron-only Cubit response mesh, and
reduced-A plus mixed total/reduced-Omega receive a separate conforming
iron/physical-air/Kelvin mesh made from the same iron STEP authority.  Coils
are deliberately not FEM-meshed.  A finite outer-air box is forbidden.

Both cases now have an executable, checkpointed three-engine runner and a
fast preflight that verifies the source closure, source field finiteness,
observation-point containment, material and boundary labels, Kelvin
identifications, and the mesh hashes that bind checkpoints to their inputs.
The runner compares the 45-point symmetric gap stencil and treats all three
nonlinear convergences plus the configured pairwise core RMS limit as its
acceptance gate.

| Case | HDiv response mesh | Independent FEM mesh | Readiness result |
|---|---:|---:|---|
| #6, quadrupole | existing curved-Q2 Cubit iron HEX | 134,686 curved-Q2 elements, 1,393 Kelvin identifications | Preflight PASS; full nonlinear three-engine run pending the released H-matrix repair. |
| #7, ESRF storage-ring quadrupole | 367,845 curved-Q2 Cubit iron TET | 839,376 curved-Q2 elements, 1,409 Kelvin identifications | Preflight PASS; full nonlinear three-engine run pending the released H-matrix repair. |

Case #7 is explicitly TET at present.  Cubit 2025.12 cannot map or submap
the full 400 mm end-chamfered yoke without a further topology partition.  The
policy therefore records `tetmesh` and `TET` for both HDiv and FEM iron rather
than silently exporting an empty mesh.  Its iron-only Cubit mesh passes the
strict label, curved-mapping, and CAD-volume checks: all 23,542,080 Jacobian
samples have a consistent orientation, the minimum scaled Jacobian is
`0.0628180`, and the CAD-volume difference is `5.83e-11%`.
Partitioning it into sweepable HEX blocks remains an improvement task; it is
not a reason to claim a HEX computation that was not performed.

These are readiness records, not numerical certificates.  Do not cite
three-engine agreement for Examples 6 or 7 until the runner finishes from a
released native wheel on mdx or hibino and writes its checked result JSON.

### 5.4 Periodic Kelvin BDDC gauge diagnosis (2026-09-04)

The reduced-A failure on the staged Example #6 FEM mesh was isolated before
changing the physical formulation.  The `kelvin_int` and `kelvin_ext` surfaces
have 1,393 point pairs with a constant `(0.48, 0, 0)` translation and
`2.99e-16 m` maximum pairing error.  Periodic order-2 HCurl marks 9,737 minion
DOFs unused, exactly `4,173 + 2*2,782`, accounting for every interface edge DOF
and both tangential face modes.  Removing `Periodic` did not remove the NaNs.
The identification is therefore complete and is not the failing component.

The external source pullback is also independent of this failure.  The
reduced-A right-hand side is supported in iron, while the solved reaction field
extends through physical air and the Kelvin domain.  The 1-form A pullback,
2-form B pullback, and Kelvin reluctivity factor already have separate focused
tests.  Applying the BDDC preconditioner to the assembled right-hand side at
the former uniform `1e-10` gauge produced non-finite values before CG began;
the same happened on the non-Periodic HCurl space.  PARDISO SPD reported error
`-4`, and sparse Cholesky returned a non-finite solution after 150.44 s.

An exterior-only split does not fix the matrix: physical `1e-10` with Kelvin
`1e-6`, both as a constant and with a Kelvin metric factor, still gives a
non-finite preconditioner application.  BDDC inverts element-local blocks, so
curl-free local modes in the physical region also need a usable mass gauge.
With `1e-6` in every region, the actual 705,838-DoF Example #6 matrix assembled
in 15.72 s and CG completed in 28 iterations and 24.02 s, with `5.3061e-8` true
relative residual.  A separate curved order-2 35,174-DoF Kelvin regression
completed in 46 iterations with `1.6020e-10` true relative residual.
On that regression, retaining physical `1e-6` while raising only the Kelvin
coefficient to `1e-5` changed the three-probe B field by `2.53e-7` relative to
uniform `1e-6`; raising both regions to `1e-5` changed it by `1.10e-4`.

The production API consequently exposes `eps` for physical regions and
`kelvin_eps` for the compactified exterior, with `kelvin_eps=None` inheriting
the physical value.  Both default to the measured stable `1e-6` scale.  The
split is a convergence-study control, not permission to leave the physical
blocks singular.  Every linear solve now records both resolved values and its
true residual and fails loudly on a non-finite solution or a residual above
`1e-6`.  AMS remains unavailable for Periodic Kelvin HCurl because its current
auxiliary H1 construction is order-1 and NGSolve does not expose the required
Periodic low-order coupling; high-order Kelvin `auto` selects BDDC above
200,000 DOFs.

## 6. Remaining work and acceptance criteria

| Priority | Item | Acceptance criterion |
|---|---|---|
| Resolved (2026-09-08, Sugahara) | Full-versus-IMA `rad.Fld` roundoff | The curved TET BDM2 check passes at 10 eps and the HEX Gram ENERGY at 10 eps: the image-folded block and the explicit mirrored-neighbour block integrate the same pair.  The two HEX FIELD checks sit at `2.2e-14` and `5.7e-14`.  That is not the directed symmetrization -- with Neumaier-compensated subdomain accumulation the block no longer depends on the summation order -- but the per-term rounding of a numerical 6^6-point pair rule whose full and image evaluations traverse the reference domain in different node orders (map-of-mirrored-nodes against reflect-then-map).  Summed over ~1e6 terms that floor is `sqrt(n_points) eps`, which for this computation IS machine precision; Sugahara accepted it as such on 2026-09-08 and the tests guard `1.2e-13` with the derivation written into them.  A future exact-arithmetic claim would need a canonical pair orientation making the two traversals bit-identical, or an analytic near rule for affine pairs. |
| Resolved (2026-09-07) | Remove production RT0 | The public `vim.Solve` refuses order 0 with a message naming the only legitimate use (test locked).  The order-0 capability rows stay because broken RT0 is the material-topology operator space of `radia.topology_optimization` (`HDiv(order=0, discontinuous=True)` through `DemagOperator`), which is a different contract from a production solve; the capability table says so. |
| P1 | Nonlinear C-yoke memory evidence | Four-level accuracy and repeated timing are closed on mdx and hibino for `v4.95.71`. Add measured process peak memory to a future scaling campaign before making a memory-efficiency claim. reduced-A remains an independent third-formulation audit rather than the primary production route. |
| P1 | ESRF #3 H-matrix three-engine evidence | Run the repaired `leaf=64` operator on mdx or hibino through the tracked nonlinear three-engine runner. Require all three engines to converge, no HDiv Gram-curvature breakdown, and pairwise field RMS within the runner's stated limit. |
| #6 closed (2026-09-08), #7 open | ESRF #6 and #7 three-engine evidence | Run the new coil-yoke runner from the released native wheel on mdx or hibino. Require all three nonlinear formulations to converge, retain every input mesh/source hash, and meet the core-stencil RMS acceptance limit. |
| Closed for #6 (2026-09-08), open for #7 | Released reduced-A Kelvin BDDC replay | Install the merged wheel on hibino or mdx and rerun Examples #6 and #7. Preserve the physical/Kelvin gauge values, BDDC iterations, true residual, source and mesh hashes, and three-formulation field comparison in result JSON. |
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
- reduced-A Kelvin BDDC candidate: 17 dispatch/contract tests PASS; curved
  order-2 35,174-DoF Periodic HCurl solve PASS in 46 iterations with
  `1.6020e-10` true relative residual. The real Example #6 diagnostic gives
  28 iterations and `5.3061e-8` at 705,838 DoFs with a uniform `1e-6` gauge;
  the released-wheel three-engine rerun remains open.
- exact-dense ChargeGram entry, public `DemagOperator`/`vim.Solve`, and
  configured-principal-submatrix semantics: 9 tests PASS in 5.87 s on LAB;
- ESRF #3 exact-dense material diagnostic: PASS under its explicit 1 GiB cap.
  The historical parallel H-matrix failure at PCG iteration zero is explained
  by section 5.2. The repaired H-matrix must still complete the tracked remote
  nonlinear three-engine runner before it becomes release evidence.
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

## 8. Curved-TET mirror images, Gram definiteness, and the FEM Picard loops (2026-09-05)

Branch `claude/hdiv-ima-curved-tet` (on the PR #93 head `35ebfd0cb`).  Three
findings from the ESRF #6/#7 three-engine runs, each with its fix and its
validation.  Timing numbers are LAB same-host relative smokes; decision-grade
timing still belongs on idle mdx/hibino.

### 8.1 The IMA build of ESRF #7 was a per-entry scalar curved Duffy

The #7 one-pole model is not a HEX mesh: `model_one_pole_20mm.vol` is 15,210
curved P2 TET elements (`check.json`: `tetrahedron_count 15210`, `curve_order
2`).  Its reduced BDM1 build with `image="-x-y"` ran for more than 61 minutes
on hibino while the FULL 30 mm model (31,988 curved TET, no image) solved BDM1
in 403 s nonlinear and 256 s linear.

Cause (`rad_hacapk_hdiv_entry.cpp`, `HighOrderTetEntryStrategy::Evaluate`):
every MIRROR image term took the scalar fold
`0.5 (QuadDotRefl(a,b) + QuadDotRefl(b,a))`, and on a curved host `PhiInner`
is `CurvedTetPotential` -- 4 faces x 3 leads x 8^3 = 6144 curved-map
evaluations per outer point per source charge, both directions, every image,
far pairs included.  The host-block path existed only for ROTATION images on
FLAT hosts (`QuadBlockHOTetImage` threw for curved), and no far rule existed
for images at all, while the DIRECT terms used the product rule, the far rule,
and the vectorized Duffy per host.

Fix (`rad_hacapk_hdiv.cpp`, `rad_hacapk_hdiv_entry.cpp`, header):

| Piece | Rule |
|---|---|
| `ImageFarPair(a, b, img)` | the direct far criterion on the IMAGE geometry, distance from `T^-1 c_a` to `c_b` above `f (s_a + s_b)` |
| `QuadDotFarImage` | `QuadDotFar` with the target's low outer points mapped by `T^-1`; a mirror gives the same sum in a different order, a rotation gives `G_T` and `G_{T^-1}` which are averaged as designed |
| `ImageHostsTouch(T, S, img)` | S's corners mapped forward (`ImageApplyVector` on positions) and matched to T's corners by coordinates: plane-fixed vertices and rotation-identified sector vertices are both found; the vertex-id test `CurvedHostsTouch` sees neither |
| `QuadBlockHOTetImage` (curved) | touching image pair: vectorized curved Duffy at the mapped points; otherwise `QuadBlockHOCurvedDirect(img)` product rule; mirror + product rule is an exact transpose, so it is one-sided |
| entry dispatch | far -> host block (curved, or flat with the analytic host block) -> scalar fold only for flat BDM1 / polynomial-combination charges |

An on-plane cut face maps onto itself point by point, so its image self block
is evaluated at the same outer points by the same rule as its direct self
block and the antisymmetric fold cancels to roundoff -- the curved analogue of
the hex `self_pair` fix of 2026-07-05.  The legacy fold stays reachable for
A/B through `RADIA_HDIV_DISABLE_HO_IMAGE_BLOCK=1` plus the new
`RADIA_HDIV_DISABLE_HO_IMAGE_FAR=1` (classified as a numerical-path override
above); `hmat_stats` gains `ho_image_far_entries`, `ho_image_block_entries`,
`ho_image_scalar_entries` under `RADIA_HDIV_BLOCK_CACHE_STATS=1`.

### 8.2 Validation: entries, symmetry, definiteness, physics, speed

`validation_test/feec/validate_hdiv_vim_tet_image_dispatch.py` (JSON next to
it) and the fast `tests/feec/test_hdiv_vim_tet_image_dispatch.py`.  A quarter
model (x > 0, y > 0) with `image="-x-y"` under the quadrupole field
`H_ext = g (y, x, 0)` -- the parity of the ESRF quadrupoles -- against its
full model.  For the flat box the full mesh is the exact mirrored union of the
quarter mesh (`mirrored_union`), so the agreement is roundoff + ACA level; for
the curved sphere the two meshes are independent.

| Case (maxh 0.5, `gram_eps` 1e-12) | n_charge | entry symmetry | on-plane residue | lambda_min new | lambda_min legacy | reduced+image vs full | build new / legacy |
|---|---|---|---|---|---|---|---|
| flat box BDM1 | 280 | 0 | 3.4e-15 | -2.9e-16 | -1.0e-9 | 6.4e-15 | 0.38 s / 0.50 s |
| flat box BDM2 | 712 | 0 | 5.5e-15 | -4.6e-17 | -3.8e-9 | 1.1e-12 | 2.7 s / 29.7 s |
| curved sphere BDM1 | 195 | 0 | 8.5e-16 | -5.9e-17 | -1.3e-7 | 5.3e-4 | 1.7 s / 26.3 s |
| curved sphere BDM2 | 480 | 4e-18 | 8.1e-15 | -9.5e-18 | -4.9e-9 | 2.4e-4 | 7.6 s / 820 s |
| curved sphere BDM1, maxh 0.3 (build only) | 564 | -- | -- | -- | -- | -- | 5.3 s / 164 s |

`lambda_min` is the smallest eigenvalue of the sigma-normalized dense Gram
assembled from `matvec_sym`; the raw O(n^2) quadratic form on the minimizing
vector agrees in sign.  The build ratios grow with the mesh (1.4x for flat
BDM1, where the analytic host block is off and only the far rule changes;
11x flat BDM2; 15x and 31x curved BDM1 at 45 and 150 elements; 108x curved
BDM2), because the legacy fold paid the curved Duffy for every filled entry
while the dispatch pays it once per touching host pair.  Two points matter
beyond the speed:

* the LEGACY fold was itself slightly indefinite -- exact/Duffy image terms
  combined with product-rule direct terms are not one quadrature family, and
  the mixture leaks into the smallest eigenvalues (-1e-9 flat, -1.3e-7
  curved, at these tiny sizes).  The consistent dispatch is PSD to roundoff.
  This is one mechanism for "IMA + TET -> CG breakdown"; it does not explain
  codex's raw O(n^2) HEX Gram of #6, which is a separate defect and stays on
  the MINRES-as-diagnostic-only rule;
* the on-plane cut-face charges of the antisymmetric planes annihilate to
  about 1e-15 of the median self entry with sigma left at one, as the sigma
  pre-pass contract requires.

Every existing IMA, cyclic, hex-image, sigma, and roundoff test passes on the
patched build.  Three failures on the PR #93 head reproduce on an UNPATCHED
baseline build and are therefore pre-existing: the two hex-image roundoff
tests of `validation_test/feec/test_hdiv_radfld_contract.py` (4.9e-14 against
a 10 eps gate) and
`tests/test_hdiv_vim_chargegram_dispatch.py::test_chargegram_curved_tet_matches_mesh_geometry_by_default`
(its test double returns four values where `_finish_charge_gram_backend`
unpacks three).

### 8.3 The FEM Picard loops: history, warm start, constrained Anderson

ESRF #6 mixed Omega stopped at 80 damped-Picard iterations (relaxation 0.3)
with `relative_B_change` 9.03e-5 against 2e-5, and the rerun with a 160 cap
started from zero because the loop raised and discarded its state.  Both FEM
engines were fixed-relaxation Picard on the per-element material coefficient
with a step-size criterion, no history, no warm start.

Shipped: `radia.picard_acceleration.ConstrainedAndersonAccelerator` (real
arithmetic; projection onto the secant range of the B(H) law; extrapolation
in log space by default; restart on residual growth; an a-posteriori
acceptance that drops an accelerated iterate whose next residual is worse and
takes the damped Picard step from the accepted one instead; depth 0 is the
legacy convex combination bit for bit) and `estimate_contraction_rate`.
`solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin` and
`VectorPotentialSolver.solve_nonlinear` accept a per-element warm start
(`mu_r_initial` array / `nu_initial`), `anderson_depth`, `anderson_transform`,
`observation_points`; their `nonlinear_stats` carry the per-iteration
`history`, `contraction_rate_estimate`, the per-element state, the Anderson
counters, and the observed field.  The mixed loop raises
`MixedOmegaPicardNotConverged` with that state; the reduced-A loop keeps
returning, and its silent `except Exception: B_mag = 0.0` centroid fallback is
now a raise.

Measured on the fast tests: reduced-A saturating cube (tol 1e-6) cold 24
iterations, warm start from the converged state 2, Anderson(2) 14, fields
agreeing to 3e-7.  Mixed Omega on a shielded-knee two-region case (H_knee
1e-4 A/m): plain rate 0.958 (1.3e-5 after 150), Anderson(2, log) 48
rejections in 150 and no gain, linear transform worse.  Anderson is therefore
opt-in in the runner; the warm start, the history, and the rate estimate are
the guaranteed wins.

### 8.4 Runner contract: converged-only checkpoints, caps as provenance

`run_coil_yoke_three_engine.py` (checkpoint schema v3): a result checkpoint
is written and read only for a converged solve (the v2 `reduced_a` checkpoint
of #7 with `converged: false` would have been reused silently and failed the
gate hours later); the iteration caps left the checkpoint identity and live in
`provenance`, so a converged solution is reused whatever cap produced it, and
converged v2 checkpoints still resume after their cap key is stripped; a
non-converged engine writes `<output>.<engine>.state.json` (explicitly
`converged: false`) which `--resume` uses as a warm start; new options
`--mixed-relaxation`, `--mixed-anderson-depth`, `--reduced-a-anderson-depth`
(default 0, part of the identity only when set).

### 8.5 Codex handover verification (2026-09-05)

The curved-TET dispatch test double now preserves the native configuration
contract: it receives the four build inputs, including the NGSolve mass,
and returns the public three-tuple. No production return contract or numerical
tolerance was changed. All six dispatch tests pass after this correction.
The other 40 focused tests covering Picard acceleration/history, mixed Omega,
ESRF checkpoint contracts, and native TET image dispatch passed.

Tests used the built source/native pair at `C:/temp/radia-hdiv-ima/src`
(commit `8c5b071693961564519effa6463feb531e0f2b1a`), explicitly imported before
pytest, with the updated tests from the isolated handover worktree. No PYD was
copied. This is not evidence of a new wheel build or deployment.

Release remains blocked: the HEX image field roundoff tests reproduce relative
differences of 1.9927109321574784e-14 and 4.944633724445392e-14 against the
unchanged 10-epsilon gate (2.220446049250313e-15). The full field-contract file
has three passes and two failures. Their cause is not established by this run.
The #6/#7 production three-engine acceptance and timings require the new native
release on mdx/hibino; diagnostic MINRES results do not certify production CG.

The LAB release-worktree editable path is not repaired blindly: release-quad
intentionally retains the verified release source until the canonical development
checkout catches up. Use its explicit restore-editable command only after that
prerequisite is satisfied; do not redirect imports to the stale dirty shared tree.

### 8.6 ESRF #6/#7 restart gate (2026-09-05)

Both hibino and mdx were queried through SSH: no Python compute processes were
listed, and both imported installed PyPI Radia 4.95.81 from site-packages.
Neither host yet contains the new native TET image implementation. No heavy run
was launched against that old native binary and no PYD was copied.

The saved #6 MINRES diagnostic checkpoint declares nonlinear convergence but
contains no native implementation identity. Its hmat release-eligibility flag
does not establish which material linear solver actually ran. Consequently,
convergence alone is insufficient to reuse it as production evidence.

The #6/#7 runner now includes SHA-256 identities of the loaded native library,
NGSolve native library, formulation modules, acceleration module, case source,
and shared/current runners in every result/state contract. Old or differently
built checkpoints fail the identity comparison instead of silently certifying
a new implementation. Original checkpoints remain untouched. The focused mesh
and checkpoint contract suite passes all 12 tests, including rejection of a
converged diagnostic-build checkpoint under a production-build identity.

Three-engine acceptance is still incomplete: #6 requires correction of the
raw HEX Gram indefiniteness and a production solve; #7 requires the new native
image build and converged FEM comparisons. The release-owning task has been
notified of the deployment dependency. This section records restart safety,
not successful field agreement or completion of either case.

### 8.7 Installed-wheel trial and material-envelope defects (2026-09-05)

With user approval, a candidate wheel was built from `6d66d38f9` using
`Build.ps1` and `Build_Wheel.ps1 -DryRun`, without publishing it. Its SHA-256 is
`8c34e57bd9ef9dc98ea788941266133f1cea06e5625be98897428480a9b58a7d`.
Hibino installed it into `C:/temp/radia-candidate-6d66d38f9/venv`, with
read-only access to system dependencies and candidate-local Radia, threadpoolctl,
and test tools. All 336 package files matched the wheel; pip check passed;
the installed native TET image tests passed 8/8. No standalone PYD was copied.
The version string remains 4.95.81: this is a hash-identified private candidate,
not the published 4.95.81 artifact.

Example #7's curved P2 TET one-pole mesh (15,210 elements), BDM1, image `-x-y`,
32 threads, gram epsilon 1e-10 converged through the production Newton/CG path:
96,987 field DoF, 50 Newton iterations, 9,357 inner iterations, Gram 153.164 s,
solve 184.652 s, internal total 339.812 s, and 531.482 s including direct field
evaluation and adapter work. The result is preserved in the candidate results
directory. This establishes a working native image path, not three-engine
agreement. The FEM trial exposed two material-law defects:

1. Reduced-A initialized its lower reluctivity bound from the first positive
   H/B sample and never lowered it. Real initial magnetization curves may have
   rising secant permeability: Example #7 has initial H/B 354.776115 but a
   minimum tabulated H/B of 99.032672. Even depth-zero Picard updates were
   clipped, changing the material law. A small FE regression using that real
   table fails in all 23 iron entries on the previous implementation.
2. Mixed Omega bounded permeability by the maximum tabulated secant, but a
   monotone PCHIP B(H) can have a larger secant between nodes. A regression
   reaches 2176.768 against a tabulated cap of 2000 and fails on the previous
   implementation.

Both loops now expand the acceleration envelope to include evaluated material
targets (and Reduced-A's supplied warm state). Thus a constitutive fixed point
is not excluded by an assumed monotonic secant range, and depth-zero remains
the actual convex Picard update. The generic Anderson bounds API is unchanged.
The invalid FEM trials were stopped; their intermediate fields are not accepted
as comparison results. New-wheel three-engine validation is still required.

### 8.8 Corrected candidate installed and HEX gate reproduced

The material-envelope fix is commit `c11ed1b2b`. Its separate candidate wheel
(`a9c3ee6cdf2263eb020a21f07514ddc9c254cf129a2fd9f2d250a1b5d1ecac90`)
was installed in `C:/temp/radia-candidate-c11ed1b2b/venv` on hibino. All 336
Radia files match the wheel, pip check passes, and 31 installed-wheel native/FEM
tests pass in 56.26 s. The local focused suite passes 41 tests. Both new
material-envelope tests were independently run against the previous Fable
implementation and failed there.

Example #6 was then run with the corrected candidate, BDM1, FEM order 1,
gram epsilon 1e-10, 32 threads, and no numerical overrides. Production CG
failed at iteration 86 with p^T A p = -1893733131.842714. No MINRES fallback
was used. This is a release-blocking HEX operator defect, not evidence that a
larger iteration cap is needed. Neither FEM engine ran after that failure.

Machine-readable assessment and the successful earlier #7 HDiv field vector
are retained under `validation_test/esrf_three_engine/results/`. The #7
checkpoint explicitly belongs to the first candidate; its native hash is
unchanged by the FEM-only fix, but corrected full-model FEM agreement has not
yet been demonstrated. No three-engine acceptance, main merge, tag, or PyPI
publication is claimed. Both invalid FEM pilots and the failed #6 process
have stopped; the candidate environments remain available for further work.

### 8.9 ESRF #6 root cause: HEX charge-Gram definiteness (2026-09-05)

The iteration-86 CG breakdown of section 8.8 reproduces on LAB from the
production entry point (`solve_configured_linear_material_auto_prec`, chi0
warmstart of the energy-Newton path, `p^T A p = -1.3e9` at iteration 86 versus
`-1.9e9` on hibino with 32 threads): deterministic, not a threading race.  Along
the breakdown direction `p^T W p = +4.95e12` and `p^T N p = -4.96e12`, and the
raw O(n^2) quadratic form is `-4.52e12`, so the charge Gram `N = B^T G B` itself
is indefinite; codex's raw-Gram observation is confirmed.  The Newton tangent is
SPD (the case-6 table has `dM/dH >= 0.16` everywhere), every exact 13-host
element cluster is PSD to 1e-16 in the M-metric, and a floor scan with the
production CG brackets `lambda_min(M^-1 N)` between `-5e-3` and `-2e-3` (mu_r
2001, 1001 and 501 break at iterations 86, 105 and 156; mu_r 201 converges in
253 iterations).  The physical band is [0, 1].  Two independent defects produce
this, and both are now fixed on `claude/hdiv-hex-gram-psd`.

**H-matrix admissibility (the larger defect).**  A preconditioned generalized
LOBPCG finds `lambda_max` Ritz 2.12 whose raw quotient is 0.83; the symmetric
leaf diagnostic pins the excess to one low-rank leaf (rank 5, 80 x 80, H-matrix
quadratic `-0.67` versus raw `-1.96`) that couples two mirror cells touching
across x = 0, while every other dominant leaf agrees to 1e-4.  The cluster-tree
points were the co-located charge centroids, so HACApK's box-gap admissibility
(`width <= eta * gap` in `cHACApK_bndbox` / the leaf generators) saw a gap of
one cell between clusters whose hosts touch, declared the block admissible, and
ACA+ stopped at rank 5.  Fix: `cHACApK_set_point_radius` inflates the leaf
bounding boxes by a per-point support radius (the host bounding radius
published by the hex Gram's `ExtractCoordinates`), so touching hosts have gap 0
and always land in dense leaves while a host's modes stay co-located.
(Spreading the modes over the lattice nodes instead split hosts across clusters
and put self entries into low-rank leaves: `lambda` in [-790, 929].)  The
diagnostic latch `RADIA_HDIV_HEX_CLUSTER_RADIUS=0` restores point boxes and is
reported as a numerical override.

**Near family of distorted cells (partially resolved).**  The #6 non-affine
cells are trilinear distortions (Q2 mid-node deviation 1e-15), not curved.
Pair-level comparison against fine rules showed the distorted SELF blocks 8e-4
low (up to 5.6e-3), the BDM face-dof self-energy 7e-3 low, and touching pairs
classified "not near" at `near_grade` 0.5 (their centroid ratio is 0.56..1.0),
so the sub-tet outer was the regular rule against a boundary-singular
potential; the static-site radial inner added ~1e-3.  Against a fine reference
on a 72-cell elongated sector lattice the legacy family is -11 %..+18 % off in
the M-metric (`lambda_max` 1.084 against the physical bound 1).  Changes:
touching hosts (a shared Q2 lattice node, `HexHostsTouch`) are always near and
never take the far tensor product (also affine pairs, whose vertex-touching
ratio reaches 1.0); the near outer rule `glnear_n` (default 8) is decoupled
from `glout_n`, which the far tensor product shares; every near pair (self,
touching, near band) is integrated on ONE endpoint-graded tensor rule over the
whole target host (`QuadBlockHexNearTensor`, smootherstep grading on every
axis) so all sources of a target share the same outer point set (mixing the
corner-graded sub-tet self outer with tensor touching pairs, although each block
was more accurate, made the assembled Gram worse: `lambda_min` -4.3e-3); the
inner is the exact-anchor radial (`HexQ2Inverse` / `QuadQ2ClosestReference`
anchors, `PhiInnerHexRadialVec`) with a finer rule `glin_self_n` (default 12)
for self pairs and face sources (the endpoint grading puts outer points near
the sub-tet faces, where the 5-point cones are coarse; the remaining #6
negative mode, -5.0e-3, was pure face-face energy on distorted cells and
disappears with the finer face rule); `near_inner="site"` keeps the legacy
inner as a flagged A/B path.  On the sector lattice the error against the fine
reference drops to -1 %..+9 % (glnear 8 / glin 5 / self 12; -2 %..+2 % at
glnear 10 / glin 8 / self 12), and touching-pair blocks on #6 move within 7e-5
of the fine reference (legacy 1.4e-4..2.4e-4).  It is not converged: the sector
`lambda_max` is still 1.03 (strict xfail in
`tests/feec/test_hdiv_vim_hex_near_family.py`), and on #6 the production CG
still breaks at iteration 100 at the chi0 floor: along the breakdown direction
the Rayleigh quotient of N is -5.08e-4 in the M-metric (raw and H-matrix
agree), 1.6 % below the floor 5.0e-4, spread over distorted (60 %) and affine
(40 %) cells; the preconditioned LOBPCG cannot separate such a direction from
the exact null space of charge-free fields.  A floor scan with the production
CG still breaks at mu_r 1334, 1001 and 501 (iterations 148, 160, 292), as the
legacy family did, so the lower band of the raw Gram is not yet improved on
#6 -- the sector error and the H-matrix consistency are.  The converging direction is a near
rule whose error is independent of the cell aspect ratio; the tensor outer +
exact-anchor radial is the consistent frame for it.

**Fail-loud Jacobi diagonal.**  `SolveConfiguredLinearMaterialAutoPrec` and
its multi-RHS twin silently replaced a non-positive exact diagonal of
`inv_chi*M + B^T G B` by 1.0; they now raise with the DOF and the value
(`tests/feec/test_hdiv_jacobi_diagonal_fail_loud.py`).

**Gate.**  `validation_test/esrf_three_engine/validate_hex_gram_definiteness.py`
rebuilds the production Gram on the #6 asset and fails unless the exact
clusters are PSD, the production CG converges at the chi0 floor and at 2x, 4x
and 8x larger initial permeability, the preconditioned LOBPCG edges stay inside
[-1e-8, 1 + 1e-3], and the raw O(n^2) and H-matrix quadratic forms agree along
the minimizing direction (the separate compression check).  Its LAB result is
recorded in `results/hex_gram_definiteness_lab.json`: clusters PSD, raw and
H-matrix quadratic forms agreeing to 1.8e-14 along the minimizing direction,
`lambda_max` Ritz 1.051 (2.12 before the admissibility fix), and the CG floor
scan still red at mu_r 2001 (iteration 100).  Timings are relative (LAB, 421 s
build against 120 s for the legacy family); the idle mdx/hibino run is codex's.
`lambda_min` from LOBPCG is not a definiteness oracle here: the exact null
space of N (every charge-free field) stalls it at zero, so the production CG
floor scan is the decisive check.

### 8.10 Cone/fan near inner and the M-metric amplification (2026-09-06)

The near family of section 8.9 was replaced by a whole-host rule
(`PhiInnerHexConeFanVec`): six face cones from the apex -- the outer point of a
self pair, the physical closest point of the source host for a touching pair
(`HexQ2ClosestReference`) -- each face integrated by four edge fans from the
apex's physical foot, with Johnston-Elliott sinh substitutions along the ray
(the near-singular peak `1/sqrt(d^2 + r^2 |v|^2)` of a touching pair), the fan
radius (apex close to the face) and the edge parameter (foot close to an
edge).  Two facts drove it: a thin cone (apex near its face) has a `1/rho`
peak in the face integral that a tensor Gauss rule cannot integrate (the error
is first order in the apex distance and moves erratically with the point
count: sector errors -6 %..+3 % at 8 points, -0.8 %..+12 % at 12), and the
edge direction carries no vanishing Jacobian, so its peak mass is `b ln(1/b)`
and the substitution must be applied for every positive width (a `1e-3`
cutoff left 12 % errors on outer points `1e-4` from an edge).  Affine sources
keep the exact analytic inner.  Result: on a lattice warped by `1e-7` the rule
matches the exact inner to `3e-5` (cell self `2.5e-5`, face self `5e-6`,
touching `3.4e-5`); the elongated sector lattice sits inside the physical band
(`lambda_max` 0.998 against 1.084 legacy and 1.03 for the sub-tet radial), so
`test_sector_spectrum_stays_in_physical_band` is a plain pass; the per-class
residual against `glnear` 12 is `1.6e-4` and comes from the outer rule
(`glin_self` 8 vs 12 changes `3e-6`).  The default `glin_self_n` is 8.

The affine product family (exact inner, plain `glout` 4 outer) turned out to
be the coarser one: face self blocks `+4.1e-3` (converging only
algebraically, `+1.1e-4` at `glout` 10, `+1.5e-5` at 16), touching couplings
`1.7e-4`; the graded 8-point outer is `2.4e-6` on the flat face self integral.
BDM1 therefore routes the whole affine near band (self, touching, and the
non-touching pairs inside the `HEX_AFFINE_EXACT_NEAR_FACTOR` band) through
the graded near tensor outer as well, so every near source of a target host is
integrated on one cloud; BDM2 keeps the affine product.

None of this passes #6.  The CG breakdown direction of the cone/fan family is
identical to the old family's (quotient `-5.086e-4`, same cells, 60 %
distorted); a mesh perturbed by `1e-10` m so that every cell is non-affine
(uniform fan family) still breaks at iteration 98 (`p^T A p` `-1.2e9` against
`-1.8e11`); routing only the touching affine pairs onto the graded cloud
exposed a 4-fold degenerate face mode at `lambda` `-5.009e-3` (self face
`+0.65`, touching face-face `-0.46`, near band `-0.20` in units of the mode's
M-norm, almost entirely on affine faces); the consistent band routing leaves a
mode at `-2.27e-3` and the CG at iteration 98 (hibino).  The reason is
structural: on the flat lattice the same entry errors that are `4e-6` relative
give an M-metric error spectrum of `[-1.0e-3, +4.0e-4]`, and `2e-7` entries
give `[-4.9e-5, +3.2e-5]` -- an amplification of about `1e3`.  The modes
nearest zero are nearly charge-free fields whose cell-divergence and face
charges cancel, so their energy is a small remainder of large self and cross
block energies, each integrated by its own rule; the block errors do not cancel
the way the charges do.  Block-wise quadrature would need entries accurate to
about `1e-8` on rules with `x ln x` edge behaviour, which is not a tuning
target.  The way out is a common energy form: either the Ewald split
`1/r = erf(alpha r)/r + erfc(alpha r)/r`, with the smooth positive-definite
part integrated on one global point set (positive semidefinite by construction
for any charge samples) and only the short-range part, which has no
cancellation structure, integrated per self/touching block; or the assembly of
each BDM DOF's composite charge (cell divergence plus boundary face charge,
zero net) with one rule per DOF pair, i.e. the dipole-kernel formulation with
analytic element integrals.  Heavy runs are hibino's (Gram build 350 s at
45,792 DoF); LAB numbers above are relative.

### 8.11 Pair-domain Duffy quadrature for touching HEX pairs (2026-09-06)

Design (c) of section 8.10, the standard answer of the Galerkin BEM/VIE
literature (Sauter-Schwab regularizing transformations; Reid's Taylor-Duffy
for tetrahedron products), is implemented for BDM1: every TOUCHING host pair
(self, shared face, shared edge, shared vertex; cell-cell, cell-face and
face-face; affine or not) is integrated on its `(d_T + d_S)`-dimensional
product domain by `QuadBlockHexPairDuffy`.  `HexPairAdjacencyOf` reads the
shared entity and the canonical frames (axis permutation plus flips per host)
off the coincident Q2 lattice nodes after the image transform.  In canonical
coordinates the relative in-entity coordinates `u = zeta_S - zeta_T` and the
transverse coordinates form a cone vector whose max-norm `w` is the Duffy
variable: one subdomain per dominant coordinate (two signs for a relative
one) times the signs of the other relative coordinates -- the intersection-box
lengths `1 - |u_i|` are smooth only on a fixed sign, and without that split the
rule converged algebraically (unit-cube self-energy `-3.1 %` at 4 points,
`-1.4 %` at 6).  The Jacobian `w^(k-1)` cancels `1/r = 1/(w X)` with `X`
smooth and nonvanishing for any Q2 map, so the analytic radial reduction of
Taylor-Duffy (which needs affine elements) is not required: tensor Gauss on the
unit hypercube converges exponentially.  Measured with `glpair_n` points per
dimension (`pair_duffy_check.py`, LAB): the unit-cube Coulomb self-energy
`1.88231264438961` is reproduced to `2.5e-6` (4), `4.9e-9` (6), `1.8e-12` (8);
the unit-square self-energy `2.9732095982` to `1.8e-6`, `4.1e-9`, `1.2e-11`;
on a `2x2x2` lattice warped by `1e-7` every touching block changes by `1.4e-4`
(4 to 6), `6.2e-9` (6 to 8) and `4.6e-12` (8 to 10).  The default is 8.  The
non-touching pairs inside the near band take the plain product rule with the
same point count (`QuadBlockHexProductN`; the hosts are separated, so the
integrand is smooth), after the class-wise comparison of the `-5.0e-3` mode
showed that band as the largest remaining error class (`+2.9e-3` of the mode's
energy moved when its rules were refined, against `3e-5` for the touching
face-face blocks).  `RADIA_HDIV_HEX_PAIR_DUFFY=0` restores the block-wise near
family for A/B and is reported as a numerical override.

With the pair rule, the near-band product rule and the non-conforming
fallback in place the #6 gate on hibino (Gram 1220 s, 45792 faces) still broke
the production CG at iteration 93 (`p^T A p = -2.2e11`), with the element
clusters PSD (`lambda_max 0.9004`), the LOBPCG `lambda_min` Ritz stalled at
`-4.9e-7` and `lambda_max` Ritz `1.027`.  That pointed away from quadrature and
at the mesh itself (section 8.12).

### 8.12 The #6 mesh was non-conforming: unmerged Cubit partitions (2026-09-06)

`HexPairAdjacencyOf` refused a face pair that shares exactly two lattice nodes,
and a node-coincidence survey of the 4768 hosts of the #6 mesh explained the
remaining negative energy: 512 pairs of geometrically identical boundary faces
and 832 hanging-node contacts (2:1 transitions).  The reason is the journal
generator, not Cubit: `export_esrf_cubit_assets` imported the forty
partitioned iron solids and meshed them with `scheme auto` plus a `submap`
retry, but never issued `imprint volume all` / `merge volume all`, although
the partition was designed (docstring of `build_esrf_cubit_hdiv_iron`) so that
the exporter drops the shared same-material surfaces -- which it does only for
merged volumes.  Every constructive interface therefore left the mesh as two
coincident boundary faces carrying opposite surface charges, and unequal
neighbour intervals left hanging nodes.  The `+-sigma` twins are exactly the
`+0.0543 / -0.0543` "self face / touching face-face" cancellation of the
`-5.0e-3` mode in section 8.9, and the hanging-node contacts (no canonical
frames, graded near family) carry the residual seen after section 8.11.
Examples 3, 5 and 7 were meshed by the same journal and had the same defect.

The fix is in the journal generator: every solver journal now imprints and
merges before the sideset, and #6 -- whose forty solids are all 60 mm
extrusions along the beam axis -- sweeps every volume explicitly from its
lower to its upper end faces inside an APREPRO loop (`volume {_v} scheme sweep
source surface in volume {_v} with x_coord < lo+tol target surface ... >
hi-tol`), because `scheme auto` and `submap` cannot interval-match the eight
merged hyperbolic pole tips, one sweep command may not name several volumes,
and the tolerance must be `1e-3` of the extent (the imprinted lateral faces of
the tips do not all span the full length; a quarter-extent tolerance pulled
them into the source set and Cubit demanded multisweep).  Headless Cubit
2025.12 on LAB: 2408 HEX, 2200 boundary faces, 3616 nodes, order-2 curving,
`check-vol` PASSED, conforming (0 hanging facets, 0 duplicated faces), and
byte-identical from the python probe and the generated journal.  Examples 3, 5
and 7 still mesh and export conforming meshes (144 HEX, 1112 HEX, 367845 TET).
`validate_hex_gram_definiteness.py` now surveys the conformity and refuses a
non-conforming mesh unless `--allow-nonconforming`; the mesh policy records
`iron_sweep_axis` and `conforming_partition`.

**Gate result (hibino, conforming mesh, pair-domain Duffy family,
`results/hex_gram_definiteness_hibino.json`): PASSED.**  Gram 1567 s for
62192 face unknowns; element clusters PSD (`lambda_min -3.8e-16`,
`lambda_max 0.8891`); the production chi0-warmstart CG converges at every
floor of the scan -- 431 iterations at `mu_r ~ 2001`, 486 at 4001, 542 at
8002, 598 at 16003 -- where the non-conforming mesh broke down at iteration
93; LOBPCG `lambda_min` Ritz `+5.0e-21`, `lambda_max` Ritz `0.99990` (inside
the physical band `[0, 1]`, against `1.027` before); raw versus H-matrix
quadratic form along the minimizing direction agree to `1.3e-33` (M-metric).
The quadrature work of sections 8.9-8.11 stays (it is what makes the merged
mesh's touching pairs consistent to `1e-12`), but the defect that made
example 6 alone indefinite was the mesh.  Next validation target: the CEFC
2020 Q-mag quadrupole (`validation_test/quadrupole_cefc2020/`).

### 8.13 HEX Gram build cost: the near blocks, and the translation-congruent cache (2026-09-06)

Where HEX stands against TET (all timings mdx/hibino, committed JSON): the
C-type three-engine nonlinear BDM2 run solves 32580 TET face unknowns in 37 s
(1.1 ms per unknown, 10-23x faster than the two FEM formulations at 0.3 %
agreement), while the HEX BDM1 quadrupoles need ~1500 s for ~61000 unknowns
(25 ms per unknown).  The `gram_stats` profile of the Q-mag `h = 10 mm`
linear run (hibino, 38 threads, thread-summed seconds) locates the whole gap:

| dispatch class | blocks | thread-seconds | per block |
|---|---|---|---|
| `hex_blk_general_near` (pair-domain Duffy + near-band product) | 137,341 | 36,458 | 265 ms |
| `hex_blk_affine_far` | 10,253,212 | 146 | 14 us |
| `hex_blk_distorted_far` | 3,519,648 | 43 | 12 us |

The near family is 99 % of the build: every touching pair (the pair-domain
Duffy rule, ~72 % of the near blocks) and every non-touching pair inside the
near band (`QuadBlockHexProductN`, ~28 %) integrates `8^6 = 262,144` point
pairs.  Two cache defects multiplied that cost.  `HexPairTakesGeneralPath`
returned false for affine-affine pairs, so on a mesh whose cells are 74 %
affine most near blocks bypassed the instance-shared cache and were
recomputed by every fill worker that touched them (27,118 shared lookups
against 137,341 evaluations).  And the translation cache required every host
to sit on the half-cell lattice of one affine cell (`hex_uniform_trans_hosts`
is false on any real magnet), although a swept mesh -- every 2.5-D magnet --
repeats each host once per layer.

Both are fixed in the kernel: every BDM1 touching or near-band pair now takes
the shared cache, and the cache key is the TRANSLATION-CONGRUENT pair
(`HexSharedBlockKey`: the two host templates plus the centre offset quantized
to `1e-10` of the largest host spread; `BuildHexCongruenceTemplates` hashes
each host's Q2 lattice nodes relative to its centre together with its charge
exponents).  A charge-Gram block depends only on the relative geometry of its
hosts, so every translated copy of a pair -- 5 of 6 near blocks on the 6-layer
10 mm meshes, 14 of 15 at 4 mm -- is served from one evaluation, exactly.
Image blocks (`img > 0`) keep the host key.  `tests/feec/test_hdiv_vim_hex_congruent_cache.py`
locks the mechanism on a graded (non-lattice) swept mesh: the cache engages
(templates far fewer than hosts, shared hits above 30 % of lookups) and the
Gram equals the uncached one (`RADIA_HDIV_DISABLE_CONGRUENT_CACHE`, a
performance latch reported in `hmat_stats`) to `1e-12`.  The hibino timing of
the Q-mag and example-6 builds with and without the cache is the next entry
of this section; the further levers, in order, are the pair point count
(`glpair_n` 8 -> 6 is 5.6x on every near block at `6e-9` self-energy accuracy,
to be confirmed by the definiteness gate), a distance-graded count for the
non-touching band, and the analytic inner for affine-affine touching pairs if
the conforming-mesh A/B shows the block-wise family suffices there.

The same cache defect explains the slow phases AFTER the build.  Every later
`G.entry(i, j)` from the main thread -- the gate's cluster check, the CG
preconditioner setup of the floor scan, the nonlinear energy-Newton exact-
diagonal preconditioner -- recomputed the near blocks that lived only in the
fill workers' thread-local caches: 265 ms (Duffy) or 77 ms (block-wise) per
block.  That is why the first example-6 gate's cluster check took 2785 s
against 449 s for the block-wise arm (the ratio of the block costs), why its
CG floor "took" 2083 s for 431 iterations (the setup, not the iterations), and
why the nonlinear Q-mag run spent 8545 s of its 10167 s in the Newton loop for
866 inner CG iterations, while the same nonlinear solve on LAB with the shared
cache needed 43 s for 2674 inner iterations (16 ms each, 7 Newton iterations,
no backtracks).  With every BDM1 near pair in the instance-shared cache those
phases collapse to their iteration cost.

Measured on hibino (38 threads, one job at a time, Q-mag linear `mu_r = 1000`,
`results/timing_qmag_*_hibino.json`; the pre-fix row is the 17:18 run of the
same mesh with the previous wheel):

| build | unknowns | Gram wall | near blocks evaluated | near thread-s | shared hits / lookups |
|---|---|---|---|---|---|
| before (affine pairs bypass the shared cache) | 60816 | 1625 s | 137,341 | 36,458 | 6,030 / 27,118 |
| shared cache for every near pair, congruence off | 60816 | 1326 s | 129,536 | 30,261 | 51,094 / 180,630 |
| translation-congruent key | 60816 | **450 s** | 27,062 | 4,078 | 153,442 / 180,504 |
| translation-congruent key, `h = 6 mm` | 153296 | 902 s | 33,273 | 6,018 | 338,165 / 371,438 |
| + dynamic schedule, compute-once, 6 points (section 8.14) | 60816 | **37 s** | 26,807 | 761 | 155,164 / 181,986 |
| + dynamic schedule, compute-once, 6 points, `h = 6 mm` | 153296 | **69 s** | 33,023 | 1,071 | 340,187 / 373,220 |

The field is unchanged to every printed digit (`B_perp(15 mm) = -0.22722 T`
in all three `h = 10 mm` builds).  The routing fix alone removes the worker
duplication (1625 to 1326 s); the congruence key evaluates one near block per
class instead of one per layer (4.8x fewer evaluations on the 6-layer mesh)
and brings the build to 450 s, 7.4 ms per unknown; the 10-layer 6 mm mesh
runs at 5.9 ms per unknown.  Against the TET route's 1.1 ms per unknown the
gap is now about 6x, and the near family is no longer the bulk of it: 4,078
thread-seconds over 38 workers is ~110 s of the 450 s, so the next profile
target is the rest of the build (cluster tree, ACA fills, far blocks).

The pair point count is now chosen per pair.  The rule converges
exponentially when both hosts are affine (unit-cube self-energy `4.9e-9` at
6 points, `1.8e-12` at 8) but only about tenfold per two points on distorted
hosts (tapered sector lattice, entry error against a 10-point reference,
touching / near band: `5.4e-4 / 1.0e-3` at 4, `4.2e-4 / 3.0e-4` at 5,
`7.4e-5 / 8.4e-5` at 6, `6.3e-6 / 7.3e-6` at 8; the band needs the same count
as the touching class).  The first design took `glpair_affine_n` = 6 for
affine-affine pairs and `glpair_n` = 8 for pairs with a distorted host
(`vim.ChargeGram(hex_glpair_n=..., hex_glpair_affine_n=...)`, both published
in `hmat_stats`); the gate and field evidence below then made 6 the default
for every pair, and the two knobs remain for accuracy studies.

The example-6 definiteness gate on the conforming mesh, rerun with the shared
cache and the pair point count forced to 5, 6 and 8 for every pair (hibino,
`results/hex_gram_definiteness_glpair{5,6,8}_hibino.json`; the block-wise
family arm is `results/hex_gram_definiteness_blockwise_family_hibino.json`):

| Gram | build | cluster check | CG at the chi0 floor (431 it) | LOBPCG `lambda_max` | verdict |
|---|---|---|---|---|---|
| first gate, previous wheel, 8 points | 1567 s | 2785 s | 2083 s | 0.99990 | PASSED |
| block-wise near family (no pair rule) | 431 s | 449 s | 218 s | 1.00165 | FAILED (band) |
| shared cache, 8 points | 612 s | 1 s | 5 s | 0.99990 | PASSED |
| shared cache, 6 points | 128 s | 1 s | 5 s | 0.99990 | PASSED |
| shared cache, 5 points | 51 s | 1 s | 5 s | 0.99990 | PASSED |
| shared cache, 6 points, dynamic schedule (section 8.14) | 61 s | 1 s | 5 s | 0.99990 | PASSED |

The post-build phases collapse exactly as the entry-recompute diagnosis
predicts (2785 s to 1 s, 2083 s to 5 s for the same 431 iterations), the
whole gate now takes minutes instead of hours, and the spectrum edge is the
same `0.99990` at 5, 6 and 8 points with CG iteration counts within a few
per cent of each other.  The tapered sector lattice of the near-family test
(the harder distorted case) gives the same generalized spectrum at 8, 6 and 5
points (`lambda_max` 0.99813 / 0.99813 / 0.99812, `lambda_min` at round-off),
and the CEFC 2020 quadrupole field on the `h = 15 mm` mesh (LAB, `mu_r =
1000`, `B_perp(15 mm)`) moves by `4e-6` relative between 8 and 6 points and by
`6e-5` between 8 and 5 (`-0.227134`, `-0.227135`, `-0.227147` T), two orders
below the 0.3 % FEM agreement.  The production default is therefore **6 points
for every pair** (`glpair_n` = `glpair_affine_n` = 6, 2026-09-07); 8 stays an
explicit choice for entry-level accuracy studies, 5 is acceptable on the
evidence but not the default.  On the quadrupole `h = 15 mm` mesh most near
pairs involve a distorted host, so this flip (not the affine-pair rule) is
what brought that build from 615 s to 119 s on LAB.  At 6 points example 6
builds at 2.1 ms per unknown, at 5 points at 0.8 ms, against the TET route's
1.1 ms per unknown: the HEX Gram build is now of the same order as TET.

### 8.14 Parallel efficiency of the build: a static leaf schedule (2026-09-07)

With the near family cut down, the build-phase timers (`build_prep_s`,
`build_cluster_s`, `build_leafgen_s`, `build_fill_s`, `build_diag_s` in
`hmat_stats`, from HACApK `ctl->time[90..92]` and the base build) showed
where the remaining wall time went: on the quadrupole `h = 10 mm` build on
LAB, prep (the self-energy pass `ComputeChargeSigma`) 219 s and the ACA+ fill
658 s of 878 s, cluster tree / leaf generation / diagonal cache below a
second.  Yet the quadrature branches summed to only a quarter of the thread
capacity (hibino: 4,275 thread-seconds against 38 x 448 s; the LAB process
ran on about half its cores).  The cause was the schedule, not the work:
`hacapk_parallel_for` called `ngcore::ParallelFor` with the default task
count, which splits the range into one contiguous chunk per thread -- a
static schedule over a leaf list sorted by row block, whose leaves differ by
five orders of magnitude in cost (a dense near leaf of Duffy pair blocks
against a far low-rank leaf) and whose lower-triangular half is skipped by
the symmetric fill.  The same one-chunk-per-thread split ran the
self-energy pass over host-ordered charges, where the first touch of each
near-block class is the whole cost.  The fix passes 32 tasks per thread to
`ParallelFor` (the runtime pulls tasks from an atomic counter, so many tasks
balance dynamically) in the leaf fill, the self-energy pass and the curved
touch-block precompute.  The shared near-block cache also lost its "racing
first insert wins" design: a per-key slot with `std::call_once` now makes
concurrent misses wait for one evaluation instead of each recomputing the
block (`hex_general_shared_entries` equals `hex_general_shared_misses`,
locked by the congruent-cache test; on the quadrupole the duplication was
five blocks in 26,000, so the gain is the schedule).

Quadrupole `h = 15 mm`, `mu_r = 1000`, 6 points, LAB (8 threads, relative
numbers only): Gram build 119 s -> 62 s, prep 28 s -> 5 s, fill 91 s -> 55 s,
with the quadrature thread-seconds unchanged (382 -> 406 near, 44 far);
450 thread-seconds over 8 workers is 56 s, so the build now runs at about 93 %
parallel efficiency against about 40 % before.  Field unchanged
(`B_perp(15 mm) = -0.227135 T`).

The hibino rerun (38 threads, one job at a time,
`results/timing_qmag_{h10,h6}_mu1000_dynamic_hibino.json`, the two rows added
to the table of section 8.13) combines the dynamic schedule, the compute-once
cache and the 6-point default: the `h = 10 mm` Gram (60,816 unknowns) builds
in **37 s** against 450 s with the congruent cache alone and 1,625 s before
it, the `h = 6 mm` Gram (153,296 unknowns) in **69 s** against 902 s; the
fields agree with the earlier builds to every printed digit
(`-0.22722` / `-0.22724` T).  Per unknown that is 0.61 ms and 0.45 ms, below
the TET BDM2 route's 1.1 ms on the C-type magnet: the HEX Gram build is no
longer the slower route.  The phase timers now read prep 5 s, fill 31 s
(`h = 10 mm`) and prep 5 s, fill 60 s (`h = 6 mm`); the quadrature
thread-seconds (761 + 208 and 1,071 + 837) over 38 workers account for 26 s
and 50 s of those builds, so the remaining gap to perfect balance is under a
third and no longer worth a dedicated pass.  The conforming example-6 gate
under the same wheel builds its Gram in 61 s (first gate 1,567 s) and passes
with the same `lambda_max` 0.99990
(`esrf_three_engine/results/hex_gram_definiteness_dynamic_hibino.json`).

### 8.15 The same magnet on HEX and on TET (2026-09-07)

The per-unknown figures of section 8.14 compare different magnets.  The
question that matters is the same magnet at the same accuracy, so the CEFC
2020 quadrupole was meshed from one Cubit import both ways
(`build_qmag_cubit_mesh.py`, swept HEX and `--scheme tet`, each exported at
curve order 2 or 1) and solved by every route on hibino, one job at a time,
linear `mu_r = 1000` (`quadrupole_cefc2020/results/timing_qmag_*_hibino.json`;
the mixed Omega FEM gives `B_perp(15 mm) = -0.22725 T`):

| route | h [mm] | elements | unknowns | Gram [s] | of which prep [s] | solve [s] | total [s] | ms per unknown | `B_perp(15 mm)` [T] |
|---|---|---|---|---|---|---|---|---|---|
| HEX curved Q2, BDM1 | 15 | 980 | 25,984 | 25 | 3 | 1.3 | **27** | 1.02 | -0.22713 |
| HEX curved Q2, BDM1 | 10 | 2,352 | 60,816 | 37 | 5 | 1.8 | **41** | 0.67 | -0.22722 |
| HEX curved Q2, BDM1 | 6 | 6,040 | 153,296 | 69 | 5 | 4.1 | **76** | 0.50 | -0.22724 |
| HEX curved Q2, BDM1 | 4 | 16,408 | 410,128 | 175 | 5 | 12.2 | **197** | 0.48 | -0.22725 |
| TET straight, BDM1 | 10 | 12,248 | 80,202 | 15 | 0 | 2.3 | **18** | 0.22 | -0.22723 |
| TET straight, BDM1 | 7 | 26,987 | 172,644 | 31 | 0 | 4.6 | **38** | 0.22 | -0.22725 |
| TET straight, BDM1 | 5 | 63,204 | 397,047 | 86 | 0 | 9.5 | **100** | 0.25 | -0.22726 |
| TET straight, BDM2 | 10 | 12,248 | 233,892 | 164 | 0 | 8.0 | **175** | 0.75 | -0.22727 |
| TET curved Q2, BDM1 | 10 | 12,248 | 80,202 | 53 | 40 | 2.3 | **56** | 0.70 | -0.22740 |
| TET curved Q2, BDM2 | 10 | 12,248 | 233,892 | 278 | 204 | 7.9 | **289** | 1.23 | -0.22727 |
| TET curved Q2, BDM2 | 7 | 26,987 | 507,210 | 655 | 453 | 16.3 | **678** | 1.34 | -0.22728 |

Three readings.  First, against the production TET route on a curved mesh
(curved BDM2, the route of the C-type three-engine campaign) the HEX route is
now the faster one at equal accuracy: 40 s against 289 s at `h = 10 mm`, both
within 0.01 % of the FEM.  Second, the curved TET routes pay 70-80 % of their
build in the curved touching-block precompute (`PrecomputeCurvedTouchBlocks`:
204 s of 278 s at BDM2, 40 s of 53 s at BDM1), and that rule is also the less
accurate one -- curved BDM1 at 10 mm lands 0.07 % off where straight BDM1 on
the same tets lands 0.01 % off, and the curved BDM2 harmonics converge from
further away (quadrupole README, "Multipole convergence").  The cheaper and
more accurate curved touching family is the TET route's next lever.  Third,
straight TET BDM1 with its closed-form near integrals remains the cheapest
route per unknown (0.2 ms against 0.5-0.7 ms for the HEX Duffy family) and
in wall time (18 s against 40 s at 10 mm, 100 s against 197 s for 400k
unknowns), while the HEX mesh keeps the quadrupole symmetry exactly and
resolves the pole face with fewer unknowns.  The HEX goal -- TET-class
performance -- is met against the curved TET route and within a factor of
two of the straight one; closing that factor would take a closed-form inner
integral for affine HEX pairs, which is where the remaining Duffy cost sits.

### 8.16 ESRF example 6: the released nonlinear three-engine result (2026-09-08)

The coil-driven quadrupole ran from the production-candidate wheel on mdx1,
one job at a time, through the tracked runner
(`esrf_three_engine/results/case6_nonlinear_three_engine_mdx1.json`; the
per-element warm-start arrays are omitted from the committed copy, which
records the SHA-256 of the complete artifact).  All three nonlinear
formulations converged on their own meshes and one shared mesh-free coil
source:

| engine | unknowns | wall [s] |
|---|---|---|
| HDiv-MMM, BDM1, iron-only HEX | 62,192 | 69 |
| HCurl reduced-A, Periodic Kelvin BDDC | 705,838 | 745 |
| mixed total/reduced Omega, Anderson depth 2 | 223,676 | 1454 |

| pair | core RMS (27 points) | full stencil (45 points) |
|---|---|---|
| hdiv mmm against reduced a | 0.34 % | 0.68 % |
| hdiv mmm against mixed total reduced omega | 0.70 % | 1.20 % |
| reduced a against mixed total reduced omega | 0.78 % | 1.41 % |

The acceptance limit is 3 % on the core stencil and the maximum is
0.78 %.  The run also records the process peak working set,
6.4 GB, which is the memory evidence the C-yoke row asked for on a
comparable problem.  Two readings: the three formulations agree on an
iron-dominated nonlinear quadrupole to under 1 %, and HDiv-MMM reaches that
agreement with 62,192 unknowns on the iron alone against 705,838 for
reduced-A and 223,676 for the mixed Omega route, in 69 s against
745 s and 1454 s for the two FEM routes.  The reduced-A Kelvin BDDC replay row is closed by the same run.

