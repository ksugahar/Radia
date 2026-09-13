# ESRF6 energy-Newton repair candidate

Status: focused tests, small native checks, and the actual ESRF6 BDM1/mass-Riesz
nonlinear residual gate passed. Other paths and three-method acceptance remain HOLD.
The preceding candidate's three-engine
`passed=true` is not accepted: 34 Armijo backtracks exhausted the line search,
but the tiny rejected step was accepted as convergence after one iteration.

This repair makes the inverse BH energy the exact antiderivative of the same
PCHIP H(M) used by the residual; its derivative supplies the radial tangent.
The saturation extension has continuous H and C1 energy, not a continuous
second derivative at its junction. Nonpositive zero-field tangents are rejected.

Material energy, weak load, and tangent now use the same NGSolve
IntegrationRuleSpace samples and integration rules. The order is
`max(6, 2*fes.globalorder+2)`; this ensures discrete derivative consistency,
not physical quadrature convergence. There is no element-average substitution
for high-order magnetization. NGSolve owns interpolation, geometry and FE assembly.

NGSolve 6.2.2606 tensor IntegrationRuleSpace evaluation in the symbolic SIMD
bilinear-form path terminated the local test process. The material integrator
therefore explicitly uses `simd_evaluate=False` and sets each element type's
integration rule. This is not a silent fallback or a claim of optimal speed.

Armijo accepts only evaluated finite descending steps, retries a stale tangent
once, and raises on exhaustion. Neither a tiny step nor a settled-step counter
certifies convergence. Every accepted stage satisfies the actual assembled
nonlinear residual divided by a fixed stage-source norm (initial residual when
the source is zero). A zero-source, zero-initial-state solve skips native warmup.

Focused checks cover PCHIP derivatives, TET/HEX/WEDGE orders 1 and 2, affine and
quadratically deformed geometry, distinct materials, zero-field Hessians, and
the production outer Newton function with a dense SPD demag test double.
Those tests do not validate the C++ inner solver or the full ESRF6 magnet.
An independently identified wheel and Hibino run are required before removing
the existing nonlinear acceptance HOLD. Mixed Omega source-load assembly
quadrature convergence remains a separate task.

## Fable 5.1 Review Coordination (2026-09-13)

Finding A independently identifies the exhausted-line-search false convergence
and the element-average energy/residual inconsistency repaired in `4d85e72cc`.
Do not build a second competing Armijo patch from the pre-repair source. Review
that repair and its tests instead. The isolated wheel used for native smoke is
identified in `results/candidate_4d85e72cc/`; the old `95721799` candidate remains
comparison evidence, not the final repair candidate.

There are 51 focused checks and 12 installed-wheel native small-mesh checks.
Six strong-drive native cases (TET/HEX/WEDGE, orders 1/2) require actual Newton
updates and finish with true relative residual at most `1.59e-5` for tolerance
`2e-5`. These are not the ESRF6 acceptance run. Keep three claims separate:
surface overlap not detected in the audited meshes; actual ESRF6 nonlinear
convergence must be established per path (BDM1/mass-Riesz now passed below);
actual ESRF6 three-method numerical acceptance HOLD.

Finding B is distinct: in the measured `4d85e72cc` wheel,
`project_source_total_hodge` assembles with default `dx`, while the mixed load
has an explicit `bonus_intorder`. A matching
quadrature contract needs a focused implementation and test. Small action-scaled
algebraic residuals do not prove physical source/load quadrature convergence.
Near-zero block RHS values must not make RHS-relative normalization the only
acceptance measure. Projection/test spaces and the gauge term also matter when
interpreting discrete orthogonality; do not assume exact zero from naming alone.

Finding C concerns diagnostic integration cost and heap behavior, not evidence
that a solver has converged. Any audit optimization must preserve the evaluated
functional and record the original error, NGSolve version, threads, and memory.

The execution task owns the fixed-input ESRF6 HDiv-only comparison after the
remaining CAD/coil/gap checks. Independent Fable review should use this baseline;
projection/load and audit changes must be coordinated before overlapping edits.
Do not use `Newton iterations > 1` as a universal acceptance rule: a good warm
start can already meet tolerance. Require the true nonlinear residual and no
accepted exhausted-line-search step. Iteration/backtracking counts remain
diagnostic evidence for this particular failing example.

### Follow-Up Implementation

The diagnostic follow-up names `tolerance`, `iteration-limit`, and
`line-search-exhausted` outcomes and reports the residual target. Exhaustion
attaches the statistics and rejected Newton correction/decrement to the error.
An undefined relative correction at zero magnetization is null, with the
absolute correction and state norms retained. No settled-step success is restored.

The shared HDiv runner requires a finite true nonlinear residual at its requested
tolerance. New checkpoint modes additionally require the recorded residual and
target; the word `tolerance` alone cannot certify a solve. Legacy checkpoints
without this mode are conservatively rejected at 33 or more aggregate backtracks.
That legacy rule is a quarantine heuristic, not a theorem that 33 accumulated
backtracks always mean exhaustion. Both committed ESRF6 and CEFC J3.0 false
convergence records are explicit regression fixtures. New successful results
are not rejected solely for a large accumulated backtrack count.

These diagnostic additions are not in the running `4d85e72cc` wheel. Its first
actual ESRF6 BDM1/mass-Riesz solve converged in eight Newton iterations, with one
backtrack and true residual 3.278340250386088e-7 at tolerance 2e-5. This isolates
the prior false-convergence issue; BDM2, Picard-route consistency and nonlinear
three-method acceptance remain separate. Gap-field agreement alone is not
nonlinear acceptance. The diagnostic runner retains all cell-average M vectors
for iron-side comparisons, not just the near-zero global average of a quadrupole.
Those averages are not a complete high-order M norm or a peak-field certificate.

### Same-Wheel FEM Rerun And Execution Correction

The execution task re-solved reduced-A and mixed Omega using the same
installed `4d85e72cc` wheel as the accepted BDM1 result. The FEM adapter uses
order 2, BDDC/CG with relaxation 0.1 for reduced-A, and mixed Omega with source
order 2, bonus 4, projected Kelvin trace, relaxation 0.3 and Anderson depth 2.
Both use tolerance 2e-5 and cap 200. The driver records full Python hashes,
native/wheel/input hashes, BH/source identity, and all sampled vector fields.
No previous FEM checkpoint certifies this rerun.

Fable should review A against `4d85e72cc` and its actual BDM1 result, and B
against the still-distinct projection and load quadrature contracts. Coordinate
overlapping repairs and remote jobs with the execution task. This document is
the shared handoff; no direct Fable tool destination is available in this task.

An earlier execution report incorrectly labelled the Windows venv launcher and
its base-Python child as duplicate jobs. The execution task stopped that child,
interrupting three Picard-warmstart attempts; those attempts are not solver
failures or successful evidence. Claims that another session auto-launched
duplicates are withdrawn. A subsequent process audit established PID 13800
(venv) -> PID 10056 (base Python) as a normal parent/child pair for the FEM run.
Future process ownership checks include PID, parent PID, creation time and run
arguments, not command-line equality alone. The completed mass-Riesz JSON
was not affected and retains SHA-256 `18a5b1c97d28f738191d4c40d1aee1a5d19997ac22824526b4fc3286d2d553ce`.

`nl_tol` is now a residual/source-load ratio for 3D energy-Newton, not relative
step size. The BDM1 result demonstrates the difference: residual 3.27834e-7
passes 2e-5, while the final relative step 2.62568e-4 would not. API documentation
and the unreleased changelog record this change. FEM/Picard stopping measures
must be interpreted under their own contracts.

The updated deep-saturation sphere regression passed on the same installed
wheel on hibino (one test, 85.48 s, two threads). It now checks the actual
nonlinear residual against 1e-6 and the final convergence flag, rather than
inferring success from returning without an exception. Its JUnit evidence is
in `S:/Radia/validation_artifacts/esrf6_mesh_audit_20260913/deep_saturation.xml`
(SHA-256 `96537d88f3f5c2cb3e98f8612335bc998f05fd2e0517fb5a8d2abef04b6f60fc`).
This is a small sphere regression, not strong-saturation ESRF6 certification.

Separate follow-up `804dcaac4` forwards the mixed `bonus_intorder` to the total
Hodge projection, records it in source diagnostics, and tests the gauge-aware
weak orthogonality with matching rules at bonuses 0 and 6. It is not installed
in the running candidate. The focused source suite passed 149 tests; this does
not establish source-quadrature convergence for the full magnet. Keep this
source repair separate from the provenance and scope of the candidate results.

### Completed Nominal BDM1 Acceptance

The same-wheel FEM rerun completed successfully. Core 27-point relative vector
RMS differences, with the left field as denominator, are HDiv/reduced-A
0.334891%, HDiv/mixed Omega 0.669467%, and reduced-A/mixed Omega 0.781827%.
All pass the predeclared 3% field-agreement gate. Reduced-A converged in 31
iterations (relative change 1.859221e-5); mixed converged in 14 (relative B
change 7.565560e-6), both at tolerance 2e-5. The repaired HDiv solve used eight
Newton iterations and an assembled residual of 3.278340e-7. FEM change criteria
are not equation residuals or true-error bounds.

See `results/candidate_4d85e72cc/README.md` for the complete scope, hashes,
recovery manifest and raw 45-point results. This closes nominal BDM1 acceptance,
not whole-method BDM2/IMA/deep-saturation qualification or release approval.
The local-built wheel must not be labelled CI-built. No tags or release are
created by the execution task. A JSON replay test checks source/BH identity,
convergence and independently recomputes the accepted field differences.
