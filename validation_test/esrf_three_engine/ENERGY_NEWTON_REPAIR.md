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

Finding B is distinct: `project_source_total_hodge` currently assembles with
default `dx`, while the mixed load has an explicit `bonus_intorder`. A matching
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
