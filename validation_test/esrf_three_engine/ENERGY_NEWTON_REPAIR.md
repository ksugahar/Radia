# ESRF6 energy-Newton repair candidate

Status: pending native validation. The preceding candidate's three-engine
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
