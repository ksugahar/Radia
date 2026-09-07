# Two-layer optimization diagnostics

`radia_mcp.optimization` owns solver-neutral guidance and diagnostics;
`topology_optimization` owns electromagnetic shape/topology knowledge and
PDE/adjoint-specific evidence. Both compose into the existing `radia-design`
optimization profile. There is no extra server command or client configuration.
The learning-only profile is unchanged; PINN/assimilation-specific extensions
remain future work.

## Implemented

- `optimization_guide`: ownership, explicit scale convention, and limitations.
- `optimization_stopping_audit`: smooth unconstrained, first-order diagnostics.
- `optimization_gradient_check`: agreement with caller-supplied independent
  gradient references at two or more distinct perturbations.
- `optimization.line_search.armijo_backtracking`: canonical home of the existing
  implementation. The old `topology_optimization.line_search` import re-exports
  the exact same function; behavior/signature and existing numerical tests stay.

Use fixed positive characteristic variable scales S_i and an objective scale q.
For x=S z and F=f/q, the dimensionless gradient is S grad(f)/q. Diagnose the
infinity norm of that gradient, max_i |delta x_i|/S_i, and |delta f|/q separately.
Choose scales before the run based on meaningful physical magnitudes, not the
current objective (which can be zero). Unit conversion must transform the
gradient by the chain rule too. The diagnostics do not infer these scales.

A small step and objective change with a large gradient is stagnation, not
convergence. All small plus non-increasing objective yields only a
`first_order_candidate`. No curvature information is supplied, so a saddle or
maximum can also satisfy the criterion. This is deliberately not a global or
local minimum certificate. Bounds, general constraints and nonsmooth problems
are rejected rather than silently judged by an unconstrained gradient.

The gradient checker does not execute a submitted objective or certify that
the caller's references are independent. For central differences evaluate at
x +/- h S_i e_i and divide by 2 h S_i, then submit the physical gradient for
each dimensionless h. Compare multiple step sizes at the same point and repeat
at nonstationary points. Truncation, cancellation and solver noise can cause
disagreement. Complex-step requires an analytic complex evaluation path; it is
invalid for clipping/abs and many branch-dependent solvers. Agreement is not a
proof of derivative correctness.

## Source and interpretation

Kanamori, Suzuki, Takeuchi and Sato, *Continuous Optimization for Machine
Learning* (Kodansha, 2016), section 3.3, printed pp. 52-54, and section 4.1.
Section 3.3 discusses stopping criteria and the effect of scaling.
The source emphasizes joint observation of gradients, iterates and objectives.
Our fixed characteristic-scale diagnostic is an engineering adaptation, NOT an
implementation of the book's exact machine-epsilon-based formula or convergence
theorem. Multi-step derivative comparison is additional verification policy.
Do not copy the book's prose or illustrations into the package. No full text is
distributed and no model-weight training is performed.

## Next slices, deliberately not claimed complete

1. Constraint-qualification and curvature evidence beyond residual diagnostics;
   keep constrained stationarity distinct from sufficiency.
2. Audit remaining global-search helper ownership without merging Bayesian,
   evolutionary and local smooth search merely because all optimize objectives.
3. Proximal methods/ADMM and structured regularization; connect verified cases
   to inverse magnetic-field/coil design rather than invent generic solvers.
4. Keep Bayesian/evolutionary search distinct; use local refinement only when
   smoothness, derivative accuracy and constraints permit it.

Public ecosystem solvers remain preferred. These MCP diagnostics are not new
Radia numerical kernels or a new MATLAB optimizer implementation. The known
legacy Armijo behavior is preserved, including returning the last rejected
trial with accepted=false; callers must not accept that trial unconditionally.

## Constrained diagnostics

`optimization_kkt_audit` handles smooth minimization with g(x)<=0 and h(x)=0.
Provide all constraint values and Jacobian **rows** at the same point, including
bounds if they exist. The Lagrangian convention is L=f+mu*g+nu*h: physical
inequality multipliers mu must be nonnegative; equality multipliers nu have no
sign restriction. The tool checks supplied multipliers, never estimates them.
Use empty lists for an absent family and at least one constraint overall.

With positive fixed scales x=S*z, F=f/q, G=g/c and H=h/d, dimensionless
multipliers are lambda=mu*c/q and eta=nu*d/q. The reported infinity norms are:

- primal: max(G,0) and H, separately;
- dual: min(lambda,0);
- complementarity: lambda*G;
- stationarity: grad_z(F)+J_G^T*lambda+J_H^T*eta.

Each check uses the supplied positive dimensionless tolerance. Infeasibility
takes status priority; all residuals within tolerance yield only
`first_order_candidate`. Constraint qualifications are explicitly
`not_verified`. Degenerate constraints may prevent KKT necessity even at an
actual minimum. Nonconvex maxima/saddles may satisfy KKT. Smoothness, constraint
completeness, Jacobian accuracy, convexity and second-order sufficiency are not
certified by caller-supplied residuals.

`optimization_projected_gradient_audit` handles smooth box-only problems.
Bounds are physical coordinates; null denotes no bound. Equal bounds fix a
variable. It reports bound infeasibility separately from
G=z-project_box(z-grad_z(F)), using a **fixed unit dimensionless step**. This
mapping has the correct boundary sign behavior even when the raw gradient is
nonzero. It is computed in displacement coordinates to avoid subtracting large,
nearly equal iterates. General constraints must use the KKT route, not this
box diagnostic. Tolerances and fixed scales still affect finite-precision
classification; a small mapping is not an accuracy or optimality certificate.

Source: Kanamori et al. (2016), section 3.2 (feasible directions and convex
sufficiency) and section 10.1 (active inequalities and KKT conditions).
Refer to these sections in the book cited above; the book is not bundled.
The dimensionless residual
contracts and box mapping are engineering adaptations, not source algorithm
transcriptions. These are read-only MCP diagnostics, not numerical solver or
MATLAB kernel additions. No new client process is required.

## Helper ownership and compatibility

The existing NumPy-only helpers now have explicit canonical homes:

| Capability | Canonical module | Compatibility module |
| --- | --- | --- |
| Nonlinear least squares (LM) | `optimization.nonlinear_lsq` | `topology_optimization.nonlinear_lsq` |
| Regularized linear inverse (TSVD/Tikhonov/L-curve) | `optimization.linear_inverse` | `topology_optimization.linear_inverse` |
| Linear CG for SPD systems | `matrix_solvers.krylov` | `topology_optimization.krylov` |

All paths are relative to `radia_mcp`. Compatibility modules re-export the
same function objects, not forked copies or alternate solvers. Signatures remain
unchanged. No new MCP tool, server, dependency or MATLAB numerical implementation
is introduced. Linear CG's algorithm is unchanged; it is not nonlinear CG.

The caller audit found direct imports in the existing helper regression tests
and references in topology application and stream-function knowledge. Those
tests deliberately retain the old imports and assert identity with the new
homes. Application guidance names the canonical paths. Field-map construction
and PDE adjoints remain domain-owned. The nonlinear multi-start evidence gate
and simplex-stationarity gate retain their existing summary contracts and tool
names in `topology_optimization`; their ownership needs a separate gate audit,
not a move bundled with these numerical helpers.
Global/evolutionary search is not migrated in this slice.

Two numerical edge cases were corrected during migration:

- LM now reports `grad_norm` at the returned iterate. `converged` is true only
  when its absolute gradient criterion holds; small steps, exhausted iterations
  and no improvement alone are not success. The additive `termination_reason`
  field distinguishes these exits. Invalid controls, callback shapes and
  nonfinite evidence fail loudly. Callbacks must be deterministic. Existing
  residual-scaling limitations still apply; gradient convergence is not a
  minimum certificate. This intentionally tightens the legacy success flag.
- Positive-lambda Tikhonov uses s/(s^2+lambda^2) without dividing by s, so exact
  null modes yield zero rather than NaN. Lambda zero uses NumPy least squares
  with `rcond=None`; filter factors for exact zero modes are zero. Negative or
  nonfinite lambda is rejected. TSVD refuses a selected exact-zero singular
  mode and asks for a smaller k. Numerical-rank/truncation choice remains the
  caller's responsibility. L-curve corner selection remains a discrete heuristic.

These compact teaching helpers are not replacements for established production
optimization libraries. Existing coil-fit, analytic and SciPy comparison tests
remain authoritative regressions; no new physics solve is needed for the move.
