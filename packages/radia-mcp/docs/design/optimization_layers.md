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
The supplied Japanese PDF's pp. 64-65 (one-based PDF numbering) explain
stopping/scaling; the rendered p. 65 was checked because OCR mangles equations.
The source emphasizes joint observation of gradients, iterates and objectives.
Our fixed characteristic-scale diagnostic is an engineering adaptation, NOT an
implementation of the book's exact machine-epsilon-based formula or convergence
theorem. Multi-step derivative comparison is additional verification policy.
Do not copy the book's prose or illustrations into the package. No full text is
distributed and no model-weight training is performed.

## Next slices, deliberately not claimed complete

1. Constraint-qualification and curvature evidence beyond residual diagnostics;
   keep constrained stationarity distinct from sufficiency.
2. Audit references/callers before relocating nonlinear least squares and linear
   inverse helpers. Linear CG belongs to linear algebra; nonlinear CG belongs
   to optimization. Do not relocate both merely because their names match.
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
sufficiency) and section 10.1 (active inequalities and KKT conditions), inspected
in the supplied PDF on PDF pages 62-63 and 170-173. The dimensionless residual
contracts and box mapping are engineering adaptations, not source algorithm
transcriptions. These are read-only MCP diagnostics, not numerical solver or
MATLAB kernel additions. No new client process is required.
