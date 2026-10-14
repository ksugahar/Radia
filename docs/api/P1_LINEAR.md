# Linear total-A magnetostatics

`radia.p1_linear.solve_p1_linear` solves a linear magnetostatic problem with
lowest-order Nedelec elements and one beta-zero AMS-PCG solve:

\[
\int_\Omega \nu\,\operatorname{curl} A\cdot\operatorname{curl} v\,dx
= \int_{\Omega_J} J\cdot v\,dx,
\qquad \nu=\frac{1}{\mu_0\mu_r}.
\]

The mesh must contain straight three-dimensional tetrahedra. Relative
permeability is a positive finite scalar constant for each named material;
materials omitted from `mu_r_dict` use vacuum permeability. The current is a
real vector coefficient function in A/m², restricted to `current_materials`.
The boundary condition is homogeneous tangential A on `dirichlet`.

```python
import ngsolve as ng
from radia.p1_linear import solve_p1_linear

ng.SetNumThreads(8)  # The caller selects the thread count.
result = solve_p1_linear(
    mesh,
    current_cf=current,
    current_materials="coil",
    mu_r_dict={"iron": 500.0},
    dirichlet="outer",
    cg_tolerance=1e-8,
)
B = result["B_cf"]
stats = result["stats"]
```

Call the function **outside** `ng.TaskManager()`. Like `solve_p1_newton`, it
constructs AMS outside a parallel region and owns scoped assembly/solve
regions. Each call constructs a new system and hierarchy. No factorization or
hierarchy is implicitly reused across calls.

The ungauged system requires a load orthogonal to the discrete gradients on
free dofs. An incompatible current raises; it is never silently projected or
regularized. `radia.meshed_current` supplies suitable closed-conductor current
sources.
The compatibility check currently includes every vertex gradient, inherited
from `p1_newton`. It supports interior current sources and may reject valid
currents reaching a Dirichlet boundary. Restricting the check to admissible
boundary gradients is a separate formulation change.
Nonconvergence raises after an independent original-equation residual
check. The check enforces the requested tolerance and the repository relative
residual limit of 1e-6. A zero source produces a zero solution.

The result contains `A`, `B_cf`, `H_cf`, elementwise `nu`, `fes` and `stats`.
`observation_points` optionally returns `observation_B_T`. Statistics include
the single linear-solve count, CG iterations, true relative residual,
gradient-load compatibility defect and phase times. Source creation before
the call is outside these times and should be measured separately for an
inclusive application time. Field observations are included in the total
when requested.

Use `solve_p1_newton` for a tabulated nonlinear B-H law. The linear entry point
does not infer linearity from a table. The small executable controls are in
`tests/test_p1_linear.py`, including independent NGSolve matrix assembly and a
gauged sparse direct field/energy comparison.

MATLAB batch access is `radia.python.p1Linear("solve_p1_linear", {mesh},
Keywords=options)`. Mesh, field and grid-function objects stay in Python;
this is an explicit Python fallback, not native MEX coverage.
