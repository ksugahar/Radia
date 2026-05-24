"""
Shape optimization for nonlinear magnetostatics — the velocity method,
shape derivative, and adjoint-state machinery for FEM constraints.
"""

SHAPE_OPT_OVERVIEW = """
# Shape optimization in computational EM: overview

The goal of shape optimization is to find the geometry Ω that
minimizes a cost functional J(Ω, u) subject to a PDE constraint
(e.g. Maxwell's equations).

## The setup

- Cost functional:        J(Ω, u) = ∫_Γ |B_r(u) - B_d|² ds
  (tracking-type: match desired B-radial profile along a circle
   Γ in the air gap of an electric motor)
- State equation:         -div(β(x, |∇u|²) ∇u) = f  in Ω, Dirichlet on ∂Ω
  (nonlinear magnetostatics; u = 3rd component of magnetic vector
   potential; β = magnetic reluctivity)
- Design domain:          Ω ⊂ Ω^ref ⊂ Ω_f^ref (part of iron core)

The optimization problem:
    minimize J(Ω, u(Ω))
    Ω ∈ O = {Ω ⊂ Ω^ref open + Lipschitz with bounded Lipschitz constant}

## Why nonlinearity matters
Iron in electric motors saturates: B-H is highly nonlinear above
~1.5 T.  The reluctivity ν(|B|²) is not constant — it depends on
the field value at each point.  This makes the shape derivative
computation MUCH harder than the linear case, because the
material-derivative of u (which classical methods need) is itself
defined by a nonlinear PDE.

## The Gangl-Sturm 2015 fix: novel Lagrangian for nonlinear problems

Instead of computing the material derivative of u directly,
introduce a Lagrangian
    G(Ω, ϕ, ψ) = J(Ω, ϕ) + ⟨state eq.(ϕ), ψ⟩
where ϕ, ψ are test functions in H¹₀(D).  The state equation is
recovered from dψG = 0 (variation w.r.t. ψ gives the constraint).
The adjoint equation comes from dϕG = 0.

Key result (Theorem 3 in Gangl 2015, due to Sturm): the shape
derivative can be computed as the **partial time derivative** of
the Lagrangian at the perturbed solution, WITHOUT ever computing
the material derivative of u.  Formally:
    dJ(Ω; V) = ∂_t G(t, u_t, p_0) |_{t=0}
where p_0 is the solution of the AVERAGED adjoint equation.

This works for nonlinear PDE constraints, generalizing the
classical (linear) Sokolowski-Zolesio framework.

## Why this is a big deal
- Before Gangl-Sturm: nonlinear shape opt either skipped the
  derivative (heuristic methods, slow convergence) or used
  formal computations without rigorous justification.
- After Gangl-Sturm: standard gradient-based optimization works
  even for nonlinear PDEs.  Convergence is fast and provable.

For motor optimization with iron-core saturation, this is the
mathematical foundation.

See `shape_derivative` topic for explicit formulas, `motor` topic
for the concrete application.
"""

SHAPE_DERIVATIVE = """
# Shape derivative via the velocity method

## The velocity method (Sokolowski-Zolesio, Allaire)

A domain Ω is deformed by the flow T_t of a velocity field V:
    d/dt x(t) = V(x(t)),  x(0) = X
    T_t(Ω) = {x(t) : X ∈ Ω}

For small t, T_t is a smooth diffeomorphism.  The shape derivative
is defined as
    dJ(Ω; V) = lim_{t → 0+} [J(T_t Ω) - J(Ω)] / t.

For a smooth velocity field V ∈ D¹(D, R²) compactly supported in D,
this limit exists and is linear in V.

## Volume vs boundary expression

The shape derivative can be expressed in two forms:

**Boundary form** (Hadamard formula): for the linear case
    dJ(Ω; V) = ∫_∂Ω g(s) · V · n ds
This requires the normal trace of V on ∂Ω.

**Volume form**: for general (including nonlinear) cases
    dJ(Ω; V) = ∫_D [a(x) · ∇V + b(x) tr(∇V)] dx
This involves volume integrals of ∇V over all of D, NOT just ∂Ω.

For numerical algorithms, the volume form is preferable: it can be
evaluated on the same FEM mesh without surface integration.

## Material derivative and shape derivative

For a field u(Ω) depending on Ω, the **material derivative** is
    ϕ' = lim_{t→0+} [u(T_t Ω) ∘ T_t - u(Ω)] / t
i.e., compare u in transported coordinates.

The **shape derivative** of u is
    u' = ϕ' - V · ∇u
(subtract the transport effect).

For shape-derivative computation of functionals, ϕ' is usually what
appears, NOT u'.  This is what makes the Lagrangian approach
attractive — ϕ' is harder to compute for nonlinear u, but the
Lagrangian sidesteps it entirely.

## Shape Lagrangian for cost-with-PDE-constraint

J(Ω) = ∫_Γ |B_r(u) - B_d|² ds  with state eq.  div(β(|∇u|²) ∇u) = f.

Define
    G(Ω, ϕ, ψ) = ∫_Γ |B_r(ϕ) - B_d|² ds + ∫_D β(|∇ϕ|²) ∇ϕ · ∇ψ dx
                  - ⟨f, ψ⟩

State equation:  dψG = 0 (gives u)
Adjoint equation: dϕG = 0 (gives adjoint variable p)

Shape derivative:
    dJ(Ω; V) = ∂_t G(t, u_t, p_0) |_{t=0}

where p_0 solves the AVERAGED adjoint equation
    ∫_0^1 dϕG(t, su_t + (1-s)u_0, p_t; ϕ̂) ds = 0  ∀ ϕ̂

The averaging over s is what allows the nonlinearity to be handled
rigorously (Sturm 2015).

## Strong-form of the adjoint equation

For the nonlinear elliptic case with β(|∇u|²):
    -div(A₁(∇u) ∇p) = 0   in Ω
    -div(A₂(∇u) ∇p) = 0   in D \\ Ω
    p = 0                  on ∂D
with appropriate transmission across Γ_f.

Here A_i is the Jacobian of β · ∇u w.r.t. ∇u:
    A_i(∇u) = β_i(|∇u|²) I + 2 β'_i(|∇u|²) (∇u) (∇u)^T

A_i is symmetric positive-definite under Assumption 1.4 of Gangl 2015
(equivalent to monotonicity of  s ↦ β(s²) s).

## In practice: gradient-based shape optimization

Algorithm:
1. Initialize Ω
2. Solve state equation for u
3. Solve averaged adjoint equation for p
4. Compute dJ(Ω; V) for a basis of velocity fields V
5. Choose descent direction V* (typically by solving an auxiliary
   variational problem to extract a smooth deformation)
6. Update Ω → T_{tV*} Ω for small step t
7. Repeat until ||dJ|| is small

The whole loop in NGSolve takes ~few minutes per iteration for a
medium-scale motor problem.
"""

NUMERICAL_ALGORITHM = """
# Numerical algorithm: gradient-based shape optimization

## Algorithm structure (Gangl 2015 §6)

```
Input: initial shape Ω_0, max iterations N
For k = 0, 1, ..., N:
  1. Mesh Ω_k (re-mesh or move existing mesh)
  2. Solve nonlinear state eq.  ν(|∇u_k|²) curl-curl u_k = f
     using Newton's method (10-30 iterations)
  3. Solve averaged adjoint eq. for p_k (linear, 1 solve)
  4. Compute shape gradient ∇J|_Ω_k as a vector field on D
  5. Smooth the gradient → descent direction V_k
  6. Line search for step size τ_k:  τ_k = argmin J(T_{τV_k} Ω_k)
  7. Update Ω_{k+1} = T_{τ_k V_k} Ω_k
  8. If ||∇J|| < tol, stop
End for
```

## Velocity-field smoothing

Direct application of the shape derivative gives V as a 1-form
(distribution on ∂Ω).  For mesh movement, V needs to be a vector
field on the WHOLE mesh (including air, magnets, etc.).  Use a
H¹-Riesz lifting:
    Find W ∈ H¹(D)²  such that  ∫_D (∇W : ∇W' + W·W') = dJ(Ω; W')
    ∀ W' ∈ H¹(D)²
This gives a smooth deformation field that vanishes far from ∂Ω.

## Convergence

Numerical experiments (Gangl 2015 §6):
- IPM motor: 50-200 iterations to converge
- Cost function reduces by 1-2 orders of magnitude
- Resulting shape has smooth boundary (no jagged contour)

Comparison with topological derivative methods:
- TD: faster initial reduction but jagged shape
- SD: slower but smoother final shape
- Best: hybrid (TD for global topology, SD for boundary smoothing)

## NGSolve implementation pattern

NGSolve has built-in support for:
- Nonlinear FEM (Newton solver for magnetostatics with B-H curve)
- Shape derivative computation via `ngsxfem` or hand-derived adjoints
- Mesh moving via `mesh.SetDeformation(gfu)`

A typical motor-opt script with NGSolve is ~200-500 lines.

## Connection to Radia / Radia C++

Radia itself doesn't have shape-opt machinery, but the FEM-side
mesh moving + NGSolve nonlinear solver + Radia C++ for the
fast-multipole field evaluation (e.g. for very large magnet
assemblies where FEM is too slow) is a natural pipeline.

A potential workflow:
1. Coarse shape opt with NGSolve FEM on a simplified motor
2. Use Radia C++ for the final field validation (rapid, no remesh)
3. Iterate
"""


def get_shape_optimization_documentation(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
      "all"
      "overview"        - Why nonlinear shape opt is hard, Gangl-Sturm 2015 fix
      "shape_derivative" - Velocity method, Hadamard formula, Lagrangian-based
                          derivation for nonlinear PDE constraints
      "algorithm"       - Gradient-based algorithm structure, NGSolve patterns
    """
    topic = topic.lower().strip()
    if topic == "overview":
        return SHAPE_OPT_OVERVIEW
    if topic in ("shape_derivative", "velocity_method", "hadamard"):
        return SHAPE_DERIVATIVE
    if topic in ("algorithm", "numerical", "ngsolve"):
        return NUMERICAL_ALGORITHM
    if topic == "all":
        return "\n\n".join([
            SHAPE_OPT_OVERVIEW,
            SHAPE_DERIVATIVE,
            NUMERICAL_ALGORITHM,
        ])
    return (
        f"Unknown topic '{topic}'. Available: "
        "all, overview, shape_derivative, algorithm."
    )
