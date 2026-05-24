"""
Applications: motor optimization, magnet design, induction heater
shape, accelerator pole face optimization.
"""

MOTOR_OPTIMIZATION = """
# IPM motor shape optimization (Gangl 2015 case study)

## Problem
Interior Permanent Magnet (IPM) brushless motor:
- Stator (outer) with coil slots
- Rotor (inner) with permanent magnets in iron cavities
- Thin air gap between them
- 2D cross-section model (typical for rotating-machine analysis)

Goal: smooth out the **cogging torque** by tuning the geometry of
a specific iron region of the rotor.

The cogging torque comes from interaction between PM stray field
and the variable reluctance as the rotor rotates.  Reducing it
improves motor smoothness and reduces audible noise.

## Mathematical formulation

Cost functional:
    J(Ω, u) = ∫_Γ₀ |B_r(u) - B_d|² ds

where Γ₀ is a circle inside the air gap, B_r is the radial
component of B = curl A (= ∇u × ẑ in 2D), and B_d is a TARGET
sinusoidal profile that gives zero cogging.

State equation:
    -div(ν(|∇u|²) ∇u) = f^PM    in D
    u = 0                        on ∂D

with f^PM the magnetization source term from the permanent magnets,
and ν the reluctivity (highly nonlinear in iron, constant in air).

Design domain Ω ⊂ Ω^ref ⊂ Ω_f^ref:  a small piece of the rotor
iron core (the user can only modify SOME of it; the rest is fixed
by manufacturing constraints).

## Constraint set

    O = {Ω ⊂ Ω^ref open, Lipschitz with bounded Lipschitz constant L_O}

Bounded Lipschitz constant prevents oscillation / chattering during
optimization.

## Existence of optimal shape

Gangl 2015 §3 proves: under Assumption 1 (β monotone + Lipschitz),
the minimization problem has a solution Ω* ∈ O.

The proof uses Hausdorff convergence of characteristic functions
and strong-H¹ convergence of u_n → u* (the latter requires a small
ε-regularity gain over H¹).

## Numerical results (Gangl 2015 §6)

Starting from an initial rotor design:
- ~100 iterations of gradient-based shape opt
- Cost J reduced by ~10× (1 → 0.1)
- Resulting Ω: smooth boundary with a characteristic "bump" that
  cancels the cogging Fourier modes

## Practical NGSolve workflow

```python
from ngsolve import *
mesh = generate_motor_mesh(...)

# Nonlinear reluctivity (B-H curve)
def nu_fn(absB):
    # Frohlich-type interpolation of B-H data
    return mur_air if absB < B_threshold else nu_curve(absB)
nu = CoefficientFunction(...)

# State equation
fes = H1(mesh, order=2, dirichlet="outer")
u = fes.TrialFunction(); v = fes.TestFunction()
a = BilinearForm(fes)
a += nu(absB) * grad(u) * grad(v) * dx
# (Newton's method for the nonlinear solve)

# Adjoint state
# (linear in p, after fixing u from the nonlinear state equation)

# Shape gradient (volume form)
# (loop over elements, accumulate ∇V contributions)

# Riesz-lift shape gradient to a smooth deformation field
W = VectorH1(mesh, order=2)
# ... solve auxiliary problem

# Mesh deformation
mesh.SetDeformation(W * step_size)
```

## Connection to Radia practice

Radia is used for:
- Accelerator magnets (where the cost is field uniformity or
  multipole content along beam axis — analogous to motor cost)
- Insertion devices (wigglers/undulators where the cost is
  spectral content of trajectory)
- Iron pole face design (cost = field flatness over a region)

For all these, the Gangl-Sturm framework applies directly.  Radia
itself doesn't run the optimization — that happens in NGSolve or
similar — but Radia gives the field-evaluation backbone for
validation of the optimized geometry.

## Connection to Radia + Mathematica + radia-mcp

Workflow for a PhD student doing magnet design:
1. Theory layer (radia_mcp.topology_optimization) — read this
   knowledge module
2. Symbolic layer (radia_mcp.mathematica) — derive shape-derivative
   for the specific cost functional in closed form
3. Discretization layer (radia_mcp.radia_ngsolve) — implement in
   NGSolve, run the optimization
4. Validation layer (Radia C++) — verify the optimized field
   against Radia's exact analytical formulas

This is the 4-layer stack from the killer demo in the README.
"""


def get_applications_documentation(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
      "all"   - All applications
      "motor" - IPM motor shape optimization (Gangl 2015 case study)
    """
    if topic.lower().strip() in ("motor", "ipm", "all"):
        return MOTOR_OPTIMIZATION
    return f"Unknown topic '{topic}'. Available: all, motor."
