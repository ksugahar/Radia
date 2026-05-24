"""
Physics-Informed Neural Networks (PINN) and Gaussian Processes (PI-GP).
"""


# Authoritative topic enum for the dispatcher tool (wired into
# `<short>_topics()` via common.register_topics_tool).
TOPICS: dict[str, str] = {
    "overview": "PI-AI vs FEM, PINN vs PI-GP comparison",
    "pi_gp": "Physics-Informed Gaussian Process foundation (Raissi 2017)",
    "pi_gp_solvers": "PI-GP generalizes classical PDE solvers (Pförtner 2023)",
    "pinn": "PINN formulation, Maxwell applications, state of the art",
    "all": "return \"\n\n\".join([",
}

PINN_OVERVIEW = """
# Physics-Informed AI for PDEs — overview

Classical FEM solves PDEs by:
1. Mesh the domain
2. Choose basis functions (Whitney etc.)
3. Assemble + solve a sparse linear system
4. Interpret the solution

Physics-Informed AI (PINN, PI-GP) replaces this with:
1. Choose a function representation (neural network OR GP prior)
2. Minimize a loss = PDE residual + boundary residual + data fit
3. Read the trained representation as the solution

## When does PIAI win over FEM?

| Situation | Win Method |
|---|---|
| Standard forward solve, simple geometry | FEM |
| Inverse problem (find unknown parameters from data) | PI-AI |
| Multi-fidelity (combine high+low fidelity data) | PI-AI |
| High-dimensional PDEs (Hamilton-Jacobi-Bellman ≥4D) | PI-AI |
| Real-time query needed after one-time training | PI-AI |
| Stiff problems, very nonlinear materials | FEM |
| Need certified error estimates | FEM |

## PINN vs PI-GP

| Aspect | PINN | PI-GP |
|---|---|---|
| Function space | Neural network (∞-dim, expressive) | RKHS (kernel-determined) |
| Training | Stochastic gradient (Adam, L-BFGS) | Maximum marginal likelihood |
| Scaling in N samples | O(N) per epoch | O(N³) without approx |
| Uncertainty | Hard to quantify (ensemble / Bayesian NN) | NATIVE Bayesian uncertainty |
| Nonlinear PDEs | Yes (loss is just the residual) | Hard (linear-only or use sampling) |
| State-of-the-art (2026) | Maxwell/NS/elasticity demos | Linear PDE inverse problems |

For computational EM (mostly linear PDEs), PI-GP has stronger
theoretical foundation (Pförtner 2023 shows it generalizes classical
linear PDE solvers).  PINN is more flexible but requires more tuning.
"""

PI_GP_FOUNDATION = """
# Physics-Informed Gaussian Process (Raissi-Karniadakis 2017)

## The construction

Goal: learn an operator L^φ_x identified by parameters φ from data
{X_u, y_u} on u(x) and {X_f, y_f} on f(x) = L^φ_x u(x).

Step 1: place a GP prior on u
    u(x) ~ GP(0, k_uu(x, x'; θ))

Step 2: linear operators preserve Gaussianity, so
    f(x) = L^φ_x u(x)  ~  GP(0, k_ff(x, x'; θ, φ))

with the KEY TRANSFORMATION:
    k_ff(x, x'; θ, φ) = L^φ_x L^φ_x' k_uu(x, x'; θ)
    k_uf(x, x'; θ, φ) = L^φ_x' k_uu(x, x'; θ)
    k_fu(x, x'; θ, φ) = L^φ_x  k_uu(x, x'; θ)

Step 3: form the joint Gaussian over (y_u, y_f) data and condition
    K = [k_uu(X_u, X_u) + σ²_n I,  k_uf(X_u, X_f)]
        [k_fu(X_f, X_u)         , k_ff(X_f, X_f) + σ²_n I]

Step 4: maximum likelihood over (θ, φ, σ²) gives the operator
parameters AND the function inference.

## Why this is powerful

- The PDE operator L is embedded in the kernel structure
- Any inverse problem (parameter ID) becomes a standard GP regression
- Solutions come with built-in Bayesian uncertainty
- Works for ODE, PDE, integro-differential, fractional operators
  uniformly

## Limitation: linear operators only

The transformation k_ff = L_x L_x' k_uu requires L to be LINEAR.
For nonlinear PDEs, need:
- Linearization around an operating point
- Or replace GP with another method (PINN)

## Example: Maxwell wave equation in 1D

L_x u = ∂²u/∂t² - c² ∂²u/∂x²

For squared-exponential prior k_uu(x, x') = σ²exp(-||x-x'||²/(2ℓ²)),
the kernel k_ff is COMPUTABLE in closed form (just second
derivatives of Gaussians).  Wave-speed c becomes a hyperparameter
learned from data.

Henderson 2023 gives the explicit kernels for the wave equation.
"""

PI_GP_PDE_SOLVERS = """
# PI-GP generalizes classical linear PDE solvers (Pförtner 2023)

A surprising 2023 result: physics-informed GP regression, in the
DETERMINISTIC LIMIT (zero prior variance, infinite data on the
boundary), recovers CLASSICAL FEM / spectral methods.

## The bridge

Consider Poisson:  -Δu = f  on Ω,  u = g on ∂Ω.

PI-GP setup:
- Prior k_uu(x, x') (kernel choice = "basis function" choice in FEM language)
- Observe f at collocation points {x_c}
- Observe g at boundary points {x_b}
- Train hyperparameters

In the limit where prior variance → 0, the GP mean satisfies:
    -Δ u(x) = f(x) at collocation points
    u(x) = g(x) at boundary points

Which is exactly the COLLOCATION METHOD for FEM.  Different kernels
give different basis functions, which give different FEM-type
methods.

## Why this matters

- Unifies "classical FEM" and "physics-informed ML" under one
  framework
- Provides a Bayesian uncertainty quantification for free in
  classical methods
- Suggests new methods: e.g. learnable kernels that adapt to
  the problem geometry

## Connection to FEEC

The kernel choice in PI-GP is analogous to the choice of finite-
element space in FEEC (Whitney elements, Nedelec, etc.).  Both
have:
- A function space (basis functions / RKHS)
- A way to enforce the PDE (Galerkin / kernel-conditioning)
- Asymptotic optimality theorems

The Whitney-element FEM has the structure-preserving guarantee
(d²=0 exact).  PI-GP doesn't automatically preserve such
structure but can be augmented to.

A natural research direction: "FEEC-Pythagorean GP" — use
Whitney-form-derived kernels in PI-GP to combine structure
preservation with Bayesian uncertainty.

Reference: Pförtner et al, "Physics-Informed Gaussian Process
Regression Generalizes Linear PDE Solvers", arXiv:2207.12382 (2023).
"""

PINN_FORMULATION = """
# Physics-Informed Neural Networks (PINN) — Raissi-Karniadakis 2019

## The construction

Approximate the PDE solution u(x) by a neural network u_NN(x; w).

Define a loss combining:
    L_residual = mean_x_collocation |L u_NN(x) - f(x)|²
    L_boundary = mean_x_boundary  |u_NN(x_b) - g(x_b)|²
    L_data     = mean_x_data      |u_NN(x_data) - y_data|²

Train the NN to minimize  L = L_residual + L_boundary + L_data.

Magic: automatic differentiation of the NN gives L u_NN(x) for
free, so we never need to discretize / mesh the PDE operator.

## Maxwell-specific PINN considerations

For Maxwell's equations:
- The CURL operator must be computed via autograd on the NN output
- Boundary conditions (perfectly conducting walls, ABC, etc.) need
  to be incorporated into L_boundary
- For coupled E-B fields, use multi-output NN
- Sharp material interfaces are hard (NN smoothness conflicts with
  jumps in B or H)

## Practical recipe

1. Choose architecture: 5-8 hidden layers, ~50 neurons each, tanh
   activation (smooth derivatives)
2. Use ADAM optimizer for ~1000 epochs, then L-BFGS for refinement
3. Sample collocation points: uniform in domain, denser near
   boundaries / shocks
4. Watch the loss curve: if L_residual plateaus, increase NN
   capacity; if L_boundary plateaus, increase boundary sampling

## State of the art (2026)

- 2D Maxwell scattering: PINN works, ~5-10 hours training, sub-1%
  L² error.  Useful for parameter sweeps after one-time training.
- 3D full-wave: still challenging.  Hybrid PINN + FEM solver
  is the trend.
- Inverse design (material parameter ID): PINN excels — single
  training run gives parameters + field simultaneously.

## Caveats vs FEM

- No mathematical proof of convergence in general
- Hyperparameters (NN size, learning rate, sampling) matter a lot
- Long training time (NN > FEM for one-off problems)
- Doesn't preserve discrete conservation (d²=0) — spurious modes
  in long-time integration

For LAB practice: use PINN for:
- Inverse problems
- Multi-fidelity (combine measurement + FEM data)
- Real-time query after offline training
Stick with FEM for:
- One-off forward solves
- Problems with sharp interfaces / shocks
- When a certified error estimate matters
"""


def get_pinn_documentation(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
      "all"
      "overview"      - PI-AI vs FEM, PINN vs PI-GP comparison
      "pi_gp"         - Physics-Informed Gaussian Process foundation (Raissi 2017)
      "pi_gp_solvers" - PI-GP generalizes classical PDE solvers (Pförtner 2023)
      "pinn"          - PINN formulation, Maxwell applications, state of the art
    """
    topic = topic.lower().strip()
    if topic in ("overview", "intro"):
        return PINN_OVERVIEW
    if topic in ("pi_gp", "gp", "raissi_2017"):
        return PI_GP_FOUNDATION
    if topic in ("pi_gp_solvers", "pde_solvers", "pfoertner", "pförtner"):
        return PI_GP_PDE_SOLVERS
    if topic in ("pinn", "neural_network", "raissi_2019"):
        return PINN_FORMULATION
    if topic == "all":
        return "\n\n".join([
            PINN_OVERVIEW, PI_GP_FOUNDATION, PI_GP_PDE_SOLVERS, PINN_FORMULATION,
        ])
    return (
        f"Unknown topic '{topic}'. Available: "
        "all, overview, pi_gp, pi_gp_solvers, pinn."
    )
