"""MCMC applications to EM engineering problems."""

DIRECT_MCMC_LAPLACE = r"""
# Direct MCMC solution of Laplace's equation (Shadare-Sadiku 2019)

Reference: Shadare, Sadiku, Musa, "Markov Chain Monte Carlo Solution
of Laplace's Equation in Axisymmetric Homogeneous Domain", Open Journal
of Modelling and Simulation 7(4):203-216, 2019.
DOI: 10.4236/ojmsi.2019.74012
Lab copy: W:/.../04_機械学習と最適化/03_最適化_モンテカルロ/

## The idea

Solve Laplace's equation ∇²V = 0 in an axisymmetric region with
inhomogeneous Dirichlet BC. The probabilistic potential theory
(Kakutani 1944, Doob, Hsu) tells us:

    V(x) = E_x[V(X_T)]

where X_T is a random walk starting at x and stopping when it first
hits the boundary; V(X_T) is the boundary value at the hit point.

→ To compute V at one interior point, run many random walks from
that point; the average boundary value is V(x).

## Algorithm (whole-field MCMC, axisymmetric)

1. **Discretize** the (ρ, z) plane on a uniform grid with spacing h.
2. **Transition probabilities** (for axisymmetric Laplace, ρ > 0):

       p_(±ρ) = 1/4 ± h/(8ρ)   (radial)
       p_(±z) = 1/4              (axial)

   These come from the discretization of `∂²V/∂ρ² + (1/ρ)∂V/∂ρ +
   ∂²V/∂z² = 0` (Equation 2 in the paper).

3. **For each interior node** x:
   - Launch N random walks. Each step picks a direction with the
     above probabilities.
   - When a walk hits a Dirichlet boundary node, record the boundary
     value V_b.
   - V(x) = mean of recorded V_b values over the N walks.

4. **Whole-field**: Repeat for every interior node. (Naive O(N_grid ×
   N_walks × steps_per_walk)). Or use the Markov chain transition
   matrix to compute hitting probabilities in one linear solve.

## Why this works

The discrete Laplace equation `V_i = sum_{j adj} p_{ij} V_j` is the
**same** equation as `V_i = sum_j P_{ij}(hitting boundary) V_b_j`
for the transition matrix P. So a random walk simulates the discrete
Laplacian.

## Comparison with FDM / FEM / BEM

| Method | Computational cost | Memory | Code complexity | Best for |
|--------|-------------------|--------|-----------------|----------|
| FDM | O(N) per point, O(N²) full field via iteration | low | low | Whole field, regular grid |
| FEM (NGSolve) | O(N) sparse solve | low | moderate | Complex geometry |
| BEM | O(N²) dense | high | high | Open boundary |
| **MCMC (this paper)** | O(N_walks × walk_length) per point | tiny | very low | Few isolated points, no iteration |

## When MCMC-Laplace pays off (Shadare's argument)

- **Single-point evaluation**: only need V at a few specific points
  → MCMC computes them directly without solving the whole field.
- **Massively parallel**: each random walk is independent → trivial
  multi-process / GPU parallelization.
- **No matrix to factor / store** → constant memory regardless of
  grid size.
- **No iteration to converge** → naturally embarrassingly parallel.

## When it FAILS to compete with FEM

- **Whole-field computation**: O(N × N_walks) >> FEM O(N) sparse
  solve. Don't use MCMC if you need the whole field at moderate
  resolution.
- **Strong material contrast / nonlinearity**: random walk
  probabilities become uneven, walks waste time.
- **High accuracy**: Monte Carlo error is O(1/sqrt(N_walks)) →
  needs huge N to match FEM. FEM gives O(h^p) deterministic.

## Lab application potential

| Problem | MCMC-Laplace works? |
|---------|---------------------|
| Coil field at single observation point | YES — better than FEM if grid is huge |
| Capacitor field around a corner | YES — adaptive walks concentrate near corner |
| Whole-domain B distribution for motor | NO — FEM dominates |
| Open-boundary problem | YES — walks naturally extend; no PML / Kelvin needed |
| **Inverse problem (boundary V from interior measurements)** | YES — recast as Bayesian inverse via MCMC |

The lab's existing Kelvin transformation in NGSolve handles open
boundary cleanly, so the unique value of MCMC-Laplace in the lab is
single-point evaluation in irregular domains.

## Python prototype

```python
import numpy as np

def mcmc_laplace_axisym(rho_z, V_boundary_fn, mask_boundary,
                         n_walks=10000, max_steps=10000, h=1e-3):
    # Estimate V at one point via axisymmetric random walks.
    # rho_z: (rho, z) tuple, the query point
    # V_boundary_fn: callable (rho, z) -> V on boundary
    # mask_boundary: callable (rho, z) -> bool (True if boundary)
    rng = np.random.default_rng()
    V_hits = []
    for _ in range(n_walks):
        rho, z = rho_z
        for _ in range(max_steps):
            if mask_boundary(rho, z):
                V_hits.append(V_boundary_fn(rho, z))
                break
            # Axisymmetric transition probabilities (Equation 2)
            p_plus_rho  = 0.25 + h / (8 * max(rho, h))
            p_minus_rho = 0.25 - h / (8 * max(rho, h))
            p_plus_z    = 0.25
            p_minus_z   = 0.25
            u = rng.uniform()
            if u < p_plus_rho:                rho += h
            elif u < p_plus_rho + p_minus_rho: rho -= h
            elif u < p_plus_rho + p_minus_rho + p_plus_z: z += h
            else:                              z -= h
            if rho < 0: rho = -rho   # reflect axis
    return np.mean(V_hits), np.std(V_hits) / np.sqrt(len(V_hits))
```

## Production lab use

For a 3D EM forward solve, FEM (NGSolve) is faster + more accurate.
But for:
- Bayesian inverse problems where you embed the forward solve in a
  posterior, MCMC-Laplace can be a cheap forward — fewer FEM calls.
- Pedagogical demonstration that EM and Bayesian inference share
  the same Markov-chain / detailed-balance mathematics.

Cross-link: `radia_mcp.mcmc('algorithms')`, `radia_mcp.radia_ngsolve`
(for the FEM alternative).
"""


BAYESIAN_INVERSE = r"""
# Bayesian Inverse Problems for EM

Given measured data `y` (B-field at sensor, terminal voltage,
inductance), infer the **posterior distribution** over the design
parameters `theta` (material properties, geometry):

    p(theta | y) propto p(y | theta) * p(theta)

The forward model `y = f(theta) + noise` provides the likelihood
`p(y | theta) = N(f(theta), sigma_noise²)`.

## Why Bayesian (not just least-squares fit)

Plain least-squares gives ONE estimate. Bayesian gives the full
posterior:
- **Uncertainty quantification**: 95% credible interval, not just
  mean
- **Multi-modality**: if the inverse is non-unique (e.g. two BH
  curves both fit), MCMC finds both
- **Identifiability check**: posterior flat in some direction = that
  param is not identifiable from the data
- **Prior incorporation**: respect physical bounds + manufacturer
  datasheets

## EM applications

### Inverse 1: BH curve identification

Given measured B(t), H(t) loops at multiple amplitudes, infer:
- Saturation flux density B_sat
- Initial permeability mu_init
- Anhysteretic curve shape parameters

```python
import pymc as pm
import numpy as np

with pm.Model() as model:
    # Priors on BH curve parameters (e.g. Frohlich-Kennelly)
    B_sat = pm.Normal("B_sat", mu=1.8, sigma=0.2)
    H_c   = pm.HalfNormal("H_c", sigma=50)         # coercivity
    mu_init = pm.HalfNormal("mu_init", sigma=5000)

    # Forward model: simulate B(H) curve
    H_grid = np.linspace(-1000, 1000, 100)
    B_pred = bh_frohlich(H_grid, B_sat, H_c, mu_init)

    # Likelihood: observed B at the same H_grid
    B_obs = pm.Normal("B_obs", mu=B_pred, sigma=0.05,
                       observed=B_meas)
    trace = pm.sample(2000, tune=1000, chains=4)
```

### Inverse 2: Play model parameter ID (lab core)

Cross-link: `radia_mcp.magnetic_materials` (Play / Energy models).

Play hysteresis model has K play operators with thresholds eta_k and
shape functions f_k(r). MCMC infers (K, eta, f_k) from measured B(H)
loop.

```python
def log_posterior(theta):
    K = int(theta[0])
    eta = np.sort(np.exp(theta[1:K+1]))   # log-uniform increments
    f_k_params = theta[K+1:]
    mat = build_play_model(K, eta, f_k_params)
    B_sim = simulate_loops(mat, H_meas)
    return -0.5 * np.sum((B_sim - B_meas)**2 / sigma_meas**2)

# Use emcee (K varies → can't NUTS) or reversible-jump MCMC
```

### Inverse 3: Permeability map from coil measurement

Measure inductance at multiple coil positions. Infer the spatially-
varying permeability mu_r(x, y, z) in an iron core.

```python
# Discretize mu_r on a coarse voxel grid
n_vox = 10**3
log_mu_r = np.array([Normal(log(100), 1) for _ in range(n_vox)])

# Forward: NGSolve FEM with mu_r(x) = exp(log_mu_r[voxel(x)])
def log_posterior(log_mu_r):
    mu_r_field = build_mu_r_grid(log_mu_r)
    L_pred = run_ngsolve_inductance(mu_r_field, coil_positions)
    return -0.5 * np.sum((L_pred - L_meas)**2 / sigma_L**2)

# 1000-dim problem → use NUTS with PINN surrogate (radia_mcp.pinn)
```

### Inverse 4: Magnetic source localization (biomagnetism)

Given B field measurements outside the body, infer the location and
strength of internal magnetic sources (current dipoles in MEG).

The Saotome 1995 Sampled Pattern Matching method (covered in
`radia_mcp.mcmc('spm')`) tackles the same problem deterministically;
MCMC adds posterior quantification.

## Tricky aspects

1. **Forward model cost** — each MCMC step is one full FEM solve.
   N_steps × time_per_FEM may be hours-to-days. Mitigations:
   - Surrogate model (PINN, GP, Bayesian NN) → `radia_mcp.pinn`
   - Adaptive MH with delayed acceptance (cheap surrogate proposes,
     full FEM only accepts/rejects)
   - Coarse FEM in early MCMC iterations, fine FEM near posterior mode

2. **Identifiability** — many EM inverse problems have non-unique
   solutions (gauge invariance, surface vs volume currents,
   sub-pixel sources). The posterior may have ridges or modes that
   indicate this. Add regularization via prior.

3. **Noise model** — Gaussian likelihood is a strong assumption.
   If measurement noise is correlated (sensor drift) or non-Gaussian
   (clipping, quantization), build a proper noise model into the
   likelihood.

## Cross-links

- `radia_mcp.magnetic_materials` — Play / Energy / BH curve models
- `radia_mcp.pinn` — surrogate models for cheap forward
- `radia_mcp.optuna` — point estimate / MAP if posterior not needed
- `radia_mcp.mcmc('spm')` — Saotome SPM (deterministic inverse)
- `radia_mcp.team_benchmark` — TEAM Problem 32 (anisotropic BH ID),
  TEAM Problem 33b (inverse identification benchmark)
"""


BAYESIAN_UQ = r"""
# Bayesian Uncertainty Quantification (UQ) for EM Design

Forward UQ: given uncertain inputs, predict the distribution of
outputs.

    p(y) = ∫ p(y | theta) p(theta) dtheta

In EM design:
- Material properties have manufacturing tolerances (mu_r ± 10%,
  sigma ± 5%)
- Geometry has fabrication tolerances (gap thickness ± 0.05 mm)
- Operating conditions vary (temperature, frequency)

UQ asks: what is the probability that the device meets spec
(e.g. P(L > 95 µH AND R < 50 mΩ))?

## Naive Monte Carlo

1. Sample theta_i ~ p(theta), i = 1..N
2. Run forward y_i = f(theta_i) for each
3. Empirical distribution of {y_i} is p(y)

Cost: N forward calls. For FEM at minutes each + N = 1000, that's
a day → use surrogate.

## Polynomial chaos expansion (PCE) as a cheaper alternative

Approximate y = f(theta) as a polynomial in theta. Compute moments
of y analytically from the polynomial coefficients.

```python
from chaospy import Normal, Uniform, generate_expansion, generate_quadrature

# Define uncertain inputs
mu_r_dist = Normal(100, 10)
sigma_dist = Uniform(1e6, 5e6)
joint = chaospy.J(mu_r_dist, sigma_dist)

# Sample at quadrature nodes
nodes, weights = generate_quadrature(4, joint, rule="gaussian")

# Run forward at each node (only ~20-50 evaluations!)
y_samples = [run_fem(*n) for n in nodes.T]

# Fit polynomial expansion
expansion = generate_expansion(4, joint)
approx = chaospy.fit_quadrature(expansion, nodes, weights, y_samples)

# Get moments
print("Mean:", chaospy.E(approx, joint))
print("Std:", chaospy.Std(approx, joint))
print("P(y > 100):", 1 - chaospy.Prob(approx, [100], joint))
```

PCE is **much** cheaper than MCMC for forward UQ when:
- Input dimension < ~10
- Output is smooth in inputs

For higher dimension or non-smooth output, fall back to MCMC /
quasi-MC.

## MCMC for INVERSE UQ

The Bayesian inverse problem IS uncertainty quantification — you
get a posterior, not a point estimate. Cross-link to
`bayesian_inverse` topic above.

## Sobol sensitivity analysis

Decompose total output variance into contributions from each input
and their interactions. Identify which input dominates the
uncertainty:

```python
from SALib.sample import sobol
from SALib.analyze import sobol as analyze_sobol

problem = {
    "num_vars": 4,
    "names": ["mu_r", "sigma", "gap", "freq"],
    "bounds": [[80, 120], [1e6, 5e6], [0.45e-3, 0.55e-3], [40e3, 60e3]],
}
param_values = sobol.sample(problem, 1024)
Y = np.array([run_fem(*p) for p in param_values])

Si = analyze_sobol.analyze(problem, Y)
print("First-order Sobol:", Si["S1"])
print("Total-order Sobol:", Si["ST"])
```

Use this BEFORE expensive UQ — fix the inputs that don't matter
(low Sobol index), reducing the dim of subsequent MCMC.

## EM-specific UQ workflow

```
1. Define uncertain inputs + bounds + distributions (manufacturing
   datasheets, calibration uncertainty)

2. Sobol sensitivity to rank inputs by importance (cheap, ~1k FEM)

3. If <5 important inputs: PCE for forward UQ (cheap, ~50 FEM)
   If >10 inputs: surrogate (PINN) + MCMC posterior

4. Compute P(spec satisfied), 95% credible intervals on outputs

5. Robust design: optimize the input MEANS to maximize
   P(spec satisfied), not just nominal output
```

## Cross-links

- `radia_mcp.pinn` — Gaussian process / Bayesian NN surrogate
- `radia_mcp.optuna` — robust design as outer BBO
- `radia_mcp.mcmc('bayesian_inverse')` — inverse UQ
"""


DESIGN_UNDER_UNCERTAINTY = r"""
# Design Under Uncertainty (Robust Optimization)

Optimize a design that performs well under input uncertainty, not
just at nominal conditions.

## Classical robust formulation

Min over design `d`:

    E_theta [ -L(d, theta) ]              (expected loss)
    + lambda * Var_theta [ L(d, theta) ]  (variance penalty)

Or chance-constrained:

    min E[L(d, theta)]
    s.t. P(g(d, theta) <= 0) >= 1 - alpha   (alpha = 0.05 means 95% feasibility)

where theta is the uncertain input vector (with known distribution).

## Why MCMC enters

Computing E and Var requires sampling theta. Naive Monte Carlo: at
each candidate `d`, draw N samples of theta, run N forward solves.
With FEM forward and a Bayesian outer loop (e.g. PyMC over `d`),
total cost = N_outer * N_inner * FEM.

Mitigations:
- **Common random numbers**: re-use the same theta samples across
  candidate `d` values
- **Surrogate**: PINN / Bayesian NN for FEM
- **Two-loop SMC** (sequential Monte Carlo): propagate uncertain
  theta and design d jointly

## EM example: IH coil tolerant design

Inputs:
- Workpiece radius R_wp uncertain: N(20mm, 0.5mm) (sample-to-sample)
- Workpiece resistivity: N(20 uOhm·cm, 5%)
- Coil-to-workpiece gap: U(2mm, 4mm) (varies in production)

Design: coil turns N, coil radius R_c, drive frequency f.

```python
import optuna
import numpy as np

def robust_objective(trial):
    # Design variables
    n_turns = trial.suggest_int("n_turns", 5, 15)
    R_c    = trial.suggest_float("R_c", 25e-3, 60e-3)
    f      = trial.suggest_categorical("f", [10e3, 25e3, 50e3])

    # Inner MC over uncertain inputs
    n_inner = 50
    P_wp_samples = []
    for _ in range(n_inner):
        R_wp = np.random.normal(20e-3, 0.5e-3)
        rho_wp = np.random.normal(20e-8, 1e-9)
        gap = np.random.uniform(2e-3, 4e-3)
        result = run_ih_panel(n_turns, R_c, R_wp, rho_wp, gap, f)
        P_wp_samples.append(result["P_wp_W"])

    P_wp_mean = np.mean(P_wp_samples)
    P_wp_std  = np.std(P_wp_samples)

    # Robust score: penalize variance
    return -(P_wp_mean - 2.0 * P_wp_std)   # maximize lower-CI bound

# Outer Optuna picks the design
study = optuna.create_study(
    direction="minimize",
    sampler=optuna.samplers.TPESampler(multivariate=True))
study.optimize(robust_objective, n_trials=200)
```

50 inner samples × 200 outer trials = 10000 FEM solves → use a
surrogate (`radia_mcp.pinn`) to scale to production sweep.

## When MCMC vs Optuna for robust design

- Use **Optuna + inner MC**: when you only need point estimates of
  robust score (mean - 2*sigma)
- Use **MCMC + variational**: when you need the full posterior of
  design given measurements (Bayesian robust design)

## Cross-links

- `radia_mcp.optuna('coil_design')` — IH design without uncertainty
- `radia_mcp.pinn` — surrogate for inner MC speedup
- `radia_mcp.mcmc('bayesian_uq')` — forward UQ (deterministic design)
"""


def get_em_applications_documentation(topic: str = "all") -> str:
    """Dispatch EM application topics.

    Topics:
        direct_mcmc_laplace      - Shadare-Sadiku 2019 MCMC for Laplace
        bayesian_inverse         - Material params / BH curve / source ID
        bayesian_uq              - Forward UQ, PCE, Sobol sensitivity
        design_under_uncertainty - Robust optimization with MCMC inner
        all                      - Everything
    """
    topic = topic.lower().strip()
    if topic in ("direct_mcmc_laplace", "laplace", "direct_mcmc",
                 "shadare", "field_solver"):
        return DIRECT_MCMC_LAPLACE
    if topic in ("bayesian_inverse", "inverse", "inverse_problem",
                 "material_id", "bh_id"):
        return BAYESIAN_INVERSE
    if topic in ("bayesian_uq", "uq", "uncertainty", "pce", "sobol"):
        return BAYESIAN_UQ
    if topic in ("design_under_uncertainty", "robust", "robust_design",
                 "uncertainty_design"):
        return DESIGN_UNDER_UNCERTAINTY
    if topic == "all":
        return "\n\n".join([DIRECT_MCMC_LAPLACE, BAYESIAN_INVERSE,
                            BAYESIAN_UQ, DESIGN_UNDER_UNCERTAINTY])
    return (f"Unknown topic '{topic}'. Available: direct_mcmc_laplace, "
            "bayesian_inverse, bayesian_uq, design_under_uncertainty, all.")
