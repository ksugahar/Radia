"""Bayesian optimization, GP regression, FMQA, surrogate model knowledge."""


# Authoritative topic enum for the dispatcher tool (wired into
# `<short>_topics()` via common.register_topics_tool).
TOPICS: dict[str, str] = {
    "gp_fundamentals": "(see knowledge file for details)",
    "bayesian_optimization": "(see knowledge file for details)",
    "fmqa": "(see knowledge file for details)",
    "physics_informed_gp": "(see knowledge file for details)",
    "surrogate_models": "(see knowledge file for details)",
    "gp_libraries": "(see knowledge file for details)",
    "llm_bo": "(see knowledge file for details)",
    "all": "return \"\n\n\".join([GP_FUNDAMENTALS, BAYESIAN_OPTIMIZATION, FMQA,",
}

GP_FUNDAMENTALS = r"""
# Gaussian Process Regression Fundamentals

A Gaussian process (GP) defines a distribution over functions. For
input `x`, the function value `f(x)` is a random variable; for any
finite set of inputs, the joint distribution is multivariate Gaussian.

## Specification

A GP is fully specified by:
- **Mean function** m(x), often 0 after centering
- **Kernel** (covariance function) k(x, x') — controls smoothness +
  characteristic length scale

The posterior at a test point `x*` given observations
`(X_train, y_train)` is Gaussian with closed-form mean and variance:

    K = k(X, X) + sigma_n^2 * I        (noise on training)
    mu(x*)  = k(x*, X) K^-1 y
    sig²(x*) = k(x*, x*) - k(x*, X) K^-1 k(X, x*)

The cubic cost `O(N^3)` for matrix inverse limits standard GP to
N ~ 1000-10000 (sparse GP / SVGP scale further).

## Common kernels

| Kernel | Form | Use case |
|--------|------|----------|
| **RBF** (squared exp) | exp(-‖x-x'‖²/(2l²)) | Default smooth functions |
| **Matern 3/2** | (1 + √3r/l) exp(-√3r/l), r=‖x-x'‖ | Rougher functions; one diff |
| **Matern 5/2** | similar with √5r/l | Standard for BO |
| **Periodic** | exp(-2 sin²(πr/p)/l²) | Periodic process (rotation angle) |
| **ARD-RBF** | exp(-Σ_i (x_i-x'_i)²/(2l_i²)) | Per-dimension length scale (auto-relevance) |
| **Linear** | x·x' + sigma² | Linear trend |
| **White noise** | sigma_n² delta | Observation noise |

Kernels are **composable** by addition and multiplication:
- k_RBF + k_Linear → smooth + linear trend
- k_RBF × k_Periodic → seasonal modulation
- This is "kernel engineering"

## ARD (Automatic Relevance Determination)

ARD-RBF has one length scale per input dimension:

    k_ARD(x, x') = sigma_f² * exp(-1/2 Σ_i (x_i - x'_i)² / l_i²)

After fitting, **large l_i means input i is irrelevant** (changes in
x_i don't affect output). Inverse `1/l_i²` gives a per-input
"relevance score".

Use ARD when you have many inputs and want feature importance for
free as a byproduct of GP fit.

Lab reference: `ARD カーネルによる非線形地震応答解析のガウス過程回帰
代替モデル構築.pdf` (public-safe curated corpus) uses ARD-RBF to
build a GP surrogate for nonlinear seismic response analysis — same
principle applies to nonlinear EM (saturating iron, BH curves).

## Hyperparameter optimization

Marginal log-likelihood:

    log p(y | X, theta) = -1/2 y^T K^-1 y - 1/2 log|K| - N/2 log(2 pi)

Maximize w.r.t. kernel params theta (length scales, sigma_f, sigma_n).
Gradient available analytically. Often multi-modal; restart from
several random initializations.

```python
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel

kernel = ConstantKernel(1.0) * RBF(length_scale=[1.0]*N_inputs) \
         + WhiteKernel(noise_level=1e-3)
gpr = GaussianProcessRegressor(kernel=kernel,
                                 n_restarts_optimizer=10,
                                 normalize_y=True)
gpr.fit(X_train, y_train)
mu, sigma = gpr.predict(X_test, return_std=True)
```

## When GP regression shines (for EM)

- Few hundred samples (N < 1000): exact GP is fast
- Smooth output (e.g. inductance vs geometry parameter)
- Need predictive uncertainty (95% CI) at test points
- Want to identify which inputs matter (ARD)
- Drop-in surrogate for FEM in optimization loop

## When GP fails

- N > 10000: cubic cost prohibitive → use sparse GP (inducing points)
  or BoTorch's variational SVGP
- High dim (D > 50): kernel breaks down ("curse of dimensionality")
  → use additive structures, deep kernels, or NN-based surrogate
- Discontinuous / non-stationary output: kernel assumes smoothness
  → use mixture of kernels, change point models, or BNN

Cross-link: `radia_mcp.pinn` for the NN-based alternative when GP
fails on dim or N.
"""


BAYESIAN_OPTIMIZATION = r"""
# Bayesian Optimization (BO)

Use GP as a surrogate for an expensive black-box function. Iteratively:
1. Fit GP to all observed (x, y) pairs
2. Define an **acquisition function** alpha(x) that scores candidates
3. Pick x_next = argmax alpha(x); evaluate f(x_next)
4. Add to data; repeat

The acquisition function balances **exploration** (sample where GP
uncertainty is high) and **exploitation** (sample where GP mean is
best).

## Acquisition functions

| Name | Formula | Behavior |
|------|---------|----------|
| **PI** (Probability of Improvement) | P(f(x) > best + xi) | Greedy; favors exploitation |
| **EI** (Expected Improvement) | E[max(f(x) - best, 0)] | Default; balanced |
| **UCB** (Upper Confidence Bound) | mu(x) + kappa * sigma(x) | Tunable via kappa |
| **KG** (Knowledge Gradient) | E[improvement of model] | Lookahead; expensive |
| **qNEI** (q-Noisy EI) | EI for batch of q points | Parallel BO |
| **qNEHVI** | qNEI for multi-objective | Pareto front BO |
| **MES** (Max-value Entropy Search) | Info gain on max | Anti-greedy |

EI is the workhorse — almost always good enough. UCB with `kappa = 2.0`
is a transparent alternative.

## Connection with Optuna

`optuna.integration.BoTorchSampler` is Optuna's BO sampler — uses
BoTorch (Meta) for GP + acquisition optimization. Use it when:
- < 100 trials affordable
- Continuous params (BoTorch struggles with categorical)
- Smooth posterior

For more trials or mixed types, TPE (Optuna default) usually wins.

Cross-link: official `optuna/optuna-mcp` or plain Optuna APIs for the
Optuna-side sampler/pruner view.

## BO libraries

| Library | Backend | Strengths |
|---------|---------|-----------|
| **BoTorch** | PyTorch | State-of-the-art acquisition, parallel q-batch, multi-fidelity, multi-objective (qNEHVI) |
| **GPyTorch** | PyTorch | GP-only; flexible kernels; large-scale via Whitening |
| **scikit-optimize** | scikit-learn | Old, simple, scipy-native |
| **GPyOpt** | GPy | Old (2019), well-documented |
| **Ax** | BoTorch + Adapt | Meta's experimentation framework on top of BoTorch |

For lab use:
- Optuna + BoTorchSampler — usual entry point
- Direct BoTorch for: parallel q-batch BO, multi-fidelity, MO

## Multi-fidelity BO

When you have a cheap-but-noisy model (coarse FEM) AND an expensive-
but-accurate one (fine FEM), multi-fidelity BO chooses between them.

```python
import botorch
from botorch.models import SingleTaskMultiFidelityGP
# Treat 'fidelity' as an extra dim; GP learns fidelity correlation
```

Lab use case: coil design where 2D FEM is 10s, 3D FEM is hours →
ramp up to 3D only near promising candidates.

## Multi-objective BO (qNEHVI)

For Pareto front (e.g. torque vs loss, P_wp vs P_coil):

```python
from botorch.acquisition.multi_objective import qNoisyExpectedHypervolumeImprovement
```

Beats NSGA-II for ~50 trials with smooth output; loses to NSGA-II
for > 500 trials.

## Robust / risk-averse BO

Optimize the lower CI bound `mu - 2*sigma` instead of the mean —
penalizes high-uncertainty regions. Useful when each trial is very
expensive and you can't afford to "explore wrong" near the end.

## Cost-aware BO

`BoTorch.acquisition.CostAwareAcquisition` divides EI by evaluation
cost (FEM time as function of mesh density). Lets BO trade off
accuracy vs speed.

## Stopping criteria

- Fixed budget: stop at N_trials (most common)
- Convergence: stop when EI < 1% of best so far
- Confidence interval: stop when 95% CI of best < tolerance
"""


FMQA = r"""
# FMQA: Factorization Machine + Quantum Annealing

A hybrid BBO algorithm for **binary / categorical** design variables
where the evaluation cost is high.

## Key idea

Build a surrogate model that is **quadratic-form-friendly** — so
the next-candidate optimization can be solved by a quantum annealer
(D-Wave) or simulated annealing.

The surrogate is a **Factorization Machine** (FM, Rendle 2010):

    f_FM(x) = w0 + Σ_i w_i x_i + Σ_{i<j} <v_i, v_j> x_i x_j

where `v_i` ∈ R^k are factor vectors of rank k. This captures
pairwise interactions with O(k * N) params (vs O(N²) for full
quadratic).

FM training: stochastic gradient descent on observed (x, y) pairs.
After training, `f_FM(x)` is a quadratic in binary x → fits the
QUBO (Quadratic Unconstrained Binary Optimization) form natively.

## Algorithm

```
1. Sample N initial random binary configs x_i; evaluate y_i = f(x_i)
2. While budget remains:
   a. Train FM on (X, y) data
   b. Extract Q matrix from FM (quadratic form coefficients)
   c. Solve QUBO: x* = argmin x^T Q x using
       - Quantum annealing (D-Wave)
       - Simulated annealing (CPU, free)
       - openjij / Fixstars Amplify (commercial SA)
   d. Evaluate y* = f(x*); add (x*, y*) to data
3. Return best (x, y) found
```

## When FMQA beats Optuna / GP-BO

- Variables are **binary / 0-1 / one-hot**
- Strong pairwise interactions matter (FM captures them)
- N > 100 dimensions where GP breaks down
- You have access to D-Wave OR simulated annealing tooling

## When FMQA fails

- Continuous variables (FM needs discretization; loses resolution)
- Few interactions (linear regression is enough)
- Very few trials (FM needs ~50 data points to fit)

## EM applications

### Metamaterial design (★ lab paper)

`量子アニーリングを用いたメタマテリアル設計.pdf` (public-safe curated corpus
04_最適化_ベイズ_BBO/). Variables: which unit cells are "filled"
or "empty" in a 2D/3D metamaterial lattice (binary). Objective:
target effective permittivity / permeability response.

```python
import numpy as np
from openjij import SQASampler

def fmqa_metamaterial(n_cells=64, target_eps_eff=2.0, budget=200):
    # n_cells binary: 1 = filled, 0 = empty
    # Initial random
    X = np.random.randint(0, 2, (20, n_cells))
    Y = np.array([cae_effective_eps(x) for x in X])

    for trial in range(20, budget):
        # 1. Fit FM
        fm = train_factorization_machine(X, Y - target_eps_eff)
        # 2. Extract Q
        Q = fm.quadratic_matrix(N=n_cells)
        # 3. Solve QUBO
        sampler = SQASampler()
        result = sampler.sample_qubo(Q, num_reads=100)
        x_new = list(result.first.sample.values())
        # 4. Evaluate
        y_new = cae_effective_eps(x_new)
        X = np.vstack([X, x_new])
        Y = np.append(Y, y_new)

    best_idx = np.argmin(np.abs(Y - target_eps_eff))
    return X[best_idx], Y[best_idx]
```

### Topology optimization (binary density)

For binary SIMP (no intermediate density), FMQA is competitive with
discrete level-set methods.

### Multi-material assembly

Picking which of M materials goes into each of N cells → discrete
NM choice; FM captures interactions between adjacent cell choices.

## libFM (the original library)

Rendle's `libFM` (Factorization Machines with libFM, JMLR 2012)
provides command-line FM training. Modern Python alternatives:
- `pyFM` (sklearn-style)
- `xLearn` (FM/FFM with GPU)
- TensorFlow / PyTorch custom FM

## QUBO solvers

| Solver | Type | Cost | Notes |
|--------|------|------|-------|
| **D-Wave Advantage** | Quantum annealer | Cloud subscription | 5000+ qubits, but limited connectivity |
| **D-Wave Hybrid** | Quantum + classical | Cloud | Larger problems |
| **openjij** | Simulated annealing (CPU) | Free | Pure Python; reasonable for N < 1000 |
| **Fixstars Amplify** | SA + GPU | Free tier + paid | Japan-based, well-documented |
| **dwave-neal** | Simulated annealing | Free | D-Wave's CPU fallback |

## Cross-references

- `radia_mcp.metamaterial` — physical metamaterial knowledge
- `radia_mcp.topology_optimization` — gradient alternative for
  binary topology
- official `optuna/optuna-mcp` — practical BO/TPE alternative when
  variables are continuous
"""


PHYSICS_INFORMED_GP = r"""
# Physics-Informed Gaussian Processes (PI-GP)

A GP surrogate that **respects a known PDE** — incorporates physics
into the kernel itself, so the model is exact for that PDE class
(no training error from physics).

## Foundational paper

Raissi-Perdikaris-Karniadakis 2017, "Machine learning of linear
differential equations using Gaussian processes", J. Comp. Physics
348:683-693. DOI: 10.1016/j.jcp.2017.07.050.

Lab copy: `public-safe curated corpus et al
- 2017 - Machine learning of linear differential equations using
Gaussian processes.pdf`

## Key construction

For a linear PDE `L u(x) = f(x)`, if `u ~ GP(0, k_u)`, then `f =
L u ~ GP(0, k_f)` with:

    k_f(x, x') = L_x L_x' k_u(x, x')

(Apply the operator on both sides of the kernel.) Now you can jointly
condition on observations of u AND of f — gives a single GP that
exactly satisfies the PDE on average.

## Why this works for EM

Linear EM PDEs (Laplace, static Maxwell, modal Helmholtz with fixed
k) are exactly the PDEs PI-GP handles. Common cases:

| PDE | Linear? | Examples |
|-----|---------|----------|
| Laplace ∇²u = 0 | yes | Electrostatic potential in dielectric |
| Poisson ∇²u = -ρ/ε | yes | Same with sources |
| Helmholtz (k fixed) | yes | Wave scattering at one frequency |
| Diffusion (linear σ, μ) | yes | Eddy current in air |
| **Nonlinear BH** | **no** | Use PINN instead |
| **Coupled circuit-FEM** | **no** | Use FEM directly |

## Pförtner 2023 advance

`Pförtner et al - 2023 - Physics-Informed Gaussian Process Regression
Generalizes Linear PDE Solvers.pdf` shows that PI-GP **generalizes
linear PDE solvers** — i.e. any classical method (FDM, FEM, BEM,
spectral) can be viewed as a special case of PI-GP with a particular
prior + observation operator.

The practical consequence:
- PI-GP gives the SAME solution as your FEM solve, PLUS uncertainty
  estimates everywhere.
- For inverse problems, PI-GP gives the posterior over the unknown
  (source distribution, boundary conditions) for free.

## Henderson 2023 (wave equation)

`Henderson et al - 2023 - Covariance models and Gaussian process
regression for the wave equation.pdf` builds wave-equation-respecting
covariance models. Useful for:
- Acoustic / EM scattering
- Time-of-flight imaging
- Antenna pattern recovery from sparse measurements

## When PI-GP beats FEM

- **Sparse measurement data** + PDE constraint: PI-GP fuses both
  cleanly; FEM doesn't directly use measurements
- **Inverse problems** where the source is unknown: marginalize over
  it in the GP posterior
- **Uncertainty quantification** on the field: free with PI-GP

## When PI-GP loses to FEM

- High-resolution forward simulation: FEM is O(N) sparse; PI-GP is
  O(N³)
- Nonlinear PDE: PI-GP requires linear L; PINN handles nonlinear
- Complex geometry with discontinuous materials: FEM's mesh refines;
  PI-GP needs custom kernels per region

## EM-flavored example: source ID from sparse field measurements

Given measurements of B field at N points around a coil, infer the
coil current distribution.

```python
import numpy as np
import gpytorch
import torch

class WaveEqKernel(gpytorch.kernels.Kernel):
    \"\"\"GP kernel satisfying Laplace equation in 2D.\"\"\"
    has_lengthscale = True
    def forward(self, x1, x2, **params):
        d = self.covar_dist(x1, x2)
        # Choose Matern 5/2 base; apply Laplacian operator
        # ... (PI-GP kernel construction)
        return ...

# Joint observe (B_field, source_proxy = 0 inside non-source region)
gp = MultitaskGPModel(wave_kernel, ...)
# Train via MLL; predict source distribution
```

## Cross-links

- `radia_mcp.pinn` — neural counterpart (nonlinear capable)
- `radia_mcp.fem` — classical FEM alternative
- `radia_mcp.bem` — BEM alternative for unbounded domains
- `radia_mcp.differential_forms` — PDE structure on manifolds
"""


SURROGATE_MODELS = r"""
# Surrogate Models for CAE (Lab Workshop Lineage)

`2024_CAE懇話会_CAEのサロゲートモデル入門/` (public-safe curated corpus
04_最適化_ベイズ_BBO/) — 5-session online workshop on CAE surrogate
modeling, June-October 2024.

## When you need a surrogate

A surrogate model replaces the FEM solve in:
- Inner loop of an optimization (BBO, variational inference)
- Real-time control / digital twin
- Design sweep (parametric study)
- Uncertainty quantification (Monte Carlo)

Rule of thumb: if you'll evaluate the FEM > 100 times in similar
input ranges, fit a surrogate once + use it many times.

## Surrogate model menu

| Model | Best for | Cost to train | Cost to predict |
|-------|----------|---------------|-----------------|
| **Linear regression** | Truly linear, baseline | seconds | microseconds |
| **Polynomial regression** | Low-degree smooth | seconds | microseconds |
| **Radial Basis Function (RBF)** | Smooth, ~100 samples | seconds | microseconds (per sample) |
| **Gaussian Process** | Smooth, with UQ | minutes (N³) | milliseconds (per pred) |
| **Random Forest** | Mixed types, robust | minutes | milliseconds |
| **Gradient Boosting (LightGBM/XGBoost)** | Tabular, fast | minutes | milliseconds |
| **MLP (Neural Net)** | Large N, complex | hours | milliseconds |
| **PINN** | Physics + small data | hours | milliseconds |
| **DeepONet** | Operator learning | hours | milliseconds |
| **Fourier Neural Operator (FNO)** | PDE solution surrogate | hours | milliseconds |

## Decision flowchart

```
Have a PDE you trust?
  YES → PI-GP (linear) or PINN (nonlinear)
  NO →
    < 100 training samples?
      YES → GP / RBF
      NO →
        < 10000 samples?
          YES → Random Forest / GBM / MLP
          NO →
            need PDE-level generalization?
              YES → DeepONet / FNO / GNN
              NO  → MLP / FNN
```

## Training data strategy

- **Design of experiments (DoE)**: LHS, Sobol, factorial designs
  give better coverage than random sampling.
- **Adaptive sampling**: use surrogate uncertainty to pick next
  training point (active learning) — much more data-efficient than
  passive sampling.
- **Multi-fidelity**: combine cheap + expensive models; surrogate
  learns the correlation.

## Lab workshop sessions (CAE懇話会 2024)

| Session | Date | Topic (inferred from folder names) |
|---------|------|------------------------------------|
| 1 | 2024-06-21 | Introduction to CAE surrogate models |
| 2 | 2024-07-26 | (GP / kernel methods) |
| 3 | 2024-08-23 | (Neural network surrogates) |
| 4 | 2024-09-27 | (Multi-fidelity / active learning) |
| 5 | 2024-10-25 | (Application + Q&A) |

Direct content from these sessions is in the lab folder; recommend
attending live for context.

## EM-flavored example: GP surrogate for IH coil inductance

```python
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel

# 1. Design of experiments: LHS over (n_turns, R_coil, gap, freq_kHz)
from scipy.stats.qmc import LatinHypercube
sampler = LatinHypercube(d=4, seed=42)
X_lhs = sampler.random(n=80)
# scale to physical ranges
X = scale_to_ranges(X_lhs, [(3, 15), (20e-3, 80e-3), (1e-3, 10e-3), (10, 200)])

# 2. Train surrogate on 80 FEM evals
y = np.array([run_inductance_fem(*xi)["L_uH"] for xi in X])
kernel = ConstantKernel(1.0) * RBF(length_scale=[1.0]*4) + WhiteKernel(1e-3)
gpr = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10,
                                 normalize_y=True)
gpr.fit(X, y)

# 3. Use surrogate inside Optuna (instead of FEM)
def fast_objective(trial):
    x = np.array([trial.suggest_float(...), ...])
    mu, sigma = gpr.predict(x[None,:], return_std=True)
    # Bonus: penalize high-uncertainty regions
    return -(mu - 2*sigma)
```

## Adaptive sampling (active learning)

```python
# After fitting GP on N0 points, pick next point by max sigma
x_grid = lhs_sample(1000)
mu, sigma = gpr.predict(x_grid, return_std=True)
x_next = x_grid[np.argmax(sigma)]
y_next = run_fem(x_next)
# Refit GP on N0+1 points; repeat
```

This is BO with the **UCB** acquisition reduced to pure exploration
(`kappa = inf`). For real BO, balance with exploitation.

## Cross-links

- official `optuna/optuna-mcp` — study/trial orchestration around
  surrogate objectives
- `radia_mcp.pinn` — neural surrogates for PDEs
- `radia_mcp.mor` — model order reduction is the "physics-aware
  surrogate" alternative (PRIMA, Cauer ladder)
- `radia_mcp.mathematica` — symbolic baseline for verification
"""


GP_LIBRARIES = r"""
# GP / BO library landscape

## Python GP libraries

| Library | Backend | Strengths | When |
|---------|---------|-----------|------|
| **scikit-learn** | numpy | Simple API, basic kernels | First experiments, baseline |
| **GPy** | numpy | Mature, many kernels, multi-output | Standalone GP research |
| **GPyTorch** | PyTorch | Scalable (Whitening), deep kernels | Large N, NN-GP hybrid |
| **GPflow** | TensorFlow | Variational, large-scale | Multi-task, latent variable |
| **Pyro / NumPyro** | PyTorch / JAX | Bayesian model integration | GP as building block |
| **BoTorch** | PyTorch | BO-focused, q-batch, multi-fidelity | Production BO |
| **Ax** | BoTorch | Experiment framework on top | Multi-experiment management |
| **Stan** | C++ | Not specifically GP, but works | Bayesian, when GP is one of many components |

## Skeleton: GPy

```python
import GPy
import numpy as np

# 1D input, 1D output
X = np.random.uniform(0, 10, (20, 1))
Y = np.sin(X) + 0.1*np.random.randn(20, 1)

kernel = GPy.kern.RBF(input_dim=1, variance=1., lengthscale=1.)
gpr = GPy.models.GPRegression(X, Y, kernel)
gpr.optimize_restarts(num_restarts=5)
gpr.plot()
```

## Skeleton: GPyTorch (scalable)

```python
import torch
import gpytorch

class ExactGPModel(gpytorch.models.ExactGP):
    def __init__(self, X_train, y_train, likelihood):
        super().__init__(X_train, y_train, likelihood)
        self.mean = gpytorch.means.ConstantMean()
        self.covar = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.RBFKernel(ard_num_dims=X_train.shape[-1])
        )
    def forward(self, x):
        return gpytorch.distributions.MultivariateNormal(
            self.mean(x), self.covar(x)
        )

likelihood = gpytorch.likelihoods.GaussianLikelihood()
model = ExactGPModel(X_train, y_train, likelihood)
mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)
# Train via Adam optimizer on -mll
```

## Skeleton: BoTorch BO

```python
import torch
import botorch
from botorch.models import SingleTaskGP
from botorch.fit import fit_gpytorch_mll
from botorch.acquisition import qExpectedImprovement
from botorch.optim import optimize_acqf

# Train GP
gp = SingleTaskGP(X_train, y_train)
fit_gpytorch_mll(gp.likelihood.mll(gp))

# Propose next point via qEI
acqf = qExpectedImprovement(gp, best_f=y_train.max())
bounds = torch.tensor([[0.]*d, [1.]*d])
candidate, _ = optimize_acqf(acqf, bounds=bounds, q=1, num_restarts=10,
                               raw_samples=512)
```

## Practical considerations

- **Always normalize inputs** to [0, 1] or zero-mean unit-var
  before GP fitting — kernel hyperparameter optimization is much
  more stable.
- **Normalize outputs** (subtract mean, divide by std).
- **Restart optimization** 5-10 times from different random kernel
  hyperparameter inits — marginal log-likelihood is multimodal.
- **Check kernel choice**: Matern 5/2 is the standard for BO; RBF
  assumes smoother functions which is too restrictive.
- **For > 1000 samples**: use SVGP (variational, GPyTorch) or sparse
  GP (inducing points).

## Cross-references

- official `optuna/optuna-mcp` — Optuna integration
- `radia_mcp.pinn` — NN-based alternatives for high-dim / nonlinear
"""


LLM_BO = r"""
# LLM-Enhanced Bayesian Optimization

`Large Language Models to Enhance Bayesian Optimization.pdf` (in
public-safe curated corpus — flagged here because BO-relevant).

## Key idea

Use an LLM (GPT-4 / Claude) as a domain-aware "guesser" for BO:

1. Sample initial random points; evaluate
2. Pass (X, y, context) to LLM as a prompt: "Given these observations
   of an electric motor design objective, where would you suggest
   trying next?"
3. LLM proposes candidate(s); evaluate; add to data
4. Mix LLM proposals with standard BO acquisition (alternate or
   ensemble)

## Why it might work

- LLM has internalized published design heuristics + similar problem
  solutions from training data
- For interpretable design variables (e.g. "n_turns", "core radius"),
  the LLM can use prior knowledge
- Especially valuable in the first 10-50 trials (cold start)

## Why it might fail

- LLM hallucinates plausible-sounding but wrong values
- Output not reproducible (temperature > 0)
- Doesn't generalize to truly novel design spaces
- No theoretical guarantee (vs BO regret bounds)

## Sane integration pattern

```python
import openai  # or anthropic

def llm_propose(X, y, context_text):
    \"\"\"Use an LLM to propose next design point.\"\"\"
    prompt = f\"\"\"{context_text}

Observed designs and outcomes:
{format_table(X, y)}

Propose ONE new design as a JSON dict matching the variable names
above. Justify briefly.\"\"\"
    response = openai.ChatCompletion.create(
        model="gpt-4", messages=[{"role": "user", "content": prompt}])
    return parse_json(response.choices[0].message.content)

# In Optuna or custom loop:
def hybrid_objective(trial):
    if trial.number % 3 == 0:
        # LLM proposal every 3rd trial
        x_dict = llm_propose(X_so_far, y_so_far, context)
        for name, val in x_dict.items():
            trial.set_user_attr(f"suggested_{name}", val)
        # Force trial to use these values
        x = {name: trial.suggest_float(name, val*0.99, val*1.01)
             for name, val in x_dict.items()}
    else:
        x = standard_suggest(trial)
    return run_fem(x)
```

## Lab-specific opportunities

- IH coil design context: LLM has seen many coil topology papers
  → reasonable starting points
- Motor design: lots of training data on PMSM/SynRM choices
- WPT: known compensation topologies are well-documented

For truly novel research (new material, new physics regime), LLM
will hallucinate and bias the search — use sparingly.

## Cross-references

- official `optuna/optuna-mcp` / plain Optuna API — enqueue_trial integration
- `radia_mcp.bayesian_opt(topic="bayesian_optimization")` — for
  pure BO
"""


def get_bayesian_opt_documentation(topic: str = "all") -> str:
    """Dispatch Bayesian opt / GP / FMQA / surrogate topics.

    Topics:
        gp_fundamentals          - GP regression, kernels, ARD, MLL
        bayesian_optimization    - BO algorithm + acquisition functions
        fmqa                     - Factorization Machines + QA for binary BO
        physics_informed_gp      - PI-GP (Raissi 2017, Pförtner 2023)
        surrogate_models         - When/how to use; CAE懇話会 workshop
        gp_libraries             - GPy/GPyTorch/BoTorch/Ax skeletons
        llm_bo                   - LLM-enhanced BO (experimental)
        all                      - Everything
    """
    topic = topic.lower().strip()
    if topic in ("gp_fundamentals", "gp", "gpr", "kernels", "ard"):
        return GP_FUNDAMENTALS
    if topic in ("bayesian_optimization", "bo", "acquisition"):
        return BAYESIAN_OPTIMIZATION
    if topic in ("fmqa", "factorization_machine", "quantum_annealing",
                 "qubo"):
        return FMQA
    if topic in ("physics_informed_gp", "pi_gp", "raissi", "pde_gp"):
        return PHYSICS_INFORMED_GP
    if topic in ("surrogate_models", "surrogate", "cae_surrogate"):
        return SURROGATE_MODELS
    if topic in ("gp_libraries", "libraries", "gpytorch", "botorch", "gpy"):
        return GP_LIBRARIES
    if topic in ("llm_bo", "llm_enhanced", "llm"):
        return LLM_BO
    if topic == "all":
        return "\n\n".join([GP_FUNDAMENTALS, BAYESIAN_OPTIMIZATION, FMQA,
                            PHYSICS_INFORMED_GP, SURROGATE_MODELS,
                            GP_LIBRARIES, LLM_BO])
    return (f"Unknown topic '{topic}'. Available: gp_fundamentals, "
            "bayesian_optimization, fmqa, physics_informed_gp, "
            "surrogate_models, gp_libraries, llm_bo, all.")
