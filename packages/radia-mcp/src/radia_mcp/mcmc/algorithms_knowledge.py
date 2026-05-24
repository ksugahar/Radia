"""MCMC algorithm theory and diagnostics."""

OVERVIEW = r"""
# MCMC Overview

## What MCMC is (and isn't)

Markov Chain Monte Carlo (MCMC) constructs a Markov chain whose
**stationary distribution** is the target distribution `pi(x)`
(typically a Bayesian posterior or a Boltzmann distribution).
Running the chain long enough produces samples from `pi(x)` without
ever computing its normalization.

- It is a **sampling** algorithm, not an optimizer. (Though the
  posterior mode IS the maximum a posteriori estimate.)
- It works in arbitrary dimension because it never tries to evaluate
  the integral over the whole space.
- It needs only `pi(x)` up to a constant: `pi(x) propto p(x)`.

## When to use MCMC

| Problem | MCMC fits when |
|---------|----------------|
| Bayesian inference | Posterior is intractable analytically |
| Inverse problem | Forward model is cheap but noisy / multi-modal |
| Design under uncertainty | Need posterior over designs given measurements |
| Material parameter ID | BH curve fit, hysteresis params with multimodal landscape |
| Bayesian UQ | Want credible intervals on L / P / B given measured data |

When NOT to use MCMC:
- Single deterministic optimization → use Optuna (`radia_mcp.optuna`)
  or gradient (`radia_mcp.topology_optimization`)
- Smooth convex objective → use scipy.optimize / CG / BFGS
- Need posterior in < seconds → use variational inference instead

## Detailed balance (the heart of MCMC)

A transition kernel `K(x -> y)` has `pi` as stationary distribution
if **detailed balance** holds:

    pi(x) K(x -> y) = pi(y) K(y -> x)

This is the sufficient condition every MCMC algorithm engineers in.
The trick is choosing `K` so that the ratio `pi(y)/pi(x)` is
computable (cancels normalization).

## The MCMC family tree

```
MCMC
├── Metropolis-Hastings (1953/1970)         ← universal, simple, slow
│   ├── Random walk MH                       ← isotropic Gaussian proposal
│   ├── Independence MH                      ← proposal independent of x
│   ├── Adaptive MH (Haario 2001)            ← auto-tune proposal scale
│   └── pCN (preconditioned Crank-Nicolson)  ← for function-space priors
│
├── Gibbs sampling (Geman-Geman 1984)        ← cycle conditional updates
│   ├── Block Gibbs                          ← update groups jointly
│   └── Collapsed Gibbs                      ← integrate out nuisance vars
│
├── Hamiltonian Monte Carlo (Duane 1987)     ← exploit gradient ∇log pi
│   ├── HMC (fixed step + step count)
│   ├── NUTS (Hoffman-Gelman 2014)           ← auto-tune trajectory
│   └── Riemannian HMC (Girolami-Calderhead) ← metric-aware
│
├── Replica methods
│   ├── Parallel tempering (Geyer 1991)      ← chains at different temperatures
│   └── Simulated tempering                  ← single chain, jumps T
│
├── Sequential Monte Carlo (SMC)              ← annealed importance sampling
│   └── Particle filters                     ← for state-space models
│
└── Slice sampling (Neal 2003)               ← rarely beats HMC/NUTS now
```

## Choosing an algorithm

| You have | Use |
|----------|-----|
| Continuous params + cheap gradient | **NUTS** (PyMC, Stan, NumPyro) |
| Continuous + black-box forward (no gradient) | **emcee** (affine-invariant) or **adaptive MH** |
| Mixed discrete + continuous | **Gibbs + MH** or **PyMC** (auto-mix) |
| Multimodal posterior | **Parallel tempering** or SMC |
| Function-space prior (PDE inverse) | **pCN** (preconditioned Crank-Nicolson) |
| Bayesian deep learning | **NUTS** in small models, **SGLD/SGHMC** in large |
| Tiny posterior (1-3 dim) | Direct grid + log-sum-exp; skip MCMC |
"""


METROPOLIS_HASTINGS = r"""
# Metropolis-Hastings (the universal MCMC algorithm)

The MH algorithm: from current state `x`, propose `y ~ q(y|x)`, accept
with probability

    alpha = min(1, [pi(y) q(x|y)] / [pi(x) q(y|x)])

If accepted, next state is `y`; otherwise stay at `x`.

When the proposal is symmetric (q(y|x) = q(x|y), e.g. isotropic
Gaussian), the q factors cancel and we get **Metropolis**:

    alpha = min(1, pi(y) / pi(x))

## Minimal Python implementation

```python
import numpy as np

def mh_sample(log_pi, x0, n_steps, proposal_sigma=0.5, rng=None):
    # Random-walk Metropolis sampler.
    # log_pi: callable returning log pi(x) up to a constant
    # x0: initial state (ndarray)
    # proposal_sigma: scale of Gaussian proposal
    rng = rng or np.random.default_rng()
    d = len(x0)
    x = np.array(x0)
    log_p_x = log_pi(x)
    samples = np.zeros((n_steps, d))
    n_accept = 0
    for t in range(n_steps):
        y = x + rng.normal(0, proposal_sigma, size=d)
        log_p_y = log_pi(y)
        if np.log(rng.uniform()) < log_p_y - log_p_x:
            x, log_p_x = y, log_p_y
            n_accept += 1
        samples[t] = x
    return samples, n_accept / n_steps
```

## Proposal scale tuning (the central practical issue)

Acceptance rate target rules of thumb (Roberts-Gelman-Gilks 1997):

| Dim | Optimal accept rate |
|-----|---------------------|
| 1 | 0.44 |
| 2 | 0.35 |
| > 5 | 0.234 (asymptotic) |

If your acceptance is 0.95 → proposal too small, chain crawls.
If 0.05 → proposal too big, chain rejects everything.

Tune by adaptive MH (Haario-Saksman-Tamminen 2001):
- Estimate empirical covariance `Sigma_t` from samples so far
- Use proposal `N(x, (2.38^2/d) * Sigma_t)`
- Optionally + small isotropic for stability

## Pitfalls

- **Burn-in**: discard first ~10-50% of samples (chain converging to
  stationary distribution). Use `pymc.sample(tune=N)` or hand-truncate.
- **Mixing**: even after burn-in, consecutive samples are correlated.
  Use **autocorrelation function (ACF)** to estimate effective sample
  size (ESS).
- **Multi-modality**: random-walk MH gets stuck in one mode.
  Use parallel tempering or HMC.
- **Step size for log-scale params**: if a param has natural log
  scale (e.g. `sigma > 0`), sample `log sigma` instead.

## EM-flavored MH example

Inverse problem: find permeability `mu_r` and conductivity `sigma`
such that simulated `L_uH` matches measured 100.3 ± 0.5 uH.

```python
import numpy as np
from radia.panels.calc_inductance import run_inductance_solve

L_meas, L_meas_std = 100.3e-6, 0.5e-6

def log_posterior(theta):
    mu_r, sigma = np.exp(theta[0]), np.exp(theta[1])
    if not (1 < mu_r < 10000): return -np.inf
    if not (1e4 < sigma < 1e8): return -np.inf
    # Prior: log-uniform within bounds (constant log_p in valid region)
    L_sim = run_inductance_solve(mu_r=mu_r, sigma=sigma)["L_uH"] * 1e-6
    log_lik = -0.5 * ((L_sim - L_meas) / L_meas_std)**2
    return log_lik

samples, acc = mh_sample(log_posterior,
                          x0=np.log([100., 1e6]),
                          n_steps=10000, proposal_sigma=0.2)
print(f"Accept rate: {acc:.2f}")
mu_r_post = np.exp(samples[5000:, 0])
print(f"Posterior mu_r: {mu_r_post.mean():.0f} +/- {mu_r_post.std():.0f}")
```

Each MH step calls one FEM solve → ~1000 trials = ~hours. For
expensive forward models, use surrogate (`radia_mcp.pinn`) inside
the inner log_posterior.
"""


HAMILTONIAN_MC = r"""
# Hamiltonian Monte Carlo (HMC) and NUTS

Exploit the gradient `∇ log pi(x)` to take large, informed steps
through high-dimensional space. The state-of-the-art MCMC for smooth,
differentiable posteriors (PyMC default, Stan default).

## Physical analogy

Augment state `x` (position) with momentum `p`. Define Hamiltonian:

    H(x, p) = -log pi(x) + 0.5 * p^T M^-1 p

The leapfrog integrator simulates Hamilton's equations:

    p(t + dt/2) = p(t) + (dt/2) * ∇ log pi(x(t))
    x(t + dt)   = x(t) + dt * M^-1 * p(t + dt/2)
    p(t + dt)   = p(t + dt/2) + (dt/2) * ∇ log pi(x(t + dt))

After L leapfrog steps, accept the new state with MH probability
using `H` change as the energy.

## Why HMC is faster than RWMH

- Random-walk MH takes O(d) steps to traverse a d-dim ball: each
  step adds O(sqrt(d)) distance.
- HMC trajectories add O(sqrt(L)) distance per leapfrog (Geometric
  drift, not random walk).
- Acceptance rate stays ~0.65-0.85 if step size is tuned.

## NUTS (No-U-Turn Sampler)

NUTS auto-tunes:
1. **Step size `dt`** — via dual averaging during warmup
2. **Trajectory length `L`** — by doubling until the trajectory
   makes a U-turn (i.e. starts coming back)

Result: zero user-tunable parameters. The reason PyMC / Stan default
to NUTS.

## When HMC/NUTS shines

- High-dimensional continuous posteriors (10-10000 dim)
- Differentiable forward model
- Smooth, unimodal-ish posterior

## When HMC fails

- Discrete params → cannot diff, need Gibbs for those
- Discontinuous posterior (e.g. hard constraints) → use reflective HMC
- Severe funnel geometry (hierarchical models) → use non-centered
  parameterization
- Multi-modal posterior → still gets stuck → use parallel tempering

## PyMC example

```python
import pymc as pm
import numpy as np

with pm.Model() as model:
    # Priors (log-uniform via log-normal proxy)
    log_mu_r = pm.Uniform("log_mu_r", lower=0, upper=10)
    log_sigma = pm.Uniform("log_sigma", lower=9, upper=18.5)
    mu_r = pm.Deterministic("mu_r", pm.math.exp(log_mu_r))
    sigma = pm.Deterministic("sigma", pm.math.exp(log_sigma))

    # Forward model (must be PyMC-differentiable: use surrogate)
    L_pred = surrogate_L(mu_r, sigma)  # JAX/Aesara graph

    # Likelihood
    L_obs = pm.Normal("L_obs", mu=L_pred, sigma=0.5e-6,
                       observed=100.3e-6)

    trace = pm.sample(2000, tune=1000, chains=4)

pm.summary(trace)  # R-hat, ESS per param
```

## For non-differentiable FEM forward model

You CANNOT use NUTS directly on a Cubit → NGSolve forward (no autograd
through subprocess). Options:

1. **Surrogate replacement**: train a PINN / GP on FEM samples, then
   NUTS on the surrogate.
2. **Adaptive MH** instead (treats FEM as a black-box log-posterior).
3. **emcee**: gradient-free, parallel ensemble; works with black-box
   FEM forward via mpiexec.

Cross-link: `radia_mcp.pinn` for surrogate construction.
"""


GIBBS = r"""
# Gibbs sampling

When the joint posterior has tractable conditional distributions, cycle
through them:

    x_1 ~ pi(x_1 | x_2, x_3, ..., x_d)
    x_2 ~ pi(x_2 | x_1, x_3, ..., x_d)
    ...
    x_d ~ pi(x_d | x_1, ..., x_{d-1})

Each step is an exact sample from the conditional, so acceptance = 1.
This is a special case of MH with proposal = conditional.

## When Gibbs works

- Conditional distributions are conjugate (Normal-Normal, Beta-
  Bernoulli, Gamma-Poisson, etc.)
- Mixture models (z labels | params have categorical conditional)
- LDA topic models

## When Gibbs fails

- Correlated parameters → slow mixing (chain crawls along ridges)
  → use **block Gibbs** (joint update of correlated groups) or HMC
- No conjugate conditionals → falls back to MH within Gibbs
- Continuous PDE-constrained → HMC/NUTS usually faster

## EM-flavored Gibbs example

Hierarchical Bayesian model: measure inductance at N temperatures.
At each T, parameters (mu_r_T, sigma_T) are drawn from a population
prior with hyperparameters mu_pop, tau_pop.

```python
# Per-temperature: mu_r_T | mu_pop, tau_pop ~ Normal
# Population:     mu_pop ~ Normal(prior_mean, prior_sd)
#                 tau_pop ~ HalfCauchy
# Cycle:
#   sample mu_r_T | rest (Normal conditional, conjugate)
#   sample mu_pop | rest (Normal conditional, conjugate)
#   sample tau_pop | rest (Inverse-Gamma conditional)
```

PyMC auto-detects and uses Gibbs for conjugate parts, NUTS for the
rest -- this is the modern "compound step" sampler.
"""


PARALLEL_TEMPERING = r"""
# Parallel Tempering (replica exchange)

For **multi-modal** posteriors that trap single chains in one mode.

## Algorithm

Run M chains at temperatures T_1=1 < T_2 < ... < T_M:

- Chain m samples `pi(x)^(1/T_m)` (flatter posterior for larger T)
- Periodically propose swapping states between adjacent chains:

      accept_swap = min(1,
          (pi(x_m) / pi(x_{m+1}))^(1/T_{m+1} - 1/T_m))

- High-T chains explore freely; low-T chains stay near modes; swaps
  let modes propagate through the ladder.

## Temperature ladder design

Common: geometric `T_m = T_1 * r^m` with r ~ 1.2-2. Tune so that
adjacent swap acceptance is ~0.2-0.3.

## When PT pays off

- Posterior has 2-10 well-separated modes
- Forward model is moderate cost (you can afford 4-8 parallel chains)
- Single-chain MH gets stuck (visible as posterior "switching" with
  random seed)

## Implementation

PyMC has `pm.sample(..., init="adapt_diag", ...)` but no built-in PT;
use **emcee** with `PTSampler` (deprecated) or **dynesty** (nested
sampling, often beats PT for moderate dim).

## EM-flavored use case

Inverse problem: find magnet pole geometry that produces a target B
field. Many geometries can produce the same field (gauge / shape
non-uniqueness) → multi-modal posterior. Parallel tempering finds
all modes; single-chain MH finds one and reports it as "the answer".
"""


DIAGNOSTICS = r"""
# MCMC Diagnostics (the "did it converge?" question)

Running MCMC without diagnostics is dangerous: the chain may not have
reached stationarity. Always inspect:

## R-hat (Gelman-Rubin)

Run M ≥ 2 chains. Compute between-chain variance B and within-chain
variance W:

    R-hat = sqrt((W + B/n) / W)

R-hat → 1 as chains converge. Rules of thumb:
- R-hat < 1.01: converged
- R-hat < 1.05: probably OK
- R-hat > 1.1: NOT converged; run longer or fix model

```python
import arviz as az
az.rhat(trace)  # per-variable R-hat
```

## Effective Sample Size (ESS)

Due to autocorrelation, N samples ≠ N independent samples. ESS
estimates the "effective" count:

    ESS = N / (1 + 2 * sum(rho_k))

where rho_k is lag-k autocorrelation. Targets:
- ESS > 400 per variable for posterior mean / variance
- ESS > 1000 for posterior quantiles (e.g. 95% CI)

```python
az.ess(trace)
```

If ESS << N, chain is sticky → tune proposal, switch to HMC, or
just run longer.

## Trace plots

Visual check:
```python
az.plot_trace(trace, var_names=["mu_r", "sigma"])
```

What to look for:
- Chains overlap (well-mixed) ✓
- No drift / trend (stationary) ✓
- Looks like white noise around a constant ✓
- Two chains visit different regions ✗ (NOT converged)
- Slow oscillation visible ✗ (under-mixed)

## Autocorrelation Function (ACF)

```python
az.plot_autocorr(trace, var_names=["mu_r"])
```

Should decay to zero within ~50 lags. Slow decay = poor mixing =
need to thin (keep every k-th sample) OR run longer.

## Posterior predictive check

Even if R-hat and ESS look good, the model may be wrong. Compare
posterior predictive samples to observed data:

```python
with model:
    ppc = pm.sample_posterior_predictive(trace)
az.plot_ppc(ppc, num_pp_samples=100)
```

## Divergences (HMC/NUTS only)

NUTS reports divergences (numerical errors during leapfrog). >1%
divergences = misspecified model or pathological geometry. Try:
- Increase `target_accept` (smaller step size)
- Non-centered parameterization for hierarchical models
- Reparameterize to remove funnels

## Typical workflow

1. Run short chain (500 samples) to debug code
2. Check R-hat / ESS / trace on this short run
3. If happy, run full chain (2000-10000 samples)
4. Re-check diagnostics
5. Posterior predictive check
6. Report mean ± SD AND 95% credible interval
"""


def get_algorithms_documentation(topic: str = "all") -> str:
    """Dispatch algorithm topics.

    Topics:
        overview            - MCMC family tree + when to use
        metropolis_hastings - Random walk MH, adaptive tuning
        hamiltonian_mc      - HMC + NUTS (default in PyMC/Stan)
        gibbs               - Gibbs sampling, conjugate updates
        parallel_tempering  - Replica exchange for multi-modal
        diagnostics         - R-hat, ESS, trace, ACF, divergences
        all                 - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "intro"):
        return OVERVIEW
    if topic in ("metropolis_hastings", "mh", "metropolis", "random_walk"):
        return METROPOLIS_HASTINGS
    if topic in ("hamiltonian_mc", "hmc", "nuts"):
        return HAMILTONIAN_MC
    if topic in ("gibbs", "gibbs_sampling"):
        return GIBBS
    if topic in ("parallel_tempering", "pt", "replica_exchange", "tempering"):
        return PARALLEL_TEMPERING
    if topic in ("diagnostics", "rhat", "ess", "convergence"):
        return DIAGNOSTICS
    if topic == "all":
        return "\n\n".join([OVERVIEW, METROPOLIS_HASTINGS, HAMILTONIAN_MC,
                            GIBBS, PARALLEL_TEMPERING, DIAGNOSTICS])
    return (f"Unknown topic '{topic}'. Available: overview, "
            "metropolis_hastings, hamiltonian_mc, gibbs, "
            "parallel_tempering, diagnostics, all.")
