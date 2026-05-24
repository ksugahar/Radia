"""MCMC libraries / frameworks landscape."""

PYMC = r"""
# PyMC (default Python Bayesian framework)

`pip install pymc` (or pymc[viz])

## Why PyMC

- Pythonic API (`with pm.Model() as model:`)
- PyTensor backend → auto-diff for NUTS
- Built-in Gibbs / NUTS / Metropolis / SMC, picks best per variable
- ArviZ for diagnostics + plots
- Used by the Suyama "ベイズ機械学習入門" book code path
  (Akagi et al, Pythonでベイズ機械学習入門 → 10_統計_ベイズ_MCMC)

## Skeleton

```python
import pymc as pm
import numpy as np

with pm.Model() as model:
    # Priors
    mu_r  = pm.LogNormal("mu_r",  mu=np.log(100),  sigma=1.0)
    sigma = pm.LogNormal("sigma", mu=np.log(1e6), sigma=1.0)

    # Forward model (must be PyTensor-traceable for NUTS)
    L_pred = surrogate(mu_r, sigma)   # see radia_mcp.pinn

    # Likelihood
    L_obs = pm.Normal("L_obs", mu=L_pred, sigma=0.5e-6,
                       observed=measurement)

    # Sample
    trace = pm.sample(2000, tune=1000, chains=4,
                      target_accept=0.9, return_inferencedata=True)

# Diagnose
import arviz as az
print(az.summary(trace))   # mean, sd, hdi, r_hat, ess_bulk
az.plot_trace(trace)
```

## Non-differentiable forward (FEM via subprocess)

PyMC requires the likelihood to be a PyTensor graph. For black-box
FEM:

```python
import pytensor.tensor as pt
from pytensor.compile.ops import as_op

@as_op(itypes=[pt.dscalar, pt.dscalar], otypes=[pt.dscalar])
def fem_forward(mu_r, sigma):
    import subprocess, json
    out = subprocess.check_output([
        "python", "calc_inductance.py",
        "--mu-r", str(float(mu_r)),
        "--sigma", str(float(sigma)),
    ])
    return float(json.loads(out)["L_uH"]) * 1e-6
```

`@as_op` makes the function gradient-free → PyMC auto-falls back to
Metropolis or SMC for that block. ~minutes-to-hours per chain.

For real production: train a surrogate (`radia_mcp.pinn`) on FEM
samples, then NUTS the surrogate.

## PyMC samplers

| Sampler | When |
|---------|------|
| NUTS | Default. Differentiable continuous |
| HamiltonianMC | NUTS without auto-tuning (rare) |
| Metropolis | Continuous non-differentiable, low dim |
| BinaryMetropolis | Boolean params |
| CategoricalGibbs | Categorical params |
| SMC | Multi-modal posterior (sequential Monte Carlo) |
| Slice | Univariate, no gradient |

Mix them per-variable:

```python
with model:
    step = [pm.NUTS([mu_r, sigma]),
            pm.Metropolis([discrete_param])]
    trace = pm.sample(2000, step=step, ...)
```

## ArviZ diagnostics

```python
import arviz as az
az.summary(trace, var_names=["mu_r", "sigma"])  # tabulated R-hat / ESS
az.plot_trace(trace)
az.plot_posterior(trace, var_names=["mu_r"], hdi_prob=0.95)
az.plot_pair(trace, var_names=["mu_r", "sigma"], kind="kde")
az.plot_ppc(ppc, num_pp_samples=100)  # posterior predictive
```
"""


STAN_PYSTAN = r"""
# Stan / pystan / cmdstanpy

`pip install cmdstanpy` (or pystan but cmdstanpy is now recommended)

## Why Stan

- Hand-tuned NUTS implementation — often beats PyMC on tricky models
- Stan modeling language → clean, declarative
- Same model file works in R / Julia / Python / Stata
- Production-grade (Bayesian standard in Bayesian academia)

## Skeleton

`inductor_model.stan`:
```stan
data {
  int<lower=1> N;
  vector[N] L_obs;
}
parameters {
  real<lower=1, upper=10> log_mu_r;
  real<lower=9, upper=18.5> log_sigma;
}
transformed parameters {
  real mu_r  = exp(log_mu_r);
  real sigma = exp(log_sigma);
  vector[N] L_pred;
  for (i in 1:N) {
    L_pred[i] = surrogate_inductance(mu_r, sigma, freq[i]);
  }
}
model {
  L_obs ~ normal(L_pred, 0.5e-6);
}
```

Python driver:
```python
from cmdstanpy import CmdStanModel
model = CmdStanModel(stan_file="inductor_model.stan")
fit = model.sample(data={"N": N, "L_obs": L_obs, "freq": freq},
                    chains=4, iter_warmup=1000, iter_sampling=2000)
print(fit.summary())
```

## Stan vs PyMC

| Feature | Stan | PyMC |
|---------|------|------|
| Speed (small) | similar | similar |
| Speed (large) | usually faster | PyMC catching up |
| Black-box forward | no | yes (`@as_op`) |
| User-defined functions | yes (`functions {}` block) | yes (PyTensor) |
| Discrete params | needs marginalization | direct support |
| Learning curve | steeper (Stan lang) | gentler (Python) |
| Debugging | tougher | better tracebacks |

If you have a smooth differentiable model and want maximum speed →
Stan. If you have black-box bits or want Python everywhere → PyMC.
"""


EMCEE = r"""
# emcee (affine-invariant ensemble sampler)

`pip install emcee`

Goodman-Weare 2010 algorithm. Maintains N "walkers" in parallel,
proposes moves based on positions of OTHER walkers. Affine-
invariant: the chain dynamics are the same under any linear
reparameterization → great for elongated / correlated posteriors
without preconditioning.

## Why emcee for EM

- **Gradient-free** → works with black-box FEM forward
- **Parallel** → 32-512 walkers per chain, embarrassingly parallel
  per step (great for cluster / mdx)
- **No tuning** → no proposal scale / step size knobs
- Default in astronomy / astrophysics Bayesian inference

## Skeleton

```python
import emcee
import numpy as np
from radia.panels.calc_inductance import run_inductance_solve

def log_prob(theta, L_obs, L_std):
    mu_r, sigma = np.exp(theta)
    if not (1 < mu_r < 10000): return -np.inf
    if not (1e4 < sigma < 1e8): return -np.inf
    L_pred = run_inductance_solve(mu_r=mu_r, sigma=sigma)["L_uH"] * 1e-6
    return -0.5 * ((L_pred - L_obs) / L_std)**2

n_walkers, n_dim = 32, 2
p0 = np.log([100, 1e6]) + 0.1 * np.random.randn(n_walkers, n_dim)

sampler = emcee.EnsembleSampler(
    n_walkers, n_dim, log_prob,
    args=(100.3e-6, 0.5e-6),
    pool=mpi_pool,  # for MPI parallel
)
sampler.run_mcmc(p0, 2000, progress=True)

# Discard burn-in, flatten
samples = sampler.get_chain(discard=1000, flat=True)
print(f"mu_r:   {np.exp(samples[:,0]).mean():.0f} +/- {np.exp(samples[:,0]).std():.0f}")
print(f"sigma:  {np.exp(samples[:,1]).mean():.0e}")
```

## When emcee shines

- 2-50 dim continuous posterior
- Cheap-to-moderate forward (32 walkers x N steps = 32N FEM solves)
- Want zero hyperparameter tuning
- Have a cluster: use `from emcee.moves import StretchMove` + MPI

## When emcee fails

- High dim (>50): autocorrelation grows; HMC dominates
- Multimodal: walkers can split between modes; use parallel
  tempering (was emcee.PTSampler, now removed → use `dynesty` or
  `ptemcee`)
- Hard constraints / boundaries: walkers cluster at edges
"""


NUMPYRO_PYRO = r"""
# NumPyro / Pyro (JAX-based)

`pip install numpyro`

NumPyro is the JAX-backend version of Pyro (Pyro is PyTorch-based).
JAX gives:
- JIT compilation (often 10-100x speedup vs Stan/PyMC)
- Auto-batching, GPU/TPU support
- Automatic differentiation through arbitrary numpy-like code

## Why NumPyro for EM

- Large NN priors (Bayesian deep learning of surrogate models)
- Many parallel chains on GPU
- Composes naturally with `jax.experimental` and `equinox`

## Skeleton

```python
import jax
import jax.numpy as jnp
import numpyro
import numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS

def model(L_obs):
    log_mu_r  = numpyro.sample("log_mu_r",
                                 dist.Uniform(0., 10.))
    log_sigma = numpyro.sample("log_sigma",
                                 dist.Uniform(9., 18.5))
    mu_r, sigma = jnp.exp(log_mu_r), jnp.exp(log_sigma)

    L_pred = jax_surrogate(mu_r, sigma)   # JAX function

    numpyro.sample("obs",
                   dist.Normal(L_pred, 0.5e-6),
                   obs=L_obs)

kernel = NUTS(model, target_accept_prob=0.9)
mcmc = MCMC(kernel, num_warmup=1000, num_samples=2000,
             num_chains=4, chain_method="parallel")
mcmc.run(jax.random.PRNGKey(0), L_obs=L_obs_jnp)
mcmc.print_summary()
samples = mcmc.get_samples()
```

## When NumPyro shines

- GPU available (chain_method="vectorized" runs all chains in parallel)
- Forward model is jax-traceable (numpy code without subprocess)
- Want to chain inference with autodiff (e.g. ELBO + MCMC hybrid)

## NumPyro vs PyMC vs Stan

| | NumPyro | PyMC | Stan |
|---|---|---|---|
| Backend | JAX | PyTensor | C++ |
| GPU | yes | partial | no |
| Stan model file | no | no | yes |
| Bayesian DL | excellent | good | poor |
| Discrete params | marginalize | direct | marginalize |
| Speed (small) | fast | fast | fast |
| Speed (GPU) | ★ | OK | n/a |
"""


BLACKJAX = r"""
# BlackJAX (research-grade composable MCMC)

`pip install blackjax`

JAX-based MCMC library by deepmind / numpyro authors. Lower-level than
NumPyro — you compose the sampler from primitives (kernel, adaptation,
warmup). Useful when you want:
- A custom MCMC variant (e.g. Riemannian HMC, pCN, MALA)
- Direct control over the warmup schedule
- Stateful samplers in JAX scan

## Skeleton

```python
import jax
import blackjax

# log-density function
def logdensity_fn(x):
    return jax.scipy.stats.norm.logpdf(x, 0, 1).sum()

# Pick an inference algorithm
nuts = blackjax.nuts(logdensity_fn, step_size=0.5,
                     inverse_mass_matrix=jnp.array([1., 1.]))
state = nuts.init(jnp.zeros(2))

# Run one step
new_state, info = nuts.step(jax.random.PRNGKey(0), state)
```

Usually you chain with `jax.lax.scan` for the full chain. See
BlackJAX docs for window-adaptation patterns.

## When BlackJAX

- Implementing or experimenting with a new MCMC algorithm
- Writing a paper that needs algorithm-level control
- Otherwise prefer NumPyro / PyMC for productivity
"""


LIBRARY_DECISION_TABLE = r"""
# Library decision table

| Need | Library |
|------|---------|
| Differentiable, mid-size | **PyMC** (Python-first) |
| Differentiable, max speed | **Stan / cmdstanpy** |
| GPU acceleration | **NumPyro** |
| Black-box FEM forward, parallel | **emcee** (gradient-free + MPI) |
| Research-level algorithm experiments | **BlackJAX** |
| Variational inference | PyMC `pm.fit()`, NumPyro `SVI` |
| Sequential Monte Carlo | PyMC `pm.smc.SMC`, BlackJAX SMC |
| Nested sampling | **dynesty** (separate; not really MCMC) |

## Lab-specific recommendations

| Lab problem | Recommended library |
|-------------|---------------------|
| BH curve / Play model param ID from measured loops | PyMC + `@as_op` or emcee |
| Bayesian UQ on FEM-computed inductance / loss | emcee (cheap surrogate) or PyMC + PINN surrogate |
| Hierarchical model: many measurement temperatures | PyMC (compound step) |
| Topology likelihood given image observation (NDT) | NumPyro on GPU + image-based likelihood |
| Direct MCMC for Laplace (Shadare-Sadiku) | Custom — see em_applications topic |
"""


def get_libraries_documentation(topic: str = "all") -> str:
    """Dispatch library topics.

    Topics:
        pymc            - PyMC (default Python Bayesian)
        stan            - Stan / cmdstanpy / pystan
        emcee           - emcee (affine-invariant, gradient-free)
        numpyro         - NumPyro / Pyro (JAX backend)
        blackjax        - BlackJAX (research composable)
        decision        - When to pick which
        all             - Everything
    """
    topic = topic.lower().strip()
    if topic in ("pymc", "pm"):
        return PYMC
    if topic in ("stan", "pystan", "cmdstanpy", "cmdstan"):
        return STAN_PYSTAN
    if topic in ("emcee", "ensemble"):
        return EMCEE
    if topic in ("numpyro", "pyro", "jax_mcmc"):
        return NUMPYRO_PYRO
    if topic in ("blackjax", "black_jax"):
        return BLACKJAX
    if topic in ("decision", "decision_table", "which", "library_choice"):
        return LIBRARY_DECISION_TABLE
    if topic == "all":
        return "\n\n".join([PYMC, STAN_PYSTAN, EMCEE, NUMPYRO_PYRO,
                            BLACKJAX, LIBRARY_DECISION_TABLE])
    return (f"Unknown topic '{topic}'. Available: pymc, stan, emcee, "
            "numpyro, blackjax, decision, all.")
