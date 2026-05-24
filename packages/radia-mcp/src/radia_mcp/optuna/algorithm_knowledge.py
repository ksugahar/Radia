"""Optuna algorithm knowledge: samplers, multi-objective, constraints,
pruning, internals."""

SAMPLERS = r"""
# Samplers (how parameters are proposed)

The sampler chooses where to evaluate next. Optuna ships several;
swap them via `sampler=...` in `create_study`.

## Sampler decision table

| Sampler | Best for | Notes |
|---------|----------|-------|
| **TPE** (default) | < ~1000 trials, mixed types, conditional space | Tree-structured Parzen estimator -- the Optuna default |
| **Random** | Baseline + ablation reference | Use first to make sure TPE is actually helping |
| **CmaEs** | Continuous, > ~100 dim | Covariance-Matrix Adaptation Evolution Strategy |
| **NSGA-II** | Multi-objective, ~100-1000 trials | Default for `directions=...` mode |
| **NSGA-III** | Multi-objective with > 3 objectives | Reference-direction selection |
| **GP / BoTorch** | < 100 expensive trials, continuous | Gaussian process surrogate (requires `botorch`) |
| **Grid** | Reproducibility, structured sweeps | Exhaustive over a finite grid |
| **Brute force** | Tiny discrete spaces | All combinations |

## TPE (default)

```python
study = optuna.create_study(
    sampler=optuna.samplers.TPESampler(
        n_startup_trials=10,      # random for first 10 trials
        n_ei_candidates=24,       # how many candidates to score EI on
        seed=42,
        multivariate=True,        # consider param correlations
        group=True,               # group conditional params by branch
    )
)
```

Key TPE knobs:
- `n_startup_trials` -- how many random trials before TPE kicks in.
  Default 10. Bump to 50-100 for high-dim problems.
- `multivariate=True` -- since Optuna 2.2; treats parameters jointly
  (covariance) rather than independently. Recommended.
- `group=True` -- for conditional search spaces, groups parameters
  by which branch they appear in. Almost always wanted.

## CMA-ES (for continuous, mid-dim)

(Textbook list 3.36)

```python
study = optuna.create_study(
    sampler=optuna.samplers.CmaEsSampler(seed=42)
)
study.optimize(objective, n_trials=5000)
```

CMA-ES (Hansen & Ostermeier 2001) is the gold standard for continuous
BBO when:
- Dimension 10-300
- No conditional structure
- > 100 trials affordable

It maintains a multivariate Gaussian proposal and adapts both mean
and covariance using rank-mu update.

**Limitation**: cannot use categorical / int parameters directly
(it casts to float, which may not respect grid). CMA-ES + TPE hybrid
sometimes used: TPE for categorical, CMA-ES for the floats inside.

## NSGA-II (default for multi-objective)

(Textbook list 3.8, 3.9, 3.11)

```python
sampler = optuna.samplers.NSGAIISampler(
    population_size=50,            # GA population
    mutation_prob=None,            # default 1/n_params
    crossover_prob=0.9,
    swapping_prob=0.5,
    crossover=optuna.samplers.nsgaii.UNDXCrossover(),  # real values
    constraints_func=lambda trial: trial.user_attrs["constraints"],
)
```

Crossover options (`optuna.samplers.nsgaii.*`):
- `UNDXCrossover` -- Unimodal Normal Distribution Crossover (real)
- `SPXCrossover` -- Simplex Crossover (real)
- `SBXCrossover` -- Simulated Binary Crossover (real)
- `BLXAlphaCrossover` -- Blend crossover
- `VSBXCrossover` -- Vector SBX
- `UniformCrossover` -- for mixed / categorical (default if not real)

For continuous-only multi-objective (e.g. ZDT benchmark), `UNDXCrossover`
gives smoother Pareto coverage than the default uniform crossover.

## Random (baseline)

```python
study = optuna.create_study(sampler=optuna.samplers.RandomSampler(seed=42))
```

Always run this first on a new problem. If TPE doesn't beat random by
2-5x on best-so-far at the same trial count, the search space is
likely too narrow / the problem is too smooth / your objective is
broken (always returning the same value).

## GP / BoTorch (for very expensive)

```python
study = optuna.create_study(
    sampler=optuna.integration.BoTorchSampler(
        n_startup_trials=10,
        seed=42,
    )
)
```

Requires `pip install botorch torch`. BoTorch sampler uses a
Gaussian Process surrogate + acquisition function (qNEI by default).
Useful when:
- Each trial is hours (CFD simulation, real experiment)
- < ~100 total trials affordable
- Continuous parameters only

The GP cubic-cost penalty limits scaling -- TPE is better past ~500
trials.

## QMC (Quasi-Monte-Carlo) sampler

```python
study = optuna.create_study(
    sampler=optuna.samplers.QMCSampler(qmc_type="sobol")
)
```

Low-discrepancy sequence (Sobol or Halton). Good for initial coverage
of a continuous space before switching to TPE / GP. Often used as
`n_startup_trials` warmup before TPE.

## Sampler swap mid-sweep

(Useful when warming up with Random then switching to TPE)

```python
# Sweep 1: 50 random trials
study.sampler = optuna.samplers.RandomSampler(seed=42)
study.optimize(objective, n_trials=50)

# Sweep 2: switch to TPE
study.sampler = optuna.samplers.TPESampler(seed=42)
study.optimize(objective, n_trials=200)
```

The TPE will use the random trials as its training set.

## Reproducibility

Always set `seed=N` on the sampler for reproducible sweeps. Optuna
also requires `study_name` to be deterministic (don't use timestamps).
Storage path should be deterministic too.

```python
sampler = optuna.samplers.TPESampler(seed=42)
study = optuna.create_study(
    study_name="reproducible_sweep_v1",   # fixed name
    storage="sqlite:///reproducible_v1.db",
    sampler=sampler,
)
```
"""


MULTI_OBJECTIVE = r"""
# Multi-objective optimization

Many EM design problems have competing objectives:
- Inductance vs DC resistance (WPT coil)
- Torque vs iron loss (motor)
- Field uniformity vs magnet volume (NMR / MRI shim)
- Heating rate vs coil current density (IH)

Optuna supports multi-objective via `directions=[...]` and returns
the Pareto front instead of a single best.

## Setup

(Textbook list 3.1, Binh-Korn benchmark)

```python
import optuna

def objective(trial):
    x = trial.suggest_float("x", 0, 5)
    y = trial.suggest_float("y", 0, 3)
    f1 = 4*x**2 + 4*y**2          # minimize
    f2 = (x-5)**2 + (y-5)**2      # minimize
    return f1, f2                  # tuple of objectives

study = optuna.create_study(
    directions=["minimize", "minimize"]
)
study.optimize(objective, n_trials=100)
```

Notice the plural `directions=[...]` (NOT `direction=`).

## Reading the Pareto front

```python
for trial in study.best_trials:
    print(f"#{trial.number} values={trial.values} params={trial.params}")
```

`study.best_trials` (plural!) returns ALL non-dominated trials.
`study.best_value` and `study.best_params` are NOT defined for MO --
use `study.best_trials` instead.

## ZDT1 benchmark (textbook list 3.8)

```python
import math, optuna

def f1(X): return X[0]

def f2(X):
    g = 1 + 9*sum(X[1:]) / (len(X) - 1)
    return g * (1 - math.sqrt(X[0]/g))

def objective(trial):
    X = [trial.suggest_float(f"x{i}", 0, 1) for i in range(30)]
    return f1(X), f2(X)

sampler = optuna.samplers.NSGAIISampler(
    crossover=optuna.samplers.nsgaii.UNDXCrossover()
)
study = optuna.create_study(
    directions=["minimize", "minimize"],
    sampler=sampler,
)
study.optimize(objective, n_trials=1000)

optuna.visualization.plot_pareto_front(study).show()
```

## Sampler choice for multi-objective

The default for `directions=[...]` is **NSGAIISampler** (not TPE).
Override with:

```python
sampler = optuna.samplers.NSGAIISampler(  # explicit
    population_size=50,
    crossover=optuna.samplers.nsgaii.UNDXCrossover(),
)
study = optuna.create_study(directions=[...], sampler=sampler)
```

For > 3 objectives, use NSGA-III:
```python
sampler = optuna.samplers.NSGAIIISampler(reference_points=...)
```

For BoTorch with qNEHVI acquisition (3 or fewer objectives):
```python
from optuna_integration import BoTorchSampler
sampler = BoTorchSampler()  # auto-picks qNEHVI for MO
```

## Hypervolume

Pareto-front quality is measured by **hypervolume** (volume dominated
by the Pareto front, w.r.t. a reference point worse than any trial).

```python
import numpy as np
from optuna.study._multi_objective import _normalize_objective_values
# Custom: extract Pareto values, compute hypervolume vs reference

best_values = np.array([t.values for t in study.best_trials])
ref = best_values.max(axis=0) * 1.1   # reference worse than any trial
# Use pymoo or pygmo for hypervolume computation
```

## Pareto-front visualization

```python
fig = optuna.visualization.plot_pareto_front(
    study,
    constraints_func=lambda trial: [trial.values[0] - 0.5],
)
fig.show()
```

For 3+ objectives, plot_pareto_front returns a parallel coordinate
or 3D scatter automatically.

## Pseudo-multi-objective via weighted sum (anti-pattern)

Avoid:
```python
def objective(trial):
    # ... returns -L + lambda*R  -- WRONG, hides Pareto trade-off
```

Multi-objective optimization exists precisely BECAUSE the right
weighting is unknown. Returning a tuple lets Optuna expose the
trade-off; you (the engineer) pick the operating point from the
Pareto front after the sweep.

## EM-flavored MO example

```python
def coil_mo(trial):
    n_t = trial.suggest_int("n_turns", 3, 15)
    R   = trial.suggest_float("R_coil", 5e-3, 50e-3)
    w   = trial.suggest_float("wire_mm", 0.5, 5.0)

    result = run_peec(n_t, R, w)
    L_uH   = result["L_uH"]
    R_mOhm = result["R_dc_mOhm"]
    return -L_uH, R_mOhm  # MAXIMIZE L (minus sign), MINIMIZE R

study = optuna.create_study(
    directions=["minimize", "minimize"],  # both as minimize
    storage="sqlite:///coil_mo.db",
    study_name="wpt_coil_pareto",
    sampler=optuna.samplers.NSGAIISampler(
        crossover=optuna.samplers.nsgaii.UNDXCrossover()),
)
study.optimize(coil_mo, n_trials=500, catch=(RuntimeError,))
```
"""


CONSTRAINTS = r"""
# Constraints

Optuna handles constraints via `constraints_func`. The function returns
a list `[c1, c2, ...]` where `ci <= 0` means feasible. NSGA-II then
uses Deb's constraint-domination rule:
- Feasible always dominates infeasible
- Among feasible: standard Pareto dominance
- Among infeasible: smaller total violation wins

## Pattern (textbook list 3.9)

```python
import optuna

def objective(trial):
    X = [trial.suggest_float(f"x{i}", 0, 1) for i in range(30)]
    v1 = X[0]
    v2 = compute_zdt1(X)

    # Store constraint values for the sampler to read
    trial.set_user_attr("constraints", [v1 - 0.5])

    return v1, v2

sampler = optuna.samplers.NSGAIISampler(
    crossover=optuna.samplers.nsgaii.UNDXCrossover(),
    constraints_func=lambda trial: trial.user_attrs["constraints"],
)

study = optuna.create_study(
    directions=["minimize", "minimize"],
    sampler=sampler,
)
study.optimize(objective, n_trials=1000)
```

The flow:
1. Objective writes constraints to `trial.user_attrs["constraints"]`
   (any key; convention is "constraints").
2. `constraints_func` reads them back when needed by the sampler.
3. NSGA-II uses constraint-domination during selection.

## Constraints in single-objective

For single-objective + TPE/CMA-ES, prefer:
- **Penalty method**: bake the constraint into the objective
  ```python
  return loss + 1e6 * max(0, c1)**2
  ```
- **Pruning**: if a constraint is detected mid-evaluation, raise
  `TrialPruned()` to abort cheaply (see `pruning` topic)

## Why constraints_func vs penalty

| Approach | Pros | Cons |
|----------|------|------|
| **constraints_func** (NSGA-II only) | Mathematically clean Pareto with feasibility | Multi-objective only; TPE can't use it directly |
| **Penalty in objective** | Works with any sampler | Penalty weight is a hyperparameter itself; biases samplers |
| **Pruning if infeasible** | Cheap (skip rest of solve) | Requires constraint to be detectable mid-trial |

## EM example: torque vs loss with airgap constraint

```python
def motor_mo(trial):
    slot_w  = trial.suggest_float("slot_w", 1e-3, 8e-3)
    yoke_h  = trial.suggest_float("yoke_h", 5e-3, 30e-3)
    airgap  = trial.suggest_float("airgap", 0.3e-3, 1.5e-3)

    # Constraint: yoke_h must accommodate slot_w
    # i.e. slot_w - 0.6*yoke_h <= 0
    trial.set_user_attr("constraints", [slot_w - 0.6*yoke_h])

    if slot_w >= 0.6*yoke_h:
        # Save the FEM solve -- prune infeasible early
        raise optuna.TrialPruned("Infeasible: slot_w too large")

    result = motor_fem(slot_w, yoke_h, airgap)
    return -result["torque"], result["iron_loss"]
```
"""


PRUNING = r"""
# Pruning (early termination of bad trials)

For iterative methods (training loop, Picard iteration, time-stepping
FEM) you can report intermediate values and let Optuna decide whether
to abort the trial early. This is a massive speedup when most trials
are bad.

## Pattern (textbook list 3.38, PyTorch CNN)

```python
import optuna

def objective(trial):
    # ... build model from trial.suggest_* ...

    for epoch in range(N_EPOCHS):
        train_one_epoch(model, train_loader)
        accuracy = evaluate(model, val_loader)

        # Step 1: report intermediate value
        trial.report(accuracy, epoch)

        # Step 2: ask the pruner
        if trial.should_prune():
            raise optuna.TrialPruned()

    return accuracy

study = optuna.create_study(
    direction="maximize",
    pruner=optuna.pruners.MedianPruner(
        n_startup_trials=5,
        n_warmup_steps=2,
        interval_steps=1,
    )
)
study.optimize(objective, n_trials=20)
```

## Pruner options

| Pruner | When to use |
|--------|-------------|
| **MedianPruner** | Default; prune if intermediate value worse than median at same step |
| **SuccessiveHalvingPruner** | Asynchronous halving (Li-Karnin); resource-efficient |
| **HyperbandPruner** | SHA wrapped in multiple brackets (best in low-FLOP regime) |
| **PercentilePruner(p=25)** | Prune below 25th percentile (more aggressive than median) |
| **ThresholdPruner(upper=0.5)** | Absolute threshold (use when you know "above 0.5 is bad") |
| **PatientPruner(inner, patience=3)** | Wait `patience` steps before pruning (gives slow-warmers a chance) |
| **NopPruner** | Disable pruning (debug) |

## MedianPruner knobs

- `n_startup_trials` -- pruning kicks in only after this many completed
  trials (need baseline). Default 5.
- `n_warmup_steps` -- per trial, pruning kicks in only after this
  many reported steps (don't prune at step 0). Default 0.
- `interval_steps` -- how often to check (default 1 = every step).

## "Pruning" as constraint-style abort (textbook list 3.11)

```python
def objective(trial):
    X = [trial.suggest_float(f"x{i}", 0, 1) for i in range(30)]
    v1 = X[0]
    if v1 > 0.5:
        raise optuna.TrialPruned(f"Too large: {v1}")  # cheap abort
    v2 = expensive_compute(X)
    return v2
```

Here `TrialPruned` is used not for "bad partial result" but for
"detected hard infeasibility before doing the FEM solve". The trial
goes into the PRUNED state and the sampler avoids that region in
future proposals.

## Pruning + multi-objective

Pruning is NOT supported for multi-objective optimization out of the
box -- TrialPruned in a MO study just marks the trial pruned but
doesn't influence the next NSGA-II generation. For MO, prefer
constraint-domination via `constraints_func`.

## EM example: Karl iteration early stop

```python
def ih_objective(trial):
    R_wp = trial.suggest_float("R_wp", 5e-3, 30e-3)
    g    = trial.suggest_float("gap", 1e-3, 10e-3, log=True)
    n_t  = trial.suggest_int("n_turns", 3, 15)

    # Karl iteration with intermediate Z_s reports
    Zs_history = []
    Zs_current = initial_Zs(R_wp)
    for karl_iter in range(MAX_KARL):
        result = run_outer_solve(R_wp, g, n_t, Zs_current)
        Zs_new = run_esim(result["H_t_rms"], R_wp)
        Zs_current = 0.5*Zs_new + 0.5*Zs_current
        Zs_history.append(abs(Zs_current))

        # Intermediate report: -L (since maximize L)
        trial.report(-result["L_uH"], step=karl_iter)
        if trial.should_prune():
            raise optuna.TrialPruned()

    return -result["L_uH"]   # minimize -> maximize L
```

50% of trials early-stop after 2 Karl iter; FEM hours -> minutes.
"""


WARM_START = r"""
# Warm starting: enqueue_trial vs add_trial

Use warm starts when you have a known-good baseline (from a previous
sweep, from literature, from intuition) and want Optuna to evaluate
it first or seed the sampler with completed examples.

## enqueue_trial (params only, will be EVALUATED)

(Textbook list 3.24)

```python
import optuna

def objective(trial):
    x = trial.suggest_float("x", -1, 1)
    y = trial.suggest_float("y", -1, 1)
    return x * y

study = optuna.create_study()

# Queue specific params; they'll be tried first
study.enqueue_trial({"x": 0.5, "y": -0.3})
study.enqueue_trial({"x": 0.9})   # partial; rest sampled

study.optimize(objective, n_trials=3)
```

Use case:
- You know the literature optimum for a benchmark; queue it as trial 0.
- You have a "production baseline" geometry; queue it to ensure the
  sweep doesn't make things worse.
- Partial params (omit some keys) -- omitted ones are sampled normally.

## add_trial (completed trial WITH value, no evaluation)

(Textbook list 3.27)

```python
import optuna

def objective(trial):
    x = trial.suggest_float("x", -1, 1)
    y = trial.suggest_float("y", -1, 1)
    return x * y

study = optuna.create_study()

search_space = {
    "x": optuna.distributions.FloatDistribution(-1, 1),
    "y": optuna.distributions.FloatDistribution(-1, 1),
}

# Inject already-known evaluations
study.add_trial(optuna.trial.create_trial(
    params={"x": 0.5, "y": -0.3},
    distributions=search_space,
    value=-0.15,
))
study.add_trial(optuna.trial.create_trial(
    params={"x": 0.1, "y": 0.1},
    distributions=search_space,
    value=0.0,
))

study.optimize(objective, n_trials=3)
```

Use case:
- Restart a sweep that was tracked outside Optuna (e.g. paper Table 1
  values) without re-running the FEM.
- Transfer learning: feed the optimum from a similar problem to
  bootstrap TPE / GP.
- Combine results from multiple sweeps into one study for analysis.

## Difference

| | enqueue_trial | add_trial |
|---|---|---|
| Triggers FEM solve? | YES | NO |
| Knows the value upfront? | NO | YES |
| Influences future samples? | Yes (after eval) | Yes (immediately) |
| Use when | "evaluate this" | "I already know this" |

## EM example: seed with current production design

```python
study = optuna.create_study(
    direction="minimize",
    study_name="ih_improve_production",
    storage="sqlite:///ih.db",
    load_if_exists=True,
)

# Try the current production design first (we know it works)
study.enqueue_trial({
    "R_wp":    0.020,
    "gap":     0.003,
    "n_turns": 8,
})

# Optuna will run the production design as trial 0, then explore
study.optimize(fem_objective, n_trials=200)
```

After the sweep, the production-design value (trial 0) is the
baseline that every "improvement" claim is measured against.

## Lab pattern: transfer between similar problems

```python
# Sweep 1: 50 mm workpiece radius
study_50mm = optuna.load_study(
    storage="sqlite:///ih.db", study_name="ih_50mm")

# Sweep 2: 60 mm workpiece radius -- use 50mm best as warm start
study_60mm = optuna.create_study(
    direction="minimize", study_name="ih_60mm",
    storage="sqlite:///ih.db")

best_50 = study_50mm.best_trial
study_60mm.enqueue_trial(best_50.params)
study_60mm.optimize(fem_objective_60mm, n_trials=200)
```
"""


PARALLELIZATION = r"""
# Parallelization (multi-worker sweeps)

## Three modes

| Mode | Storage | Best for |
|------|---------|----------|
| **In-process threads** | in-memory (default) | I/O-bound, GIL-released compute (numpy, FEM C++) |
| **Multiprocess on one machine** | SQLite (lock-prone) or MySQL | CPU-bound Python, many cores |
| **Distributed across machines** | MySQL / PostgreSQL | Production HPC sweeps |

## In-process threads (`n_jobs=N`)

```python
study.optimize(objective, n_trials=100, n_jobs=4)
```

Uses Python `threading`. OK if the objective releases the GIL
(numpy / scipy / pytorch / C-extension FEM calls do). Bad for pure
Python loops.

## Multiprocess via RDB (textbook list 3.30)

Worker scripts (one per process / per shell):

```python
# worker.py -- run N copies in parallel
import optuna

def objective(trial):
    # ... your expensive function ...
    return value

study = optuna.load_study(
    study_name="ch3-parallel",
    storage="mysql+pymysql://user:test@localhost:3306/optunatest",
)
study.optimize(objective, n_trials=100)  # each worker does 100
```

Launch:
```bash
# Terminal 1
python worker.py
# Terminal 2
python worker.py
# Terminal 3
python worker.py
```

The three workers all coordinate via the shared MySQL. Each acquires
a fresh trial number atomically. Total trials = 300 (each requests
100; they share the queue).

For "total 100 across all workers" use:
```python
study.optimize(objective, n_trials=100 // N_WORKERS)
```
or stop manually via `timeout=`.

## Why SQLite is bad for parallel

SQLite locks the entire database on each write. With 4+ workers you
get `database is locked` errors every few seconds.

Use SQLite only for:
- Single-worker local development
- Read-only post-hoc analysis (dashboard, plot scripts)

Use MySQL / PostgreSQL for:
- Any parallel workers
- Production / shared studies
- > 10k trials

## Crash recovery (textbook list 3.35)

```python
from optuna.storages import RetryFailedTrialCallback, RDBStorage

storage = RDBStorage(
    url="mysql+pymysql://user:test@localhost:3306/optunatest",
    heartbeat_interval=60,   # workers ping every 60 s
    grace_period=120,        # if no ping in 120 s, mark trial as stale
    failed_trial_callback=RetryFailedTrialCallback(max_retry=3),
)

study = optuna.load_study(study_name="ch3-parallel", storage=storage)
study.optimize(objective, n_trials=100)
```

Behavior:
- Worker pings DB every 60 s while a trial is running.
- If a worker crashes (OOM, segfault, killed), its trial stays
  RUNNING in the DB.
- After 120 s of no pings, the trial is marked FAIL.
- `RetryFailedTrialCallback(max_retry=3)` then re-enqueues the same
  params for up to 3 retries (different worker picks it up).

**Use this for any overnight FEM sweep**. Without it, a single
license error during a trial kills the whole study.

## Distributed across machines

Same as above but the MySQL is on a shared server:

```python
storage = "mysql+pymysql://user:pw@gpu-server.lab.local:3306/optuna"
```

Workers can be on any machine that can reach `gpu-server.lab.local`.
Mix CPU machines (cheap, slow) and GPU machines (expensive, fast) --
Optuna doesn't care.

## Lab pattern: LAB + 100号機 hybrid sweep

```python
# Setup: optuna RDB on LAB (192.168.11.100 = NAS, shared)
storage = "mysql+pymysql://optuna:pw@192.168.11.100:3306/lab_optuna"

# LAB: run worker
# python sweep_worker.py --study coil_design --n-trials 100

# 100号機: run worker
# python sweep_worker.py --study coil_design --n-trials 100

# mdx (cluster): run 32 workers via SLURM
# srun -n 32 python sweep_worker.py --study coil_design --n-trials 10
```

Use `RetryFailedTrialCallback` because mdx jobs hit walltime limits.
"""


INTERNALS = r"""
# How Optuna works internally (textbook Ch 5)

Understanding the internals helps when you need to:
- Debug a stuck study
- Optimize sampler performance
- Implement a custom sampler
- Migrate / inspect storage

## Database schema

```
studies:
  study_id (PK)
  study_name
  directions (JSON array)

trials:
  trial_id (PK)
  study_id (FK)
  number (within-study counter)
  state (RUNNING, COMPLETE, PRUNED, FAIL, WAITING)
  value (single-objective)
  datetime_start, datetime_complete

trial_values:                                 # multi-objective
  trial_id (FK), value_index, value

trial_params:
  trial_id (FK), param_name, param_value (numeric internal),
  distribution_json (the suggest_* spec)

trial_intermediate_values:                    # for pruning
  trial_id (FK), step, intermediate_value

trial_user_attrs / trial_system_attrs:
  trial_id (FK), key, value_json
```

Param values are stored as floats internally:
- `suggest_float` -- direct
- `suggest_int` -- float (rounded on read)
- `suggest_categorical("x", ["a", "b", "c"])` -- internal int 0/1/2,
  mapped via the distribution

## Trial state machine

```
WAITING (enqueued via enqueue_trial)
  |
  v
RUNNING (worker picked it up)
  |
  +-- objective returns ----> COMPLETE
  +-- TrialPruned raised --> PRUNED
  +-- other exception -----> FAIL
  +-- heartbeat lost ------> FAIL (then maybe retried -> WAITING)
```

## TPE algorithm in one paragraph

For each parameter independently (or jointly with `multivariate=True`):
- Split observed trials into "good" (best gamma quantile, default 25%)
  and "bad" (the rest).
- Fit Parzen-window kernel density `l(x)` to good values,
  `g(x)` to bad values.
- Sample N candidates from `l(x)` (default N=24).
- Pick the candidate maximizing the ratio `l(x) / g(x)` (Expected
  Improvement proxy).
- That candidate is the next param value.

For categorical: replace kernels with weighted histograms.

## TPE with multivariate (the default since 2.2)

Joint density via Gaussian Process Mixture or copula. Captures
correlations between parameters. Almost always better than the
univariate default; turn it on:

```python
sampler = optuna.samplers.TPESampler(multivariate=True, group=True)
```

`group=True` -- for conditional search spaces, parameters appearing
in the same branch are kept together (don't try to fit a joint
density across branches).

## NSGA-II in 4 steps

1. Maintain population of `population_size` (default 50) parents.
2. Generate offspring via crossover (UNDX / SBX / Uniform) +
   mutation (default Gaussian noise on a random param).
3. Combine parents + offspring; sort by:
   - Pareto rank (non-dominated front first)
   - Crowding distance (within same rank, prefer less crowded)
4. Top `population_size` become next generation's parents.

`constraints_func` adds Deb's rule: feasible always dominates
infeasible regardless of objectives.

## Where the work happens

`study.optimize(objective, n_trials=N)` does:

```python
for trial_number in range(N):
    # Atomic in RDB:
    trial = storage.create_new_trial(study_id)
    trial_obj = Trial(study, trial._trial_id)
    try:
        value = objective(trial_obj)         # YOUR FUNCTION
        storage.set_trial_state_values(trial._trial_id,
                                        TrialState.COMPLETE,
                                        values=[value])
    except TrialPruned:
        storage.set_trial_state_values(trial._trial_id,
                                        TrialState.PRUNED)
    except Exception:
        storage.set_trial_state_values(trial._trial_id,
                                        TrialState.FAIL)
        raise   # unless catch=()

    # Sampler reads ALL completed trials to propose next:
    next_params = sampler.sample_independent(study, trial, ...)
```

The sampler reads the entire study history each iteration; this is
the parallelization bottleneck. With > 10k trials, samplers can be
seconds per call -- consider `CmaEsSampler` (constant memory) or
batch sampling APIs.

## Custom sampler

```python
class MyConstantSampler(optuna.samplers.BaseSampler):
    def sample_relative(self, study, trial, search_space):
        return {}   # no joint sampling

    def sample_independent(self, study, trial, name, distribution):
        return distribution.low   # always pick the minimum :)

study = optuna.create_study(sampler=MyConstantSampler())
```

Useful for debugging or for embedding domain knowledge (e.g. always
respect a physical constraint between params).
"""


def get_algorithm_documentation(topic: str = "all") -> str:
    """Dispatch algorithm topics.

    Topics:
        samplers         - TPE / Random / CmaEs / NSGA-II / GP
        multi_objective  - Pareto fronts, NSGA-II, ZDT benchmark
        constraints      - constraints_func + Deb domination
        pruning          - trial.report, MedianPruner, Hyperband
        warm_start       - enqueue_trial vs add_trial
        parallelization  - in-process / multiprocess / distributed
        internals        - DB schema, TPE/NSGA-II algorithms
        all              - Everything
    """
    topic = topic.lower().strip()
    if topic in ("samplers", "sampler", "tpe", "cmaes", "nsga"):
        return SAMPLERS
    if topic in ("multi_objective", "multi", "mo", "pareto"):
        return MULTI_OBJECTIVE
    if topic in ("constraints", "constraint", "feasibility"):
        return CONSTRAINTS
    if topic in ("pruning", "prune", "early_stop", "early_stopping"):
        return PRUNING
    if topic in ("warm_start", "warmstart", "enqueue", "add_trial"):
        return WARM_START
    if topic in ("parallelization", "parallel", "distributed",
                 "multi_worker"):
        return PARALLELIZATION
    if topic in ("internals", "internal", "schema", "how"):
        return INTERNALS
    if topic == "all":
        return "\n\n".join([SAMPLERS, MULTI_OBJECTIVE, CONSTRAINTS,
                            PRUNING, WARM_START, PARALLELIZATION,
                            INTERNALS])
    return (f"Unknown topic '{topic}'. Available: samplers, "
            "multi_objective, constraints, pruning, warm_start, "
            "parallelization, internals, all.")
