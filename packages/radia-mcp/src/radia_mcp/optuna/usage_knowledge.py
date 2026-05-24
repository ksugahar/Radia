"""Optuna usage knowledge: basics, storage, visualization."""

OVERVIEW = r"""
# Optuna for Radia / EM analysis

## When Optuna is the right tool

Optuna is a **derivative-free / black-box** Bayesian hyperparameter
optimization framework. Use it when:

- The objective is a **non-differentiable / non-analytic** function of
  design variables -- e.g. a full FEM solve with a complicated mesh /
  material / circuit topology.
- The objective is **expensive** (one trial = one FEM solve = minutes
  to hours) so naive grid / random search is wasteful.
- The design space includes **categorical / conditional** variables
  (e.g. winding topology choice + per-topology parameters).
- You want **multi-objective Pareto** trade-offs (e.g. torque vs loss).
- The hyperparameter being tuned is the **physical design** itself
  (coil geometry, lamination shape, material grade, ...).

Do NOT use Optuna when:

- A gradient is cheaply available -- use `radia_mcp.topology_optimization`
  (SIMP, level-set, MMA), which is far more efficient.
- The objective is a closed-form formula -- use scipy.optimize directly.
- The problem is purely combinatorial with structure -- use a SAT /
  MILP solver instead.

## Position in the radia-mcp universe

| Approach | Tool | Best for |
|----------|------|----------|
| Gradient-based topology opt | `radia_mcp.topology_optimization` (SIMP/LS) | Density field over fixed mesh |
| Shape sensitivity (cheap grad) | `radia_mcp.topology_optimization` (shape_derivative) | Smooth boundary deformation |
| Black-box BBO (this MCP) | `radia_mcp.optuna` | Categorical + mixed + expensive |
| Multi-objective gradient | `radia_mcp.topology_optimization` (MMA) | Many local constraints |
| Multi-objective BBO | `radia_mcp.optuna` (NSGA-II) | Pareto front, no gradient |
| PINN / surrogate-based BBO | `radia_mcp.pinn` + `radia_mcp.optuna` | Many parameters, replace FEM |

## Textbook reference (canonical lab source)

佐野正太郎・秋葉拓哉・今村秀明・太田健・水野尚人・柳瀬利彦,
『Optuna によるブラックボックス最適化』, オーム社, 2023.2
ISBN 9784274230103
(The authors are the Preferred Networks Optuna core development team;
Akiba is the original Optuna author.)

Lab copy: `W:/03_文献・論文/04_機械学習と最適化/00_教科書/
Optunaによるブラックボックス最適化/`

The textbook covers:
- Ch 1: BBO fundamentals (no code)
- Ch 2: Getting started with Optuna (`basic_usage`, `storage`)
- Ch 3: Fluent usage (`samplers`, `multi_objective`, `constraints`,
        `pruning`, `parallelization`, `warm_start`, `visualization`)
- Ch 4: Application -- Bayesian cookie recipe (real-world demo)
- Ch 5: How Optuna works internally (`internals`)
- Ch 6: BBO algorithm theory (no code)

## Optuna's core design (from the KDD 2019 paper)

1. **Define-by-run API** -- the search space is built by calling
   `trial.suggest_*` inside the objective function, NOT declared
   upfront. This gives full Python control flow including conditional
   branches (e.g. `if clf == "RF": suggest RF params else suggest GB`).

2. **Storage abstraction** -- in-memory (default), SQLite, MySQL,
   PostgreSQL. Same study can be resumed / shared across processes by
   pointing at the same database.

3. **Sampler / Pruner separation** -- samplers propose parameters,
   pruners early-stop bad trials based on intermediate values.

4. **Distributed via storage** -- N workers all call
   `optuna.load_study(storage=...)` then `study.optimize(...)`; they
   coordinate through the shared RDB. No special framework needed.

## Why Optuna won (vs Hyperopt / Vizier / SMAC)

- Hyperopt's MongoDB-only storage was a pain.
- Vizier was Google-internal until very recently.
- SMAC is excellent but pre-deep-learning era; less Pythonic.
- Optuna grew up alongside PyTorch / LightGBM and naturally inherited
  modern Python idioms (context manager, async / parallel, plotly viz).
"""


BASIC_USAGE = r"""
# Basic usage pattern

The mental model: write an `objective(trial)` function that returns a
scalar; Optuna calls it `n_trials` times, each time choosing parameter
values via `trial.suggest_*`.

## Minimal single-objective example (Beale test function)

(Textbook list 2.4)

```python
import optuna

def objective(trial):
    x = trial.suggest_float("x", -4.5, 4.5)
    y = trial.suggest_float("y", -4.5, 4.5)
    return (1.5 - x + x*y)**2 + \
           (2.25 - x + x*y**2)**2 + \
           (2.625 - x + x*y**3)**2

study = optuna.create_study(direction="minimize")
study.optimize(objective, n_trials=1000)

print(f"Best value: {study.best_value}")
print(f"Best params: {study.best_params}")
```

## Suggest API

| Call | Type | Notes |
|------|------|-------|
| `suggest_float(name, low, high)` | continuous float | `step=` for grid, `log=True` for log-scale |
| `suggest_int(name, low, high)` | discrete int | `step=` integer grid |
| `suggest_categorical(name, [a, b, c])` | choice | values can be int / float / str / bool |
| `suggest_float(name, low, high, log=True)` | log-uniform | for learning rate, regularization, ...|
| `suggest_discrete_uniform` | (deprecated) | use `suggest_float(..., step=...)` |

## Conditional search space

The killer feature of Optuna's define-by-run is **branching**:

```python
def objective(trial):
    clf_name = trial.suggest_categorical("clf", ("RF", "GB"))
    if clf_name == "RF":
        max_depth = trial.suggest_int("rf_max_depth", 2, 32)
        min_split = trial.suggest_float("rf_min_samples_split", 0, 1)
        clf = RandomForestClassifier(max_depth=max_depth,
                                      min_samples_split=min_split)
    else:
        max_depth = trial.suggest_int("gb_max_depth", 2, 32)
        min_split = trial.suggest_float("gb_min_samples_split", 0, 1)
        clf = GradientBoostingClassifier(max_depth=max_depth,
                                          min_samples_split=min_split)
    return cross_val_score(clf, X, y, cv=3).mean()
```

Naming convention: prefix params with the branch tag (`rf_*`, `gb_*`)
so the TPE sampler treats them as separate distributions.

## ML hyperparameter example

(Textbook list 2.14, with sklearn cross_val_score)

```python
import optuna
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

def objective(trial):
    clf = RandomForestClassifier(
        max_depth=trial.suggest_int("max_depth", 2, 32),
        min_samples_split=trial.suggest_float("min_samples_split", 0, 1),
    )
    return cross_val_score(clf, X, y, cv=3).mean()

study = optuna.create_study(direction="maximize")
study.optimize(objective, n_trials=100)
```

## Direction

- `direction="minimize"` -- objective is loss / energy / time
- `direction="maximize"` -- objective is accuracy / efficiency
- For multi-objective: `directions=["minimize", "minimize"]` (see
  `multi_objective` topic)

## Stopping conditions

`optimize()` accepts EITHER `n_trials=N` OR `timeout=seconds` (or
both). Use `timeout` for production sweeps when you care about wall
clock more than trial count.

```python
study.optimize(objective, n_trials=1000, timeout=3600)  # whichever first
```

## Catch-and-continue

By default a failed trial (raised exception) aborts `optimize()`.
For long sweeps with intermittent FEM failures:

```python
study.optimize(objective, n_trials=1000,
               catch=(ValueError, RuntimeError))
```

Failed trials are stored with `state=FAIL` and ignored by samplers.

## EM-flavored example skeleton

```python
import optuna
import radia as rad
# ... your radia / NGSolve solve function

def fem_objective(trial):
    R_wp = trial.suggest_float("R_workpiece", 5e-3, 30e-3)
    g    = trial.suggest_float("gap", 1e-3, 10e-3, log=True)
    n_t  = trial.suggest_int("n_turns", 3, 15)
    freq = trial.suggest_categorical("freq_kHz", [10, 50, 200])

    result = run_ih_panel(R_wp, g, n_t, freq*1e3)  # one FEM solve
    return -result["L_total_uH"]   # maximize L

study = optuna.create_study(
    direction="minimize",
    study_name="ih_coil_design",
    storage="sqlite:///ih_design.db",
    load_if_exists=True,
)
study.optimize(fem_objective, n_trials=200, catch=(RuntimeError,))
```
"""


STORAGE = r"""
# Persistent storage and study resume

## Why storage matters

Without storage, a study lives in RAM and is lost when the Python
process exits. With storage you can:

- **Resume** a long sweep after a crash / reboot
- **Parallelize** across multiple processes (each calls `load_study`)
- **Inspect history** later (`optuna-dashboard`, custom plots)
- **Extend** the search space mid-sweep (Ch 5 use case)

## SQLite (default for local single-user)

```python
study = optuna.create_study(
    study_name="my_study",
    storage="sqlite:///my_study.db",
    direction="minimize",
    load_if_exists=True,   # reopen existing study
)
study.optimize(objective, n_trials=100)
```

`load_if_exists=True` is the safe default for "resume if there,
otherwise start fresh".

## Loading an existing study

(Textbook list 2.19)

```python
import optuna

study = optuna.load_study(
    storage="sqlite:///optuna.db",
    study_name="ch2-conditional",
)
print(f"Best value: {study.best_value}")
print(f"Best params: {study.best_params}")

# Inspect all trials
for trial in study.trials:
    print(trial.number, trial.value, trial.params)
```

## MySQL / PostgreSQL (for parallel sweeps)

(Textbook list 3.30)

SQLite locks the whole database on write -- it does NOT scale to >1
concurrent worker. For parallel sweeps use MySQL or PostgreSQL.

```python
# Worker 1 .. N (run in parallel):
study = optuna.load_study(
    study_name="ch3-parallel",
    storage="mysql+pymysql://user:test@localhost:3306/optunatest",
)
study.optimize(objective, n_trials=100)
```

Setup:
```sql
CREATE DATABASE optunatest;
CREATE USER user IDENTIFIED BY 'test';
GRANT ALL ON optunatest.* TO 'user';
```

Then one-time-only:
```python
optuna.create_study(
    study_name="ch3-parallel",
    storage="mysql+pymysql://user:test@localhost:3306/optunatest",
)
```

## Heartbeat + RetryFailedTrialCallback (crash recovery)

(Textbook list 3.35)

Long EM solves crash sometimes (OOM, license error, segfault).
Optuna can detect a stale trial (no heartbeat for `grace_period`
seconds) and retry it.

```python
from optuna.storages import RetryFailedTrialCallback, RDBStorage

storage = RDBStorage(
    url="mysql+pymysql://user:test@localhost:3306/optunatest",
    heartbeat_interval=60,          # ping every 60 s
    grace_period=120,               # mark stale if no ping for 120 s
    failed_trial_callback=RetryFailedTrialCallback(max_retry=3),
)

study = optuna.load_study(study_name="ch3-parallel", storage=storage)
study.optimize(objective, n_trials=100)
```

## Extending the search space mid-sweep

(Textbook Ch 5)

The define-by-run API lets you add new parameters in a later session
without invalidating older trials.

```python
# Session 1: just (max_depth, min_samples_split)
def objective_v1(trial):
    clf = RandomForestClassifier(
        max_depth=trial.suggest_int("max_depth", 2, 32),
        min_samples_split=trial.suggest_float("min_samples_split", 0, 1),
    )
    return cross_val_score(clf, X, y, cv=3).mean()

study = optuna.create_study(direction="maximize",
                             storage="sqlite:///optuna.db",
                             study_name="ch5-rf")
study.optimize(objective_v1, n_trials=100)
```

Then later:

```python
# Session 2: add n_estimators
def objective_v2(trial):
    clf = RandomForestClassifier(
        max_depth=trial.suggest_int("max_depth", 2, 32),
        min_samples_split=trial.suggest_float("min_samples_split", 0, 1),
        n_estimators=trial.suggest_int("n_estimators", 10, 200),
    )
    return cross_val_score(clf, X, y, cv=3).mean()

# Reuse the study -- old trials become a "warm cache" for the new sampler
study = optuna.load_study(storage="sqlite:///optuna.db",
                          study_name="ch5-rf")
study.optimize(objective_v2, n_trials=100)
```

The TPE sampler treats `n_estimators` as a new dimension and uses
the older trials' values to bootstrap conditional density estimates.

## Storage layout (Ch 5 internals)

Optuna RDB schema (simplified):

```
studies(study_id, study_name, direction)
trials(trial_id, study_id, number, state, value, datetime_*)
trial_values(trial_id, value)                  # multi-objective
trial_params(trial_id, param_name, param_value, distribution_json)
trial_intermediate_values(trial_id, step, intermediate_value)
trial_user_attributes(trial_id, key, value_json)
trial_system_attributes(trial_id, key, value_json)
```

Knowing this helps when:
- Custom queries (`select * from trials where state='COMPLETE'`)
- Cross-study comparison (read 2 study DBs, join in pandas)
- Migration / archival
"""


VISUALIZATION = r"""
# Visualization (Plotly-based)

All Optuna visualizers return a Plotly Figure -- call `.show()` to
open in browser, `.write_image()` for PNG, `.write_html()` for
interactive HTML.

## Optimization history

```python
import optuna.visualization as vis

fig = vis.plot_optimization_history(study)
fig.write_image("history.png")
```

Shows trial number vs value (with best-so-far envelope). First sanity
check after any sweep -- if it's still descending at trial N, run more.

## Parameter importance

(Textbook list 3.18)

```python
fig = vis.plot_param_importances(study)
fig.write_image("importances.png")
```

Uses a random forest internally to estimate per-parameter importance.
Tells you which knobs matter and which can be safely fixed.

## Slice plot

(Textbook list 3.19, 3.6)

```python
fig = vis.plot_slice(study, params=["x", "y", "z"])
fig.write_image("slice.png")
```

Per-parameter 1D scatter (value vs that param holding others arbitrary)
-- shows whether the search range was too wide / too narrow.

If the optimum is at a boundary -- the range was clipped too tightly.
If the value barely changes with the parameter -- it doesn't matter.

## Contour plot

```python
fig = vis.plot_contour(study, params=["x", "y"])
```

2D contour of objective over a pair of parameters -- shows interaction
effects.

## Pareto front (multi-objective)

(Textbook list 3.5)

```python
fig = vis.plot_pareto_front(
    study,
    constraints_func=lambda trial: [trial.values[0] - 0.5],  # optional
)
fig.show()
```

Plots dominated vs non-dominated trials. With `constraints_func`,
infeasible trials are shown separately.

## Intermediate values (pruning visualization)

(Textbook list 3.38)

```python
fig = vis.plot_intermediate_values(study)
```

For pruning studies: shows the per-trial learning curves and which
were pruned at which step.

## Parallel coordinate plot

```python
fig = vis.plot_parallel_coordinate(study, params=["x", "y", "z"])
```

Each trial = one polyline across the parameter axes + final value axis.
Useful for spotting which parameter combinations correlate with low /
high objective.

## EDF (empirical distribution function)

```python
fig = vis.plot_edf([study1, study2, study3])
```

Compare cumulative distribution of best-so-far values across multiple
studies (e.g. compare TPE vs Random vs CmaEs).

## Headless / batch use

For HPC / CI sweeps without a display:

```python
import plotly
plotly.io.kaleido.scope.default_format = "png"
fig.write_image("results.png", width=1200, height=800)
```

Requires `pip install plotly kaleido`.

## Optuna Dashboard

For live monitoring of a long sweep:

```bash
pip install optuna-dashboard
optuna-dashboard sqlite:///my_study.db
```

Opens a web UI at http://localhost:8080 with all of the above plots
auto-updating as trials complete. Useful when running an overnight
FEM sweep -- you can check progress from a laptop / phone.

## EM-flavored visualization pattern

For an IH coil design sweep:

```python
# After 200-trial sweep:
import optuna.visualization as vis

vis.plot_optimization_history(study).write_image("01_history.png")
vis.plot_param_importances(study).write_image("02_importances.png")
vis.plot_slice(study, params=["R_workpiece", "gap"]).write_image(
    "03_slice_R_gap.png")
vis.plot_contour(study, params=["R_workpiece", "gap"]).write_image(
    "04_contour_R_gap.png")
```

Drop the 4 PNGs in a paper figure. The contour plot is especially
useful for showing reviewers that the design trade-off was explored.
"""


def get_usage_documentation(topic: str = "all") -> str:
    """Dispatch usage topics.

    Topics:
        overview      - Position in radia-mcp universe, textbook ref
        basic_usage   - objective function pattern, suggest API,
                        conditional search space, ML/EM examples
        storage       - SQLite/MySQL/PostgreSQL, resume, heartbeat
        visualization - Plotly history / importance / slice / Pareto
        all           - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "intro"):
        return OVERVIEW
    if topic in ("basic_usage", "basic", "usage", "suggest"):
        return BASIC_USAGE
    if topic in ("storage", "rdb", "sqlite", "mysql", "resume"):
        return STORAGE
    if topic in ("visualization", "viz", "plot", "dashboard"):
        return VISUALIZATION
    if topic == "all":
        return "\n\n".join([OVERVIEW, BASIC_USAGE, STORAGE, VISUALIZATION])
    return (f"Unknown topic '{topic}'. Available: overview, basic_usage, "
            "storage, visualization, all.")
