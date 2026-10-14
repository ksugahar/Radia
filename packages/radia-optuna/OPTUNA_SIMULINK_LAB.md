# Radia Optuna Simulink laboratory

This laboratory is designed for iteration, not for drawing an optimizer out
of signal lines. The public `Optuna Study` block has two inputs, four scalar
convenience outputs, and one typed monitor Bus; the complete
six-input/eighteen-output runtime stays inside the masked subsystem.

| Port | Meaning |
|---|---|
| input `start` | Rising edge starts or resumes the saved study. |
| input `cancel` | Rising edge cancels the current run. |
| output `best` | Best objective value so far. |
| output `status` | `0` idle, `1` complete, `2` running, `3` cancelled, `4` paused, `-1` failure. |
| output `progress` | Number of attempted trials. |
| output `best trial` | Zero-based trial number of the current best result. |
| output `monitor` | Fixed-schema `OptunaMonitorBusV1` for Bus Selector, Scope, SDI, and logging. |

Sampler, seed, pruner, bounds, total trial budget, objective, and MAT study
path are mask parameters. Pareto points, failed/pruned states, constraints,
intermediate values, and all parameter values remain in the normalized study
tables instead of becoming permanent top-level wiring. The monitor Bus carries
only fixed scalar telemetry, so increasing the trial budget or Pareto-front
size never changes model topology.

## 1. Open the known-optimum exercise

```matlab
addpath("matlab")
radia.simulink.buildOptunaTeachingModel(Exercise="quadratic")
open_system("radia_optuna_teaching")
```

Run the model. The constrained quadratic has its unconstrained minimum at
`x = 0.25`. The four visible scalar signals are sufficient to tell whether the
run completed and which trial is currently best. Use a Bus Selector on
`monitor` when you also want attempted/completed/pruned/failed counts, the
current trial, the most recent value, or Pareto revisions.

Double-click `Optuna Study` and press **Review saved study**, or run:

```matlab
review = radia.simulink.reviewOptunaStudy( ...
    "radia_optuna_teaching/Optuna Study");
review.trials
review.parameters
```

The mask button also assigns `radia_optuna_review` and
`radia_optuna_trials` in the base workspace.

## 2. Perform a fair comparison without rewiring

A different sampler, seed, pruner, or search space is a different experiment.
Give it a new **Study MAT file** so its history is not mixed with the first
experiment. Change exactly one factor, run again, and compare the two tables.

For example, change sampler from `tpe` to `qmc`, seed from `20260829` to `17`,
and use a new file only if that combination is the experiment being tested.
For a one-factor comparison, keep every other field unchanged.

To add trials to the same experiment, keep its MAT file and configuration,
increase **Total trial budget**, and run again. The completed history is
preserved and only the additional trials execute. Neither operation changes a
signal line.

An empty seed draws fresh private entropy for a new study, matching upstream
`seed=None`. The resolved seed and sampler state are stored in the MAT study;
extending that same study reuses them, so the continuation is exact. Changing
to a different explicit seed requires a new MAT study and fails loudly if the
old path is reused.

Record for every comparison:

- MATLAB and `optuna_mex` versions;
- sampler, seed, pruner, ordered search space, and total budget;
- study MAT path;
- proposal order, trial states, best-so-far history, and selected result.

An equal seed implies an equal sequential proposal sequence only when sampler
options, ordered search space, and prior history are also equal.

## 3. Select and apply a result

Set **Trial to apply** in the block mask (`-1` means the best completed
feasible trial) and press **Apply selected trial**. The same operation is
available from MATLAB:

```matlab
values = radia.simulink.applyOptunaTrial( ...
    "radia_optuna_teaching/Optuna Study", -1);
sim("radia_optuna_teaching")
```

The selected parameter values are applied to the model workspace when the
variables exist there, otherwise to the base workspace. Applying a trial is
explicit; reviewing or monitoring never changes the model.

## 4. Pareto selection

```matlab
radia.simulink.buildOptunaTeachingModel(Exercise="pareto")
```

The objectives are `x^2` and `(x-1)^2`. Inspect `review.pareto` and the
parameter table, choose a compromise trial, apply it from the same mask, and
simulate again. No Pareto-vector signal bundle is required at the top level.

## 5. Pruned and failed trials

```matlab
radia.simulink.buildOptunaTeachingModel(Exercise="reliability")
```

This five-point brute-force exercise deliberately produces a failed region,
a pruned region, and completed trials. Inspect `review.trials` and
`review.intermediate`. Explain why failed, pruned, infeasible, and cancelled
trials are distinct states and why none should be rewritten as a poor
objective value.

## 6. Advanced control boundary

`radia.simulink.buildOptunaBlock` remains the advanced runtime interface for
models that genuinely need wired pause/resume, selected-trial/apply commands,
full Pareto telemetry, failure codes, or checkpoint revisions. It retains the
stable six-input/eighteen-output ABI. Students and ordinary application models
should start with `buildOptunaStudyBlock`; advanced wiring is an explicit
opt-in, not the default teaching surface.

Parallel steady-state runs may be faster, but their proposal order is
schedule-dependent and is not reproducible from the seed alone.
