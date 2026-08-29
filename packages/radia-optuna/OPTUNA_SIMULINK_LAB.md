# Radia Optuna Simulink laboratory

This worksheet uses the same MATLAB-native `OptimizationSession` and
`optuna_mex` backend as the production block. No Global Optimization Toolbox,
Simulink Design Optimization, or Python process is required during a trial.

The interface deliberately follows the familiar MathWorks workflow: define
tunable parameters and bounds, choose an optimization policy, run the model,
inspect trials, and apply a selected design. Optuna 4.9.0 remains the oracle
for shared sampler and study behavior.

## 1. Known optimum

```matlab
addpath("matlab")
radia.simulink.buildOptunaTeachingModel(Exercise="quadratic")
open_system("radia_optuna_teaching")
```

Run the model. The objective is a constrained quadratic with its unconstrained
minimum at `x = 0.25`. Inspect `teaching_best`, `teaching_attempted`, and the
saved trial table in `C:\temp\radia_optuna_teaching_quadratic.mat`.

Repeat the study after changing exactly one of sampler, seed, bounds, or trial
budget. Record the change and compare proposal order, best-so-far history, and
the final selected parameter. Equal explicit seeds are meaningful only when
sampler options, ordered search space, and history are also equal.

## 2. Pareto selection

```matlab
radia.simulink.buildOptunaTeachingModel(Exercise="pareto")
```

The objectives are `x^2` and `(x-1)^2`. Neither endpoint dominates the other.
Use `pareto_count` and the trial table to choose a compromise, enter its trial
number at `Selected Trial`, pulse `Apply`, and run the model again. Monitoring
does not apply a trial automatically.

## 3. Pruned and failed trials

```matlab
radia.simulink.buildOptunaTeachingModel(Exercise="reliability")
```

This exercise uses a five-point brute-force grid. One region deliberately
fails, another is deliberately pruned after an intermediate observation, and
the remaining trials complete. Explain why failed, pruned, infeasible, and
cancelled trials are different states and why none should be rewritten as a
poor objective value.

## 4. Session controls

- `start`, `cancel`, `pause`, `resume`, and `apply` are rising-edge commands.
- `selected trial` is a zero-based completed-trial number.
- status is `0` idle, `1` complete, `2` running, `3` cancelled, `4` paused,
  and `-1` infrastructure failure.
- `checkpoint` changes whenever session state is persisted.

Before submitting the exercise, save the MAT trial table and note MATLAB,
`optuna_mex`, sampler, pruner, seed, and search-space ordering. Parallel
steady-state runs may be faster, but their proposal order is schedule-dependent
and is not reproducible from the seed alone.
