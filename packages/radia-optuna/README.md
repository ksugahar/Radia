# radia-optuna

`radia-optuna` is the separately distributed MATLAB optimization component of
the Radia monorepo. It installs the `radia.optuna` MATLAB namespace, the
independent 20-command `optuna_mex`, and the checked Optuna 4.9.0 compatibility
contracts. It does not install or load the Radia solver, NGSolve, oneMKL, or
Cubit.

Install it by itself:

```powershell
python -m pip install radia-optuna
radia-optuna-path
radia-optuna-doctor --json
```

Or install the independently versioned release validated by Radia through:

```powershell
python -m pip install "radia[optuna]"
```

`radia[optuna]` installs the native MATLAB/Simulink package without heavy
Python numerical dependencies. Use `radia[optuna-upstream]` when GP,
scrambled-QMC, or importance parity also needs the pinned upstream Python,
SciPy, and PyTorch stack.

Add the printed directory to MATLAB, then use the upstream-shaped API:

```matlab
addpath("<output of radia-optuna-path>")
study = radia.optuna.create_study( ...
    sampler=radia.optuna.TPESampler(Seed=42));
study.optimize(@(trial) (trial.suggest_float("x", -2, 2) - 0.25)^2, 100);
```

Generic Simulink optimization is part of the standalone contract.
`radia.optuna.SimulinkRunner` configures `Simulink.SimulationInput` objects,
runs the model, extracts objectives and constraints, classifies failed trials,
and records reproducible execution metadata without loading Radia. Radia-owned
electromagnetic models and application blocks remain in the main distribution.

The wheel also ships the generic `radia.simulink.buildOptunaBlock` Level-2
MATLAB S-Function block and `radia.simulink.addOptunaMonitor`. The block runs
one trial per sample, persists the normalized study tables after every state
transition, and exposes best value, trial counts, status, Pareto points, and
failure telemetry as ordinary Simulink signals. The monitor uses Simulink Scope
and XY Graph blocks; it does not require a browser or the Radia solver.

The distribution is Windows x64 because the current native artifact is
`optuna_mex.mexw64`. Native Random/TPE/evolutionary/pruner workflows do not
start Python. Install `radia-optuna[upstream]` for features intentionally
executed through pinned upstream Python packages, including checked GP
acquisition, scrambled QMC, and parameter importance.

`LTspiceRunner`, `SheetMetalRunner`, and `internal.runLTspiceTrial` are shipped
as explicitly classified Radia integration adapters. The standalone core does
not call them; choosing one requires the full `radia` installation.

## Release boundary

CI builds and verifies a distinct `radia-optuna-wheel` artifact. A
`radia-optuna-v<version>` tag may publish only that exact CI artifact, after
rechecking its tag, version, API inventory, MEX inventory, and dependency
boundary. The first PyPI publication additionally requires the repository's
trusted publisher to be registered for the new `radia-optuna` project; until
that external registration and release tag exist, install the verified wheel
artifact rather than assuming the PyPI name is already live.
