# radia-optuna

`radia-optuna` is an independent, separately distributed MATLAB optimization
component from the Radia monorepo. It installs the `radia.optuna` MATLAB
namespace, the 20-command `optuna_mex`, and checked compatibility contracts
whose behavioral oracle is Optuna 4.9.0. It does not install or load the Radia
solver, NGSolve, oneMKL, or Cubit.

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

## Upstream attribution

This is an independent, unofficial project. It is not affiliated with,
sponsored by, or endorsed by Preferred Networks, Inc. or the Optuna project.
It does not use the Optuna logo or present itself as an official Optuna
distribution.

Optuna, the Optuna logo and any related marks are trademarks of Preferred Networks, Inc.

Optuna and the official `optuna/optuna-mcp` server are separate MIT-licensed
upstream projects and are not bundled in this wheel. See
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for their copyright and
license notices. Shared MCP Study/Trial/visualization operations remain owned
by the official server; `radia-mcp` covers only MATLAB/Simulink differences.

## Release boundary

CI builds and verifies a distinct `radia-optuna-wheel` artifact. A
`radia-optuna-v<version>` tag may publish only that exact CI artifact, after
rechecking its tag, version, API inventory, MEX inventory, and dependency
boundary, including the bundled third-party notice.

Publication also uses the standalone four-machine release-quad lane. After the
successful `main` CI run, execute
`python tools/release_quad.py optuna-candidate --ci-run-id <id> --target all`,
then pass the retained wheel to
`python tools/release_quad.py optuna-done --wheel <path>`. Only after that gate
passes may the matching commit be tagged. The manual release workflow requires
the same CI run ID and the candidate SHA256, publishes that exact wheel to
PyPI, and attaches it to the matching GitHub Release.
