# MATLAB Optuna component

This directory is installed by the `radia-optuna` distribution and is the path
to add in MATLAB.

```matlab
addpath("<this-directory>")
status = radia.optuna.nativeStatus();
assert(status.gateway == "optuna_mex")
```

The package contains a separately compiled 21-command MEX gateway. It has no
Radia solver, NGSolve, oneMKL, Cubit, or Python binary dependency. Missing or
incompatible native gateways fail loudly and never redirect through
`radia_mex`.

`radia.optuna.SimulinkRunner` is part of this standalone distribution. It uses
`Simulink.SimulationInput` and the public `sim`/`parsim` APIs without loading
Radia. Radia electromagnetic models and application-specific adapters are not
part of that generic contract.

The checked Optuna 4.9.0 public inventory is closed: all 816 inventoried
symbols and public class members are present and differential-oracle mapped.
See `optuna_upstream_compatibility.json` and `optuna49_api_coverage.json` for
the exact machine-readable boundary. Optuna 4.9.0 is the common algorithmic
source of truth, while native MATLAB/MEX vectorization may execute that
algorithm without Python. This is a MATLAB API, not a Python binary drop-in.
MATLAB-only parallel
execution, MAT/table storage, Simulink operation, and Radia adapters are
extensions rather than Optuna parity evidence.

Native execution covers concurrent-RUNNING constant-liar TPE, advanced CMA-ES
modes, and deterministic unscrambled Sobol generation through 21,201
dimensions. The Sobol table is bundled data and does not require Python or
SciPy at MATLAB runtime.

Use `radia.optuna.export_study` and `radia.optuna.import_study` for an explicit
JSON handoff to the `radia-optuna-bridge` CLI and a real upstream Optuna
storage. The handoff preserves original names, constraints, metric names,
attributes, distributions, and trial state; it is not a runtime fallback.

`Study.trials_dataframe` follows the Optuna 4.9.0 `attrs` expansion and column
ordering while returning a native MATLAB `table`:

```matlab
frame = study.trials_dataframe( ...
    attrs=["number","value","params","user_attrs","state"]);
```

With `multi_index=true`, the flattened variable names remain convenient MATLAB
identifiers and the exact pandas-style two-level labels are available in
`frame.Properties.UserData.column_levels`. Multi-objective metric names are
preserved and sorted using the upstream column contract.

This is an independent, unofficial project and is not affiliated with,
sponsored by, or endorsed by Preferred Networks, Inc. or the Optuna project.
Optuna, the Optuna logo and any related marks are trademarks of Preferred Networks, Inc.
The Optuna logo is not used. See `THIRD_PARTY_NOTICES.md` in
this directory for upstream MIT license notices and provenance links.
