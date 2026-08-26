# MATLAB Optuna component

This directory is installed by the `radia-optuna` distribution and is the path
to add in MATLAB.

```matlab
addpath("<this-directory>")
status = radia.optuna.nativeStatus();
assert(status.gateway == "optuna_mex")
```

The package contains a separately compiled 20-command MEX gateway. It has no
Radia solver, NGSolve, oneMKL, Cubit, or Python binary dependency. Missing or
incompatible native gateways fail loudly and never redirect through
`radia_mex`.

`radia.optuna.SimulinkRunner` is part of this standalone distribution. It uses
`Simulink.SimulationInput` and the public `sim`/`parsim` APIs without loading
Radia. Radia electromagnetic models and application-specific adapters are not
part of that generic contract.

The MATLAB API is a differentially verified subset of upstream Optuna 4.9.0,
not a claim of complete package compatibility. See
`optuna_upstream_compatibility.json` and `optuna49_api_coverage.json` for the
machine-readable boundary. MATLAB-only parallel execution, MAT/table storage,
Simulink operation, and Radia adapters are extensions rather than Optuna parity
evidence.

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
