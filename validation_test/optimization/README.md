# Optimization validation workflows

Validation-class optimization workflows for CAE-style objectives. These
workflows use external Optuna, but radia-mcp itself does not vendor or serve
Optuna; install the official public `optuna` / `optuna-mcp` packages separately
when needed.

| Example | Shows | Capabilities used |
|---|---|---|
| [`validation_optuna_waveguide_slab.py`](validation_optuna_waveguide_slab.py) | Validation-class Optuna study for a waveguide dielectric slab: enqueue analytic half-wave candidates, minimize `|S11|`, record feasible trials | `optuna`, `waveguide_dielectric_slab_sparams`, `best_feasible_record`, `constraint_violation` |
| [`validation_optuna_waveguide_bragg_filter.py`](validation_optuna_waveguide_bragg_filter.py) | Validation-class Optuna study for a TE10 Bragg stopband filter: enqueue quarter-wave candidates, minimize `|S21|`, record feasible trials | `optuna`, `waveguide_cascade_sparams`, `best_feasible_record`, `constraint_violation` |
| [`validate_matlab_optuna_quality.m`](validate_matlab_optuna_quality.m) | Equal-budget, fixed-seed quality and ask/tell cost for MATLAB Random/TPE/CMA-ES and Random/MOTPE/NSGA-II | `radia.optuna`, Branin, correlated valley, ZDT1 |

```powershell
python validation_test\optimization\validation_optuna_waveguide_slab.py
python validation_test\optimization\validation_optuna_waveguide_bragg_filter.py
matlab -batch "addpath('validation_test/optimization'); validate_matlab_optuna_quality"
```

Run the MATLAB quality benchmark on a designated compute host. Its timing is
the sampler, table-backed ask/tell lifecycle, and analytic objective together;
it is not a field-solver benchmark. The generated JSON records the host,
runtime, trial budgets, all seeded outcomes, median regret, ZDT1 front error,
front coverage, and milliseconds per trial. Lower regret/front error and higher
coverage are better. This benchmark supports sampler selection but does not
claim bit-for-bit or sample-efficiency parity with Python Optuna.

These validation workflows contain only open, analytic waveguide formulas and
optional plain Optuna usage. Study/trial/dashboard MCP operation belongs to the
official external `optuna/optuna-mcp` server.
