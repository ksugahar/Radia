# Optimization validation workflows

Validation-class optimization workflows for CAE-style objectives. These
workflows use external Optuna, but radia-mcp itself does not vendor or serve
Optuna; install the official public `optuna` / `optuna-mcp` packages separately
when needed.

| Example | Shows | Capabilities used |
|---|---|---|
| [`validation_optuna_waveguide_slab.py`](validation_optuna_waveguide_slab.py) | Validation-class Optuna study for a waveguide dielectric slab: enqueue analytic half-wave candidates, minimize `|S11|`, record feasible trials | `optuna`, `waveguide_dielectric_slab_sparams`, `best_feasible_record`, `constraint_violation` |
| [`validation_optuna_waveguide_bragg_filter.py`](validation_optuna_waveguide_bragg_filter.py) | Validation-class Optuna study for a TE10 Bragg stopband filter: enqueue quarter-wave candidates, minimize `|S21|`, record feasible trials | `optuna`, `waveguide_cascade_sparams`, `best_feasible_record`, `constraint_violation` |
| [`validate_matlab_optuna_quality.m`](validate_matlab_optuna_quality.m) | Equal-budget, fixed-seed quality gates and ask/tell cost for MATLAB Random/TPE/CMA-ES and Random/MOTPE/NSGA-II | `radia.optuna`, Branin, correlated valley, ZDT1, automatic TPE intersection, full-covariance CMA-ES |
| [`benchmark_optuna49_python.py`](benchmark_optuna49_python.py) + [`benchmark_matlab_optuna49.m`](benchmark_matlab_optuna49.m) | Same-host warmed throughput comparison with exact seeded proposal checksums | upstream `optuna==4.9.0`, scalar TPE, incremental-history grouped conditional TPE, fused native proposal kernels |
| [`benchmark_optuna_mex_cold_start.ps1`](benchmark_optuna_mex_cold_start.ps1) | Fresh-process first-call cost of the independent optimization gateway | `optuna_mex`, binary size/hash, seven-process median |
| [`validate_matlab_adjoint_quality.m`](validate_matlab_adjoint_quality.m) | Analytic state/adjoint derivative QA plus constrained MMA/SQP cross-check on a material-dependent field operator | `radia.topopt.optimizeAdjoint`, directional derivative, volume constraint, MMA, SQP |

```powershell
python validation_test\optimization\validation_optuna_waveguide_slab.py
python validation_test\optimization\validation_optuna_waveguide_bragg_filter.py
matlab -batch "addpath('validation_test/optimization'); validate_matlab_optuna_quality"
python validation_test\optimization\benchmark_optuna49_python.py --output C:\temp\optuna49_python.json
matlab -batch "addpath('validation_test/optimization'); benchmark_matlab_optuna49('C:/temp/optuna49_matlab.json')"
pwsh -ExecutionPolicy Bypass -File validation_test\optimization\benchmark_optuna_mex_cold_start.ps1 -Output C:\temp\optuna_mex_first_call.json
matlab -batch "addpath('validation_test/optimization'); validate_matlab_adjoint_quality"
```

Run the MATLAB quality benchmark on a designated compute host. Its timing is
the sampler, table-backed ask/tell lifecycle, and analytic objective together;
it is not a field-solver benchmark. The generated JSON records the host,
runtime, trial budgets, all seeded outcomes, median regret, ZDT1 front error,
front coverage, and milliseconds per trial. Lower regret/front error and higher
coverage are better. This benchmark supports sampler selection but does not
claim bit-for-bit or sample-efficiency parity with Python Optuna.

The paired Optuna 4.9 performance benchmark is a warmed, sequential,
same-host differential gate. Run both commands while the machine is otherwise
idle and compare the medians only when the host and workload settings match.
Each runner fails if its explicit-seed proposal checksum changes. The checked
result JSON records the environment, raw medians, throughput ratios, and the
claim boundary; persistence, parallel scheduling, objective cost, and cold
process startup are deliberately reported outside this shared-behavior gate.
The checked result also records the separate `optuna_mex` operational boundary:
20 commands, binary size and hash, imported libraries, and seven fresh-process
first-call measurements. This proves that optimization code is not loaded by a
non-optimization run; it does not include MATLAB executable startup time.

The adjoint validation is a fast numerical correctness workflow rather than a
performance benchmark. It checks the complete directional derivative and then
requires MMA and SQP to converge to the same feasible material design using
only the supplied analytic state/adjoint sensitivities.

These validation workflows contain only open, analytic waveguide formulas and
optional plain Optuna usage. Study/trial/dashboard MCP operation belongs to the
official external `optuna/optuna-mcp` server.
