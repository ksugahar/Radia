# radia-optuna 5.0 compatibility review

This review supersedes the earlier Claude Code Opus 5 review of the 4.x
implementation. Its useful defect findings remain incorporated, but its API
counts, version pins, sampler boundary, and release conclusions are historical.

- Reviewed: 2026-09-08
- Behavioral and algorithm oracle: `optuna==5.0.0`
- Official MCP stable contract: `optuna-mcp==0.2.0`
- Observed upstream MCP source: `0.3.0.dev` (not claimed as an installed release)
- SciPy data/runtime pin used by Sobol fixtures: `scipy==1.17.1`
- Release candidate version: `radia-optuna==0.2.0`

## Verdict

The active implementation now targets the Optuna 5 design directly. It does
not retain a 4.x namespace, fixture lane, public multi-objective TPE sampler,
removed integration exports, public system-attribute API, categorical-distance
option, baseline-quantile option, or old sampler-state restore shim.

The generated inventory reports 812/812 Optuna 5 public entries present. Of
these, 749 are backed by executable upstream evidence and 63 wider
Python-language or bridge entries are explicit assertions. All 401 required
MATLAB entries are executable-evidence mapped; none is closed by an assertion.
There are zero partial, unmapped, or missing entries.

That is a checked compatibility claim, not a Python binary-drop-in claim.
MATLAB objects, MAT/table persistence, Simulink blocks, native MEX execution,
and parallel scheduling remain explicit MATLAB extensions.

## Optuna 5 replacements

- `TPESampler` is the sole public TPE sampler for scalar and multi-objective
  studies. Its private multi-objective path uses uniform good-trial weights,
  the Optuna 5 multi-objective gamma rule, and the shared Parzen estimator.
- Omitted `Multivariate` uses the Optuna 5 automatic policy, `ConstantLiar`
  defaults to true, numerical bandwidth uses nearest-neighbor distances, and
  the removed categorical-distance option is absent.
- Constraints are named dictionaries set by `Trial.set_constraint`. Missing
  and empty dictionaries are feasible; trials may use different names and
  different numbers of constraints. The deprecated `ConstraintsFcn` path
  remains only with the upstream 5.0-to-7.0 future warning.
- NSGA-II and NSGA-III accept `BaseMutation`; `PolynomialMutation` implements
  the upstream seeded formula and default distribution index.
- QMC includes categorical coordinates, uses the union of concurrent RUNNING
  trial spaces before the first result, warns for conditional pending spaces,
  and freezes to the first COMPLETE/PRUNED trial afterwards.
- Parameter importance defaults to PED-ANOVA. Metric-name columns retain
  declaration order.
- Removed integrations and removed public APIs fail as absent rather than
  being simulated by MATLAB compatibility shims.

## Numeric and performance work retained from the earlier review

The earlier review correctly identified decimal-step adjustment, PRUNED TPE
history, arbitrary-name collision, next-down boundaries, ties-to-even rounding,
sparse trial lookup, PRUNED CMA-ES participation, and multi-objective history as
important compatibility risks. These remain covered by the 5.0 oracle suite.

The required MEX is intentionally not optional. The native gateway fails
loudly when absent or incompatible; there is no silent `radia_mex` or MATLAB
algorithm substitution. For the Optuna 5 TPE bandwidth change, the MEX also
reproduces NumPy 2.5.2's short-array unstable argsort network so tied
observations consume the same seeded proposal sequence.

Long timing and scaling work belongs under `validation_test/optimization`.
The active paired lane is `benchmark_optuna50_python.py` plus
`benchmark_matlab_optuna50.m`. Earlier 4.x mdx JSON remains useful historical
engineering evidence but cannot satisfy the 0.2.0 release gate.

The scalar automatic-multivariate path now updates its intersection search
space only for newly finished trials instead of rebuilding every historical
intersection on every `ask`. A same-host LAB development run on 2026-09-08
(100 trials, 11 repeats, 3 warm-ups, explicit study names) measured 783 versus
818 scalar trials/s for MATLAB versus upstream Python, 336 versus 286 grouped
conditional trials/s, and 241k versus 199k `trials_dataframe` rows/s. The
proposal checksums were identical. These are development measurements; the
published release result still requires the mdx validation lane.

## Distribution, MCP, and licensing boundary

The wheel contains the MATLAB namespace, the 21-command standalone
`optuna_mex`, the audited generic Simulink subset, and required notices. It does
not require the Radia solver, NGSolve, oneMKL, or Cubit. Three explicitly named
Radia adapters stay classified in the manifest and require Radia only when
chosen.

Shared Study/Trial/visualization MCP operation belongs to the official
`optuna/optuna-mcp` server. `radia-mcp` owns only MATLAB/Simulink/MEX health,
code generation, differential-oracle planning, performance evidence, and the
MATLAB-specific release gate. The project remains independent and unofficial,
does not use the Optuna logo, and carries the upstream MIT and SciPy/Joe--Kuo
notices.

## Test and release policy

Fast deterministic regressions live under `tests`. Long performance, scaling,
parallel-efficiency, and maximum-dimension runs live under `validation_test`.
Fixtures and coverage are regenerated from the pinned 5.0 environment and must
be byte-stable. A release requires the full short suite, isolated installed-
wheel checks, fresh 5.0 performance evidence, CI, merge, tag, PyPI publication,
and the four-machine release-quad gate.
