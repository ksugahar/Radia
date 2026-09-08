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

The scalar automatic-multivariate path updates its intersection only for newly
finished trials and reuses identical encoded distributions. The existing native
history kernel now also serves the default sequential TPE configuration. A sole
new RUNNING trial has no observations to contribute to ConstantLiar; constraints,
PRUNED history, concurrent RUNNING trials, custom weights/gamma, and persistence
retain their general paths. Reseeding invalidates the native history. Metadata
comparison uses `isequaln` so unspecified numeric steps do not invalidate every
cached search space. Empty constraint tables avoid unnecessary column access.

The [paired LAB result](../../validation_test/optimization/results_optuna50_paired_lab_20260908.json)
records 974 versus 797 scalar trials/s for MATLAB versus Python (1.22x) and
539 versus 247 grouped conditional trials/s (2.18x). Both scripts now execute
11 prewarm workloads before the 11 measured repeats; the first three measured
repeats are discarded. This avoids mixing MATLAB JIT compilation into the
warmed claim. All proposal checksums agree within 1e-12; the executable oracle
suite passes all 74 tests. Table export measured 6.336 ms versus 5.616 ms for
1,000 rows (MATLAB/Python throughput ratio 0.886), so this is specifically a
TPE throughput improvement, not a universal speed claim.

Post-change validation: 74 upstream-oracle tests plus 76 MATLAB integration
tests passed (150 total, zero failures/incomplete), including table persistence,
parallel execution, session resume, Simulink blocks, and the teaching model.
All 11 package Python tests passed with pinned Optuna 5.0.0. The rebuilt 0.2.0
wheel passed strict source fidelity for 222 MATLAB files and 21 MEX commands.

Follow-up performance work updates the native non-group intersection at completed
trial insertion, skips `makeValidName` for already valid parameter names, and
constructs exported table columns in one call. A new upstream-generated seeded
case exercises changed bounds, removed/reintroduced parameters, and an empty
completed trial; it guards against stale intersection caches.
That case exposed three additional upstream differences, now corrected:
compatible relative proposals are reused when contained in changed bounds,
historical observations outside the current bounds are retained, and trial
ranking precedes per-parameter observation selection. The last point reuses
the existing global trial split instead of ranking only trials containing a
parameter. All 75 upstream-oracle tests pass after these corrections, including
the 24-trial changing-space sequence. The fixture is byte-stable across two
regenerations with pinned Optuna 5.0.0. The 24 table/reliability/core MATLAB
regressions and all 11 package Python tests also pass.

The [follow-up LAB record](../../validation_test/optimization/results_optuna50_followup_lab_20260908.json)
retains both fresh-Engine measurements preceding the three compatibility
corrections: scalar 1,281/1,280 trials/s, grouped
616/666 trials/s, and 1,000-row export 4.080/4.485 ms. The intervening Python
measurement is 880 scalar and 275 grouped trials/s, with 5.223 ms export.
All checksums remain within 1e-12. Python also improved relative to the preceding
record, so the apparent before/after gain includes host-load variation.
Historical 4.9 LAB scalar time was 72.717 ms per 100 trials (about 1,375 trials/s),
versus 78.040/78.105 ms now. This is a regression-investigation baseline, not a
controlled version comparison: prewarm settings differ and grouped proposals
changed in 5.0. Recovering historical performance still requires matched-host
validation; exceeding Python alone is not sufficient evidence.

The final-code rerun (after all three compatibility corrections) recorded
1,144 scalar and 473 grouped trials/s versus Python 874/282. Table export was
5.191 ms versus Python 4.885 ms. These slower MATLAB timings are retained in
the same record with separate final-source hashes; earlier peak throughput
must not be presented as the final-code acceptance result. Neither stable
grouped improvement nor historical 4.9 performance recovery is established.

MATLAB Engine 26.1 is installed on mdx2. The official dedicated Engine
startup/calculation/shutdown diagnostic passed on both mdx runner accounts in
[run 34210024491](https://github.com/ksugahar/Radia/actions/runs/34210024491).
SSH startup timed out with both installed and bundled Engine, and owned probe
processes were reaped. There is no mdx2 MATLAB timing result from these attempts;
fresh mdx performance validation remains a release requirement.

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
