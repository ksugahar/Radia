# Optuna 5.0 binary64 differential precision audit

## Conclusion

The exercised TPE and CMA-ES proposals are numerically very close to the
upstream oracle, but are **not bit-identical and not uniformly within 1 ULP**.
The largest absolute difference observed was 4.440892098500626e-15 (advanced
CMA-ES). Do not replace this measured statement with a universal machine-
precision compatibility claim. Existing 5e-12 acceptance tolerances were not
changed by this audit, and are not the measured errors.

| Exercised family | Maximum absolute difference | Maximum reference-relative difference | Maximum ULP distance |
|---|---:|---:|---:|
| TPE, including grouped, constrained/running, pruned and multi-objective cases | 2.8866e-15 | 1.2525e-13 | 962 |
| CMA-ES, including advanced modes | 4.4409e-15 | 1.8025e-13 | 1280 |
| NumPy RNG seed contract and RandomSampler sequence | 0 | 0 | 0 |
| NSGA-II/III baseline seeded proposal sequences | 0 | 0 | 0 |
| NSGA-II crossover cases | 6.9389e-18 | 3.0366e-16 | 2 |
| QMC cases | 0 | 0 | 0 |

Each column is an independent maximum; it need not come from the same pair.
For example, the largest CMA-ES ULP distance compares approximately
0.0123184109424432 (MATLAB) with 0.0123184109424410 (oracle): absolute
difference 2.2204e-15, relative difference 1.8025e-13, 1280 ULP. Small absolute
differences can span many representable numbers at a small magnitude. This
does not demonstrate a causal source of the discrepancy; arithmetic-order,
library, oracle serialization and algorithm-path attribution require separate
investigation. No optimizer changes were made here.

## Measurement scope and limitations

The audit intercepts real-double, same-shaped `verifyEqual` assertions in the
76-test upstream suite while preserving the exact original verification and
tolerance arguments. It records all calls, not only failures or nonzero errors.
It includes numeric metadata and internal consistency checks as well as direct
fixture comparisons. Nested struct/cell values, complex/single/integer arrays,
and inequalities are not traversed. Passing string/state assertions remain
ordinary tests; they are not floating-point measurements.

The run records 1,890 numeric comparison calls with 2,641 finite elements,
including repeated persisted/resumed cases, not 2,641 independent scenarios.
532 finite comparisons differ numerically; 2,109 are exactly equal under
numeric equality. Nonfinite mismatches and zero-reference mismatches are zero.
Signed zeros are treated as equal. NaN pairs are recorded separately from
finite distances, not presented as bit equality.

GP's default `Backend="upstream-python"` also gives zero difference in the
exercised proposal tests. That is Python delegation/bridge evidence, **not**
an independent MATLAB GP numerical implementation result. Similarly, do not
generalize a family table to untested options, dimensions, histories, hosts,
parallel schedules or native versus delegated backends.

Relative error is `abs(actual-expected)/abs(expected)`. Reference-zero cases
are excluded from the relative maximum and counted separately if unequal.
ULP is exact binary64 representable-value distance using ordered uint64 keys;
the full integer distance is stored as decimal text to avoid JSON's 53-bit
numeric limit. Every maximum records its actual/expected pair and source
test/line. The diagnostic self-test covers adjacent positive/negative doubles,
subnormals across zero, nonfinite values, zero references, and distances above
2^53. It is a MATLAB-only diagnostic test, not an optimization oracle.

## Evidence and reproduction

Raw results: `results_optuna_precision_lab_20260910.json`.
MATLAB R2026a Update 3; Python 3.12.10, Optuna 5.0.0, NumPy 2.5.2,
SciPy 1.17.1, cmaes 0.13.1, PyTorch 2.11.0+cu126. Python packages were
queried from the configured pinned oracle environment. Fixture SHA256:
`1b8640596d0d3b20abe8b34e621738abefc9b7e9f81e7c85e7cb8a7b0ff7a2b3`.
The JSON retains the fixture dependency versions and measured MATLAB host/runtime.

Configure MATLAB `pyenv` to the pinned Optuna 5.0 oracle environment, add
`validation_test/optimization` to the MATLAB path, then run:

```matlab
report = validate_optuna_precision("C:/temp/optuna_precision.json");
```

This runs the diagnostic self-test and the full upstream suite, writes test
pass/fail/incomplete counts into the report, and errors if either gate fails.
The optional recorder is inactive during ordinary tests unless
`RADIA_OPTUNA_PRECISION_OUTPUT` is set. The validator restores that environment
setting on exit; it does not change the configured Python environment or
optimizer behavior. Use a fresh owned Engine and quit it after completion.

Acceptance: diagnostic self-test and all 76 upstream tests pass; 11 standalone
package tests pass. An initial development run found an out-of-date test
classification manifest while the new diagnostic test was being registered;
the regenerated manifest passes the final audit. Failed development output is
not used as the accepted evidence. This is a finite test-corpus precision audit,
not a proof for every possible input or a performance measurement.
