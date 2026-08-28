# radia-optuna compatibility review

Review of the MATLAB Optuna component after reconciling the earlier Claude
Code Opus 5 review (`00c9576bc`, `09f4608c3`, `cbc029319`) with the current
monorepo implementation.

- Reviewed: 2026-08-28
- Behavioral and algorithm oracle: `optuna==4.9.0`
- SciPy data/runtime pin used by Sobol fixtures: `scipy==1.17.1`
- Base package version during review: `radia-optuna==0.1.4`
- Working branch: `codex/complete-radia-optuna`

## Verdict

The earlier **not complete** verdict is stale for the declared Optuna 4.9.0
compatibility surface. The generated inventory now reports **816 present, 816
oracle-mapped, zero partial, zero unmapped, zero missing**, and
`full_compatibility_complete=true`. The complete MATLAB Optuna regression set
passes **127/127**.

That conclusion has a precise scope. `radia-optuna` is a MATLAB implementation
of the Optuna API and algorithms, with explicit calls to pinned upstream Python
for the families whose implementation remains upstream-owned. It is not a
Python binary drop-in and does not claim that MATLAB-only parallel execution,
MAT/table storage, Simulink telemetry, or Radia adapters are upstream behavior.

This review does not turn a worktree result into a release claim. Full CI,
merge, versioning, tag, PyPI publication, and release-quad remain separate
gates.

## Current checked state

| Contract | Current result |
|---|---:|
| Optuna 4.9 public entries | 816 / 816 present and mapped |
| MATLAB Optuna suite | 127 / 127 passed |
| Upstream Python differential tests | 77 |
| Official upstream MCP tests | 3 |
| MATLAB integration tests | 47 |
| Standalone native gateway | 21 commands |
| Native unscrambled Sobol limit | 21,201 dimensions |
| Packaged `radia.optuna` MATLAB files | 213 |

Additional release-candidate checks performed in this review:

| Check | Result |
|---|---:|
| `packages/radia-optuna` plus focused `radia-mcp` pytest | 23 / 23 passed |
| Fresh wheel archive verification | PASS, 213 MATLAB files / 21 commands |
| Isolated installed-wheel MATLAB smoke | PASS, no repository on MATLAB path |
| Installed-wheel seeded TPE smoke | 8 / 8 trials completed |
| Installed-wheel Sobol 64 x 4,096 checksum | 131,040 |
| Native Sobol maximum-dimension validation | 21,201 dimensions, PASS |

The oracle JSON regenerates byte-for-byte with SHA-256
`212737E09EA0CBD5FD3617033AB20DF174680ED7ACDA30755F3B04F022C1774D`.
The generated API coverage, test manifest, and Sobol binary also regenerate
byte-stably.

## Disposition of the earlier findings

The decimal step handling, PRUNED TPE history, arbitrary-name collision,
`nextDown`, ties-to-even rounding, sparse trial lookup, editable layout,
PRUNED CMA-ES, and PRUNED multi-objective TPE findings are retained and covered
by the current oracle suite.

One recommendation was deliberately not imported: `NativeKernels.has` must not
make a missing or incompatible MEX look optional. The standalone optimizer
gateway is a required, exact contract and fails loudly. This matches the
published no-fallback policy and avoids silently changing numerical backends.

## Improvements made from this review

### Seeded TPE tie resolution

The final full-suite run exposed one real seeded difference in univariate
constant-liar TPE. The MT19937 state, candidate stream, good/bad split, and
categorical probabilities all matched upstream. The difference was two ULPs
introduced by the MSVC scalar `log`/`exp` path in the MEX categorical density:
an exact NumPy acquisition tie between choices A and B was no longer a tie, so
MATLAB selected a later category.

Categorical acquisition now uses the same vectorized log-sum-exp evaluation
order as the upstream NumPy oracle. Numerical Parzen sampling and density work
remain native. The concurrent constant-liar fixture now matches all 32 rows
(16 univariate and 16 multivariate), and the complete 127-test suite passes.
This is why seeded compatibility must compare proposal sequences, not merely
distribution moments or objective quality.

### History storage

The useful performance work from `09f4608c3` was ported onto the current
implementation without replacing its newer API surface:

- `freezeTrials` converts timestamp columns once per batch.
- trials without reports reuse the shared empty intermediate table.
- `TrialState.toStorage` short-circuits canonical strings.
- `IntermediateTable` is a lazy materialized view over typed columns.
- `TrialRowIndex` provides rebuildable per-trial row buckets for parameters,
  objectives, intermediate values, attributes, and constraints.

The index is a cache, never the authority. A row-count mismatch rebuilds from
the source column, so a missed append notification can cost speed but cannot
change results.

The durable long benchmark is
`validation_test/optimization/benchmark_optuna_history_store.m`. On LAB, the
250/500/1000/2000-trial validation produced a freeze exponent of 0.928; at
2,000 trials the indexed probe was 5.0 times faster than the scan reference.
All indexed results matched the scan reference. These LAB numbers are useful
for regression only; release performance should still be measured on an idle
compute host.

### Explicit Optuna storage handoff

The bridge idea from `cbc029319` was retained as an explicit separation, not a
runtime fallback:

```text
MATLAB Study -> export_study -> study-export.v1 -> radia-optuna-bridge
     ^                                                   |
     +------------- import_study <- upstream storage ----+
```

The original draft lost constraint vectors and study-level system attributes.
The revised schema preserves those, plus original parameter and attribute
names, distributions, intermediate values, timestamps, metric names, and
trial states. The Python bridge refuses any Optuna version other than 4.9.0.
It passes a real SQLite round trip; the MATLAB round trip also preserves names
that are not valid MATLAB field names.

### Previously named sampler gaps

- Constant-liar TPE now has an oracle fixture with multiple genuinely
  concurrent `RUNNING` trials.
- CMA-ES covers source-trial warm start, separable CMA, margin CMA,
  learning-rate adaptation, pruned-trial consideration, warnings, and invalid
  option combinations.
- Unscrambled Sobol uses the SciPy Joe--Kuo criterion-6 table through dimension
  21,201. A standalone MEX command generates point batches without Python.

For contiguous Sobol batches the MEX follows the single Gray-code bit change
between consecutive samples. On the same LAB host, 64 dimensions by 4,096
points measured 0.7587 ms median through the MATLAB API versus 0.8898 ms for
SciPy 1.17.1, while preserving the upstream point sequence. The separate
maximum-dimension validation generated the first point in all 21,201
dimensions successfully (0.208 s on this run).

## Distribution and licensing boundary

The fresh local `0.1.4` `py3-none-win_amd64` wheel was checked byte-for-byte
against the monorepo sources, then installed into an isolated venv. This is
verification of the worktree against the current base version, not authority to
republish that already assigned version. `radia-optuna-doctor`
resolved the wheel layout, all 213 MATLAB files, all eight standalone Simulink
entry points, the 21-command `optuna_mex`, and the third-party notices without
requiring the Radia solver package. MATLAB then loaded that installed tree with
the repository absent from its path and executed seeded TPE and native Sobol.

The distribution retains the Optuna, SciPy, and Joe--Kuo notices, identifies
itself as independent and unofficial, does not use the Optuna logo, and pins
`optuna==4.9.0` whenever an explicit upstream handoff is requested.

## Policy conclusion

Optuna 4.9.0 remains the common algorithmic source of truth: equations,
transforms, state transitions, boundary behavior, and seeded random
consumption order are shared. MATLAB vectorization, MEX kernels, batching,
parallel scheduling, MAT/table persistence, Simulink telemetry, and the
explicit storage bridge may improve performance and workflow without defining
an alternative compatibility truth.

Short regression tests remain under `tests`; long scaling and maximum-dimension
work remains under `validation_test`.

## Remaining release gates

- Run the repository CI matrix on the final commit.
- Choose the next independent `radia-optuna` version (do not republish 0.1.4)
  and synchronize the root extras, package metadata, and manifest.
- Merge, tag, publish the wheel to PyPI, and complete release-quad.

Until those gates finish, the correct statement is: **the reviewed worktree
closes the declared Optuna 4.9.0 compatibility inventory and passes its local
differential, integration, wheel, and maximum-dimension checks; it is not yet a
published release.**
