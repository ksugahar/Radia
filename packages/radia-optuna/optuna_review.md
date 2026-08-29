# radia-optuna compatibility review

Review of the MATLAB Optuna component after reconciling the earlier Claude
Code Opus 5 review (`00c9576bc`, `09f4608c3`, `cbc029319`) with the current
monorepo implementation.

- Reviewed: 2026-08-29
- Behavioral and algorithm oracle: `optuna==4.9.0`
- SciPy data/runtime pin used by Sobol fixtures: `scipy==1.17.1`
- Base package version during review: `radia-optuna==0.1.4`
- Working branch: `codex/optuna-simulink-student-workflow`

## Verdict

The earlier **not complete** verdict is stale for the required Optuna 4.9.0
compatibility scope. The generated evidence ledger reports **816/816 API
entries present**, with **748 backed by executable upstream evidence** and
**68 explicitly marked as inventory assertions**. There are zero partial,
unmapped, or missing entries. More importantly, all **400/400 required entries
have executable upstream evidence**; none of that required scope is closed by
an assertion. The ledger therefore reports
`full_compatibility_complete=true` under its checked scope rule. The complete
fast MATLAB Optuna regression set passes **148/148**.

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
| Optuna 4.9 public entries | 816 / 816 present |
| Executable upstream evidence | 748 |
| Explicit wider-inventory assertions | 68 |
| Required compatibility scope | 400 / 400 executable-evidence mapped |
| MATLAB Optuna suite | 148 / 148 passed |
| Upstream Python differential tests | 77 |
| Official upstream MCP tests | 3 |
| MATLAB integration tests | 68 |
| Standalone native gateway | 21 commands |
| Native unscrambled Sobol limit | 21,201 dimensions |
| Packaged `radia.optuna` MATLAB files | 226 |
| Standalone Simulink entry points | 10 |

Additional release-candidate checks performed in this review:

| Check | Result |
|---|---:|
| `packages/radia-optuna` plus focused `radia-mcp` pytest | 28 / 28 passed |
| Fresh wheel archive verification | PASS, 226 MATLAB files / 21 commands / 10 Simulink entries |
| Isolated installed-wheel MATLAB E2E | PASS, no repository or Radia on MATLAB path |
| Installed-wheel `OptimizationSession` save/resume | PASS, 4 / 4 trials completed |
| Installed-wheel student Simulink model | PASS, 12 / 12 trials attempted |
| Installed-wheel table resume | PASS, all seven typed tables restored |
| mdx warmed MATLAB/Python ratios | 0.670 scalar / 0.491 grouped / 0.630 table |
| mdx deterministic 4-worker batch | 2.357x sequential throughput |
| mdx 4,000-trial indexed history lookup | 5.683x scan reference |
| mdx first MEX call | 11.44 ms median over seven fresh MATLAB processes |
| radia-mcp release evidence gate | PASS, `status=ready` |
| Native Sobol maximum-dimension validation | 21,201 dimensions, PASS |
| Full Radia MEX/Simulink provenance regeneration | HIBINO, 86 / 86 passed |

The oracle JSON regenerates byte-for-byte with SHA-256
`7A561D22470DF24F8B62E7797F47A3011CD858C3CF6691B736E6C1BA4BB5FA1F`.
The API coverage SHA-256 is
`43A5F72D2961FF1797C7E9FEC537D5308FD162744D2C99118A3AC64523C05BFC`;
the test-manifest SHA-256 is
`3663ABEFEC833BD28928F0DE036DCAB8E8D09E74A4AE56B7282B5042F7282972`.
The generated coverage and manifest regenerate byte-stably.

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
(16 univariate and 16 multivariate), and the complete 148-test suite passes.
This is why seeded compatibility must compare proposal sequences, not merely
distribution moments or objective quality.

### MATLAB and Simulink student workflow

The MATLAB surface now has a toolbox-shaped entry layer:
`OptimizationParameter`, `OptimizeOptions`, `optimoptions`, `optimize`, and
`getParameterFromModel`. These names provide a familiar MATLAB workflow but do
not depend on Global Optimization Toolbox or Simulink Design Optimization.
Sampler, pruner, seed, search-space, callback, persistence, and stopping
configuration delegate to the same `radia.optuna` implementation used by the
lower-level API.

`OptimizationSession` owns the explicit configured/running/paused/completed/
cancelled lifecycle, checkpoint/restore, stale-RUNNING recovery, trial
selection, and application of selected parameters. The version-2 Level-2
MATLAB S-Function retains the original 14 output positions and appends selected
trial, pruned-trial, current-trial, and checkpoint telemetry. Its six inputs
are start, cancel, pause, resume, selected trial, and apply.

The tracked `radia_optuna_teaching.slx` and
`OPTUNA_SIMULINK_LAB.md` cover a known quadratic optimum, a biobjective Pareto
exercise, and deterministic complete/pruned/failed behavior. Both the teaching
model and the production library passed the required official-agent
read/edit/check/save/reopen lane, clean-path reopen, full-window visual QA, and
embedded-SLX scans with zero U+FFFD or suspicious `???` runs.

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
`validation_test/optimization/benchmark_optuna_history_store.m`. The mdx run
covered 250/500/1000/2000/4000 trials, produced a freeze exponent of 0.862,
and made the 4,000-trial indexed probe 5.683 times faster than the scan
reference. Every indexed result matched the independent scan result.

### mdx release-candidate performance

Pinned upstream Python and MATLAB were run consecutively on mdx with identical
100-trial workloads and 11 repeats, discarding the first three. MATLAB/Python
warmed-time ratios were 0.670 for scalar TPE, 0.491 for grouped conditional
TPE, and 0.630 for a 1,000-row `trials_dataframe`; lower is faster. Seeded
checksums and table shape matched.

The MATLAB-only deterministic batch benchmark froze RandomSampler proposals
before objective evaluation. With 64 trials, a calibrated scalar objective,
and four process workers, the warmed median improved from 5.256 s sequential
to 2.230 s parallel: 2.357x speedup and 58.9% worker efficiency. This measures
scheduler/worker throughput and does not predict the speedup of a particular
CAE solver. A 5 ms smoke was slower in parallel, documenting that cheap
objectives should remain sequential.

The required MEX first call had an 11.44 ms median over seven fresh MATLAB
processes. A separate unrelated `measured_jointconvex_ridge.py` process
remained active on one core; it was not stopped. Pre-run total CPU was
2.72--4.73% with 51,918 MiB free memory. The raw evidence records this load,
and `radia-mcp.matlab_optuna_release_gate` returned `status=ready` with no
errors.

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
republish that already assigned version. `radia-optuna-doctor` resolved the
wheel layout, all 226 MATLAB files, all ten standalone Simulink entry points,
the 21-command `optuna_mex`, and the third-party notices without requiring the
Radia solver package. MATLAB then loaded that installed tree with the
repository and Radia absent from its path and exercised Study/SimulinkRunner,
`OptimizationSession` checkpoint/resume, the version-2 optimization block, the
tracked student workflow, and table persistence.

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
closes the required Optuna 4.9.0 compatibility scope with executable upstream
evidence, exposes the entire generated surface with assertions identified
separately, and passes its differential, integration, installed-wheel, mdx
performance, and maximum-dimension checks; it is not yet a published release.**
