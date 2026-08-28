# radia-optuna review

Review of the MATLAB Optuna 4.9.0 component (`matlab/+radia/+optuna`, the
`radia-optuna` distribution, and their fixtures), plus the work done to
resolve what it found.

- First reviewed: 2026-08-28
- Re-reviewed: 2026-08-29, on `codex/claude-optuna-20260828`
- Upstream oracle: `optuna==4.9.0` (installed and read directly, not quoted
  from documentation)
- Verification host for timings: **hibino** (idle, MATLAB R2026a, 38 cores).
  LAB timings are polluted by concurrent builds and are never quoted.

## Verdict

**Complete as an Optuna-compatible MATLAB implementation**, with the
compatibility claim now resting on derived evidence rather than assertion.

`matlab_optuna_health` reports `ok: true`. It reported `ok: false` with five
errors at the first review; all five are resolved:

| Gate error (first review) | Resolution |
|---|---|
| upstream license/trademark notices are incomplete | `THIRD_PARTY_NOTICES.md` exists |
| public API coverage is not completely oracle mapped | 400/400 required entries present, mapped and evidence-backed |
| compatibility contract does not declare complete API closure | closure now true, and it now means something checkable |
| upstream oracle SHA256 differs from the coverage record | was an editable-install artifact, not a repository defect |
| public API inventory SHA256 differs from the coverage record | same artifact |

What "complete" does **not** cover: the toolbox-shaped surface described
below is a MATLAB extension with no upstream counterpart, and 68 entries in
the bridged and out-of-scope scopes are `asserted`, not verified. Both are
visible in the ledger rather than folded into the claim.

## Current state

```
scope       required 400 | bridged 261 | out-of-scope 85 | python-language 46 | replaced 24
required    present 400 (100%)  missing 0  oracle-mapped 400  asserted 0
surface     816 entries: 748 verified (each naming its oracle sections), 68 asserted
tests       MATLAB 137;  Python 10;  radia-mcp 188
closure     full_compatibility_complete = true
```

Only entries scoped `required` count towards closure. Every other scope must
name what discharges it, so the ledger cannot quietly excuse an
unimplemented feature.

## Findings and disposition

The first review raised nine findings; eight were real and one was withdrawn.
Auditing the remaining samplers found two more of the same class. The
re-review found one further defect, and one measurement defect in the ledger
itself.

| # | Finding | Status |
|---|---|---|
| 1 | `suggestFloat` decided step divisibility in binary floating point, so `[0, 0.7]` step `0.1` warned and moved `high` above the requested bound. Upstream uses `decimal.Decimal(str(x))`. | Fixed; locked by `numeric_untransform` |
| 2 | Univariate TPE dropped PRUNED trials from both the startup gate and the observation set. Upstream reads `states=(COMPLETE, PRUNED)`. | Fixed; locked by `tpe_pruned_history_seed_113` |
| 3 | Claimed CmaEs divided by `high-low` for single-valued distributions. | **Withdrawn** — `IntersectionSearchSpace.ExcludeSingle` already defaults to true, so the NaN path is unreachable. The related real defect (a single-valued distribution reaching the sampler and consuming a random number) was fixed at the `Trial` level. |
| 4 | `Trial` silently renamed colliding parameter names while `Study.freezeTrial` recomputed keys with `makeValidName`, losing a distribution. | Fixed; locked by a collision test |
| 5 | Float untransform clamped to `high`; upstream returns `np.nextafter(high, -inf)`. | Fixed |
| 6 | MATLAB `round` (half away from zero) used where upstream uses `np.round` (half to even). | Fixed; single shared implementation |
| 7 | Vectorized observation lookup assumed a dense trial table and silently dropped rows. | Fixed with a search fallback |
| 8 | `NativeKernels.has` threw instead of returning false, making every fallback branch unreachable. | Fixed; missing gateway warns once |
| 9 | `radia_optuna.matlab_path` reported a bare FAIL for an editable install. | Fixed; `layout()` reports which tree resolved |
| 10 | `CmaEsSampler` ignored `ConsiderPrunedTrials` in the startup gate and never told the optimizer a pruned solution. | Fixed; locked by `cmaes_pruned_sampler_seed_31` |
| 11 | Multi-objective TPE used a per-parameter COMPLETE-only startup count and dropped PRUNED trials. | Fixed; locked by `multiobjective_tpe_pruned_seed_41` |
| 12 | **The ledger asserted `verified` instead of deriving it.** An entry was verified when it was present in MATLAB and its name appeared in a hand-maintained Python set — and a *module* was verified merely by being present, with no list membership at all. The test guarding the ledger demanded that all 816 entries be verified, which locked the overclaim in place. | Fixed; see below |
| 13 | The merge left `test_optuna_table.m` with two copies of the collision test. MATLAB refuses to build a suite from a file that declares the same local function twice, so the whole `test_optuna*` run aborted rather than reporting a failure. | Fixed; the surviving copy also asserts the collision precondition |

NSGA-II was audited for the same PRUNED class and is correct: upstream
`BaseGASampler` also uses `states=[COMPLETE]`.

Every new upstream fixture was confirmed to **fail** against the pre-fix code
before being accepted, so each one locks a real defect rather than describing
current behaviour.

### Finding 12: what "verified" now means

`build_oracle()` maps each section of `optuna49_oracle.json` to the function
that produces it, so the upstream names that function references are the
names that section exercises. An entry may be `verified` only when such a
section exists; one that is present and listed without any is reported as
`asserted`, and every entry records the sections behind it, so the claim is
auditable from the JSON rather than from the generator's source.

The link is a *necessary* condition, not a sufficient one — it proves an
upstream artifact covers the name, not that the artifact is exhaustive. That
is still strictly stronger than a hand-kept list, and it fails loudly when a
name is claimed with no upstream evidence at all.

Measured against the previous ledger, **28 of the 400 required entries had no
upstream artifact of any kind**:

- `Study.best_params`
- the four trial classes' v3.0-deprecated `suggest_uniform` /
  `suggest_loguniform` / `suggest_discrete_uniform` aliases
- `BaseGASampler.get_trial_generation`, on NSGA-II and NSGA-III
- `StudyDirection.NOT_SET` and `StudySummary`
- the `terminator` base classes and `should_terminate`

All 28 exist in Optuna 4.9.0 and all 28 were already implemented in MATLAB,
so the gap was evidence, not capability. Five new contracts cover them, and
MATLAB reproduces each: the deprecated aliases return their modern
equivalents' values from the same seeded stream, the generation sequence is
the same step function of `population_size`, and `should_terminate` flips on
the same trial.

Removing a single contract registration flips `full_compatibility_complete`
to false and fails the ledger test — which is what the old rule could not do.

## Global Optimization Toolbox surface

Optuna is configured in Optuna's vocabulary and the Global Optimization
Toolbox in MATLAB's; a Simulink user coming from `ga` or `sdo.optimize` had
to learn the second one to reach a study. `radia.optuna.optimize` adds the
toolbox idiom on top of the existing study.

```matlab
p = radia.optuna.getParameterFromModel("plant", ["Kp" "Ki"]);
p(1).Minimum = 0;    p(1).Maximum = 10;
p(2).Minimum = 1e-3; p(2).Maximum = 1;  p(2).Transform = "log";

opts = radia.optuna.optimoptions(MaxTrials=200, Sampler="tpe", ...
    Display="iter", OutputFcn=@myOutputFcn);

[x, fval, exitflag, output] = radia.optuna.optimize(@cost, p, opts);
```

| Toolbox concept | Here |
|---|---|
| `optimoptions(...)` | `radia.optuna.optimoptions` — bare Name=Value, a leading `"optuna"`, or an options object to copy |
| `MaxTime` / `FunctionTolerance` / `Display` / `OutputFcn` / `PlotFcn` / `UseParallel` | same names, same meanings |
| `MaxGenerations` / `MaxStallGenerations` | `MaxTrials` / `MaxStallTrials` — a trial is the unit the study records |
| `[x, fval, exitflag, output]` | same, with ga's exitflag codes |
| `output.funccount` / `.message` / `.rngstate` | same, plus `study`, `bestTrialNumber`, `prunedcount`, `failedcount` |
| `stop = fcn(x, optimValues, state)` | same, with `"init"` / `"iter"` / `"done"` |
| `sdo.getParameterFromModel` | `radia.optuna.getParameterFromModel` |
| — | `Sampler`, `Pruner`, `Seed`, `StoragePath`, `Directions` keep their Optuna names |

Exit flags follow `ga`: `1` stall, `0` MaxTrials, `-1` stopped by a hook,
`-2` no feasible point, `-5` MaxTime.

Two deliberate departures from Optuna's own shape:

- **The objective takes values, not a `Trial`.** That is the toolbox's
  contract, and it is what makes `UseParallel` honest: suggesting mutates the
  study and has to stay serial, so only a values objective can be evaluated
  on workers without the sampler's sequence depending on scheduling.
  `UseParallel` asks for a batch, evaluates it on the pool, then tells the
  results; with no pool open it errors rather than quietly running serially,
  which would make a timing comparison meaningless.
- **Bounds are never invented.** `getParameterFromModel` returns the model's
  values with `Minimum`/`Maximum` left infinite, exactly as `sdo` does, and
  `optimize` refuses to search a parameter whose range the user has not set.

The short-name-to-sampler mapping moved out of `optunaSFunction`, where
nothing outside Simulink could reach it and the seed was pinned to 0, into
`radia.optuna.internal.samplerFromName`. Both callers share it, and it
refuses `"auto"` rather than guessing without the study's shape.

This surface is a **MATLAB extension**. Upstream Optuna has no
`optimoptions` and no exitflag, so it is classified `matlab-integration` and
must never be cited as evidence of Optuna compatibility.

## History storage

The initial hypothesis — that per-trial scanning dominated — was wrong, and
the measurement said so. Profiling `get_trials()` over 2000 trials on hibino
showed `datetime` conversion at 44% of the time and history lookup nowhere in
the top fourteen entries.

What actually helped was reducing call counts, not table management:

1. `freezeTrials` converts each serial timestamp column once and slices it.
   `datetime()`'s per-call overhead is roughly 250x its per-element cost
   (6000 scalar calls 0.749 s versus one vector call 0.0030 s).
2. Trials that never reported reuse a shared empty table template
   (2002 table constructions down to 667).
3. `TrialState.toStorage` short-circuits an already-canonical string.

Result: **`get_trials()` at 2000 trials fell from 2.600 s to 1.000 s**, and
`datetime` left the profile entirely.

A rejected alternative worth recording: replacing
`datetime(s,'ConvertFrom','datenum','TimeZone','local')` with
`epoch + days(s)` is **not** equivalent — measured divergence up to 1.1e6 ms.
Constructing unzoned and then assigning `.TimeZone` is equivalent.

`IntermediateTable` also became a column store behind a materialized view,
and every store gained a trial-number bucket index
(`radia.optuna.internal.TrialRowIndex`). The index is a cache, never an
authority: each lookup receives the store's own key column and rebuilds when
the two disagree, so a missed notification costs speed, not correctness.

### Measured scaling (hibino, idle)

| trials | write (s) | freeze (s) | scan/indexed lookup |
|---|---|---|---|
| 500 | 1.35 | 0.30 | 1.6x |
| 1000 | 1.77 | 0.57 | 3.8x |
| 2000 | 3.39 | 1.01 | 4.0x |
| 4000 | 6.41 | 1.94 | 5.5x |
| 8000 | 13.56 | 3.88 | 10.8x |

Fitted exponents: write 0.85, freeze 0.92, indexed lookup **-0.18** (flat),
scan lookup 0.43.

The index only shows up when the probe count is held fixed while the store
grows. Measured through `get_trials()` as a whole it is invisible, because
per-trial construction dominates. Reproduce with
`validation_test/optimization/benchmark_optuna_history_store.m` and its
committed results JSON.

## Storage bridge

MATLAB's table/MAT storage is not an `optuna.storages` backend, so "if MATLAB
is missing a feature, just run Python Optuna" was not true — the study could
not be opened from Python at all. That gap is closed by an explicit handoff,
which is also what makes scoping `optuna.storages` as `bridged` honest rather
than an excuse.

```
MATLAB Study --export_study--> study-export.v1 --radia_optuna.bridge--> Optuna storage
      ^                                                                        |
      +----------- import_study <-- bridge.from_study <-------------------------+
```

- Names travel as records, not JSON object keys, because a parameter or
  attribute name need not be a valid MATLAB field name. `Study.addTrial`
  gained an `OriginalNames` map so an imported `x-1` does not come home as
  `x_1`.
- Distributions travel as Optuna's own `distribution_to_json` strings, so
  neither side reinterprets them.
- Verified end to end: MATLAB to JSON to a sqlite Optuna storage and back
  preserves states, parameters, objectives and attributes exactly.
- CLI: `radia-optuna-bridge load|dump`.

## Deliberate omissions

Recorded here because a plausible wrong answer is worse than a missing method.

- **`BaseSampler.sample_relative`** — the base returns the trial's
  already-computed relative params filtered to the requested search space,
  rather than recomputing them. (An earlier draft of this review argued the
  method should be omitted entirely; that was wrong. Upstream declares it
  abstract, and reading back what the sampler's own relative phase produced
  is a real answer, not a guess.)
- **`BaseSampler.reseed_rng`** — a documented no-op, which is exactly what
  upstream's base class does (`pass`); samplers that own an RNG override it.
  (Also corrected from the first review, which argued for omission.)
- **`gamma(n)` under a dynamic search space** — clamped to the observation
  count. Identical to upstream whenever every trial carries the parameter.

## Remaining work

1. **68 asserted entries** in the bridged and out-of-scope scopes. They are
   discharged by the storage bridge or by an explicit out-of-scope decision,
   not by a differential test, and the ledger says so. Turning any of them
   into `verified` requires a real upstream contract, not a list entry.
2. **The derived link is per-name, not per-behaviour.** It proves an oracle
   section exercises the name; it does not prove the section covers every
   behaviour of that name. Tightening it would mean attributing individual
   assertions to individual API entries.
3. **`CmaEsSampler`** restart, margin, separable, warm-start and
   learning-rate modes remain explicit gaps, recorded as limitations.

## Reproducing

```
python tests/matlab/fixtures/generate_optuna49_oracle.py
python tests/matlab/fixtures/generate_optuna49_api_coverage.py
python tests/matlab/fixtures/generate_optuna_test_manifest.py
python -m pytest packages/radia-optuna/tests -q
python -m pytest packages/radia-mcp/tests -q -k "optuna or matlab"
```

MATLAB suites: every `tests/matlab/test_optuna*.m` (137 tests). Timing work
runs on mdx or hibino only; hibino has no NAS mount in an SSH session, so
stage the tree to `C:\temp` with tar and run it there.
