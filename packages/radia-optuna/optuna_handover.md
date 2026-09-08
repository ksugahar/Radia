# radia-optuna Simulink optimization — implementation handover

Rewritten 2026-08-29 and revised 2026-09-08 for the Optuna 5 migration. This is
the canonical product and migration handover; the unrelated EQNEDT64 document
is not an Optuna authority.

## 1. Product goal

The product is a MATLAB-native Optuna implementation that lets a student learn
optimization by changing a Simulink model, defining design variables and
requirements, running trials, inspecting why they succeeded or failed, and
applying a selected trial back to the model.

The intended experience is informed by two MathWorks products:

- Global Optimization Toolbox supplies the MATLAB solver idiom: parameter
  bounds, an `optimoptions`-style options object, solver progress callbacks,
  exit flags, and a result structure.
- Simulink Design Optimization Response Optimizer supplies the interactive
  workflow: obtain parameters from a model, define goals and constraints, run
  and stop an optimization, inspect iteration and Pareto plots, save a session,
  select an outcome, and update the model.

Those products are references for the MATLAB and Simulink operating surface.
They are not runtime dependencies and do not replace Optuna's algorithms.

The algorithmic source of truth is pinned upstream `optuna==5.0.0`.
MATLAB vectorization, MEX kernels, deterministic batching, table/MAT
persistence, Simulink signals, and session tooling may improve performance and
teaching value, but they must not silently change a shared Optuna algorithm.

## 2. Scope and non-goals

### In scope

- A toolbox-shaped MATLAB entry point over the existing Study/Trial/Sampler
  implementation.
- A persistent optimization session with explicit lifecycle state.
- A generic masked Optuna block in the single Radia Simulink library.
- Model design variables, scalar or multiobjective values, constraints,
  intermediate metrics, pruning, failure classification, and telemetry.
- Seeded reproducibility and explicit parallel scheduling semantics.
- Trial history and Pareto inspection while a run is in progress.
- Pause, resume, cancel, checkpoint, restore, and apply-selected-trial actions.
- Short regression tests under `tests` and long performance/scaling checks
  under `validation_test`.
- Standalone `radia-optuna` packaging without Radia, NGSolve, oneMKL, or Cubit.
- Radia-specific MCP tools only for MATLAB/Simulink behavior that the official
  Optuna MCP server does not own.

### Out of scope

- Reimplementing a Simulink solver inside the block.
- Calling Python once per simulation step.
- Treating Global Optimization Toolbox or Simulink Design Optimization as a
  required dependency.
- Claiming MATLAB-only table storage, parallel scheduling, or Simulink
  telemetry as upstream Optuna behavior.
- Hiding a missing or incompatible `optuna_mex` behind a fallback.
- Automatically inventing bounds, objective directions, constraints, or
  parameter scales from a current model value.
- Equating API-name presence with verified upstream compatibility.

## 3. Authority hierarchy

When references disagree, use this order:

1. Upstream Optuna 5.0.0 is the behavioral oracle for shared Study, Trial,
   sampler, pruner, distribution, storage-state, and seeded random behavior.
2. Global Optimization Toolbox and Simulink Design Optimization define familiar
   MATLAB workflow conventions.
3. This document defines MATLAB-only session, signal, storage, and teaching
   behavior.
4. Existing MATLAB output is never the compatibility oracle.

Shared algorithm tests must derive expectations by executing pinned upstream
Optuna. MATLAB-only behavior must be marked `matlab-integration` and must not be
presented as evidence of upstream parity.

## 4. Migration baseline and current state

The migration started from a completed Optuna 4.9 baseline:

- `Study`, `Trial`, samplers, pruners, distributions, storage, visualization,
  integrations, and the required standalone `optuna_mex`.
- Seeded Random/TPE/evolutionary behavior, advanced CMA-ES modes, concurrent
  RUNNING constant-liar TPE, and native unscrambled Sobol through 21,201
  dimensions.
- `SimulinkRunner` with `SimulationInput`, fast restart, serial execution,
  deterministic parallel batches, schedule-dependent steady-state execution,
  failure classification, and model/provenance hashing.
- A Level-2 MATLAB S-Function that runs one trial per sample, persists normalized
  tables, and emits numerical telemetry.
- An 816-entry generated upstream public-surface inventory.

That baseline was migration input, not the target compatibility claim. The
active branch now contains only the 5.0 fixture names and bridge pin, uses
unified public TPE, named constraints, Optuna 5 defaults and mutation points,
and has deleted removed integrations and APIs. The local 150-test fast suite,
wheel verification, installed-wheel Simulink E2E, and paired LAB development
benchmark pass. Release status remains pending until the fresh long mdx
performance lane, CI, merge, tag, publication, and release-quad gates pass.

The follow-up scalar performance change enables the completed-history MEX for
default sequential TPE, reuses unchanged distributions and NaN-bearing metadata,
and invalidates native history on reseed. The paired prewarmed LAB measurement
in `validation_test/optimization/results_optuna50_paired_lab_20260908.json`
reports MATLAB/Python throughput ratios of 1.222 scalar and 2.184 grouped TPE.
The separate table-export ratio is 0.886; no universal speed claim is made.
All 150 MATLAB tests (74 shared-oracle and 76 integration) and 11 Python package
tests pass after this change; the rebuilt wheel passes strict source fidelity.
Engine startup/calculation/
shutdown passed on both mdx CI runner accounts in run 34210024491. mdx2 has
Engine 26.1 installed, but its SSH startup still timed out; this does not
provide a MATLAB performance result on mdx2.

The prior coverage audit also remains binding: verified entries must be derived
from the test manifest and oracle fixture provenance. A maintained allow-list
may describe a mapping, but it cannot turn that mapping into verification.

Therefore the public claim is:

> The required declared Optuna 5.0.0 compatibility scope is oracle-covered only
> after the 5.0 inventory, direct oracle fixture, real-transport MCP fixture,
> MATLAB fast suite, and evidence ledger all agree. Historical 4.9 evidence does
> not satisfy this claim.

Do not restore the older blanket “816 verified” wording. The generated 5.0
ledger is 812 present, 749 oracle-verified, 63 explicitly asserted, and all
401 required entries oracle-mapped.

## 5. Architecture

`OptimizationParameter` and `OptimizeOptions` form the user contract.
`OptimizationSession` owns lifecycle, history, persistence, and selection.
`Study` owns Optuna state. `SimulinkRunner` evaluates model trials.
`optuna_mex` owns performance-critical native kernels. The Level-2 MATLAB
S-Function is a readable Simulink adapter, not the optimizer implementation.

```text
masked Radia Optuna block
        |
        v
OptimizationSession ---- checkpoint/session metadata
        |
        +---- Study / Trial / Sampler / Pruner
        |             |
        |             +---- optuna_mex
        |
        +---- SimulinkRunner ---- SimulationInput / sim / parsim
        |
        +---- normalized trial tables + telemetry + plots
```

There is one canonical short-name sampler factory. Both the toolbox-shaped API
and the Simulink path call it. The block must never carry a separate sampler
mapping or a hard-coded seed.

## 6. Student workflow

A complete exercise follows this visible sequence:

1. Open or build a model and run a baseline simulation.
2. Select scalar tunable parameters from model or base workspace.
3. Set initial value, lower/upper bounds, type, scale/transform, and optional
   allowed categorical values.
4. Select objective and constraint signals or provide a tested objective
   function.
5. Select sampler, pruner, explicit seed, trial budget, time budget, and
   execution mode.
6. Start the session from an explicit trigger.
7. Observe current trial, best value, completed/pruned/failed counts, elapsed
   time, latest failure code, and Pareto revision.
8. Pause or cancel without corrupting completed history.
9. Save and restore the session.
10. Select any completed or Pareto trial, inspect its parameters and metrics,
    apply it to the model, and simulate again.
11. Compare two saved sessions that differ in search space, sampler, seed, or
    pruning policy.

The interface must expose the experiment, not only the winning number. A
student must be able to answer:

- Which parameter values were proposed?
- Which objective or constraint caused rejection?
- Was the trial complete, pruned, failed, or cancelled?
- Which sampler and seed generated it?
- What changed when another optimization policy was selected?
- Can the selected trial be reproduced from the saved session?

## 7. MATLAB public API

### Parameters

`OptimizationParameter` represents one scalar design variable:

- `Name` and current `Value`
- `Minimum` and `Maximum`
- `Free`
- `Type`: continuous, integer, or categorical
- `Transform`: linear or log
- `Step` for quantized numeric variables
- `Choices` for categorical variables
- optional model/workspace binding metadata

`getParameterFromModel(model,names)` reads current values but leaves bounds
unbounded. Optimization refuses a free unbounded parameter rather than
inventing a range.

### Options

`radia.optuna.optimoptions` accepts bare name/value arguments, a leading
`"optuna"` solver name, or an existing options object to copy and modify.

The initial contract includes:

- `MaxTrials`, `MaxTime`, `MaxStallTrials`, `FunctionTolerance`
- `Sampler`, `Pruner`, `Seed`
- `Direction` or `Directions`
- `Display`, `OutputFcn`, `PlotFcn`
- `UseParallel`, `ParallelMode`, `BatchSize`
- `StoragePath`, `StudyName`, `Resume`
- `CatchObjectiveErrors`

Unknown or invalid options fail loudly.

### Solver contract

```matlab
[x,fval,exitflag,output] = radia.optuna.optimize(fun,parameters,options)
```

The returned `x` contains optimized and fixed parameters. `output` contains the
Study, trial counts, elapsed time, best/Pareto information, stop reason, and
provenance. Output and plot callbacks use `"init"`, `"iter"`, and `"done"`
states.

This is a MATLAB integration API. It does not pretend upstream Optuna has
`optimoptions` or MATLAB exit flags.

## 8. OptimizationSession contract

`OptimizationSession` is the unit that both MATLAB scripts and Simulink operate.
It has explicit states:

```text
configured -> running <-> paused -> completed
                  |          |
                  +----------+-> cancelled
                  +------------> failed
```

Required operations:

- `configure` or constructor validation
- `start`
- `runNext` for one trial
- `run` for a bounded batch
- `pause` and `resume`
- `cancel`
- `save` and static `load`
- `selectTrial`
- `selectedParameters`
- `applySelectedToModel`
- `snapshot` for telemetry and plots

Every state transition is persisted when a storage path is configured. A
process interruption may lose the currently RUNNING trial but must not lose
completed history. On restore, an orphaned RUNNING trial is converted to a
documented failed/stale state before new work begins.

## 9. Objective, constraints, intermediate values, and pruning

The objective boundary supports two forms:

- Upstream-shaped callback: `fun(trial)`, which may call `suggest_*`,
  `report`, and `should_prune`.
- Model-shaped callback: `fun(values,context)`, where `values` is the parameter
  struct and `context` can evaluate a Simulink model, record intermediate
  metrics, and return objective/constraint values.

A model result uses a documented structure:

```matlab
result.Objectives
result.InequalityConstraints   % <= 0 is feasible
result.EqualityConstraints
result.IntermediateValues
result.UserData
```

Signal-sink blocks may publish objective, inequality, equality, and intermediate
metric values into the session runtime store. The optimizer translates them
into the same Study/Trial contract. Pruning occurs only at declared observation
steps and records the last intermediate value.

## 10. Simulink block v2

The production interface is deliberately split. `Optuna Study` is a masked
two-input/five-output subsystem for students and ordinary models: four scalar
convenience signals plus one fixed-schema `OptunaMonitorBusV1`. It contains
the generic Level-2 MATLAB S-Function as an advanced internal runtime, so the
algorithm and checkpoint behavior are unchanged while top-level wiring stays
small. Variable-length history and Pareto data remain in normalized MAT tables,
so trial-budget changes never alter Bus width. The standalone advanced builder retains the stable six-input/eighteen-
output ABI for models that truly consume every command and telemetry signal.

### Mask parameters

- objective source/function
- parameter specification or session/config file
- trial count and optional time limit
- direction(s)
- sampler and pruner
- explicit seed; an empty seed means fresh private entropy
- storage/checkpoint path and resume policy
- trial sample time
- execution mode and batch size
- live monitor choice
- failure policy

### Compact public interface

- inputs: start and cancel;
- outputs: best objective, status, attempted-trial progress, and best trial;
- mask actions: review saved study and apply a selected trial (`-1` = best);
- mask experiment controls: objective, sampler, seed, pruner, ordered search
  space, total budget, storage, direction, sample time, and target model.

A new sampler/seed/search-space comparison receives a new MAT study path. A
budget extension keeps the same path and resumes the existing history. Both
workflows leave model wiring unchanged.

### Advanced input commands

The advanced runtime retains start, cancel, pause, resume, selected trial, and
apply. Commands are edge-triggered and idempotent.

### Output telemetry

Keep the existing numeric outputs for backward compatibility:

- best value and best trial number
- session status
- completed trial count
- latest value and elapsed time
- best-update pulse
- Pareto count, first two objective coordinates, and revision
- failed and attempted counts
- latest failure code

Add selected-trial number, pruned count, current trial, and checkpoint revision
only through a versioned, documented extension. Do not silently reorder old
ports.

### Monitor

Scopes and XY Graphs remain valid lightweight sinks. A richer MATLAB monitor may
show:

- best and current objective history
- complete/pruned/failed state timeline
- parameter trajectories
- constraint violation
- Pareto scatter with trial selection
- session and seed metadata

Monitoring is read-only. Applying a trial is an explicit operation.

## 11. Reproducibility and parallel execution

An explicit equal seed, equal sampler options, equal ordered search space, equal
trial history, and equal objective values must consume randomness in the same
order as upstream Optuna for shared sequential behavior.

`seed=[]` or `None` is nondeterministic constructor behavior. Each sampler draws
private entropy and must not mutate MATLAB's global random stream.

For a persisted Simulink study, the resolved private seed and sampler state are
part of the checkpoint: a budget-only extension restores both even when the
mask seed remains empty. A changed sampler, seed, pruner, or direction must use
a new MAT study path and fails loudly against incompatible saved metadata.

Parallel modes are different products and must be named:

- `sequential`: upstream proposal-sequence parity lane.
- `parallel-batch`: a reproducible batch of suggestions is frozen before model
  evaluation. Results are committed in trial-number order.
- `parallel-steady-state`: workers are refilled as they complete. This can be
  faster but proposal order is schedule-dependent.

No benchmark or documentation may call steady-state execution reproducible from
seed alone.

## 12. Storage and artifacts

Normalized Study tables remain the durable MATLAB source of truth:

- trials
- parameters
- objective values
- intermediate values
- user attributes
- system attributes
- constraints

Session metadata adds:

- schema version
- model path and SHA
- sampler/pruner configuration and seed
- parameter definitions and ordering
- objective/constraint definitions
- execution mode and runtime versions
- selected trial and UI state
- checkpoint revision and stop reason

The explicit `study-export.v2` bridge to an upstream Optuna storage remains a
batch handoff, never a runtime fallback.

## 13. Performance contract

Measure end-to-end performance on identical models and trials:

- cold MATLAB/Python startup separately
- warmed steady-state median over repeated runs
- suggestion time
- model evaluation time
- history freeze/query time
- persistence time
- parallel throughput and worker utilization
- data-transfer cost

The MATLAB version should meet or exceed upstream throughput where native MEX,
vectorization, and Simulink integration provide an advantage. It must never
trade away seeded algorithm parity on the sequential lane merely to win a
benchmark.

The automatic-multivariate scalar TPE path maintains an incremental
intersection-search-space cache keyed by the ordered finished-trial prefix.
Appending a trial filters the current intersection once; a changed or restored
history invalidates the prefix and rebuilds from source data. This removes the
former quadratic history scan while preserving the upstream proposal sequence.

Long scaling benchmarks belong in `validation_test/optimization`. Fast tests
may assert results and conservative non-regression bounds; they must not encode
machine-specific timing as correctness.

## 14. Test and oracle plan

### Fast `tests`

- Every shared behavior reads an upstream-generated fixture.
- Every MATLAB test function appears in the checked oracle manifest.
- Toolbox-shaped API tests are explicitly `matlab-integration`.
- Session lifecycle tests cover pause/resume/cancel/save/load/select/apply.
- Simulink tests cover clean open, initialization, command edges, telemetry,
  persistence, typed failure, repeat runs, and teardown.
- Installed-wheel tests run without the repository or Radia on the MATLAB path.
- The coverage ledger is regenerated from evidence and is byte-stable.

### Long `validation_test`

- high-dimensional Sobol and scaling
- large trial-history query/freeze scaling
- multiworker batch and steady-state throughput
- long Simulink fast-restart and lifecycle stability
- performance comparison with pinned upstream Optuna
- large multiobjective/Pareto sessions

Every long run writes machine, MATLAB, Python, Optuna, NumPy, SciPy, PyTorch,
`cmaes`, compiler, and MEX versions with the measured result.

## 15. Simulink artifact editing rule

The tracked `matlab/radia_simulink_library.slx` is production source. Structural
changes must use MathWorks' official Simulink Agentic Toolkit:

1. close scratch and orphaned models;
2. `model_read` the exact tracked file;
3. inspect the declared customer-library block knowledge;
4. `model_edit` the model;
5. `model_check` it;
6. save, close, and reopen the exact tracked path;
7. verify `FileName` and visually inspect the complete application window.

Do not patch SLX ZIP/XML. A MATLAB builder is useful for temporary reconstruction
and regression tests, but raw builder output is not promoted directly to the
tracked production library.

## 16. Packaging, MCP, and licensing

`radia-optuna` remains independently installable. Generic MATLAB Optuna,
`SimulinkRunner`, the generic block, monitor, and `optuna_mex` belong in that
wheel. Radia-specific electromagnetic adapters remain declared adapters and may
require the main Radia package.

The official `optuna/optuna-mcp` server owns shared public Study/Trial operations.
`radia-mcp.matlab` owns:

- MATLAB installation and MEX health
- upstream-oracle fixture/manifest plans
- installed-wheel Simulink checks
- session and performance evidence
- Radia-only adapter validation

Keep the independent/unofficial notice, Optuna/SciPy/Joe--Kuo notices, and pinned
upstream handoff version. Do not use the Optuna logo or imply endorsement.

As of 2026-09-08 the latest released `optuna-mcp` is 0.2.0. Upstream `main`
declares 0.3.0.dev, but its server tool implementation is unchanged from 0.2.0;
the observed changes are packaging, container publication, CI hardening, Python
3.13 coverage, action pinning, attestations, and development dependencies. Use
released 0.2.0 with `optuna==5.0.0` for the real-transport fixture. Record the
server version and transport in the fixture. Do not label an unreleased
0.3.0.dev checkout as the supported server.

## 17. Optuna 5 migration and deletion plan

This is a replacement migration, not a side-by-side compatibility feature.
There will be one MATLAB implementation and one active upstream oracle.

### Stage A — freeze and inventory

- Preserve the last green Optuna 4.9 commit/tag as history; do not create a
  `V49` namespace, compatibility mode, or runtime version switch.
- Inventory 5.0 public symbols/signatures and classify added, removed, and
  behavior-changed entries.
- Mark the working branch and package documentation as migration in progress.

Exit: the 5.0 API delta and every active 4.9 reference are mechanically listed.

### Stage B — establish the 5.0 oracle before changing MATLAB behavior

- Generate direct fixtures in an isolated environment pinned to
  `optuna==5.0.0`, recording Python, NumPy, SciPy, PyTorch, and `cmaes`.
- Generate the MCP fixture through a real stdio session using released
  `optuna-mcp==0.2.0` running against Optuna 5.0.0.
- Run each generator twice and require byte-identical JSON.
- Add explicit coverage for TPE defaults and bandwidth behavior,
  multi-objective default TPE, named constraints and duplicate/NaN behavior,
  new NSGA-II mutation types, and removed APIs.

Exit: the complete 5.0 fixtures are deterministic and regeneration fails with
any other Optuna version.

### Stage C — replace behavior in place

- Change the existing MATLAB classes and MEX kernels directly. Do not retain
  parallel 4.9 implementations.
- Make TPE multivariate and constant-liar behavior follow 5.0 defaults and port
  the 5.0 bandwidth rule from the upstream algorithm.
- Make TPE the default for both single- and multi-objective studies.
- Replace positional constraint storage at the public boundary with named
  constraints while keeping an efficient internal numeric view for samplers.
- Add `BaseMutation` and `PolynomialMutation` and connect them to NSGA-II/III.
- Remove `categorical_distance_func`, public Study/Trial system-attribute APIs,
  and the removed StudySummary field instead of emulating them.

Exit: focused MATLAB tests pass against the 5.0 fixture and no production path
selects old behavior.

### Stage D — atomic oracle cutover and 4.9 deletion

The first commit that declares 5.0 compatibility must do all of the following
together:

- switch every active test, package verifier, manifest, bridge, README, policy,
  and CI command from `optuna49`/4.9.0 to `optuna50`/5.0.0;
- delete the 4.9 oracle JSON, public-API inventory, coverage ledger, MCP fixture,
  and their 4.9 generator scripts;
- delete removed API code and constructor options, not merely stop testing them;
- regenerate the distribution manifest and installed-wheel expectations;
- fail a repository scan if active source or tests still mention `optuna49`,
  `optuna==4.9.0`, or the removed compatibility members.

Historical CHANGELOG entries, tagged release evidence, and named archived
performance result files may retain 4.9 text. They are immutable provenance and
must be clearly described as historical; they are never loaded by active tests.

Deletion happens here—not before Stage B, because that would remove the last
working oracle before its replacement is proven, and not after release, because
shipping both active oracle generations would leave ambiguous truth and dead
code.

Exit: exactly one active oracle generation exists and the fast suite passes
from a clean checkout and an isolated installed wheel.

### Stage E — validation and release

- Run long seeded parity and performance/scaling jobs under `validation_test`.
- Compare warmed throughput without weakening sequential seeded parity.
- Review licenses and third-party notices, merge, tag, publish the four
  distributions, and run release-quad only from the same verified commit.

Exit: release evidence names Optuna 5.0.0 and optuna-mcp 0.2.0 explicitly, and
no release gate consumes an Optuna 4.9 active artifact.

## 18. Historical implementation record

Historical Optuna 4.9 release-candidate record (2026-08-29; not valid 5.0
release evidence):

- fast MATLAB Optuna suite: 150/150 passed, 0 failed, 0 incomplete;
- oracle ledger: 816/816 present, 748 executable-evidence verified,
  68 explicitly asserted, 0 partial/unmapped/missing;
- required compatibility scope: 400/400 executable-evidence mapped,
  0 asserted or unmapped;
- focused package/radia-mcp Python tests: 28/28 passed;
- isolated installed-wheel verification: 226 MATLAB files, 21 MEX commands,
  15 Simulink entries, seed-less four-to-six-trial block continuation,
  best-trial model application with unchanged topology, four-trial session
  checkpoint/resume, twelve-trial teaching-model run, and seven-table restore
  all passed;
- tracked production library and teaching model: official-agent
  read/edit/check/save/reopen, clean-path reopen, full-window visual QA, and
  embedded-text scans passed.
- mdx long validation: MATLAB/Python warmed-time ratios were 0.670 scalar,
  0.491 grouped-conditional, and 0.630 table export; four-worker deterministic
  batch evaluation was 2.357x sequential throughput; 4,000-trial indexed
  lookup was 5.683x the scan reference with freeze exponent 0.862; first MEX
  call median was 11.44 ms; the radia-mcp release gate returned `ready`.

### Phase 0 — specification and truthful evidence

- [x] Rewrite this handover before implementation.
- [x] Replace allow-list “verified” assertions with a derived evidence ledger.
- [x] Add the new MATLAB-integration tests to the oracle manifest.
- [x] Keep required upstream scope closure distinct from wider API presence.

### Phase 1 — toolbox-shaped MATLAB API

- [x] `OptimizationParameter`
- [x] `OptimizeOptions` and `optimoptions`
- [x] shared `samplerFromName`
- [x] `optimize` with output callbacks and exit contract
- [x] `getParameterFromModel`

### Phase 2 — persistent session

- [x] `OptimizationSession` lifecycle and snapshots
- [x] checkpoint/restore and orphan RUNNING handling
- [x] trial selection and application to model
- [x] model-shaped objective/constraint context

### Phase 3 — Simulink block v2

- [x] use the shared sampler/pruner/seed configuration
- [x] bind the block to `OptimizationSession`
- [x] add pause/resume/select/apply without breaking existing ports
- [x] place the full runtime behind a two-input/five-output student facade
- [x] expose fixed scalar telemetry through `OptunaMonitorBusV1` without
  putting variable-length history or Pareto arrays on the Bus
- [x] prove compare, budget-extension, review, and apply without rewiring
- [x] update the single Radia library through the official toolkit
- [x] clean-open, check, save, reopen, and visual acceptance

### Phase 4 — teaching and validation

- [x] a fast, toolbox-free teaching model with a known optimum
- [x] a multiobjective model with a selectable Pareto trial
- [x] pruning and failed-trial exercise
- [x] student worksheet/notebook with saved results
- [x] long mdx performance and parallel validation

### Phase 5 — distribution and release

- [x] installed-wheel API/session/Simulink proof
- [x] radia-mcp MATLAB difference-gate update
- [ ] version, CI, review, merge, tag, PyPI, and release-quad

## 19. Definition of done

The implementation is complete only when:

- a student can configure, run, compare, extend, inspect, select, and apply an
  optimization from the compact Simulink block without changing signal lines;
- wired pause/resume and full telemetry remain available only through the
  explicit advanced interface;
- equal seeded sequential shared behavior matches pinned Optuna fixtures;
- MATLAB-only extensions are separately classified;
- the evidence ledger contains no assertion disguised as verification;
- the production library passes official clean-open/edit/check/save/reopen
  handling and visual QA;
- fast tests and long validation tests are correctly separated;
- installed-wheel tests pass with no repository or Radia on the path;
- performance evidence reports cold and warm results on identical workloads;
- documentation states precisely which compatibility scope is proven.

The implementation may be better suited to MATLAB and Simulink than upstream
Optuna without becoming a different optimizer. That separation—shared
algorithm, MATLAB-native execution, explicit teaching workflow—is the design.
