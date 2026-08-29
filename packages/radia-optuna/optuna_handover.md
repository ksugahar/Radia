# radia-optuna Simulink optimization — implementation handover

Rewritten 2026-08-29. No tracked `optuna_handover.md` existed on the current
main branch, so this is the new canonical handover rather than a continuation
of the unrelated EQNEDT64 document.

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

The algorithmic source of truth remains pinned upstream `optuna==4.9.0`.
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

1. Upstream Optuna 4.9.0 is the behavioral oracle for shared Study, Trial,
   sampler, pruner, distribution, storage-state, and seeded random behavior.
2. Global Optimization Toolbox and Simulink Design Optimization define familiar
   MATLAB workflow conventions.
3. This document defines MATLAB-only session, signal, storage, and teaching
   behavior.
4. Existing MATLAB output is never the compatibility oracle.

Shared algorithm tests must derive expectations by executing pinned upstream
Optuna. MATLAB-only behavior must be marked `matlab-integration` and must not be
presented as evidence of upstream parity.

## 4. Current baseline and required correction

The current main branch already has:

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
- A 816-entry generated upstream public-surface inventory.

One audit correction must land before expanding the API. The current coverage
generator calls some entries verified because they appear in a maintained
allow-list. Verification must instead be derived from the test manifest and
oracle fixture provenance. The reconciled ledger distinguishes 748
evidence-derived entries from 68 asserted mappings, while the required shared
scope remains 400/400 directly mapped with no asserted required entry.

Therefore the public claim is:

> The required declared Optuna 4.9.0 compatibility scope is oracle-covered.
> The wider MATLAB surface is present and mapped, but the coverage ledger must
> distinguish evidence from assertion until every wider entry has direct
> evidence.

Do not restore the older blanket “816 verified” wording unless the generator
can derive that number from actual oracle tests.

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

The production block remains a generic Level-2 MATLAB S-Function in
`matlab/radia_simulink_library.slx`. It delegates numerical work to the tested
MATLAB API and standalone MEX.

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

### Input commands

The existing start and cancel inputs remain supported. The v2 command contract
adds pause/resume and apply-selected-trial through either typed command inputs
or explicit mask callbacks. Commands must be edge-triggered and idempotent.

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

The explicit `study-export.v1` bridge to an upstream Optuna storage remains a
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

## 17. Implementation order

### Phase 0 — specification and truthful evidence

- [x] Rewrite this handover before implementation.
- [ ] Replace allow-list “verified” assertions with a derived evidence ledger.
- [ ] Add the new MATLAB-integration tests to the oracle manifest.
- [ ] Keep required upstream scope closure distinct from wider API presence.

### Phase 1 — toolbox-shaped MATLAB API

- [ ] `OptimizationParameter`
- [ ] `OptimizeOptions` and `optimoptions`
- [ ] shared `samplerFromName`
- [ ] `optimize` with output callbacks and exit contract
- [ ] `getParameterFromModel`

### Phase 2 — persistent session

- [ ] `OptimizationSession` lifecycle and snapshots
- [ ] checkpoint/restore and orphan RUNNING handling
- [ ] trial selection and application to model
- [ ] model-shaped objective/constraint context

### Phase 3 — Simulink block v2

- [ ] use the shared sampler/pruner/seed configuration
- [ ] bind the block to `OptimizationSession`
- [ ] add pause/resume/select/apply without breaking existing ports
- [ ] update the single Radia library through the official toolkit
- [ ] clean-open, check, save, reopen, and visual acceptance

### Phase 4 — teaching and validation

- [ ] a fast, toolbox-free teaching model with a known optimum
- [ ] a multiobjective model with a selectable Pareto trial
- [ ] pruning and failed-trial exercise
- [ ] student worksheet/notebook with saved results
- [ ] long mdx performance and parallel validation

### Phase 5 — distribution and release

- [ ] installed-wheel API/session/Simulink proof
- [ ] radia-mcp MATLAB difference-gate update
- [ ] version, CI, review, merge, tag, PyPI, and release-quad

## 18. Definition of done

The implementation is complete only when:

- a student can configure, run, pause, resume, cancel, inspect, save, restore,
  select, and apply an optimization from Simulink;
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
