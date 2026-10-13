# radia-optuna

The native gateway is built with MATLAB `mexcuda`. TPE Parzen density,
candidate scoring and first-maximum selection can run on CPU or GPU in double
precision. Select the execution explicitly in each MATLAB process:

```matlab
radia.optuna.executionMode("cpu"); % OpenMP host kernels, same mexcuda binary
radia.optuna.executionMode("gpu"); % CUDA kernels; the initial/default mode
mode = radia.optuna.executionMode();
```

CPU/GPU selection is a development addition after 0.3.0; the published 0.3.0
wheel retains its original GPU-only TPE execution. See the
[dual-route validation record](../../validation_test/radia_optuna/results/execution-choice-20261013.json).

CPU mode does not initialize CUDA, allocate device memory or transfer data.
GPU mode requires Parallel Computing Toolbox and a supported NVIDIA GPU;
device errors raise and never switch to CPU. This is a process-wide execution
setting, independent of Study storage and RNG state. Set it before a run and
in each parallel worker; reloading `optuna_mex` resets it to GPU. Simulink
uses the setting in its MATLAB process. No automatic size-based dispatch is
provided. `api.info` reports `execution_mode`, CPU call counters and CUDA
launch/transfer counters so the chosen path can be checked.

Study/Trial management, estimator construction, seeded random consumption
and candidate generation remain host work in both modes. CPU host parallel
kernels use OpenMP in the mexcuda build; ordered random consumption remains
serial. Ordinary `mex` compilation remains prohibited. GPU users do not need
Python; optional upstream-backed features retain their requirements below.

Study history remains fully in host memory. Public tables are generated on
demand from column stores; TPE does not materialize tables on every proposal.
Independent sampling uses the original column scan for histories below 512
trials. Larger histories index parameter rows by their exact names, append
new rows to that index, and preserve chronological observation order.
Single-parameter histories use the existing columns directly. For a single
objective with at least 512 completed trials and multiple parameters,
without pruning, constraints or a custom Gamma function, TPE also
reuses the global good-trial split between parameter requests until the
trial count or states change. Other splits retain their existing computation.
These caches do not discard trials, change RNG consumption or select a
different CPU/GPU backend. They are rebuilt after loading saved history.

An empty `StoragePath` performs no persistence. To retain explicit checkpoints
without writing a complete MAT snapshot at each mutation, configure
`StoragePath` with `AutoSave=false`, then call `study.save()` at the required
checkpoint boundaries. Unsaved progress is lost if the process exits; the
default automatic persistence contract is unchanged. Additional RAM retains
more history, but does not remove serialization or synchronous storage costs.

`validation_test/optimization/benchmark_optuna_parameter_history.m` measures
independent TPE with all history retained and CPU or GPU execution in the mexcuda
module. Its JSON includes module initialization, first-trial and warmed
end-to-end times, MATLAB memory, and proposal arrays for differential
comparison. Optional MAT fixtures reuse setup inputs only; sampling storage
remains in memory. This benchmark does not compare Python/DB or Rustuna,
or measure the default multivariate native-history route.
See the [history-cache validation record](../../validation_test/radia_optuna/results/memory-history-20261013.json).

The [history-cache-only A100 measurements](../../validation_test/optimization/results_optuna_memory_performance_a100_20261013.json)
precede CPU OpenMP dispatch tuning. Their 1.52x GPU/CPU ratio used a serial CPU
density kernel at the default 24 candidates and is superseded by the comparison below.

CPU density now parallelizes independent candidate rows when there are at least
256 candidates or `candidates * mixture_components * dimensions >= 16384`,
provided more than one OpenMP worker is available. The team is bounded by the
candidate count. Component summation and first-maximum selection remain ordered.
Set `OMP_NUM_THREADS` before starting MATLAB; `api.info` reports
`cpu_parallel_density_calls` and `cpu_max_density_threads` to verify actual use.
CPU execution still requires only the same mexcuda binary and performs no CUDA
initialization or transfers.

[Formal OpenMP A100 measurements](../../validation_test/optimization/results_optuna_openmp_performance_a100_20261013.json)
use the same in-memory history implementation for both variants, default 24
candidates, fresh-session ABBA order, 1 and 16 configured threads, and 18 warmed
samples per variant/case. Times include `ask`, all suggestions, recording and
`tell`, with complete history retained and persistence disabled. The 16-thread
independent-TPE comparison is:

| History trials | Parameters | Previous CPU (ms) | Tuned CPU (ms) | GPU (ms) |
| ---: | ---: | ---: | ---: | ---: |
| 24 | 1 | 1.097 | 1.092 | 1.408 |
| 256 | 1 | 1.325 | 1.312 | 1.357 |
| 1000 | 1 | 2.541 | 1.452 | 1.797 |
| 5000 | 1 | 9.101 | 3.548 | 4.220 |
| 24 | 40 | 26.640 | 26.049 | 39.778 |
| 256 | 40 | 47.736 | 47.728 | 50.215 |
| 1000 | 40 | 93.983 | 46.427 | 60.641 |
| 5000 | 40 | 389.318 | 152.884 | 177.174 |

The threshold change speeds up the large independent cases by 1.75-2.57x.
Tuned CPU is faster than GPU in all eight measured independent cases; small-case
changes are close to baseline. All 17,056 paired proposal values match exactly,
and CPU/GPU candidate proposals also match. Earlier interrupted runs are excluded.
These are development-source results, not published 0.3.0 wheel timings or a
comparison with Python storage/Rustuna.

[Joint-TPE A100 measurements](../../validation_test/optimization/results_optuna_openmp_joint_a100_20261013.json)
use `Multivariate=true`, 40 parameters, the same default 24 candidates and
complete histories. CPU 1/4/8/16 and GPU 16 run in palindromic order, with a
separate two-session CPU-12 check because 24 candidate rows divide evenly into
12 workers. Each configuration/case has 18 warmed samples. The best measured
CPU configuration is selected separately for each case:

| History trials | CPU threads | Tuned CPU (ms) | GPU (ms) | GPU speedup |
| ---: | ---: | ---: | ---: | ---: |
| 24 | 8 | 11.792 | 12.000 | 0.98x |
| 1000 | 16 | 17.826 | 16.533 | 1.08x |
| 5000 | 16 | 53.613 | 43.897 | 1.22x |

CPU parallelism reduces the 5000-trial joint case from 209.264 to 53.613 ms
(3.90x versus serial CPU). GPU retains a modest 1.22x advantage in that case,
while its 1.08x advantage at 1000 trials is small and it loses at 24 trials.
All proposal arrays agree exactly across the measured CPU thread counts and GPU.
The recommendation is CPU for routine use, with GPU retained as an explicit
option for large joint histories whose end-to-end measurements justify it.
These results do not support promoting GPU as a general TPE speed advantage.
The process default remains GPU for compatibility; explicitly select
`radia.optuna.executionMode("cpu")` to use the recommended CPU route.

Those comparisons probe fixed histories. The
[closed-loop inference measurements](../../validation_test/optimization/results_optuna_resident_inference_a100_20261013.json)
instead optimize Branin and a 40-dimensional correlated quadratic valley for
800 trials, seeds 11/37, default 24 candidates and `Multivariate=true`.
Twelve fresh idle-A100 MATLAB Engine sessions run CPU / previous GPU / resident
GPU / resident GPU / previous GPU / CPU, separately for batch sizes 1 and 8.
CPU uses 16 configured threads. The following medians pool both seeds and
both sessions, trials 513-800; inference includes `ask` and suggestions,
with objective evaluation and `tell` timed separately:

| Problem | Batch | CPU inference (ms) | Previous GPU (ms) | Resident GPU (ms) | GPU / CPU speedup |
| --- | ---: | ---: | ---: | ---: | ---: |
| Branin, 2D | 1 | 1.146 | 1.441 | 1.138 | 1.01x |
| Correlated valley, 40D | 1 | 13.016 | 12.351 | 12.000 | 1.08x |
| Branin, 2D | 8 | 1.758 | 2.075 | 1.797 | 0.98x |
| Correlated valley, 40D | 8 | 23.363 | 23.632 | 23.414 | 1.00x |

Full proposal arrays, objective values and convergence curves agree exactly
across backends within each seed/problem/batch configuration. The retained
workspace improves warmed inference over the previous GPU by 1.01-1.27x,
but establishes no broad advantage over tuned CPU. Cold-inclusive totals are
also recorded: first-use GPU initialization can reverse a warmed advantage.
The default `ConstantLiar=true` includes pending observations, changing models
during batch proposals; these runs avoid allocations but save no input-transfer
bytes. GPU residency alone cannot eliminate host model construction.
In the sequential 40D case, `ask` contains multivariate model preparation and
candidate sampling through `beforeTrial`, plus Study setup. Its median is
CPU 5.492 ms versus resident GPU 4.559 ms (1.20x). The subsequent 40 suggestions,
which retrieve relative parameters and record them, take 7.516/7.435 ms.
These API phases show where the gain appears; neither is a kernel-only timer.

A separate six-session
[model-reuse comparison](../../validation_test/optimization/results_optuna_resident_reuse_a100_20261013.json)
explicitly sets `ConstantLiar=false`, with batch size 8 and the same cases,
seeds and measurements. This preserves that option's algorithm on all three
backends; its trajectories must not be compared with `ConstantLiar=true`
as if they were the same optimization configuration. Input transfers fall
by 85.8%, confirming that model data remain on the GPU between proposals.
Warmed Branin inference is CPU 1.786 ms, previous GPU 2.072 ms, resident GPU
1.783 ms. For the 40D valley the corresponding times are 23.245, 23.713 and
23.740 ms. Thus transfer reuse does not establish an overall GPU speed
advantage: workspace reuse helps small GPU calls, while the higher-dimensional
inference remains close to tuned CPU and its previous GPU implementation.
Host estimator construction and packing remain in this implementation;
isolated phase profiling is needed before attributing the remaining cost.
The GPU route remains optional, with no general acceleration claim.

GPU times include allocation, transfers, density and selection; isolated
transfer seconds are not measured. JSON includes first-trial and module/Engine
startup times, before/after MATLAB memory and native temporary memory counters.
Whole-runtime peak memory and objective-function execution time are not measured.
The speedups concern proposal processing, not the full optimization including
an expensive objective. CPU thread screening is limited to the configurations
shown and does not establish a universal optimum.

Use the same TPE from MATLAB or the included masked Simulink
study block and `SimulinkRunner`. Simulink is optional for MATLAB-only use.

Build on the designated 100 host with `Build.ps1 -OptunaMexOnly -Test`.
MATLAB selects its supported CUDA toolkit and Visual Studio compiler.
The focused acceptance saves `validation_test/optimization/results_optuna_cuda_100.json`,
including upstream 5.0.0 density/selection checks, exact random consumption,
CUDA launch counters, runtime versions, transfers and temporary device memory.
This is correctness evidence; formal speed comparisons use the idle A100 host.

`radia-optuna` is an independent, separately distributed MATLAB optimization
component from the Radia monorepo. It installs the `radia.optuna` MATLAB
namespace, the 22-command `optuna_mex`, and checked compatibility contracts
whose behavioral oracle is Optuna 5.0.0. Its checked public inventory contains
812 present entries: 108 are oracle-verified in the checked coverage ledger
and 704 remain assertion-mapped. All 401 required entries are present;
79 have recorded oracle mappings and 322 still require that evidence.
Full API compatibility verification is therefore incomplete. Focused TPE
oracle tests certify their exercised cases, rather than the entire API. It
does not install or load the Radia solver, NGSolve, oneMKL, or Cubit.

Optuna 5.0.0 is also the algorithmic source of truth: native MATLAB and MEX
paths preserve its equations, transforms, state updates, boundary handling,
defaults, extension points, and seeded random-consumption order. The active
implementation uses the 5.0 design directly; removed 4.x APIs and state shims
are not kept. MATLAB vectorization, parallel trial
scheduling, table/MAT persistence, and Simulink telemetry are performance and
workflow extensions around that common algorithm rather than alternative
default optimizers.

Install it by itself:

```powershell
python -m pip install radia-optuna
radia-optuna-path
radia-optuna-doctor --json
```

For the optional upstream Python dependencies of this MATLAB distribution:

```powershell
python -m pip install "radia-optuna[upstream]"
```

The CPU package's `radia[optuna]` extra installs upstream Python Optuna;
install `radia-optuna` separately for MATLAB/Simulink. Use
`radia-optuna[upstream]` when GP,
scrambled-QMC, importance, advanced terminators, storage transports,
visualization, or lazy integration exports need the pinned upstream Python,
SciPy, and PyTorch stack. Individual visualization and integration targets
retain the same optional third-party dependencies as upstream Optuna.

Add the printed directory to MATLAB, then use the upstream-shaped API:

```matlab
addpath("<output of radia-optuna-path>")
study = radia.optuna.create_study( ...
    sampler=radia.optuna.TPESampler(Seed=42));
study.optimize(@(trial) (trial.suggest_float("x", -2, 2) - 0.25)^2, 100);
```

For a Global Optimization Toolbox-shaped MATLAB workflow, use parameter
objects and a persistent session:

```matlab
p = radia.optuna.OptimizationParameter("Kp", ...
    Value=1, Minimum=0, Maximum=10);
options = radia.optuna.optimoptions("optuna", ...
    Sampler="tpe", Pruner="median", Seed=42, MaxTrials=50);
[x,fval,exitflag,output] = radia.optuna.optimize(@myObjective,p,options);

session = radia.optuna.OptimizationSession(@myObjective,p,options);
session.start(); session.runNext(); session.pause();
session.resume(); session.run();
```

Generic Simulink optimization is part of the standalone contract.
`radia.optuna.SimulinkRunner` configures `Simulink.SimulationInput` objects,
runs the model, extracts objectives and constraints, classifies failed trials,
and records reproducible execution metadata without loading Radia. Radia-owned
electromagnetic models and application blocks remain in the main distribution.

The wheel ships two deliberately separated Simulink interfaces. The default
`radia.simulink.buildOptunaStudyBlock` facade exposes only `start`, `cancel`,
`best`, `status`, `progress`, and `best trial`. Sampler, seed, pruner, bounds,
budget, storage, study review, and trial application live in its mask, so a
student can compare configurations, extend a saved study, and apply a result
without rewiring the model. The complete six-input/eighteen-output
`radia.simulink.buildOptunaBlock` Level-2 runtime remains available as an
advanced interface. `radia.simulink.addOptunaMonitor` remains available where
wired Scope/XY telemetry is genuinely useful.

Both block interfaces require MATLAB/MEX sampling at step time. Their `auto`
route uses TPE for small fixed-numeric budgets where batch routing selects GP.
Explicit GP and scrambled QMC (the short-name `qmc` route) are rejected, including
legacy blocks and resumed GP/QMC studies: they call Python per trial. Use those
samplers in a batch study instead. The block does not silently substitute
unscrambled QMC or the non-oracle MATLAB GP backend. User objective callbacks
remain responsible for obeying the same step-time boundary.

`radia_optuna_teaching.slx` and
`radia.simulink.buildOptunaTeachingModel` provide known-optimum, Pareto, and
pruned/failed student exercises. See
[`OPTUNA_SIMULINK_validation host.md`](OPTUNA_SIMULINK_validation host.md). These exercises do not
require Global Optimization Toolbox or Simulink Design Optimization.

The distribution is Windows x64 because the current native artifact is
`optuna_mex.mexw64`. Native Random/TPE/evolutionary/pruner workflows do not
start Python. Install `radia-optuna[upstream]` for features intentionally
executed through pinned upstream Python packages, including checked GP
acquisition, scrambled QMC, and parameter importance.

The native sampler surface includes unified scalar and multi-objective TPE,
automatic multivariate selection, concurrent-RUNNING constant liar, named
constraints, polynomial mutation, source-trial/separable/margin/learning-rate
CMA-ES modes, and deterministic unscrambled Sobol generation through 21,201
dimensions. The Sobol path uses a
checked binary conversion of SciPy 1.17.1's Joe--Kuo criterion-6 direction
numbers and does not import Python or SciPy at MATLAB runtime.

MAT/table persistence remains a MATLAB extension. An explicit, versioned
handoff is available when the same completed history must be opened by an
upstream Optuna storage:

```matlab
radia.optuna.export_study(study, Path="study.json")
```

```powershell
radia-optuna-bridge load study.json --storage sqlite:///study.db
radia-optuna-bridge dump returned.json --storage sqlite:///study.db --study-name <name>
```

`radia.optuna.import_study("returned.json")` restores trial states, original
parameter and attribute names, distributions, intermediate values,
constraints, metric names, and study attributes. This is an explicit handoff,
not a per-trial Python fallback.

`LTspiceRunner`, `SheetMetalRunner`, and `internal.runLTspiceTrial` are shipped
as explicitly classified Radia integration adapters. The standalone core does
not call them; choosing one requires the full `radia` installation.

## Upstream attribution

This is an independent, unofficial project. It is not affiliated with,
sponsored by, or endorsed by Preferred Networks, Inc. or the Optuna project.
It does not use the Optuna logo or present itself as an official Optuna
distribution.

Optuna, the Optuna logo and any related marks are trademarks of Preferred Networks, Inc.

Optuna and the official `optuna/optuna-mcp` server are separate MIT-licensed
upstream projects and are not bundled in this wheel. See
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for their copyright and
license notices and for the SciPy Sobol direction-data notice. Shared MCP
Study/Trial/visualization operations remain owned by the official server;
`radia-mcp` covers only MATLAB/Simulink differences.

## Release boundary

CI builds and verifies a distinct `radia-optuna-wheel` artifact. A
`radia-optuna-v<version>` tag may publish only that exact CI artifact, after
rechecking its tag, version, API inventory, MEX inventory, and dependency
boundary, including the bundled third-party notice. The verifier byte-matches
every staged Python, MATLAB, Simulink, MEX, contract, license, and notice
payload against the checked monorepo source; matching version and file counts
alone are not sufficient.

`radia-mcp.matlab` owns the executable MATLAB difference gate. Use
`matlab_optuna_health`, `matlab_optuna_oracle_plan`, and
`matlab_optuna_benchmark_plan`, then submit the resulting evidence to
`matlab_optuna_release_gate`. Installed-wheel evidence can be captured with:

```powershell
pwsh -File packages/radia-optuna/tests/run_installed_wheel_simulink.ps1 `
  -Wheel <wheel> -EvidenceOutput private-runtime-path
```

MATLAB runs through the MATLAB Engine for Python. The runner installs the
`matlabengine` release matching that MATLAB's version into its isolated venv,
checks that the Engine binds to that MATLAB root, and records both in the
evidence; the acceptance runs in a function workspace. With
`-EngineSession <name>` the runner reuses that shared session after a PID check
and restores its path, folder and environment without quitting it; it refuses a
session with loaded diagrams or Radia/Optuna MEX. Without it, the runner starts
one dedicated session only when the host has no MATLAB process and no shared
Engine, and fails before launch otherwise. `release_quad optuna-candidate`
accepts the same choice per target as `--engine-session HOST=NAME`.

This evidence includes source-fidelity verification, the installed doctor,
standalone Simulink success and typed-failure paths, and all seven persisted
tables after MAT reload. Shared Study/Trial MCP operations remain upstream.

Acceptance requires a GPU architecture included in the compiled artifact. Additional device acceptance remains pending until the same artifact is verified there.

### Reproducible performance evidence

The 2026-10-12 CUDA gateway packed inputs into one host-to-device transfer and
used one per-call GPU allocation, including joint TPE evaluation. It retained ordered
random consumption and the existing GPU kernels; it never chooses a CPU
fallback for a declared GPU operation.

The development gateway now retains a MATLAB-owned GPU workspace, up to 64 MiB
per device. An exact host byte shadow identifies unchanged input ranges;
only the contiguous changed range is uploaded. This reuses model data when
successive calls share it without caching a stale estimator. MATLAB GPU-data
liveness detects device resets, including reused allocation addresses.
Switching to CPU or unloading the module releases retained workspaces.
Larger calls use temporary allocations; `allocator` describes that temporary
allocation route. `api.info` additionally reports `device_allocation_calls`,
`workspace_reuses`, `avoided_upload_bytes` and `cached_device_bytes`.
Host estimator construction and ordered random consumption remain unchanged.

`benchmark_optuna_inference.m` measures a complete optimization loop on Branin
and a 40-dimensional correlated quadratic valley. It separates inference
(`ask` plus suggestions, including preparation/recording), objective evaluation,
and history updates. `BatchSize=8` proposes eight trials before updating their
observations, with serial objective evaluation. These public optimization
benchmarks do not certify engineering CAE objectives or change TPE into a GP
surrogate. CPU/GPU comparisons must preserve the full proposal and convergence
trajectories for each seed and batch size.

An idle, user-selected A100 40 GB host was used for the 2026-10-12 formal
comparison because the original timing host had no NVIDIA GPU. MATLAB R2026a
used one host thread. Four fresh MATLAB sessions ran baseline / optimized /
optimized / baseline, with 31 repetitions per shape per session. The TPE call
includes candidate sampling, input handling, transfers and GPU selection;
random-handle creation/destruction is outside the timer.

| Candidates | Mixture components | Dimensions | Baseline TPE (us) | Optimized TPE (us) | Speedup |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 24 | 32 | 1 | 373.2 | 315.3 | 1.18x |
| 24 | 32 | 5 | 407.8 | 367.4 | 1.11x |
| 256 | 128 | 1 | 468.1 | 420.3 | 1.11x |
| 256 | 128 | 5 | 625.7 | 595.5 | 1.05x |
| 4096 | 512 | 1 | 3020.7 | 2768.8 | 1.09x |
| 4096 | 512 | 5 | 5624.1 | 5450.8 | 1.03x |

Small density-only calls remain slower than the vectorized MATLAB CPU
reference. These results compare two GPU implementations and do not establish
a speedup over Optuna Python or a CPU C++ implementation. Selected candidates
match across variants; both the RTX 5000 build host and A100 pass the Optuna
5.0.0 oracle checks. See the [machine-readable results](../../validation_test/optimization/results_optuna_cuda_performance_a100_20261012.json)
for raw observations, hashes, uncertainty and limitations.

To reproduce, use the same hash-verified binaries built on the build host and
run `benchmark_optuna_cuda(binaryDirectory, outputPath, variant, 31)` from
`validation_test/optimization` in a fresh owned MATLAB session for each variant.
Set `OMP_NUM_THREADS=1` before starting MATLAB. Check idle GPU/CPU state,
running jobs, queues and available memory before each session; never run
formal timing alongside other computation.

Long benchmarks live under `validation_test/optimization`. The Optuna 5 lane
uses `benchmark_optuna50_python.py` and `benchmark_matlab_optuna50.m`, records
the complete machine/runtime/load state, and checks seeded proposal checksums
before comparing warmed medians. Results produced against 4.x are historical
and cannot satisfy the 0.3.0 release gate. Parallel results remain calibrated
scheduler evidence rather than a promise that every CAE objective scales by
the same factor; cheap objectives should remain sequential.

Publication also uses the standalone exact-wheel CUDA lane on the supported GPU acceptance host. After the
successful `main` CI run, execute
`python tools/release_quad.py optuna-candidate --ci-run-id <id> --target all`,
then pass the retained wheel to
`python tools/release_quad.py optuna-done --wheel <path>`. Only after that gate
passes may the matching commit be tagged. The manual release workflow requires
the same CI run ID and the candidate SHA256, publishes that exact wheel to
PyPI, and attaches it to the matching GitHub Release.


## CUDA execution boundary (2026-10-11)

Build `optuna_mex` with `Build.ps1 -OptunaMexOnly`, which calls MATLAB
`mexcuda` with its supported compiler configuration. A CUDA-capable NVIDIA
host and Parallel Computing Toolbox are required for GPU-mode density evaluations.
No NGSolve, TaskManager, Radia CPU solver, or CPU-only MEX build is linked.

Parzen numerical/categorical mixture densities and univariate/joint/grouped
TPE candidate densities and selection execute in FP64 CUDA in GPU mode or
CPU host kernels in CPU mode. Study/Trial management,
stateful random-number consumption, candidate generation and Pareto/Sobol
remain host operations. Selection tie-breaking stays ordered on both routes. Device errors raise; they do not
switch to a CPU density implementation. RNG state is restored on failed
candidate evaluation. `api.info` reports this boundary and the build route;
MATLAB rejects old CPU-only MEX binaries even if their command list matches.

The oracle remains Optuna 5.0.0. Numerical parity and transfer-inclusive
measurements are required before claiming acceleration. Building with
mexcuda alone does not certify speedup or every sampler's GPU acceleration.
