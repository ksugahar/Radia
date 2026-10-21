# Claude Code - Radia Project Policy

Investigations remain in Git.
## Mission
Radia is AI-native, AI-operable electromagnetic CAE through MCP. Users set analysis plans, physical assumptions and acceptance criteria; AI assists tool operation and execution.
Do not market autonomous AI design. Extend NGSolve only for missing engineering capability.

Scope:

1. Complete and validate HDiv-MMM topology optimization.
2. Improve repository and `radia-mcp` quality, maintenance, validation,
   build, CI, packaging, and release operation.

Do not add unrelated features. No backward compatibility: keep only the latest design; delete legacy paths, aliases, flags, shims and old releases; repository health comes first.
## Repository Boundaries
The monorepo independently releases `radia`, `radia-gpu`, `cubit-mesh-export`, `radia-mcp`, `radia-optuna` and `eqnedit64`.
Scope commits to the owner. Each PyPI distribution owns its CI lane in `tools/ci_distribution_manifest.json`; shared files trigger multiple lanes only for real shared ABI, build or integration changes. Migration exceptions only shrink.
`radia`: CPU Python/pybind11 only. `packages/radia-gpu`: all MATLAB/Simulink/LTspice source, dependencies, build and distribution. `radia-mcp`: cross-package knowledge/tools, optional GPU backend. Keep `radia-optuna` independent.
CPU `radia` releases exclude mexcuda, MATLAB, LTspice and Simulink artifacts and their CI gates. These belong to independent distribution releases; shared and CPU numerical gates remain mandatory.

- `src/`, `packages/`: production; integrations belong to their package owner.
- `tests/`: fast deterministic bug/contract protection; every PR.
- `validation_test/`: correctness gates (physics, CPU/GPU parity, platform agreement); nightly and on the exact release artifact.
- `validation_records/`: once-passed evidence (speedups, one-offs): script, source SHA, artifact hash, result JSON; CI checks schema. Timings never gate.
- `docs/**/*.ipynb`: executed, saved feature introductions and presentation-ready
  narratives linking those samples; not production interfaces or benchmarks.
- `C:\temp`: disposable prototypes and generated work.
- `examples/`: retired; never add files.

Do not commit logs, transient solver output, binary backups, lock copies, generated inventories or root scratch.
Do not keep two tests with the same purpose and failure signal.  Commercial converters stay in independent private repositories, separate from MCP servers.
## Engineering Architecture
Prefer established public abstractions over proprietary plumbing.

- NGSolve owns spaces, orientation, Piola maps, curved geometry, quadrature,
  assembly, `CoefficientFunction`, and `GridFunction` behavior.
- Radia owns electromagnetic formulations, open-boundary operators,
  constitutive models, HDiv-MMM/VIM, PEEC, SIBC/ESIM, Kelvin/DtN, stream
  functions, application coupling, and validated native kernels.
- Cubit and build123d own CAD authoring. Netgen/NGSolve and Cubit own solver
  mesh generation. Gmsh is Radia's post-processing target. Electromagnetic modeling APIs use neutral domain names; do not introduce reference commercial-solver brands into repository identifiers, documentation, tests, or commit messages, or claim vendor API compatibility.
- Build thin Python/MATLAB/MCP adapters around tested domain workflows; do not duplicate solver logic or expose every helper.
  Consolidate with preserved contracts per `packages/radia-mcp/CONTRIBUTING.md`.
- Keep two genuinely independent analysis routes for important models when
  feasible.
### Numerical Rules
- Radia/radia-mcp must not include implementations infringing CLN patents.
  Preserve citations and independent methods; renaming is not removal.
  Use measured operating-band error to select reduced order.
  Follow [CLN implementation policy](docs/CLN_POLICY.md).
- Call eigenmode-bulk + surface-impedance coupling **Modal-bulk SIBC** (formerly Foster + SIBC), never "mixed Galerkin".

- Use SI units; magnetization is A/m.
- PARDISO is not used in principle: Python/MEX direct FE solves/subsolves use `sparsecholesky`, never a PARDISO fallback.
  MKL serves dense BLAS/LAPACK/FFT only.
  Linear true residual <= RELATIVE_LIMIT (1e-6), else raise. Nonlinear/validation gates stay unchanged.
  Residual checks cannot prevent allocation/index overflow.
- Compare vector fields with `norm(B1 - B2)`; ESRF coil-yoke HDiv/FEM validation defaults to 1 % relative RMS on its declared core stencil, separately from solver convergence.
- Radia C++ matrix storage is row-major unless an external API owns the layout.
- Radia core Green functions remain Laplace/MQS/Darwin kernels.
- Use HACApK for Radia core compact interactions; do not restore core ExaFMM.
- NGSolve BEM is supported. Its native multipole/FMM facilities may be used for
  smooth free-space BEM and Biot-Savart work.
- Python FE work follows caller-owned `ngsolve.TaskManager()`. Exception: `p1_newton.solve_p1_newton` and `p1_linear.solve_p1_linear` own regions and reject active TaskManager for AMS setup; helpers stay caller-owned.
- Fail loudly on unsupported geometry, labels, ABI, convergence, or backend
  state. Never silently substitute a numerically different route.
## Interfaces
Python/MCP is the first-class AI interface. Operational workflows use MCP
or human–AI collaboration through Cubit for CAD/mesh and masked blocks in the
single Radia Simulink library for application models. Implement and study both
MCP+LLM and Simulink+MCP workflows; their relative effectiveness is an active
research question.

Standalone PySide/PyQt Radia panels and notebook workbenches are retired.
Coreform Cubit's private PySide6 is allowed only inside Cubit for the
`cubit-mesh-export` toolbar. Normal Radia Python must not depend on Qt. LLM/MCP Cubit execution is always batch/nographics and fails rather than launching or attaching to a GUI; the only GUI run is the separately scoped human-facing release test, while journal/checkpoint handoff remains artifact-only.

### Shared MCP Runtime Ownership
- radia-mcp releases independently of Radia: targets build (editable) and NGSolve compute validation (wheel). Never deploy radia-mcp to the JMAG evaluation host or retired execution targets; isolated CI tests are not deployment. Do not wait for a solver release.
- Routine editable updates: update the usual source, reconnect, check the live source and one harmless affected tool. Pure Python edits need no reinstall; relocation or dependency/metadata changes do. No new daemon or mandatory ledger.
- MCP updates are forward-only: fix forward, never revert/reinstall/redistribute an old version or keep rollback copies; remove obsolete copies after consumer and unique-work checks.
- MCP is experimental development tooling, not a numerical solver release. Developers may edit live MCP source and change its editable source with
  `pip install -e`; no dedicated branch, frozen snapshot or separate deployment
  approval is required for routine MCP experiments.
- Coordinate overlapping edits and environment changes; preserve others' WIP and active CAD/MATLAB jobs. Never mass-kill processes to refresh MCP.
- Verify the interpreter and actual import path after repointing. Report each
  client's live loaded source separately; reload/reconnect affected clients as
  needed. Unknown or mixed evidence must stay unverified.
- Keep solver release tests isolated from development; numerical gates are unchanged.
  See [Shared MCP runtime policy](packages/radia-mcp/docs/operations/mcp-runtime-policy.md).

### Canonical Bibliography

- The single parent is `packages/radia-mcp/src/radia_mcp/bibliography/data/references.bib`.
  Resolve it with `bibliography_canonical_path`; correct verified entries there.
- Manuscript folders contain generated `.bbl` only, never local `.bib` copies.
  Use `bibliography_make_bbl`, regenerate after changes, and check citations.

### CPU / GPU Execution Policy (2026-10-11)
- CPU: Python (pybind11)/NGSolve, built with Visual Studio 2022 Build Tools on Windows; NGSolve TaskManager owns CPU parallel execution. CPU MEX is unsupported.
- GPU: `radia-gpu` uses MATLAB mexcuda (optional Simulink/LTspice adapters, no external products). GPU setup, preconditioners and execution stay in MATLAB/mexcuda without Python adapters or snapshots. NGSolve-inherited DLLs are permitted. Python provides independent CPU accuracy validation, never a GPU-user prerequisite.
- mexcuda modules may retain CPU host control. Never include, link, or use NGSolve TaskManager in mexcuda builds, except a read-only guard that fails if a foreign TaskManager is active. MATLAB owns host scheduling and CUDA owns device execution; evaluate bounded OpenMP host work using numerical parity and transfer-inclusive timing.
- All builds: VS2022/MSVC14.44, mexcuda-selected. CUDA: MATLAB-bundled12.8 only; no external CUDA/unsupported flags. Order: MATLAB 'a' release (no 'b'), supported compiler/CUDA, all hosts; upgrade yearly together, same Update, and rebuild. Verify interop and AMS/IH on both GPUs.
- No implicit CPU/GPU/Python fallback. When both routes exist, select explicitly (`Backend="gpu"|"cpu"`), record the backend and raise on failure. Internal CPU/GPU agreement does not replace Python CPU validation. Measure host setup and CPU solves separately.
- GPU development/builds/validation: A100, MATLAB-bundled CUDA only. Keep NGSolve off. 100 prioritizes CI/debt and same-wheel RTX 5000 release acceptance, with independent review; no GPU dev on 100. No GPU on 102.
- Record parity, artifact hashes and validation JSON. Poisson does not certify ICCG/AMS/HDiv-MMM.

### GPU Kernel Research Strategy (2026-10-11)

- Use mexcuda for GPU kernels with explicit assumptions. Implement missing Radia singular/near-field kernels in CUDA; preserve NGSolve FE ownership and HACApK policy.
- Data: HDF5 only (schema/version, SI units, index base, labels, dtypes, UTF-8); keep NGSolve `.vol`/`.sol`. No MAT or MAT fallback; migrate existing MAT per package, bit-exact with hashes. Never a mandatory GPU adapter; validate readers/writers above 2 GB, never split.
- Reuse gpuArray/BLAS where they implement the declared formulation. Keep GPU kernels independently callable. Validate iterative solvers for the actual operator.
- Establish CPU/GPU agreement and performance on small publishable cases on the local RTX GPU before scaling; save timings, transfers, memory, versions and true residuals in JSON. Backend names are not speedup evidence.
- Humans alone start/stop A100. Bounded jobs require resource admission and evidence recovery. Offline: explicit failure, no host fallback or indefinite PR queue.
- Target development and compute GPU architectures with validated PTX compatibility; rebuild per OS/MATLAB ABI and verify before claiming portability or scaling.

### MATLAB And Simulink

- Use MathWorks' official MATLAB MCP Server and Simulink Agentic Toolkit.
  Record stable tool versions and regression-test upgrades. Use official
  MATLAB Engine for SSH/batch or shared sessions when appropriate; record
  interpreter, MATLAB version and route. Engine does not waive Toolkit
  model-edit, check, save/reopen or diagnostic gates. Close only owned sessions.
- A tracked production `.slx` passes read, edit, check, save, close, and reopen
  on the exact path. Never patch SLX ZIP/XML directly.
- Public model-authored UI text is English. Mojibake, replacement glyphs,
  suspicious `???`, broken wiring, or unresolved blocks block release.
- Production blocks use readable Level-2 MATLAB S-Functions. Reusable numerical
  kernels remain independently callable standalone MEX functions.
- Native objects use checked `uint64` handles with type, generation,
  ownership, and liveness validation; never expose raw pointers.
- Explicit Python CPU workflows are allowed at initialization, update, artifact
  or batch boundaries, never per Simulink step or as implicit CUDA fallback.
- Classify capabilities as Python CPU, MATLAB CUDA GPU, or both in the parity
  manifest. Do not duplicate CPU MEX for parity. Share domain code and verify
  physical contracts without duplicating NGSolve FE semantics.
- Measure cold start, warmed median, transfer cost, memory, versions, and
  reliability; do not assume MEX or pybind11 is faster.
- Level-2 S-Function templates are user-customizable. Do not bury a reusable
  kernel exclusively in a MEX S-Function.

Production models require typed ports, masks, sample-time semantics, dependency
checks, lifecycle tests, numerical checks, and durable `run.log` /
`result.json`. Spatial runs emit checked Gmsh `.msh` v4.1 artifacts.

## Mesh And CAD

SAT is important for Cubit's ACIS workflow; STEP is the portable standard. The
solver boundary is a checked `.vol` regardless of the creation route.

- Cubit-generated `.vol` files belong to `cubit-mesh-export` on a licensed
  Cubit machine; without Cubit, use build123d + Netgen. radia CI consumes committed
  fixtures. A missing required fixture fails the test; only explicitly optional
  validation inputs may produce a visible skip. Heavier `.vol` work
  belongs to `validation_test/`, and `docs/**/*.ipynb` may show the Cubit
  generation step.
- Drive Cubit headless (APREPRO/Python CLI, never GUI; only `cubit-mesh-export` owns GUI release tests) for production/validation CAD and meshes, exporting through `cubit-mesh-export`; consider Sculpt for suitable HEX domains; build123d/Netgen alternatives must be explicit.
- `cubit_mesh_export.mcp` owns Cubit MCP/APIs without Radia; each distribution owns its runtime helpers and requires neither the other MCP package nor a shared private distribution.
  Radia may depend on Cubit for applicable workflows and prefers it when licensed; build123d + Netgen is the explicit alternative. MCP handoff is artifact-only, never GUI control.
- Every solver-bound `.vol` passes `check-vol` with its versioned label
  contract before solver or Simulink initialization.
- Label checks validate topology/naming; DesignSpec validates physical data.
- Cubit export stays in the C++ `.ccm` backend. Its embedded PySide6 toolbar
  collects options and invokes backend commands.
- The only Radia-owned VTK export surface is `cubit-mesh-export`. NGSolve may
  keep native `VTKOutput`; Radia application post-processing uses Gmsh.
- Preserve canonical CAD, mesh, contracts, and generation sources. Remove only
  verified scratch or superseded output.

## Evidence And CI

### CI Execution, Validation Evidence, and Notebook Policy (2026-09-03)

**POLICY**: Shared CPU CI/builds use 100号機 with the `100` runner label; GPU development/builds use A100. Stage3: fast-contracts uses `ci-light`, no `100`; heavy jobs keep `100`. Final: `100` + one capability.
100号機 permits isolated heavy validation.
Isolate CI from user runtimes. NGSolve compute validation consumes radia/radia-mcp wheels, never Cubit.
Build runs SSH preflight. Use the on-demand GPU worker only for mexcuda GPU work; Retired execution targets are excluded from deployment, CI and computation.
Release-quad: same SHA on build/validation; separate NGSolve timing.
CPU: build editable, NGSolve compute validation wheels. GPU: isolated reviewed source and verified artifacts.
`cubit-mesh-export`: build only, never compute validation; solver releases do not manage it.
NGSolve hosts start from checked exported .vol files, without Cubit; `A100` validation runs without Cubit. Never install/run Cubit there; checker unit tests may use isolated CI; no silent uninstalls.

MATLAB-capable runners use a MATLAB-authenticated account and pass an Engine
startup/calculation/shutdown check in the actual runner context. SYSTEM and SSH
success are not substitutes. Manual diagnostics; authenticated interactive `run.cmd`, never a MATLAB runner service.

CI scope begins at the independently released distribution boundary. The
checked `radia_mcp.meta` catalog and each server's live `tools/list` response
are the tool-discovery source of truth. Generated tool inventory snapshots are
local diagnostics, not committed artifacts or CI oracles. normal pull-request
and main-push CI runs a stable compact contract set, tests selected from the
changed source/test paths, and only the affected server selftests. The complete
package pytest suite, all-server selftests, and live-catalog audit belong to an
explicit full-audit workflow.

normal CI optimizes for fast, high-signal feedback and does not automatically
rerun a failed deterministic test. Do not keep two tests whose purpose and
failure signal are the same. When CI exposes numerical uncertainty, run only
the relevant validation lane and retain its result JSON.

`docs/**/*.ipynb` introduces features through saved calculation records. Saved output is sufficient;
an adjacent JSON and a runtime gate are not required. A docs-only contract lane
parses changed notebooks. Public examples include saved parameterized WebGUI
geometry/mesh and primary-field scenes.

Pre-push runs only impact-scoped build-target preflight. Release
workflows, never developer hooks, publish immutable artifacts.

No workflow selects a LAB runner. EqnEdit64 signed executables reach GitHub-hosted CI
through the eqnedit64-staging release, staged by sync_to_o.ps1 while updating O:.
O: is the human hand-test entry point, not a CI input; service accounts cannot use
LAB-owned mapped drives or the workgroup share. Do not restore that dependency.

### Compute Host Routing

`jmag`: JMAG evaluation only. Keep Python 3.12.10, SSH service/config/keys, VC++ runtime and `inetpub` for control from 100. No Radia deployment, builds, CI/preflight, solver runs, validation or timing.

**POLICY**: 100 owns CPU CI/debt; A100 owns GPU development (mexcuda only). Prioritize HDiv-MMM, then AMS. Compare AMS/ICCG on the same operator/residual gate, including setup/transfers. Retired targets are excluded from deployment, CI and computation.
Before A100 work, verify `ssh -o ConnectTimeout=6 -o BatchMode=yes A100 hostname`, jobs and memory. One heavy job at a time; never stop others' jobs.
Compare CPU/GPU on the same host with fixed, recorded threads (16-thread baseline where supported). Measure speed, scaling, setup, transfers and total time; save hardware, memory and versions in JSON.
100: CPU development/builds/CI and GPU acceptance; 100/LAB timings are diagnostic. Gmsh needs a desktop.
Measurements remain valid provenance: publish neutral roles, not hosts/IPs/paths; keep maps/originals/hashes private. Never relabel hosts. Recover unique work before authorized VM/disk deletion.

**POLICY**: 全てのベンチマークスクリプトは機械可読な JSON 結果を保存すること。

### Compute Scratch Cleanup

**POLICY (2026-10-13)**: Remote compute jobs, including NGSolve and retirement cleanup on jmag, finish only after verified recovery and cleanup of their job-owned `C:\temp` inputs, outputs, staging, environments and helpers.
- Recover JSON, logs, commands, runtime/source identity, hashes and non-reconstructible inputs
  to LAB durable storage. Numerical evidence goes in `validation_test/`, public demonstrations
  in executed `docs/**/*.ipynb`; raw operational logs/inventories stay outside tracked source.
- Verify destination hashes and commit intended evidence before deleting its remote copy,
  including failed-run diagnostics. Verify deletion and report locations and cleanup outcome.
- Never blanket-delete `C:\temp`: preserve active/queued jobs and CI; resolve approved absolute
  targets, reject escaping links/junctions and use literal-path operations.
- Retention requires an explicit owner, reason and cleanup trigger for unfinished work.
  Never silently retain completed work or reuse stale inputs. Compute hosts are not archives.

## Build And Release

- Build ownership: CPU artifacts on 100; GPU artifacts on A100, using
  reviewed isolated source and current `Build.ps1`/CMake; 100 owns local source.
  Consumers use hash-verified artifacts. No CPU builds on LAB or
  compute hosts; preserve others' runtimes. Reuse immutable verified-source evidence.
- Upgrade tools/dependencies and validate; no backward compatibility. Retire old builds/fallbacks/toolchains after consumer/unique-work checks; no rollback reserves. History is provenance, not runtime. Pins: metadata/CI.
- Install MKL from its supported package dependency. Do not bundle MKL or an
  Intel-specific NumPy build; preserve user NumPy ownership.
- `radia_motor_rom.dll` is Radia-owned. Release wheels contain no third-party DLLs
  except dependencies explicitly permitted by packaging policy.
- Build with current `Build.ps1`/CMake. Remove superseded recipes.
- Standalone MEX and Level-2 S-Function layers have separate API, error,
  lifecycle, repeated-run, and performance responsibilities.

Use `tools/release_quad.py` and the `release-quad` skill. Publish only when
CI, exact package hashes, native/MEX/SLX checks, and required machine gates pass
for the same commit. MCP: compute validation wheels; Cubit: build only. Verify fresh imports
and live sources.
Before tagging, dispatch `Radia Native Release` on the exact release SHA; `ci-verify` requires its successful native check.

## Optuna
`radia-optuna` uses GPU numerical kernels compiled with MATLAB mexcuda. Keep management, sequential RNG and work without a speed requirement on CPU. Choose GPU work using measured transfer-inclusive cost; compilation alone is not GPU execution or speedup. No CPU-only MEX build or implicit fallback for declared GPU operations.
Pinned Optuna 5.0.0 is the sole oracle for shared MATLAB behavior: seed, options,
parameter/search-space order, history, named constraints, values, states, warnings,
defaults, extension points, and random consumption match upstream-generated fixtures.
Handwritten MATLAB output is not compatibility truth.
Adopt the Optuna 5 design in place: delete active optuna49 fixtures/pins, removed APIs/options,
legacy sampler-state restore and public multi-objective TPE; unified TPESampler owns both objectives.
Table/HDF5 storage, Simulink, parallelism and MEX are extensions, not algorithm changes.
Native acceleration requires correct differential results; maintain API/oracle manifests and fail loudly.
Released optuna-mcp==0.2.0 owns generic Study/Trial MCP for Optuna 5.0.0; 0.3.0.dev is not a stable release.
radia-mcp owns MATLAB, Simulink, MEX, differential-oracle, performance-gate and Radia composition.

## Git And Agents

- Never reset, clean, stash, delete, or overwrite another session's WIP.
- Commit only reviewed files owned by the task; use clean integration worktrees.
- Inspect branches, PRs and remote main before merge or rebase.
- Do not finish while required tests, CI, build, or release commands run.
- Only Codex has task sessions, directed by its レポジトリ管理 hub.
  Claude Code has one レポジトリ管理 session, with read-only reviewer agents
  for deep reviews. Cross-platform exchange is only between these hubs.

Codex implements and validates. Claude owns independent review, GitHub
integration/CI/merge, release/deploy, `C:\temp` deletion and user decisions.
Domain knowledge lives in repo knowledge/skills, not session memory.
User assignments prevail.

## Detailed Guidance

Details: focused skills and `packages/radia-mcp/src/radia_mcp/**/knowledge/`. Fix conflicts with focused regressions. Do not expand this file into a second manual.

## Private Commercial Reference Comparisons (2026-10-07)

Commercial comparisons, including all inputs, outputs, values, figures, logs, reviews and provenance, are private lab data. Exclude them from public commits, docs, packages, releases, CI, papers, presentations and links; anonymization never permits publication.
Transfer only general algorithms without comparison data, validated independently on analytical, public or publishable cases. Publication needs an explicit user policy change and license clearance. Applies to Codex and Claude; never a reason to delete, move or rewrite existing artifacts or history.
