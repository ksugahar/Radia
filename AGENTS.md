# Codex - Radia Project Policy

Investigations remain in Git.
## Mission
Radia is AI-native, AI-operable electromagnetic CAE through MCP. Users set analysis plans, physical assumptions and acceptance criteria; AI assists tool operation and execution. Do not market autonomous AI design. Extend NGSolve only for missing engineering capability.
Scope: (1) complete and validate HDiv-MMM topology optimization; (2) improve repository and `radia-mcp` quality, validation, build, CI, packaging and release operation.
Do not add unrelated features. No backward compatibility: keep only the latest design; delete legacy paths, aliases, flags, shims and old releases; repository health comes first.
## Repository Boundaries
The monorepo independently releases `radia`, `radia-simulink`, `cubit-mesh-export`, `radia-mcp` and `radia-optuna`; EqnEdit64 lives in its own repository. Scope commits and CI to the owner; each PyPI distribution owns its CI lane in `tools/ci_distribution_manifest.json`, and shared files trigger multiple lanes only for real shared ABI, build or integration changes.
`radia`: Python/pybind11 CPU and GPU, Windows and Linux; its code and tests never import `radia_simulink`. `packages/radia-simulink`: all MATLAB/Simulink/LTspice source, build and distribution, depending on `radia` one-way (the retired `radia-gpu` name is not published). `radia-mcp`: cross-package knowledge/tools, optional GPU backend. `radia-optuna` stays independent. `radia` releases exclude mexcuda, MATLAB, LTspice and Simulink artifacts and their CI gates; shared and CPU numerical gates remain mandatory.

- `src/`, `packages/`: production; integrations belong to their package owner.
- `tests/`: protect against breakage by future changes, every push, impact-selected: surrogate tests (an adopted case's coarsest convergence level vs its stored baseline; ~1 % off analytic is fine: it guards code, not accuracy) and lint-like contracts. Docstrings state `Guards:`/`Assumes:`; never independent validation.
- `validation_test/`: formulation soundness: p-convergence against cited, CAS-checked analytic solutions, independent solvers or measurements, plus CPU/GPU and platform parity; runs on triggers (covered code or dependency change, surrogate-test failure, baseline regeneration, adoption) and every adopted case gates release on the exact wheel (failed or incomplete blocks it; release notes mark validated vs experimental features).
- `validation_records/`: adoption decisions (user-approved criterion, rationale, evidence hashes) and champion data (finest/highest-order runs, speed and scale records) with version, source SHA and conditions; marked stale, never deleted, when covered code changes; CI checks schema only. Timings never gate.
- `docs/**/*.ipynb`: executed notebooks explaining Radia's implementation; not production interfaces or benchmarks.
- Job-owned scratch (`C:\temp`): disposable prototypes and generated work. `examples/`: retired; never add files.

Do not commit logs, transient solver output, binary backups, lock copies, generated inventories or root scratch. Commercial converters stay in independent private repositories, separate from MCP servers.
## Engineering Architecture
Prefer established public abstractions over proprietary plumbing.
- NGSolve owns spaces, orientation, Piola maps, curved geometry, quadrature, assembly, `CoefficientFunction` and `GridFunction` behavior.
- Radia owns electromagnetic formulations, open-boundary operators, constitutive models, HDiv-MMM/VIM, PEEC, SIBC/ESIM, Kelvin/DtN, stream functions, application coupling and validated native kernels.
- Cubit and build123d own CAD authoring; Netgen/NGSolve and Cubit own solver meshing; Gmsh is the post-processing target. Modeling APIs use neutral domain names: no reference commercial-solver brands in identifiers, docs, tests or commit messages, and no vendor API compatibility claims.
- Build thin Python/MATLAB/MCP adapters around tested domain workflows; do not duplicate solver logic or expose every helper (contracts per `packages/radia-mcp/CONTRIBUTING.md`). Keep two genuinely independent analysis routes for important models when feasible.
### Numerical Rules
- Radia/radia-mcp must not include implementations infringing CLN patents; preserve citations and independent methods, renaming is not removal, and select reduced order by measured operating-band error. Follow [CLN implementation policy](docs/CLN_POLICY.md).
- Call eigenmode-bulk + surface-impedance coupling **Modal-bulk SIBC** (formerly Foster + SIBC), never "mixed Galerkin".
- Use SI units; magnetization is A/m. Radia C++ matrix storage is row-major unless an external API owns the layout.
- PARDISO is not used in principle: direct FE solves use `sparsecholesky`, never a PARDISO fallback; MKL serves dense BLAS/LAPACK/FFT only. Linear true residual <= RELATIVE_LIMIT (1e-6), else raise; nonlinear/validation gates stay unchanged.
- Compare vector fields with `norm(B1 - B2)`; ESRF coil-yoke HDiv/FEM validation defaults to 1 % relative RMS on its declared core stencil, separately from solver convergence.
- Radia core Green functions remain Laplace/MQS/Darwin kernels. Use HACApK for core compact interactions; do not restore core ExaFMM. NGSolve BEM and its multipole/FMM facilities may serve smooth free-space BEM and Biot-Savart work.
- Python FE work follows caller-owned `ngsolve.TaskManager()`. Exception: `p1_newton.solve_p1_newton` and `p1_linear.solve_p1_linear` own regions and reject active TaskManager for AMS setup; helpers stay caller-owned.
- Fail loudly on unsupported geometry, labels, ABI, convergence or backend state. Never silently substitute a numerically different route.
## Interfaces
Python/MCP is the first-class AI interface. Operational workflows use MCP or human-AI collaboration through Cubit for CAD/mesh and masked blocks in the single Radia Simulink library; the effectiveness of MCP+LLM vs Simulink+MCP is an open research question. Standalone PySide/PyQt panels and notebook workbenches are retired; Cubit's private PySide6 serves only the `cubit-mesh-export` toolbar, and normal Radia Python must not depend on Qt.

### Shared MCP Runtime Ownership
- radia-mcp releases independently of Radia: it deploys to the build host (editable) and the validation host (verified wheel), never to the JMAG evaluation host or retired execution targets; isolated CI tests are not deployment.
- Routine editable updates: update the usual source, reconnect, check the live source and one harmless affected tool. Pure Python edits need no reinstall; relocation or dependency/metadata changes do.
- MCP updates are forward-only: fix forward, never revert/reinstall/redistribute an old version or keep rollback copies; remove obsolete copies after consumer and unique-work checks.
- MCP is experimental development tooling, not a numerical solver release. Developers may edit live MCP source and change its editable source with `pip install -e`; no dedicated branch, frozen snapshot or separate deployment approval is required for routine MCP experiments.
- Coordinate overlapping edits and environment changes; preserve others' WIP and active CAD/MATLAB jobs. Never mass-kill processes to refresh MCP.
- Verify the interpreter and actual import path after repointing. Report each client's live loaded source separately; reconnect as needed. Unknown or mixed evidence must stay unverified.
- Keep solver release tests isolated from development. See [Shared MCP runtime policy](packages/radia-mcp/docs/operations/mcp-runtime-policy.md).

### Canonical Bibliography
The single parent is `packages/radia-mcp/src/radia_mcp/bibliography/data/references.bib`; resolve it with `bibliography_canonical_path` and correct verified entries there. Manuscript folders contain generated `.bbl` only, never local `.bib` copies: use `bibliography_make_bbl`, regenerate after changes and check citations.

### CPU / GPU Execution Policy (2026-10-11)
- CPU: Python (pybind11)/NGSolve (Windows: Visual Studio Build Tools); NGSolve TaskManager owns CPU parallel execution under the caller-owned region policy. A mexcuda-built CPU route is permitted only when selected explicitly and proven to create no CUDA context; no separate CPU-only MEX build.
- radia first: functions are developed and validated in `radia` (CPU, Python/pybind11); `radia-simulink` exposes them to MATLAB/Simulink through thin MEX gateways and blocks. `radia-simulink` kernels are C++/CUDA built with MATLAB mexcuda (NGSolve-inherited DLLs allowed), developed and tested first through a pybind11 binding; Equivalence: (1) pybind11 vs MEX on the same kernel (tests, bitwise); (2) radia CPU reference vs GPU kernel (validation_test, tolerance). Exceptions (Simulink-only behavior such as sample time, ports, S-Function lifecycle, LTspice; GPU-only kernels without a CPU counterpart) are declared in the parity manifest with a reason and validated separately. Python CPU accuracy validation is never a GPU-user prerequisite.
- mexcuda modules may retain CPU host control. Never include, link, or use NGSolve TaskManager in mexcuda builds. MATLAB owns host scheduling, CUDA owns device execution; judge bounded OpenMP host work by numerical parity and transfer-inclusive timing.
- Windows builds: VS Build Tools 2022/MSVC 14.44 (mexcuda-compatible), MATLAB-bundled CUDA 12.8; Linux `radia` GPU: nvcc 12.8; no other CUDA/compiler overrides. BLAS: MKL (x86), NVPL (aarch64). Python GPU runtime (2026-10-21): exact-version pip NVIDIA CUDA 12.8 runtime packages, loaded lazily on GPU paths only; wheels do not redistribute the CUDA runtime.
- No implicit CPU/GPU/Python fallback: select explicitly (`Backend="gpu"|"cpu"`), record the backend and raise on failure. Record parity, artifact hashes and validation JSON; Poisson does not certify ICCG/AMS/HDiv-MMM.
- `radia-simulink` (CPU/GPU) is Windows-only: built on the build host, accepted on its RTX GPU with independent review, GPU-deployed to a Windows A100 platform. Prove CPU/GPU agreement on small publishable cases before scaling; backend names are not speedup evidence. Data is HDF5 only (schema/version, SI units, index base, labels, dtypes, UTF-8) plus NGSolve `.vol`/`.sol`; no MAT. Validate readers/writers above 2 GB, never split. Claim GPU portability only after validated architecture/PTX compatibility and a rebuild per OS/MATLAB ABI.

### MATLAB And Simulink
- Use MathWorks' official MATLAB MCP Server and Simulink Agentic Toolkit; record stable tool versions and regression-test upgrades. MATLAB Engine sessions never waive Toolkit model-edit, check, save/reopen or diagnostic gates. Close only owned sessions.
- A tracked production `.slx` passes read, edit, check, save, close and reopen on the exact path; never patch SLX ZIP/XML directly. Public model-authored UI text is English; mojibake, broken wiring or unresolved blocks block release.
- Production blocks are readable, user-customizable Level-2 MATLAB S-Functions; reusable kernels remain independently callable standalone MEX functions. Native objects use checked `uint64` handles; never expose raw pointers.
- Explicit Python CPU workflows run at initialization, update, artifact or batch boundaries, never per Simulink step or as implicit CUDA fallback. Measure cold start, warmed median, transfers and memory; do not assume MEX or pybind11 is faster.
- Production models require typed ports, masks, sample-time semantics, dependency checks, lifecycle tests, numerical checks and durable `run.log`/`result.json`. Spatial runs emit checked Gmsh `.msh` v4.2 artifacts.

## Mesh And CAD
STEP is the portable CAD standard (SAT for Cubit's ACIS workflow); the solver boundary is a checked `.vol` regardless of the creation route.
- Cubit-generated `.vol` files belong to `cubit-mesh-export` on a licensed Cubit machine; build123d + Netgen is the explicit alternative. radia CI consumes committed fixtures: a missing required fixture fails the test; only explicitly optional validation inputs may skip visibly.
- Drive Cubit headless (APREPRO/Python CLI) and export through `cubit-mesh-export`. LLM/MCP Cubit execution is always batch/nographics and fails rather than launching or attaching to a GUI; only `cubit-mesh-export` owns the human-facing GUI release test, and journal/checkpoint handoff is artifact-only.
- `cubit_mesh_export.mcp` owns Cubit MCP/APIs without Radia; each distribution owns its runtime helpers and requires neither the other MCP package nor a shared private distribution. Cubit export stays in the C++ `.ccm` backend; the only Radia-owned VTK export surface is `cubit-mesh-export`.
- Every solver-bound `.vol` passes `check-vol` with its versioned label contract before solver or Simulink initialization. Label checks validate topology/naming; DesignSpec validates physical data. Preserve canonical CAD, mesh, contracts and generation sources; remove only verified scratch or superseded output.

## Evidence And CI

### CI Execution, Validation Evidence, and Notebook Policy (2026-09-03)

**POLICY**: CI and builds run on the build host with the `100` runner label, isolated from user runtimes; `ci-light` serves fast contracts, `100` heavy CPU jobs, plus one GPU label. Use the on-demand GPU worker only for mexcuda GPU work.
Deployment (2026-10-30): normal Radia operation uses explicitly bound editable sources on both build and validation roles. Exact-wheel release QA remains isolated on both roles. radia-mcp independently remains editable on the build role and a verified wheel on the validation role; cubit-mesh-export stays build-role-only. Solver releases never install, uninstall, repoint or version-gate those independent packages. Retired targets remain excluded.
The build host owns preflight and native artifact builds. Each normal editable runtime has an explicit private hostname/interpreter/source binding and is separate from exact-wheel QA. The validation role runs without Cubit/cme. NGSolve uses the version pinned in package metadata (6.2.2608).
NGSolve hosts start from checked exported .vol files, without Cubit; the validation host runs without Cubit.
MATLAB-capable runners use a MATLAB-authenticated interactive account and pass an Engine startup/calculation/shutdown check in the actual runner context; never a MATLAB runner service.

CI scope begins at the independently released distribution boundary. The checked `radia_mcp.meta` catalog and each server's live `tools/list` response are the tool-discovery source of truth. Generated tool inventory snapshots are local diagnostics, not committed artifacts or CI oracles. normal pull-request and main-push CI runs a stable compact contract set, tests selected from the changed source/test paths, and only the affected server selftests. The complete package pytest suite, all-server selftests, and live-catalog audit belong to an explicit full-audit workflow.

normal CI optimizes for fast, high-signal feedback and does not automatically rerun a failed deterministic test. Do not keep two tests whose purpose and failure signal are the same. When CI exposes numerical uncertainty, run only the relevant validation lane and retain its result JSON.

`docs/**/*.ipynb` introduces features through saved calculation records. Saved output is sufficient; an adjacent JSON and a runtime gate are not required. A docs-only contract lane parses changed notebooks. Public examples include saved parameterized WebGUI geometry/mesh and primary-field scenes.

Pre-push runs only impact-scoped build-target preflight. Release workflows, never developer hooks, publish immutable artifacts.
No workflow selects a retired or user-desktop runner.

### Tests For AI-Authored Code

Radia is implemented only by AI agents; humans do not write code. Tests therefore target AI-typical failures, not human slips:
- Acceptance criteria (validation thresholds, reference values, tolerances) are human-owned. Changing them requires explicit user approval; agents may only propose, with evidence.
- Lifecycle: computation time is scale. A feature is adopted when its user-approved criterion passes in `validation_test/`; the coarsest level becomes a `tests/` surrogate test calling the same production entries (<= 30 s target; runtime recorded); at creation and regeneration a temporary production mutation per guarded feature must fail it and pass once restored. Lint-like contracts: policy, packaging, API, fail-loud. Development unit tests are triaged at PR time: contract, validation or deletion. Tests are fluid, not a constitution: a new test may retire those it supersedes, a suite that only grows is a defect, overlap across paths is welcome; only full duplicates (same purpose, quantity, model and failure signal) are retired. Baselines regenerate only after validation re-certification, with review (reason, diff). Explanatory milestones go to `docs/**/*.ipynb`, validated workflows to radia-mcp.
- Verified-reference rule: test a new implementation on identical inputs against a validated program, not its own output (MEX/GPU vs `radia` Python CPU; radia-optuna vs pinned upstream Optuna), in CI on every build. Its author never certifies expected values or baselines; reciprocal review checks oracle independence.
- Assertions report actual and expected values. Prefer behavioral tests and generated single sources of truth over string-pinned text checks.
- Never claim results that were not executed; evidence carries the source SHA and artifact hashes. Keep files small enough for an agent's context.

Every AI output passes at least one independent check: an agent that is about 90 % right becomes about 99 % right after one. A check counts only if independent (a different model or agent, or independent evidence), actually executed (not a summary read) and able to fail. Scale checks by impact: routine changes take one check (reciprocal review plus CI); releases, deletions and acceptance-criteria changes take at least two, and acceptance-criteria changes also need user approval.
Reliability comes from cheap independent checks on every change, not exhaustive test enumeration. Aim for the necessary, sufficient and minimal set: two independent detectors where code changes often, one where it is mature; no just-in-case variants.

### Compute Host Routing

Roles: the build host owns CI, release artifacts, isolated numerical validation and CPU/GPU development; the validation role operates an editable Radia source with separate exact-wheel QA and diagnostic timing, never Cubit/cme; the JMAG evaluation host runs no Radia work. CPU-compute and GB200 GPU-compute release cases are build-only unless separately authorized for execution. Timings are diagnostic.
Before remote work, probe reachability, running jobs and available memory. Use an idle host and one heavy job at a time; preserve active timing jobs and never stop another session's jobs.
Heavy jobs (2026-10-24): development venv read-only, never the release runtime; a job venv only for a differing pin, admitted and deleted; use only owner-pool cores, never CI runner cores; obtain manager resource admission before start and record it in `owner.json` in the job root; record topology, affinity, NUMA node and thread settings in result JSON. Runs longer than about 30 minutes are never PR gates; one-off comparisons are not repeated.
Compare CPU/GPU on the same host with fixed, recorded threads (16-thread baseline where supported). Measure speed, scaling, setup, transfers and total time; save hardware, memory and versions in JSON.
Historical measurements remain valid provenance: public evidence uses neutral roles; role-to-host maps, originals and hashes stay in private storage. Never relabel to another host. Keep hosts, IPs and paths private.

**POLICY**: Every benchmark script saves a machine-readable JSON result.

### Compute Scratch Cleanup

**POLICY (2026-10-13)**: Remote compute jobs finish only after verified recovery and cleanup of their job-owned scratch. Recover JSON, logs, commands, identities, hashes and non-reconstructible inputs to private storage; verify hashes before deleting the remote copy, then verify and report the deletion. Never blanket-delete scratch: preserve active/queued jobs and CI, resolve approved absolute targets and reject escaping links/junctions. Retention needs an explicit owner, reason and cleanup trigger; compute hosts are not archives.

## Build And Release

- Build ownership (2026-10-06): the build host builds solver candidates and release artifacts from an isolated reviewed checkout using current `Build.ps1`/CMake; remove superseded recipes. Other compute hosts consume the same hash-verified artifacts; never replace another session's runtime.
- Dependency versions live in package metadata and CI; do not duplicate pins here. Retire old builds/fallbacks/toolchains after consumer/unique-work checks; no rollback reserves. History is provenance, not runtime.
- Install MKL from its supported package dependency; do not bundle MKL or an Intel-specific NumPy build. `radia_motor_rom.dll` is Radia-owned. Release wheels contain no third-party DLLs except dependencies explicitly permitted by packaging policy.
- Use the `release-quad` wrapper for four source-bound cases: build-role editable, validation-role editable, CPU-compute build and GB200 GPU-compute build. Host-specific bindings and skills stay private. The controller remains `tools/release_quad.py`; isolated exact-wheel and scientific gates are separate from normal editable deployment. Publish only when CI, exact package hashes, the owning distribution's native/artifact checks and required machine gates pass for the same commit on the build host; on-demand validation-host timing is separate. Verify fresh imports and live sources after deployment. Before tagging, dispatch `Radia Native Release` on the exact release SHA; `ci-verify` requires its successful native check.

## Optuna
`radia-optuna` uses GPU numerical kernels compiled with MATLAB mexcuda. Keep management, sequential RNG and work without a speed requirement on CPU. Choose GPU work by measured transfer-inclusive cost; no implicit fallback for declared GPU operations.
Pinned Optuna 5.0.0 is the sole oracle for shared MATLAB behavior: seed, options,
parameter/search-space order, history, named constraints, values, states, warnings,
defaults, extension points, and random consumption match upstream-generated fixtures.
Handwritten MATLAB output is not compatibility truth.
Adopt the Optuna 5 design in place: delete active optuna49 fixtures/pins, removed APIs/options,
legacy sampler-state restore and public multi-objective TPE; unified TPESampler owns both objectives.
Table/HDF5 storage, Simulink, parallelism and MEX are extensions, not algorithm changes.
Maintain API/oracle manifests and fail loudly.
Released optuna-mcp==0.2.0 owns generic Study/Trial MCP for Optuna 5.0.0.
radia-mcp owns MATLAB, Simulink, MEX, differential-oracle, performance-gate and Radia composition.

## Git And Agents

- Never reset, clean, stash, delete, or overwrite another session's WIP.
- Commit only reviewed files owned by the task; use clean integration worktrees. Inspect branches, PRs and remote main before merge or rebase.
- Do not finish while required tests, CI, build, or release commands run.
- Every candidate is reviewed by a different agent before merge; each review runs the shared CI-equivalent precheck (`tools/ci_precheck.py`), not a diff read alone.

Domain knowledge lives in repo knowledge/skills, not session memory. Session roles, review transport and host maps are private local configuration, not repository policy. User assignments prevail.

## Detailed Guidance

Details: focused skills and `packages/radia-mcp/src/radia_mcp/**/knowledge/`. Fix conflicts with focused regressions. Do not expand this file into a second manual.

## Private Commercial Reference Comparisons (2026-10-07)

Commercial comparisons, including all inputs, outputs, values, figures, logs, reviews and provenance, are private lab data. Exclude them from public commits, docs, packages, releases, CI, papers, presentations and links; anonymization never permits publication. Transfer only general algorithms without comparison data, validated independently on analytical, public or publishable cases. Publication needs an explicit user policy change and license clearance. Applies to Codex and Claude; never a reason to delete, move or rewrite existing artifacts or history.
