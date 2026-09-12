# Claude Code - Radia Project Policy

Keep active rules here; details belong in source, tests, package docs, skills,
and `radia-mcp` knowledge. Historical investigations remain in Git history.

## Mission

Radia is an AI-native electromagnetic CAE platform, not another standalone
solver. AI designs; Radia provides the engineering platform. Respect NGSolve as
the numerical foundation and extend it only for missing engineering capability.

Until explicitly closed, work is limited to:

1. Complete and validate HDiv-MMM topology optimization.
2. Improve repository and `radia-mcp` quality, maintenance, validation,
   build, CI, packaging, and release operation.

Do not add unrelated features during this maintenance program.

## Repository Boundaries

The monorepo contains independently released `radia`, `cubit-mesh-export`,
`radia-mcp`, `radia-optuna`, and `eqnedit64`. Keep commits and CI scoped
to the owning distribution. A shared file may trigger multiple lanes only when
it changes a real shared ABI, build, or integration contract.

- `src/`, `matlab/`, `packages/`: production implementation.
- `tests/`: fast deterministic bug and contract protection.
- `validation_test/`: numerical, performance, native, GUI, and multi-machine
  evidence with machine-readable result JSON.
- `docs/**/*.ipynb`: executed, result-bearing public demonstrations and
  presentation-ready narratives, not production interfaces or benchmarks.
- `C:\temp`: disposable prototypes and generated work.
- `examples/`: retired; never add files.

Do not commit logs, transient solver output, binary backups, lock copies,
generated inventories, or root scratch. Do not keep two tests with the same
purpose and failure signal.

## Engineering Architecture

Prefer established public abstractions over proprietary plumbing.

- NGSolve owns spaces, orientation, Piola maps, curved geometry, quadrature,
  assembly, `CoefficientFunction`, and `GridFunction` behavior.
- Radia owns electromagnetic formulations, open-boundary operators,
  constitutive models, HDiv-MMM/VIM, PEEC, SIBC/ESIM, Kelvin/DtN, stream
  functions, application coupling, and validated native kernels.
- Cubit and build123d own CAD authoring. Netgen/NGSolve and Cubit own solver
  mesh generation. Gmsh is Radia's post-processing target.
- Build thin Python/MATLAB/MCP adapters around tested domain workflows; do not duplicate solver logic or expose every helper.
  Consolidate with preserved contracts per `packages/radia-mcp/CONTRIBUTING.md`.
- Keep two genuinely independent analysis routes for important models when
  feasible.

### Numerical Rules

- Name eigenmode-bulk/surface coupling **Foster + SIBC**; use **CLN + SIBC**
  for an actual CLN/Krylov bulk basis. Do not call these scalar enriched-space
  models "mixed Galerkin". Legacy API/path identifiers remain compatible.
- Foster + SIBC is the production default for this bulk/surface method.
  CLN remains a comparison route, not a required migration. Decide retained
  mode counts from error/convergence in the operating band, not naming alone.

- Use SI units; magnetization is A/m.
- Compare vector fields with `norm(B1 - B2)`.
- Radia C++ matrix storage is row-major unless an external API owns the layout.
- Radia core Green functions remain Laplace/MQS/Darwin kernels.
- Use HACApK for Radia core compact interactions; do not restore core ExaFMM.
- NGSolve BEM is supported. Its native multipole/FMM facilities may be used for
  smooth free-space BEM and Biot-Savart work.
- Python FE work follows caller-owned `ngsolve.TaskManager()`. MATLAB
  parallel work follows MATLAB's runtime.
- Fail loudly on unsupported geometry, labels, ABI, convergence, or backend
  state. Never silently substitute a numerically different route.

## Interfaces

Python/MCP is the first-class AI interface. Masked blocks in the single Radia
Simulink library are the human production interface. Implement and study both
MCP+LLM and Simulink+MCP workflows; their relative effectiveness is an active
research question.

Standalone PySide/PyQt Radia panels and notebook workbenches are retired.
Coreform Cubit's private PySide6 is allowed only inside Cubit for the
`cubit-mesh-export` toolbar. Normal Radia Python must not depend on Qt.

### Shared MCP Runtime Ownership

- MCP is experimental development tooling, not a numerical solver release.
  Developers may edit live MCP source and change its editable source with
  `pip install -e`; no dedicated branch, frozen snapshot or separate deployment
  approval is required for routine MCP experiments.
- Coordinate overlapping edits and environment changes; preserve others' WIP
  and active CAD/MATLAB jobs. Never mass-kill processes to refresh MCP.
- Verify the interpreter and actual import path after repointing. Report each
  client's live loaded source separately; reload/reconnect affected clients as
  needed. Unknown or mixed evidence must stay unverified.
- Solver/native numerical acceptance is unchanged. Keep release installation
  tests isolated from development; see the
  [Shared MCP runtime policy](packages/radia-mcp/docs/operations/mcp-runtime-policy.md).

### Canonical Bibliography

- The single parent is `packages/radia-mcp/src/radia_mcp/bibliography/data/references.bib`.
  Resolve it with `bibliography_canonical_path`; correct verified entries there.
- Manuscript folders contain generated `.bbl` only, never local `.bib` copies.
  Use `bibliography_make_bbl`, regenerate after changes, and check citations.

### MATLAB And Simulink

- Use MathWorks' official MATLAB MCP Server and Simulink Agentic Toolkit.
  The official MCP Server remains the standard MATLAB operation foundation.
  LAB/100 track compatible stable releases: record versions before release,
  review upgrade notes and regression-test affected Radia composition contracts.
  Use the official MATLAB Engine for Python when available and appropriate,
  particularly for SSH batch execution and shared-session connections. Keep
  Engine installed and importable on mdx1/mdx2; record the interpreter, MATLAB
  version, and execution route. Reuse Radia MATLAB entry points through either
  route. Engine execution does not waive the official Toolkit model-edit,
  check, save/reopen, and diagnostic requirements for production SLX files.
  Close only sessions owned by the operation and retain caller-owned sessions.
- A tracked production `.slx` passes read, edit, check, save, close, and reopen
  on the exact path. Never patch SLX ZIP/XML directly.
- Public model-authored UI text is English. Mojibake, replacement glyphs,
  suspicious `???`, broken wiring, or unresolved blocks block release.
- Production blocks use readable Level-2 MATLAB S-Functions. Reusable numerical
  kernels remain independently callable standalone MEX functions.
- Native objects use checked `uint64` handles with type, generation,
  ownership, and liveness validation; never expose raw pointers.
- Python fallback is allowed at initialization, explicit update, artifact
  generation, or batch-solve boundaries, never per Simulink time step.
- Every public Python capability has a named MATLAB entry point or a checked
  parity-manifest classification. Prefer one C++ source of numerical truth with
  thin pybind11 and MEX adapters.
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

- Producing a `.vol` is `cubit-mesh-export`'s responsibility and runs on a
  licensed Cubit machine. radia CI never generates one: it consumes committed
  fixtures. A missing required fixture fails the test; only explicitly optional
  validation inputs may produce a visible skip. Heavier `.vol` work
  belongs to `validation_test/`, and `docs/**/*.ipynb` may show the Cubit
  generation step.
- Only `cubit-mesh-export` launches the Cubit GUI and owns Cubit GUI tests on licensed hosts;
  Radia solver/application/validation/CI lanes must not launch it or duplicate those tests.
- Prefer Cubit APREPRO/Python batch for production/validation CAD and meshes;
  consider Sculpt for suitable HEX domains; build123d/Netgen alternatives must be explicit.
- Radia normally reads checked `.vol` files. Generation uses APREPRO or Cubit's Python API
  in batch/headless mode with `cubit-mesh-export` owning export; CI remains fixture-only.
- `radia-mcp.cubit` supports human-AI collaboration through the `cubit-mesh-export` GUI;
  it does not authorize solver-side GUI launch or interruption of human-owned sessions.
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

**POLICY**: **mdx1 and mdx2** are Radia's self-hosted CI and preflight pool. LAB and
100号機 are development machines. Both mdx hosts give CI and preflight priority.
GitHub Actions uses the shared `mdx` label and assigns jobs to an available runner.
Release-quad requires LAB, 100号機, mdx1, and mdx2 for the same release commit.
LAB/100号機 retain verified editable installs; mdx1/mdx2 consume release wheels.
Release-quad deploys `cubit-mesh-export` and its Cubit plugin/toolbar only to
LAB and 100号機. Do not install or run them on mdx1/mdx2 in the release lane;
their Cubit version/compatibility probe fields are explicitly not applicable.
This deployment boundary does not prohibit Cubit-independent checker unit tests
in isolated CI environments. Existing installations are not silently removed.
hibino remains a computation host and is not a release-quad acceptance target.
Long solver work should use hibino first when it is available and may use mdx
only when the mdx CI queue is idle.

MATLAB-capable runner services must use a MATLAB-authenticated account and pass
an Engine startup/calculation/shutdown check in the actual service context;
SYSTEM and SSH success are not substitutes for that acceptance. Keep the
diagnostic manual-only. mdx1 is scheduled for retirement in March 2027; migrate
CI routing and release-quad targets together before retirement, following the
`release-quad` skill and a user-approved replacement plan.

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

`docs/**/*.ipynb` is the public calculation record. Saved output is sufficient;
an adjacent JSON and a runtime gate are not required. A docs-only contract lane
parses changed notebooks. Public examples include saved parameterized WebGUI
geometry/mesh and primary-field scenes.

Developer pre-push hooks run only the impact-scoped mdx preflight. Release
workflows, never developer hooks, publish immutable artifacts.

No workflow selects a LAB runner. The EqnEdit64 signed-standalone release job
once did, to read the LAB-owned OneDrive release manifest, and that was never
reachable: a runner service executes as NETWORK SERVICE, which sees no per-user
mapped drive, and these machines are the workgroup HFEM rather than a domain,
so the laboratory share refuses its machine account outright. Making LAB a
runner would also hand it Radia's heavy build and test work. The signed
executable now reaches CI as an asset of the eqnedit64-staging release, staged
by sync_to_o.ps1 in the same transaction that updates O:, and the whole release
lane runs GitHub-hosted. O: stays the human hand-test entry point, not a CI
input.

### Compute Host Routing

**POLICY (2026-09-08)**: Run solver-heavy validation, optimization, scaling,
memory, and timing work on hibino when it is **already running and idle**.
hibino is a SPOT instance: starting it is a human action, not an agent one.
Probe with `ssh -o ConnectTimeout=6 -o BatchMode=yes hibino hostname` — ICMP is
blocked, so `ping` reports a false "down" — then check for a running python
job, because hibino takes one heavy job at a time. Otherwise use whichever of
mdx1/mdx2 is idle, after checking both its CI runner and its job queue are idle.
Compute work must never delay or destabilize CI/preflight.
Historical mdx measurements remain valid provenance. Record host, runtime,
versions, and measured quantities in validation JSON.

Core count does not imply speed; settle the setting before the host. One #6
Gram build moved from 51 s to 1567 s on a single host by quadrature rule alone
(`validation_test/esrf_three_engine/results/hex_gram_definiteness_*.json`).
hibino's real advantage is memory — 230 GB with no pagefile — rather than its
76 logical cores.

**POLICY**: 全てのベンチマークスクリプトは機械可読な JSON 結果を保存すること。

## Build And Release

- Dependency versions live in package metadata and CI; do not duplicate pins
  here.
- Install MKL from its supported package dependency. Do not bundle MKL or an
  Intel-specific NumPy build; preserve user NumPy ownership.
- `radia_motor_rom.dll` is Radia-owned. Release wheels contain no third-party DLLs
  except dependencies explicitly permitted by packaging policy.
- Build with current `Build.ps1`/CMake. Remove superseded recipes.
- Standalone MEX and Level-2 S-Function layers have separate API, error,
  lifecycle, repeated-run, and performance responsibilities.

Use `tools/release_quad.py` and the `release-quad` skill. Publish only when
CI, exact package hashes, native/MEX/SLX checks, and required machine gates pass
for the same commit. LAB and 100号機 retain approved, verified editable sources
after release; source changes follow Shared MCP Runtime Ownership, not an automatic reset.
Before tagging, dispatch `Radia Native Release` on the exact release SHA; `ci-verify` requires its successful native check.

## Optuna

Pinned Optuna 5.0.0 is the sole oracle for shared MATLAB behavior: seed, options,
parameter/search-space order, history, named constraints, values, states, warnings,
defaults, extension points, and random consumption match upstream-generated fixtures.
Handwritten MATLAB output is not compatibility truth.
Adopt the Optuna 5 design in place: delete active optuna49 fixtures/pins, removed APIs/options,
legacy sampler-state restore and public multi-objective TPE; unified TPESampler owns both objectives.
Table/MAT storage, Simulink, parallelism and MEX are extensions, not algorithm changes.
Native acceleration requires correct differential results; maintain API/oracle manifests and fail loudly.
Fast deterministic tests belong in tests; long performance/scaling/parallel/dimension work in validation_test.
Released optuna-mcp==0.2.0 owns generic Study/Trial MCP for Optuna 5.0.0; 0.3.0.dev is not a stable release.
radia-mcp owns MATLAB, Simulink, MEX, differential-oracle, performance-gate and Radia composition.

## Git And Agents

- Never reset, clean, stash, delete, or overwrite another session's WIP.
- Commit only reviewed files owned by the task; use clean integration worktrees.
- Inspect branches, PRs, and remote main before merge or rebase.
- Do not finish while required tests, CI, build, or release commands run.

Claude normally stops after a tested local commit and reports branch/SHA. Codex
owns push, CI, fix-forward, tags, publication, deployment, and `release-quad`,
unless the user explicitly assigns that work to Claude for the specific task.

## Detailed Guidance

Use focused skills and `packages/radia-mcp/src/radia_mcp/**/knowledge/` for detail; correct policy conflicts there with a focused regression.
Do not expand this file into a second manual.
