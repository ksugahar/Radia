# Codex - Radia Project Policy

This file contains active decision rules only. Implementation details belong
in source, tests, package documentation, focused skills, and `radia-mcp`
knowledge. Historical investigations remain in Git history and must not be
copied back here.

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
- Build thin Python, MATLAB, and MCP adapters around coarse tested workflows.
  Do not expose every helper as a tool or duplicate solver logic in wrappers.
- Keep two genuinely independent analysis routes for important models when
  feasible.

### Numerical Rules

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

### MATLAB And Simulink

- Use MathWorks' official MATLAB MCP Server and Simulink Agentic Toolkit.
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

**POLICY**: **mdx** is Radia's self-hosted CI and preflight host. LAB and
100号機 are development machines. mdx gives CI and preflight work priority.
Long solver work should use hibino first when it is available and may use mdx
only when the mdx CI queue is idle.

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

The EqnEdit64 signed-standalone release job is the sole LAB runner exception:
it verifies the LAB-owned signing certificate and OneDrive release manifest
after successful EqnEdit64 tag CI. It is a release signing gate, not regular
repository CI, and must not acquire Radia test or build work.

### Compute Host Routing

**POLICY (hibino-first; mdx CI-first, 2026-09-03)**: Run solver-heavy
validation, optimization, scaling, memory, and timing work on hibino first.
Use mdx only when hibino is unavailable and both the mdx CI runner and its job
queue are idle. Compute work must never delay or destabilize CI/preflight.
Historical mdx measurements remain valid provenance. Record host, runtime,
versions, and measured quantities in validation JSON.

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
for the same commit. LAB and 100号機 return to verified canonical editable
installs after release.

## Optuna

Pinned upstream Optuna is the oracle for shared MATLAB behavior. Seed, options,
search-space order, history, constraints, values, states, warnings, and random
consumption must match upstream fixtures. Handwritten MATLAB output is not
compatibility truth.

MATLAB table/MAT storage, Simulink monitoring, parallel execution, and MEX are
extensions, not permission to alter the compatible algorithm. Keep API coverage
and oracle manifests current. Unsupported behavior fails loudly. Official
`optuna/optuna-mcp` owns generic Study/Trial MCP; `radia-mcp` owns MATLAB,
Simulink, MEX, and Radia-domain composition.

## Git And Agents

- Never reset, clean, stash, delete, or overwrite another session's WIP.
- Commit only reviewed files owned by the task; use clean integration worktrees.
- Inspect branches, PRs, and remote main before merge or rebase.
- Do not finish while required tests, CI, build, or release commands run.

Claude normally stops after a tested local commit and reports branch/SHA. Codex
owns push, CI, fix-forward, tags, publication, deployment, and `release-quad`,
unless the user explicitly assigns that work to Claude for the specific task.

## Detailed Guidance

Use focused skills and `packages/radia-mcp/src/radia_mcp/**/knowledge/` for
operational detail. In particular see `release-quad`, `simulink-app-health`,
`verify-deploy`, `gmsh-verify`, and `api-inventory`. If detailed guidance
conflicts with this file, update the stale detail and its focused regression.
Do not expand this file into a second manual.
