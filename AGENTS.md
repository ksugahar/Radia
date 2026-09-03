# Codex - Radia Project Development Guidelines

This document contains development guidelines and policies for the Radia project when working with Codex.

---

## Monorepo Structure

Radia is a **monorepo** containing all components. Only 2 repositories exist:

| Repository | Purpose |
|-----------|---------|
| **ksugahar/Radia** | Everything (C++, Python, MCP servers, sparsesolv, Cubit plugin) |
| **ksugahar/netgen** | Netgen fork (historical, archived — all features merged into official 6.2.2603) |

```
S:\Radia\01_GitHub\
  src/radia/              # Python package (pip install radia)
    _radia_pybind.pyd     # C++ extension
    mcp/                  # MCP servers (radia, ngsolve, cubit)
      radia/server.py     # mcp-server-radia
      ngsolve/server.py   # mcp-server-ngsolve
      cubit/server.py     # mcp-server-cubit
    panels/               # Headless calc_*.py, samples, artifact contracts, Cubit toolbar
    *.py                  # Python modules
  src/core/               # C++ source
  src/ext/
    HACApK_LH-Cimplm/    # H-matrix library (MIT)
    sparsesolv/           # Compact AMS/COCR (built into radia wheel, exposed as radia.sparsesolv_ngsolve)
  packages/
    cubit-mesh-export/    # Independent PyPI package (pip install cubit-mesh-export)
      src/cubit_mesh_export/
        check.py          # check-vol CLI + check_consistency() API
        cubit_mesh_curver.pyd  # C++ pybind11 module (bundled)
    radia-mcp/            # Independent PyPI package (pip install radia-mcp)
    radia-optuna/         # Standalone MATLAB Optuna (pip install radia-optuna)
    eqnedit64/            # Windows equation editor + Python API (pip install eqnedit64)
  tests/                  # Radia tests + tests/mcp/
  examples/               # retired; do not add new files
  docs/
  Build.ps1               # MSVC + MKL build
```

**PyPI packages** (5 independent distributions in the same monorepo):

| Package | Install | Purpose |
|---------|---------|---------|
| **radia** | `pip install radia` | C++ core + Python (HDiv-VIM/PEEC, panels, MCP) |
| **cubit-mesh-export** | `pip install cubit-mesh-export` | High-order curved mesh export from Cubit (does NOT require radia) |
| **radia-mcp** | `pip install radia-mcp` | MCP servers + skills for AI-assisted workflows |
| **radia-optuna** | `pip install radia-optuna` | Standalone MATLAB Optuna namespace + lightweight `optuna_mex`; no Radia solver/NGSolve/MKL dependency |
| **eqnedit64** | `pip install eqnedit64` | Windows equation editor, structural Python API, native TeX rendering/clipboard backend, and Web assets |

**Installation**:
```bash
pip install radia               # Python package and native solver kernels
pip install radia[cubit]        # Also installs cubit-mesh-export
pip install radia-optuna        # MATLAB Optuna only
pip install radia[optuna]       # Radia + independently versioned, validated radia-optuna
pip install eqnedit64           # Equation editor Python/API installation (Windows x64)
cubit-plugin-install            # Deploy Cubit .ccm backend + PySide6 toolbar (skip if no Cubit)
```

`radia-optuna` is the explicit distribution exception to the rule that Radia
numerical solver kernels stay in the main `radia` wheel. Optuna is a generic,
independently useful MATLAB compatibility component rather than a Radia physical
solver. Its canonical optimization sources remain under `matlab/+radia/+optuna`;
the separate wheel stages that tree plus the 21-command `optuna_mex`. The audited
generic Simulink subset (`buildOptunaBlock`, its Level-2 MATLAB S-Function/runtime
store, and `addOptunaMonitor`) is staged from `matlab/+radia/+simulink` too. It must
not acquire Radia-core, NGSolve, or MKL dependencies. Three named LTspice/sheet-
metal adapters may require Radia and must remain declared in the checked
distribution manifest. Generic `SimulinkRunner` and block operation belong to the
standalone package and must be tested from an installed wheel without the
repository or Radia on the MATLAB path. The block must persist trial tables and
expose progress/failure telemetry as Simulink signals; `optuna_mex` is required,
not an optional missing-MEX fallback.

**Deleted repositories** (integrated into Radia):
- ~~ksugahar/mcp-server-cae-ai~~ → `src/radia/mcp_server/`
- ~~ksugahar/ngsolve-sparsesolv~~ → `src/ext/sparsesolv/` (source only, build is separate)

---

## Repository Layering & Development Policy

Radia is a layered CAE stack.  When choosing where new work belongs, preserve
these boundaries instead of turning the repository into one generic solver or
one generic GUI.

| Layer | Canonical role |
|-------|----------------|
| Differential geometry | de Rham-complex concepts via NGSolve / Mathematica, with Radia collocation formulations kept compatible with that language. |
| Analysis methods | FEM / BEM through NGSolve and ngsolve.bem; Radia-owned magnetic moment, multipole moment, PEEC, and source-provider methods. |
| Linear algebra | Reuse proven solvers and compression tools such as shifted ICCG, AMS, ACA / TSVD, HACApK, BiCGSTAB, and BDDC. |
| Physical methods | Primary development focus: ESIM / SIBC, reduced potentials, Kelvin transforms, CLN, stream functions, and related open-region physics. |
| Application examples | Induction heating, MagLev, electromagnets, printed circuit boards, motors, and similar concrete engineering workflows. |
| Interfaces | Simulink application blocks for human production, Python/MCP for AI, result-bearing docs notebooks, and the Cubit toolbar. |

**Policy**:
- Follow NGSolve's Python API design wherever possible.  If NGSolve already has
  the right abstraction, extend around it instead of inventing a parallel Radia
  vocabulary.
- Expose Radia C++ functionality to Python with pybind11 and keep the Python
  surface idiomatic for NGSolve / NumPy users.
- Do not reinvent wheels.  Use public ecosystem components for CAD, FEM, BEM,
  meshing, visualization, and linear algebra whenever they are fit for purpose;
  implement only the missing Radia-specific electromagnetic / multiphysics
  pieces.
- Focus repository effort on the physical-method layer rather than on generic
  infrastructure.
- For an important model, prefer two or more independent analysis routes when
  feasible, so cross-validation is possible without relying on one formulation.
- Human production interfaces belong to application-specific Simulink blocks.
  Reusable capability stays in C++/Python/MATLAB APIs, headless CLI tools,
  result-bearing notebooks, validation tests, and MCP servers.
- MCP servers are part of the development loop: they encode executable
  knowledge, support autonomous validation / self-learning, and should reflect
  the same layer boundaries as the source tree.

### Active Two-Track Program (2026-09-01)

**POLICY**: Until this policy is explicitly superseded, active Radia delivery is
limited to these two outcomes:

1. Complete **HDiv-MMM topology optimization** as a usable, reproducible
   engineering capability.
2. Make the GitHub repository and **radia-mcp** a high-quality public software
   project.

The HDiv-MMM topology-optimization definition of done is an end-to-end,
named engineering design workflow, not merely an optimizer or a field-solve
prototype.  It must have one canonical solver/API route; constrained design
variables and objectives; reproducible input, output, and provenance artifacts;
an independent evaluation of the selected design; a recorded initial-versus-final
comparison; focused regression tests; and an executable public Python/MCP entry
point with documentation that can reproduce the evidence.

The GitHub/MCP quality definition of done includes a discoverable public API and
README, a single authoritative MCP catalog, documented tool contracts, focused
unit and real-protocol integration coverage, deterministic release/version
metadata, and a repository that contains neither obsolete installation paths nor
machine-specific execution dependencies.  Keep fast gates in `tests/`; place
expensive numerical evidence in `validation_test/` with checked provenance.

New application domains, exploratory solvers, UI work, packaging variants, and
unrelated refactors are deferred.  A change outside the two tracks is permitted
only when it fixes a release-blocking defect or directly unblocks one of their
acceptance checks.  Record other ideas in the backlog or memory rather than
expanding the public surface or CI matrix.

### Thin MCP Adapter and Coarse Workflow Policy (2026-09-02)

**POLICY**: An MCP server is a thin discovery and orchestration adapter over
canonical package APIs and durable artifacts. It must not become a second
implementation of the solver, a numerical-results database, or a flat catalog
of every test helper.

- Keep directly registered MCP tools for primary user workflows, status,
  topics, and coarse domain operations. Families of three or more one-purpose
  validation, identity, and fixture gates stay ordinary importable/testable
  Python functions and are exposed through the server's searchable
  `<server>_validation_catalog` plus `<server>_validation_run` pair. A smaller
  family remains direct because replacing one operation with two discovery
  tools would enlarge the public surface.
- The default `core` profile is the production public surface. The `full`
  profile, selected with `RADIA_MCP_TOOL_PROFILE=full` or
  `--tool-profile full`, temporarily restores legacy individual validation
  tools for migration, debugging, and release comparison; new clients must not
  depend on that flat surface.
- Server import must not eagerly load RAG/Chroma, optional numerical stacks, or
  validation-only modules. Resolve those dependencies on the first operation
  that needs them. Keep hot reload effective by resolving a lazy target's
  current module attribute at call time rather than retaining an old function
  object.
- MCP knowledge explains ownership, selection, workflow, and artifact
  interpretation. Numerical values and benchmark truth remain in versioned
  `validation_test/` JSON artifacts; MCP reads or points to those artifacts
  instead of copying their values into prose.
- Treat tool count, serialized `tools/list` size, cold import time, and first
  useful call latency as validation metrics. Record performance studies under
  `validation_test/`, while `tests/` retain fast contracts for profile
  selection, primary-tool visibility, lazy loading, catalog dispatch, status
  metadata, and real MCP transport.

### NGSolve-Native Discretization Policy (2026-07-16)

**POLICY**: Let NGSolve own finite-element plumbing.  Radia supplies the
electromagnetic method, open-boundary operator, constitutive model, and coupling;
NGSolve remains the source of truth for element orientation, local/global DOF
transforms, Piola maps, curved geometry, quadrature, weak-form assembly, and
`CoefficientFunction` / `GridFunction` evaluation.

- Build FE loads and couplings with NGSolve spaces, forms, grid functions, and
  mapped evaluation APIs.  Do not reconstruct high-order physical basis values
  from `CalcShape` and `GetDofNrs` in Python when NGSolve can evaluate or assemble
  the same quantity.
- This is especially strict for HDiv/HCurl HEX, WEDGE, curved, and high-order
  elements: their orientation and element transformations include NGSolve-owned
  local DOF transforms that are not a public Python reimplementation contract.
- Hand-derived reference-element algebra is allowed inside a tested Radia C++
  kernel when it is the new physical method itself.  Its Python boundary must
  still consume NGSolve meshes/spaces and be cross-checked against an independent
  NGSolve weak-form or `GridFunction` evaluation.
- Keep the caller-owned `with ngsolve.TaskManager():` convention for Python FE
  work; C++ kernels use the repository's TaskManager self-wrap policy.

### EarlyTimes Curvilinear Field-Map and Loft-Chain Policy (2026-08-15)

**POLICY**: EarlyTimes consumes an NGSolve HCurl finite-element A field on a
design-orbit-centred Bishop/RMF loft chain.  HDiv-MMM owns the continuous source
fields, NGSolve owns the conforming A projection, and the Lie-map and independent
Runge--Kutta routes remain separate enough to cross-check one another.

- The HDiv-MMM boundary supplies vector-potential and magnetic-flux-density
  `CoefficientFunction` objects.  Project A to an `HCurl` `GridFunction`
  before EarlyTimes Lie or canonical A-map tracking.  Independent Cartesian
  B-map Runge--Kutta may evaluate the original HDiv-MMM B
  `CoefficientFunction` directly; an `HDiv(order=4)` `GridFunction` remains
  an optional conforming-projection cross-check, not the source-field truth.
  The fourth-order Lie-map baseline is `HCurl(order=5)` for A.
- Treat measured median-plane B only as an HDiv-MMM inverse-design objective.
  Build observation rows at the physical probe locations and vary the pole
  topology/shape until complete three-dimensional field re-solves match the
  measured components within their bands.  Never use B-splines or a polynomial
  continuation of measured plane data to fabricate off-plane B or an A-map.
  Generate the HCurl A-map and independent B-map from the accepted physical
  magnet solution; unused measurement samples remain validation observations.
- The Lie route accepts only the constrained HCurl A-map and inserts A itself
  into the canonical Hamiltonian; it never substitutes `curl(A)`.  Canonical
  A-map Runge--Kutta evaluates `A_s,A_y` from that same HCurl field directly.
  Independent Cartesian B-map Runge--Kutta evaluates the HDiv-MMM B source
  directly.  Treat order 6 as a convergence check, not as the truth field or a
  substitute for the independent B-map route.
- Build each transverse section from four orthogonal quadrilateral strips split
  only in x, with symmetric nodes such as `[-a, -c, 0, c, a]`; each strip spans
  the full `[-b, b]` y interval.  Loft matching section nodes along the
  Bishop/RMF frame to form four HEX elements per longitudinal cell.  Thus the
  median plane `y=0` lies inside the elements, while the design orbit lies on the
  central `x=0` face.  On that face, `A_s` and `A_y` are tangential HCurl traces
  and therefore single-valued; `A_x` is normal and is removed by the gauge.
- Do not use a 2-by-2 quadrant split that places `y=0` on an internal face when
  the direct median-plane value of `A_y` is required.  In particular, never
  define `A_y` by naively averaging two HCurl normal traces.  If vertical
  refinement is needed, add symmetric y layers while keeping `y=0` inside a
  central layer.
- Use the full upper and lower volume in the source symmetrisation and conforming
  projection.  In a median-plane-symmetric gauge, enforce the reflection parity
  `A_x,A_s` even and `A_y` odd; correspondingly, `B_y` is even and `B_x,B_s` are
  odd.  Fit or extract `A_s` with even-y terms and `A_y` with odd-y terms so that
  off-plane derivative information is retained.  This parity-aware projection,
  not an upper/lower trace average, is the accuracy-improving operation.
- Apply the local gauge before projection: require `A_x=0` and
  `A_s=A_y=0` on the design orbit, using a gauge scalar even in y so the symmetry
  parity is preserved.  For a complete closed orbit, explicitly check the
  longitudinal circulation obstruction before claiming that a single-valued
  global gauge can set `A_s=0` around the entire ring.
- Keep the section axes `(e_x,e_y,t)` orthonormal and use
  `R(x,y,s)=r_0(s)+x e_x(s)+y e_y(s)`.  Use curved geometry of sufficient order
  for the order-5 field and refine the s subdivision until both geometry and
  field-map observables converge; do not select the longitudinal count from a
  fixed rule of thumb alone.
- Certify the usable order-5 Lie-map aperture against the independent HDiv-MMM
  B `CoefficientFunction` Cartesian Runge--Kutta result, with an HDiv
  GridFunction as an optional projection check.  Report separately the Lie
  truncation error
  (Lie versus unexpanded canonical A-map RK), the A/B field-route discrepancy,
  and the total Lie-versus-B-map discrepancy.

### MATLAB and Simulink Production Interface Policy (2026-07-20)

**POLICY**: Use MathWorks' official MATLAB MCP Server as the MATLAB execution
foundation and the Simulink Agentic Toolkit for generic model operations.
Radia-specific MATLAB, LTspice, NGSolve, optimization, and Simulink behavior is
the domain layer above those official foundations.  The final human-facing
production interface for Radia applications is a masked block in the single
**Radia** Simulink Library Browser entry, not a Jupyter notebook panel.

| Layer | Responsibility |
|-------|----------------|
| MATLAB MCP Server | Start/connect MATLAB, evaluate code, run `.m` files and tests, detect toolboxes, and run Code Analyzer checks. |
| Simulink Agentic Toolkit | Read, edit, check, test, and query generic Simulink models and Model-Based Design settings. |
| Radia C++ / Python / MATLAB APIs | Own numerical methods, NGSolve integration, artifact schemas, and the headless computation contract. |
| Radia Simulink blocks | Provide the final application-specific human operating surface through masks, typed ports, initialization, diagnostics, and result artifacts. |
| Python / MCP | Provide the first-class AI operating surface over the same numerical and artifact contracts. |

**Rules**:
- Keep this entire policy section synchronized verbatim with the corresponding
  section in the other agent-policy file in the same change.
- A tracked production `.slx` MUST be created, structurally edited, checked,
  and saved through MathWorks' official Simulink Agentic Toolkit. The required
  sequence is `model_read` before editing, `model_edit` for model operations,
  `model_check` after editing, then save and reopen the exact tracked path.
  Direct ZIP/XML patching is forbidden. A MATLAB `new_system` / `add_block` /
  `save_system` builder remains useful for temporary reconstruction and
  regression tests, but its raw output MUST NOT be promoted directly to the
  tracked production `.slx` without the official-agent open/check/save lane.
- A production-model acceptance check MUST start by closing every scratch,
  harness, and deleted-on-disk test model. Reopen the exact tracked `.slx` in a
  clean model state, verify its resolved `FileName`, and exercise the same
  open path a user gets by double-clicking the file. Tests and agent trials that
  create a temporary model MUST close it with `onCleanup` / teardown before
  deleting its file, including error and missing-toolbox paths. Never leave an
  orphaned in-memory model visible to the user.
- Simulink visual QA MUST inspect the complete application window, not only a
  cropped model canvas. The title bar, left and right dock strips, status bar,
  block labels, masks, and confirmation/error dialogs are all in scope. Public
  model-authored UI text is English. Any replacement glyph, square-box text,
  mojibake, or suspicious `???` is a release-blocking failure: discard or
  regenerate the affected artifact through the official agent; do not blame
  Qt, patch fonts, or repair the `.slx` archive by hand. Also scan embedded SLX
  XML for U+FFFD and suspicious question-mark runs before acceptance.
- Publish application blocks through one **Radia** library. The production
  application grouping covers Electromagnet, PCB/PEEC, Motor, Stream Function,
  and Induction Heating. Additional validated MEX/ROM, LTspice, optimization,
  or application subsystems may coexist without changing this contract.
- All Jupyter notebook panels and workbenches, including the former `radia-ih`
  comparison workbench, are retired as production interfaces. Remove their
  launchers and adapters; do not retain notebook-only production features or
  create new notebook panels. Every Radia application uses its masked Simulink
  block as the final human operating surface.
- `docs/**/*.ipynb` is the canonical public example, derivation, reproduction,
  and validation layer, not a production GUI. Every repository-published CAE
  example MUST be an executed, result-bearing notebook with narrative, code,
  saved output, and saved `ngsolve.webgui.Draw` scenes. An adjacent JSON is not
  required. Draw the geometry
  or mesh and the primary computed field when one exists; use
  `netgen.webgui.Draw` for pre-mesh CAD geometry. A source catalog, migration
  archive, static PNG, or script excerpt alone does not qualify as an example.
- Once a historical migration is complete, remove its source-catalog or
  migration-ledger notebook and JSON archive. Git history is the archive;
  `docs/` keeps maintained public notebooks, not migration bookkeeping.
- A notebook field scene MUST pass the field and mesh explicitly and name the
  view, for example `Draw(field, mesh, name="B_magnitude", ...)`. Keep useful
  display choices such as volume/surface visibility, autoscaling, clipping,
  vectors, and widget size as explicit arguments. A bare `Draw(field)` is not
  the public field-visualization contract.
- Public example notebooks set `metadata.radia.notebook_role="example"` and
  `metadata.radia.webgui_required=true`. The document-meta notebook audit must
  fail an example that lacks an executed `Draw` cell with saved rich output.
  When a primary computed field exists, also set
  `metadata.radia.webgui_field_required=true`; the audit then requires an
  executed parameterized field scene with saved rich output.
- A tracked native `.slx` interface sample may accompany the library when it has a
  canonical `.m` generator, a direct load/update regression test, and an
  explicit backend declaration. It demonstrates block composition and signal
  wiring; it does not replace the result-bearing docs notebook that owns the
  numerical evidence. A production IH sample must use readable Level-2 MATLAB
  Eddy/Thermal S-Functions backed by native object-handle MEX kernels and must
  not be a LUT or lumped-state-space surrogate.
- Blocks must delegate to tested `radia.*` APIs, MEX/ROM handles, or validated
  headless application entry points. They must not reimplement solver logic.
- Every operation exposed through Radia's Python API MUST also be possible from
  MATLAB with the same numerical meaning, state transition, error behavior, and
  artifact contract. MATLAB may reach the shared C++ kernel through a standalone
  MEX ABI or use an explicit Python fallback at a non-step-time boundary, but it
  must not provide a reduced or behaviorally different substitute. When users
  configure, generate, inject, or observe time-domain waveforms, prefer a
  Level-2 MATLAB S-Function so the signal flow, sampling, state lifecycle, and
  diagnostics remain visible and editable in Simulink.
- Radia follows the same high-level-orchestration / narrow-native-kernel design
  used by NGSolve's Python / pybind11 boundary. The default Simulink execution
  layer is a readable Level-2 MATLAB S-Function; Level-1 S-Functions are
  forbidden. Its `setup`, ports, sample time, state, Outputs/Update separation,
  termination, and diagnostics remain inspectable MATLAB code. Performance-
  critical numerical work is exposed through independently callable MEX
  function ABIs operating on numeric arrays, sparse matrices, structs, files,
  and checked `uint64` handles. This M-file S-Function plus standalone-MEX
  composition is the production default for LTspice, NGSolve, Radia solvers,
  optimization, and other application blocks, not merely a prototype fallback.
  Promote the complete wrapper to a Level-2 C/C++ MEX S-Function only when the
  application requires hard/soft real-time determinism, deployable generated
  code, zero-copy native state, native resource ownership tied directly to the
  Simulink engine, or reproducible end-to-end measurements demonstrate a
  material benefit that cannot be obtained through the standalone MEX ABI.
  IH Eddy/Thermal follows this default through separate native object handles.
  Every block declares its backend and passes lifecycle and numerical-parity tests;
  MATLAB and C++ implementations must never silently substitute for one another.
- When conversion or reconstruction of a native solver object would materially
  affect performance, use an explicit handle ABI: `create` returns an opaque,
  checked `uint64` handle; `step` / `evaluate` / `update` accept that handle;
  and `destroy` releases it. The native registry MUST validate type, generation,
  ownership, and liveness on every call so stale or cross-kernel handles fail
  loudly. A MATLAB handle class or Level-2 MATLAB S-Function owns the token and
  guarantees cleanup through `delete` / `Terminate` / `onCleanup`; never expose
  a raw pointer as the public ABI or rely only on `mexLock` for lifetime safety.
- IH production uses Level-2 MATLAB S-Functions for ports and lifecycle, with
  independent native C/C++ MEX object handles for Eddy and Thermal numerical
  state. Eddy receives current, angle, and temperature distribution and emits
  heat distribution; Thermal receives that heat distribution, ambient
  temperature, and angle and emits the accepted temperature distribution.
  Python fallback is not permitted per time step.
- The initial application-block backend may be the validated Python/headless
  CLI, launched only on an explicit trigger or update. MEX/ROM is an optional
  later backend promotion, not a prerequisite for the Simulink interface. Move
  a path to MEX only after numerical parity, error propagation, lifecycle, and
  long-run stability are independently tested; keep the block contract stable.
- Expensive CAD, mesh, basis construction, and field solves occur during model
  initialization, explicit update commands, or artifact generation. Per-step
  simulation uses checked native MEX/ROM state where available and must not
  spawn Python once per time step.
- Each production block must have a mask parameter contract, typed and
  documented ports, sample-time semantics, fail-fast dependency checks,
  `run.log` / `result.json` compatible provenance, a MATLAB test entry point,
  and a numerical cross-check against the corresponding headless golden.
- Every solver-bound `.vol` MUST pass `check-vol` after mesh export and before
  solver or Simulink initialization. A production application/mode owns a
  versioned `radia.vol-label-contract.v1` file and enables strict labels; keep
  the `cubit-mesh-export.vol-check.v1` JSON report in the run directory and
  index it from `result.json`. A sibling `.vol.json` CAD reference is used when
  present but is not required for the structural/label gate. The checker does
  not infer material constants from labels; DesignSpec/configuration validates
  conductivity, permeability, BH data, frequency, and other physical values.
- Every Simulink application run that computes a spatial field MUST write a
  checked GMSH `.msh v4.1` post-processing artifact inside its run directory
  and list it in `result.json`. The runner owns the output path; this is not a
  user mesh-generation setting. Scalar/circuit-only modes with no spatial
  field must declare the artifact not applicable and must not fabricate one.
  GMSH remains a visualization/post-processing target, never the solver mesh
  generator or the NGSolve interchange format.
- `cubit-mesh-export` is the repository's only supported VTK producer. All
  Radia solver, validation, docs, and MCP spatial-field output uses checked
  GMSH `.msh v4.1`. Outside the Cubit export component, direct
  `ngsolve.VTKOutput` calls and custom `.vtk`, `.vtu`, or `.vts` writers are
  forbidden. Negative tests and documentation of external formats may mention
  VTK, but they must not provide a Radia-owned VTK output path. NGSolve's
  upstream `VTKOutput` remains a valid NGSolve API and may be documented as
  such; Radia workflows choose GMSH as their post-processing contract.
- Spatial GMSH figures, including contours, flux lines, streamlines, and
  sections, use a physical 1:1:1 axis scale by default. A deliberately
  exaggerated axis is allowed only when the scale factor is explicit in the
  render configuration and stated in the figure caption; silent distortion is
  not a valid production visualization.
- Treat `radia-mcp` as the canonical executable manual for Radia-specific
  MATLAB and Simulink workflows. Do not duplicate generic MathWorks guidance.
- Long or solver-heavy MATLAB validation runs execute on hibino first. Use mdx
  only when hibino is unavailable and the mdx CI queue is idle; LAB and 100号機
  remain development and fast-test hosts.

### Python-to-MATLAB Capability Parity and Fallback Policy (2026-07-21)

**POLICY**: Every user-facing capability available through Radia's Python API,
including pybind11-backed kernels, Python numerical modules, and NGSolve
integrations, MUST have a named MATLAB entry point. Native C++/MEX is the
preferred implementation and MUST be attempted on a best-effort basis. When a
complete native port is not practical, the MATLAB entry point MUST remain
available through an explicit Python fallback; unsupported behavior must never
be silently omitted or replaced by a numerically different MATLAB algorithm.

**Rules**:
- Keep C++ as the single numerical source of truth wherever a capability can be
  expressed with numeric arrays, sparse matrices, files, or checked native
  handles. Expose that source through thin pybind11 and MEX adapters. Do not
  maintain independent Python and MATLAB copies of the same numerical kernel.
- Every production `.py` module under `src/radia/` that exposes user-facing
  behavior MUST be represented by a corresponding `.m` entry point under
  `matlab/+radia/`, or by an explicit entry in the checked parity manifest. The
  manifest classifies each module as `native-mex`, `matlab-native`,
  `python-fallback`, or `private/not-applicable` and names the owning `.m` file.
  A generic command runner by itself is not a module-level MATLAB counterpart.
- Prefer MEX for pybind11-backed methods, repeated numerical kernels, large
  array operations, and Simulink step-time code. Translate Python objects into
  MATLAB-friendly numeric arrays, structs, sparse matrices, file/path-based
  meshes, and checked `uint64` handles rather than trying to share Python object
  identity.
- A standalone MEX function ABI is the standard native deliverable. It provides
  the numerical kernel used by the Level-2 MATLAB S-Function, a low-level MATLAB
  entry point, a reproducible debugging surface, and a boundary that can be
  tested independently of Simulink. Do not bury a reusable numerical kernel
  exclusively inside a MEX S-Function.
- Standalone MEX and S-Function layers have separate release responsibilities.
  Standalone MEX must pass API, numerical, error-propagation, and performance
  tests; production Level-2 MATLAB S-Functions must pass initialization,
  step-time, termination, repeated-run, and Simulink integration tests. When an
  approved native MEX S-Function exception exists, it must pass both categories
  plus its generated-code, real-time, zero-copy, or native-lifecycle requirement.
- Use an explicit Python fallback when the capability fundamentally depends on
  Python callbacks, Python-only ecosystem objects, or an NGSolve operation that
  has no stable native C++ boundary. The owning `.m` function may use MATLAB's
  `py.*` / `pyrun*` interface and therefore depend on the configured
  `python312.dll`. It MUST check `pyenv`, fail fast with installation guidance,
  convert inputs and outputs through documented numeric/struct/file contracts,
  and expose which backend actually ran.
- A Python fallback is allowed for initialization, explicit update, artifact
  generation, and batch solves. It MUST NOT be invoked once per Simulink time
  step. Per-step execution requires a checked MEX/ROM/state-space path.
- A new public Python feature is not production-complete until its `.m` entry
  point, parity-manifest entry, documentation, and focused MATLAB regression
  test exist. Intentional limitations must be recorded in the manifest and
  exercised by a fail-loud test.
- Backend promotion requires an accuracy check and a reproducible performance
  comparison on identical inputs. Measure cold startup separately from warmed
  steady-state execution, report repeated median timings and data-transfer
  costs, and record the machine/runtime versions. Neither pybind11 nor MEX is
  assumed faster a priori; select the production backend from measured end-to-
  end latency, throughput, memory behavior, and operational reliability.
- Keep the official MATLAB MCP Server as the execution foundation; this policy
  governs Radia's domain API and fallback behavior above that execution layer.
- Treat these weak MATLAB families as the active native-promotion backlog, in
  the checked order recorded by `matlab/python_api_parity_manifest.json`:
  axifem; high-level VIM / ESIM / IH; high-level BEM / PEEC / SIBC;
  Kelvin / DtN; acoustic CQ-BEM / FSI; Motor / MagLev; and coil CAD.
  Each family entry MUST name its current native surface, next stable numeric or
  artifact boundary, deliberately retained Python/ecosystem boundary, and test
  gates. Landing one focused MEX command advances the family but does not mark
  the whole family complete.

### CI Execution, Validation Evidence, and Notebook Policy (2026-09-03)

**POLICY**: **mdx** is Radia's self-hosted CI and preflight host. LAB and
100号機 are development machines and MUST NOT execute pull-request, push,
scheduled, or release CI jobs. A pre-push check sends the unpushed candidate
commit to mdx and runs the same fast contract lane there; it is not a
best-effort substitute executed on the developer desktop.

mdx gives CI and preflight work priority over ad-hoc compute. Long-running
optimization, parameter sweeps, solver validation, and benchmark jobs use
hibino first when it is available. They may use mdx only when hibino is
unavailable and the mdx CI queue is idle; a compute job must not delay, starve,
or destabilize the self-hosted CI runner.

CI scope begins at the independently released distribution boundary. Keep a
commit focused enough that `radia`, `cubit-mesh-export`, `radia-mcp`,
`radia-optuna`, and `eqnedit64` changes can be identified without rebuilding or
retesting the other distributions. Pull-request and main-push workflows run
only the owning distribution lanes plus genuinely shared repository contracts.
Within `radia-mcp`, classify changes as docs-only, packaging/version metadata,
or implementation/tests. Docs-only changes run relevant document checks;
generated-inventory checks run when the inventory or its generator changes;
metadata-only changes build and import-check one wheel; implementation/tests
run the supported-Python compile/import matrix for the affected package
families. On the newest supported Python, normal pull-request and main-push CI
runs a stable compact contract set, tests selected from the changed source/test
paths, and only the affected server selftests. A shared `common` change selects
all server selftests. The complete package pytest suite, all-server selftests,
and full generated-inventory audit run only through the explicit full-audit
workflow or a named pre-tag release-candidate audit on the same commit. A
release-tag workflow builds and verifies the exact tagged artifact without
repeating the full tests already accepted for that source commit. Full
numerical, GUI, performance, and multi-machine evidence remains an explicit
validation or release gate rather than an automatic response to every commit.
Because AI agents are expected to diagnose and repair failures interactively,
normal CI optimizes for fast, high-signal feedback rather than exhaustive proof;
strict breadth belongs to explicit validation and release acceptance. Normal CI
does not automatically rerun a failed deterministic test; it returns the first
failure promptly for diagnosis.
Developer pre-push hooks run only the impact-scoped mdx preflight. They MUST
NOT upload native binaries or mutate GitHub Releases; exact native artifacts
are published by release workflows from the accepted CI run.

The repository has three complementary, non-duplicative evidence surfaces:

- `tests/` contains deterministic unit and contract tests that fit the fast CI
  budget. A test belongs here only when it is self-contained and protects a
  focused public or source-level invariant against a concrete implementation
  bug. Do not keep two tests whose purpose and failure signal are the same.
- `validation_test/` contains solver-heavy numerical evidence, long native
  builds, GUI/Simulink checks, golden comparisons, performance measurements,
  and multi-machine studies. When normal CI exposes a defect or uncertainty,
  run only the relevant validation lane and retain its result JSON. Validation
  artifacts require machine-readable JSON with the checked values, runtime,
  versions, and execution host. Validation runs explicitly on hibino first, or
  on mdx only when hibino is unavailable and the mdx CI queue is idle, through
  a scheduled validation lane or named release gate; it does not run on every
  pull request.
- `docs/**/*.ipynb` is the public calculation record: derivation, executed
  evidence, figures, and a presentation-ready narrative. It is not a retired
  workbench. A notebook should be promotable directly into a talk or paper
  figure sequence without reimplementing the calculation elsewhere. The saved
  notebook output is sufficient; an adjacent JSON and a runtime gate are not
  required.

MCP and notebooks serve the same artifact contract from different directions:
MCP is the executable, AI-facing discovery/orchestration surface, while the
notebook is the human-readable, result-bearing record of that execution.
Neither may carry a private numerical implementation that the other cannot
reproduce through the shared Python/C++ API and durable artifacts. A docs-only
contract lane parses changed notebooks and rejects malformed JSON or
replacement characters; a named validation lane re-executes the notebook sets
named by a changed method or release gate on the permitted compute host.

Every new test declares its tier from measured runtime and dependency scope.
When a formerly fast test grows beyond the CI budget, move its evidence and
artifacts to `validation_test/` while retaining a focused contract in `tests/`.
Every new public notebook example retains executed outputs and the required
parameterized WebGUI scene; a talk reuses those saved results rather than
creating a second calculation path.

### MATLAB Optuna Upstream Differential-Oracle Policy (2026-08-21)

**POLICY**: Upstream Optuna 4.9.0 is the sole behavioral oracle for every
MATLAB Optuna test that exercises behavior shared with Optuna. A MATLAB
implementation detail, historical Radia result, or handwritten expected value
must never define compatibility. MATLAB-only storage, Simulink, native-MEX, and
parallel-execution checks may remain as integration tests, but they are not
evidence of Optuna parity.

- Generate shared-behavior expectations by executing pinned
  `optuna==4.9.0`. Use the same explicit sampler seed, sampler options, search
  space and parameter order, trial numbers, completed/pruned/failed history,
  constraints, and objective values on both sides. Run proposal-sequence
  comparisons sequentially so scheduling cannot change random consumption.
- Treat the upstream Optuna algorithm as the common algorithmic core, even when
  MATLAB or MEX executes an independent native implementation. Preserve the
  upstream equations, transforms, state updates, boundary handling, and seeded
  random-consumption order; do not silently replace them with a different
  optimizer or approximation. Native vectorization, MEX kernels, batched work,
  parallel trial execution, table/MAT persistence, and Simulink telemetry may
  improve performance around that core. A deliberate MATLAB-only algorithmic
  extension must use a distinct option or type, be documented as an extension,
  and must not change the upstream-compatible default path.
- Record and check every backend version that can affect the result: Python,
  NumPy, SciPy, PyTorch, and `cmaes` where applicable. Regenerators must fail
  when the pinned Optuna version is not present and must produce byte-stable
  JSON on repeated runs.
- Use direct upstream Optuna execution for seeded numerical and state-machine
  comparisons. Use the official `optuna/optuna-mcp` server over a real MCP
  transport for the public MCP Study/Trial tool contract. The MCP server is not
  a seeded numerical oracle while its `set_sampler` tool does not expose a
  seed.
- radia-mcp supports the MATLAB distribution without becoming a second Optuna
  MCP server. Every shared operation present in the official
  `optuna/optuna-mcp` live `tools/list` MUST be executed by that upstream
  server. `mcp-server-radia-matlab` owns only MATLAB differences: table/MAT
  persistence, Simulink monitoring and failure telemetry, MATLAB parallel
  execution, the standalone `optuna_mex`, and Radia CAE artifact adapters. It
  MUST NOT proxy, rename, or reimplement an upstream Optuna MCP tool.
- Preserve the applicable Optuna and `optuna-mcp` MIT copyright and permission
  notices whenever upstream source or a substantial portion is copied or
  redistributed, and record copied/adapted provenance. `radia-optuna` is an
  independent, unofficial project; it is not affiliated with, sponsored by,
  or endorsed by Preferred Networks, Inc. or the Optuna project. Public docs
  MUST include exactly: "Optuna, the Optuna logo and any related marks are
  trademarks of Preferred Networks, Inc." Do not use the Optuna logo or imply
  official status. Differential MCP regeneration uses local stdio and a fresh
  per-run temporary SQLite database, never shared/production storage; routine
  tests use checked fixtures, never launch Dashboard, and never automatically
  create upstream issues or pull requests.
- A test of shared behavior must read its expected values, states, ordering,
  warnings/errors, or proposal sequence from an upstream-generated fixture.
  Do not add handcrafted sampler sequences, quality bands, dominance orders,
  pruning decisions, default values, or private-state assumptions as the
  compatibility truth.
- Prefer public observable behavior. An internal MATLAB test is allowed only
  when it protects a MATLAB-specific invariant that upstream Optuna cannot
  express; classify it explicitly as `matlab-integration` and state the
  boundary. Such a test must not be cited in docs, changelogs, or releases as
  evidence of upstream parity.
- Every `tests/matlab/test_optuna*.m` test function must appear in the checked
  Optuna oracle manifest as `upstream-python`, `upstream-mcp`, or
  `matlab-integration`. Shared-behavior tests classified as
  `matlab-integration` are a policy failure.
- Generate the complete pinned public-API inventory, including exported symbols
  and public class members, directly from `optuna==4.9.0`. Compare it with the
  MATLAB surface in a checked coverage file. A complete-compatibility claim is
  forbidden until every required entry is present and exactly oracle-mapped;
  partial family evidence is not closure. MATLAB-only parallel execution and
  table/MAT storage are extensions and do not waive any shared API behavior.
- Treat upstream `seed=None` as nondeterministic constructor behavior, not an
  exact proposal sequence. Verify that the upstream default is `None`, that
  separately constructed MATLAB samplers draw fresh private entropy, and that
  this does not mutate MATLAB's global RNG. Exact sequence parity still uses an
  equal explicit seed on both sides.
- When upstream and MATLAB disagree, first add or regenerate the upstream
  fixture, then fix MATLAB. Never update the expected result from the MATLAB
  output. A deliberate unsupported upstream feature must fail loudly and be
  recorded as a limitation, not normalized into a passing alternative.

---

## Critical Policies

### Self-Driving Loop Discipline (2026-06-24)

Drive an autonomous / self-paced loop by **task completion, not a clock**: finish
one verified step, then immediately start the next — do not insert idle interval
ticks. Every iteration must produce **concrete, verified progress** (a number that
came from running the code, never an estimate), recorded to the project's internal
validation notes. **Never overclaim** — a "pass" must trace to a checked value. Stop
and ask only where the decision is genuinely the user's.

### Validation-Class Examples (2026-06-24)

Do not force long-running or solver-heavy validation problems into `tests/`.
Keep CI tests small and fast enough to catch regressions, but still run heavier
research validation when it teaches something new.  New development trials and
scratch scripts live under `C:\temp`, not in the tracked source tree.  Once a
trial becomes worth keeping, promote it to `validation_test/<topic>/`, a
result-bearing `docs/<topic>/*.ipynb`, a notebook-local `docs/<topic>/*.py`
helper, or `src/` according to the classification below. A masked Simulink
block is the final human operating surface promoted from a mature headless
contract, not from a notebook workbench. A validation-class problem may be too slow for pytest and still
be mandatory learning material.

When a validation-class example is being promoted, the executable verification
lane is the repository's actual `validation_test/` directory, not `docs/`
(and not a differently named `tests_validation/` tree).  A result-bearing
`docs/<topic>/*.ipynb` may render the theory,
tables, and plots for humans, but it does not replace the runnable validation
surface.  If historical material has `validation_*.py`, `validate_*.py`,
`*_summary.json`, or references from `validation_test/`, classify it first as
`validation_test` / protected-validation-corpus material; add a docs notebook
only as an executed showcase layer when it helps readers.

Heavy `validation_test/` runs are compute-host-only. Use LAB/100号機 for fast
import/path smoke checks, small correctness probes, and build checks, but run
solver-heavy sweeps, timing claims, memory/scaling runs, and research-grade
validation on hibino first. Use mdx only when hibino is unavailable and the mdx
CI queue is idle. Confirm the selected host is idle before launch, and record
the hostname and runtime in the result JSON/log. LAB/100号機 timings are smoke
observations only and must not be presented as benchmark data in docs, MCP, or
papers.

### Retired Examples / Promotion Triage (2026-06-28, updated 2026-07-04)

**POLICY**: `examples/` is retired and must not be recreated.  It is neither a
scratch area nor a teaching tier.  New development experiments run in `C:\temp`;
tracked work enters the repository only after promotion into one of the durable
lanes below.  Historical `examples/` references are migration blockers.

| Class | Destination | Rule |
|-------|-------------|------|
| Development-in-progress / superseded / failed iteration | keep in `C:\temp` until distilled, then delete | Preserve the lesson in `memory/<topic>.md` or a short docs note when it matters. Git history and `C:\temp` are the scratch/archive path, not `examples/`. |
| Reusable computation, parser, mesh reader, solver helper, formula, or API surface | `src/` | Promote to a named public or internal API and add focused tests. Do not keep it as a loose example helper. |
| Fast implementation regression or minimal fixture | `tests/` | Keep small enough for CI / developer feedback. |
| Important numerical verification, benchmark, golden lock, convergence sweep, or regression corpus | `validation_test/<topic>/` plus optional docs notebook | The executable check and required result JSON live in `validation_test/`; docs may render selected theory, tables, and plots for humans. |
| User-facing explanation, tutorial, or method showcase | `docs/<topic>/*.ipynb` | Notebook must be executed and result-bearing. Integrate Markdown explanation, code, saved results, and WebGUI where applicable; no adjacent JSON is required. |
| Notebook-only subroutine / local renderer / catalog helper | `docs/<topic>/*.py` | Allowed only when tightly coupled to the notebook. If another topic, panel, MCP server, or validation uses it, promote to `src/` instead. |
| Mesh/CAD/journal/result assets | keep until owning script/notebook is migrated | Mesh definitions, Cubit `.jou`, tracked `.msh`, figures, and JSON results are protected by preservation/reproducibility policy. |

The migration order is strict: inventory and reference search first; create or
refresh the executed docs notebook if the result is user-facing; move reusable or
validation code to `src/` or `validation_test/`; update docs/MCP/panel
references; then delete the historical source.  Never leave two live copies of
the same implementation after API promotion is complete, and never add new
long-lived references to `examples/`.

`protected_*` / "保護参照あり" is a temporary blocker, not a destination.
If a docs notebook, validation test, panel sample, MCP knowledge file, or
README still references `examples/<topic>`, record the blocker and the
`target_after_unblock` (`docs`, `src`, `validation_test`, or distill-delete),
then migrate the reference.  Public docs may refer to other `docs/` artifacts,
and code may refer to `src/` APIs, but new long-lived references to
`examples/` should not be introduced.

### Publish Boundary: No Validation Provenance in Public Artifacts (2026-06-24)

Public artifacts (this repo, PyPI packages, public docs) **lead with analytic
solutions**. Do NOT put into them: internal absolute paths or local working-directory
names; a third-party tool's **benchmark numbers** used as a validation basis; or
"verified / validated against <external tool>" attributions, including the names of
that tool's source files. Citing a *published* convention or a peer-reviewed paper is
fine — attributing validation to an internal or third-party reference is not. Keep
cross-validation provenance in machine-local notes only; stored regression-reference
values may remain, but **unattributed**.

### No Development Cruft in SOURCE — Distill the Lesson to memory/ (2026-06-26)

**POLICY**: Development-in-progress iterations MUST NOT accumulate in the
tracked SOURCE tree. The SOURCE keeps only the **final / canonical** version of
a given piece of code; superseded snapshots, abandoned formulations, and debug
stepping-stones are removed once they are superseded. The **lesson** (so the
same rut is not re-walked) is distilled into the memory system — a
`memory/<topic>.md` file plus a one-line `MEMORY.md` index entry — NOT left as
dead code. Pruning a dead branch and recording *why* it died grows the
repository's knowledge; leaving the dead branch in the tree does not.

**Concrete rules**:

1. **Canonical-only in SOURCE.** "Canonical" = the single version reachable
   from the current committed tree that the build / tests / goldens exercise.
   Superseded iteration snapshots identified by ad-hoc suffixes
   (`_v2`/`_v3`/`_v9`/`_old`/`_new`/`_corrected`/`_fixed`/`_tmp`/`_backup`/`_wip`
   or duplicate-with-numeral filenames) are NOT canonical and do not belong in
   the tracked tree once the canonical version lands.

2. **Distill, then delete — in that order.** When retiring a superseded
   iteration, FIRST write the lesson to `memory/<topic>.md` and add a one-line
   `MEMORY.md` index entry (per "Promotion-After-Verify Policy"; keep the index
   line under ~200 chars — `MEMORY.md` is already over its length limit, so push
   detail into the topic file). THEN remove the dead snapshot. The knowledge
   must survive the deletion. Recover the old code from git history if ever
   needed — it is not lost, only un-tracked.

3. **Removal targets dead implementation code only.** In scope: renamed
   snapshots (`rad_*_v2.cpp`, `calc_*_old.py`), abandoned-formulation source,
   and debug stepping-stone modules that no live path imports. A negative result
   that warrants a minimal reproducing test/fixture keeps that test — prune the
   dead snapshot, not the test that documents why the approach failed.

4. **Genuinely co-valid alternatives are NOT iteration snapshots.** Two
   implementations that are BOTH live and user-selectable (e.g.
   `--coil-solver peec|bem-a`, or independent FEM/BEM validation paths) are separate
   canonical artifacts. Do not delete one as "superseded" — neither superseded
   the other.

**Exceptions** (this policy prunes superseded ITERATION CODE only — it NEVER
deletes protected data, assets, or fixtures):

- **Golden test fixtures** under `tests/**/fixtures/` and golden-band lock files
  are protected (per the promotion ladder). Distinct tier artifacts of one
  geometry (a `tests/` fixture vs a `docs/` result notebook vs a
  `src/radia/panels/samples/` `.jou`) are SEPARATE canonical artifacts, not duplicate
  snapshots — keep all.
- A **single, explicitly-named frozen reference baseline** (e.g. a
  `src/radia/panels/samples/*` loft baseline kept deliberately for regression comparison)
  is a protected canonical artifact, not cruft — keep it even though a newer
  variant exists.
- **Committed figure/table/published-result-backing data** — `.json`/`.csv`/
  `.png`/`.pdf` files committed next to their script or figure to make a
  committed figure/table/published result reproducible — are protected results,
  NEVER iteration snapshots. Removing such data would make the figure
  non-regenerable; do not delete it under this policy.
- **Tracked mesh definitions and mesh-gen assets** — mesh files
  (`.bdf`/`.nas`/`.msh`, plus Cubit-export `.vtk` fixtures), Cubit `.jou` journals, and
  **mesh-generation scripts** — are protected by "Mesh File Preservation" and
  are NEVER deleted by this policy, even when a `_v2`/`_old` name makes them look
  like an iteration. Historical tracked meshes formerly under `examples/` should
  be migrated to their owning `docs/`, `validation_test/`, or
  `src/radia/panels/` lane. A
  file that is BOTH a `_v2`-named snapshot AND a protected asset is governed by
  the preservation policy, not this one.
- **LAB-local, non-committed authoring references** (e.g. a `.nb` kept beside its
  canonical `.wls`) are out of scope — this policy governs the TRACKED source
  tree only.

**Why**: dead snapshot code rots — it gets imported by accident, copied as a
template, or mistaken for the live path, producing wrong numbers nobody can
trace (the same failure class as silent fallbacks). The lesson, by contrast,
compounds when it lives in `memory/` where the next session reads it before
re-deriving the dead end. A clean tree where every tracked file is the canonical
one, with the rationale for each removed branch preserved in memory, is the
repository actually being tended.

### Green's Function: Laplace Kernel Only (MQS/Darwin)

**POLICY**: Radia uses **Laplace kernel only**: $G(r) = 1/(4\pi r)$. Target regime is MQS (Magneto-Quasi-Static) to Darwin approximation.

**Do NOT**:
- Add Helmholtz kernel ($e^{-jkr}/r$) to any Green's function
- Use wave number $k$ in field calculations (except for skin depth)
- Implement full-wave EFIE or MFIE formulations

Skin depth is computed from frequency for SIBC, but field propagation uses quasi-static approximation.

**Affected Components**: `rad_green_fullwave.h/cpp`, `rad_conductor.cpp` (`GreenFunction()`), `rad_hacapk.cpp`.

### Matrix Storage: Row-Major (C-style)

**POLICY**: All interaction matrices use **row-major [target][source] format**.
- `A[i][j]` stored at `i * stride + j`; represents effect ON target i FROM source j
- All BLAS calls use `CblasRowMajor`
- Python interface returns NumPy C-contiguous (row-major) arrays

**Source Files**: `rad_interaction.cpp`, `rad_relaxation_methods.cpp`, `rad_hacapk.cpp`.

### Binary File Policy

**POLICY**: Native build outputs (`.pyd`, `.dll`, `.so`, `.lib`, `.exe`,
`.mex*`) are release artifacts, not source files. They are built in isolated
CI/release environments and delivered through wheels or versioned GitHub
Release packages. Do not upload binaries from a pre-push hook and do not add
source-tree download/bootstrap scripts for them.

- The reviewed `cubit_mesh_export.ccm` distributed by the independent
  `cubit-mesh-export` package is the explicit tracked-plugin exception.
- A release workflow consumes the artifact built and accepted for that exact
  commit; development machines must not copy native outputs directly to mdx.
- Tracked `.slx`, result-bearing notebook output, `.png`, and `.pdf` are allowed
  when required by their documented interface or publication role.
- Generated solver meshes and field outputs (`.msh`, `.vol`)
  stay out of source control unless a focused fixture policy explicitly owns
  them.

### GitHub Release Publication Gate (2026-07-23)

**POLICY**: A Radia Simulink library release, including
`radia_simulink_library.slx`, MATLAB support files, and MEX assets, MUST NOT
be published to GitHub Releases until the complete `release-quad` four-machine
verification has passed. The four machines are LAB, 100号機, mdx, and hibino.

- The release candidate is assembled and tested before publication.
- The `release-quad` `done` gate is the authoritative publication decision.
- A failed, partial, or manually waived machine check is not a release pass.
- GitHub Release assets must include the versioned Simulink package,
  `manifest.json`, and `SHA256SUMS.txt` when the candidate contains them.
- The same gate applies to later revisions of the library, not only the first
  Simulink publication.

### File Placement Policy

**POLICY**: Development scratch outputs belong in `C:\temp`.  Committed output
files (`.png`, `.msh`, `.vol`, validation JSON) must be placed next to
their owning `tests/`, `validation_test/`, `docs/`, or `src/radia/panels/`
driver.
- Do NOT place generated files at the repository root
- Build output goes to `build*/` or `dist/` (both gitignored)
- `docs/<topic>/` MAY contain `.py` helpers imported by that topic's
  result-bearing notebook. Reusable behavior belongs in a `src/` API.
- Every `docs/<topic>/*.ipynb` method/showcase notebook must be result-saving:
  execute it before committing so code-cell outputs, figures, and tables are
  embedded. It is a public demonstration, not benchmark evidence, and does not
  require an adjacent JSON or runtime threshold. Numerical truth, timing,
  convergence, and publication evidence belong in `validation_test/<topic>/`
  with a machine-readable result JSON that records versions, host, and runtime.
- Every public CAE example notebook must also save an interactive WebGUI scene:
  use `ngsolve.webgui.Draw` for the generated mesh and primary
  `CoefficientFunction`/`GridFunction`, or `netgen.webgui.Draw` for CAD geometry
  before meshing. Execute the `Draw` cell and preserve its rich output. Mark the
  notebook metadata with `radia.notebook_role="example"` and
  `radia.webgui_required=true`; the notebook audit enforces this contract.
- `src/radia/panels/` owns headless `calc_*.py` application entry points,
  samples, and artifact schemas. It contains no application notebook
  workbench and is not the final human UI location. Do not create a repo-root
  `panels/` tree.
- `matlab/+radia/+simulink/` owns Simulink builders and runtime adapters;
  `matlab/radia_simulink_library.slx` is the distributable human-facing block
  library. Native block kernels and build sources may live under
  `src/radia/simulink/`.

### Promotion Ladder: C:\temp → tests / validation_test / docs / src / Simulink (2026-07-20)

**POLICY**: New exploratory scripts start outside the repository in `C:\temp`.
`examples/` is retired forever.  A file enters the source tree only when it has
a durable role:

| Lane | Purpose (intent) | Audience | Ships in wheel? |
|------|------------------|----------|-----------------|
| `tests/**` | **実装バグの検出** — small deterministic regression, fixture, and API/error contract. Do not duplicate an existing test purpose. | CI / Codex / developer | No |
| `validation_test/<topic>/` | **重要な検証・ベンチ・golden lock** — numerical truth with required result JSON; large runs use hibino first and mdx only behind an idle CI queue. | developer / agent / research validation | No |
| `docs/<topic>/*.ipynb` | **ユーザーに理論と代表結果を同時に見せる** — executed, output-bearing demonstration; no JSON sidecar required. | users / collaborators / future agents | Docs |
| `docs/<topic>/*.py` | Notebook-local helper only. | notebook readers / MCP if local | Docs |
| `src/` | Reusable API, parser, formula, solver helper, and computation kernel. | package users / validation / MCP / blocks | Yes |
| `src/radia/panels/` | Validated headless application CLI, `DesignSpec`, samples, and artifact contracts. | blocks / AI / validation | Yes |
| `matlab/+radia/+simulink/` + `matlab/radia_simulink_library.slx` | Final application-specific human operating surface. | end users | MATLAB distribution |

### Test Runtime Placement Policy (2026-08-28)

**POLICY**: `tests/` is the short default CI/debug suite; `validation_test/`
(the repository's canonical singular directory name) owns long-running and
environment-heavy checks. Do not keep a long test under `tests/` and hide it
with `@pytest.mark.slow` or a node-ID side list.

- A test that is at least 10 seconds in two comparable successful CI runs is
  moved to `validation_test/<topic>/` and marked `slow`. Re-run a one-off timing
  spike before reclassifying it.
- Real Office/Cubit GUI startup, licensed external applications, solver
  convergence studies, long golden/reference runs, benchmarks, and publication
  validations belong to `validation_test/` regardless of a lucky short run.
- Fast API/shape/error contracts remain in `tests/`. If a module mixes fast and
  slow coverage, keep reusable setup as a non-collected helper and expose the
  long assertion through a test collected from `validation_test/`.
- `validation_test/slow_nodeids.txt` may classify older measured validation
  tests. It must contain only `validation_test/` node IDs; there is no
  corresponding slow-node list under `tests/`.

### Independent Validation Oracle Policy (2026-09-03)

**POLICY**: An analytical solution, manufactured solution, symbolic derivation,
or high-accuracy reference whose purpose is to judge a production solver stays
outside the production C++ kernels. Do not port a validation-only oracle to C++
merely to make validation faster.

- Keep validation oracles readable and independent under `validation_test/`, or
  in a tracked Python, MATLAB, or Mathematica reference used from that lane.
  Solver-heavy validation runs on hibino first, with mdx allowed only when
  hibino is unavailable and the mdx CI queue is idle; oracle speed alone is not
  a reason to reduce implementation independence.
- An oracle must not call the same production function, pybind11/MEX binding,
  or native kernel whose numerical behavior it certifies. Do not expose a
  validation-only analytical reference through the production pybind11 or MEX
  surface.
- A closed-form or analytical algorithm may still be implemented in C++ when
  it is a genuine shipped capability used by applications at production scale.
  That C++ implementation is then part of the product under test, not an
  independent reference; retain a separate readable oracle for validation.
- Store validation results and provenance in the owning `validation_test/`
  JSON. Record the reference implementation and runtime so a passing result is
  not mistaken for an implementation self-comparison.

### Maintenance Stabilization Window (2026-09-03 to 2026-10-03)

**POLICY**: During this window, do not add new user-facing features. Spend
repository effort on completing, reviewing, simplifying, documenting, and
validating capabilities that already exist.

This section has precedence during the window. It overrides standing native-
promotion, parity-expansion, application-promotion, and development-priority
guidance whenever following that guidance would create a new supported surface.
The end date does not automatically resume feature development: review the
maintenance results and explicitly set the next policy first.

- Allowed work includes bug fixes, API and policy consistency, removal of
  retired paths, dependency and build cleanup, CI/runtime reduction, test and
  validation classification, documentation repair, and completion of an
  already-claimed capability whose current implementation is incomplete.
- Do not introduce a new solver family, application, public API family, MCP
  tool family, Simulink product block, or optional dependency during the
  window. Record promising ideas for later instead of implementing them now.
- Security fixes, upstream compatibility fixes, and the smallest change needed
  to prevent data loss or restore a broken released workflow are exceptions.
  Keep any exception narrowly scoped and state why it could not wait.
- Maintenance changes must reduce or leave unchanged the long-term support
  surface. A rewrite that merely moves unfinished behavior behind a new API is
  a feature change, not maintenance.

**Scope boundary**:

- A new solver or formulation, application, public Python/MATLAB/MEX API
  family, MCP tool family, Simulink product block, artifact schema, or runtime
  dependency is a new feature and is forbidden during the window.
- Filling a recorded gap in an already-public and already-claimed contract may
  be maintenance when it adds no new numerical meaning or workflow. Examples
  include a missing MEX mapping, parity-manifest classification, packaging
  metadata, error propagation, or a broken released entry point.
- Refactoring for reliability or measured performance is maintenance only when
  the public contract and artifact meaning stay unchanged. Do not use the
  window to redesign a working API under a new name.
- Reference formulas used only as numerical oracles stay in readable Python,
  MATLAB, or validation code. Do not promote analytical validation references
  into C++, pybind11, or MEX merely to make them look production-ready.

**Maintenance baseline and active backlog**:

- Completed baseline: canonical `radia-optuna` source selection is
  deterministic; the public accelerator-field, lamination, and MMM-topology
  MATLAB entries are classified; and `ObjTetrahedronCurrent` plus the Kelvin
  field/potential surfaces have explicit standalone-MEX coverage or boundary
  classifications. HDiv-MMM environment controls are finite and exhaustively
  classified, and native stats record the effective release-relevant path.
- Active work: continue CI, build, dependency, and repository cleanup. Remove
  duplicate tests, retired interfaces, and stale generated files while keeping
  fast regression coverage for implementation bugs and JSON-backed numerical
  evidence in `validation_test/`.
- Treat any regression in the completed baseline as a maintenance bug, not as
  a reason to reopen an alternative source-root, fallback, or duplicate API.

Work on one maintenance family per commit. Rebase it onto the latest clean
`main`, run its focused fast tests, and report unrelated global failures
separately. Never make a gate pass by deleting an independent oracle, weakening
a tolerance without evidence, skipping the affected lane, or reading expected
results from the implementation under test.

**Promotion gates**:

- **C:\temp → tests/**: the behavior is small, deterministic, and useful for
  fast regression.
- **C:\temp → validation_test/**: the run is a numerical validation,
  benchmark, convergence sweep, golden lock, or regression corpus; heavy runs
  use hibino first, or mdx only when hibino is unavailable and the mdx CI queue
  is idle, and are labelled with the actual validation host.
- **C:\temp → docs/**: the result teaches a method or workflow to humans; the
  notebook must be executed, output-bearing, and Markdown-integrated. It has no
  mandatory JSON sidecar. A public CAE example also includes saved WebGUI
  `Draw` scenes for its geometry/mesh and primary field.
- **C:\temp → src/**: the code is reusable by more than one documentation,
  validation, application-block, or MCP path.
- **Mature docs/headless contract → Simulink**: publish a masked block only
  after the headless command and golden are stable. The block delegates to that
  contract and adds no numerical implementation.

Same geometry may exist in multiple lanes only when the artifacts have distinct
roles: e.g. a minimal fixture in `tests/`, a heavy truth run in
`validation_test/`, a human-facing result notebook in `docs/`, and a packaged
headless sample in `src/radia/panels/`. No lane points back to `examples/`.

### Application Interface Promotion Policy (2026-07-20)

**POLICY**: Application interfaces are built in **four gated stages**. The
production human surface is a Simulink block; a documentation notebook is not
promoted into a notebook GUI.

**Stage 1 — Enumerate the app-specific variables.**
Write down every knob the user of this specific application might want
to change.  This is a list, not code.  Pin the solver-specific variables
(mesh size, frequency, material, source current, ...) and the
**solver-switch variable** itself (e.g. `--impedance-model linear|esim`,
`--solver pardiso|ams`).

**Stage 2 — CLI Python script (`src/radia/panels/calc_*.py`).**
Turn the Stage-1 list into an argparse-driven Python script.
Computation only, no GUI.  JSON on stdout.  The solver switch **must**
also be a CLI flag so the same script can drive any supported backend.
Map those arguments into a small UI-neutral `DesignSpec` dataclass so Python,
MCP, MATLAB, and Simulink share one settings-to-command contract.

Stage 2 is validated by running the application mode end-to-end against its
sample input and comparing the resulting scalar against a golden band in
`validation_test/panels/test_*_golden.py`, with the evidence recorded in JSON.
Two or more solver choices must be
exercised and produce numerically consistent results (within the
mode's documented tolerance).

Stage 2 is considered **合格 (pass)** when:
- `python calc_<mode>.py --help` exits 0 and prints all knobs
- running against the sample with **each** supported solver switch
  produces JSON whose key numbers are inside the golden band
- `validation_test/panels/test_<mode>_golden.py` locks the result and writes or
  checks its owning result JSON

**Stage 3 — result-bearing documentation notebook (`docs/<topic>/*.ipynb`).**
Explain the validated method, inputs, outputs, equations, and representative
result. Save outputs in the notebook; no adjacent JSON is required. For a CAE example,
include executed `ngsolve.webgui.Draw` cells for the mesh and primary result
field and save their rich output. The notebook may call the Stage-2 contract
but may not become an alternate production workbench.

**Stage 4 — masked Simulink application block.**
Publish the human interface in the single Radia library. The block exposes a
stable mask contract, typed ports, explicit sample-time/trigger semantics,
fail-fast dependency checks, and `run.log` / `result.json` artifacts. It
delegates to the Stage-2 contract. The initial backend may launch the validated
Python CLI once per explicit trigger; it must never launch Python every time
step. A later MEX/ROM backend keeps the same block contract and is promoted only
after independent parity and long-run tests.

Stage 4 is considered **合格 (pass)** when:
- the block is present in `matlab/radia_simulink_library.slx` with an
  application-specific mask and documented typed ports
- a MATLAB test loads, updates, and executes the block's supported path
- the same sample produces the same golden-band quantities as Stage 2
- success, timeout, dependency failure, and solver failure leave inspectable
  result/log diagnostics

No Radia application has a notebook-workbench stage. The Cubit Export Mesh
toolbar remains the separate allowed PySide6 surface inside Coreform Cubit's
embedded Python.

**Why**:
- Forces the hard thinking about *what is changeable* before any widget
  code is written (Stage 1 is where over-scoping is cheapest to cut).
- The solver switch being a Stage-2 argument keeps backend bugs visible to the
  same sample and golden.
- Stage 3 preserves readable theory and results without making Jupyter an
  operational dependency.
- Stage 4 gives human operation a stable Model-Based Design surface while
  Python/MCP remains the first-class AI surface.

Related:
- "Panel Samples Quality Policy" above — Stage 2's validation relies
  on trustworthy samples.
- "Cubit Panel Architecture" below (§ 4-Layer) — Simulink and Cubit remain
  separate interface processes over the same headless artifacts.

### Panel Visualization Routing Policy (2026-06-26)

**POLICY**: Panel-facing visualization follows the existing ecosystem
instead of custom Radia viewers.

- **GUI / notebook route (human-facing)**: use `netgen.webgui`
  (`Draw(...)`, browser/Jupyter scene widgets) for in-panel or
  notebook visual inspection.
- **LLM / headless route (automation-facing)**: export `.msh v4.1`
  and use `gmsh` / `GmshPostExport` for file-based inspection,
  screenshots, geometry/field artifacts, and scripted validation.
- **Notebook input IO files**: `.vol` and `.sol` are user-facing
  notebook inputs/outputs.  Double-click viewing should open the plain
  Netgen viewer, not a Radia-specific viewer or a GMSH desktop window.

Do not invert these roles. Notebook or panel UX should not depend on a
GMSH desktop window for ordinary interactive viewing, and LLM workflows
should not depend on a transient `netgen.webgui` browser state when a
durable `.msh`/`.json` artifact can be produced.

Windows note: the pip-installed Netgen launcher is `netgen.exe`, backed
by the `netgen.__main__` module.  As of this memo, `netgen.__main__`
ships file handlers for `.py`, `.geo`, `.step`, and `.stl`, but not
`.vol`/`.sol`.  If `.vol` double-click does not open, the fix belongs
in Netgen's Python package startup (`netgen/__main__.py`, with any DLL
setup still handled by `netgen/__init__.py`), adding a `.vol` handler
using Netgen's native loader (`Ng_LoadMesh`) and then associating
`.vol`/`.sol` with `netgen.exe "%1"`.  Do not treat the old
`radia-vol-viewer --register` association as the default notebook IO
route; keep it only as a legacy/helper path when a custom `.sol`
companion-mesh heuristic is explicitly needed.

### Verify-First Policy: FES inspection before physics solve (2026-04-25)

**POLICY**: When debugging FEM setup (Periodic BC, Dirichlet, material
labels, mesh.Curve, Kelvin identification, etc.), check the **finite
element space** first -- BEFORE running any physics solve.  FES checks
are sub-second; a botched Kelvin/Periodic setup wasted on a 5-minute
solve is unacceptable.

The minimum verify trio:

1. **Materials / boundaries spelled as expected** -- one-line check:

   ```python
   print(mesh.GetMaterials(), mesh.GetBoundaries())
   ```

2. **Periodic / Dirichlet actually constrains DOFs** -- compare
   `H1(...)` vs the constrained variant (`Periodic`, `dirichlet=...`):

   ```python
   fb = H1(mesh, order=p, dirichlet="GND")
   fp = Periodic(fb)
   slaved = sum(fb.FreeDofs()) - sum(fp.FreeDofs())   # must be > 0
   ```

3. **Functional boundary test for Periodic** -- set 1.0 on slave bnd,
   integrate on master bnd, ratio must be 1.0:

   ```python
   gfu = GridFunction(fp); gfu.vec[:] = 0
   gfu.Set(1.0, definedon=mesh.Boundaries("kelvin_int"))
   r = Integrate(gfu*gfu, mesh, definedon=mesh.Boundaries("kelvin_ext")) \
       / Integrate(gfu*gfu, mesh, definedon=mesh.Boundaries("kelvin_int"))
   ```

If any of (1)-(3) fails, the issue is in geometry / Identify / labels --
fix that BEFORE solving.  Do NOT iterate solver runs to discover an
FES-level bug; reading the wrong number from a wrong solve is the
fastest way to mis-diagnose the next layer.

**Real example (2026-04-25)**: commit 3297b5c9 reverted the EM panel
default to `fes_order=1` based on a 5-minute physics-solve discrepancy
("p=2 gives -138 mT vs ELF -228 mT"), claiming the C++ Kelvin pair
identification only handled linear vertices.  An FES inspection done
directly on `em_elf_quarter.vol` showed `slaved=6085` at p=2 with
`Set(1)|kelvin_int -> kelvin_ext ratio=1.0`: the Kelvin BC was already
correct at p=2.  The actual cause of the discrepancy is a separate
Periodic-Omega-reduced single-space p-stability issue.  The expensive
solve cost ~10 minutes; the FES check costs <1 second.

### AI-Driven Cubit: Probe, Don't Guess (2026-04-25)

**POLICY**: When AI is authoring Cubit Python that identifies entities
(volumes / surfaces / curves / vertices) by geometric properties, it
MUST first **probe Cubit** (`cubit.parse_cubit_list("...", "all")` +
per-entity centroid / bbox / area) and **print the actual values**
before writing the classification logic.  Do not derive a centroid
filter "by hand" from the .jou source -- Cubit's webcut, subtract,
imprint, and merge operations create surprising intermediate volumes,
renumber entities, and (with `subtract ... keep`) leave artifacts.

**Why**: AI does not have the geometry in its head.  An educated-guess
centroid bound (e.g. "the kelvin outer cap has cx > kelvin_offset, the
cut face has cx == kelvin_offset") is wrong as often as it is right --
e.g. for an offset 1/4-sphere octant the outer cap's x-centroid is
ALSO at the sphere center because of y-z mirror symmetry of the 1/4
patch.  Same lesson at copy-mesh anchor selection: 1/8 spherical caps
have 3 equal-length boundary arcs, so `max(curves, key=length)` ties
non-deterministically and Cubit's listing order can pick different
arcs on source vs target -- see
`memory/feedback_kelvin_1_8_blocker.md` and the deterministic
`min(curves, key=(centroid_z, y, x))` fix in
`_add_kelvin_cubit_reduction`.  Running with assertion failures and patching by inspection is
slow and noisy; running ONE probe pass first is fast and final.

**Pattern**:

```python
# Step 1: probe -- print everything we might filter on.
for vid in cubit.parse_cubit_list("volume", "all"):
    v = cubit.volume(vid)
    c = v.centroid()
    bb = v.bounding_box()
    print(f"vol {vid}: c=({c[0]:.3f},{c[1]:.3f},{c[2]:.3f}), "
          f"extent=({bb[3]:.3f},{bb[4]:.3f},{bb[5]:.3f}), "
          f"vol={v.volume():.3e}")
for sid in cubit.parse_cubit_list("surface", "all"):
    s = cubit.surface(sid)
    cx, cy, cz = s.center_point()         # NOT centroid() -- Surface API
    bb = s.bounding_box()
    print(f"surf {sid}: c=({cx:.3f},{cy:.3f},{cz:.3f}), "
          f"area={s.area():.3e}, extent=({bb[3]:.3f},{bb[4]:.3f},{bb[5]:.3f})")

# Step 2: classify based on observed numbers.
mag_vol = next(v for v in volumes if cubit.volume(v).volume() < 1e-4)
# ...
```

**Specifics**:
- `cubit.volume(vid).centroid()` exists (returns 3-tuple).
- `cubit.surface(sid).center_point()` is the equivalent on Surface
  (Cubit does NOT expose `.centroid()` on Surface).
- `cubit.volume(vid).bounding_box()[3:6]` and `cubit.surface(sid).bounding_box()[3:6]`
  are the (extent_x, extent_y, extent_z) tuple -- a flat cut face
  has zero extent in its cut direction.
- `cubit.volume(vid).volume()` returns the measured volume; for a
  1/8 octant of radius R it is `(4*pi*R^3/3) / 8`.
- After `subtract A from B keep`, expect the SUBTRAHEND (`A`, the
  thing being removed) to remain at its original ID, the MINUEND
  (`B`) to be modified to (B - A), AND a duplicate of B to also
  appear -- delete the duplicate explicitly via centroid + volume
  inspection (do not assume it does not exist).
- After `webcut volume all with plane <plane> offset 0`, all volumes
  intersecting the plane are split; volumes entirely on one side are
  NOT split.  Do NOT assume "N volumes -> 2N volumes".

This POLICY tightens "Journal File Portability Policy" (no hardcoded
IDs, identify by geometric properties) for AI-authored Cubit scripts:
the geometric predicate must be derived from a printed probe, not
from an a-priori derivation.  Apply when:
- A Cubit Python script fails an assertion on entity identification.
- The script uses webcut, subtract, intersect, sweep, or imprint+merge.
- The next-step user is the AI itself (i.e. you are about to debug
  your own Cubit script).

### Promotion-After-Verify Policy (2026-04-25)

**POLICY**: When a debug session ends with a verified result, propagate
the knowledge to ALL three layers BEFORE moving on:

1. **`memory/<topic>.md` + `MEMORY.md` index entry** -- the lesson is
   useless if the next conversation re-discovers it from scratch.
2. **Owning durable artifact** -- record the validated result where it now
   lives: `validation_test/<topic>/` JSON/README for executable truth, or an
   executed result-bearing `docs/<topic>/*.ipynb` for
   user-facing explanation, or MCP knowledge for agent-operational rules.
   Include the specific checked value (e.g. "VERIFIED p=2: slaved=8914
   DOFs, ratio=1.0") so future contributors know the result is golden,
   not "should work".
3. **`AGENTS.md`** if the lesson is a method (e.g. "always check FES
   first") that applies to future debugging.

A debug session that produced a fix but did not propagate is
**incomplete**.  The repository must grow knowledge along with code, or
the next contributor (including future-you) will repeat the same
multi-hour investigation.

### MCP Knowledge Placement Policy (2026-04-21, updated 2026-04-24)

**POLICY**: All MCP knowledge ships from the single Radia monorepo.
LAB-private mcp-server packages have been retired.

| Where | What |
|-------|------|
| `S:\Radia\01_GitHub\packages\radia-mcp\src\radia_mcp\<topic>\` (**public**) | All MCP knowledge — general FEM/BEM/Kelvin **and** application-specific (induction heating, electromagnet, peec). Ships to PyPI as `radia-mcp`. Each topic is a subpackage with its own `server.py`. |
| ~~`S:\mcp-server\mcp-server-ih\`~~ (**retired 2026-04-24**) | Promoted into `radia_mcp.ih` subpackage. Old path no longer exists. |

**Subpackage layout** (one server per concern):

| Subpackage | Scope | Promoted from |
|------------|-------|---------------|
| `radia_mcp.cubit` | Cubit scripting / scaffolding / hex-mesh | n/a |
| `radia_mcp.build123d` | build123d STEP authoring | n/a |
| `radia_mcp.gmsh` | GMSH MSH v4.1 inspect/validate/convert | n/a |
| `radia_mcp.elf` | Historical ELF reference | n/a |
| `radia_mcp.interop` | Cross-CAD interop | n/a |
| `radia_mcp.radia_ngsolve` | **General** Radia + NGSolve (FEM/BEM/Kelvin/PEEC inductance/sparsesolv/MSH post) | n/a |
| `radia_mcp.ih` | **IH-specific**: induction heating workflow, ESIM cell problem, workpiece SIBC, Karl iteration, screening physics | `s:\mcp-server\mcp-server-ih\` (2026-04-24) |

**Splitting rule** (general vs application):
- If the topic is generally useful for FEM/BEM/Kelvin/.vol pipeline — `radia_ngsolve`
- If the topic only makes sense in a specific application context (induction heating, accelerator magnets, PCB) — that application's subpackage

**New research topics**: in-flight WIP lives in `C:\temp` until stable.
Promotion into a `radia_mcp.<topic>` subpackage requires: feature committed,
promoted to its durable lane (`src/`, `validation_test/`, `docs/`, or
`src/radia/panels/`), deploy-verified, golden-tested, and knowledge stops referencing
unpublished scratch files. There is **no longer** a separate
`S:\mcp-server\mcp-server-*\` tree — promote directly into the
public subpackage when ready.

**Past examples** (historical):
- PEEC-inductance: `mcp-server-ih` → `radia_mcp.radia_ngsolve.peec_inductance_knowledge` (2026-04-21, general technique)
- Induction Heating: `mcp-server-ih` → `radia_mcp.ih` (2026-04-24, application-specific subpackage)
- Analytical Formulas: in `radia_mcp.radia_ngsolve.analytical_formulas` (2026-05-01, extended same day) — closed-form reference layer (Wakao-Igarashi-Fujiwara-Kameari Part 1-9). Group B+C (radia 4.20.0): ellipsoid demag/torque, AC vector locus, magnetic shielding, 2D rectangular bar, thin-plate eddy current, Fabri solenoid, three-phase line, K(k)/E(k) Hastings, Gauss-Legendre. Group D (radia 4.21.0, Part 6/8/9): plate Joule dissipation, AC thin-shell shielding, magnetic-shell interior fields, planar surface impedance, full Bessel cylindrical-conductor AC impedance, Gauss-Patterson nested quadrature, cuboid average B (numerical-integration path). **radia 4.22.0**: cuboid_average_field closed-form C++ kernel shipped — sympy-derived G1, G2 antiderivatives + 64-corner inclusion-exclusion sum (~40 µs/call, 817× faster than Gauss-Legendre baseline); `method="numerical"` fallback retained for cross-checks and the V_T ≪ V_S ULP-cancellation regime. The originally-cited Stafl 1967 §3.4 was confirmed unrelated (2D rectangular conductor, not 3D cuboid magnetisation). MCP tool `analytical_formulas(topic)` exposes 11 topics including a `validation_use_cases` mapping that says "given analysis X, which closed form is the trusted reference?". Use it as the FIRST QUESTION when validating any new analysis result.

### Axisymmetric FE: Henrotte for Magnetic, Standard H1 for Scalar (FEMM-Canonical)

**POLICY (2026-05-10, refined from earlier "all axisym Henrotte"
framing)**: Axisymmetric FE convention follows the FEMM 4.2 split:

| Physics | Basis | Reason |
|---------|-------|--------|
| **Magnetic A_phi (curl-curl)** | **Henrotte** `{1, r^2, z}` (`radia.radia_axifem`) | The cylindrical curl operator `B_z = (1/r) d(r A_phi)/dr` produces a `1/r` integrand that standard FE Gauss quadrature cannot integrate accurately near the axis.  Henrotte's `s = r^2` substitution gives clean closed-form integration. |
| **Scalar T / phi (Laplacian)** | **Standard NGSolve `H1`** + `2 pi r` weighting | The weak form `int k grad T . grad v . 2 pi r dr dz` has `2 pi r` as a **smooth Jacobian** (not a `1/r` integrand).  Standard FE handles this fine; no axis-special treatment is needed. |

This follows the documented FEMM 4.2 axisymmetric convention:

- **magnetic** uses the Henrotte `{1, r^2, z}` basis
- **heat** uses a **standard P1 triangle** with `2 pi r` evaluated at the
  element centroid — **no** Henrotte basis, **no** `s = r^2` substitution

**Why not "all axisym Henrotte"**: Henrotte basis IS the natural
function space for axisymmetric scalars (the parity / even-function
argument is mathematically correct).  But the practical accuracy
benefit is small for scalar Laplacians because the `2 pi r` Jacobian
suppresses spurious odd-r modes automatically.  FEMM ships
production-grade thermal accuracy with standard P1; we follow that
proven convention.

**API for magnetic axisym**:

```python
import radia.radia_axifem as ax

mesh = Mesh(...)                                  # axis-aligned (r, z) mesh
fes  = ax.H1Henrotte(mesh, order=p)               # p = 1 (Q1) or p = 2 (Q2)

a_mag = BilinearForm(fes, symmetric=True)
a_mag += ax.AxiHenrotteStiffnessBFI(mu_cf)
a_mag += ax.AxiHenrotteSigmaMassBFI(sigma_cf)     # eddy-current term
```

**API for scalar axisym (heat, electric potential, diffusion)**:

```python
from ngsolve import H1, BilinearForm, x as r_coord, dx, ds, grad, InnerProduct
import math

mesh = Mesh(...)
fes  = H1(mesh, order=p)                          # standard NGSolve H1

weight = 2 * math.pi * r_coord                    # axisym Jacobian
a_heat = BilinearForm(fes, symmetric=True)
a_heat += k_cf * InnerProduct(grad(u), grad(v)) * weight * dx
a_heat += h_conv * v * u * weight * ds(surface_label)   # Robin
```

**Optional Henrotte heat infrastructure**: The
`radia.radia_axifem.AxiHenrotteHeat{Stiffness,Mass}BFI` classes
(added in radia 4.31.0) and the `H1Henrotte` BND DiffOp (radia
4.32.0) are kept in the codebase as parity-conscious infrastructure
for research / publication uses (e.g. comparing convergence rates of
Henrotte vs standard H1 on a scalar problem).  They are NOT used by
production heat solvers and are NOT required.

**Near-axis order rule for axisym SCALAR solves (2026-09-03 study)**:
standard P1/Q1 cannot reproduce a nonconstant even quadratic radial
profile while also representing `dT/dr = 0` at the `r = 0` axis.  The
near-axis profile therefore shows a piecewise-linear cusp.  Production
axisym heat (`calc_heat_axisym.py`) defaults to `--fes-order 2`; the 3D
thermal default remains order 1.  On the committed uniform-flux cylinder
validation, order 2 reduces the spurious radial axis slope from about
139 C/m to below 0.3 C/m.  This refines, not replaces, the standard-H1
convention above -- a Henrotte switch is still not needed for scalars.
Initialize a uniform order>=2 state with `gfT.Set(CF(T0))`, never by
filling its hierarchical coefficient vector.  Report extrema from
physical GridFunction evaluations over vertices plus deterministic
volume and boundary samples; neither raw coefficient extrema nor
vertex-only extrema are valid for a general order>=2 field.

**Reference**: see
[`docs/axifem/FORMULATION.md`](docs/axifem/FORMULATION.md)
sections 5-6 (Henrotte basis derivation for magnetic) and 10b/10c
(optional heat BFIs).

### Unit System Policy

**POLICY**: Radia always uses **meters**. There is no unit conversion in C++. All coordinates are in meters, all current densities in A/m^2.

**`FldUnits` is removed**: Do NOT call `rad.FldUnits()` in any code. Radia always uses meters with no configuration needed.

```python
# CORRECT
magnet = rad.ObjHexahedron(vertices, [0, 0, 954930])  # meters, A/m

# WRONG - hard-coded conversion
x_mm = x_m * 1000.0  # DO NOT DO THIS
```

**Radia Units** (always meters, no conversion):
- All coordinates in meters
- B in Tesla, H in A/m, A in T*m
- Current density J in A/m^2
- Physical constants in `rad_constants.h`: `MU_0_OVER_FOUR_PI = 1e-7`, `INV_FOUR_PI = 1/(4*pi)`

### Magnetization Units: A/m (NOT Tesla)

**POLICY**: Radia uses **M in A/m**. Common conversion: `M = Br / mu_0` (e.g., Br=1.2T -> M=954930 A/m).

Do NOT confuse M (A/m) with J (magnetic polarization, Tesla): J = mu_0 * M.

### Windows Console Encoding (cp932)

**POLICY**: NEVER use Unicode mathematical symbols in print statements. Use ASCII equivalents: `^2` not `²`, `->` not `→`, `<=` not `≤`, etc. Windows console defaults to cp932 in Japanese environments.

### Repository Language: English

**POLICY**: All source code, documentation (`docs/**/*.md`), comments, commit messages, and docstrings in the Radia repository MUST be written in **English**. Japanese text is NOT allowed in tracked files. Exception: `AGENTS.md` may contain Japanese policy descriptions. Conversation with the user may be in Japanese, but repository content must remain English-only.

### Naming Policy: External Project References

**POLICY**: Do NOT use "ELF" or "ELF_MAGIC" in Radia source code, documentation, or comments. Radia is an independent project. Academic citations are allowed.

### No Console Output from C++ Code

**POLICY**: No `printf`/`cout`/`cerr` in C++ code for logging. All user-facing output through Python. Allowed: error messages via `Send.ErrorMessage(...)` and `#ifdef DEBUG_...` guards.

### No Fallbacks — Fail Fast, Fail Loud

**POLICY**: Do NOT write fallback chains (`try API_A except: try API_B except: try API_C`). Pick the **one** API that works for the project's target environment (Cubit 2025.12, NGSolve 6.2.2606, Python 3.12) and commit to it. If the chosen API stops working, fix the call site or raise — never bury the breakage under another path.

**Design philosophy**: this is the **fail-fast** principle (Jim Shore, 2004) combined with **"errors should never pass silently"** (PEP 20, Zen of Python) and **"explicit is better than implicit"**. The deeper rationale is *user agency*: a fallback path performs a computation the user did not ask for and cannot inspect. Silent fallback violates the **Principle of Least Astonishment** — the user gets a number, has no way to know which code path produced it, and trusts it. That trust is unrecoverable when the result turns out to be from the wrong path.

**Erroring out is FAR better than silently producing a sloppy / wrong result.** Compare:

| Behavior | What user sees | What user can do |
|---|---|---|
| Raise with available labels | "label 'src' not found; available: [source, sink, ...]" | Fix the .jou (the actual root cause) in 30 seconds |
| Silent fuzzy match → wrong label | A plausible-looking number | Notice the bug 6 months later, after publishing |
| Silent default constant (e.g. `H_t_rms = 5.0`) | A plausible-looking number | Never notice |

The raise is a **feature**, not a defect.

**Why fallbacks are harmful**:
- They hide bugs. When the primary call silently fails, the fallback masks it; the next person sees "it works" and never learns the primary is broken.
- They obscure intent. Readers cannot tell which path is the supported one.
- They defeat root-cause debugging. The exception that would have pointed at the real problem gets swallowed.
- They violate user agency. The user neither requested nor can audit the fallback path.
- Wrong numbers from a fallback are worse than no numbers — they get put in papers and presentations.

Note: the **Robustness Principle** ("be liberal in what you accept") is now widely regarded as harmful for both protocols and scientific code. Prefer **strict in what you accept**.

**How to apply**:
- Cubit Python: use `get_block_id_list()` / `get_sideset_id_list()` directly. Do NOT add `parse_cubit_list("sideset", "all")` as a fallback.
- File loading: pick one format (e.g. `.vol`), not "try .vol, fall back to .cub5".
- Solver: pick one preconditioner (e.g. CompactAMS), not "try AMS, fall back to BoomerAMG".
- Surface mesh extraction: use the in-memory NGSolve extractor (`_extract_surface_mesh_filtered`), not "try in-memory, fall back to .vol-text rewrite".
- Material/boundary lookup: if the expected label is missing, raise with the available label list — do not guess by case-insensitive prefix matching, kind-of-name fallback, or numeric→string fallback.
- Numerical defaults: NEVER substitute a "reasonable-looking" constant (e.g. `H_t_rms = 5.0` if Karl iteration didn't converge). Raise.
- GND/source identification: if the labelled GND vertex doesn't exist, raise. Do NOT pick "the closest vertex to the origin" as a substitute — it changes silently when the mesh is regenerated.
- Periodic boundary identification: if the C++ identification path fails, raise with a clear message. Do NOT silently re-do the work in Python with different tolerances.

**Allowed (these are not fallbacks)**:
- A single `try/except` that converts a system error into a clean user message — but it must `raise` or `return {"error": ...}`, never silently try another path and continue.
- Two-path code where the user **explicitly chose** the path (e.g. `--mode bem` vs `--mode fem`). The user is in control.
- Optional features genuinely not needed for the result (e.g. "if matplotlib is installed, also save a PNG").
- Default *parameter values* declared at the function signature — these are the user's contract, not a runtime fallback.

### STEP-Only Centerline: Auto-Detect or Fail (2026-05-15)

**POLICY**: PEEC coil filament generation extracts the centerline
**exclusively from the STEP file** via
`coil_from_cad.extract_centerline_from_step`.  There is no
caller-supplied centerline override (`path_points_m` was removed in
v4.48.1) and no JSON ingestion path.  If auto-detection cannot
recover a centerline that covers the conductor's bounding box, the
panel / CLI **raises** -- the user fixes the CAD rather than papering
over the breakage with a hand-crafted polyline.

**Reasoning**: A user-supplied centerline JSON is unauditable -- the
user types numbers into a file, PEEC consumes them, and the only
sanity check is the user re-reading their own JSON.  Real-world
incidents (e.g. keiko 2026-05-15 `1turn_coil_loft_outsideline.step`)
show that the auto-detect failure mode is the geometry being
ambiguous -- not the auto-detect being wrong -- so the right fix is
to make the geometry unambiguous, not to bypass detection.

**How it is enforced**:

- `extract_centerline_from_step` uses **classification-based single
  dispatch** (5 positive-match predicates: multi-station loft / united
  multi-turn / revolution+plane / OPEN / CLOSED), no `try/except`
  cascade across paths.
- `_centerline_from_topology_spine` is **CLOSED-only**: it raises if
  called with `topo.is_open=True` (programming-error indicator, not a
  soft-fallback signal).
- `_find_lateral_surface` checks UV-closure inline -- when it returns
  a face, downstream UV sampling MUST succeed (no `try/except` in
  `filaments_from_step` Path 1).
- `_check_filaments_cover_solid_bbox` raises with a HINT pointing at
  CAD regeneration or BEM-A switch -- it never suggests
  "pass --path-points-m" because that escape hatch no longer exists.

**Failure mode users should expect**: a hard `ValueError` with a
diagnostic that names the axis / overshoot / gap.  Users address it
by (a) regenerating the STEP with cleaner loft vertex alignment so
the lateral surface is a single dominant BSPLINE / TORUS, OR (b)
switching the panel to `--coil-solver bem-a --coil-vol <pre-meshed.vol>`
which bypasses spine extraction entirely (meshed conductor in, no
filament topology needed).

### Field Comparison: Vector Difference

**POLICY**: Compare magnetic fields using **vector difference** `norm(B1 - B2)`, not scalar magnitude difference `abs(|B1| - |B2|)`. Magnetic field is a vector quantity.

### Scattered-Field Robin RHS: Removed (2026-04-24)

**POLICY**: The PEEC Biot-Savart scattered path
(`solve_fem_biot_savart`, `--source-mode scattered` in
`calc_fem_kelvin.py`) was **removed 2026-04-24** after giving a
~3.4x P_wp under-prediction that could not be traced to a clean
formulation fix (the empirical correction factor landed at ~0.78,
not the derivation's +1.0; remainder untraced).

Production paths for PEEC + FEM-SIBC:
- **P_wp**: `calc_inductance.py --coil-solver peec --vol <wp>`
  (PEEC+BEM weak coupling; emits L_coil + ΔL_telegen + P_wp) or
  `calc_fem_coilmesh.py` (full FEM with volumetric coil).
- **L_total**: `calc_fem_coilmesh.py` (volumetric coil + workpiece
  SIBC + Kelvin, intrinsic back-reaction) OR `calc_fem_kelvin.py`
  with the remaining total-field line-integral RHS.

`calc_peec_bem.py` was unified into `calc_inductance.py` in v4.25.0
(2026-05).  References to the old script name in older release notes
or knowledge files predate that refactor.

`kelvin_source.biot_savart_A_cf` / `biot_savart_B_cf` helpers
remain available for research uses (line-integral back-reaction,
validation harnesses). See `memory/scattered_robin_rhs_bug.md` for
the full investigation record.

### FMM (Fast Multipole Method): Removed (2026-03-06)

**ExaFMM-t was removed from the repository**. Do NOT re-implement FMM acceleration.

**Why FMM failed for Radia**:

1. **Point multipole approximations are poor for extended element sources**: distributed face/volume sources cannot be collapsed to a single dipole at intermediate distances (r ~ 2-5 element sizes) without unacceptable engineering error.

2. **FMM Solve (Method 3) was useless**: Compact geometries (C-type magnets, iron yokes) have 87% near-field pairs. Near-field correction memory equals the full dense matrix, eliminating FMM's O(N log N) advantage. HACApK (H-matrix, Method 2) is 10-100x faster because ACA+ compression works on the same near-field blocks.

3. **FMM field evaluation had no benefit over direct**: For typical Radia models (N < 10,000 elements), direct B_genComp with TaskManager parallelization is fast enough. FMM overhead (tree build, M2L translation) exceeds direct computation time for these sizes.

4. **HACApK covers all large-scale needs**: H-matrix acceleration (ACA+) provides O(N log N) memory and O(N log^2 N) MatVec for the interaction matrix, which is the actual bottleneck.

**Lesson**: FMM is effective for point charges/dipoles in unbounded space (N-body) and for smooth BEM kernels. It is not Radia's default for compact core operators with extended element sources and near-field-heavy geometry.

### GmshBuilder: Removed (2026-03-13)

**GmshBuilder was removed from the repository**. Do NOT re-implement GMSH-based mesh generation.

**POLICY**: GMSH is used for **visualization and post-processing only**, NOT for mesh generation.

**Mesh generation is 2-path only**:
1. **STEP -> Netgen** (via NGSolve OCC): For tet meshes with `mesh.Curve(order)` support
2. **STEP -> Cubit** (Coreform Cubit): For structured hex meshes and complex topology

**Do NOT**:
- Use GMSH Python API (`gmsh.model.occ.*`) for geometry or mesh creation
- Import `from radia.gmsh_builder import GmshBuilder` (removed)
- Write new GMSH mesh generation scripts

**GMSH is allowed for**:
- Opening and visualizing `.msh` files (GMSH GUI)
- Post-processing field data (GMSH views)
- Reading `.msh` file format for visualization verification

### GMSH .msh Format Version Policy (2026-04-15 update)

**POLICY**: **全リポジトリで GMSH .msh v4.1 のみ**。v2.2 は全廃。
netgen の I/O は常に **`.vol` 経由** を研究室の正式な運用プロセスとする。

| Component | Output | Purpose |
|-----------|--------|---------|
| **Cubit plugin** (`export gmsh`) | `.msh v4.1` | Mesh export → GMSH viewer |
| **Cubit plugin** (`export netgen`) | `.vol` | NGSolve mesh interchange |
| **Radia post** (`GmshPostExport` / `vol2msh`) | `.msh v4.1` | Field post-processing |

**Shared routine layout** (both mesh_export and post_export emit the same v4.1 structure):
1. `$MeshFormat 4.1 0 8`
2. `$PhysicalNames` (blocks → material names, sidesets → boundary names)
3. `$Entities` (one per physical group, linked via `physicalTag`)
4. `$Nodes` (block-structured, one block per entity)
5. `$Elements` (block-structured, one block per entity × element type)
6. `$NodeData` / `$ElementData` (post-processing only)

**Netgen interchange**:
- Cubit → NGSolve: `.vol` only (never `.msh`).  No `ReadGmsh` path.
- NGSolve → GMSH view: `.vol` + `.sol` → `vol2msh()` → `.msh v4.1`.
- `.msh v2.2` support has been removed from `GmshPostExport` and
  `ExportGmshCommand`.  The `version` keyword on `export gmsh` is
  accepted for back-compat with old `.jou` files but is ignored
  (always emits v4.1 with a warning if `version 2` is passed).

### Mesh Export Consistency Check Policy

**POLICY**: Run the canonical `check-vol` gate after every solver-bound `.vol`
is exported/generated and before solver or Simulink initialization.

| Gate | Always available | Failure condition |
|------|------------------|-------------------|
| NGSolve structure | Yes | Cannot load, invalid counts/order, or malformed topology |
| Label contract | Yes | Missing/unexpected/generated/case-colliding labels or invalid Kelvin/source/symmetry relations |
| Curved mapping | Default | Non-positive Jacobian, scaled-Jacobian threshold, wrong curve order, or element-family mismatch |
| CAD consistency | When `.vol.json` exists | Per-label volume/area or total edge-length error above threshold; count/order mismatch |

Production application/mode meshes use
`--contract <labels.json> --strict-labels`; strict naming permits lower snake case,
`sym_bn=0_<axis>` / `sym_ht=0_<axis>`, and the reserved `GND` point label.
The default CAD threshold is 1%, and errors are displayed in scientific
notation. The sidecar is auto-discovered when present; explicitly passing
`--json` makes it mandatory.

```bash
check-vol model.vol --contract labels.json --strict-labels \
  --report-json run/vol_check.json
```

The checker validates mesh topology and semantic region names. Material
constants remain in the selected DesignSpec/configuration and are checked
there; never infer conductivity, permeability, BH data, or frequency from a
`.vol` label.

### GMSH API Node Ordering Verification Policy

**POLICY**: All high-order mesh exports (.msh, .bdf, .vtk) MUST be verified using the **GMSH API** (`getJacobians`). Do NOT rely on custom parsers or ReadGmsh for HO verification.

**Why**: GMSH's `getJacobians()` computes Jacobian determinants using its own isoparametric basis functions. Negative determinants indicate inverted elements caused by incorrect node ordering. This is the authoritative test because it verifies that GMSH itself correctly interprets the exported file.

**Test procedure** (`tests/cubit/test_ho_volume_all_formats.py`):

```python
import gmsh
gmsh.initialize()
gmsh.open("exported.msh")

# Get integration points and Jacobians
etypes, etags, ntags = gmsh.model.mesh.getElements(dim=3)
for et in etypes:
    local_coords, weights = gmsh.model.mesh.getIntegrationPoints(int(et), "Gauss4")
    jac, det, pts = gmsh.model.mesh.getJacobians(int(et), local_coords)
    # Volume = sum(det * weight), negative det = inverted element
```

**Pass criteria**:
1. **No negative Jacobian determinants** (node ordering correct)
2. **Volume error < 1%** vs analytical (mid-node positions correct)
3. **Cross-format consistency** (all formats produce same volume within 0.01%)

**Netgen reference coordinate mapping** (root cause of past bugs):
- `ref(0,0,0)` -> `el.vertices[3]` (NOT `el.vertices[0]`)
- `ref(1,0,0)` -> `el.vertices[0]`
- `ref(0,1,0)` -> `el.vertices[1]`
- `ref(0,0,1)` -> `el.vertices[2]`
- Vertex `el[i]` -> `lam[(i+1)%4]` (barycentric index shifted by +1 mod 4)

### Cubit Learn Edition 50k Element Cap: IGNORE — exceeding 50k is OK

**POLICY**: **Exceeding 50,000 elements is OK and expected.** Tune the
mesh to whatever density the physics needs (typically the BEM-validated
reference density) and let the element count fall where it may.

Cubit Learn Edition prints
`ERROR: Coreform Cubit - Learn Edition restricts export to models with
less than 50k elements.` whenever it sees a model larger than 50k. This
ERROR is **harmless for the Radia workflow**:

- The Radia in-tree `export netgen` plugin **bypasses the cap**
  and writes the .vol successfully regardless of the warning. Verified
  2026-04-12 with a 147,234-element coil model:

      *****ERROR: Coreform Cubit - Learn Edition restricts export to
      models with less than 50k elements.
      Phase1f: 147234 volume elements added
      Exported Netgen Vol (order 1): ih_fem_sample.vol
                                    (25301 nodes, 147234 elements)

  The "ERROR" line and the "Exported" line BOTH appear in the same run.
  The export succeeded.

- The cap only applies to Cubit's own built-in exporters
  (`export gmsh` / `export vtk` / `export exo`), which the Radia
  workflow does not use.

- The deployed LAB / 100号機 / mdx machines run Cubit Pro and never see
  the warning anyway. hibino is a PyPI package / MCP consumer; Cubit is
  optional there and the QUAD deploy skips the Cubit plugin smoke when
  Coreform Cubit 2025.12+ is absent.

**Do NOT** coarsen sample .jou meshes to "fit under" the 50k cap, and
**do NOT** treat the Learn Edition ERROR line as a failure when the
"Exported" line is present in the same Cubit run.

### Cubit Block/Sideset Label Convention

**POLICY**: Separate blocks for material and boundary labels. Do NOT mix volume elements and surface elements in the same block.

```python
# CORRECT: separate blocks
block 1 add volume 1           # material block
block 1 name "iron"
block 2 add tri in surface 1   # boundary block  
block 2 name "source"

# ALSO CORRECT: use sideset for boundaries (preferred for FEM)
sideset 1 add surface 1
sideset 1 name "source"

# WRONG: mixed volume + surface in one block
block 1 add volume 1
block 1 add tri in surface 1   # tris invisible via get_block_tris when type=TETRA
block 1 name "mixed"           # boundary label LOST
```

**Label priority** (Netgen Vol export):
- Material: block (volume membership) > entity name > `volume_N`
- Boundary: sideset > block (tri/face overlap) > entity name > `surface_N`
- Edge (BBND): sideset on curve > entity name (TODO: SetCD2Name crashes, needs investigation)

### Journal File Portability Policy

**POLICY**: `.jou` and Cubit `.py` files MUST NOT hardcode entity IDs.
Humans pick IDs by visual inspection in the GUI; LLMs and automation MUST
identify entities from **geometric properties** (area, volume, centroid,
distance) and material/sideset names.

**Why**: IDs change between Cubit versions, after imprint/merge, after
edit-rebuild cycles, and depending on operation order. Hardcoded IDs make
scripts silently produce the wrong geometry (e.g., a 50% half-coil instead
of a 355° gapped torus).

**Guidelines**:
- **No hardcoded volume/surface/curve/vertex IDs.** Use `cubit.get_last_id()`
  immediately after the `create` call, then thread the variable through.
- **Identify entities by geometric predicates**:
  - "the surface whose area ≈ π·a²" -> coil gap face
  - "the volume whose centroid x ≈ offset_x" -> Kelvin sphere
  - "the surface adjacent to material 'kelvin'" -> Kelvin boundary
- **Prefer named blocks/sidesets** over IDs in downstream commands
  (`block 3 name "kelvin"` -> let the C++ exporter look up by name).
- **C++ export plugin auto-detection**: Prefer boundary detection from
  material topology (e.g., Kelvin inner/outer detected from block names
  "air"/"kelvin") over manual sideset assignment.

```python
# CORRECT: capture IDs immediately, identify by geometry
cubit.cmd('create surface curve 1')
cubit.cmd('sweep surface 1 axis 0 0 0 0 0 1 angle 355')
coil_vid = cubit.get_last_id("volume")
cubit.cmd('block 1 add volume %d' % coil_vid)
cubit.cmd('block 1 name "coil"')

# CORRECT: detect gap faces by area
A_gap = math.pi * a_coil**2
gap_faces = [s for s in cubit.parse_cubit_list("surface", "in volume %d" % coil_vid)
             if abs(cubit.surface(s).area() - A_gap) / A_gap < 0.05]

# WRONG: hardcoded IDs
block 1 add volume 1                # may be a half-coil after webcut!
sideset 1 add surface 2             # surface 2 may not be the gap face
```

---

## Architecture Overview

### Terminology: HDiv-VIM / BEM / PEEC

**POLICY**: Radia's production soft-iron demagnetization route is
**HDiv-VIM**. Do not revive retired collocation demag backend names or aliases.
Use precise terminology:

| Term | Method | Library | Description |
|------|--------|---------|-------------|
| **HDiv-VIM** | FEEC flux / charge-Gram volume integral method | Radia C++ + `radia.vim` | Soft-iron demag path, NGSolve mesh/FES coupled |
| **BEM** | Boundary Element Method | **ngsolve.bem** | EFIE/MFIE, HDivSurface, Maxwell/Laplace kernels |
| **PEEC** | Partial Element Equivalent Circuit | Radia Python + C++ | Loop-Star, circuit extraction (L,R,C,M) |

**Do NOT** blur these:
- BEM (ngsolve.bem): surface integral equations (EFIE/MFIE) on conductor/dielectric boundaries
- HDiv-VIM: volume/mesh-based soft iron demag, charge-Gram/HACApK acceleration,
  reduced-FEM/NGSolve coupling
- PEEC: conductor circuit extraction and MQS/Darwin surface-impedance analysis

**Decision (2026-07-05, Sugahara): Radia soft-iron demag is HDiv-VIM only.**
`demag_backend="auto"` and `"hdiv"` are the supported names; `"auto"` selects
the HDiv-VIM path for soft iron. Mesh-backed operation goes through `.vol` ->
NGSolve `Mesh` -> `radia.vim.soft_iron_from_mesh` (or
`soft_iron_from_vol(...)`). Hand-built element primitives are an internal
representation detail, not a parallel user-facing demag backend.

**When to use which**:
- Soft iron, nonlinear demag, reduced-FEM coupling -> **HDiv-VIM** (Radia)
- Permanent-magnet source fields -> **rad.Fld / analytical Radia field kernels**
- Eddy currents, shielding, impedance extraction → **BEM** (ngsolve.bem) or **PEEC** (Radia)
- High-frequency scattering → **BEM** (ngsolve.bem, Helmholtz kernel)

### Development Strategy: Complement NGSolve

Radia's role is to **complement NGSolve**, not compete with it. Focus on areas where FEM is weak.

```
┌─────────────────────────────────────────────────────────────────┐
│                    Electromagnetic Analysis                      │
├─────────────────────────────────────────────────────────────────┤
│  NGSolve (FEM)              │  Radia (HDiv-VIM/PEEC/field)      │
│  ───────────────────────────│──────────────────────────────────│
│  OK: Bounded domains        │  OK: Unbounded domains (open BC) │
│  OK: Complex geometry       │  OK: Permanent magnets (no mesh) │
│  OK: Nonlinear materials    │  OK: Thin conductors (PEEC)      │
│  OK: Transient analysis     │  OK: SPICE circuit extraction    │
│  OK: Multi-physics coupling │  OK: Model order reduction (MOR) │
│  WEAK: Open boundary (PML)  │  OK: Natural open boundary       │
│  WEAK: Thin structures      │  OK: Surface impedance (SIBC)    │
│  WEAK: Circuit parameters   │  OK: L, R, C, M extraction       │
└─────────────────────────────────────────────────────────────────┘
```

### Accelerator Magnet Solver Architecture

The complete pipeline for accelerator electromagnet analysis:

```
┌─────────────────────────────────────────────────────────────────┐
│  CoilBuilder                                                     │
│  ─────────────────────────────────────────────────────────────  │
│  add_straight() / add_arc() → continuous path (beam optics style)│
│  close() → multi-variable optimization for loop closure          │
│  mirror() / rotate_copies() → symmetry (dipole, quadrupole)     │
│  to_radia() → Biot-Savart source Hs (NO coil mesh needed)       │
│  write_step() → GMSH visualization                              │
└──────────────────────────┬──────────────────────────────────────┘
                           │ Hs (analytical)
┌──────────────────────────┼──────────────────────────────────────┐
│  Cubit                   │                                       │
│  ─────────────────────────────────────────────────────────────  │
│  Iron yoke + air + Kelvin domain → hex sweep mesh                │
│  export_NGSolveCurvedMesh(order=N) → arbitrary-order hex mesh    │
│  CallbackGeometry → ACIS direct curving (NO STEP/OCC)            │
└──────────────────────────┬──────────────────────────────────────┘
                           │ curved hex mesh
┌──────────────────────────┼──────────────────────────────────────┐
│  NGSolve FEM             │                                       │
│  ─────────────────────────────────────────────────────────────  │
│  Omega-reduced Omega (2-scalar, fastest formulation)             │
│  Kelvin transformation (open boundary, no PML)                   │
│  Energy-based B-input Play model (nonlinear hysteresis)          │
│    - Reversible/irreversible separation → convex energy          │
│    - LU factored ONCE, back-substitution iteration               │
│    - Fast inverse: Picard 2-3 iter (vs Newton ~100 iter)         │
│  GmshPostExport → GMSH visualization (+ coil STEP overlay)      │
└─────────────────────────────────────────────────────────────────┘
```

**Unique capabilities** (no other software provides all of these):
- Coil mesh NOT needed (Biot-Savart analytical source)
- STEP/OCC NOT needed (ACIS CallbackGeometry)
- Hex mesh with arbitrary-order curving
- Energy-based hysteresis with fast inverse
- Open source, `pip install radia`

**Do NOT Implement** (use existing libraries):
- FEM solvers (use NGSolve)
- General sparse solvers (use MKL/MUMPS)
- Full-wave BEM (use ngsolve.bem for high frequency)
- CAD geometry kernels (use OpenCASCADE via NGSolve)
- Mesh generation wrappers (use Netgen or Cubit directly, NOT GMSH)
- PEEC from scratch (use PAMELA)
- Custom H-matrix algorithms (use HACApK)

**Radia C++ Core** (maintain and enhance):
1. HDiv-VIM soft-iron demag and charge-Gram kernels
2. Field computation - B, H, A, Phi in unbounded domains
3. PEEC / BEM Laplace-kernel matrix assembly where Radia owns the workflow
4. NGSolve integration - RadiaField CoefficientFunction

### Soft-Iron Solver Method: HDiv-VIM

HDiv-VIM is the supported soft-iron demagnetization method. It is built around
NGSolve mesh/FES concepts, HDiv flux continuity, and Radia's charge-Gram /
HACApK acceleration path. Retired collocation demag names are not supported as
public backends.

NGSolve family selection is explicit in documentation: `HDiv(mesh, order=p)`
is BDM, while Raviart--Thomas requires `HDiv(mesh, order=p, RT=True)`.
Radia's established `vim.Solve`, `PlanarDemagBody`, `MagnetizationSource`,
charge-Gram, nonlinear, IMA, and `rad.Fld` production paths use BDM1/BDM2.
Do not call these RT1/RT2.  The actual RT family is an explicit comparison or
research path until separately promoted.

### Unified Field Computation Architecture

**POLICY**: All field computation MUST use `rad_field_unified.h/cpp`.

```
┌─────────────────────────────────────────────────────────────────┐
│                    rad_field_unified.h/cpp                       │
│  ─────────────────────────────────────────────────────────────  │
│  ComputeFieldSingle()     - Single point, static field          │
│  ComputeFieldBatch()      - Batch points, TaskManager parallelized │
│  ComputeComplexFieldSingle() - Complex (AC) field               │
│  ComputeComplexFieldBatch()  - Complex batch with TaskManager   │
│  IsPointInsideAnyElement() - Inside/outside classification      │
│  ComputeBFromMagnetization() - Dipole field from M (complex)    │
└─────────────────────────────────────────────────────────────────┘
                              │
           ┌──────────────────┼──────────────────┐
           ▼                  ▼                  ▼
    ┌─────────────┐    ┌─────────────┐    ┌─────────────┐
    └─────────────┘    └─────────────┘    └─────────────┘
```


**Key Features**: Inside/outside classification, TaskManager parallelized batch, complex field support for PEEC/SIBC workflows.

### Field Calculation: Surface Current vs Surface Charge

- **ObjRecMag**: Surface current model (rectangular blocks). 8-corner BufVect formula, efficient and non-cancelling on symmetry axes.
- **ObjHexahedron/ObjTetrahedron**: Surface charge model (general polyhedra). Face-based solid angle integration. A field may be zero on symmetry axes (mathematical cancellation, not a bug).

**rad.Fld() inside materials**: Prefer external field probes or HDiv/NGSolve
state variables for validation. Interior material fields depend on the active
formulation and should not be compared against retired collocation conventions.

### Vector Potential A Field

A field is **implemented** for all element types using face integration (Wilton et al. formula). Formula: `A = (mu_0/4pi) * (M x BufVect)`. Satisfies `B = curl(A)` (verified numerically). Verification script: `validation_test/ngsolve_integration/verify_curl_A_equals_B/`.

### User-Facing Element APIs

- `rad.ObjRecMag(center, dimensions, magnetization)` -- Rectangular magnets (optimized formulas)
- `rad.ObjHexahedron(vertices, magnetization)` -- Arbitrary hexahedra (8 vertices)
- `rad.ObjTetrahedron(vertices, magnetization)` -- Tetrahedra (4 vertices)
- `rad.ObjWedge(vertices, magnetization)` -- Wedges (6 vertices)
- Mesh import functions (`netgen_mesh_to_radia`) for complex geometries

---

## API Guardrails

### Common Mistakes Checklist

**1. ObjBckg Requires Callable (CRITICAL)**
```python
bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])  # CORRECT
bkg = rad.ObjBckg([0, 0, 0.1])             # WRONG - not a callable
```

**2. UtiDelAll() Cleanup**: Every script must call `rad.UtiDelAll()` before exiting.

**3. Relative Path Imports**:
```python
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))  # CORRECT
sys.path.insert(0, r'S:\Radia\01_GitHub\src\radia')  # WRONG - machine-specific
```

**4. MatLin Usage**: For isotropic materials, ALWAYS use single-argument form `MatLin(mu_r)`. MatLin is for soft magnetic materials only -- permanent magnets specify magnetization directly in `ObjHexahedron(vertices, [Mx, My, Mz])`.

**5. Docstring Units**: Use "in constructor length units", not "in mm".

**6. State Mutation**: Computation methods must NOT leave object state inconsistent on exception.

### Background Field API

```python
bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])      # Uniform 0.1T in z
bkg = rad.ObjBckg(quadrupole_field_function)    # Spatially varying
container = rad.ObjCnt([mag_obj, bkg])
rad.Solve(container, 0.0001, 1000, 1)
```
Legacy `ObjBckg([Bx, By, Bz])` array form is NOT supported. Callback receives `[x, y, z]` in current units and returns `[Bx, By, Bz]` in Tesla.

### Memory Management

```cpp
// Exception-safe pattern
Type* ptr = nullptr;
try {
    ptr = new Type(...);
    Handle h(ptr);
    ptr = nullptr;  // Ownership transferred
} catch(...) {
    if(ptr) delete ptr;
    Initialize();
    return 0;
}
```
Prefer RAII containers (`std::vector`) over manual `new`/`delete`.

### Deprecated Relaxation API

| Deprecated | Replacement |
|------------|-------------|
| `RlxPre()`, `RlxMan()`, `RlxAuto()` | `rad.Solve(obj, prec, maxiter, method)` |
| `RlxUpdSrc()`, `SetRelaxSubInterval()` | `rad.Solve()` |

---

## Build & Release

### Target Versions

| Component | Version | Notes |
|-----------|---------|-------|
| **Python** | 3.12.10 | System Python for Radia/NGSolve. Cubit toolbar launches notebook/headless workflows via subprocess. |
| **Coreform Cubit** | 2025.12 | Embedded Python 3.10 + PySide6. Cannot import NGSolve/Radia directly. |
| **NGSolve** | 6.2.2606 | pinned by `pyproject.toml`; BEM correctness, thread-local TaskManager/LocalHeap, curved mesh support |

**Cubit panel subprocess constraint**: Cubit embeds Python 3.10; Radia/NGSolve require 3.12. Same-process import is impossible. All computation runs via `subprocess.run([python3.12, calc_*.py])` with JSON output.

### Build: MSVC + Intel MKL

**POLICY**: Use **MSVC** compiler with **Intel MKL**. Intel oneAPI compiler (icx-cl) is NOT compatible with NGSolve linking.

```powershell
pwsh -NoProfile -ExecutionPolicy Bypass -File "Build.ps1"
pwsh -NoProfile -ExecutionPolicy Bypass -File "Build.ps1" -Rebuild  # Clean rebuild
```

**Required Software**: Visual Studio 2022 (MSVC), Python 3.12, pinned
NGSolve/Netgen, pybind11, and pip `mkl-devel`.

### BLAS/LAPACK Runtime Boundary

**POLICY**: NGSolve owns its packaged OpenBLAS runtime. Radia's HACApK,
PARDISO, and dense native kernels use Intel MKL from the selected Python
environment's pip `mkl-devel`/`mkl` packages. These runtimes may coexist, but
Radia must not replace NGSolve's BLAS or bundle an alternative NumPy build.
`MKLROOT` is an explicit controlled fallback, not a machine-wide default.

### Parallelization: NGSolve TaskManager

**POLICY**: Follow the **NGSolve-native execution model**. Use **NGSolve TaskManager** for
thread-level parallelization, NOT raw OpenMP parallel regions, `std::thread`, `std::async`,
or a second project-local thread pool.

NGSolve's TaskManager provides work-stealing task-based parallelism and is the single shared
threading substrate for Radia + NGSolve workflows. All new Radia parallel code should use
`ngcore::ParallelFor` / `ParallelForRange`; long C++ entry points that can be called directly
from Python should stand up or reuse an `ngcore::RegionTaskManager` so a bare `rad.Solve(...)`
does not accidentally run TaskManager loops serially. Python/NGSolve assembly code follows
NGSolve convention and is caller-wrapped with `with ngsolve.TaskManager():`.

External threaded kernels are the exception, not an alternative Radia threading model. Dense
LU (`dgesv_`) uses MKL's internal threading under `radia::MKLThreadGuard` while
`ngcore::SuspendTaskManager` prevents nested TaskManager/MKL oversubscription. Iterative
Radia solvers (BiCGSTAB, HACApK/method 2, HDiv VIM) stay TaskManager-native.

```cpp
// CORRECT: NGSolve TaskManager
#include <core/taskmanager.hpp>
ngcore::RegionTaskManager rtm(std::max(1, ngcore::TaskManager::GetMaxThreads()));
ngcore::ParallelFor(ngcore::IntRange(n), [&](size_t i) {
    // compute...
});

// AVOID: raw OpenMP / private thread pools
#pragma omp parallel for
for (int i = 0; i < n; i++) { ... }
```

**When to use TaskManager**:
- Field computation loops (ComputeFieldBatch)
- Interaction matrix assembly
- BiCGSTAB vector operations and matrix-vector products
- HACApK H-matrix build/matvec/solve loops
- HDiv VIM C++ solve loops
- Any embarrassingly parallel loop

**When non-TaskManager threading is acceptable**:
- MKL internal threading for dense BLAS/LAPACK/PARDISO calls, guarded by
  `SuspendTaskManager` + `MKLThreadGuard` where Radia controls the call
- Legacy code only until it is migrated; do not add new OpenMP regions

### PyPI Release Workflow (Automated via GitHub Actions)

**POLICY**: Build and publish one immutable commit. Main-push CI provides fast
source and contract feedback; a release tag starts the native mdx artifact
lane. GitHub Release publication still requires `release-quad done`.

**Release Flow**:
1. Make a distribution-focused change; synchronize version files and changelog
   when the distribution version changes.
2. Run the impact-scoped mdx preflight and focused local tests.
3. Review/rebase, merge through a pull request, and push `main`.
4. Require the exact main SHA's check-runs to pass with
   `python tools/release_quad.py ci-verify`.
5. Tag that exact SHA; the native `build-test.yml` lane builds in an isolated
   mdx environment and retains its wheel/MEX/SLX artifacts.
6. Publish from the accepted artifact, then run the four-machine deployment and
   `release-quad done` before publishing the GitHub Release:
   ```bash
   git tag vX.Y.Z
   git push origin vX.Y.Z
   ```
7. Monitor with `python tools/check_ci.py --sha <SHA> --watch`.

**General User Install** (after PyPI publish):
```bash
pip install radia[cubit]
cubit-plugin-install       # Cubit .ccm backend + PySide6 toolbar
```

**CI/CD Pipeline** (`.github/workflows/`):
```
main / PR -> impact-scoped source and package contracts
v* tag    -> mdx native build -> retained exact-SHA artifact
accepted artifact -> PyPI OIDC publication -> release-quad -> GitHub Release
```

**No API tokens stored**. Uses PyPI OIDC Trusted Publishers (id-token: write).

### GitHub Automation Tooling

Committed automation prefers the token-aware REST helpers in `tools/gh_api.py`
so CI and fresh clones do not require an interactive GitHub CLI login. `gh` is
supported on LAB and 100号機 for operator PR/Actions work and may be used by a
GitHub Actions job with `GITHUB_TOKEN`; mdx runtime jobs must not depend on a
developer login. `release_quad.py ci-verify` reads SHA-bound check-runs through
the REST helper rather than inspecting a runner's local process or files.

### Distribution Test Policy (2026-04-24, updated 2026-05-19)

**POLICY (上位ルール、2026-05-19 追加)**: **LAB は開発マシン。 菅原研究室で開発する
パッケージは全て editable install を default とする** — 特定の 4 パッケージに限らない
**一般原則**。 source edit が即座にランタイムへ反映される dev loop を維持するため。
- `pip show <pkg>` → `Editable project location:` が `S:\Radia\01_GitHub\...` を指す
  ことを確認できればよい。他の場所 (site-packages snapshot、CI runner clone、PyPI
  install) を指していたら **発見次第 LAB source に戻す**:
  ```
  pip install -e S:\Radia\01_GitHub\packages\<pkg> --no-deps --no-cache-dir
  ```
- LAB 上で `pip install --upgrade <Sugahara-lab-package>` は **禁止** — editable を
  silently 上書きして dev loop を破壊する (2026-04-28 incident、2026-05-19 再発)。
- Upgrade route は **mdx / hibino 側** (PyPI consumer)。 LAB と 100号機 は release 後も
  metadata 同期のため `pip install -e <path> --no-deps --no-cache-dir` で editable を維持する。
- CI/CD 環境 (e.g. `C:\actions-runner\_work\Radia\Radia\...`) は別管理 (NETWORK
  SERVICE 所有)。 LAB の editable pointer がそちらに drift していたら戻す。

**POLICY (2026-09-03 update)**: **QUAD 配布**: LAB と 100号機 は editable
(NAS source、developer/user feedback loop)、mdx と hibino は PyPI install。
mdx は self-hosted CI/preflight と exact-artifact verification を優先し、
`radia` + `cubit-mesh-export` を検証する。machine-wide `radia-mcp` は不要。
hibino は重い validation/benchmark の第一選択であり、PyPI 経由の MCP
consumer として `radia-mcp` も入れる。

**QUAD 配布モデル (2026-06-25)**:

| Stage | マシン | install 形態 | 目的 |
|-------|--------|-----------|------|
| 1 | LAB | `pip install -e .` + `pip install -e packages/cubit-mesh-export` + `pip install -e packages/radia-mcp` | 開発者ループ。最速フィードバック (NAS source 直接編集) |
| 1 | 100号機 | `pip install -e \\192.168.11.100\work\00_CAE\Radia\01_GitHub` + `pip install -e ...\packages\cubit-mesh-export` + `pip install -e ...\packages\radia-mcp` | 共有ユーザ環境も editable。radia と MCP の学習/修正が即反映される |
| 2 | mdx | `pip install radia / cubit-mesh-export` (PyPI) + `cubit-plugin-install --all-users` | self-hosted CI/preflight と exact-artifact verification。compute は hibino unavailable かつ CI idle 時のみ。`radia-mcp` は不要 |
| 2 | hibino | `pip install radia / radia-mcp / cubit-mesh-export` (PyPI) + `cubit-plugin-install --all-users` | 重い validation/benchmark と PyPI 経由の MCP consumer verification |

**変更点 (2026-06-25)**:
- 旧: LAB editable / 100号機 + mdx 両方 PyPI (2-tier).
- 新: LAB + 100号機 editable / mdx + hibino PyPI (QUAD).
- mdx は CI/preflight 優先で `radia-mcp` 不要。hibino は重い計算と PyPI
  経由の MCP consumer verification を担う。

**LAB / 100号機 editable パッケージ**:
- `radia` (LAB: `S:\Radia\01_GitHub`, 100号機: `\\192.168.11.100\work\00_CAE\Radia\01_GitHub`)
- `cubit-mesh-export` (`packages\cubit-mesh-export`)
- `radia-mcp` (`packages\radia-mcp`)
- `mcp-server-document` (LAB: `S:\mcp-server`) -- LAB-private (PyPI 配布なし)

共有worktreeに並行WIPがあるreleaseでは、そのWIPをstash/clean/resetしてはならない。
代わりに同一NAS上へexact SHAのclean release worktreeを作り、
`RADIA_RELEASE_EDITABLE_REPO_LAB` と `RADIA_RELEASE_EDITABLE_REPO_100` で
LAB/100号機から見える各pathを `release_quad all` と `done` の両方へ渡す。
QUADはprocess停止やinstallより前にSHA一致とtracked-cleanを強制する。公開完了まで
release worktreeを保持し、並行WIPがmainへ着地した後に通常のeditable pointerへ戻す。

LAB / 100号機で `pip install --upgrade <pkg>` を流すと editable が静かに上書きされて壊れるので注意。release 後の LAB / 100号機 側 metadata 同期は `pip install -e <path> --no-deps --no-cache-dir` で再 editable 化。`pip install --upgrade` は **mdx / hibino 用** (PyPI から通常通り upgrade).

**mdx / hibino 全ユーザー PyPI install**: `C:\Program Files\Python312`
の machine-wide site-packages に PyPI install。mdx は
`pip install --upgrade radia==X.Y.Z cubit-mesh-export==X.Y.Z`、hibino は
`pip install --upgrade radia==X.Y.Z radia-mcp==X.Y.Z cubit-mesh-export==X.Y.Z`
+ `cubit-plugin-install --all-users` を実行。

**mdx / hibino Cubit plugin (regular file)**:
- `<Cubit>\bin\plugins\cubit_mesh_export.ccm` (regular file from PyPI wheel)
- `<Cubit>\bin\plugins\cubit_mesh_curver.cp312-win_amd64.pyd` (regular file from PyPI wheel)

LAB の `Build.ps1` 出力は **NAS の `S:\Radia\01_GitHub` に書かれるため、LAB / 100号機
editable には反映される**。mdx / hibino の PyPI install には反映されないので、C++/plugin
変更を mdx / hibino で試すには PyPI release を切るのが正規ルート。

### CI Testing Policy

**POLICY**: CI/CD のテストだけでは不十分。Cubit が必要な機能（`export_curved`, Cubit toolbar, BEM extractor）は **Cubit 環境でのローカルテストが必須**。CI は C++ ビルドと基本テスト（Cubit 不要なもの）のみ。

**リリース前の必須テスト**:
1. CI: C++ ビルド + pytest（Cubit 不要テスト）
2. ローカル: Cubit + system Python で `export_curved` テスト（球、トーラス）
3. ローカル: Cubit toolbar と Simulink application block の動作確認

CI が通っても Cubit テストに通らなければリリースしない。

### Cubit Batch Self-Testing Policy

**POLICY**: Codex は Cubit を **完全ヘッドレス** (`-batch -nographics -nojournal`) で起動し、自力で機能試験を走らせること。GUI や人間の操作は不要。

**起動方法** (既存テストのパターン):
```python
import cubit
cubit.init(['cubit', '-nojournal', '-batch', '-nographics',
            '-commandplugindir', <plugin_dir>])
cubit.cmd("create sphere radius 0.05")
cubit.cmd("mesh volume 1")
cubit.cmd('export femeem "C:\\temp\\cub" overwrite')
```

**対象**: `radia_export {gmsh,netgen,nastran,vtk,femeem}`、Cubit toolbar の非 GUI ロジック、BEM extractor、`export_curved`。
Cubit toolbar の実画面確認だけをヘッドレス試験の例外とする。

**前提**: `cubit` は Python API import (`tools/find_cubit.ps1` でインストールを検出)。バッチ起動でライセンス消費あり。

**特記**: FEMEEM エクスポートの出力パスは **40 文字以下** にすること。`inpin.f90::chkinib(filename*40)` が長い Python `tempfile.TemporaryDirectory()` パスを切り詰めて `forrtl severe (29)` を起こす。`C:\temp\<short>\` 等を使う。

### Cubit Driving Policy: APREPRO/Python Headless is Primary; GUI is a User Debugging Aid (2026-08-05)

**POLICY** (Sugahara): AI-agent-driven Cubit work runs through **APREPRO
commands and Python** on the **headless/batch route** — `.jou` playback,
`cubit.cmd(...)` scripting, the batch stdio daemon, and headless dry-run
tools (`cubit_batch_try`, `cubit_mesh_auto`, mesh races). This is the 本命
(primary) path for automation, testing, CI, and validation. The **GUI
session** (`cubit_show`, the persistent GUI daemon, `cubit_snapshot`) exists
**for the USER's visual debugging** — a human watching the model, inspecting
a mesh, capturing a figure — not as the agent's default execution surface.

- Agent workflows (mesh generation, exports, gates, validation runs) default
  to batch/headless; they must never REQUIRE a GUI window to function.
- GUI-session features are maintained and tested as the **user-debugging
  surface** (that is why `validation_test/radia_mcp/test_cubit_session_e2e.py`
  exercises the GUI path: to protect the user-facing debug surface, incl.
  `cubit_snapshot`, which needs a rendering window — batch reports an honest
  ok=false there).
- This restates and extends the "Cubit Batch Self-Testing Policy" above to
  the MCP-server era: batch first, GUI when the user wants to see something.

**Wheel Verification** (automated by Build_Wheel.ps1, also manual):
```python
import zipfile
whl = zipfile.ZipFile('dist/radia-X.Y.Z-cp312-cp312-win_amd64.whl')
for info in whl.infolist():
    if info.filename.endswith('.pyd'):
        print(f'{info.filename}: {info.file_size} bytes')
# Must contain radia/_radia_pybind.pyd (> 2 MB)
# Must NOT contain any .dll files (MKL policy)
```

### MKL DLL Policy: Do NOT Bundle

**POLICY**: PyPI packages MUST NOT bundle Intel MKL DLLs. `pyproject.toml` declares `mkl>=2024.2.0` as dependency; pip installs MKL DLLs to `{sys.prefix}/Library/bin/`. `__init__.py` adds the path via `os.add_dll_directory()`.

**Do NOT**: Copy MKL/Intel OpenMP DLLs into `src/radia/` or include `*.dll` in `package_data`.

### Package Structure

```
src/radia/
  __init__.py           # DLL path setup + re-export from C++ module
  _radia_pybind.pyd     # Main C++ extension (includes RadiaField CoefficientFunction)
  cln_core.pyd          # CLN transient solver
  peec_matrices.pyd     # PEEC matrix assembly
  *.py                  # Python utility modules
  # NO .dll files
```

**Always use `Build.ps1`** for building. Never use manual cmake commands -- the script handles CMake configure + build + `.pyd` copy to `src/radia/`.

---

## Mesh & NGSolve Integration

### NGSolve Version Requirement

**CRITICAL**: Use NGSolve **6.2.2606** as pinned by `pyproject.toml`. Required for the validated BEM fixes, thread-local TaskManager/LocalHeap behavior, curvedelements .vol Load, hex/prism curving, and Periodic BC fix.

Reference: https://forum.ngsolve.org/t/ngsolve-periodic-boundary-condition-regression-bug-report/3805

Official PyPI ngsolve **6.2.2606** includes the Periodic BC fix,
**curvedelements Save/Load**, **p-version hex/prism curving**, and the
current `ngsolve.bem` APIs used by Radia. Radia links Intel MKL for its own
BLAS/LAPACK/PARDISO calls; NGSolve's wheel dependencies are managed by PyPI.

**Netgen fork is no longer required.** The ksugahar/netgen repository is historical only.
All curvedelements, CallbackGeometry, and curving features are now in the official release.

```bash
pip install radia[cubit]       # Installs Radia + cubit-mesh-export dependency
cubit-plugin-install           # Deploys Cubit .ccm backend + PySide6 toolbar
```

### SetGeomInfo API (Netgen PR#232)

SetGeomInfo is no longer needed for typical workflows. The Cubit plugin handles
UV coordinates and geometry projection internally via CallbackGeometry (embedded in C++).
PR: https://github.com/NGSolve/netgen/pull/232 (historical reference)

### NGSolve Recommended Configuration

```python
fes = HDiv(mesh, order=2)  # BDM2 (NGSolve default)
B_gf = GridFunction(fes)
B_gf.Set(rad.RadiaField(radia_obj, 'b'))  # C++ CoefficientFunction in _radia_pybind.pyd
```

- Evaluate GridFunction at distances > 1 mesh cell from magnet surface
- Use CoefficientFunction directly for maximum accuracy near boundaries
- Avoid GridFunction evaluation within 1 mesh cell of magnet surface

### NGSolve Magnetization → Radia Open Boundary Field Evaluation

NGSolve FEM solves M(x) inside bounded domains but struggles with open boundary (PML needed). Radia provides natural open boundary evaluation using **exact analytical formulas** (NOT dipole approximation).

```
NGSolve FEM Solve → M per element → netgen_mesh_to_radia() → Radia objects → rad.Fld()
```

**POLICY**: Do NOT use dipole approximation (m=M*V) for NGSolve → Radia pipeline. Register elements as proper Radia ObjHexahedron/ObjTetrahedron with solved magnetization. Radia's surface charge/surface current analytical formulas are exact for constant M per element, with no approximation error at any distance.

**Use cases**:
- External field from FEM-solved nonlinear iron core (no PML needed)
- Stray field evaluation at large distances (exact, not approximate)
- Particle trajectory through FEM-solved magnet assembly
- NGSolve CoefficientFunction for coupling back into FEM

**Workflow**:
```python
import radia as rad
from ngsolve import *
from radia.netgen_mesh_import import netgen_mesh_to_radia

rad.UtiDelAll()

# 1. NGSolve solves nonlinear problem → M per element
# (user's FEM solve code here)

# 2. Convert mesh to Radia objects with per-element magnetization
def material_from_ngsolve(el_idx):
    M = get_element_magnetization(gf_M, mesh, el_idx)  # user function
    return {'magnetization': M.tolist()}

container = netgen_mesh_to_radia(mesh, material=material_from_ngsolve, units='m')
# No Solve() needed - M is already known from NGSolve

# 3. Evaluate field at arbitrary external points (exact analytical formulas)
B = rad.Fld(container, 'b', [0, 0, 0.1])          # single point (shape (3,))
B_batch = rad.Fld(container, 'b', obs_points)      # batch (shape (N,3))
```

**Why Radia objects, not dipoles**:
- Surface charge model: exact for constant M, zero approximation error
- Near-field: no distance limitation (dipoles fail at r < 2*element_size)
- `netgen_mesh_to_radia()` already supports per-element material via callable

### NGSolve Mesh Access Policy

**POLICY**: All mesh access MUST use functions from `src/radia/netgen_mesh_import.py`. NEVER directly access `mesh.ngmesh.Points()` or `el.vertices[].nr` -- NGSolve has two indexing schemes (0-indexed vs 1-indexed) that cause off-by-one errors.

```python
# CORRECT
from netgen_mesh_import import netgen_mesh_to_radia, extract_elements
radia_obj = netgen_mesh_to_radia(mesh, material={'magnetization': [0,0,0]}, units='m')

# WRONG - index confusion
pt = mesh.ngmesh.Points()[v.nr]  # Off-by-one!
```

### Mesh Generation Policy

**POLICY**: Mesh generation uses **2 paths only**. GMSH is NOT used for mesh generation.

| Path | Workflow | Element Types | Use Case |
|------|----------|---------------|----------|
| **STEP -> Netgen** | STEP -> NGSolve OCC -> `Mesh()` | Tet4 (+ `mesh.Curve(order)`) | General purpose, curved boundaries |
| **STEP -> Cubit** | STEP -> Coreform Cubit -> `.msh` export | Hex8, Wedge6, Tet4 | Structured hex, complex topology |

**Radia supports 1st order only** (Tet4, Hex8, Wedge6). 2nd order planned.

**CRITICAL**: For curved geometries, use `mesh.Curve(3)` after Netgen meshing. Without it, polygon approximation of circles loses ~2% area -> ~9% force error.

### Mesh Import Paths

Both paths produce `.vol` files consumed by `Mesh("model.vol")`:

```
Path A: Cubit (recommended for hex)
  STEP -> Cubit -> export netgen "model.vol" order N -> Mesh("model.vol") -> Radia

Path B: OCC (recommended for tet)
  STEP -> NGSolve OCC -> Mesh() -> netgen_mesh_import.py -> Radia
```

**Key import functions**:

| Module | Function | Purpose |
|--------|----------|---------|
| `netgen_mesh_import` | `netgen_mesh_to_radia(mesh, ...)` | NGSolve mesh -> Radia (recommended) |
| `netgen_mesh_import` | `create_hex_mesh_grid(...)` | Structured hex grid (no external tool) |

### Cubit Mesh Export (cubit-mesh-export)

For high-order curved mesh export from Coreform Cubit, use the **cubit-mesh-export** package.
Mesh export is C++ only (`export netgen` APREPRO command in the Cubit plugin).

**Install**: `pip install cubit-mesh-export` (or `pip install radia[cubit]`)
**Source**: `packages/cubit-mesh-export/` in the Radia monorepo

**Consistency checking** (does NOT require Cubit):
```bash
check-vol model.vol                         # Sidecar auto-discovered if present
check-vol model.vol --contract labels.json --strict-labels \
  --report-json run/vol_check.json
```
```python
from cubit_mesh_export.check import check_consistency, check_label_contract
```

**Module names**:
- `cubit_mesh_export` — canonical Python package (PyPI: cubit-mesh-export)
- `cubit_mesh_curver` — C++ pybind11 module (bundled in cubit_mesh_export)
- `check_vol_consistency` — thin backward-compat re-export in `src/radia/panels/` (imports from cubit_mesh_export.check)

Cubit workflow for journal files: define blocks before export, use the Cubit plugin commands (`cubit.cmd('export gmsh/jmag_nastran/vtk ...')`). Requires `CUBIT_PLUGIN_DIR` environment variable (set by `cubit-plugin-install`).

### PEEC Conductor Mesh

PEEC conductors use **surface mesh only** (SIBC handles skin effect). Generate surface meshes via Netgen or Cubit. Supported: Tri3, Quad4 (1st order), Tri6, Quad8/9 (2nd order).

### Nastran Format: REMOVED

Nastran BDF support is **REMOVED**. Use Cubit -> `.msh` export or Netgen direct. Cubit can read legacy `.bdf` files if needed.

### Mesh Operations: Dropped APIs

`ObjDivMag`, `ObjDivMagPln`, `ObjCutMag` are NOT supported. All mesh operations use external tools (Netgen, Cubit).

### Mesh File Preservation

**NEVER DELETE** canonical mesh files (`.bdf`, `.nas`, `.msh`, and Cubit-owned
`.vtk` fixtures), Cubit journal files (`.jou`), or mesh generation scripts.
These are difficult to recreate. This protection does not preserve superseded
Radia VTK field output outside the Cubit component.

### Available Mesh Access Functions

From `src/radia/netgen_mesh_import.py`:
- `netgen_mesh_to_radia()` -- Convert entire mesh to Radia (recommended)
- `extract_elements()` -- Extract element data for custom processing
- `compute_element_centroid()` -- Centroid from vertex list
- `create_radia_tetrahedron()` / `create_radia_hexahedron()` -- Single elements
- `create_hex_mesh_grid()` -- Structured hex grid (no external tool)
- Constants: `TETRA_FACES`, `HEX_FACES`, `WEDGE_FACES`, `PYRAMID_FACES` (1-indexed face topology)

---

## H-Matrix Acceleration (HACApK)

### Policy: Use HACApK Only

**POLICY**: Do NOT implement custom H-matrix algorithms. Use the HACApK library at `src/ext/HACApK_LH-Cimplm/` (MIT license).

**Solver Methods**:

| Method | Name | Use Case |
|--------|------|----------|
| 0 | LU | Small problems (N < 500), guaranteed convergence |
| 1 | BiCGSTAB | General purpose, medium problems |
| 2 | HACApK | Large problems (N > 1000), O(N log N) memory |

**ソルバー選択ガイドライン**:
- **小規模 (N<500)**: LU推奨 (確実な収束)
- **中規模 (500<N<2000)**: BiCGSTAB推奨 (最速)
- **大規模 (N>2000)**: HACApK推奨 (メモリ効率)

### Solver Configuration (Unified API)

```python
rad.SolverConfig(hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)
rad.SolverConfig(bicgstab_tol=1e-4, relax_param=0.3, newton_method=True)
config = rad.GetSolverConfig()  # Returns dict with all settings
```

| Keyword | Default | Description |
|---------|---------|-------------|
| `hacapk_eps` | 1e-4 | ACA tolerance (1e-6 to 1e-2) |
| `hacapk_leaf` | 10 | Minimum cluster size for H-matrix clustering |
| `hacapk_eta` | 2.0 | Admissibility parameter |
| `bicgstab_tol` | 1e-4 | BiCGSTAB convergence tolerance |
| `relax_param` | 0.0 | Under-relaxation (0=full step, <1=damped) |
| `newton_method` | False | True=Newton-Raphson, False=Picard |
| `newton_damping` | True | Enable Newton line search damping |

See `docs/HMATRIX_EVALUATION.md` for full evaluation report.

### Sign Convention: +N (Physical) Everywhere

**POLICY**: All kernel computation functions return **+N** (positive physical quantity). The sign flip to **-N** for the system matrix happens in **ONE place only**: `ComputeEntry()` in `rad_hacapk.cpp`.

```
Compute*BlockFast() returns +N (physical demagnetization tensor)
       ↓
GetInteractionMatrixElement() returns +N
       ↓
ComputeEntry(): A_val = -N_val + delta_ij * inv_chi[i]
       ↓
H-matrix stores system matrix A = -N + diag(1/chi)
```

**Sign convention applies uniformly to all supported interaction blocks.** No
DOF-type-specific sign conditionals.

| Layer | Sign | Description |
|-------|------|-------------|
| `Compute*BlockFast` | **+N** | Physical quantity (rad_interaction.cpp) |
| `m_flatInteractMatrix` | **+N** | Flat storage for LU/BiCGSTAB |
| `m_flat_N_data` | **+N** | HACApK pre-computed flat (rad_hacapk.cpp) |
| `ComputeEntry` callback | **-N+1/chi** | System matrix for H-matrix fill |
| LU/BiCGSTAB MatVec | **alpha=-1.0** | BLAS negates +N to get -N |
| `UpdateDiagonal` | **-diag_N+inv_chi** | Diagonal update after chi change |

**Do NOT** add sign flips in `GetCached*Element` or `GetCachedMixedElement`. These must return +N.

### 1/(4pi) Factor Convention

**POLICY**: Two field computation functions exist with DIFFERENT 1/(4pi) conventions. Do NOT mix them.

| Function | Location | 1/(4pi) | Use in |
|----------|----------|---------|--------|
| `RadFieldFromTriangleFaceWithBasis` | `rad_poly_analytical.cpp` | **Included** (in weight W) | `RadHACApKManager::Compute3x3BlockFast` |
| `FieldFromChargedTriangleLocal` | `rad_interaction.cpp` | **NOT included** | `radTInteraction::Compute3x3BlockFast`, `ComputeMixedBlockFast` |

- Functions using `FieldFromChargedTriangleLocal` MUST multiply by `RadConst::INV_FOUR_PI` at output
- Functions using `RadFieldFromTriangleFaceWithBasis` must NOT multiply again
- Both produce the same final result (+N with 1/(4pi)), but the intermediate values differ
- `B_comp()` (PreRelax mode) includes 1/(4pi) internally — do NOT scale its output

**Geometry indexing**: Tet, hex, and wedge all use **type-specific indices** (via `m_tetraElemIndices`, `m_hexaElemIndices`, `m_wedgeElemIndices`). Convert global element indices using `m_globalToTetraIdx` (O(1) lookup) or linear search for hex/wedge.

### Under-Relaxation for Nonlinear Problems

```python
rad.SolverConfig(relax_param=0.3)  # 30% damping (0.0 = full step)
rad.Solve(container, 0.0001, 1000, 1)
rad.SolverConfig(relax_param=0.0)  # Reset to full step
```

### Hantila Polarization Method

Hantila (1975) splits the constitutive relation into constant linear part + residual:

```
B = mu_0*(1+alpha)*H + mu_0*R    where R = M - alpha*H
```

Historically, Hantila polarization split the material law so the geometry-only
demag operator could be reused across nonlinear iterations:

```
H = H_ext + N*M
Substituting M = alpha*H + R:
(I - alpha*N)*H = H_ext + N*R    <- constant LHS, LU factored ONCE
```

**Advantages over Picard/Newton**:

| Feature | Picard (rad.Solve) | Newton | Hantila |
|---------|-------------------|--------|---------|
| Matrix factorization | Every iteration | Every iteration | **Once** |
| Jacobian needed | No | Yes (dM/dH) | **No** |
| BH curves | Yes | Yes | **Yes** |
| Hysteresis | No | No | **Yes** |
| Cost per iteration | O(N^3) LU | O(N^3) LU | **O(N^2) back-sub** |

**Current policy**: keep Hantila as solver-design background only. The
production soft-iron route is HDiv-VIM plus its current nonlinear iteration and
preconditioner stack; do not add a separate public Hantila demag backend.

**Usage**:
```python
from radia.hantila_solver import solve_hantila

# BH curve case
result = solve_hantila(iron_container, source=coil,
                       bh_data=BH_DATA, alpha=500.0, tol=1e-4)

# Hysteresis case (per-element material handles)
result = solve_hantila(iron_container, source=coil,
                       mat_handles=handles, alpha=500.0, relax=0.5)

# Result: M and H per element, convergence info
M = result['M']  # (n_elem, 3) in A/m
B = rad.Fld(iron_container, 'b', [0, 0, 0.05])  # Field evaluation
```

Reference: F.I. Hantila, Rev. Roum. Sci. Techn. - Electrotechn. et Energ., 1975.

---

## Compact HX Preconditioner (radia.sparsesolv_ngsolve)

Compact AMS/AMG/COCR types live in the `radia.sparsesolv_ngsolve` submodule.
Source: `src/ext/sparsesolv/` (monorepo integrated).
The `.pyd` is built by Build.ps1 via `add_ngsolve_python_module(sparsesolv_ngsolve ...)`
and shipped inside the radia wheel at `src/radia/sparsesolv_ngsolve.pyd`.

Import:
```python
import radia.sparsesolv_ngsolve as ssn
from radia.sparsesolv_ngsolve import (
    CompactAMSPreconditioner,
    ComplexCompactAMSPreconditioner,
    COCRSolver,
    SparseSolvSolver,
)
```

(History: 2026-05 the standalone `ngsolve-sparsesolv` PyPI package was
retired and absorbed into the Radia wheel. The top-level
`import sparsesolv_ngsolve` no longer resolves — always use
`radia.sparsesolv_ngsolve`. An earlier AGENTS.md draft described these
symbols as living in `ngsolve.la`; that integration was aspirational
and never landed.)

### Policy: Compact HX for HCurl Problems

**POLICY**: Use **Compact HX** (Compact Hiptmair-Xu) as the default preconditioner for HCurl curl-curl + mass systems. Compact HX is a HYPRE-free, TaskManager-native AMS implementation available via `radia.sparsesolv_ngsolve`.

**Name origin**: HX = Hiptmair-Xu (2007), "Nodal auxiliary space preconditioning in H(curl) and H(div) spaces", SIAM J. Numer. Anal. 45(6). "Compact" = lightweight, HYPRE-free, TaskManager-native.

**Configuration** (validated on complex eddy current @ 30 kHz, 155k-1.44M DOFs):

| Parameter | Value | Description |
|-----------|-------|-------------|
| Cycle type | 1 (01210) | pre-smooth, G-correct, Pi-correct, G-correct, post-smooth |
| Outer solver | BiCGStab | Non-symmetric Krylov solver |
| Fine smoother | l1-Jacobi | Fully parallel (TaskManager) |
| Subspace solver | CompactAMG | PMIS + classical interp + l1-Jacobi V-cycle |
| Pi mode | Separate Pix/Piy/Piz | Multiplicative correction |
| AMG theta | 0.25 | Strength-of-connection threshold |
| Correction weight | 1.0 | No damping |

**Performance** (mesh1_3.5T, 197k DOFs, BiCGStab, tol=1e-10):
- Compact HX + CompactAMG: 25 iterations (matches HYPRE AMS)
- HYPRE AMS + BoomerAMG: 25 iterations (reference)

**Source files** (`src/ext/sparsesolv/`):

| File | Description |
|------|-------------|
| `compact_amg.hpp` | Algebraic multigrid (PMIS, classical interp, l1-Jacobi) |
| `compact_ams.hpp` | AMS cycle (Pi, G subspace corrections, l1-Jacobi smoother) |
| `complex_compact_ams.hpp` | Complex Re/Im parallel wrapper (TaskManager) |

**Do NOT**:
- Add HYPRE dependency for new AMS features (use CompactAMG)
- Use sequential Gauss-Seidel in the fine smoother (breaks TaskManager parallelism)
- Use combined Pi with CompactAMG (combined Pi requires BoomerAMG num_functions=3)

**HYPRE option**: BoomerAMG subspace solver remains available behind `#ifdef SPARSESOLV_USE_HYPRE` for comparison benchmarks (subspace_solver=2).

### Shifted Preconditioner for Air+Conductor Problems

**POLICY**: For HCurl eddy current with air regions (σ=0), use **Shifted Preconditioner** instead of system regularization. Add ε·mass to the preconditioner only, not the system matrix.

```python
# Preconditioner: shifted (non-singular)
a_shifted += eps * nu * u * v * dx   # eps = 1e-6 * nu

# System: original (singular in air, but physically correct)
a += nu * curl(u) * curl(v) * dx
a += 1j * omega * sigma_cf * u * v * dx("cond")  # no eps here

# Solve original system with shifted preconditioner
c = Preconditioner(a_shifted, "bddc")
inv = CGSolver(a.mat, c.mat, ...)
```

**Verified**: eps value does NOT affect the solution (1e-4 to 1e-8 give identical ||B||²). Without shift: diverges (nan). See `src/ext/sparsesolv/examples/hiruma/shifted_ams_experiment.py`.

---

## IMA (Image Method of Analysis)

### IMA Sign Selection Policy

| Field vs Mirror Plane | IMA Sign |
|----------------------|----------|
| Field **parallel** to mirror | **+** (symmetric) |
| Field **perpendicular** to mirror | **-** (antisymmetric) |

```python
# Z-field, X-Z quarter model
rad.Solve(container, 0.0001, 100, 0, image='+x-z')  # Bz parallel to X-mirror, perp to Z-mirror

# X-field, X-Z quarter model
rad.Solve(container, 0.0001, 100, 0, image='-x+z')
```

### IMA Boundary Element Limitation

IMA produces incorrect results for **boundary elements** (faces ON symmetry plane) when observation points are **also on the symmetry plane** (~0.5x magnitude). Off-plane observation points work correctly (fixed 2026-02-04).

**Workaround**: Use explicit element duplication for models with boundary elements.

**When IMA is Safe**:
1. Non-boundary elements only (offset from symmetry planes)
2. Observation points off-plane
3. Explicit full-model solve -- use when the symmetry cut is not covered by IMA

---

## PEEC & Conductor Solver

### Architecture Overview

**Approach**: PEEC (Partial Element Equivalent Circuit) with SIBC (Surface Impedance Boundary Condition) and ESIM (Effective Surface Impedance Method).

**Target**: Induction heating (1-500 kHz), WPT (6.78/13.56 MHz), power electronics (DC-1 MHz).

See `docs/` for detailed PEEC documentation.

### Filament-Panel Architecture (FastImp Style)

```
Surface Mesh -> Face -> Panel (Star: charge)  -> P matrix
             -> Edge -> Filament (Loop: current) -> L, R matrices
```

Loop-Star basis transformation is NOT needed -- filaments and panels are inherently separate in PEEC.

### PEEC System Equation

```
[R + jwL + Zs    jwM_LS  ] [I_filament]   [V]
[jwM_LS^T        P/(jw)  ] [Q_panel   ] = [0]
```

### Node-Segment Topology API

```python
from peec_matrices import PyPEECBuilder
from peec_topology import PEECCircuitSolver

builder = PyPEECBuilder()
n1 = builder.add_node_at(0, 0, 0)
n2 = builder.add_node_at(0.1, 0, 0)
builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
builder.add_port(n1, n2)
topo = builder.build_topology()

solver = PEECCircuitSolver(topo)
Z = solver.compute_port_impedance(freq=1e6)
```

### Multi-Filament (nwinc/nhinc)

Use `nwinc`/`nhinc` parameters to subdivide conductor cross-sections for skin/proximity effect:
```python
builder.add_connected_segment(n1, n2, 3e-3, 3e-3, sigma=5.8e7, nwinc=3, nhinc=3)
```

Guidelines: DC=1x1, moderate skin (d/delta~2-5)=3x3, strong skin (d/delta>5)=5x5+.

### FastHenry .inp Parser

```python
from fasthenry_parser import FastHenryParser
parser = FastHenryParser()
parser.parse_file('inductor.inp')
result = parser.solve()
```

Supports: `.Units`, `N`/`E` definitions, `.external`, `.freq`, `.default`, `.equiv`, `.magnetic` blocks, line continuation `+`.

### Coupled PEEC + HDiv-VIM

```python
from peec_coupled import CoupledPEECSolver
solver = CoupledPEECSolver(topology_dict, magnetic_objects=[core_id])
solver.compute_coupling_matrix()  # N_seg Radia Solve calls
Z = solver.compute_port_impedance(freq)
Z_sweep = solver.frequency_sweep(freqs)
L_total = solver.get_effective_inductance()  # L_air + Delta_L
```

For linear materials, `Delta_L` is frequency-independent (computed once).

**Physics**: `Z_eff(f) = diag(R + Zs(f)) + jw * (L_air + Delta_L)` where Delta_L comes from Biot-Savart -> `rad.ObjBckg()` + `rad.Solve()` -> vector potential A -> mutual inductance.

**FastHenry .magnetic Block** for coupled simulations:
```
.magnetic
  type=box
  center=0.05,0.01,0.0
  size=0.06,0.01,0.01
  mu_r=1000
.endmagnetic
```

### SIBC Implementations

| Conductor Type | Method | Library |
|---------------|--------|---------|
| Circular | Bessel I0/I1 | `scipy.special.iv` |
| Rectangular (d << w) | Dowell formula | C++ rad_peec_surface_impedance.cpp |
| Nonlinear magnetic | ESIM cell problem | `esim_cell_problem.py` |

**Bessel**: Use `scipy.special.iv` (modified Bessel), NOT `jv` (regular Bessel). MKL does not provide Bessel functions.

### ESIM (Effective Surface Impedance Method)

ESIM solves 1D cell problem for H-dependent surface impedance: `d/dz[(1/mu(z)) * dH/dz] = jw*sigma*H`.

Supports complex permeability: `mu = mu' - j*mu"` for magnetic hysteresis/grain eddy current losses.

Use for: induction heating workpieces, nonlinear iron cores, lossy ferrite at high frequency.

Reference: `src/radia/esim_cell_problem.py`, `src/radia/esim_coupled_solver.py`.

### Deleted Legacy PEEC APIs

The following C++ APIs are **REMOVED**: `CndLoop`, `CndRecBlock`, `CndLoopFromHelix`, `CplMagCreate`, `CplMagSolve`, `CplMagSetFrequency`, `CndHexahedron`, `CndWire`, `CndSpiral`, `MatSIBC`. Use `PEECBuilder` and `CoupledPEECSolver` instead.

### PRIMA Model Order Reduction

**POLICY**: Use PRIMA (not CLN/Cauer) terminology. Both use Lanczos tridiagonalization; PRIMA (1998, IEEE TCAD) is the standard reference.

Key classes: `SPICEExtractionConfig`, `PRIMASchurExtractor`, `LoopStarMagneticCoupled` in `lanczos_reduction.py`.

### ngsolve.bem Integration

Radia PEEC works alongside ngsolve.bem:

| Range | Solver |
|-------|--------|
| DC - 1 MHz | Radia PEEC + SIBC |
| DC - 1 MHz | ngsolve.bem (Weggler EFIE, low-freq stable) |
| 1 MHz - GHz | ngsolve.bem (Helmholtz) |

Radia PEEC unique features: direct circuit extraction (L, R, C), native SPICE netlist, Lanczos MOR, and HDiv-VIM core coupling.

### Integration Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│  NGSolve Geometry & Mesh                                        │
└──────────────────────────┬──────────────────────────────────────┘
                           │
           ┌───────────────┴───────────────┐
           ▼                               ▼
┌─────────────────────┐         ┌─────────────────────┐
│  Radia PEEC         │         │  ngsolve.bem        │
│  - Loop-Star        │         │  - EFIE/MFIE        │
│  - SIBC/ESIM        │  <--->  │  - H-matrix         │
│  - Lanczos MOR      │ coupling│  - Helmholtz/Laplace│
│  - SPICE output     │         │  - Low-freq Weggler │
└─────────────────────┘         └─────────────────────┘
           │                               │
           └───────────────┬───────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│  Unified Solution (NGSolve GridFunction)                        │
└─────────────────────────────────────────────────────────────────┘
```

### PEEC Source Files

**C++ Core**:
- `src/core/rad_peec_matrices.h/cpp` -- PEECSegment, PEECPort, PEECMatrices, MutualInductance
- `src/core/rad_peec_surface_impedance.cpp` -- Dowell formula
- `src/lib/rad_peec_matrices_api.cpp` -- pybind11 bindings

**Python**:
- `src/radia/peec_topology.py` -- PEECCircuitSolver (MNA nodal admittance)
- `src/radia/peec_coupled.py` -- CoupledPEECSolver
- `src/radia/fasthenry_parser.py` -- FastHenry .inp parser
- `src/radia/esim_cell_problem.py` -- ESIM cell problem solver
- `src/radia/lanczos_reduction.py` -- PRIMA model order reduction

---

## Material Specification

### MatLin - Linear Materials

```python
mat = rad.MatLin(mu_r)                       # Isotropic (preferred)
mat = rad.MatLin([mu_r_par, mu_r_perp], [ex, ey, ez])  # Anisotropic
```

For isotropic materials, ALWAYS use single-argument form. MatLin is for soft magnetic materials only.

### MatSatIsoTab - Nonlinear (B-H Curve)

```python
BH_DATA = [[0.0, 0.0], [100.0, 0.1], [1000.0, 1.2], [50000.0, 2.0]]
mat = rad.MatSatIsoTab(BH_DATA)  # [[H(A/m), B(T)], ...]
```

### Permanent Magnets

For fixed magnetization PM, specify directly -- no `Solve()` needed:
```python
pm = rad.ObjHexahedron(vertices, [0, 0, 954930])  # M in A/m
B = rad.Fld(pm, 'b', [0, 0, 0.1])
```

Call `Solve()` only when soft iron is present alongside permanent magnets.

PM material classes (`MatMagFixed`, `MatMagLinear`, `MatMagCurve`) are available but currently all behave as fixed magnetization. Full demagnetization is planned.

See `docs/ELF_CONVENTIONS.md` for detailed unit system documentation.

### Hysteresis Materials (Play and Energy Models)

Two B-input play hysteresis models are available. The Play model is recommended (faster, no sign constraints).

```python
# Play model (recommended): B-input, direct Forward O(K)
from radia.hysteresis_io import load_hys_file
K, eta, f_k_tables = load_hys_file('material.hys')
mat = rad.MatPlayHysteresis(K, eta, f_k_tables)
# K: number of play operators
# eta: ndarray[K], play thresholds in Tesla
# f_k_tables: list of (r_array, f_array) tuples (shape functions)

# Energy model: B-input, Egger Schur complement Newton
mat = rad.MatEnergyHysteresis(K, eta, f_k_tables, eps=1e-6)
# Same parameters + eps convergence tolerance
# Requires non-negative, monotonically increasing shape functions (convex U_k)
```

**State management** (works for both Energy and Play models):
```python
rad.MatApl(iron, mat)
# ... solve quasi-static step ...
state = rad.MatHysSaveState(mat)     # Save state (ndarray, length K*9)
rad.MatHysRestoreState(mat, state)   # Restore state
rad.MatHysCommitState(mat)           # Commit converged state for next step
```

**Play vs Energy model comparison**:

| Feature | Play Model | Energy Model |
|---------|-----------|--------------|
| Forward (B->H) | O(K) direct | Newton (100 iter) |
| Inverse (H->B) | Newton + analytical Jacobian | K independent Newton |
| Shape functions | No sign constraint (negative OK) | Must be non-negative |
| Speed | 4-9 us/eval (Forward) | 100-500 us/eval |

**MatMvsH** - Query M(H) for any material:
```python
M = rad.MatMvsH(mat, [Hx, Hy, Hz])  # Returns [Mx, My, Mz] in A/m
```

### Permanent Magnet + Soft Iron Interaction

When combining PM with soft iron, use `Solve()`:
```python
pm = rad.ObjHexahedron(pm_vertices, [0, 0, 954930])  # Fixed PM
iron = rad.ObjHexahedron(iron_vertices, [0, 0, 0])    # Zero initial M
mat_iron = rad.MatLin(1000)
rad.MatApl(iron, mat_iron)
assembly = rad.ObjCnt([pm, iron])
result = rad.Solve(assembly, 0.0001, 1000, 0)  # LU solver
B = rad.Fld(assembly, 'b', [0, 0, 0.1])
```

---

## File & Naming Conventions

### Python Script Path Import

```python
# Prefer installed/editable package imports.  Only add a repo-local src path
# for tests/validation/docs helpers that must run before installation.
from pathlib import Path
repo = next(p for p in Path(__file__).resolve().parents if (p / "src" / "radia").exists())
sys.path.insert(0, str(repo / "src"))
```

Import from the `radia` package or the repo-local `src` tree (not build
directories).  Do not add new `examples/` import patterns.

### Script Naming Convention

Use **snake_case** with functional prefixes:

| Prefix | Purpose | Example |
|--------|---------|---------|
| `demo_` | Educational demonstration | `demo_batch_evaluation.py` |
| `benchmark_` | Performance measurement | `benchmark_solver_scaling.py` |
| `verify_` | Correctness verification | `verify_curl_A_equals_B.py` |
| `compare_` | Method comparison | `compare_radia_ngsolve.py` |
| (none) | Physical model name | `sphere_in_quadrupole.py` |

### VTK Export


### Error Display Policy: Scientific Notation

**POLICY**: Volume error and area error MUST use **scientific notation** (e.g., `-8.24e-02%`, `+2.35e-04%`). Do NOT use fixed-point notation like `-0.08%` or `+0.0002%` — this hides significant digits at small values and makes p-convergence trends unreadable.

Apply this to: Mesh Evaluation tables, Joachim correspondence, benchmark results, test output.

### Benchmark Policy

**POLICY (hibino-first; mdx CI-first, 2026-09-03)**: Run solver-heavy
validation, optimization, scaling, memory, and timing work on hibino first.
Use mdx only when hibino is unavailable and both the mdx CI runner and its job
queue are idle. Compute work must never delay or destabilize CI/preflight.
- Check the selected host's active processes and load before launch. For an mdx
  fallback, also check the CI runner/service and queued or active jobs.
- If hibino is unavailable and mdx is busy, LAB may run correctness/smoke probes
  only; defer publication or decision-grade timing until a permitted compute
  host is idle.
- Record hostname, start time, runtime, memory conditions, and runtime versions
  in the result JSON or log. Historical mdx measurements remain valid
  provenance, but they do not define the current routing policy.
- Keep this policy synchronized between AGENTS.md and CLAUDE.md.

**POLICY**: 全てのベンチマークスクリプト (`bench_*.py`) は JSON 形式の結果ファイルを出力すること。

**実行ルール**:
1. 1ケース毎に実行（並列実行しない、メモリ測定の正確性のため）
2. メモリ使用量を記録（`psutil` を使用、`tracemalloc` はC++メモリを追跡しない）
3. 結果JSONファイル名: `results_{benchmark_name}.json` (スクリプトと同じディレクトリ)

**JSON必須フィールド**:

| フィールド | 型 | 説明 |
|-----------|------|------|
| `peak_memory_mb` | `float` | ピークメモリ使用量 (psutil peak_wset/rss) |
| `t_setup` | `float` | 前処理セットアップ時間 (秒) |
| `t_solve` | `float` | 線形ソルバー実行時間 (秒) |
| `iterations` | `int` | 反復回数 |
| `converged` | `bool` | 収束判定 |

**JSONメタデータ** (トップレベル):

| フィールド | 型 | 説明 |
|-----------|------|------|
| `timestamp` | `str` | ISO 8601形式 |
| `hostname` | `str` | `platform.node()` |
| `benchmark` | `str` | ベンチマーク名 |
| `problem` | `dict` | 問題パラメータ (ndof, ne, order等) |
| `results` | `list[dict]` | 各ケースの結果 |

```python
import json, os, platform, psutil, time
from datetime import datetime

def get_peak_memory_mb():
    mem = psutil.Process(os.getpid()).memory_info()
    return mem.peak_wset / (1024 * 1024) if hasattr(mem, 'peak_wset') else mem.rss / (1024 * 1024)

def save_benchmark_results(filename, benchmark_name, problem, results):
    data = {
        "timestamp": datetime.now().isoformat(),
        "hostname": platform.node(),
        "benchmark": benchmark_name,
        "problem": problem,
        "results": results,
    }
    with open(filename, "w") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"Results saved to {filename}")
```

---

## .vol Pipeline: 2-Path Generation, 1-Path Computation

### File Management: Reproducible CAD and Explicit Exports

**POLICY**: A loaded Cubit model exports directly to a user-selected path. The
toolbar never forces a journal save and never couples output names to a `.jou`
basename. A known source journal may supply an editable initial directory and
basename only. The toolbar does not save `.cub5` automatically.

| File | Role |
|------|------|
| `.sat` / `.step` / `.jou` | Optional CAD sources. SAT preserves Cubit's ACIS-native geometry, STEP is the portable exchange format, and tracked journals are replayable meshing/CAD recipes. |
| `.vol` | **Computation interface**. Sole solver-mesh interface between Cubit and NGSolve. No Cubit ABI dependency. |

**Design**:
- A loaded Cubit session can export without a current `.jou`.
- If no model is loaded, the toolbar may offer to replay a `.jou`; this is a
  model-loading convenience, not an output-path contract.
- If Cubit reports a current journal, its path is only an initial editable hint.
- The export dialog owns the output directory and filename and remembers the
  per-format directory. Cancelling the dialog performs no save or export.
- Journal capture is an explicit reproducibility action. It is never an export
  side effect and temporary journals do not belong at the repository root.
- `.vol` is the branch point into NGSolve/Radia computation; `.msh` is the GMSH
  post-processing artifact where required.

### Cubit/NGSolve Complete Separation Policy

**POLICY**: Radia-NGSolve computation scripts (`calc_*.py`, panels) must **NEVER `import cubit`**. The `.vol` file is the **sole interface** between Cubit and NGSolve.

**Why**: Coreform Cubit is expensive commercial software (annual license). NGSolve/Radia computation must work without Cubit. `.vol` files can also be generated by Netgen standalone (STEP -> Netgen -> `.vol`), so the computation pipeline must not assume Cubit is available.

**Mesh export is C++ only** (`export netgen` APREPRO command):
- Cubit -> `export netgen "model.vol" order N` -> `.vol` with labels + curving

**1-Path Computation** (`.vol` only, no Cubit dependency):
```python
# calc_*.py — NGSolve computation script
# Accepts --vol only. NEVER imports cubit.
from ngsolve import Mesh, ...
mesh = Mesh("model.vol")   # labels + curving loaded from .vol
# ... FEM solve, post-processing ...
```

**calc_*.py accepts `--vol` only** (no `.cub5`):
- `calc_volume.py --vol model.vol` — volume/area integration
- `calc_peec.py --step coil.step` — PEEC filament inductance (no mesh needed for coil)
- `calc_fem_kelvin.py --vol model.vol` — FEM Kelvin + SIBC (IH workpiece)
- `calc_verify_vol.py --vol model.vol` — consistency check vs companion JSON
- `calc_mesh_eval.py --vol-base model` — p-convergence (`_p1.vol` ... `_p5.vol`, C++ exports)

### IH: No Source/Sink Sidesets Required (PEEC+FEM)

**POLICY**: The IH production path (PEEC+FEM) does NOT require source/sink sidesets.
Coil current is defined by filament topology (STEP -> centerline -> PEECBuilder ports).
FEM workpiece uses Biot-Savart H_s from PEEC filaments as the incident field +
SIBC Robin BC on the workpiece surface (no current injection faces needed).

**IH has 2 paths** (both source/sink-free):
1. **PEEC+FEM** (production): STEP -> filaments -> PEEC coil L,R + workpiece FEM-SIBC+Kelvin
2. **FEM** (reference): full volume mesh with coil included + Kelvin + SIBC/ESIM

Note: Omega-reduced H_s formulation is for the **accelerator magnet panel**, not IH.
IH uses Biot-Savart from filaments (PEEC path) or volume mesh coil (FEM path).

**Workpiece .vol needs only material blocks**:
- `workpiece` (sigma, mu_r)
- `air`
- `kelvin`

No sidesets needed. No source/sink labels.

**BEM (legacy)**: reusable BEM solver modules live in `src/radia` as
`radia.bem_inductance`, `radia.bem_coupled_solver`, and `radia.ngsbem_*`.
Executable reference scripts and sweep results live under
`validation_test/induction_heating/bem_reference/`. BEM knowledge is in
`mcp-server-radia-ngsolve` (ngsbem_inductance topic).

**References**:
- Djordjevic & Notaros, "Double higher order MoM", IEEE TAP 2004 (geometry/basis independence)
- Marussig et al., "Fast Isogeometric BEM based on Independent Field Approximation", arXiv 2014
- Dolz et al., "Bembel: Fast Isogeometric BEM", arXiv 2019

**`.vol` Must Be Self-Contained**:
- Material labels: `SetMaterial()` -> `materials` section
- Boundary labels: `SetBCName()` -> `bcnames` section
- High-order curving: `curvedelements` text section (upstream Netgen master feature)
- No external STEP/geometry file needed for computation

**Cubit Plugin Responsibility**: The `export netgen` C++ command handles all label + curving embedding into `.vol`. Higher maintenance cost is acceptable for complete separation.

---

## Cubit Panel Architecture

### 4-Layer Architecture

**POLICY**: Cubit, Simulink, and Radia-NGSolve computation are separate
processes. They exchange checked files and artifacts; no interface process
imports another process's private runtime.

```
┌─────────────────────────────────────────────────────────────────┐
│  Layer 1: C++ .ccm APREPRO backend (no Qt)                      │
│  ─────────────────────────────────────────────────────────────  │
│  Export Mesh commands (GMSH/Nastran/VTK/Netgen Vol/FEMEEM/MEG)  │
│  Mesh evaluation is a docs/notebook workflow, not a toolbar item │
│  Explicit output path; no implicit .jou/.cub5 save              │
│  export netgen/gmsh/jmag_nastran/vtk (APREPRO commands)         │
└─────────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────────┐
│  Layer 2: Cubit GUI Python (Python 3.10 + PySide6, Cubit 2025.12)│
│  ─────────────────────────────────────────────────────────────  │
│  register_toolbar.py -> Solve menu management                   │
│    Radia-NGSolve / Generate Coil / Kelvin / Reload / Verify     │
│  import cubit OK (same process). import radia/ngsolve FORBIDDEN │
│  Exports checked .vol/.sol assets for Layers 3 and 4            │
└─────────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────────┐
│  Layer 3: Radia Simulink application blocks (MATLAB/Simulink)   │
│  ─────────────────────────────────────────────────────────────  │
│  matlab/radia_simulink_library.slx + radia.simulink adapters    │
│  Explicit trigger launches Layer 4 and reads result/log files.  │
│  Native MEX/ROM is optional and admitted only after parity tests.│
└─────────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────────┐
│  Layer 4: Computation (Python 3.12, no GUI)                     │
│  ─────────────────────────────────────────────────────────────  │
│  calc_peec.py (PEEC coil) / calc_fem_kelvin.py (FEM workpiece)  │
│  import cubit FORBIDDEN. import PySide6 FORBIDDEN.              │
│  NGSolve + GMSH only. .vol input, JSON stdout output.           │
└─────────────────────────────────────────────────────────────────┘
```

**Interfaces between layers**: self-contained `.vol` / `.sol` inputs plus
versioned configuration, `run.log`, solver JSON, and `result.json` artifacts.

**Filename convention**: The user or application run owns the output path.
Related artifacts use the selected/run basename (`{base}.vol`, `{base}.msh`,
`{base}_J.sol`, `{base}_q.sol`) without requiring a same-named journal.

### Layer Isolation Rules

| Rule | Rationale |
|------|-----------|
| Layer 4 must NOT `import cubit` | Cubit is expensive commercial software. Computation must work without it. |
| Layer 4 must NOT import PySide6/PyQt5 | Headless computation. JSON stdout only. |
| Layer 3 must NOT load Cubit's Python runtime | Simulink is a separate MATLAB process; it uses files and the validated headless runner. |
| Layer 2 must NOT `import radia` or `import ngsolve` | DLL conflicts (Cubit bundles its own numpy/scipy). |
| Layer 1 (C++) has no Python/Qt dependency | `.ccm` links the Cubit C++ API directly. The Qt5 `.ccl` GUI was removed in radia 4.80.0. |

### Application Interface Files

| File | Layer | Purpose |
|------|-------|---------|
| `src/radia/panels/radia_export_menu.py` | 2 (PySide6) | Export Mesh menu + dialogs inside Cubit |
| `src/radia/panels/register_toolbar.py` | 2 (Cubit Python) | Solve/export menu integration; no notebook launcher |
| `matlab/radia_simulink_library.slx` | 3 (Simulink) | Final EM, PCB, Motor, Stream Function, and IH blocks |
| `matlab/+radia/+simulink/` | 3 (MATLAB) | Block builders, explicit runner, masks, and runtime adapters |
| `src/radia/simulink/application.py` | 3→4 bridge | Python-backed explicit-trigger artifact runner |
| `ih_design.py` | 4 contract | Shared UI-neutral IH settings and command mapping |
| `*_design.py` | 4 contract | UI-neutral settings-to-command mapping shared by blocks and AI |
| `src/radia/panels/calc_*.py` | 4 (no GUI) | Headless application computations |
| `src/radia/panels/calc_mesh_eval.py` | 4 (no GUI) | p-convergence + format QA, called from docs/notebooks |

### Cubit Plugin: C++ First, No Python ABI Dependency

**POLICY**: Cubit plugin functionality MUST be implemented in C++ to avoid Python ABI mismatch. Cubit embeds Python 3.10; NGSolve/Radia use Python 3.12. Sharing Python objects between them causes segfaults and DLL conflicts.

- `.ccm`: Link Cubit C++ API (cubiti, cubit_util) directly -- no Python or Qt dependency
- `cubit_mesh_curver.pyd`: pybind11 for Python 3.12 -- does NOT link Cubit C++ libraries
- Netgen `SetNCD2Names()` is not exposed to Python -- call from C++ side in `NetgenCurverPure`
- Interface between Cubit and NGSolve: **.vol file** (text format, no ABI dependency)
- Export Mesh computation/export stays in the C++ `.ccm` backend, while the
  supported Cubit-facing GUI is the Python/PySide6 toolbar/menu. Do not
  reintroduce the Qt5 `.ccl` / `RadiaComp.cpp` path.

---

## Visualization Policy

### NGSolve-Through Principle

**POLICY**: NGSolve を経由すれば済む問題は、NGSolve を経由するワークフローとする。Radia 独自のポスト処理機能は作らない。

- **幾何形状**: STEP 出力 → GMSH GUI で読み込み
- **磁場ポスト**: NGSolve GridFunction → GmshPostExport (.msh v4.1) → GMSH GUI
- **rad.Fld()**: デバッグ・確認用（ポストには使わない）

**Why**: NGSolve 経由なら構造格子に限定されない。非構造メッシュ上で高次要素の精度を保ったまま場を評価・出力できる。

| 用途 | ツール | 出力 |
|------|--------|------|
| **幾何形状** | CoilBuilder.write_step(), OCC shapes | **.step** → GMSH |
| **磁場ポスト** | NGSolve GridFunction → **GmshPostExport** | **.msh v4.1** → GMSH |
| 確認用 | `rad.Fld()` (点評価) | スクリプト内のみ |

**Do NOT** implement custom visualization in Radia C++ code.

**Removed APIs**: `rad.ObjDrwVTK()`, `exportGeometryToVTK()`, `radia_pyvista_viewer.py`.

### Standard Output Format: GMSH .msh v4.1

GMSH を可視化ツールとして使用する理由:
- **高次要素ネイティブ対応** (Tri6, Tri10, Tri15, ..., arbitrary p)
- **STEP ファイルを直接読み込み** → 幾何形状と磁場を重ねて表示
- Per-material Physical Groups で材料別表示
- **.msh v4.1 only** (lab-wide standard; netgen I/O は常に .vol 経由)

### GMSH Invocation Policy: Python API Only

**POLICY**: GMSH は **常に pip-gmsh の Python API 経由** で呼ぶ. 単独
``gmsh.exe`` を探したり (CST 同梱のものを含む) 別 install するスクリプトを書かない.

理由:
- pip-gmsh (PyPI ``gmsh``) は OCC + FLTK GUI を含む完全な Windows
  バイナリ. `gmsh.fltk.run()` は blocking で window が安定して開く.
- 仮に Python プロセスから起動した GMSH window が「すぐ消える」と
  いう報告があっても、それは **呼び出し側で process kill / timeout
  をしている** ことが原因 (実例: 2026-05-02 Codex が `timeout 4`
  + `taskkill` で background process を殺し、ユーザに「GUI が瞬殺
  される」と誤って報告).  pip-gmsh GUI 自体は問題ない.
- ``C:\Program Files\CST Studio Suite ...\gmsh.exe`` のような他社
  バンドル版を呼ばない. version 不一致 (CST: 4.11.1, pip: 4.15+) や
  別途 install を強制するワークフローは避ける.

```python
# CORRECT
import gmsh
gmsh.initialize()
gmsh.open(msh_path)
# ... display options ...
gmsh.fltk.run()       # blocking GUI; user closes the window to exit
gmsh.finalize()
```

```python
# WRONG -- do not do this
subprocess.Popen([r"C:\Program Files\CST...\gmsh.exe", msh_path])  # external binary
subprocess.Popen([r"C:\Tools\gmsh\gmsh.exe", msh_path])             # ad-hoc install path
```

If `fltk.run()` appears to flash and exit on a user's machine, debug:
1. Confirm the script reaches `gmsh.fltk.run()` (add a print before).
2. Check that no caller wraps the launch in a `timeout` or kills the
   Python process. background launches (`run_in_background=True`) MUST
   NOT be later `taskkill`'d if you want the user to keep interacting
   with the window.
3. If the GUI genuinely doesn't display, check pip-gmsh installation
   (`pip show gmsh`) -- a corrupted install, not a missing executable,
   is the failure mode.

### GmshPostExport: High-Order Field Visualization

```python
from radia.gmsh_post_export import GmshPostExport

# BEM/FEM surface visualization (arbitrary order curved elements)
post = GmshPostExport(mesh, boundary=True)  # boundary=True for BND from volume mesh
post.add_field("|J|", node_J, ncomp=1)      # per-vertex scalar
post.add_vector_field("J", gf_J)            # vector field
post.write("results.msh")
# -> GMSH renders Tri6/Tri10/Tri15/... with correct curved interpolation
```

**Key features**:
- `boundary=True`: exports BND surface elements from a volume mesh (BEM use case)
- **Arbitrary order** support: Curve(p) → Tri type auto-selected (p=2: Tri6, p=3: Tri10, p=4: Tri15, p=5: Tri21)
- High-order nodes extracted via **GetTrafo + GMSH reference coordinates** (exact curved positions)
- Per-material Physical Groups for selective field display in GMSH GUI
- NodeData and ElementData support
- 2-phase output: mesh first (before solve), field added after solve

**Supported GMSH triangle types**:

| Order | GMSH Type | Nodes | Nodes/edge | Interior |
|-------|-----------|-------|------------|----------|
| 1 | 2 (Tri3) | 3 | 0 | 0 |
| 2 | 9 (Tri6) | 6 | 1 | 0 |
| 3 | 21 (Tri10) | 10 | 2 | 1 |
| 4 | 23 (Tri15) | 15 | 3 | 3 |
| 5 | 25 (Tri21) | 21 | 4 | 6 |

**Implementation note**: High-order node positions are extracted via `mesh.GetTrafo(el)` evaluated at GMSH reference coordinates (obtained from `gmsh.model.mesh.getElementProperties()`). Each BND element's transformation is evaluated at equidistant reference points matching GMSH's Lagrange node layout. Edge nodes are cached across shared edges with direction correction. H1 order=p GridFunction approach was found **unreliable for p>=4** (`Set()` L2 projection + averaging corrupts coordinates).

**GMSH display setting**: `Mesh.NumSubEdges = 4` required to render curved surfaces (default=1 draws straight lines). Set via GMSH console: `Mesh.NumSubEdges = 4;`

### Visualization Workflow

```
GMSH GUI:
  Merge "coil.step"        ← CoilBuilder.write_step()
  Merge "magnet.step"      ← OCC shape → STEP
  Merge "field.msh"        ← NGSolve → GmshPostExport (.msh v4.1)
  → 幾何形状 + 磁場を重ねて可視化
```

**Design principle**: Radia C++ に可視化コードを持たない。NGSolve + GMSH に任せる。少人数で最大の成果を出すため。

**Panel-specific split**: when the visualization is part of a
human-facing panel/notebook workflow, prefer `netgen.webgui`; when the
consumer is an LLM or a headless validation run, prefer durable
GMSH `.msh v4.1` artifacts.

For notebook IO, `.vol` and `.sol` double-click behavior should be the
plain Netgen viewer.  If double-click fails, remember that pip Netgen's
argument dispatch is in `netgen.__main__` and may need a `.vol`/`.sol`
file handler; `netgen.__init__` is only the package startup/DLL setup
side.

---

## Universal Relaxation Network (URN)

URN material now lives under `docs/universal_relaxation_network/` for the
showcase / paper / result artifacts and `src/radia/urn/` for reusable API code.
Do not recreate `examples/universal_relaxation_network/`.

**Policy**:
- Synthetic data MUST be clearly marked as synthetic
- Real-world datasets MUST include license and citation info
- All paper results reproducible from scripts in this directory

---

---

## Cubit Mesh Export Module

**POLICY**: Export Mesh computation is **C++ backend first**. Mesh export
functionality lives in the Cubit `.ccm` command plugin
(`cubit_mesh_export.ccm`), while the supported Cubit-facing GUI is the
Python/PySide6 toolbar/menu (`src/radia/panels/radia_export_menu.py`).
Do NOT reintroduce the retired Qt5 `.ccl` / `RadiaComp.cpp` GUI path.

### C++ Plugin Architecture

| Component | File | Purpose |
|-----------|------|---------|
| `.ccm` (plugins/) | `cubit_mesh_export.ccm` | APREPRO commands: `export gmsh/netgen/vtk/femeem/meg`, `export jmag_nastran` |
| PySide6 toolbar | `src/radia/panels/radia_export_menu.py` | Export Mesh menu + dialog inside Cubit |
| `.pyd` (Python 3.12) | `cubit_mesh_curver.pyd` | pybind11: Cubit-free mesh curving |

**Export formats** (all in C++, ACIS geometry projection for curving):

| Format | Command | Max Order | Notes |
|--------|---------|-----------|-------|
| Netgen Vol | `export netgen "f.vol" order 3` | 1-5 | Primary format for NGSolve FEM |
| GMSH v4.1 | `export gmsh "f.msh"`           | 1-3 | Lab-wide standard; structured entity blocks |
| Nastran BDF | `export jmag_nastran "f.bdf"` | 1-2 | CTETRA/CTETRA(10), nopyramid option |
| VTK | `export vtk "f.vtk"` | 1-2 | Legacy format, cell types 10/24 |

**GMSH order limit**: Order 4-5 is an error (not fallback). NetgenCurver face/volume
interior node extraction is unreliable at p>=4 (linear interpolation fallback causes
negative Jacobians in GMSH). Use `export netgen` for order 4-5.

**High-order mesh curving** (order >= 2):
- `NetgenCurver` (compact_netgen, static link): `CallbackGeometry` + `BuildCurvedElements` for order 1-5
- **No fallback**: `HighOrderMesh` is removed. NetgenCurver failure = error.
- ACIS surface projection via `closest_point_uv_guess` (UV-guided Newton, falls back to `closest_point_trimmed`)
- Surface elements: `parse_cubit_list("tri/face") + get_connectivity` (shared FEM nodes)
- Segment elements on curves: `parse_cubit_list("edge", "in curve N")` for PointBetweenEdge
- Edge projection callback: `RefEdge::get_curve_ptr()->closest_point_trimmed()`
- Requires `cubit_geom.dll` (DELAYLOAD)
- **ACIS is NOT thread-safe**: OpenMP parallelization of BuildCurvedElements is impossible

### Mesh Export Policy

**POLICY**: Mesh export uses `export netgen` C++ command only. Pure Python reference (`cub5_to_vol.py`) is maintained in the netgen fork, not in Radia. Run `test_vol_multi_geometry.py` (10 shapes) after any NetgenCurver change.

### Companion JSON (.vol.json)

`export netgen` writes a companion JSON alongside the .vol file:
```json
{
  "materials": {"sphere": 5.235988e-04},
  "boundaries": {"surface_1": 3.141593e-02},
  "edges": {"curve_1": 3.141593e-01},
  "n_elements": 10359, "n_points": 2071, "order": 3
}
```
- CAD volume per material (RefVolume::measure)
- CAD area per boundary (RefFace::area)
- CAD length per edge (RefEdge::measure)
- calc_verify_vol.py reads .vol.json for consistency checks

### Source Files

| File | Purpose |
|------|---------|
| `src/cubit_plugin/ExportGmshCommand.cpp` | GMSH v4.1 writer |
| `src/cubit_plugin/ExportNastranCommand.cpp` | Nastran BDF writer |
| `src/cubit_plugin/ExportVtkCommand.cpp` | VTK Legacy writer |
| `src/cubit_plugin/ExportNetgenCommand.cpp` | Netgen .vol writer + companion JSON |
| `src/cubit_plugin/MeshData.cpp` | Shared mesh extraction from Cubit |
| `src/cubit_plugin/NetgenCurver.cpp` | Order 1-5 curving via compact_netgen |
| `src/cubit_plugin/callbackgeom.cpp` | ACIS projection callbacks for CallbackGeometry |
| `packages/cubit-mesh-export/src/cubit_mesh_export/cubit_mesh_export.ccm` | Packaged Cubit APREPRO command plugin |
| `packages/cubit-mesh-export/src/cubit_mesh_export/cubit_mesh_curver.pyd` | Packaged Python 3.12 curving helper |

### Build

Compact Netgen is statically linked — no nglib.dll/ngcore.dll needed for the Cubit plugin.

### Testing

```bash
python tests/cubit/test_export_combinations.py   # All format x option combinations
python tests/cubit/test_ho_volume_all_formats.py  # Order=2 volume accuracy (sphere)
```

---

## Agent Division of Labor: Claude Stops at the Commit; Codex Owns Everything After (2026-06-21, scope widened 2026-08-28)

**POLICY**: Split of work between AI agents.

**Scope**: this applies to **EVERY Claude session — every project, every
worktree, every concurrent session** — not only the session that happens to be
reading this file. The same policy is kept machine-wide in the user-level
`~/.claude/CLAUDE.md`; keep this section synchronized verbatim with the
corresponding section in the other agent-policy file in the same change.

- **Claude's responsibility ENDS at the local `git commit`.** Claude implements,
  tests locally, and commits (this-session files BY NAME, per the existing commit
  hygiene), then STOPS and reports the commit SHA and branch.
- **Everything after the commit is codex's job — NOT Claude's.** codex pushes,
  watches CI to green, fixes CI infrastructure, cuts tags, runs the PyPI
  publish, and deploys (100号機 / mdx / hibino).

Claude does **NOT**: `git push` (for release or otherwise), monitor/poll GitHub
Actions CI, wait for CI-green, run `tools/ci_preflight.py` or
`tools/check_ci.py` watch-loops, push tags, invoke the `release-qud` flow, or
publish to PyPI. If asked to "release", Claude prepares and commits the work,
then hands off to codex.

**Why**: CI monitoring and release driving (queue watching, tag/publish, remote
deploy) is long-running, polling-heavy work that does not need Claude's
reasoning and burns Claude turns. Keep Claude on implement → test → commit;
codex owns everything downstream.

**Exception**: only if the user EXPLICITLY asks Claude to push / handle CI /
release **in that specific task** does Claude do it. That permission does not
carry over to the next task, and it is **never inherited from another session**
or from a previous turn. The default is **commit-and-stop**.

**Codex-side consequence**: a Claude session hands off a local commit SHA and
branch. Taking it from there — `tools/ci_preflight.py`, push, CI-green, tag,
publish, `release-qud` — is codex's lane, and codex should not wait for Claude
to do any of it.

---

**Last Updated**: 2026-04-02
**For**: Codex AI Assistant
**Project**: Radia Magnetic Field Computation
