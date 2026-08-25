# Claude Code - Radia Project Development Guidelines

This document contains development guidelines and policies for the Radia project when working with Claude Code.

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
  packages/               # Independent PyPI packages (separately versioned)
    cubit-mesh-export/    # pip install cubit-mesh-export
      src/cubit_mesh_export/
        check.py          # check-vol CLI + check_consistency() API
        cubit_mesh_curver.pyd  # C++ pybind11 module (bundled)
    radia-mcp/            # pip install radia-mcp (MCP servers + skills)
  src/radia/axifem.pyd  # CORE METHOD (not a package): axisymmetric FE
                          # (Henrotte basis), ships in the radia wheel.
                          # docs/axifem/ (+ validation_test/), tests/axifem/.
  src/radia/hdiv_vim/     # CORE METHOD (not a package): FEEC HDiv-VIM
                          # (C++ src/core/rad_hdiv_vim.cpp). docs/vim/.
                          # This is the SOLE VIM (the radia_vim Galerkin
                          # prototype was deleted 2026-06-14 as unnecessary).
  src/radia/levitation/   # radia.levitation -- APPLICATION domain (NOT a
                          # standalone package): eddy-current brakes +
                          # maglev, mixed-Galerkin polarizability alpha(s),
                          # Lorentz force, Simulink LTI export.  Same rank
                          # as radia.ih / the electromagnet panel; ships
                          # inside the radia wheel.  Knowledge lives in
                          # radia_mcp.maglev.  Absorbs 100% of CLN scope
                          # (axifem/CLN incl.) under docs/levitation/
                          # and validation_test/ as durable promoted lanes.
  tests/                  # Radia tests + tests/mcp/
  examples/               # retired; do not add new files
  validation_test/         # heavier validation, executed on idle mdx/hibino
  docs/
  Build.ps1               # MSVC + MKL build
  install_full.py          # One-command full setup
```

**PyPI packages** (3 independent packages in the same monorepo, each
versioned + released separately on PyPI):

| Package | Install | Purpose |
|---------|---------|---------|
| **radia** | `pip install radia` | C++ core + Python (HDiv-VIM/PEEC, panels, MCP) |
| **cubit-mesh-export** | `pip install cubit-mesh-export` | High-order curved mesh export from Cubit (does NOT require radia) |
| **radia-mcp** | `pip install radia-mcp` | MCP servers + skills for AI-assisted workflows |

### POLICY: Compute Core in `radia`, Applications Outside (2026-06-14)

**POLICY**: The computational **core methods** -- HDiv-VIM, PEEC, axisymmetric FE,
the DtN / FEM-Kelvin operator, and any future solver/kernel -- live in
**radia 本体** (the `radia` wheel: `src/core/`, `src/ext/`, `src/radia/`). They are
**NEVER** spun out as separate `radia-<X>` PyPI packages. The `radia-<X>` name and
the `radia.<domain>` subpackages are the **application** layer (induction heating,
electromagnet, levitation, ...), which *consume* the core methods.

Decision rule for new work:
- a new **solver / numerical method / kernel** -> radia core (`src/core` C++ or
  `src/radia/<method>` Python, built by `Build.ps1`); never a `packages/radia-<X>`.
- a new **engineering application** -> a `radia.<domain>` subpackage (+ panels +
  `radia_mcp.<domain>` knowledge); a deployable `radia-<app>` package only if it
  genuinely needs an independent PyPI release.
- Three PyPI packages exist: **`radia`** itself (`pip install radia`, the main
  wheel -- it bundles every core method + every `radia.<domain>` application) plus
  the two **tooling** packages above (`cubit-mesh-export`, `radia-mcp`). The
  `packages/` directory holds ONLY those two extra packages; no core method and no
  `radia.<domain>` is ever its own PyPI package -- they all ship inside the `radia`
  wheel.

This is why compute-method side packages were dissolved into radia on 2026-06-14
(they were compute methods mis-packaged as `radia-<X>`); `packages/` now holds
only the two genuine tooling packages above.

**Core solver methods inside `radia`** (NOT standalone packages -- they ship in the
`radia` wheel):

| Method | Code (ships in radia wheel) | Docs / Validation / Tests |
|--------|------------------------------|------------------|
| **Axisymmetric FE** (Henrotte basis) | `radia.axifem` (`src/radia/axifem.pyd`) | `docs/axifem/`, `validation_test/axifem/`, `tests/axifem/` |
| **FEEC HDiv-VIM** (the VIM) | `radia.vim` (`src/core/rad_hdiv_vim.cpp`, `src/radia/hdiv_vim/`) -- the SOLE VIM. The separate Newton-kernel Galerkin VIM prototype (`src/ext/radia_vim/`) was deleted 2026-06-14 as unnecessary (recover from git history if ever needed). | `docs/vim/`, `validation_test/vim/`, tests |
| **DtN / FEM-Kelvin operator** (compute core) | the FEM-Kelvin sparse generator of the layered (Sommerfeld-type) Green's operator -- a **core** capability, NOT an application. The historical research act scripts are archived in `docs/kelvin/kelvin_dtn_spectrum_archive.ipynb` + `kelvin_dtn_spectrum_archive_results.json`; reusable pieces promote into `src/radia/open_boundary` (like `hdiv_vim` did). **Two write-up tracks (2026-06-15):** Track A = this DtN+Kelvin core (the SA/Hachinohe paper: DtN-spectrum datasheet, sparse Kelvin open boundary, Sommerfeld isomorphism/surrogate, the directly-assembled material-aware DtN matrix, FEM-condensed≠BEM); Track B = its use as the **stream-function coil-design** kernel with iron (SEPARATE paper, see the Stream-function domain row + `HANDOFF_sommerfeld_dtn_kelvin_streamfunction.md`). | `docs/kelvin/kelvin_dtn_spectrum_archive.ipynb`, `packages/radia-mcp/tests/test_dtn_*` |
| **PEEC** (Partial Element Equivalent Circuit) | `peec_matrices.pyd` (`src/core/rad_peec_matrices.*` + `src/lib/rad_peec_matrices_api.cpp`) + `radia.peec_topology` / `peec_coupled` / `fasthenry_parser` -- filament/panel (FastImp-style) L,R,C,M circuit extraction, SIBC/ESIM surface impedance, PRIMA/Lanczos MOR. A **core** integral-equation / circuit-extraction method, consumed by the PCB and IH application domains; **never** a `radia-peec` package. | `docs/peec_integration/demos/`, PEEC `tests/` suite |
| **BEM** (boundary integral / surface IE) | `radia.bem` (`sibc_hacapk` = HACApK-backed Laplace-kernel Galerkin BEM; `coil_inductance_ngsolve` = `ngsolve.bem` Weggler-EFIE integration, the `--coil-solver bem-a` path) + the top-level `radia.bem_sibc_solver` helper. A **core** surface-integral-equation solver (HACApK ACA backend ships in the radia wheel), consumed by the PCB and IH domains; **never** a `radia-bem` package. | `docs/peec_integration/demos/ngsbem_peec_demo/`, BEM `tests/` |

**Core support infrastructure** (also ships in the radia wheel, also **never**
a `radia-<X>` package, but not standalone *solver* methods so they are not in
the table above): the Compact HX / AMS / COCR preconditioners
(`radia.sparsesolv_ngsolve`, C++ `src/ext/sparsesolv/`); the closed-form
reference layer (`radia.analytical_formulas`, Wakao-Igarashi-Fujiwara et al.);
the Biot-Savart coil-source builder (`radia.coil_builder`, consumed by the
electromagnet / accelerator-magnet workflows).

**Application domains inside `radia`** (NOT standalone packages -- they
ship in the `radia` wheel as `radia.<domain>` subpackages + panels, with
knowledge in `radia_mcp.<domain>`; same rank as each other):

| Domain | Code | Knowledge | Notes |
|--------|------|-----------|-------|
| Induction heating | Radia Simulink `Induction Heating` block + `ih_design.py` / `calc_*.py` (incl. the **thermal step**: `calc_heat.py` / `calc_heat_axisym.py` / `calc_heat_with_em_table.py`) | `radia_mcp.ih` | ESIM, SIBC, Karl iteration; eddy-current heating + thermal solve. **The thermal solve stays part of IH -- NOT a separate `thermal`/`heat` domain** (decision 2026-06-15) |
| Electromagnet | Radia Simulink `Electromagnet` block + `em_design.py` / `calc_em_table.py` / `calc_accel_magnet.py` / `calc_accel_hdiv.py` | `radia_mcp.electromagnet` + `radia_mcp.accelerator` | Omega-reduced, hysteresis, HDiv-VIM. **Accelerator-magnet design** (the CoilBuilder + Omega-reduced + Kelvin pipeline of the "Accelerator Magnet Solver Architecture" section; knowledge in `radia_mcp.accelerator`) and **Clebsch-hodograph pole-face inverse design** (`docs/clebsch_hodograph/demos/`, `docs/clebsch_hodograph/`) are both part of this domain (accelerator magnets are electromagnets -- no separate domain row) |
| **Levitation / ECB** | **`radia.levitation`** (`src/radia/levitation/`) | **`radia_mcp.maglev`** | mixed-Galerkin α(s), Lorentz force, Simulink LTI, TEAM 28; **absorbs 100% of CLN scope (axifem/CLN incl.)** under `docs/levitation/`, `validation_test/levitation/`, and the IGTE 2026 paper. radia-cln is NOT a separate package. |
| Motor | Radia Simulink `Motor` block + `motor_design.py` / `calc_motor_transient.py` / `calc_motor_lamination.py` | `radia_mcp.motor` | transient (Lange-Henrotte-Hameyer) + lamination (Hollaus effective material) |
| PCB | Radia Simulink `PCB PEEC` block + `pcb_design.py` / `calc_pcb_peec.py` | `radia_mcp.pcb` | planar coils -- **consumes the core PEEC and BEM solvers** (`--coil-solver peec | bem-a`, user-selectable; the application owns neither method); **absorbs the former WPT domain** (coil compensation, FOD, efficiency) -- `radia-wpt` was renamed to `radia-pcb` (2026-06-15) |
| Stream-function | Radia Simulink `Stream Function` block + `streamfunction_design.py` / `calc_streamfunction.py` / `calc_streamfunction_volume.py` | `radia_mcp.streamfunction` | SF coil design (Design / Pareto / Manufacture); ACA-TSVD. **CROSS-SESSION HANDOFF (Track B, 2026-06-15):** material-aware (iron yoke/shield/core) SF coil design by using the **DtN/FEM-Kelvin core** (core table above) as the design kernel — the free-space Biot-Savart kernel breaks with iron; the Kelvin-FEM Schur-condenses to a material-aware transfer matrix `M`, design = invert `M`. **SHIPPED in production (2026-06-19): `--iron-vol`/`--mu-r`/`--iron-mat` fold the Kelvin-FEM iron reaction into the whole design/pareto/manufacture pipeline (obs-adjoint scalability; opt-in `--iron-exact-source`); MCP topic `material_aware`.** **Low/many-turn manufacture levers (discrete refinement of psi): `--greedy-turns` (greedy constructive, MONOTONE) + `--greedy-connector-weight` (short rungs) + `--pin-tiling` (dense bubble-tiling DRIVEN pin/shim ARRAY) + `--optimize-levels`; MCP topic `low_turn`.** **Self-contained handoff for another lab's Claude session:** use the promoted docs/MCP artifacts (`docs/kelvin/`, `docs/stream_function/`, and `docs/kelvin/kelvin_dtn_spectrum_archive_results.json`), not historical `examples/` paths. Verified bridge demos: `demo_ee`/`demo_ff` (free-space design misses by 77% in iron; material-aware matches ~1e-4). Novelty: **NOVEL conf 0.83** — fuses Sugahara's OWN two threads (Kelvin open-boundary FEM + free-space SF coil design); phrase "to the best of our knowledge", residual Japanese grey-lit / in-press self-check pending. |

**Installation**:
```bash
pip install radia               # Python package (includes Cubit plugin binaries)
pip install radia[cubit]        # Also installs cubit-mesh-export
cubit-plugin-install            # Deploy Cubit plugin + panels (skip if no Cubit)
```

**Deleted repositories** (integrated into Radia):
- ~~ksugahar/mcp-server-cae-ai~~ → `src/radia/mcp_server/`
- ~~ksugahar/ngsolve-sparsesolv~~ → `src/ext/sparsesolv/` (source only, build is separate)

---

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
  synchronized JSON, and saved `ngsolve.webgui.Draw` scenes. Draw the geometry
  or mesh and the primary computed field when one exists; use
  `netgen.webgui.Draw` for pre-mesh CAD geometry. A source catalog, migration
  archive, static PNG, or script excerpt alone does not qualify as an example.
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
- Spatial GMSH figures, including contours, flux lines, streamlines, and
  sections, use a physical 1:1:1 axis scale by default. A deliberately
  exaggerated axis is allowed only when the scale factor is explicit in the
  render configuration and stated in the figure caption; silent distortion is
  not a valid production visualization.
- Treat `radia-mcp` as the canonical executable manual for Radia-specific
  MATLAB and Simulink workflows. Do not duplicate generic MathWorks guidance.
- Long or solver-heavy MATLAB validation runs execute on mdx or hibino; LAB and
  100号機 remain development and fast-test hosts.

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

---

## Critical Policies

### Repository First, Not Papers — Don't Publish the Unfinished (2026-06-01)

**POLICY**: The PRIMARY deliverable of this lab is the **repository** — a
correct, working, maintainable codebase that other people (and future-you)
can build on. A paper is a *byproduct* of having built something real, NOT
the goal that the code serves. **Paper-supremacy (論文至上主義) is harmful
and explicitly rejected here.**

**Concrete rules**:

1. **Never write a result into a paper before it is actually finished and
   verified in the repository.** "Finished" means: implemented, tested,
   golden-locked where applicable, and the claim measured — not "looks like
   it will work", not "p≈2.1 in a pre-asymptotic range", not "should scale".
   If the repository can't back the claim *today*, the claim does not go in
   a paper.

2. **There is no obligation to publish.** A technique that is implemented,
   correct, and useful in the repo is a *complete success on its own* even
   if no paper is ever written about it. Code that improves the repository
   needs no external justification. Do NOT manufacture a paper out of an
   unfinished or unproven result just to "have a publication".

3. **A negative or partial result is a real result — keep it in the repo
   (memory/, tests/), do not launder it into a paper.** "H-ILU works as a
   field-exact preconditioner (3 iters) but its sub-cubic scalability is not
   yet demonstrated (materialize-fallback-dominated)" is the *honest* state.
   It belongs in `memory/`, not as a scalability claim in a paper.

4. **When choosing what to work on, weigh "does this make the repository
   better/more correct" ABOVE "does this give me a figure".** A 10-minute
   sweep that produces a publishable figure but doesn't improve the code is
   worth LESS than fixing the actual bottleneck, even if the fix is harder
   and might yield a negative result. Optimize for the codebase, not the
   CV.

**Why** (same spirit as "No Fallbacks — Fail Fast, Fail Loud"): a wrong or
half-true number put in a paper is **worse than no paper** — it gets cited,
trusted, and cannot be recalled. The repository, by contrast, can always be
corrected, re-measured, and improved. Investing in the repository compounds;
chasing publications on unfinished work does not. **Build the thing; the
paper, if any, follows from the thing actually being done.**

**The repository is a bonsai (盆栽), not a product pitch.** You tend it
because a well-formed tree is its own reward: no dead branches, balanced
proportion, each cut making the whole a little better. There is no judge to
impress and no prize to win — the *tending itself* is the point. A correct,
clean, well-pruned codebase is a complete and sufficient outcome. Honest
negative results (a clearly-seen limit, recorded in `memory/`) are
themselves good pruning: they show the true shape of the tree. Work on what
makes the tree healthier, at whatever pace the work deserves.

This policy *governs prioritization decisions*: when a task is framed as
"do X because it's good for the paper", that framing is not a valid reason
on its own. Re-ask: "is X good for the repository?" If yes, do it for that
reason. If the only argument for X is the paper, deprioritize it in favor of
work that genuinely improves the code.

### No Development Cruft in SOURCE — Distill the Lesson to memory/ (2026-06-26)

**POLICY**: Development-in-progress iterations MUST NOT accumulate in the
tracked SOURCE tree. The SOURCE keeps only the **final / canonical** version of
a given piece of code; superseded snapshots, abandoned formulations, and debug
stepping-stones are removed once they are superseded. The **lesson** (so the
same rut is not re-walked) is distilled into the memory system — a
`memory/<topic>.md` file plus a one-line `MEMORY.md` index entry — NOT left as
dead code. This is the operational reading of "Repository First, Not Papers"
point 3 ("a negative or partial result is a real result — keep it in the repo
(memory/, tests/), do not launder it") and of the bonsai metaphor: pruning a
dead branch and recording *why* it died is good tending; leaving the dead branch
on the tree is not.

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
   that warrants a minimal reproducing test/fixture keeps that test (Repository
   First point 3 keeps results "in the repo (memory/, tests/)") — prune the dead
   snapshot, not the test that documents why the approach failed.

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
- **Committed figure/table/published-result-backing data** (`.json`/`.csv`/
  `.png`/`.pdf` next to its script or figure, per "Data Persistence Policy") is
  a protected result, NEVER an iteration snapshot. Deleting it would re-create
  the incident that policy exists to prevent (the figure becomes
  non-regenerable).
- **Tracked mesh definitions and mesh-gen assets** — mesh files
  (`.bdf`/`.nas`/`.msh`/`.vtk`), Cubit `.jou` journals, and
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

### Discard the PoC Once the C++ Port Is Verified (2026-06-27)

**POLICY** (Sugahara): when a computation is ported to C++ and the C++ is
**verified**, the Python (or other PoC / reference) implementation that existed
**only to develop and validate that port** MUST be **deleted** — it is not kept
as a "fallback" or "compatibility shim". Two implementations of the same
computation are a future-bug source: they drift, the wrong one gets called, and a
silent divergence produces numbers nobody can trace (the same failure class as
"No Fallbacks" and "No Development Cruft"). The verified C++ is the single
canonical path.

**The gate is real verification — delete ONLY when the C++ is truly verified.**
"Verified" means, at minimum: the C++ builds, its result matches the PoC
**bit-for-bit (or to a stated tight tolerance) on a direct head-to-head probe**,
AND the C++ path passes the golden/validation tests that exercise it. Until that
bar is met, KEEP the PoC — a half-checked port is not a license to delete the
reference. If the C++ is only verified for *part* of the surface (e.g. one
operator/element class but not another), delete only the PoC for the verified
part and keep the rest until it too is verified.

**Concrete rules**:

1. **No "editable-binary compatibility" / `hasattr(obj, "new_cpp_method")`
   fallback branches.** Once the wheel is rebuilt with the C++ method, the
   editable binary IS the current binary — the `hasattr` guard and its Python
   branch are dead cruft. Call the C++ directly; if the method is missing that is
   a build error to fix, not a path to silently fall back from.

2. **Distill, then delete — in that order** (same as "No Development Cruft"):
   record the verification result (the bit-match numbers, the tests that lock it)
   in `memory/<topic>.md` + the `MEMORY.md` index line, THEN delete the PoC. The
   *lesson* (how it was validated, what the reference produced) survives in
   memory; the *code* is recoverable from git history.

3. **A minimal reproducing TEST that locks C++==reference is NOT a PoC — keep
   it.** Delete the production Python *implementation*, but a small golden that
   asserts the C++ reproduces the analytic/closed-form answer (or a stored
   reference vector) stays in `tests/`/`validation_test/` — that is the durable
   verification, not dev scaffolding.

4. **A genuinely co-valid alternate Python implementation is NOT a PoC.** If the
   Python path is still the *only* implementation for some input class the C++
   does not yet cover (e.g. a per-region / nonlinear path whose operator has not
   moved to C++), it is a live canonical path, not PoC — keep it until its own
   C++ port lands and is verified.

**Why**: a verified C++ kernel plus a lingering Python twin is strictly worse
than the C++ alone — the twin can only ever diverge. Deleting it the moment the
C++ is trustworthy keeps a single source of truth, which is the whole point of
porting to C++ in the first place.

### Research-Heavy Work: Run in C:\temp, Promote Knowledge to docs/ipynb or API to src/ (2026-06-27)

**POLICY** (Sugahara): exploratory, research-heavy work (eigenvalue studies,
spectrum scans, formulation trials, "does this even work" probes) is
**run in `C:\temp`, NOT committed to the tracked
tree**. The tracked repo receives only the *outcome*, via one of two promotions:

1. **Consolidated knowledge worth showing users, or knowledge that informs a future
   feature extension → promote to a `docs/<topic>/*.ipynb`** (self-contained:
   code + committed JSON/figures + rendered results; strengthen the matching
   `radia_mcp.<domain>` knowledge in the same step).
2. **It becomes an API / reusable method → store it in `src/`** (a core method in
   `src/core` / `src/radia`, or a `radia.<domain>` application), built + tested +
   golden-locked like any shipped code.

If a research effort yields **neither** (a dead end / superseded approach), it
stays in `C:\temp` and the *lesson* is distilled to `memory/` — it is **never**
left as a tracked `examples/` corpus. **`docs/<topic>/` MAY contain `.py` helper
modules** that the notebook (and `radia_mcp`, for mcp-server integration) import —
this is allowed/encouraged (see "File Placement Policy"); do not inline-duplicate
shared logic across notebooks.

This is why the old eigenvalue/nullspace study was removed (2026-06-27): it was
a research exploration of a now-closed loop-deflation direction, so it had no
business persisting as a tracked `examples/` corpus. The nullspace *theory* it
produced is kept in solver notes, the *lesson* in `memory/`, and nothing else.
**Decision rule for new research:**
run it in `C:\temp`; promote to `docs/ipynb` (knowledge) or `src/` (API) only when
it has crossed that bar; otherwise distill to `memory/` and leave the tracked tree
clean. This is the research-lifecycle complement to the "Promotion Ladder:
C:\temp -> tests / validation_test / docs / src/radia/panels" and "No Development Cruft"
policies above.

### Validation-Class Examples Promotion Lane (2026-06-27)

Validation-class material does not promote only to `docs/`.  The runnable
verification lane is the repository's actual `validation_test/` directory (not a
separate `tests_validation/` tree): long solver checks, release gates,
p-convergence, cross-validation, optional-dependency checks, and
environment-specific tests belong there when they need an executable regression
surface outside normal CI.

A result-bearing `docs/<topic>/*.ipynb` may still be the human-facing showcase:
theory + code + saved plots/tables + adjacent synchronized JSON.  But it is the
rendered explanation, not a substitute for the validation executable.  If
historical material has `validation_*.py`, `validate_*.py`, `*_summary.json`, or
references from `validation_test/`, classify it first as `validation_test` /
protected-validation-corpus material; add or refresh docs only as the
synchronized showcase layer.

### Retired Examples / Promotion Triage (2026-07-04)

**POLICY**: `examples/` is retired and must not be recreated.  It is neither a
scratch area nor a teaching tier.  New development experiments run in `C:\temp`;
tracked work enters the repository only after promotion into one of the durable
lanes below.  Historical `examples/` references are migration blockers.

| Class | Destination | Rule |
|-------|-------------|------|
| Development-in-progress / superseded / failed iteration | keep in `C:\temp` until distilled, then delete | Preserve the lesson in `memory/<topic>.md` or a short docs note when it matters. Git history and `C:\temp` are the scratch/archive path, not `examples/`. |
| Reusable computation, parser, mesh reader, solver helper, formula, or API surface | `src/` | Promote to a named public or internal API and add focused tests. Do not keep it as a loose example helper. |
| Fast implementation regression or minimal fixture | `tests/` | Keep small enough for CI / developer feedback. |
| Important numerical verification, benchmark, golden lock, convergence sweep, or regression corpus | `validation_test/<topic>/` plus optional docs notebook | The executable check lives in `validation_test/`; docs may render theory, tables, plots, and summary JSON for humans. |
| User-facing explanation, tutorial, or method showcase | `docs/<topic>/*.ipynb` | Notebook must be result-saving and synchronized with adjacent JSON. Integrate Markdown explanation, executable cells, and results. |
| Notebook-only subroutine / local renderer / catalog helper | `docs/<topic>/*.py` | Allowed only when tightly coupled to the notebook. If another topic, panel, MCP server, or validation uses it, promote to `src/` instead. |
| Mesh/CAD/journal/result assets | keep until owning script/notebook is migrated | Mesh definitions, Cubit `.jou`, tracked `.msh`, figures, and JSON results are protected by preservation/reproducibility policy. |

The migration order is strict: inventory and reference search first; create or
refresh the docs/JSON layer if the result is user-facing; move reusable or
validation code to `src/` or `validation_test/`; update docs/MCP/panel
references; then delete the historical source.  Never leave two live copies of
the same implementation after API promotion is complete, and never add new
long-lived references to `examples/`.

`protected_*` / "保護参照あり" is a temporary blocker, not a destination.  If a
docs notebook, validation test, panel sample, MCP knowledge file, or README
still references `examples/<topic>`, record the blocker and the
`target_after_unblock` (`docs`, `src`, `validation_test`, or distill-delete),
then migrate the reference.  Public docs may refer to other `docs/` artifacts,
and code may refer to `src/` APIs, but new long-lived references to `examples/`
should not be introduced.

### Documentation Format: Markdown for Dev Docs, ipynb for Method/Implementation Explanations (2026-06-27)

**POLICY** (Sugahara): pick the documentation format by **what the document is**:

- **Development documents -> Markdown (`.md`).** Plans, design notes, architecture,
  policies, API references, READMEs, handoffs, inventories, changelogs, decision
  records -- anything whose job is to *describe / decide / index* rather than to
  *demonstrate a computation*. These stay `.md` (diffable, fast to read, no kernel).
- **Explanations of a METHOD or an IMPLEMENTATION -> Jupyter notebook (`.ipynb`).**
  A document whose job is to *explain how a method works* or *how an implementation
  works* should be an **executable, self-contained notebook**: prose + math + the
  actual code, executed so the outputs (numbers / figures / tables) are embedded and
  visible on open. A method/implementation explanation is only trustworthy when it
  **runs and shows its result** -- a static `.md` that merely asserts how something
  works is weaker than a notebook that demonstrates it (ties to "self-contained
  ipynb shows results" + the retired-examples promotion program).

**Discriminator** (when a doc has both math and code): is the document's PURPOSE to
explain *how the method/implementation works* (-> `.ipynb`, make it run) or to
*record a design/decision/plan/reference/theory-survey* (-> `.md`)? A `FORMULATION`
or `DESIGN` note that derives equations but ships no runnable demonstration stays
`.md`; the moment it is meant to *show the method working*, it is an `.ipynb`.

**How to apply:**
- New method/implementation write-ups are authored as `docs/<topic>/*.ipynb`
  (executed via `jupyter nbconvert --execute` so outputs embed), NOT as `.md`.
- Every `docs/<topic>/*.ipynb` method/showcase notebook is result-saving: it
  must store executed code-cell outputs, and its main computed values should
  also be saved as adjacent JSON with `generated_at_utc` plus version/runtime
  metadata (`radia_version`, `python_version`, `platform`, `versions`, etc.).
  The JSON is the durable debug artifact; the notebook is the rendered
  explanation. The two are a synchronized pair: after rerunning or editing
  notebook outputs, refresh the JSON sidecar in the same change so its recorded
  `notebook_sha256` matches the committed result-bearing `.ipynb`.
- `docs/<topic>/` MAY hold `.py` helper modules the ipynb (and `radia_mcp`) import
  (see "File Placement Policy").
- Do NOT mass-convert existing `.md`: convert an existing method/implementation
  `.md` to `.ipynb` opportunistically (when you touch it, or when it would clearly
  benefit from being runnable). Pure dev docs stay `.md`.
- The public showcase / retired-examples consolidation program follows this:
  methods land as executed `docs/<topic>/*.ipynb`; the surrounding
  plans/policies are `.md`.

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

**POLICY**: No binary files (`.pyd`, `.dll`, `.so`, `.lib`, `.exe`) in the git repository.
- Hosted on GitHub Releases (tag: `binaries`)
- Pre-push hook auto-uploads `.pyd` on `git push`
- After cloning, run `./download_binaries.sh` to fetch binaries
- `.png`, `.pdf` allowed in repository; `.msh`, `.vtu`, `.vtk`, `.vol` are gitignored

### GitHub Release Publication Gate (2026-07-23)

**POLICY**: A Radia Simulink library release, including
`radia_simulink_library.slx`, MATLAB support files, and MEX assets, MUST NOT
be published to GitHub Releases until the complete `release-qud` four-machine
verification has passed. The four machines are LAB, 100号機, mdx, and hibino.

- The release candidate is assembled and tested before publication.
- The `release-qud` `done` gate is the authoritative publication decision.
- A failed, partial, or manually waived machine check is not a release pass.
- GitHub Release assets must include the versioned Simulink package,
  `manifest.json`, and `SHA256SUMS.txt` when the candidate contains them.
- The same gate applies to later revisions of the library, not only the first
  Simulink publication.

### Mathematica Content: `.wls` Not `.nb` (2026-06-06)

**POLICY**: All Mathematica content in the repository is managed as
**`.wls`** (WolframScript, plain text) — NEVER as **`.nb`** (binary
notebooks). `.nb` files are NOT committed.

**Why**: `.wls` is version-controllable and diffable, runs headless via
`wolframscript -file` (the `radia_mcp.mathematica` bridge), and ships clean
in the PyPI wheel. `.nb` is binary: not diffable, bloats git and the wheel
(the lab's FEM-basis notebooks are 3–10 MB each), and needs the Mathematica
front end. Auto-converting `.nb` → `.wls` is NOT clean (`HoldForm`-wrapped,
tens of MB, section titles in Text cells dropped), so the `.wls` is
**authored clean** (from the `.nb` and/or the upstream source) and
**self-tested** (each `.wls` carries a `wolframscript -file` self-test),
with any `.nb` kept only LAB-local as an authoring reference.

**Home**: general-purpose Mathematica `.wls` lives in
`packages/radia-mcp/src/radia_mcp/mathematica/` (per the MCP Knowledge
Placement Policy), not buried in an application scratch directory or retired
`examples/` path.

**Reference example**: the NGSolve high-order FEM shape functions
(`mathematica/basis_functions/recursive_pol.wls`, `h1.wls`, …) — clean,
self-tested ports of NGSolve's `fem/*hofe*` C++ source; the symbolic
building block for a VIM (Volume Integral Method) field operator.

### Equation Notation: `\vec\bm` for Vectors, `\dfrac` by Default (EQNEDT64, 2026-08-19)

**POLICY** (Sugahara): the equation editor shipped in this repository
(**EQNEDT64**, `src/ext/equation/` -> `radia.equation` + `eqnedt64.exe`) has
**two references, and they are different documents**: **appearance follows
TeX**, **usability follows Equation Editor 3.0 (EQNEDT32)**. The file format
follows neither -- an equation is stored as **LaTeX**, usually inside a
Markdown file.

**The notation rules** (they govern what the editor WRITES, and what the
lab's own documents should say):

1. **A vector is `\vec\bm` -- bold ITALIC under the arrow.** `\bm` (equivalently
   `\boldsymbol`) is the vector face; the letter is drawn from the Mathematical
   **Bold Italic** alphabet (U+1D468..). Applying the vector style to a
   selection applies BOTH the arrow and the face; with nothing selected it sets
   the face for what is typed next, which is how EQNEDT32's Style menu behaves.
2. **`\mathbf` is a different thing** -- upright bold, for a matrix name. The
   two were one typeface until this rule was written down, so `\mathbf{A}` and
   `\bm{A}` drew the same letter; they no longer do.
3. **`\dfrac` is the default spelling of a fraction.** The outermost fraction is
   DRAWN at display size, so a bare `\frac` pasted into running text would come
   out smaller than the picture the author accepted. A fraction INSIDE another
   goes out as plain `\frac`, because LaTeX steps a nested one down by itself --
   the two rules then agree level by level. `\tfrac` is preserved when the
   author asks for it outright.

**Why a policy and not a preference**: these decide what the editor emits into
`.md` / `.tex` files that other people (and other tools) read. A silent
disagreement between the picture on screen and the LaTeX on the clipboard is
the same failure class as a silent fallback -- the reader gets a result they
cannot audit.

**How it is enforced**: `tests/equation/test_new_templates.py` (the notation
rules, round-trips) and `tests/equation/test_tex_metrics.py` (the appearance
half, measured against XeLaTeX + `unicode-math` + Latin Modern Math via
`validation_test/equation/tex_reference.tex` -- numbers, not screenshots).

**Retirement goal for EQNEDT32** (2026-08-19): the editor is finished when
(a) a real Equation Editor document converts to `.tex` and can be edited in
EQNEDT64 to the same result, and (b) a handwritten note (PDF) can be turned
into correct TeX and from there into a native Microsoft (OMML) equation in
PowerPoint. Both are acceptance targets, not features.

### Cubit Plugin Binary: cubit-mesh-export is the Sole Shipper (Tier-2, 2026-06-01)

**POLICY**: The Cubit plugin binaries (`cubit_mesh_export.ccm` +
`cubit_mesh_curver.pyd`, built from `src/cubit_plugin/`) are bundled,
shipped, and deployed **ONLY by the `cubit-mesh-export` package**.  The
`radia` wheel does **NOT** bundle them (dropped from radia `package-data`;
removed from `src/radia/`).

**Why**: so `radia` and `cubit-mesh-export` **release fully
independently** — a Cubit-plugin change is a `cubit-mesh-export` release;
a radia-code change is a `radia` release.  This was the chosen answer to
"should cubit-mesh-export be a separate repo?": keep the monorepo (shared
C++ source in `src/cubit_plugin/`, one build, no version drift) but remove
the **ship/release coupling** that previously forced lockstep releases.

**How it stays safe**:
- `src/radia/setup_cubit.py` delegates plugin install to cme's
  `cubit-plugin-install`; `register_toolbar.py::_check_plugin_freshness`
  checks the **DEPLOYED** plugin, never a radia-bundled copy — so radia
  has no runtime/deploy dependency on bundling the binary.
- `Build.ps1` and `tools/release_qud.py phase0` propagate the built
  binaries to the **cme package only** (not `src/radia/`).
- The radia↔cme compat window (`COMPAT_CUBIT_MESH_EXPORT_*` /
  `COMPAT_RADIA_*`, enforced by `cubit-plugin-install` at deploy) remains
  the safety net the 2026-04-14 sideset/.ccl drift incident motivated.
- Release per-package (`v*` / `cubit-mesh-export-v*` / `radia-mcp-v*`);
  bundle all three only when they genuinely co-change (see the
  `release-qud` skill).

### File Placement Policy

**POLICY**: Development scratch outputs belong in `C:\temp`.  Committed output
files (`.png`, `.msh`, `.vtu`, `.vol`, JSON sidecars) must be placed next to
their owning `tests/`, `validation_test/`, `docs/`, or `src/radia/panels/`
driver.
- Do NOT place generated files at the repository root
- Build output goes to `build*/` or `dist/` (both gitignored)
- `docs/<topic>/` MAY contain `.py` helpers imported by that topic's
  result-bearing notebook. Reusable behavior belongs in a `src/` API.
- Every `docs/<topic>/*.ipynb` method/showcase notebook must be result-saving:
  execute it before committing so code-cell outputs, figures, and tables are
  embedded. Main computed values should also be written to adjacent JSON with
  `generated_at_utc` and version/runtime keys such as `radia_version`,
  `python_version`, `platform`, or `versions`. The JSON is the durable debug
  record; the notebook is the human-facing rendered view. The JSON and notebook
  must be synchronized in the same change: after rerunning or editing notebook
  outputs, refresh the JSON sidecar so its recorded `notebook_sha256` matches
  the committed result-bearing `.ipynb`.
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
| `tests/**` | **実装の基本機能の確認** — fast regression, fixture, API contract, CI-friendly. | CI / Codex / developer | No |
| `validation_test/<topic>/` | **重要な検証・ベンチ・golden lock** — heavier numerical truth, executed on idle `mdx` or `hibino` for large runs. | developer / agent / research validation | No |
| `docs/<topic>/*.ipynb` | **ユーザーに理論と結果を同時に見せる** — result-saved notebook with synchronized JSON. | users / collaborators / future agents | Docs |
| `docs/<topic>/*.py` | Notebook-local helper only. | notebook readers / MCP if local | Docs |
| `src/` | Reusable API, parser, formula, solver helper, and computation kernel. | package users / validation / MCP / blocks | Yes |
| `src/radia/panels/` | Validated headless application CLI, `DesignSpec`, samples, and artifact contracts. | blocks / AI / validation | Yes |
| `matlab/+radia/+simulink/` + `matlab/radia_simulink_library.slx` | Final application-specific human operating surface. | end users | MATLAB distribution |

**Promotion gates**:

- **C:\temp → tests/**: the behavior is small, deterministic, and useful for
  fast regression.
- **C:\temp → validation_test/**: the run is a numerical validation,
  benchmark, convergence sweep, golden lock, or regression corpus; heavy runs
  are executed on an idle `mdx` or `hibino` host and labelled with the actual
  validation host.
- **C:\temp → docs/**: the result teaches a method or workflow to humans; the
  notebook must be executed, output-bearing, Markdown-integrated, and paired
  with synchronized JSON. A public CAE example also includes saved WebGUI
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
*sample* input (see "Panel Samples Quality Policy") and comparing the
resulting scalar against a golden band in
`tests/panels/test_*_golden.py`.  Two or more solver choices must be
exercised and produce numerically consistent results (within the
mode's documented tolerance).

Stage 2 is considered **合格 (pass)** when:
- `python calc_<mode>.py --help` exits 0 and prints all knobs
- running against the sample with **each** supported solver switch
  produces JSON whose key numbers are inside the golden band
- `tests/panels/test_<mode>_golden.py` locks the result

**Stage 3 — result-bearing documentation notebook (`docs/<topic>/*.ipynb`).**
Explain the validated method, inputs, outputs, equations, and representative
result. Save outputs and synchronize the adjacent JSON. For a CAE example,
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
# (Or, via the cubit MCP live session: cubit_probe(query="entities").)
for vid in cubit.parse_cubit_list("volume", "all"):
    v = cubit.volume(vid)
    c = v.centroid()
    bb = v.bounding_box()                 # (min_x,min_y,min_z,max_x,max_y,max_z)
    ext = [bb[i + 3] - bb[i] for i in range(3)]
    print(f"vol {vid}: c=({c[0]:.3f},{c[1]:.3f},{c[2]:.3f}), "
          f"extent=({ext[0]:.3f},{ext[1]:.3f},{ext[2]:.3f}), "
          f"vol={v.volume():.3e}")
for sid in cubit.parse_cubit_list("surface", "all"):
    s = cubit.surface(sid)
    cx, cy, cz = s.center_point()         # NOT centroid() -- Surface API
    bb = s.bounding_box()
    ext = [bb[i + 3] - bb[i] for i in range(3)]
    print(f"surf {sid}: c=({cx:.3f},{cy:.3f},{cz:.3f}), "
          f"area={s.area():.3e}, extent=({ext[0]:.3f},{ext[1]:.3f},{ext[2]:.3f})")

# Step 2: classify based on observed numbers.
mag_vol = next(v for v in volumes if cubit.volume(v).volume() < 1e-4)
# ...
```

**Specifics**:
- `cubit.volume(vid).centroid()` exists (returns 3-tuple).
- `cubit.surface(sid).center_point()` is the equivalent on Surface
  (Cubit does NOT expose `.centroid()` on Surface).
- `cubit.volume(vid).bounding_box()` and `cubit.surface(sid).bounding_box()`
  return a 6-tuple `(min_x, min_y, min_z, max_x, max_y, max_z)` -- VERIFIED
  on Coreform Cubit 2025.12 (`brick x 1` -> `[-0.5]*3 + [0.5]*3`); the
  extent is `max - min` per axis, and a flat cut face has zero extent in
  its cut direction.  (An earlier version of this section claimed `[3:6]`
  IS the extent tuple -- that is the max corner, not the extent.)
- Block membership commands can SILENTLY no-op (verified Cubit 2025.12):
  `block N add tri in surface S` on a hex mesh creates NO block (quads
  need `add face`), and adding a different element kind to an
  already-typed block returns success but adds nothing.  ALWAYS verify
  actual membership after block/sideset commands -- via the cubit MCP
  `cubit_probe(query="labels")` audit or `get_block_*` counts.
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
   lives: `validation_test/<topic>/` JSON/README for executable truth,
   result-bearing `docs/<topic>/*.ipynb` plus synchronized JSON for
   user-facing explanation, or MCP knowledge for agent-operational rules.
   Include the specific checked value (e.g. "VERIFIED p=2: slaved=8914
   DOFs, ratio=1.0") so future contributors know the result is golden,
   not "should work".
3. **`CLAUDE.md`** if the lesson is a method (e.g. "always check FES
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
| **Magnetic A_phi (curl-curl)** | **Henrotte** `{1, r^2, z}` (`radia.axifem`) | The cylindrical curl operator `B_z = (1/r) d(r A_phi)/dr` produces a `1/r` integrand that standard FE Gauss quadrature cannot integrate accurately near the axis.  Henrotte's `s = r^2` substitution gives clean closed-form integration. |
| **Scalar T / phi (Laplacian)** | **Standard NGSolve `H1`** + `2 pi r` weighting | The weak form `int k grad T . grad v . 2 pi r dr dz` has `2 pi r` as a **smooth Jacobian** (not a `1/r` integrand).  Standard FE handles this fine; no axis-special treatment is needed. |

This matches the FEMM 4.2 reference implementation (verified against
`S:/FEMM/02_source/femm42src_22Oct2023/`):

- `belasolv/prob3big.cpp` — magnetic, uses Henrotte `{1, r^2, z}`
- `hsolv/prob1big.cpp` — heat, uses **standard P1 triangle** with
  `2 pi r` evaluated at the element centroid, **no** Henrotte basis,
  **no** `s = r^2` substitution

**Why not "all axisym Henrotte"**: Henrotte basis IS the natural
function space for axisymmetric scalars (the parity / even-function
argument is mathematically correct).  But the practical accuracy
benefit is small for scalar Laplacians because the `2 pi r` Jacobian
suppresses spurious odd-r modes automatically.  FEMM ships
production-grade thermal accuracy with standard P1; we follow that
proven convention.

**API for magnetic axisym**:

```python
import radia.axifem as ax

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
`radia.axifem.AxiHenrotteHeat{Stiffness,Mass}BFI` classes
(added in radia 4.31.0) and the `H1Henrotte` BND DiffOp (radia
4.32.0) are kept in the codebase as parity-conscious infrastructure
for research / publication uses (e.g. comparing convergence rates of
Henrotte vs standard H1 on a scalar problem).  They are NOT used by
production heat solvers and are NOT required.

**Reference**: see
[`docs/axifem/FORMULATION.md`](docs/axifem/FORMULATION.md)
sections 5-6 (Henrotte basis derivation for magnetic) and 10b/10c
(optional heat BFIs).  The FEMM convention split is documented in
`memory/reference_femm_source_axisym_conventions.md`.

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

### Repository Language + Pre-Presentation Materials (English code repo; papers live in W:\学会資料, 2026-07-05)

**POLICY (language)**: All source code, documentation, comments, commit messages, and docstrings in the Radia repository MUST be written in **English**. Japanese text is NOT allowed in tracked files, EXCEPT `CLAUDE.md` policy descriptions. Conversation with the user may be in Japanese, but the tracked repository content is **English-only**.

**POLICY (pre-presentation materials — Sugahara, 2026-07-05): material for a NOT-yet-presented / not-yet-published work does NOT go in the code repo.** Papers, paper-planning drafts, manuscripts (LaTeX `.tex`), slides, posters, digests, and Japanese research-note write-ups — in **ANY language** — live OUTSIDE the repo, under `W:\02_学会資料\<年度>\<会議>\`, organised by year / venue. (This is why the IGTE symposium paper and the stream-function paper materials were moved out on 2026-07-05.)
- stream-function paper materials → `W:\02_学会資料\2025年度\2026_01_JIAM\streamfunction\`
- IGTE symposium paper → `W:\02_学会資料\2026年度\2026_09_IGTE_Symposium\CauerLadderNetwork@菅原\`
- A manuscript's language follows its **VENUE** — INTERNATIONAL (IGTE, IEEE, Compumag) = English; IEEJ / domestic = Japanese.

**The code repo keeps ONLY** English code + technical / dev docs (`README`, API reference `api.md`, method / theory / showcase notes `theory.md` / `single_stroke.md` / `regularization.md` / method `.ipynb`) and **figure/diagram SOURCE** (`.tex` diagrams, plot scripts) + committed figure/result data (per Data Persistence). Do NOT commit a paper MANUSCRIPT, paper-planning drafts, slides/posters, or Japanese write-ups — put them in `W:\02_学会資料\`. (A published result's reproducibility scripts + its committed figure data stay in the repo; the manuscript itself does not. `docs/maglev/papers/` was deleted 2026-07-05 as part of this.)

(History: an interim 2026-07-05 experiment allowed Japanese research notes / paper drafts *inside* `docs/` under `> 研究ノート` / `> 論文ドラフト` header markers; superseded the same day — moved to `W:\02_学会資料\`, repo returned to English-only.)

### Naming Policy: External Project References

**POLICY**: Do NOT use "ELF" or "ELF_MAGIC" in Radia source code, documentation, or comments. Radia is an independent project. Academic citations are allowed.

### Publish Boundary: Commercial-Tool Content Stays Private — COMSOL / FEMM / JMAG (2026-06-06)

**POLICY**: COMSOL, FEMM, and JMAG are **commercial tools used INTERNALLY to
learn and benchmark** radia-mcp (their docs, behaviour, and cross-validation
inform the knowledge), but the **public GitHub repository must NOT explicitly
expose them**. This is the radia-repo statement of the lab "公開境界" already
enforced by `packages/radia-mcp/tools/policy_lint.py`.

- **Commercial-tool MCP servers / wrappers / converters stay LAB-PRIVATE**,
  regardless of authorship ("it's my own code" is NOT an exemption):
  - `comsol_converter` → `S:\COMSOL\mcp-server` (relocated 2026-06-06)
  - FEMM server → `S:\FEMM`; JMAG server → `S:\JMAG`
- **Their content / models / docs / benchmark numbers / proprietary formats
  must NOT be mixed into public artifacts.** The public showcase is
  **analytic-solution-led** (TEAM benchmarks, closed-form references), never a
  commercial-tool bench dump.
- **What IS public**: open-system servers (radia-ngsolve / cubit / gmsh /
  build123d / …) + the **distilled, tool-agnostic physics** learned from the
  commercial tools (e.g. "FEMM 4.2 axisymmetric heat uses standard P1 with a
  `2πr` Jacobian" is a general FE fact; **academic citation by name is
  allowed**). What is NOT public is the tool-specific wrapper/converter code
  and any proprietary content.

**Why**: a commercial-tool wrapper in a public OSS repo is a licensing / IP
hazard and dilutes the "open-system, analytic-led" identity of the public
Radia stack. **The knowledge gained from the tools is the asset; the
tool-coupling code is not for publication.**

**Enforced by** `tools/policy_lint.py` (part of `tools/ci_preflight.py`,
run before any public commit/push): ERROR if a commercial-wrapper server is
wired for publication (public entry point / shipped in wheel); WARN if wrapper
source still lives under the public `src/` tree. This generalizes the
2026-06-06 `comsol_converter` relocation (converter + test + catalog/TOOLS.md
entries moved to `S:\COMSOL\mcp-server`) to FEMM and JMAG.

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

### FMM (Fast Multipole Method): Removed from Radia core (2026-03-06)

**ExaFMM-t was removed from the repository**. Do NOT re-implement FMM
acceleration in the Radia C++ core (rad_*.cpp).

**SCOPE CLARIFICATION (2026-05-30, updated 2026-07-05)**: This policy targets
Radia-owned compact Laplace-kernel interaction matrices in the core. It does
NOT forbid the use of FMM-based libraries at OTHER LAYERS of the
stack -- specifically:

  - **ngsolve.bem** (NGSolve's BEM module) ships an FMM-style
    hierarchical Biot-Savart / Laplace / Helmholtz backend
    (``BiotSavartCF``, ``BiotSavartRegularMLCF``,
    ``RegularMLExpansion``, ``SphericalHarmonicsCF``,
    ``IntegralOperator.NearFieldMatrix`` / ``CalcSubMatrix``,
    NGSolve 6.2.2604+).  This is permitted and recommended for
    SF coil design on smooth surfaces where Biot-Savart is the
    natural kernel (free-space, far-field-dominated) -- exactly
    the geometry class where FMM math works well.
  - **HACApK ACA+** remains the choice for Radia's own HDiv-VIM charge-Gram and
    PEEC/BEM interaction matrices when repeated matvecs or compression matter.

The two libraries live at different layers and serve different
geometry classes; using them in their natural domains is not a
policy contradiction.

**WAIT for Hlib (H-matrix) to land in `ngsolve.bem` -- do NOT
roll our own H-matrix at the ngsolve.bem layer (2026-06-19)**:
`ngsolve.bem` today ships only the **FMM** multipole backend
(``BiotSavartCF``, ``*MLCF``, ``RegularMLExpansion``, ...).  An
**H-matrix (Hlib / ACA) backend is expected to land in
`ngsolve.bem`** in an upcoming NGSolve release.  **Until it does,
do NOT bolt a custom H-matrix / ACA path onto the ngsolve.bem
layer** (e.g. wiring Radia's HACApK or
``radia.stream_function.aca_tsvd`` into the SF-coil / BEM
Biot-Savart assembly).  Use the existing FMM backend for the
free-space coil source (apply-once, far-field -- see
archived DtN-spectrum demo_qq/rr sources), and
when a genuine **repeated-apply** H-matrix need arises at that
layer, **prefer the forthcoming native `ngsolve.bem` Hlib over a
hand-rolled integration**.

*Why wait*: a custom H-matrix shim at the ngsolve.bem layer would
(a) duplicate what NGSolve is about to provide natively, (b)
create a maintenance + version-drift burden the moment Hlib lands,
and (c) violate the "No Fallbacks / one supported path" and
"complement NGSolve, do not reimplement" principles.  The
benchmark backing this decision --
archived `bench_fmm_vs_aca_biotsavart.py` (full source: `docs/kelvin/kelvin_dtn_spectrum_archive_results.json`)
(FMM vs ACA+ vs direct on the free-space Biot-Savart coil source,
N=Q 500..8000) -- shows the split the native backends will serve:
**FMM wins one-shot field eval** (low O(N) setup; total 0.63 s vs
direct 4.24 s at N=8000), while an **H-matrix wins repeated apply**
(near-free matvec, 0.015 s vs FMM eval 0.12 s at 24k DoF; rank
~900 constant in N, storage O(N*r) = 14x less than dense).  This
does NOT change the Radia core: **HACApK ACA+ stays Radia's own
HDiv-VIM / PEEC / BEM H-matrix route** per the scope
clarification above; the wait-for-Hlib rule is only about the
`ngsolve.bem` (free-space BEM / SF-coil) layer.

**Why FMM was removed from Radia (the original reasoning, still binding for
compact Radia core matrices)**:

1. **Point multipole approximations are poor for extended element sources**:
   distributed face/volume sources cannot be collapsed to a single dipole at
   intermediate distances (r ~ 2-5 element sizes) without unacceptable
   engineering error.

2. **FMM Solve (Method 3) was useless**: Compact geometries (C-type magnets, iron yokes) have 87% near-field pairs. Near-field correction memory equals the full dense matrix, eliminating FMM's O(N log N) advantage. HACApK (H-matrix, Method 2) is 10-100x faster because ACA+ compression works on the same near-field blocks.

3. **FMM field evaluation had no benefit over direct**: For typical Radia models (N < 10,000 elements), direct B_genComp with TaskManager parallelization is fast enough. FMM overhead (tree build, M2L translation) exceeds direct computation time for these sizes.

4. **HACApK covers all large-scale needs**: H-matrix acceleration (ACA+) provides O(N log N) memory and O(N log^2 N) MatVec for the interaction matrix, which is the actual bottleneck.

**Lesson**: FMM is effective for point charges/dipoles in unbounded
space (N-body), and for smooth surface BEM kernels (= what
ngsolve.bem targets).  It is not Radia's default for compact core
operators with extended element sources and near-field-heavy geometry.
Choose the acceleration
library by GEOMETRY CLASS, not by reflex.

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

### CAD Authoring: build123d or Cubit, NOT Ad-Hoc netgen.occ (2026-07-05)

**POLICY** (Sugahara): for CAD AUTHORING -- parametric solids, booleans,
sweeps, any real geometry construction -- use **build123d** (the lab's hardened
`radia_mcp.build123d` parametric library + the `radia-build123d` MCP) or
**Cubit** -- a FIRST-CLASS CAD authoring environment in its own right (the
robust **ACIS** solid kernel: `create`, `unite`/`subtract`/`intersect`,
`webcut`, `imprint`/`merge`, `sweep`/`tweak` via `.jou` + the lab plugin; NOT
just a mesher).  Do NOT author CAD with ad-hoc **`netgen.occ`** scripting.

`netgen.occ` is a MESH library with only a THIN OCC binding -- fine for the
STEP -> mesh path (read a STEP, mesh it, write `.vol`), but NOT a CAD authoring
kernel.  Its OCC ops are FRAGILE and mislead (all measured 2026-07-05 building
the SF printable-former; see `memory/sf_printable_former_cad_status.md`):
- OCCT `Pipe`/sweep over a self-crossing or sharp-cornered path SEGFAULTS
  (uncatchable);
- `Glue` of overlapping solids SILENTLY returns the input UNCUT (it is
  shared-face assembly, NOT a boolean union);
- `Fuse` of many solids is slow (~O(N^2)) and goes DEGENERATE on self-crossing
  tools;
- a small profile swept along a large path yields a degenerate SHELL
  (volume ~0), not a solid.

**Decision rule**:
- new parametric CAD (magnets, yokes, coils, formers, primitives + booleans)
  -> the DEDICATED **`radia_mcp.build123d`** API: `modeling` verbs (`swept`,
  `coil`, `revolved`, `lofted`, `tube`, `racetrack_coil`, `path_array`, boolean
  `+`/`-`) + `archetypes` (`c_core`, `halbach_ring`, `cos_theta_dipole`,
  `solenoid`, `litz_wire`, ...), run via `execute_build123d`; tested +
  Netgen-meshable + region-labelled.  Do NOT hand-roll raw build123d
  (a missing `.solid()` / wrong profile plane silently yields a degenerate
  shell -- use the helpers, which do `Plane(origin=path@0, z_dir=path%0) *
  profile` then `sweep(...).solid()`).
- CAD authoring is EQUALLY well served by **Cubit** (`.jou` + the lab plugin
  `export netgen`/`gmsh`) -- and it is the route for hex meshing / complex
  topology.  Its **ACIS** kernel is a DIFFERENT, mature solid modeller (NOT
  OCCT), so tricky booleans / webcuts / sweeps that are FRAGILE in OCCT
  (build123d / netgen.occ) may just work in Cubit -- reach for Cubit on hard
  CAD, not only for meshing.  (Pick by fit: build123d for Python-parametric
  families + Netgen-tet; Cubit for ACIS-robust construction + hex.)
- reading/meshing an EXISTING STEP -> `netgen.occ` / OCC is fine (that is I/O,
  not authoring).

**Caveat -- build123d AND netgen.occ BOTH wrap the OCCT kernel**: the dedicated
`modeling.swept`/`coil` work for NORMAL wire/coil ratios (helices, pipes), but
an arbitrary-path SWEEP of THIN-profile / large-path / SELF-CROSSING geometry is
fragile in OCCT ITSELF -- verified 2026-07-05: even `radia_mcp.build123d.swept`
gives a DEGENERATE volume on the SF former's 0.8 mm-wire / 0.15 m-coil arc and
RAISES `StdFail_NotDone BRep_API` on the self-crossing single-stroke wire (the
lab's build123d knowledge already flags "closed-path sweep is fragile in the
OCCT kernel; use `cos_theta_dipole`").  For a ROBUST boolean of THAT
self-crossing geometry (the SF printable-former channel) use a **MESH boolean**
(`trimesh` + `manifold3d`) -> **STL** -- watertight in ~1 s where OCCT
segfaults / degenerates, and STL is the native 3D-print format (shipped as
`calc_streamfunction.py --former-stl`).

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
  `ExportGmshCommand`.  The `version` keyword on `export gmsh` was
  also removed (radia 4.80.0); the current syntax is
  `export gmsh "f" order N [dimension D]` and always emits v4.1.
  An old `.jou` (or test) that still passes `version N` now ERRORS
  ("Unrecognized Identifier: 'version'") -- drop the keyword.  **DECISION
  (2026-05-29): the `version` keyword is ABOLISHED and v2.2-route NGSolve
  loading is NOT supported.**  Do NOT restore accept-and-ignore; v4.1 is the
  only emitted `.msh` format and NGSolve interchange is `.vol`-only (no
  `ReadGmsh`).  (The gmsh-based cohomology path was abolished 2026-06-13:
  `cohomology_cut.py` no longer writes any temporary `.msh` -- the T-Omega
  cohomology CUT is computed gmsh-free by the pure-Python `radia.cohomology`
  engine.  GMSH stays allowed only for visualization / `.msh v4.1` post.)

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

- The deployed LAB / 100号機 / mdx machines all run Cubit Pro and never
  see the warning anyway.

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
  reduced-FEM/NGSolve coupling.
- PEEC: conductor circuit extraction and MQS/Darwin surface-impedance analysis.

**Decision (2026-07-05, Sugahara): Radia soft-iron demag is HDiv-VIM only.**
The public API must not preserve retired collocation or yano backend aliases.
`demag_backend="auto"` and `"hdiv"` are the supported names;
`"auto"` selects the HDiv-VIM path for soft iron. Mesh-backed operation goes
through `.vol` -> NGSolve `Mesh` -> `radia.vim.soft_iron_from_mesh` (or
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

### NGSolve-Native Discretization Policy (2026-07-16)

**POLICY**: Radia owns the missing electromagnetic method; NGSolve owns finite-
element plumbing.  Treat NGSolve as the source of truth for element orientation,
local/global DOF transforms, Piola mappings, curved geometry, quadrature,
weak-form assembly, and `CoefficientFunction` / `GridFunction` evaluation.

- Prefer NGSolve spaces, forms, mapped evaluations, and grid functions over a
  Python reconstruction of physical high-order basis values.
- Never assume `CalcShape` plus `GetDofNrs` is the complete physical transform
  for HDiv/HCurl HEX, WEDGE, curved, or high-order elements; NGSolve may apply
  additional local DOF orientation transforms.
- A Radia C++ reference-element kernel is appropriate when it implements the
  new physical method, but it must take NGSolve-owned geometry/space data and be
  checked independently against an NGSolve weak form or grid-function value.
- Python FE callers own `with ngsolve.TaskManager():`; C++ kernels follow the
  TaskManager self-wrap policy below.

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

### Reduce Proprietary API Surface — Plumbing to netgen/ngsolve, Methods Stay (2026-06-19)

**POLICY**: Extend "Complement NGSolve" to the **API surface itself**:
**aggressively delete Radia-proprietary APIs that duplicate netgen/ngsolve
(the plumbing), and lean on netgen/ngsolve for those parts.**  BUT this is a
**plumbing-vs-method** line, NOT a "uses-ngsolve-vs-not" line — the genuine
numerical METHODS stay (even when built on ngsolve, like `axifem` or the
HDiv-VIM), because they are exactly the value NGSolve does not provide.

**The decision rule (apply every time):**
> "Does netgen/ngsolve already provide this?"
> - **YES (plumbing)** → delete Radia's version, delegate to netgen/ngsolve.
> - **NO (a method NGSolve lacks: open-boundary analytic field, HDiv-VIM,
>   PEEC, Henrotte axisymmetric FE, ...)** → KEEP it (and over time DEMOTE it
>   from a user-facing API to an internal C++/representation detail — keep ≠
>   expose; the un-pybind path is gradual, see the 2-layer API below).

**DELETE / delegate (plumbing — netgen/ngsolve own it):**
- Mesh generation & representation → netgen mesh / `.vol` (e.g. retire
  `create_hex_mesh_grid`; use `MakeStructured3DMesh` / OCC / Cubit).
- Mesh I/O, geometry/CAD kernels → NGSolve / OCC.
- Visualization & mesh export → GMSH / the Cubit plugin (mostly already removed:
  `ObjDrwVTK`, pyvista, `ObjDivMag`…).
- Generic linear algebra → MKL / NGSolve.
- **The geometry PRIMITIVES (`ObjHexahedron`/`ObjRecMag`/…) as the user's
  hand-built-mesh API** → replaced by `.vol` → `soft_iron_from_mesh`/`_from_vol`
  + the intent objects below.  (The primitives are NOT deleted — see KEEP — they
  are demoted to the internal representation.)

**KEEP (genuine methods — Radia's reason to exist; never delete, even on ngsolve):**
the whole `core` set — `rad.Fld` **analytic open-boundary field** (the crown
jewel), **HDiv-VIM**, **`axifem`** (Henrotte
axisymmetric FE — NGSolve has no native Henrotte magnetic basis), **DtN/FEM-Kelvin
operator**, **PEEC**, **BEM** (`sibc_hacapk`…), **sparsesolv** (Compact AMS/AMG/COCR),
**HACApK**, **analytical_formulas**, **coil_builder** (mesh-free Biot-Savart source),
and the application methods (**levitation/ECB**, **stream-function**).  These are
KEPT as internal kernels even as their direct user-facing pybind surface shrinks.

**Target = a 2-layer API:**
- **User layer (pybind, intent-based):** `radia.SoftIron("yoke.vol", mu_r=)` /
  `Magnet(...)` / `CoilBuilder(...)` — "place a magnet / soft iron / coil",
  geometry from `.vol` or simple-shape constructors; backend selected, not
  re-APIed.  Built on NGSolve mesh/geometry where possible.
- **Internal layer (shrinking pybind surface → C++-internal):** the proprietary
  field kernels, demag solvers, and element primitives (`ObjHexahedron`…).  Move
  builders to C++ and un-pybind the primitives **gradually** (CoilBuilder,
  panels, docs notebooks, validation lanes, or tests may depend on them today —
  demote first, remove after migration).

**Soft-iron path**: the SAME `.vol` -> mesh -> HDiv-VIM object -> `rad.Fld`
workflow is the supported path. Backend selectors are limited to the current
HDiv-VIM names and must not carry retired compatibility aliases.

### Maglev Analysis: Radia + NGSolve, Not FEM Alone

**POLICY**: Magnetic levitation (maglev) analysis -- the eddy-current
electromagnetic FORCE between a permanent magnet / coil and a (moving)
conductor -- is solved with **Radia analytic field kernels + NGSolve (FEM) weak
coupling**, NOT with standalone FEM.

**Why pure FEM is the wrong tool for maglev**:
- The PM<->conductor air gap is large and must be meshed; magnet MOTION
  forces re-meshing every step (the dominant cost).
- Air-region discretisation error degrades the force accuracy.
- Open-boundary truncation needs PML or a large air box.

**The Radia + NGSolve method** (lab research line, Yano & Sugahara,
CAE-AI; see `radia_mcp.maglev` topics `radia_iem_fem` / `cln_mor_control`):
- Radia IEM computes the open-boundary external field (A_ext, H_ext)
  ANALYTICALLY -- no air mesh, exact open boundary.
- NGSolve reduced-potential FEM (A-phi or T-Omega) solves ONLY the eddy
  reaction field in the conductor; the Radia field is the source term.
- Weak (sequential) coupling, fed back to demagnetisation.
- **Magnet motion needs only an external-field UPDATE -> NO re-meshing.**
- Optionally, Cauer Ladder Network (CLN) model-order reduction compresses
  the eddy-current FEM into an equivalent circuit for real-time
  control-coupled simulation (~1/500 full-FEM time; TEAM 28).

Validated on the standard eddy-current levitation benchmarks (TEAM
Problem 7 = eddy-current "asymmetrical conductor with a hole" validating
the force/loss solver; TEAM Problem 28 = the electrodynamic levitation
device). Refs: Chadebec 2006 (IEM open boundary), Biro 2000 (reduced
potential), Kameari-Ebrahimi-Sugahara-Shindo-Matsuo 2018 (CLN, IEEE TMag
54(3):7201804).

This is the maglev-specific instance of the "Complement NGSolve" strategy
above: FEM is weak at open boundary + moving magnets + thin conductors;
Radia supplies exactly those (analytic open-boundary field, no air mesh)
and NGSolve supplies the eddy-current solve. Do NOT solve maglev with a
full-FEM air-box model when the Radia + NGSolve weak coupling applies.

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
| **Coreform Cubit** | 2025.12 | Embedded Python 3.10 + PySide6 (Qt6). No Qt5/PyQt5. Cannot import NGSolve/Radia directly. |
| **NGSolve** | 6.2.2606 | BEM correctness, thread-local TaskManager/LocalHeap, curvedelements Load, hex/prism curving, Periodic BC fix |

**Cubit panel subprocess constraint**: Cubit embeds Python 3.10; Radia/NGSolve require 3.12. Same-process import is impossible. All computation runs via `subprocess.run([python3.12, calc_*.py])` with JSON output.

### Build: MSVC + Intel MKL

**POLICY**: Use **MSVC** compiler with **Intel MKL**. Intel oneAPI compiler (icx-cl) is NOT compatible with NGSolve linking.

```powershell
powershell.exe -ExecutionPolicy Bypass -File "Build.ps1"
powershell.exe -ExecutionPolicy Bypass -File "Build.ps1" -Rebuild  # Clean rebuild
```

**Required Software**: Visual Studio 2022 (MSVC), Intel oneAPI Base Toolkit (MKL only, NOT the compiler).

### BLAS/LAPACK: Intel MKL Only

**POLICY**: OpenBLAS is NOT supported. MKL provides optimized BLAS/LAPACK. MKL internally uses Intel OpenMP (`libiomp5md.dll`) for its own threading, but Radia no longer links it directly.

**Required MKL DLLs** (loaded at runtime via pip dependency): `mkl_rt.*.dll`, `mkl_core.*.dll`, `mkl_intel_thread.*.dll`, `mkl_def.*.dll`, `mkl_avx2.*.dll`, `mkl_vml_*.dll`, `libiomp5md.dll` (MKL dependency), `libmmd.dll`, `svml_dispmd.dll`.

### Parallelization: NGSolve TaskManager

**POLICY**: Use **NGSolve TaskManager** for thread-level parallelization, NOT raw OpenMP parallel regions.

NGSolve's TaskManager provides work-stealing task-based parallelism that integrates with MKL and avoids nested OpenMP issues. All new parallel code in Radia should use TaskManager.

```cpp
// CORRECT: NGSolve TaskManager
#include <ngstd.hpp>
TaskManager::CreateJob([&](const TaskInfo& ti) {
    // work-stealing parallel loop
    for (size_t i = ti.task_nr; i < n; i += ti.ntasks) {
        // compute...
    }
});

// AVOID: Raw OpenMP parallel for (legacy code only)
#pragma omp parallel for
for (int i = 0; i < n; i++) { ... }
```

**When to use TaskManager**:
- Field computation loops (ComputeFieldBatch)
- Interaction matrix assembly
- Any embarrassingly parallel loop

**When OpenMP is acceptable**:
- MKL internal threading (controlled by `mkl_set_num_threads`)
- Legacy code not yet migrated

### TaskManager Wrap Policy: Caller Wraps, Helper Does NOT (2026-05-27)

**POLICY**: The TaskManager wrap is the **caller's responsibility**.
Helper functions / solver classes / library modules that perform
NGSolve operations (`.Assemble()`, `.Inverse(...)`, `mesh.Curve(p)`,
`Integrate()`, `GridFunction.Set(...)`) MUST NOT include `with
TaskManager():` internally.  The **caller** (a `calc_*.py` panel
script, a `validation_test/**.py` driver, a test, or a notebook cell) is
the one that opens the `with TaskManager():` context once and lets
all helper calls inside the region run in parallel.

**Rationale**:

1. **Composability**: a caller that runs `helper_a()` then
   `helper_b()` then `helper_c()` in sequence should open ONE
   `TaskManager` region for the whole batch, not pay the
   start/stop overhead three times.  When helpers wrap internally,
   the caller cannot achieve this.
2. **Predictability**: a reader of `calc_*.py` can see the
   parallelism intent at the call site, not buried inside a helper
   they have to chase.
3. **Single audit point**: the parallel-correctness audit becomes
   "grep `with TaskManager():` in `calc_*.py`, `validation_test/**.py`,
   `tests/**.py`, and docs notebook helpers".
   No need to chase helper modules.
4. **Removes silent double-wrap noise**: `with TaskManager(): with
   TaskManager(): ...` is a no-op for the inner context (NGSolve
   detects + reuses), but it is visual noise that obscures the
   actual intent.

**Concrete rules**:

- Helper modules under `src/radia/**.py` (NOT `panels/calc_*.py`)
  MUST NOT contain `with TaskManager():`.  If they need to assert
  the caller has wrapped, use the diagnostic at the top of the
  helper (see "Helper diagnostic" below).
- Caller modules MUST wrap.  `tools/audit_taskmanager.py` enforces
  this check repo-wide; a violation = parallelism bug.
- Custom C++ kernels using `ngcore::ParallelFor` (e.g.
  `rad_equivalence_source.cpp`) are NOT helpers in this sense —
  they are leaf parallel kernels and naturally honour the caller's
  TaskManager context via the runtime.
- `validation_test/**.py`, `tests/**.py`, docs notebook helpers, and panel
  `calc_*.py` follow the same rule as callers: every script
  that does `.Assemble()` / `.Inverse(inverse=...)` / `mesh.Curve(p)`
  MUST wrap in `with TaskManager():`.

**Helper diagnostic** (optional belt-and-suspenders pattern that
some helpers may use to fail loudly when the caller forgot to
wrap):

```python
def _check_task_manager_active():
    """Yell if the caller forgot to open a TaskManager context.
    NGSolve's actual API does not expose an "am I in TM" flag, so
    this is a best-effort timing/affinity check; helpers can SKIP
    this if benchmarking the helper itself never runs serial."""
    import ngsolve
    # If global threads is 1 (serial), the caller may have set
    # SetNumThreads(1) for debug -- DO NOT raise.  Only raise if
    # threads >1 but the user is clearly running outside a TM block
    # (heuristic: ngsolve.ngsglobals.task_manager is None).
    ...  # implementation-specific
```

**Audit tool**:

```bash
python tools/audit_taskmanager.py
# Exit 0 = clean; exit non-zero = list of violations.
```

The audit checks:
1. Helper modules (anything in `src/radia/` except `panels/calc_*.py`)
   have ZERO `with TaskManager():`.
2. Caller modules (`panels/calc_*.py`, `validation_test/**.py`,
   `tests/**.py`, docs notebook helpers) that
   contain `.Assemble()` / `.Inverse(.*inverse=` / `mesh.Curve(`
   DO have `with TaskManager():` at the top of the containing
   function.

Run after editing any solver helper, validation driver, panel calc script, or
docs notebook helper.  The check is
fast (~1 s repo-wide); add to pre-commit if desired.

**Migration note (2026-05-27)**: The 7 helpers under `src/radia/`
that originally had internal `with TaskManager():` were converted
to no-internal-wrap.  Historical example callers were swept to add
the caller-side wrap before promotion/deletion.  See
`taskmanager('helper_vs_caller')` in
the radia-mcp MCP knowledge for the full policy + audit history.

### TaskManager-Only Policy: Align with NGSolve, No Alternatives (2026-05-27)

**POLICY**: TaskManager is an **NGSolve-side mechanism**; Radia
aligns with it as the **sole** parallelization path for
NGSolve-driven computation.  Do NOT introduce any other
parallelization mechanism for code that touches NGSolve --
no raw `#pragma omp parallel`, no `threading.Thread`, no
`multiprocessing.Pool`, no `concurrent.futures` for numerics.
NGSolve's `with TaskManager():` is the canonical entry point;
Radia C++ kernels reach the same threadpool via `ngcore::ParallelFor`
(which honours the active TaskManager context).

**Even ops NGSolve has not yet internally parallelized go through
TaskManager.**  Examples: `BilinearForm(...)`, `LinearForm(...)`,
`H1(mesh, ...)`, `HCurl(...)`, `HDiv(...)`, `Periodic(...)`,
`GridFunction(fes)`, `Mesh(geo.GenerateMesh(...))`, `Integrate(...)`.
Some of these are mostly serial today; future NGSolve releases that
add parallelism to FES construction, mesh generation, or quadrature
will benefit automatically without code changes on Radia's side.

**Why this expansion (beyond just `.Assemble()` / `.Inverse()` /
`mesh.Curve()` / `gf.Set(cf)`)**:

1. **One coherent computation block** — readers see one
   `with TaskManager():` that brackets "the NGSolve work", not a
   patchwork of micro-wraps around isolated parallel calls.
2. **Future-proof** — when NGSolve internally parallelizes a new
   op, callers get the speed-up the day they `pip install --upgrade
   ngsolve`, without a Radia-side audit pass.
3. **No "is this op parallel?" gotcha** — contributors do not need
   to memorize which NGSolve calls happen to be parallel today.
4. **Aligns with NGSolve's own docs/examples** — the canonical
   NGSolve tutorial pattern wraps the whole solve block, not
   individual lines.

**Concrete shape** (the canonical script template):

```python
from ngsolve import *
ngsolve.SetNumThreads(args.nthreads)   # process-wide thread cap
with TaskManager():
    # ----- mesh -----
    mesh = Mesh(geo.GenerateMesh(maxh=h))
    mesh.Curve(p)

    # ----- FES + forms -----
    fes  = H1(mesh, order=p, dirichlet="outer")
    u, v = fes.TnT()
    a    = BilinearForm(fes); a += grad(u)*grad(v)*dx; a.Assemble()
    f    = LinearForm(fes);   f += rhs*v*dx;           f.Assemble()

    # ----- solve -----
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    # ----- post -----
    err = Integrate(InnerProduct(grad(gfu)-grad_u_exact,
                                 grad(gfu)-grad_u_exact), mesh)
    print(f"H1 error = {sqrt(err):.3e}")
```

The wrap covers FES construction and post-processing too, even
though those are not currently parallel-hot.  This is intentional
per the "future-proof" rationale above.

**Audit expansion**: `tools/audit_taskmanager.py` now also flags
`Integrate(`, `BilinearForm(`, `LinearForm(`, FES constructors
(`H1`, `HCurl`, `HDiv`, `L2`, `VectorH1`, `Periodic`), and
`Mesh(geo.GenerateMesh(...))` outside a TaskManager region.

**Not affected** (these are NOT NGSolve calls and do NOT need
wrapping):

- Pure numpy / scipy operations
- `radia.*` functions (Radia's own C++ kernels honour the active
  TaskManager context via `ngcore::ParallelFor`; you still wrap
  them, but they are not what triggers the audit)
- Pure I/O (file read/write, json, ...)
- argparse, logging, print

**Forbidden** (use TaskManager instead):

- `import threading` in any numeric loop
- `import multiprocessing` for solving / assembly
- `concurrent.futures.ThreadPoolExecutor` for NGSolve work
- `#pragma omp parallel` in new C++ that calls NGSolve

### C++ HACApK Self-Wrap Policy: every BUILD and SOLVE-LOOP stands up a RegionTaskManager (2026-06-23)

**POLICY**: The "Caller Wraps, Helper Does NOT" rule above governs **Python**
helpers.  The **C++ HACApK kernels are the opposite**: the H-matrix leaf fill and
the H-matvec call `ngcore::ParallelFor`, which **silently falls back to
single-threaded when NO `RegionTaskManager` region is active**.  A non-panel
caller of `rad.Solve(..., method=2)` (or a bare `hdiv_demag_solve`) does NOT open
a `with TaskManager()` region, so without a C++ self-wrap the **entire HACApK
build + Krylov solve runs serial** -- the exact failure mode that made a retired
demag method-2 path silently single-threaded off the panel path (2026-06-23).

**RULE**: every C++ entry point that drives an HACApK `ngcore::ParallelFor` MUST
self-wrap **exactly once, at the right granularity**:

- **BUILD** -- `RadHACApKBase::BuildHMatrix` stands up the region at the top of
  its body.  This is the **single central** protection for EVERY build site
  (HDiv ChargeGram, PEEC, BEM, the diagnostic densify/smoke entries).  Do NOT
  duplicate it per build caller.
- **SOLVE LOOP** -- each Krylov / iterative loop (BiCGSTAB, CG, MINRES, Picard)
  that repeatedly calls `MatVec` stands up **ONE** region around the whole loop.
  **NEVER wrap inside `RadHACApKBase::MatVec` itself**: a per-matvec region would
  stand up / tear down the threadpool on every iteration (the very case we are
  fixing), far slower than wrapping the loop once.
- **One-shot `MatVec`** (single apply / smoke test) needs no wrap -- the build it
  follows is already covered, and a lone matvec gains nothing from a persistent
  pool.

**The reference pattern** (use verbatim):
```cpp
ngcore::RegionTaskManager rtm(std::max(1, ngcore::TaskManager::GetMaxThreads()));
```

**Nesting is always safe**: `RegionTaskManager` **reuses the caller's pool when
one is already active** (no-op), so a panel that already opened `with
TaskManager()` is unaffected, and a build's self-wrap nesting inside a solver's
self-wrap is a no-op.  The self-wrap is belt-and-suspenders, never a conflict.
This does NOT contradict "Caller Wraps, Helper Does NOT" (that rule is about
**Python** helpers); C++ HACApK kernels self-wrapping with `RegionTaskManager`
is the established C++ pattern -- already used by
`RadHACApKChargeGram::SolveMaterialMINRES` (`rad_hacapk_hdiv.cpp`) and the H-LU
bridge `chacapk_par_region` (`cHACApK_harith_par.cpp`).

**Self-wrapped surface (2026-06-23 sweep)** -- BUILD: `RadHACApKBase::BuildHMatrix`;
SOLVE LOOPS: `RadHACApKChargeGram::SolveLinearMaterial` / `SolveNonlinearPicard` (HDiv;
`SolveMaterialMINRES` already had it), and the `HMatrixDensify` densify loop.
**Any NEW HACApK solver loop or build path MUST add the self-wrap.**  (The stale
`OpenMP`-worded comments in the HACApK callback state were corrected to
`TaskManager` in the same sweep -- HACApK has zero OpenMP; its parallelism IS
TaskManager via `hacapk_parallel_for` = `ngcore::ParallelFor`.)

### In-Job Code Must Not Nest CreateJob; Evaluation Must Not Mutate Shared State (2026-07-11)

**POLICY**: the self-wrap pattern above is ONLY valid when no TaskManager job
is running.  `ngcore::TaskManager::CreateJob` keeps STATIC job state, so a
`ParallelFor` issued from INSIDE a running job (e.g. a CoefficientFunction's
`Evaluate` called by NGSolve assembly workers) corrupts the running job --
hard process death (0xC0000374 / 0xC0000005), no Python traceback.  Three
rules, learned from the 2026-07-10/11 RadiaField heap-corruption incident:

1. **C++ code that can execute inside a TaskManager job must use SERIAL
   loops** and let the outer job own the parallelism.  The RadiaField CF
   evaluates via `RadFldBatchSerial` / `RadFldPhiSerial` / `RadFldASerial`
   (`radTApplication::Compute*Batch(..., parallel=false)`); the parallel
   self-wrapped entries (`RadFldBatch` etc.) are for callers with no active
   job (the Python `rad.Fld` path).  A new batch C API needs both variants.
2. **CoefficientFunction evaluation must be GIL-free**: no
   `py::gil_scoped_acquire` / Python round-trip inside `Evaluate` --
   assembly workers + GIL save/restore corrupt the interpreter.  Call the
   C API directly (single-point `RadFld` is also banned there: it
   round-trips results through the non-thread-safe global `ioBuffer`).
3. **Field evaluation must never mutate shared element state**: the
   mutating `radTGroup::B_genComp` / `radTPolyhedron::B_genComp` overrides
   were removed (groups/polyhedra inherit the transform-via-point
   `radTg3d::B_genComp`; the polyhedron override had also silently IGNORED
   TrfOrnt on hex/tet/wedge elements), and the `radTHandle` refcount is now
   `std::atomic<int>` (concurrent handle copies during parallel evaluation
   lost increments and deleted live objects mid-assembly).

Regression lock: `tests/test_radiafield_transformed_container.py` (crash
paths exercised in subprocesses; values locked vs the rotated-reference
solution and pointwise `rad.Fld`).

### PyPI Release Workflow (Automated via GitHub Actions)

**POLICY**: PyPI publishing is automatic. Push a version tag (`v*`) and CI/CD handles the rest.

**Release Flow**:
1. Bump version in `pyproject.toml` AND `src/radia/__init__.py` (must match)
2. Update `CHANGELOG.md`
3. `git commit` (do NOT push yet)
4. `/deploy` — build wheel, deploy to 100号機 (WinRM) & mdx (SSH)
5. Test on remote machines (Cubit toolbar, Simulink application blocks, result-bearing docs notebooks, Mesh Evaluation, etc.)
6. If tests pass: `git push origin main`
7. **Confirm main CI is GREEN before tagging** (gh-free; `gh` is not on LAB):
   `python tools/release_qud.py ci-verify` — waits for the self-hosted
   runner job to finish, then checks the workspace junit XMLs
   (failures=errors=0). Tag CI = the same `build-test.yml` on the same commit,
   so a green main CI guarantees a green tag CI — and avoids burning a version
   number on a broken commit (the v4.80.0→v4.80.5 saga). If RED, fix-forward
   on main and re-run.
8. Tag and push **only after step 7 is green** (triggers PyPI publish):
   ```bash
   git tag vX.Y.Z
   git push origin vX.Y.Z
   ```
9. Monitor PyPI propagation (gh-free): `pip index versions radia` (and
   `radia-mcp`, `cubit-mesh-export`). The tag CI just re-runs the
   already-green main commit, so it is expected green.

**General User Install** (after PyPI publish):
```bash
pip install radia[cubit]
cubit-plugin-install       # Cubit plugin + panels
```

**CI/CD Pipeline** (`.github/workflows/`):
```
git push v* tag
  -> CI (build-test.yml): Build.ps1 -> pytest -> Build_Wheel.ps1 -DryRun -> upload artifacts
  -> Release (release.yml): download wheel artifact -> pypa/gh-action-pypi-publish (OIDC Trusted Publishers)
```

**No API tokens stored**. Uses PyPI OIDC Trusted Publishers (id-token: write).

**NGSolve on CI runner**: The self-hosted runner (NETWORK SERVICE) cannot access S: drive. NGSolve must be copied locally:
```powershell
robocopy S:\NGSolve\01_GitHub\install_ngsolve C:\NGSolve /MIR
```

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

**POLICY (2026-06-25 update)**: **QUD 配布**: LAB と 100号機 は editable
(NAS source、developer/user feedback loop)、mdx と hibino は PyPI install。
mdx は compute/Cubit verification point なので `radia` + `cubit-mesh-export`
のみでよく、`radia-mcp` は不要。hibino は PyPI 経由の MCP consumer として
`radia-mcp` も入れる。

**QUD 配布モデル (2026-06-25)**:

| Stage | マシン | install 形態 | 目的 |
|-------|--------|-----------|------|
| 1 | LAB | `pip install -e .` + `pip install -e packages/cubit-mesh-export` + `pip install -e packages/radia-mcp` | 開発者ループ。最速フィードバック (NAS source 直接編集) |
| 1 | 100号機 | `pip install -e \\192.168.11.100\work\00_CAE\Radia\01_GitHub` + `pip install -e ...\packages\cubit-mesh-export` + `pip install -e ...\packages\radia-mcp` | 共有ユーザ環境も editable。radia と MCP の学習/修正が即反映される |
| 2 | mdx | `pip install radia / cubit-mesh-export` (PyPI) + `cubit-plugin-install --all-users` | PyPI wheel + Cubit plugin の compute verification。`radia-mcp` は不要 |
| 2 | hibino | `pip install radia / radia-mcp / cubit-mesh-export` (PyPI) + `cubit-plugin-install --all-users` | PyPI 経由の MCP consumer verification |

**変更点 (2026-06-25)**:
- 旧: LAB editable / 100号機 + mdx 両方 PyPI (2-tier).
- 新: LAB + 100号機 editable / mdx + hibino PyPI (QUD).
- mdx は `radia-mcp` 不要。hibino を PyPI 経由の MCP consumer として追加。

**LAB / 100号機 editable パッケージ**:
- `radia` (LAB: `S:\Radia\01_GitHub`, 100号機: `\\192.168.11.100\work\00_CAE\Radia\01_GitHub`)
- `cubit-mesh-export` (`packages\cubit-mesh-export`)
- `radia-mcp` (`packages\radia-mcp`)
- `mcp-server-document` (LAB: `S:\mcp-server`) -- LAB-private (PyPI 配布なし)

共有worktreeに並行WIPがあるreleaseでは、そのWIPをstash/clean/resetしてはならない。
代わりに同一NAS上へexact SHAのclean release worktreeを作り、
`RADIA_RELEASE_EDITABLE_REPO_LAB` と `RADIA_RELEASE_EDITABLE_REPO_100` で
LAB/100号機から見える各pathを `release_qud all` と `done` の両方へ渡す。
QUDはprocess停止やinstallより前にSHA一致とtracked-cleanを強制する。公開完了まで
release worktreeを保持し、並行WIPがmainへ着地した後に通常のeditable pointerへ戻す。

LAB / 100号機で `pip install --upgrade <pkg>` を流すと editable が静かに上書きされて壊れるので注意。release 後の LAB / 100号機 側 metadata 同期は `pip install -e <path> --no-deps --no-cache-dir` で再 editable 化。`pip install --upgrade` は **mdx / hibino 用** (PyPI から通常通り upgrade).

**POLICY (2026-05-27 追加): release 後の LAB editable 再確認**

PyPI release (tag push → CI publish) 後、**LAB の editable pointer
が drift していないか必ず確認する**。これは `release-qud` skill
の Definition of Done の暗黙の前提条件であり、CLAUDE.md「LAB editable
default」原則を実運用で守るチェック。

**いつ実施するか** (再発した historical incidents):

| 日付 | 何が起きたか | 検出経緯 |
|---|---|---|
| 2026-04-28 | `pip install --upgrade` で 3 パッケージ全部の editable が上書き | dev loop が機能しないことで発覚 |
| 2026-05-19 | 同じ pattern 再発 | LAB から「source 編集が反映されない」報告で発覚 |
| 2026-05-26 | release 直後の CI runner clone (`C:\actions-runner\_work\...`) に `radia-mcp` editable が drift | LAB 側で knowledge file の編集が radia-mcp 経由で見えないことで発覚 |
| 2026-05-27 | 同上、再発 | MCP tool 経由で新規 topic が dispatch されない (Unknown topic) で発覚 |
| 2026-08-06 | v4.95.47 release 後、radia + cubit-mesh-export + radia-mcp の 3 つ全部が **release-qud worktree** (`release-qud/Radia-v4.95.47-<sha>`) に drift(新パターン: CI runner ではなくリリース用worktreeで editable install が走った) | release 直後の本チェックで発覚。**重要教訓: 同バージョンへの `pip install -e` 再実行は `.pth` を書き換えない no-op になり得る — `pip show` の Editable location が直っても import は旧パスのまま**(実測: show=canonical tree なのに `__editable__.*.pth` は release-qud、radia.gmsh_post_export が旧ソースを解決)。修復は `pip uninstall -y <pkgs>` → 再 `-e`。**検証は `pip show` でなく site-packages の `__editable__*.pth` の中身と `<pkg>.__file__` で行うこと** |
| 2026-08-07 | v4.95.48 release 後、同じ release-qud クローンパターンで 3 つとも再 drift(`Radia-v4.95.48-a6a5ddde`) | **ソース修正が実行時に反映されない**ことで発覚(coil_builder の OCC pose 修正を入れたのに `write_step` の出力が 1 バイトも変わらなかった)。**教訓: `sys.path.insert` するテストは自分のツリーを見るので緑のまま、素の `import radia` を使うデモ/スクリプトだけが旧ソースを掴む** — テストが通っていることは editable が正しい証拠にならない。デモの結果が「直したはずなのに変わらない」ときは真っ先に `<pkg>.__file__` を見る |

**再発するため、release 後の確認を policy 化**。

**チェック手順** (release-qud Phase 8 / Phase 9 の直後に流す):

```powershell
# 4 パッケージ全部の editable pointer が LAB source を指しているか確認
pip show radia cubit-mesh-export radia-mcp mcp-server-document |
  Select-String -Pattern "^(Name|Version|Editable)"

# 期待:
#   radia               -> S:\Radia\01_GitHub
#   cubit-mesh-export   -> S:\Radia\01_GitHub\packages\cubit-mesh-export
#   radia-mcp           -> S:\Radia\01_GitHub\packages\radia-mcp
#   mcp-server-document -> S:\mcp-server
```

**Drift 検出時の修復手順**:

```powershell
# 1. mcp-server-* プロセスを止める (editable uninstall の file lock を解放)
Get-Process | Where-Object { $_.Name -like "mcp-server*" } |
  Stop-Process -Force
Start-Sleep -Seconds 3

# 2. drift しているパッケージを uninstall + LAB source で再 editable 化
pip uninstall -y <pkg>
pip install -e S:\Radia\01_GitHub\packages\<pkg> --no-deps --no-cache-dir

# 3. 確認
pip show <pkg> | Select-String "Editable project location"
```

**自動化候補** (TODO): `tools/verify_lab_editable.py` を作って
`release-qud done` から呼び出し、drift があれば exit non-zero。
今は手動チェックで運用。

**Drift する原因** (analysis):

1. **CI runner と LAB が同じ NAS source を共有**: `\\192.168.11.100\work\00_CAE\Radia\01_GitHub` を `S:\` で参照 (LAB) または UNC で参照 (CI runner)。CI が editable install を走らせると、CI 側の pip metadata が LAB の Python の site-packages にも書き戻されるケースがある (NAS-mounted Python env でない限り通常起こらないが、`pip install -e .` を CI で実行すると `.egg-info` 等が source tree に書かれ、その後の `pip show` 解決順序を狂わせる)。
2. **`pip install --upgrade <lab-pkg>` を LAB で実行**: editable 上書き。**禁止**。
3. **release-qud Phase 8 の PyPI install command**: mdx / hibino で実行されるべきコマンドを誤って LAB / 100号機 shell で実行。

**予防策**:

- LAB shell では決して `pip install <lab-pkg>` (non-editable) / `pip install --upgrade <lab-pkg>` を打たない
- release-qud skill の deploy commands は対象マシンの tier に従って実行する。LAB / 100号機は editable、mdx / hibino は PyPI。
- 不安な場合は release 直後に上記チェック手順を流す

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

**POLICY**: CI/CD のテストだけでは不十分。Cubit が必要な機能（`export_curved`, panels, BEM extractor）は **Cubit 環境でのローカルテストが必須**。CI は C++ ビルドと基本テスト（Cubit 不要なもの）のみ。

**リリース前の必須テスト**:
1. CI: C++ ビルド + pytest（Cubit 不要テスト）
2. ローカル: Cubit + system Python で `export_curved` テスト（球、トーラス）
3. ローカル: Cubit パネルの動作確認

CI が通っても Cubit テストに通らなければリリースしない。

### CI Preflight Policy: commit → CI check → push (2026-06-05; ownership 2026-08-17)

**OWNERSHIP (2026-08-17): every agent that pushes runs `ci_preflight` first —
Claude included.**  Per *Agent Division of Labor: Claude Owns Push + CI; Codex
Owns Release*, Claude implements → tests locally → commits this-session files by
name → runs `tools/ci_preflight.py` → pushes → watches CI to green.  The
pre-push hook runs the preflight anyway; running it deliberately means finding
the breakage before the push rather than being stopped by it.

**POLICY**: Run the CI gates **LOCALLY before every push to `main`**, not
after.  CI must not be the FIRST place a catchable error surfaces.  The
single command is:

```bash
python tools/ci_preflight.py          # ~2-3 min; exit 0 = safe to push
python tools/ci_preflight.py --fix    # also auto-regenerate a stale TOOLS.md
python tools/ci_preflight.py --full   # also run the full top-level pytest (slow)
```

**Why** (empirical, 2026-06-05 analysis of the last 80 failed GitHub
Actions runs): the failures are concentrated and *all locally
detectable before push* —

| count | workflow / step | class |
|---|---|---|
| 30x | radia-mcp matrix `Pytest` | heavy-import collection / meta-health / version |
| 16x | CI `Run basic tests` | top-level pytest import / golden / flaky |
| 9x | radia-mcp matrix `TOOLS.md drift gate` | committed TOOLS.md ≠ regenerated |
| 7x | Policy Lint `Policy 4 CblasColMajor` | C++ allowlist miss |
| 3x | radia-mcp matrix `Meta health` | catalog import / links / tags |

`tools/ci_preflight.py` mirrors the three CI workflows
(`policy-lint.yml`, `radia-mcp-matrix.yml`, `build-test.yml` "Run basic
tests"), fast-first: policy-lint (7 policies) → version consistency →
TOOLS.md drift (**WIP-aware**: warns when `radia_mcp/src` has uncommitted
`.py` changes that would contaminate the regenerated inventory) →
**radia-mcp matrix under minimal-dep simulation** (`RADIA_MCP_FORCE_MINIMAL=1`,
which reproduces the ubuntu collection on a full-env box — this catches the
ngsolve/netgen module-import collection break) → top-level collect-only.

**Enforcement (pre-push hook)**: `python tools/install_git_hooks.py`
installs `tools/git-hooks/pre-push`, which runs `ci_preflight --since
<remote-sha>` on every push to `main` (path-aware, so a non-radia-mcp
push stays fast) and **aborts the push if a gate is red**.  Emergency
bypass: `CI_PREFLIGHT_SKIP=1 git push`.  Run the installer once per clone
(it is worktree-safe and idempotent).

**The recurring CI failure classes are cataloged** in `bug_patterns.py`
(`bug_patterns_lookup(topic="ci")`): `tools-md-drift-wip-contamination`,
`heavy-import-collection-break-minimal-dep-matrix`,
`init-py-version-mismatch-vs-pyproject`, etc. — each names ci_preflight as
the detection tool.

**Release flow** (release-qud) keeps its own release gates, but
ci_preflight is the everyday "before any push" gate — run it even for a
non-release push to `main`.

### No GitHub CLI (gh) Policy (2026-06-05)

**POLICY**: Do NOT use the GitHub CLI (`gh`) anywhere in the repo's
committed tooling, scripts, workflows, or knowledge recipes.  `gh` is
NOT installed on LAB (retired 2026-05-24) and must not be assumed
present on any machine (fresh clone, CI runner, 100号機, mdx).  Use the
**raw GitHub REST API** (via `urllib`, no extra dependency) instead.

**Established gh-free replacements** (use these, never `gh`):

| `gh` command | gh-free replacement |
|---|---|
| `gh run list` / `gh run view` / `gh run watch` | `python tools/check_ci.py [--sha X] [--branch main] [--watch]` |
| `gh release download -p PAT` | `python tools/download_release_asset.py --pattern PAT --dest D` |
| `gh release download NAME`   | `python tools/download_release_asset.py --name NAME --dest D` |
| `gh release upload`          | `python tools/upload_release_asset.py --tag T --file F` |
| `gh api <path>`              | `tools/gh_api.py::gh_get(path)` (token-aware) |
| `gh auth token`              | read `$GH_TOKEN` / `$GITHUB_TOKEN` / `~/.radia/gh_token` directly |

**Auth**: the REST helpers (`tools/gh_api.py`, `check_ci.py`,
`download/upload_release_asset.py`) read a Personal Access Token from
`$GH_TOKEN` / `$GITHUB_TOKEN` (or `~/.radia/gh_token`) to get the
authenticated **5000 req/hr** limit; with no token they fall back to the
anonymous **60 req/hr** (enough for occasional checks on this public
repo).  Setting the token does NOT require installing `gh` — a PAT is
just a string.  `python tools/check_ci.py --rate` shows the current
limit + whether a token was found.

**Escape hatch**: if `gh` is ever genuinely unavoidable for a one-off
local task, it MAY be installed — but the **committed** tooling MUST
stay gh-free so a fresh clone / CI runner never depends on it.  When you
catch yourself reaching for `gh`, add the missing capability to the REST
helpers above instead.

### Cubit Batch Self-Testing Policy

**POLICY**: Claude は Cubit を **完全ヘッドレス** (`-batch -nographics -nojournal`) で起動し、自力で機能試験を走らせること。GUI や人間の操作は不要。

**起動方法** (既存テストのパターン):
```python
import cubit
cubit.init(['cubit', '-nojournal', '-batch', '-nographics',
            '-commandplugindir', <plugin_dir>])
cubit.cmd("create sphere radius 0.05")
cubit.cmd("mesh volume 1")
cubit.cmd('export femeem "C:\\tmp\\cub" overwrite')
```

**対象**: `export {gmsh,netgen,jmag_nastran,vtk,femeem}`、panels の非 GUI ロジック、BEM extractor、`export_curved`。
GUI が絶対必要なもの (panel dialog のレンダリング) のみ例外。

**前提**: `cubit` は Python API import (`S:/Radia/01_GitHub/src/radia/install_panels.py` の `find_cubit_bin()` で自動検出可)。バッチ起動でライセンス消費あり。

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

**CRITICAL**: Use NGSolve **6.2.2606** as pinned by `pyproject.toml`. Required for the validated BEM fixes, thread-local TaskManager/LocalHeap behavior, curvedelements .vol Load, hex/prism curving, Periodic BC fix, and the `ngsolve.bem` FMM-based hierarchical BEM APIs.

Reference: https://forum.ngsolve.org/t/ngsolve-periodic-boundary-condition-regression-bug-report/3805

Official PyPI ngsolve **6.2.2606** includes the Periodic BC fix,
**curvedelements Save/Load**, **p-version hex/prism curving**, and an **FMM-style hierarchical Biot-Savart / Laplace / Helmholtz backend in `ngsolve.bem`**. The wheel uses `ngsolve-openblas`; optional PARDISO support is selected when the separate `mkl` package is installed. The FMM backend is appropriate for SF coil design on smooth surfaces (free-space Biot-Savart on plane / cylinder / sphere); see the "FMM Removed from Radia core (2026-03-06)" policy section for the scope clarification -- Radia's own HDiv-VIM / PEEC / BEM repeated-apply matrices remain on HACApK ACA+ where compression is the right tool.

Installed on LAB 2026-08-25: `pip show ngsolve` -> `Version: 6.2.2606`, `Location: C:\Program Files\Python312\Lib\site-packages`. Keep the exact `ngsolve==6.2.2606` and `netgen-mesher==6.2.2606` pins synchronized with every native rebuild.

**Netgen fork is no longer required.** The ksugahar/netgen repository is historical only.
All curvedelements, CallbackGeometry, and curving features are now in the official release.

```bash
pip install radia[cubit]       # Installs everything (NGSolve, MKL, Cubit plugin binaries)
cubit-plugin-install           # Deploys Cubit plugin + panels
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
- `cubit_mesh_curver` — C++ pybind11 module (bundled in cubit_mesh_export, unchanged)
- `check_vol_consistency` — thin backward-compat re-export in `src/radia/panels/` (imports from cubit_mesh_export.check)

Cubit workflow for journal files: define blocks before export, use the Cubit plugin commands (`cubit.cmd('export gmsh/jmag_nastran/vtk ...')`). Requires `CUBIT_PLUGIN_DIR` environment variable (set by `cubit-plugin-install`).

### PEEC Conductor Mesh

PEEC conductors use **surface mesh only** (SIBC handles skin effect). Generate surface meshes via Netgen or Cubit. Supported: Tri3, Quad4 (1st order), Tri6, Quad8/9 (2nd order).

### Nastran Format: REMOVED

Nastran BDF support is **REMOVED**. Use Cubit -> `.msh` export or Netgen direct. Cubit can read legacy `.bdf` files if needed.

### Mesh Operations: Dropped APIs

`ObjDivMag`, `ObjDivMagPln`, `ObjCutMag` are NOT supported. All mesh operations use external tools (Netgen, Cubit).

### Mesh File Preservation

**NEVER DELETE** mesh files (`.bdf`, `.nas`, `.msh`, `.vtk`), Cubit journal files (`.jou`), or mesh generation scripts. These are difficult to recreate.

### Available Mesh Access Functions

From `src/radia/netgen_mesh_import.py`:
- `netgen_mesh_to_radia()` -- Convert entire mesh to Radia (recommended)
- `extract_elements()` -- Extract element data for custom processing
- `compute_element_centroid()` -- Centroid from vertex list
- `create_radia_tetrahedron()` / `create_radia_hexahedron()` -- Single elements
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

### H-Matrix Route Policy: HDiv/BEM/PEEC Use HACApK (2026-07-05)

HACApK is canonical for the Radia routes whose operators are designed as
H-matrix problems. Retired collocation demag H-matrix adapters are not a supported
surface and must not be used as a backend-selection reason.

HACApK is canonical for the routes where the operator is designed as an H-matrix problem:

1. **HDiv-VIM**: charge-Gram H-matrix and HDiv system managers are the production
   large-scale / loop-free path.  This is the route for high-mu, loop-heavy, and
   FEEC-coupled demag problems.
2. **BEM and PEEC**: the scalar Galerkin BEM and PEEC managers remain maintained
   HACApK adapters.
3. **Point-kernel / Gauss tools**: point-kernel HACApK support is kept as a
   building block for the Gauss-point charge-Gram direction.

**HACApK-connected surface currently kept**:

| Subclass / adapter | Role | Python entry / test surface |
|---|---|---|
| `RadHACApKBEMManager` | BEM / SIBC Laplace Galerkin | `radia.bem_sibc_solver` and BEM HACApK matvec validation |
| `RadHACApKChargeGram` | production HDiv-VIM charge-Gram H-matrix | HDiv linear/nonlinear charge-Gram validation |
| `RadHACApKPointKernel` | Gauss-point Laplace kernel building block | charge-Gauss operator path |
| `RadHACApKHDivSystemTet` | HDiv-VIM tetrahedral system | HDiv tet/system validation |
| `RadHACApKPEECManager` | PEEC circuit extraction | `radia.peec_hacapk_solver` smoke/solver validation |
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

### Nonlinear Iteration: Anderson+Picard Default (Hantila Historical)

**POLICY (2026-07-03, Sugahara "Anderson+Picardで行きます")**: the production
nonlinear iteration for the moment path (ALL soft iron -- tet/wedge/pyramid/hex
face-charge elements) is **Picard + safeguarded Anderson(1)**, default ON
(`moment_anderson_depth = 1`, gated on ANY moment solve: method 0 LU, method 1
dense-K, method 2 H-matrix; hysteresis uses the B-input moment Picard).  The
safeguard accepts the accelerated iterate only when it reduces the residual, so
linear and well-behaved solves are unaffected.  `rad.SolverConfig(
moment_anderson_depth=0)` opts out to plain Picard.

**Why**: plain Picard DIVERGES on strongly-coupled hysteresis blocks at the
descending-branch steep-slope steps (measured 2026-07-03, 4x4x4 hex block;
`relax_param=0.3` does NOT rescue); Anderson(1) completes the full loop at
~4.5 iters/step (method 2) / ~10 (method 0), methods agreeing to ~6e-5.
Golden: `validation_test/hysteresis/test_binput_moment.py::
test_E_coupled_block_loop_default_anderson` (also locks the default).

**Hantila polarization method (historical)**: Hantila (1975) splits
`B = mu_0*(1+alpha)*H + mu_0*R` (R = M - alpha*H), giving a constant LHS.
This remains useful as background for hysteresis solver design, but the
production soft-iron route is HDiv-VIM plus its current nonlinear iteration and
preconditioner stack. Do not add a separate public Hantila demag backend.

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
`radia.sparsesolv_ngsolve`. An earlier CLAUDE.md draft described these
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

Permanent magnets are specified by direct magnetization (`ObjHexahedron(verts, [Mx, My, Mz])`) or `MatPM(Br, Hc, axis)`. (The `MatMagFixed` / `MatMagLinear` / `MatMagCurve` skeleton trio was removed 2026-06-26 -- it duplicated direct-M / MatPM and never implemented real demagnetization; full PM demagnetization is planned via `MatPM`.)

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

### Data Persistence Policy: Always Save Data to Committed JSON

**POLICY**: 計算で得たデータは**常に `.json` に保存せよ**。とくに
**図・表・公表結果の裏付けとなるデータ**は、`.json` として保存し、その
`.json` を**版管理された永続的な場所**（図・スクリプトと同じディレクトリ、
commit 済み）に置くこと。そのデータを **`C:/temp/...` や `%TEMP%`、未コミット
の sweep ディレクトリにのみ置いてはならない**（transient で消える）。

**Why**: commit 済みの図の元データが `C:/temp` にしか無いと、temp が
クリアされた瞬間に**図が再生成不能**になる。実例 (2026-05, IGTE ESIM
digest): `sweep_heatmap.png` は commit されていたが、その
`sweep_results.json` は `C:/temp/igte_bench/` にしか無く消失 — ~2 時間の
32 ケース sweep を再走しない限り図を再生成できなくなった。

**Rules**:
- commit 済みの図 (`.png`/`.pdf`) は、その元データ `.json` を**同じ
  ディレクトリに commit** すること。「図は git、データは temp」は禁止。
- プロット/解析スクリプトの**入力**データパスは、既定で commit 済みの場所
  （スクリプト隣）を指す。`--out-dir` 等で scratch を指すのは使い捨て実行の
  時だけで、正式な run は `.json` をリポジトリ内に書く。
- これは下記 **Benchmark Policy** と上記 **File Placement Policy** の拡張：
  `.json` は「`.py`／図の隣に置くデータ」であり、`.png` と同じ扱い。
- 対象: sweep (`sweep_*.py`)、benchmark (`bench_*.py`)、出力がプロット/
  公表される全ての `calc_*.py`。

### Benchmark Policy

**POLICY (MUST: 重い検証・計算時間の測定は mdx または hibino で、2026-07-14 Sugahara)**:
**壁時計 / タイミング / スケーラビリティ計測と solver-heavy な `validation_test/` の実行は、
アイドル状態の `mdx` または `hibino` で行うことが MUST**。`mdx`を既定の静音計算ホストとし、
MATLAB、大規模メモリ、長時間ジョブ、またはmdx混雑時は`hibino`を使う。LAB での timing は codex の
並列 build / pytest / 他計算に汚染され無意味なので、**論文・docs・意思決定に用いる時間データは LAB で
測ってはならない**。LAB で許されるのは correctness / smoke(数値一致・収束確認)のみ。mdx が塞がって
いれば `hibino`へ振り分けるか、両方が塞がっていれば timing は延期する(下記)。この MUST は benchmark script (`bench_*.py`) だけでなく、
ad-hoc な timing 計測・scaling sweep・build-time 測定すべてに適用される。

**POLICY (mdx/hibino = 静音計算ホスト、他ジョブ終了後にのみ走らせる、2026-07-14)**: `mdx` と
`hibino` は研究室の **計算用・静音マシン**。壁時計 / タイミング計測および重い計算ジョブは、
**他のプロセス（別の計算ジョブ・build・CI・pytest・他ユーザ / codex の計算）が終わってから**＝
選択したホストが **アイドルのときだけ** 開始する。実行中の別ジョブと **並走させない** — 並走は
計算ホストの
"静音で再現可能" という唯一の価値を壊し、かつ他人のジョブを汚染する。
- **開始前に必ず選択したホストの稼働状況を確認**する：`ssh mdx pwsh` または `ssh hibino pwsh` で
  `Get-Process python` の本数 / CPU負荷 / build・CIの有無を見る。重いジョブが走っていれば **待つ**
  （横入りしない）。
- `mdx` と `hibino` がともに塞がっているとき：**正しさ照合 (correctness / smoke — 数値が一致するか・
  収束するか) は LAB で可**（LAB は codex 競合下でも一致確認は問題ない）。**タイミング計測は
  いずれかがアイドルになるまで延期**する。LAB のタイミングは codex の並列 build / pytest に汚染されて
  無意味（crash / hang / SIGKILL / noise）なので信用しない。
- 重い `validation_test/` は実行したホスト名、開始時刻、実行時間、メモリ条件を結果JSONまたはログに残す。
- これは codex↔claude の **共有ポリシー**（AGENTS.md の Benchmark Policy にも同文を置く）。
  背景と過去インシデントは memory `benchmark_on_mdx_quiet_machine.md`。

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

### File Management: .jou and .vol Only

**POLICY**: `ensure_jou_path()` (Layer 2, PySide6 `radia_export_menu.py`) saves `.jou` before any export or evaluation. `.jou` basename determines all output filenames. `.cub5` is NOT saved by the pipeline.

| File | Role |
|------|------|
| `.jou` | **Single source of truth**. Text, diff-able, version-controllable, reproducible across machines and Cubit versions |
| `.vol` | **Computation interface**. Sole interface between Cubit and NGSolve. No ABI dependency |

**Design — .jou loaded or saved before proceeding**:
- Every Export Mesh operation calls `ensure_jou_path()` first
- `ensure_jou_path()` resolves in 3 steps:
  1. `.jou` already loaded (via `play` or `get_current_journal_file`) → use it
  2. `.jou` saved earlier in this session (`s_lastJouPath`) → use it
  3. Neither → prompt user to save `.jou` now (QFileDialog) → `save journal`
- **No operation proceeds without a known `.jou` path**
- All output filenames derive from `.jou` basename: `{base}.vol`, `{base}.msh`, `{base}_J.sol`, etc.
- `.cub5` is NOT saved by the pipeline. `.vol` is the sole computation interface

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

### Cubit GUI: PySide6-Only -- No Qt5 / PyQt5 (radia 4.80.0+)

**POLICY**: The Cubit panel GUI is **PySide6 (Qt6) only**.  Do NOT use Qt5,
PyQt5, or the old C++ Qt5 `.ccl` Claro component anywhere.

- **Target**: Coreform Cubit **2025.12**, which bundles **PySide6** in its
  embedded Python 3.10.  It does NOT ship PyQt5, and cannot load the old Qt5
  `.ccl` (missing Qt5 DLLs) -- which is why the GUI is now PySide6.
- **Cubit GUI = Python PySide6**, not C++ Qt:
  - Layer 2 (in-Cubit, Python 3.10): `src/radia/panels/register_toolbar.py` (Solve
    menu) + `src/radia/panels/radia_export_menu.py` (Export Mesh menu + dialogs).
    `from PySide6.QtWidgets import ...`; note `QAction`
    is in `PySide6.QtGui` (was `QtWidgets` in Qt5).
  - Layer 3 (MATLAB/Simulink): application blocks in the single Radia library,
    backed by UI-neutral `*_design.py` and headless `calc_*.py` contracts.
    Normal Radia Python does **not** depend on PySide6.
- **No fallback** (per "No Fallbacks -- Fail Fast"): never
  `try PySide6 except PyQt5`.  An old Cubit without PySide6 must raise the
  ImportError loudly so the operator fixes the environment.
- **C++ plugin is Qt-free**: `cubit_mesh_export.ccm` (APREPRO `export`
  commands) + `cubit_mesh_curver.pyd` link only the Cubit C++ API +
  statically-linked netgen -- NO Qt (verified: 0 Qt5 imports).  No Qt SDK
  is needed to build the plugin.
- **Removed in radia 4.80.0**: the Qt5 `.ccl` GUI (`RadiaComp.cpp`,
  `HighOrderMesh.cpp`) and the `cubit_mesh_export_ccl` build target.  A leftover
  `<Cubit>/bin/cubit_mesh_export.ccl` is a stale artifact (Cubit 2025.12 ignores
  it; remove on redeploy).

### Simulink Block Layout Policy (2026-07-20)

**POLICY**: Current Radia analysis interfaces are masked Simulink blocks, not
desktop PySide6 windows or Jupyter workbenches. Settings remain in UI-neutral
`DesignSpec` contracts; a block mask exposes the application configuration,
explicit execution semantics, timeout/backend choice, and artifact location
without hiding headless flags in callback-only state. Ports use fixed types and
stable dimensions, and the block icon remains compact and readable in a normal
engineering model.

No Radia application has a notebook-workbench exception. Cubit may use PySide6
for its toolbar/menu, but normal Radia applications must not add new
`ModePanel`, `AnalysisWindow`, or notebook-workbench surfaces.

### Result Output Policy: ne / DoF / time + analysis integral quantities (2026-05-29)

**POLICY**: Every application block's result (and the `--output` JSON) MUST surface
the analysis's key reporting quantities: **element count, DoF, and a
breakdown of compute time** (every `t_*_s`), plus the **important
integral quantities of that analysis**.  For induction-heating (radia-ih)
this includes **heat (P_wp, the workpiece power dissipation)** and, for
the thermal step, **temperature reported as mean (volume-averaged
`∫T dV / ∫dV`), max, and min** -- not a single peak value.

**How it is enforced**:
- Simulink result adapters summarize ne / DoF /
  all `t_*_s` timings / heat `P_*` / temperature mean-max-min from the JSON
  result artifact, keyed
  on the ACTUAL emitted names (`wp_ndof`/`ndof`,
  `t_bem_solve_s`/`t_solve_s`, `P_wp_W`, `T_mean_C`/`T_max_C`/`T_min_C`).
  Do NOT regress to a fixed-spelling per-solver cascade -- it silently
  showed nothing when a calc script used a different key.
- `calc_*.py` emit the integral quantities in the result dict;
  `calc_common.calc_main` writes the full dict to the `--output` JSON.
  A new analysis MUST add its element count, DoF, timing breakdown and
  the physically meaningful integral quantities (total power, total
  energy, mean/max/min of the primary field) to its result dict.

### Bug-Pattern Catalog Policy (2026-05-31)

**POLICY**: The lab keeps a **learned bug-pattern catalog** in
`packages/radia-mcp/src/radia_mcp/meta/bug_patterns.py`, exposed
through MCP tools `bug_patterns_lookup(...)` + `bug_patterns_stats()`
on the `mcp-server-radia-meta` server.  Every bug class that bites
in a real incident gets one entry with: `id` / `title` / `topics` /
`severity` / `first_seen` / `last_seen` / `what` / `root_cause` /
`detection` / `prevention` / `related`.

**Workflow**:

1. **BEFORE writing new code** in an affected area, call
   `bug_patterns_lookup(topic="<area>")` (panel / release / cubit /
   cubit-license / build / ngsolve etc.).  Read every entry's
   `prevention` field; those are the rules to follow.
2. **WHEN a new bug class fires**, add a new entry to
   `PATTERNS` in `bug_patterns.py`.  Entry must point at the test /
   audit / skill that catches the regression, so the pattern isn't
   just a memory aid -- it's anchored to enforced infrastructure.
3. **WHEN an existing bug class re-fires**, bump `last_seen` to today's
   date.  Repeated occurrences are a signal the prevention step isn't
   strong enough -- consider tightening the static gate.

**Why this exists** (2026-05-31): the session-after-session repeat of
the same bug classes (TaskManager late-import UnboundLocalError, .log
truncated by super-then-append, init.py vs pyproject.toml version
mismatch, phantom block/sideset/nodeset, Cubit 2025.8+ logout-only-
local) showed that memory entries + skill docs are not enough on
their own.  Putting the catalog behind an MCP tool means Claude
encounters it as a first-class capability the moment it picks
`mcp-server-radia-meta` -- a query is a natural step in the
diagnosis flow, not a hidden discipline.

### New Application Block Contract Policy (2026-07-20)

**POLICY**: A new human-facing analysis application **MUST** follow the
Simulink-block recipe in
[`docs/panels/ADDING_NEW_PANEL.md`](docs/panels/ADDING_NEW_PANEL.md):

| File | Role |
|---|---|
| `src/radia/panels/calc_<topic>.py` | headless CLI: `build_argparser()` + `run(args) -> dict` + `calc_main(run, parser)` |
| `src/radia/<topic>_design.py` | UI-neutral `DesignSpec` dataclass; maps settings to calc CLI argv |
| `matlab/+radia/+simulink/` | block builder/runtime adapter over the shared application contract |
| `matlab/radia_simulink_library.slx` | masked application block in the single Radia library |
| `tests/panels/test_<topic>_golden.py` or `validation_test/panels/` | golden-band/API lock on the canonical sample |
| `tests/matlab/test_simulink_workflow.m` | library, mask, port, failure, and execution contract |

**Why this matters**: the recipe makes argparse/`DesignSpec` the single source
of truth for Python, MCP, MATLAB, and Simulink:

- the block cannot silently invent a solver flag outside `DesignSpec`
- the numerical implementation remains in the tested API/CLI, not a mask callback
- the explicit runner writes one versioned result/log artifact contract
- a future MEX backend can replace the Python backend without changing ports

**Enforcement**: `python tools/audit_new_panel_contract.py` is the
static gate for the calc/DesignSpec side, `tests/test_simulink_application.py`
locks the Python-backed runner, and `tests/matlab/test_simulink_workflow.m`
locks the library. Do not add desktop PySide6 or notebook panels.

**Output convention**: GMSH `.msh` for visualization + the text `.log`
file from Persistence Policy.  Do NOT add new bespoke output formats
unless the analysis genuinely cannot fit those two channels.  The
result `dict` returned by `run(args)` is what the application runner records;
keep the keys aligned with the Result Output Policy table (2026-05-29).

### Result Output Persistence Policy (2026-05-30)

**POLICY**: Every application Run produces a **persistent artifact pair**.
Headless CLI output may stay beside the input `.vol`/`.step`; a Simulink block
run stores the pair in its configured run root:

| File | Content | Producer |
|---|---|---|
| `<base><suffix>.json` | Structured result dict (per [Result Output Policy 2026-05-29](#result-output-policy-ne--dof--time--analysis-integral-quantities-2026-05-29)) | `calc_*.py --output` |
| `<base><suffix>.log` or `<run>/run.log` | Verbatim command/stdout/stderr/failure record | Simulink application runner |

- **Naming**: `<base>` from input path basename (`json_output()` /
  `msh_output()` convention), `<suffix>` from the panel-mode suffix
  (e.g. `_peec_bem`, `_fem_kelvin`, `_omega_reduced`, `_peec_ind`).
  The `.log` path is derived as `os.path.splitext(json_path)[0] + ".log"`.
- **Overwrite**: `.log` is **overwritten** on each Run (user-confirmed
  2026-05-30). Older artifacts do NOT accumulate -- one Run = one
  artifact pair. To preserve a history, copy the `.log` out before
  re-running.
- **Write timing**: At Run end, regardless of exit code -- a failed Run
  also leaves a `.log` capturing the failure tail. This is required
  for triage of incident reports ("panel crashed" / "wrong number")
  without asking the user to copy-paste from the Output box.
- **Scope**: All application blocks uniformly.
  The Simulink contract writes `command.txt`, `solver_result.json`, `run.log`,
  and versioned `result.json` in one run directory. An application promoted to
  MEX/ROM must preserve equivalent provenance and diagnostics.

**How it is enforced**:
- Simulink blocks build commands from `DesignSpec.build_command()` through
  `radia.simulink.application` and persist `result.json` / `run.log` artifacts.
- Python and MATLAB tests lock result creation, run-log capture, application
  identity, mask wiring, and failure behavior.

**Why**:
- 1 Run = 1 artifact pair beside the geometry → audit / reproducibility
  / paper figures all flow from this convention without per-panel
  custom code.
- A user incident report can be reduced to "attach the `.log`" and
  the diagnostic is complete (no copy-paste, no truncation, no
  encoding loss).
- This complements `C:/radia_panel_log.txt` (Cubit session-wide
  Claude debug log) which is per-session rather than per-Run; the
  `.log` is the per-Run audit artifact.

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
│  ensure_jou_path(): .jou save -> basename for all output files  │
│  export netgen/gmsh/jmag_nastran/vtk (APREPRO commands)  │
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

**Filename convention**: `ensure_jou_path()` (Layer 2, PySide6 `radia_export_menu.py`) saves `.jou` first. All output files derive basename from `.jou`: `{base}.vol`, `{base}.msh`, `{base}_J.sol`, `{base}_q.sol`, etc.

### Layer Isolation Rules

| Rule | Rationale |
|------|-----------|
| Layer 4 must NOT `import cubit` | Cubit is expensive commercial software. Computation must work without it. |
| Layer 4 must NOT import PySide6/PyQt5 | Headless computation. JSON stdout only. |
| Layer 3 must NOT load Cubit's Python runtime | Simulink is a separate MATLAB process; it uses files and the validated headless runner. |
| Layer 2 must NOT `import radia` or `import ngsolve` | DLL conflicts (Cubit bundles its own numpy/scipy). |
| Layer 1 (C++) has no Python/Qt dependency | `.ccm` links the Cubit C++ API directly. The Qt5 `.ccl` was removed in radia 4.80.0. |

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
| `src/radia/panels/calc_mesh_eval.py` | 4 (no GUI) | p-convergence + format QA called from docs/notebooks |

### Cubit Plugin: C++ First, No Python ABI Dependency

**POLICY**: Cubit plugin functionality MUST be implemented in C++ to avoid Python ABI mismatch. Cubit embeds Python 3.10; NGSolve/Radia use Python 3.12. Sharing Python objects between them causes segfaults and DLL conflicts.

- `.ccm`: Link Cubit C++ API (cubiti, cubit_util) directly -- no Python, no Qt (the Qt5 `.ccl` GUI was removed in radia 4.80.0; the GUI is the PySide6 toolbar in Layer 2)
- `cubit_mesh_curver.pyd`: pybind11 for Python 3.12 -- does NOT link Cubit C++ libraries
- Netgen `SetNCD2Names()` is not exposed to Python -- call from C++ side in `NetgenCurverPure`
- Interface between Cubit and NGSolve: **.vol file** (text format, no ABI dependency)
- Export Mesh BACKEND (`export` APREPRO commands) is C++ only (see Cubit Mesh Export Module below). The Export Mesh GUI menu IS Python/PySide6 (Layer 2, `radia_export_menu.py`) -- that is the supported GUI since radia 4.80.0, not a forbidden one.

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
  をしている** ことが原因 (実例: 2026-05-02 Claude が `timeout 4`
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

**POLICY**: The Export Mesh **backend** is **C++ only** -- all mesh extraction / curving / file writing lives in the `cubit_mesh_export.ccm` APREPRO commands (`export ...`).  The Export Mesh **GUI** is the PySide6 toolbar (`src/radia/panels/radia_export_menu.py`, Layer 2), which only collects options and calls the C++ `export` command via `cubit.cmd`.  Do NOT re-implement export logic in Python, and do NOT add a second GUI.  (The Qt5 `.ccl` GUI was removed in radia 4.80.0.)

### C++ Plugin Architecture

| Component | File | Purpose |
|-----------|------|---------|
| `.ccm` (plugins/) | `cubit_mesh_export.ccm` | APREPRO commands: `export gmsh/jmag_nastran/vtk/netgen` |
| GUI (Layer 2) | `src/radia/panels/radia_export_menu.py` | PySide6 Export Mesh menu + dialogs (replaced the Qt5 `.ccl`, removed in radia 4.80.0) |
| `.pyd` (plugins/) | `cubit_mesh_curver.pyd` | pybind11: Cubit-free mesh curving |

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
| `src/radia/panels/radia_export_menu.py` | PySide6 Export Mesh GUI (replaced RadiaComp.cpp `.ccl`, removed in radia 4.80.0) |
| `src/cubit_plugin/RadiaPlugin.cpp` | Command plugin registration (.ccm) |
| `src/cubit_plugin/cubit_mesh_export_pybind.cpp` | pybind11 module (.pyd) |

### Build

Compact Netgen is statically linked — no nglib.dll/ngcore.dll needed for the Cubit plugin.

### Testing

```bash
python tests/cubit/test_export_combinations.py   # All format x option combinations
python tests/cubit/test_ho_volume_all_formats.py  # Order=2 volume accuracy (sphere)
```

---

## Agent Division of Labor: Claude Owns Push + CI; Codex Owns Release (2026-08-17)

**POLICY**: Split of work between AI agents on this repo.

- **Claude carries its own work to green CI on `main`.** Claude implements,
  tests locally, commits (this-session files BY NAME, per the existing commit
  hygiene), runs `python tools/ci_preflight.py`, pushes to `main`, and then
  watches GitHub Actions to green with `python tools/check_ci.py --watch`. A red
  CI on Claude's own commit is Claude's to fix-forward, not to hand over.
- **Release is codex's job.** codex cuts tags, runs the PyPI publish, drives the
  `release-qud` four-machine verification, and deploys (100号機 / mdx / hibino).

Claude does **NOT**: push tags, invoke the `release-qud` flow, publish to PyPI,
or deploy to remote machines. If asked to "release", Claude gets the work
committed, pushed and CI-green, then hands off to codex for the tag and publish.

**Why the line moved (2026-08-17, Sugahara)**: an author who stops at the local
commit cannot tell whether the work actually builds anywhere but their own
machine, and a handoff at that point makes the person who fixes the breakage
someone who did not write it. Push and CI are part of finishing a change.
Release is genuinely different work — version numbers, four machines, PyPI
propagation — and stays where it was.

**What this does not license**: pushing to `main` while a red preflight is
outstanding, force-pushing, or pushing another agent's uncommitted work. The
codex↔Claude source-edit mutex still applies: when the shared tree has
in-flight work from another agent, do the change in a worktree and push the
branch rather than reaching into `main`.

---

**Last Updated**: 2026-04-02
**For**: Claude Code AI Assistant
**Project**: Radia Magnetic Field Computation
