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
    panels/               # Cubit GUI panels
    *.py                  # Python modules
  src/core/               # C++ source
  src/ext/
    HACApK_LH-Cimplm/    # H-matrix library (MIT)
    sparsesolv/           # Compact AMS/COCR (built into radia wheel, exposed as radia.sparsesolv_ngsolve)
  packages/
    cubit-mesh-export/    # Independent PyPI package (pip install cubit-mesh-export)
      src/cubit_mesh_export/
        check.py          # check-vol CLI + check_consistency() API
        radia_cubit_mesh.pyd  # C++ pybind11 module (bundled)
  tests/                  # Radia tests + tests/mcp/
  examples/
  docs/
  Build.ps1               # MSVC + MKL build
  install_full.py          # One-command full setup
```

**PyPI packages** (2 independent packages in same monorepo):

| Package | Install | Purpose |
|---------|---------|---------|
| **radia** | `pip install radia` | C++ core + Python (MMM/MSC/PEEC, panels, MCP) |
| **cubit-mesh-export** | `pip install cubit-mesh-export` | High-order curved mesh export from Cubit (does NOT require radia) |

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

## Critical Policies

### Self-Driving Loop Discipline (2026-06-24)

Drive an autonomous / self-paced loop by **task completion, not a clock**: finish
one verified step, then immediately start the next — do not insert idle interval
ticks. Every iteration must produce **concrete, verified progress** (a number that
came from running the code, never an estimate), recorded to the project's internal
validation notes. **Never overclaim** — a "pass" must trace to a checked value. Stop
and ask only where the decision is genuinely the user's.

### Publish Boundary: No Validation Provenance in Public Artifacts (2026-06-24)

Public artifacts (this repo, PyPI packages, public docs) **lead with analytic
solutions**. Do NOT put into them: internal absolute paths or local working-directory
names; a third-party tool's **benchmark numbers** used as a validation basis; or
"verified / validated against <external tool>" attributions, including the names of
that tool's source files. Citing a *published* convention or a peer-reviewed paper is
fine — attributing validation to an internal or third-party reference is not. Keep
cross-validation provenance in machine-local notes only; stored regression-reference
values may remain, but **unattributed**.

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

### File Placement Policy

**POLICY**: Generated output files (`.png`, `.msh`, `.vtu`, `.vol`) must be placed **next to their corresponding `.py` script**.
- Example outputs belong in `examples/<category>/` alongside their script
- Do NOT place generated files at the repository root
- `.msh` files in `examples/**/gmsh_models/` are tracked (pre-generated mesh definitions)
- Build output goes to `build*/` or `dist/` (both gitignored)

### Sample Promotion Ladder: tests → examples → panels (2026-05-02)

**POLICY**: Every sample lives in exactly ONE of three tiers.  The tiers
have distinct purposes and a strict promotion ladder:

| Tier | Purpose (intent) | Audience | Ships in wheel? |
|------|------------------|----------|-----------------|
| `tests/**/fixtures/` | **実装の基本機能の確認** — golden test fixture, machine-readable, minimal. | CI / Codex / developer | No |
| `examples/<topic>/` | **研究的側面も含む例題の提供** — research-oriented demonstration, exploratory geometry, README-backed. | researchers, contributors | No (gitignored outputs OK) |
| `src/radia/panels/samples/` | **工学的実問題を動く形で提供** — engineering real-problem solution, runs end-to-end through the panel. | end users (Cubit panel Browse dialog) | Yes (package-data) |

**Promotion gates** (each is a hard gate; no skipping):

- **tests/ → examples/**: golden test for the geometry passes
  (`tests/panels/test_*_golden.py` JSON inside hard band) AND a
  human-readable `README.md` is added that explains the physics /
  research question / expected ballpark numbers AND the example runs
  standalone (`python <example>.py` without the panel UI).
- **examples/ → panels/samples/**: the example runs end-to-end through
  the **panel UI** (Layer 3, Cubit `play <sample>.jou` → `radia_export`
  → panel Run button) on the actual engineering geometry — not a toy
  proxy.  Sample listed in `pyproject.toml` package-data.  Wheel
  manifest audit (`deploy` skill, L0) clean.

**Why three tiers, not two**: a sample that locks numerical correctness
(tests) is a different artifact from a sample that teaches a research
concept (examples) which is again different from a sample that solves
an engineer's actual problem on the panel (panels).  Conflating them
either ships incomplete work to end users (panels = examples) or
buries production-ready engineering examples in tests/.

**Concrete consequences**:

- A new geometry STARTS at `tests/**/fixtures/` with a golden lock.
  Promote up only after the next-tier gate is met.
- Demoting (panels → examples → tests) is allowed: e.g. an engineering
  sample that turns out to need rework can move back to `examples/`
  while the issue is investigated, with the wheel package-data line
  removed in the same commit.
- Same geometry MAY exist at multiple tiers if the artifacts differ in
  scope (a minimal fixture in tests/, a richer commented version in
  examples/, a panel-friendly .jou in panels/samples/).  Each lives
  separately and is maintained separately.
- A sample at `panels/samples/` MUST also have a passing golden test
  in `tests/panels/test_*_golden.py` — the upper tier inherits the
  correctness gate of the lower.

**Why end users can't tell broken samples from user error**: they try
the sample, it fails, and the panel looks broken.  One broken sample
discredits the whole panel.  The tier discipline above prevents this.

### Panel Design Workflow Policy (2026-04-23)

**POLICY**: Panels are built in **three strict stages**, each gated by
validation of the previous stage.  Do NOT jump straight to PySide.

**Stage 1 — Enumerate the app-specific variables.**
Write down every knob the user of this specific application might want
to change.  This is a list, not code.  Pin the solver-specific variables
(mesh size, frequency, material, source current, ...) and the
**solver-switch variable** itself (e.g. `--impedance-model linear|esim`,
`--solver pardiso|ams`).

**Stage 2 — CLI Python script (`calc_*.py` under `src/radia/panels/`).**
Turn the Stage-1 list into an argparse-driven Python script.
Computation only, no GUI.  JSON on stdout.  The solver switch **must**
also be a CLI flag so the same script can drive any supported backend.

Stage 2 is validated by running the panel mode end-to-end against its
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

**Stage 3 — PySide panel (`radia_ih.py` / `radia_em.py` / ... under
`src/radia/`).**
Wrap the **validated** Stage-2 script with a PySide `AnalysisWindow`
widget.  The panel launches the CLI via `subprocess.Popen` (per the
4-Layer Architecture) and is forbidden from re-implementing any
computation.  Stage 3 ships only after Stage 2 passes its golden test.

Only Stage-3-ready panels go into `panel_registry.json` and the
`radia_*.py` auto-discovery under `src/radia/`.  Stage-2-only panels
live as CLI scripts and wait for the Stage-3 promotion gate.

**Why**:
- Forces the hard thinking about *what is changeable* before any widget
  code is written (Stage 1 is where over-scoping is cheapest to cut).
- The solver switch being a Stage-2 argument means we catch
  solver-specific bugs with the same sample + golden test; the panel
  UI does not hide them.
- Stage-3 promotion requires a passing golden — stops the historical
  failure mode of shipping a panel whose Run button produces a wrong
  number that nobody notices until a user publishes it.

Related:
- "Panel Samples Quality Policy" above — Stage 2's validation relies
  on trustworthy samples.
- "Cubit Panel Architecture" below (§ 4-Layer) — Stage 3 corresponds
  to Layer 3 (PySide6 window), Stage 2 corresponds to Layer 4 (headless
  calc_*.py).

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
2. **`examples/.../README.md`** mark the validated sample as
   **VERIFIED** with the specific check (e.g. "VERIFIED p=2: slaved=8914
   DOFs, ratio=1.0") so future contributors know the sample is golden,
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

**New research topics**: in flight WIP can live in `examples/` or
`docs/research/` (`.gitignored` for LAB-only) until stable. Promotion
into a `radia_mcp.<topic>` subpackage requires: feature committed +
deploy-verified + golden-tested + knowledge stops referencing
unpublished files. There is **no longer** a separate
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
| **Magnetic A_phi (curl-curl)** | **Henrotte** `{1, r^2, z}` (`radia.radia_axifemm`) | The cylindrical curl operator `B_z = (1/r) d(r A_phi)/dr` produces a `1/r` integrand that standard FE Gauss quadrature cannot integrate accurately near the axis.  Henrotte's `s = r^2` substitution gives clean closed-form integration. |
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
import radia.radia_axifemm as ax

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
`radia.radia_axifemm.AxiHenrotteHeat{Stiffness,Mass}BFI` classes
(added in radia 4.31.0) and the `H1Henrotte` BND DiffOp (radia
4.32.0) are kept in the codebase as parity-conscious infrastructure
for research / publication uses (e.g. comparing convergence rates of
Henrotte vs standard H1 on a scalar problem).  They are NOT used by
production heat solvers and are NOT required.

**Reference**: see
[`docs/axifemm/FORMULATION.md`](docs/axifemm/FORMULATION.md)
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

**POLICY**: Do NOT write fallback chains (`try API_A except: try API_B except: try API_C`). Pick the **one** API that works for the project's target environment (Cubit 2025.3, NGSolve 6.2.2603+, Python 3.12) and commit to it. If the chosen API stops working, fix the call site or raise — never bury the breakage under another path.

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

1. **Dipole approximation accuracy is poor for MSC elements**: MSC (surface charge) elements have distributed charge on 4-8 faces. A single dipole m=M*V approximates this poorly at intermediate distances (r ~ 2-5 element sizes). The O((a/r)^2) error is unacceptable for engineering accuracy.

2. **FMM Solve (Method 3) was useless**: Compact geometries (C-type magnets, iron yokes) have 87% near-field pairs. Near-field correction memory equals the full dense matrix, eliminating FMM's O(N log N) advantage. HACApK (H-matrix, Method 2) is 10-100x faster because ACA+ compression works on the same near-field blocks.

3. **FMM field evaluation had no benefit over direct**: For typical Radia models (N < 10,000 elements), direct B_genComp with TaskManager parallelization is fast enough. FMM overhead (tree build, M2L translation) exceeds direct computation time for these sizes.

4. **HACApK covers all large-scale needs**: H-matrix acceleration (ACA+) provides O(N log N) memory and O(N log^2 N) MatVec for the interaction matrix, which is the actual bottleneck.

**Lesson**: FMM is effective for point charges/dipoles in unbounded space (N-body). It is NOT effective for MSC where source distributions are extended (face integrals) and geometries are compact.

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
| **Cubit plugin** (`radia_export gmsh`) | `.msh v4.1` | Mesh export → GMSH viewer |
| **Cubit plugin** (`radia_export netgen`) | `.vol` | NGSolve mesh interchange |
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
  `ExportGmshCommand`.  The `version` keyword on `radia_export gmsh` is
  accepted for back-compat with old `.jou` files but is ignored
  (always emits v4.1 with a warning if `version 2` is passed).

### Mesh Export Consistency Check Policy

**POLICY**: Before exporting .vol, verify mesh correctness by comparing NGSolve integration values against Cubit ACIS CAD values, per label:

| Check | NGSolve | Cubit CAD | Threshold |
|-------|---------|-----------|-----------|
| **Volume** (per material) | `Integrate(CF(1), mesh, definedon=mesh.Materials(mat))` | `cubit.volume(vid).volume()` | > 1% warning |
| **Area** (per boundary) | `Integrate(CF(1), mesh, BND, definedon=mesh.Boundaries(bnd))` | `cubit.surface(sid).area()` | > 1% warning |
| **Length** (per BBND) | `Integrate(CF(1), mesh, BBND, definedon=mesh.BBoundaries(bnd))` | `cubit.curve(cid).length()` | > 1% warning |

All three checks (Volume, Area, Length) passing per label confirms the mesh is geometrically correct.

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

- The Radia in-tree `radia_export netgen` plugin **bypasses the cap**
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

### Terminology: MMM/MSC vs BEM

**POLICY**: Radia's core solvers (MMM, MSC) are **NOT** BEM (Boundary Element Method). Use precise terminology:

| Term | Method | Library | Description |
|------|--------|---------|-------------|
| **MMM** | Magnetic Moment Method | Radia C++ | Volume magnetization M as DOF, dipole interaction |
| **MSC** | Magnetic Surface Charge | Radia C++ | Surface charge sigma as DOF, solid angle integration |
| **BEM** | Boundary Element Method | **ngsolve.bem** | EFIE/MFIE, HDivSurface, Maxwell/Laplace kernels |
| **PEEC** | Partial Element Equivalent Circuit | Radia Python + C++ | Loop-Star, circuit extraction (L,R,C,M) |

**Do NOT** call MMM/MSC "BEM". They are integral equation methods but with different formulations:
- BEM (ngsolve.bem): surface integral equations (EFIE/MFIE) on conductor/dielectric boundaries
- MMM: volume integral equation for magnetization, solved element-by-element
- MSC: surface charge on element faces, solved via solid angle kernel

**When to use which**:
- Permanent magnets, soft iron → **MMM/MSC** (Radia)
- Eddy currents, shielding, impedance extraction → **BEM** (ngsolve.bem) or **PEEC** (Radia)
- High-frequency scattering → **BEM** (ngsolve.bem, Helmholtz kernel)

### Development Strategy: Complement NGSolve

Radia's role is to **complement NGSolve**, not compete with it. Focus on areas where FEM is weak.

```
┌─────────────────────────────────────────────────────────────────┐
│                    Electromagnetic Analysis                      │
├─────────────────────────────────────────────────────────────────┤
│  NGSolve (FEM)              │  Radia (MMM/MSC/PEEC)             │
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
1. MMM - Magnetic Moment Method for permanent magnets and soft iron
2. MSC - Magnetic Surface Charge for hexahedra/tetrahedra
3. Field computation - B, H, A, Phi in unbounded domains
4. NGSolve integration - RadiaField CoefficientFunction

### Solver Methods: MMM and MSC

| Method | Element | DOF | Description |
|--------|---------|-----|-------------|
| **MMM** | Tetrahedra (4 faces) | 3 (Mx, My, Mz) | Magnetic dipole distributions |
| **MSC** | Hexahedra (6 faces) | 6 (sigma/face) | Surface charge solid angle integration |
| **MSC** | Wedges (5 faces) | 5 (sigma/face) | Transition elements |

**Mixed Element Support**: All solvers (LU, BiCGSTAB, HACApK) support mixed hex+wedge+tet meshes. Variable DOF offset arrays: `m_elemDOF`, `m_elemDOFOffset`, `m_totalDOF`.

**BiCGSTAB Block Jacobi**: Automatically switches to block Jacobi preconditioner when diagonal ratio > 10 or min dominance < 0.1 (distorted elements). Uses LAPACK `dgetrf_`/`dgetri_` for block inversion.

**Interaction Matrix Blocks** (mixed elements):
- **3x3** (tet-tet), **5x5** (wedge-wedge), **6x6** (hex-hex)
- **5x6 / 6x5** (wedge-hex cross), **3x6 / 3x5 / 6x3 / 5x3** (tet-hex/wedge cross)
- Implementation: `SetupInteractMatrix_VariableDOF()`, compile flag `RADIA_MSC_SUPPORT`

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


**Key Features**: Inside/outside classification (solid angle method), TaskManager parallelized batch, complex field support (PEEC+MMM AC).

### Field Calculation: Surface Current vs Surface Charge

- **ObjRecMag**: Surface current model (rectangular blocks). 8-corner BufVect formula, efficient and non-cancelling on symmetry axes.
- **ObjHexahedron/ObjTetrahedron**: Surface charge model (general polyhedra). Face-based solid angle integration. A field may be zero on symmetry axes (mathematical cancellation, not a bug).

**rad.Fld() inside materials**: MMM gives dipole approximations inside materials; MSC gives uniform field per element. For validation, compare sigma values or external field points, not internal fields.

### Vector Potential A Field

A field is **implemented** for all element types using face integration (Wilton et al. formula). Formula: `A = (mu_0/4pi) * (M x BufVect)`. Satisfies `B = curl(A)` (verified numerically). Verification script: `examples/ngsolve_integration/verify_curl_A_equals_B/`.

### User-Facing Element APIs

- `rad.ObjRecMag(center, dimensions, magnetization)` -- Rectangular magnets (optimized formulas)
- `rad.ObjHexahedron(vertices, magnetization)` -- Arbitrary hexahedra (8 vertices)
- `rad.ObjTetrahedron(vertices, magnetization)` -- Tetrahedra (4 vertices)
- `rad.ObjWedge(vertices, magnetization)` -- Wedges (6 vertices)
- Mesh import functions (`netgen_mesh_to_radia`) for complex geometries

### EIEM2 Evaluation Point Convention

**POLICY**: The MSC interaction matrix evaluation point for face `i` is:
```cpp
EvalPt = 0.5 * (FaceCenter[i] + ElementCenter)
```
Do NOT change this. This matches ELF's EIEM2 convention exactly.

**MSC Source Files**: `rad_polyhedron.cpp` (element dispatch), `rad_poly_analytical.cpp` (triangle/quad integration), `rad_interaction.cpp` (interaction matrix, `PrecomputeHexaGeometry()`).

See `docs/MSC_QUICK_START.md` for quick start guide.

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
| **Python** | 3.12.10 | System Python for Radia/NGSolve. Cubit panels call via subprocess. |
| **Coreform Cubit** | 2025.3 | Embedded Python 3.10 + PySide6. Cannot import NGSolve/Radia directly. |
| **NGSolve** | 6.2.2603+ | curvedelements Load, hex/prism curving, Periodic BC fix |

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
Radia solvers (BiCGSTAB, HACApK/method 2, HDiv VIM, moment-yano) stay TaskManager-native.

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
- HDiv VIM and moment-yano C++ solve loops
- Any embarrassingly parallel loop

**When non-TaskManager threading is acceptable**:
- MKL internal threading for dense BLAS/LAPACK/PARDISO calls, guarded by
  `SuspendTaskManager` + `MKLThreadGuard` where Radia controls the call
- Legacy code only until it is migrated; do not add new OpenMP regions

### PyPI Release Workflow (Automated via GitHub Actions)

**POLICY**: PyPI publishing is automatic. Push a version tag (`v*`) and CI/CD handles the rest.

**Release Flow**:
1. Bump version in `pyproject.toml` AND `src/radia/__init__.py` (must match)
2. Update `CHANGELOG.md`
3. `git commit` (do NOT push yet)
4. `/deploy` — build wheel, deploy to 100号機 (WinRM) & mdx (SSH)
5. Test on remote machines (Cubit panels, Mesh Evaluation, etc.)
6. If tests pass: `git push origin main`
7. Wait for CI to pass: `gh run list --limit 3`
8. Tag and push (triggers PyPI publish):
   ```bash
   git tag vX.Y.Z
   git push origin vX.Y.Z
   ```
9. Monitor: `gh run list --workflow release.yml --limit 3`

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
- Upgrade route は **100号機 / mdx 側** (PyPI consumer)。 LAB は release 後に
  metadata 同期のため `pip install -e <path> --no-deps --no-cache-dir` で再 editable 化。
- CI/CD 環境 (e.g. `C:\actions-runner\_work\Radia\Radia\...`) は別管理 (NETWORK
  SERVICE 所有)。 LAB の editable pointer がそちらに drift していたら戻す。

**POLICY (2026-05-02 update)**: **2-tier 配布**: LAB は editable (NAS source、developer
loop)、100号機 と mdx は両方 PyPI install (`pip install radia[cubit] radia-mcp` +
`cubit-plugin-install --all-users`). リリース wheel + Cubit plugin の end-to-end
検証点を 100号機 と mdx の 2 マシンで二重化。mdx editable は 2026-05-02 に retire
(`tools/push_pyds_to_mdx.py` は branch-test 用にのみ残置).

**2 ステージ配布モデル (2026-05-02 simplified)**:

| Stage | マシン | install 形態 | 目的 |
|-------|--------|-----------|------|
| 1 | LAB | `pip install -e .` + `pip install -e packages/cubit-mesh-export` + `pip install -e packages/radia-mcp` | 開発者ループ。最速フィードバック (NAS source 直接編集) |
| 2 | 100号機 / mdx | `pip install radia / radia-mcp / cubit-mesh-export` (PyPI) + `cubit-plugin-install --all-users` (regular-file deploy) | PyPI wheel + Cubit plugin の end-to-end 検証点。100号機 = 21 ユーザの本番、mdx = 別マシンでの cross-machine consistency probe (release-triple Phase 9). Stage-2 が両機で通るまで release OK を宣言しない。 |

**変更点 (2026-05-02)**:
- 旧: LAB editable / 100号機 PyPI / mdx editable (3-tier).
- 新: LAB editable / 100号機 + mdx 両方 PyPI (2-tier).
- 理由: mdx editable は (a) gh CLI 不在で `download_binaries.sh` 不可, (b) legacy site-packages shadow の手動削除が必要, (c) `.pyd` を base64-over-ssh で push する `tools/push_pyds_to_mdx.py` が必須, など落とし穴が多い割に PyPI 検証の代替価値が小さかった。100号機 と同じ "PyPI が動くか" の純粋な検証点に統一。

**LAB のみ editable な 4 パッケージ**:
- `radia` (LAB: `S:\Radia\01_GitHub`)
- `cubit-mesh-export` (LAB: `S:\Radia\01_GitHub\packages\cubit-mesh-export`)
- `radia-mcp` (LAB: `S:\Radia\01_GitHub\packages\radia-mcp`)
- `mcp-server-document` (LAB: `S:\mcp-server`) -- LAB-private (PyPI 配布なし)

LAB で `pip install --upgrade <pkg>` を流すと editable が静かに上書きされて壊れるので注意 (2026-04-28 incident)。release 後の LAB 側 metadata 同期は `pip install -e <path> --no-deps --no-cache-dir` で再 editable 化。`pip install --upgrade` は **100号機 / mdx 用** (PyPI から通常通り upgrade).

**100号機 / mdx 全ユーザー PyPI install**: `C:\Program Files\Python312`
の machine-wide site-packages に PyPI install。リリース毎に admin が
`pip install --upgrade radia==X.Y.Z radia-mcp==X.Y.Z cubit-mesh-export==X.Y.Z`
+ `cubit-plugin-install --all-users` を実行。

**100号機 / mdx Cubit plugin (regular file)**:
- `<Cubit>\bin\radia_cubit.ccl` (regular file from PyPI wheel)
- `<Cubit>\bin\plugins\radia_cubit.ccm` (regular file from PyPI wheel)
- `<Cubit>\bin\plugins\radia_cubit_mesh.cp312-win_amd64.pyd` (regular file from PyPI wheel)

LAB の `Build.ps1` 出力は **NAS の `S:\Radia\01_GitHub` に書かれるが、100号機 / mdx の
PyPI install には反映されない**。C++ 変更を 100号機 / mdx で試すには PyPI release を
切るのが正規ルート。緊急時のみ `tools/push_pyds_to_mdx.py` (mdx) や
`pip install --force-reinstall --no-cache-dir //192.168.11.100/work/00_CAE/Radia/01_GitHub`
(100号機 NAS source override) を使う。通常運用は **PyPI release → 100号機 / mdx で
`pip install --upgrade`**。

### CI Testing Policy

**POLICY**: CI/CD のテストだけでは不十分。Cubit が必要な機能（`export_curved`, panels, BEM extractor）は **Cubit 環境でのローカルテストが必須**。CI は C++ ビルドと基本テスト（Cubit 不要なもの）のみ。

**リリース前の必須テスト**:
1. CI: C++ ビルド + pytest（Cubit 不要テスト）
2. ローカル: Cubit + system Python で `export_curved` テスト（球、トーラス）
3. ローカル: Cubit パネルの動作確認

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
cubit.cmd('radia_export femeem "C:\\tmp\\cub" overwrite')
```

**対象**: `radia_export {gmsh,netgen,nastran,vtk,femeem}`、panels の非 GUI ロジック、BEM extractor、`export_curved`。
GUI が絶対必要なもの (panel dialog のレンダリング) のみ例外。

**前提**: `cubit` は Python API import (`S:/Radia/01_GitHub/src/radia/install_panels.py` の `find_cubit_bin()` で自動検出可)。バッチ起動でライセンス消費あり。

**特記**: FEMEEM エクスポートの出力パスは **40 文字以下** にすること。`inpin.f90::chkinib(filename*40)` が長い Python `tempfile.TemporaryDirectory()` パスを切り詰めて `forrtl severe (29)` を起こす。`C:\temp\<short>\` 等を使う。

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
  mmm_core.pyd          # MMM solver
  peec_matrices.pyd     # PEEC matrix assembly
  *.py                  # Python utility modules
  # NO .dll files
```

**Always use `Build.ps1`** for building. Never use manual cmake commands -- the script handles CMake configure + build + `.pyd` copy to `src/radia/`.

---

## Mesh & NGSolve Integration

### NGSolve Version Requirement

**CRITICAL**: Use NGSolve **6.2.2603** or later. Required for curvedelements .vol Load, hex/prism curving, and Periodic BC fix.

Reference: https://forum.ngsolve.org/t/ngsolve-periodic-boundary-condition-regression-bug-report/3805

Official PyPI ngsolve **6.2.2603**+ includes: **MKL**, **PARDISO**, Periodic BC fix,
**curvedelements Save/Load**, **p-version hex/prism curving**.

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
fes = HDiv(mesh, order=2)  # Best accuracy
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
  STEP -> Cubit -> radia_export netgen "model.vol" order N -> Mesh("model.vol") -> Radia

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
Mesh export is C++ only (`radia_export netgen` APREPRO command in the Cubit plugin).

**Install**: `pip install cubit-mesh-export` (or `pip install radia[cubit]`)
**Source**: `packages/cubit-mesh-export/` in the Radia monorepo

**Consistency checking** (does NOT require Cubit):
```bash
check-vol model.vol                         # CLI (installed with cubit-mesh-export)
check-vol model.vol --json model.vol.json   # With companion JSON from export
```
```python
from cubit_mesh_export.check import check_consistency  # API
```

**Module names**:
- `cubit_mesh_export` — canonical Python package (PyPI: cubit-mesh-export)
- `radia_cubit_mesh` — C++ pybind11 module (bundled in cubit_mesh_export, unchanged)
- `check_vol_consistency` — thin backward-compat re-export in `src/radia/panels/` (imports from cubit_mesh_export.check)

Cubit workflow for journal files: define blocks before export, use the Cubit plugin commands (`cubit.cmd('radia_export gmsh/nastran/vtk ...')`). Requires `CUBIT_PLUGIN_DIR` environment variable (set by `cubit-plugin-install`).

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
| `hacapk_leaf` | 10 | Minimum cluster size (elements). 10 for MSC 6DOF hex (~66 DOF/leaf) |
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

**Sign convention applies uniformly to ALL DOF types** (MMM 3DOF, MSC 5/6DOF, mixed). No DOF-type-specific sign conditionals.

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

For Radia MMM, the interaction matrix N maps M -> H_demag (constant, geometry-only):

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

**Current limitation**: MMM (tetrahedra, 3 DOF) only. MSC (hexahedra, 6 DOF) requires sigma-M conversion (future).

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
3. MMM (tetrahedra) -- no limitation

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

### Coupled PEEC + MMM

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

Radia PEEC unique features: direct circuit extraction (L, R, C), native SPICE netlist, Lanczos MOR, MMM coupling.

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
# From examples/: use ../../src/radia
# From tests/: use ../src/radia
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
```

Import from `src/radia` package (not build directories).

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

**POLICY**: `ensure_jou_path()` (C++ .ccl) saves `.jou` before any export or evaluation. `.jou` basename determines all output filenames. `.cub5` is NOT saved by the pipeline.

| File | Role |
|------|------|
| `.jou` | **Single source of truth**. Text, diff-able, version-controllable, reproducible across machines and Cubit versions |
| `.vol` | **Computation interface**. Sole interface between Cubit and NGSolve. No ABI dependency |

**Design — .jou loaded or saved before proceeding**:
- Every Export Mesh / Mesh Evaluation operation calls `ensure_jou_path()` first
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

**Mesh export is C++ only** (`radia_export netgen` APREPRO command):
- Cubit -> `radia_export netgen "model.vol" order N` -> `.vol` with labels + curving

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

**BEM (legacy)**: BEM solver modules are in `examples/induction_heating/bem_reference/`
for research reference. BEM knowledge is in `mcp-server-radia-ngsolve` (ngsbem_inductance topic).

**References**:
- Djordjevic & Notaros, "Double higher order MoM", IEEE TAP 2004 (geometry/basis independence)
- Marussig et al., "Fast Isogeometric BEM based on Independent Field Approximation", arXiv 2014
- Dolz et al., "Bembel: Fast Isogeometric BEM", arXiv 2019

**`.vol` Must Be Self-Contained**:
- Material labels: `SetMaterial()` -> `materials` section
- Boundary labels: `SetBCName()` -> `bcnames` section
- High-order curving: `curvedelements` text section (upstream Netgen master feature)
- No external STEP/geometry file needed for computation

**Cubit Plugin Responsibility**: The `radia_export netgen` C++ command handles all label + curving embedding into `.vol`. Higher maintenance cost is acceptable for complete separation.

---

## Cubit Panel Architecture

### 4-Layer Architecture

**POLICY**: Cubit, Radia-NGSolve, and computation are **3 separate processes**. No layer imports another layer's libraries.

```
┌─────────────────────────────────────────────────────────────────┐
│  Layer 1: C++ Qt5 (.ccl)                                        │
│  ─────────────────────────────────────────────────────────────  │
│  Export Mesh menu (GMSH/Nastran/VTK/Netgen Vol/FEMEEM/MEG)      │
│  Mesh Evaluation (_p1.vol ... _p5.vol + format QA exports)      │
│  ensure_jou_path(): .jou save -> basename for all output files  │
│  radia_export netgen/gmsh/nastran/vtk (APREPRO commands)        │
└─────────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────────┐
│  Layer 2: Cubit GUI Python (Python 3.10 + PySide6/Qt5)          │
│  ─────────────────────────────────────────────────────────────  │
│  register_toolbar.py -> Solve menu management                   │
│    Radia-NGSolve / Generate Coil / Kelvin / Reload / Verify     │
│  import cubit OK (same process). import radia/ngsolve FORBIDDEN │
│  Launches Layer 3 via subprocess.Popen (detached)               │
└─────────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────────┐
│  Layer 3: Radia-NGSolve PySide6 (Python 3.12 + PySide6)        │
│  ─────────────────────────────────────────────────────────────  │
│  radia_ih.py (IHWindow) — standalone PySide6 application        │
│  Separate process from Cubit. import cubit FORBIDDEN.           │
│  Receives .vol path as CLI argument.                            │
│  Launches Layer 4 via subprocess for computation.               │
└─────────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────────┐
│  Layer 4: Computation (Python 3.12, no GUI)                     │
│  ─────────────────────────────────────────────────────────────  │
│  calc_peec.py (PEEC coil) / calc_fem_kelvin.py (FEM workpiece)  │
│  import cubit FORBIDDEN. import PySide6 FORBIDDEN.              │
│  NGSolve + GMSH only. .vol input, JSON stdout output.           │
└─────────────────────────────────────────────────────────────────┘
```

**Interface between layers**: `.vol` file (text format, no ABI dependency)

**Filename convention**: `ensure_jou_path()` (Layer 1, C++) saves `.jou` first. All output files derive basename from `.jou`: `{base}.vol`, `{base}.msh`, `{base}_J.sol`, `{base}_q.sol`, etc.

### Layer Isolation Rules

| Rule | Rationale |
|------|-----------|
| Layer 4 must NOT `import cubit` | Cubit is expensive commercial software. Computation must work without it. |
| Layer 4 must NOT import PySide6/PyQt5 | Headless computation. JSON stdout only. |
| Layer 3 must NOT `import cubit` | Separate process. Cubit embeds Python 3.10; Layer 3 is Python 3.12. |
| Layer 2 must NOT `import radia` or `import ngsolve` | DLL conflicts (Cubit bundles its own numpy/scipy). |
| Layer 1 (C++) has no Python dependency | `.ccm`/`.ccl` link Cubit C++ API directly. |

### Panel Files

| File | Layer | Purpose |
|------|-------|---------|
| `RadiaComp.cpp` (.ccl) | 1 (C++ Qt5) | Export Mesh menu + Mesh Evaluation |
| `panels/register_toolbar.py` | 2 (Cubit Python) | Solve menu + Radia-NGSolve launcher |
| `radia_ih.py` | 3 (PySide6) | IH analysis window (PEEC+FEM / FEM) |
| `panels/calc_peec.py` | 4 (no GUI) | PEEC filament coil inductance |
| `panels/calc_fem_kelvin.py` | 4 (no GUI) | FEM Kelvin + SIBC (IH workpiece) |
| `panels/calc_mesh_eval.py` | 4 (no GUI) | p-convergence + format QA |

### Cubit Plugin: C++ First, No Python ABI Dependency

**POLICY**: Cubit plugin functionality MUST be implemented in C++ to avoid Python ABI mismatch. Cubit embeds Python 3.10; NGSolve/Radia use Python 3.12. Sharing Python objects between them causes segfaults and DLL conflicts.

- `.ccm`/`.ccl`: Link Cubit C++ API (cubiti, cubit_util) directly -- no Python dependency
- `radia_cubit_mesh.pyd`: pybind11 for Python 3.12 -- does NOT link Cubit C++ libraries
- Netgen `SetNCD2Names()` is not exposed to Python -- call from C++ side in `NetgenCurverPure`
- Interface between Cubit and NGSolve: **.vol file** (text format, no ABI dependency)
- Export Mesh is C++ only (see Cubit Mesh Export Module section below). Do NOT add Python export dialogs or panels.

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

---

## Universal Relaxation Network (URN)

All URN examples, data, and scripts in `examples/universal_relaxation_network/`.

**Policy**:
- Synthetic data MUST be clearly marked as synthetic
- Real-world datasets MUST include license and citation info
- All paper results reproducible from scripts in this directory

---

---

## Cubit Mesh Export Module

**POLICY**: Export Mesh is **C++ only**. All mesh export functionality is in the C++ plugin (`radia_cubit.ccm` + `radia_cubit.ccl`). Do NOT add Python export dialogs, panels, or scripts.

### C++ Plugin Architecture

| Component | File | Purpose |
|-----------|------|---------|
| `.ccm` (plugins/) | `radia_cubit.ccm` | APREPRO commands: `radia_export gmsh/nastran/vtk/netgen` |
| `.ccl` (bin/) | `radia_cubit.ccl` | Qt5 GUI: Export Mesh menu + dialog |
| `.pyd` (plugins/) | `radia_cubit_mesh.pyd` | pybind11: Cubit-free mesh curving |

**Export formats** (all in C++, ACIS geometry projection for curving):

| Format | Command | Max Order | Notes |
|--------|---------|-----------|-------|
| Netgen Vol | `radia_export netgen "f.vol" order 3` | 1-5 | Primary format for NGSolve FEM |
| GMSH v4.1 | `radia_export gmsh "f.msh"`           | 1-3 | Lab-wide standard; structured entity blocks |
| Nastran BDF | `radia_export nastran "f.bdf"` | 1-2 | CTETRA/CTETRA(10), nopyramid option |
| VTK | `radia_export vtk "f.vtk"` | 1-2 | Legacy format, cell types 10/24 |

**GMSH order limit**: Order 4-5 is an error (not fallback). NetgenCurver face/volume
interior node extraction is unreliable at p>=4 (linear interpolation fallback causes
negative Jacobians in GMSH). Use `radia_export netgen` for order 4-5.

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

**POLICY**: Mesh export uses `radia_export netgen` C++ command only. Pure Python reference (`cub5_to_vol.py`) is maintained in the netgen fork, not in Radia. Run `test_vol_multi_geometry.py` (10 shapes) after any NetgenCurver change.

### Companion JSON (.vol.json)

`radia_export netgen` writes a companion JSON alongside the .vol file:
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
| `src/cubit_plugin/RadiaComp.cpp` | Qt5 GUI component (.ccl) |
| `src/cubit_plugin/RadiaPlugin.cpp` | Command plugin registration (.ccm) |
| `src/cubit_plugin/radia_cubit_pybind.cpp` | pybind11 module (.pyd) |

### Build

Compact Netgen is statically linked — no nglib.dll/ngcore.dll needed for the Cubit plugin.

### Testing

```bash
python tests/cubit/test_export_combinations.py   # All format x option combinations
python tests/cubit/test_ho_volume_all_formats.py  # Order=2 volume accuracy (sphere)
```

---

**Last Updated**: 2026-04-02
**For**: Codex AI Assistant
**Project**: Radia Magnetic Field Computation
