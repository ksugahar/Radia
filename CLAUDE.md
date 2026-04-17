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
    panels/               # Cubit GUI panels
    *.py                  # Python modules
  src/core/               # C++ source
  src/ext/
    HACApK_LH-Cimplm/    # H-matrix library (MIT)
    sparsesolv/           # Compact AMS/COCR (source, build is separate PyPI)
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

**POLICY**: All source code, documentation (`docs/**/*.md`), comments, commit messages, and docstrings in the Radia repository MUST be written in **English**. Japanese text is NOT allowed in tracked files. Exception: `CLAUDE.md` may contain Japanese policy descriptions. Conversation with the user may be in Japanese, but repository content must remain English-only.

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

### Field Comparison: Vector Difference

**POLICY**: Compare magnetic fields using **vector difference** `norm(B1 - B2)`, not scalar magnitude difference `abs(|B1| - |B2|)`. Magnetic field is a vector quantity.

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

### CI Testing Policy

**POLICY**: CI/CD のテストだけでは不十分。Cubit が必要な機能（`export_curved`, panels, BEM extractor）は **Cubit 環境でのローカルテストが必須**。CI は C++ ビルドと基本テスト（Cubit 不要なもの）のみ。

**リリース前の必須テスト**:
1. CI: C++ ビルド + pytest（Cubit 不要テスト）
2. ローカル: Cubit + system Python で `export_curved` テスト（球、トーラス）
3. ローカル: Cubit パネルの動作確認

CI が通っても Cubit テストに通らなければリリースしない。

### Cubit Batch Self-Testing Policy

**POLICY**: Claude は Cubit を **完全ヘッドレス** (`-batch -nographics -nojournal`) で起動し、自力で機能試験を走らせること。GUI や人間の操作は不要。

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

**特記**: FEMEEM エクスポートの出力パスは **40 文字以下** にすること。`inpin.f90::chkinib(filename*40)` が長い Python `tempfile.TemporaryDirectory()` パスを切り詰めて `forrtl severe (29)` を起こす。`C:\tmp\<short>\` 等を使う。

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

## Compact HX Preconditioner (ngsolve.la)

Compact AMS/AMG/COCR types are in `ngsolve.la`. Source: `src/ext/sparsesolv/` (monorepo integrated).
Import: `from ngsolve.la import CompactAMSPreconditioner, COCRSolver`

### Policy: Compact HX for HCurl Problems

**POLICY**: Use **Compact HX** (Compact Hiptmair-Xu) as the default preconditioner for HCurl curl-curl + mass systems. Compact HX is a HYPRE-free, TaskManager-native AMS implementation available via `ngsolve.la`.

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
**For**: Claude Code AI Assistant
**Project**: Radia Magnetic Field Computation
