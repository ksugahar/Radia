"""
Accelerator electromagnet analysis knowledge base for the Radia MCP server.

Covers: complete pipeline for accelerator magnet design and analysis using
Radia (CoilBuilder, HDiv-VIM, IMA) + Cubit (hex mesh) + NGSolve (FEM, Kelvin).

Sources:
  - CLAUDE.md "Accelerator Magnet Solver Architecture"
  - CLAUDE.md "Hantila Polarization Method"
  - src/radia/panels/samples/em/c_type_dipole/ (rescued panel-only STEP fixtures)
  - src/radia/panels/samples/em/ (legacy packaged EM panel samples)
  - src/radia/coil_builder.py
  - src/radia/coil_geometry.py
"""


# ============================================================
# Topic: Overview
# ============================================================

# Authoritative topic enum for the dispatcher tool (wired into
# `<short>_topics()` via common.register_topics_tool).
TOPICS: dict[str, str] = {
    "overview":             "Accelerator electromagnet analysis pipeline architecture",
    "coilbuilder":          "CoilBuilder API: add_straight/add_arc + close optimization, "
                             "mirror/rotate_copies for symmetry, to_radia for Biot-Savart Hs",
    "kelvin_workflow":      "Kelvin transform open-boundary workflow: Cubit hex mesh + "
                             "ACIS callbackgeometry curving + NGSolve Omega-reduced",
    "hantila":              "Hantila polarization method: B = mu_0(1+alpha)H + mu_0*R, "
                             "LU factored ONCE, back-substitution iteration",
    "hysteresis":           "Energy-based / Play B-input hysteresis: convex U_k, "
                             "reversible + irreversible separation, fast Picard inverse",
    "ima":                  "Image Method of Analysis: parallel-to-mirror (+), "
                             "perp-to-mirror (-), boundary element on-plane limitation",
    "harmonics":            "Magnetic field harmonic analysis on circular bore "
                             "(dipole/quad/sext quality, multipole expansion)",
    "clebsch_hodograph":    "radia-em Clebsch hodograph mode and accelerator "
                             "pole-face inverse-design references",
    "symmetry_reductions":  "1/2, 1/4, 1/8 sector reductions: identify pairing + "
                             "anti-pairing planes, periodic master/slave BC setup",
    "1_8":                  "(alias) 1/8 sector reduction",
    "eighth":               "(alias) 1/8 sector reduction",
    "c_yoke_symmetry":      "(alias) C-yoke symmetry reduction",
    "kelvin_benchmark_vs_em": "(alias) Kelvin benchmark vs accelerator EM",
    "all":                  "Concatenated dump of all topics above",
}

OVERVIEW = """
# Accelerator Electromagnet Analysis Pipeline

## Architecture

The complete pipeline for accelerator electromagnet analysis:

```
CoilBuilder (Radia)
  add_straight() / add_arc() -> continuous path (beam optics style)
  close() -> multi-variable optimization for loop closure
  mirror() / rotate_copies() -> symmetry (dipole, quadrupole)
  to_radia() -> Biot-Savart source Hs (NO coil mesh needed)
  write_step() -> GMSH visualization
       |
       | Hs (analytical, Biot-Savart)
       v
Cubit (hex mesh)
  Iron yoke + air + Kelvin domain -> hex sweep mesh
  export netgen "model.vol" order N -> arbitrary-order hex mesh
  CallbackGeometry -> ACIS direct curving (NO STEP/OCC)
       |
       | curved hex mesh (.vol)
       v
NGSolve FEM
  Omega-reduced scalar potential (2-scalar, fastest formulation)
  Kelvin transformation (open boundary, no PML)
  Energy-based B-input Play model (nonlinear hysteresis)
    - Reversible/irreversible separation -> convex energy
    - LU factored ONCE, back-substitution iteration
    - Fast inverse: Picard 2-3 iter (vs Newton ~100 iter)
  GmshPostExport -> GMSH visualization (+ coil STEP overlay)
```

## Unique Capabilities (no other software provides all of these)

- Coil mesh NOT needed (Biot-Savart analytical source)
- STEP/OCC NOT needed for meshing (ACIS CallbackGeometry)
- Hex mesh with arbitrary-order curving (order 1-5)
- Energy-based hysteresis with fast inverse
- Open source: `pip install radia`

## Coordinate Convention (beam optics)

The coordinate system follows accelerator convention:
- x = horizontal (beam orbit plane)
- y = vertical (field direction for dipole)
- z = beam direction (longitudinal)

Cubit journal files use the SAME coordinate system. No rotation needed.

## Electromagnet Simulink block modes

`radia_simulink_library/Applications/Electromagnet` / `EMDesignSpec` exposes
five Method choices. The non-IH notebook workbench is retired:

| Method | Layer-4 script | Purpose |
|--------|----------------|---------|
| Omega | `calc_accel_magnet.py --formulation omega` | Reduced scalar-potential FEM with Kelvin open boundary |
| A-Phi | `calc_accel_magnet.py --formulation a` | Vector-potential FEM path |
| HDiv-VIM | `calc_accel_hdiv.py` | Mesh-backed HDiv-VIM path for iron-yoke magnets |
| Kelvin Benchmark | `calc_kelvin_benchmark.py` | Analytical sphere-in-uniform-field check for the Kelvin pipeline |
| Clebsch hodograph | `calc_clebsch_hodograph.py` | Self-meshed reduced-potential demonstration; forward field quality and Clebsch potentials |

The Clebsch hodograph mode is the block-facing, verified forward demo.  The
broader accelerator pole-face inverse-design examples and saturation/end-pack
studies live under `docs/clebsch_hodograph/`,
`docs/clebsch_hodograph/demos/`, and the `mcp-server-accelerator`
`accelerator("two_plane_design")` / `accelerator("endpack_cobake")` topics.

## Typical Parameters

| Magnet Type | NI (A-turns) | Gap (mm) | mu_r | B_gap (T) |
|-------------|-------------|----------|------|-----------|
| Dipole      | 2,000-20,000 | 20-80  | 1000-5000 | 0.2-2.0 |
| Quadrupole  | 1,000-10,000 | 10-40  | 1000  | gradient 5-20 T/m |
| Sextupole   | 500-5,000   | 20-60  | 1000  | ... |

## Required Packages

```bash
pip install radia[cubit]     # Radia + Cubit plugin + NGSolve + MKL
cubit-plugin-install         # Deploy Cubit plugin + panels
```
"""


# ============================================================
# Topic: CoilBuilder
# ============================================================
COILBUILDER = """
# CoilBuilder: Fluent Interface for Coil Geometry

## Overview

CoilBuilder provides a beam-optics-inspired path builder for multi-segment coils.
Each segment's end state (position, orientation) automatically becomes the next
segment's start state, enabling continuous coil path definition without manual
coordinate tracking.

## Basic Usage

```python
from coil_builder import CoilBuilder

mm = 1e-3
coil = (CoilBuilder(current=2000)
    .set_start([0, 0, 0])
    .set_cross_section(width=35*mm, height=105*mm)
    .add_straight(62.5*mm)
    .add_arc(radius=22.5*mm, arc_angle=90)
    .add_straight(50*mm)
    .add_arc(radius=22.5*mm, arc_angle=90)
    .add_straight(62.5*mm)
    .add_arc(radius=22.5*mm, arc_angle=90)
    .add_straight(50*mm)
    .add_arc(radius=22.5*mm, arc_angle=90))
```

## Key Methods

| Method | Description |
|--------|-------------|
| `set_start([x, y, z])` | Set initial position (meters) |
| `set_cross_section(width, height)` | Set cross-section (meters) |
| `add_straight(length)` | Add straight segment along current direction |
| `add_arc(radius, arc_angle, tilt=0)` | Add circular arc (degrees) |
| `close()` | Multi-variable optimization to close loop |
| `mirror(plane)` | Mirror coil across plane ('xy', 'xz', 'yz') |
| `rotate_copies(n, axis)` | Create n rotated copies |

## Outputs

| Method | Returns | Use Case |
|--------|---------|----------|
| `to_radia()` | List of Radia objects | Biot-Savart field source |
| `to_occ()` | OCC shape | STEP export, visualization |
| `write_step(path)` | None (writes file) | GMSH overlay |
| `to_wire_segments()` | (segments, current) | Wire model for panels |
| `combined_occ([other_coils])` | Fused OCC shape | Multi-coil STEP |

See also: `docs/complex_coil_geometry/complex_coil.ipynb` -- 8-segment beam-steering
coil showcase using CoilBuilder add_straight/add_arc with a Biot-Savart field map.

## Dipole Example (C-type magnet)

```python
from coil_builder import CoilBuilder
import radia as rad

mm = 1e-3
NI = 2000

# Upper racetrack coil around main leg
upper = (CoilBuilder(current=NI)
    .set_start([47.5*mm, 100*mm, 0])
    .set_cross_section(width=35*mm, height=105*mm)
    .add_straight(62.5*mm)
    .add_arc(radius=22.5*mm, arc_angle=90)
    .add_straight(50*mm)
    .add_arc(radius=22.5*mm, arc_angle=90)
    .add_straight(62.5*mm)
    .add_arc(radius=22.5*mm, arc_angle=90)
    .add_straight(50*mm)
    .add_arc(radius=22.5*mm, arc_angle=90))

# Lower coil: true geometric mirror across the XZ plane (same current;
# the mirrored circulation reverses as a pseudovector, m -> -M @ m)
lower = upper.mirror('xz')

# Radia field source
rad.UtiDelAll()
objs = upper.to_radia() + lower.to_radia()
all_coils = rad.ObjCnt(objs)

# Source field CoefficientFunction for NGSolve
Hs = rad.RadiaField(all_coils, 'h')

# STEP export for GMSH overlay
upper.combined_occ([lower]).WriteStep("dipole_coils.step")
```

## Cubit Panel Integration

For the Accelerator Magnet panel, create a coil script with `build_coil()`:

```python
# coil_dipole.py - loaded by the panel
from coil_builder import CoilBuilder

def build_coil():
    mm = 1e-3
    coil = (CoilBuilder(current=2000)
        .set_start([47.5*mm, 100*mm, 0])
        .set_cross_section(width=35*mm, height=105*mm)
        .add_straight(62.5*mm)
        .add_arc(radius=22.5*mm, arc_angle=90)
        .add_straight(50*mm)
        .add_arc(radius=22.5*mm, arc_angle=90)
        .add_straight(62.5*mm)
        .add_arc(radius=22.5*mm, arc_angle=90)
        .add_straight(50*mm)
        .add_arc(radius=22.5*mm, arc_angle=90))
    return coil
```

## CRITICAL Notes

- All coordinates in meters (Radia unit policy)
- Coil is NOT meshed -- it is a Biot-Savart field source only
- `mirror('xy')` reverses current direction automatically
- `arc_angle` is in degrees, positive = left turn in current frame
- `tilt` parameter rotates the bending plane (default 0 = horizontal bend)
"""


# ============================================================
# Topic: Kelvin Workflow
# ============================================================
KELVIN_WORKFLOW = """
# Cubit Hex Mesh -> NGSolve Curved -> Kelvin Transform Pipeline

## Overview

For accelerator magnets, the workflow uses:
1. Cubit for structured hex mesh (yoke + air + Kelvin exterior domain)
2. NGSolve for FEM solve with Kelvin transformation (open boundary)
3. Radia CoilBuilder for analytical coil source field

## Kelvin Transformation (Open Boundary)

The Kelvin transformation maps infinite exterior to a bounded shell:

```
Physical domain: r in [0, a]     (yoke + air)
Kelvin domain:   r' in [a, R]    maps to r = a^2/r' in [a, inf)
```

Key property: (a/r')^2 weight factor for ALL material constants.
This preserves the Laplace equation's conformal symmetry.

## Cubit Journal Structure

```python
# 1. Create yoke geometry
brick x {W} y {D} z {H}
# ... boolean operations for C-shape ...

# 2. Create air region (sphere or box)
sphere radius {a}
subtract volume {yoke_id} from volume {air_id}

# 3. Create Kelvin exterior domain
sphere radius {R}     # R > a
subtract volume {air_id} from volume {kelvin_id}

# 4. Mesh
volume all scheme tetmesh   # or scheme sweep for hex
mesh volume all

# 5. Labels
block 1 add volume {yoke_id}
block 1 name "yoke"
block 2 add volume {air_id}
block 2 name "air"
block 3 add volume {kelvin_id}
block 3 name "kelvin"

# Periodic BC on Kelvin inner/outer surfaces
sideset 1 add surface {kelvin_inner}
sideset 1 name "kelvin_int"
sideset 2 add surface {kelvin_outer}
sideset 2 name "kelvin_ext"

# 6. Export
export netgen "magnet.vol" order 2 overwrite
```

## NGSolve FEM Solve (Omega-Reduced Scalar Potential)

```python
from ngsolve import *
import radia as rad

mesh = Mesh("magnet.vol")

# Kelvin weight CF
a_kelvin = 0.15  # Kelvin radius [m]
kc = [0, 0, 0]   # Kelvin center
rp_sq = (x-kc[0])**2 + (y-kc[1])**2 + (z-kc[2])**2 + 1e-20
kelvin_fac = a_kelvin**2 / rp_sq

kelvin_weight = mesh.MaterialCF({
    "kelvin": kelvin_fac,
}, default=CF(1.0))

# Material properties
MU_0 = 4e-7 * 3.141592653589793
mu_cf = mesh.MaterialCF({
    "yoke": CF(MU_0 * mu_r),
}, default=CF(MU_0))

coeff = mu_cf * kelvin_weight

# Coil source field (Biot-Savart, no mesh needed)
rad.UtiDelAll()
coil_objs = coil_builder.to_radia()
Hs = rad.RadiaField(rad.ObjCnt(coil_objs), 'h')

# Omega-reduced formulation
fes = Periodic(H1(mesh, order=2, dirichlet="kelvin_ext"))
u, v = fes.TnT()

a = BilinearForm(fes)
a += coeff * grad(u) * grad(v) * dx(bonus_intorder=4)
a.Assemble()

f = LinearForm(fes)
f += coeff * Hs * grad(v) * dx(bonus_intorder=4)
f.Assemble()

phi = GridFunction(fes)
phi.vec.data = a.mat.Inverse(fes.FreeDofs(), "pardiso") * f.vec

# Total field
H_total = Hs - grad(phi)
B_total = mu_cf * H_total
```

## Kelvin Weight Rule

**CRITICAL**: The Kelvin weight (a/r')^2 applies to ALL material constants:
- mu in the Kelvin domain
- Source terms (Hs) in the Kelvin domain
- Both bilinear and linear forms

Reference: Sugahara, IEEE Trans. Magn. 2022.

## Visualization

```python
from radia.gmsh_post_export import GmshPostExport

post = GmshPostExport(mesh)
post.add_field("|B|", CF(B_total.Norm()), ncomp=1)
post.write("magnet_field.msh")

# Overlay with coil STEP in GMSH:
# gmsh magnet_field.msh dipole_coils.step
```

## For More Detail

Use the `kelvin_transformation` tool from mcp-server-ngsolve for the full
Kelvin transformation mathematics and implementation details.
"""


# ============================================================
# Topic: Hantila
# ============================================================
HANTILA = """
# Hantila Polarization Method

## Overview

Hantila (1975) splits the constitutive relation into a constant linear part
plus a nonlinear residual:

```
B = mu_0*(1+alpha)*H + mu_0*R    where R = M - alpha*H
```

For HDiv-VIM, the charge-Gram demag operator maps M -> H_demag (constant, geometry-only):

```
H = H_ext + N*M
Substituting M = alpha*H + R:
(I - alpha*N)*H = H_ext + N*R    <- constant LHS, LU factored ONCE
```

## Key Advantage: LU Factored Once

The system matrix `(I - alpha*N)` is constant because:
- N depends only on geometry (fixed)
- alpha is a chosen constant (fixed)
- Only R changes each iteration (material nonlinearity)

This means:
1. Factor `(I - alpha*N)` with LU decomposition ONCE
2. Each nonlinear iteration = cheap back-substitution only

## Comparison with Standard Methods

| Feature | Picard (rad.Solve) | Newton | Hantila |
|---------|-------------------|--------|---------|
| Matrix factorization | Every iteration | Every iteration | **Once** |
| Jacobian needed | No | Yes (dM/dH) | **No** |
| BH curves | Yes | Yes | **Yes** |
| Hysteresis | No | No | **Yes** |
| Cost per iteration | O(N^3) LU | O(N^3) LU | **O(N^2) back-sub** |

## Usage

```python
from radia.hantila_solver import solve_hantila

# BH curve (nonlinear, no hysteresis)
result = solve_hantila(
    iron_container, source=coil,
    bh_data=BH_DATA,
    alpha=500.0,    # linear part slope (tune for convergence)
    tol=1e-4
)

# Hysteresis (per-element material handles)
result = solve_hantila(
    iron_container, source=coil,
    mat_handles=handles,
    alpha=500.0,
    relax=0.5       # under-relaxation for stability
)

# Result
M = result['M']       # (n_elem, 3) in A/m
H = result['H']       # (n_elem, 3) in A/m
iters = result['iterations']
B = rad.Fld(iron_container, 'b', [0, 0, 0.05])
```

## Alpha Selection

- alpha should approximate the average differential susceptibility
- Too small: slow convergence
- Too large: oscillation (use relax parameter)
- Typical: alpha = 500-2000 for silicon steel

## Current Limitation

The production soft-iron route is HDiv-VIM.  Keep new accelerator
magnet work on mesh-backed HDiv spaces rather than resurrecting the
old moment/surface-charge soft-iron paths.

## Reference

F.I. Hantila, Rev. Roum. Sci. Techn. - Electrotechn. et Energ., 1975.
"""


# ============================================================
# Topic: Hysteresis
# ============================================================
HYSTERESIS = """
# B-input Play and Energy Hysteresis Models

## Overview

Radia provides two B-input hysteresis models for quasi-static electromagnet
simulation. The Play model is recommended (faster, no sign constraints).

## Play Model (Recommended)

```python
from radia.hysteresis_io import load_hys_file

K, eta, f_k_tables = load_hys_file('material.hys')
mat = rad.MatPlayHysteresis(K, eta, f_k_tables)
# K: number of play operators
# eta: ndarray[K], play thresholds in Tesla
# f_k_tables: list of (r_array, f_array) tuples (shape functions)
```

## Energy Model

```python
mat = rad.MatEnergyHysteresis(K, eta, f_k_tables, eps=1e-6)
# Same parameters + eps convergence tolerance
# Requires non-negative, monotonically increasing shape functions
```

## Comparison

| Feature | Play Model | Energy Model |
|---------|-----------|--------------|
| Forward (B->H) | O(K) direct | Newton (~100 iter) |
| Inverse (H->B) | Newton + analytical Jacobian | K independent Newton |
| Shape functions | No sign constraint | Must be non-negative |
| Speed | 4-9 us/eval | 100-500 us/eval |

## State Management

Both models support state save/restore for quasi-static stepping:

```python
rad.MatApl(iron, mat)

# After each converged time step:
rad.MatHysCommitState(mat)           # Commit state for next step

# Save/restore for branching or checkpointing:
state = rad.MatHysSaveState(mat)     # Save (ndarray, length K*9)
rad.MatHysRestoreState(mat, state)   # Restore
```

## Integration with Hantila Solver

The Hantila polarization method naturally supports hysteresis:
- Standard Picard/Newton: cannot handle path-dependent materials
- Hantila: R = M - alpha*H captures the nonlinear residual including
  hysteresis history, updated each iteration with the current state

## Quasi-Static Stepping Example

```python
for step, H_ext in enumerate(excitation_sequence):
    # Solve with current hysteresis state
    result = solve_hantila(iron, source=coil_at(H_ext),
                           mat_handles=handles, alpha=500)

    # Commit converged state
    for h in handles:
        rad.MatHysCommitState(h)

    # Extract B at gap center
    B_gap = rad.Fld(iron, 'b', [0, 0, 0])
    print(f"Step {step}: B_gap = {B_gap[2]:.4f} T")
```
"""


# ============================================================
# Topic: IMA
# ============================================================
IMA = """
# IMA (Image Method of Analysis) for Accelerator Magnets

## Overview

IMA exploits geometric symmetry to reduce problem size. For a quarter model
of a C-type dipole, IMA replaces 3/4 of the elements with image contributions.

## Sign Selection Rule

| Field vs Mirror Plane | IMA Sign |
|----------------------|----------|
| Field **parallel** to mirror | **+** (symmetric) |
| Field **perpendicular** to mirror | **-** (antisymmetric) |

## Examples

```python
# Dipole: Bz field, X-Z quarter model
# Bz is parallel to X-mirror, perpendicular to Z-mirror
rad.Solve(container, 0.0001, 100, 0, image='+x-z')

# Quadrupole: Bx field, X-Z quarter model
# Bx is perpendicular to X-mirror, parallel to Z-mirror
rad.Solve(container, 0.0001, 100, 0, image='-x+z')
```

## How to Determine Signs

1. Identify the dominant field component at the observation point
2. For each symmetry plane:
   - If the field component is PARALLEL to the plane -> **+**
   - If the field component is PERPENDICULAR to the plane -> **-**

## Quarter Model Example (C-type dipole)

```
Full model:          Quarter model (x>0, z>0):
  +---------+         +----+
  |  yoke   |         |yoke|
  |  [gap]  |    ->   |[  ]|  + IMA '+x-z'
  |  yoke   |         +----+
  +---------+
```

## Boundary Element Limitation

IMA produces incorrect results (~0.5x magnitude) for **boundary elements**
(faces ON the symmetry plane) when observation points are also on the plane.

**When IMA is Safe**:
1. Elements offset from symmetry planes (not touching)
2. Observation points off-plane
3. Mesh-backed HDiv-VIM with explicit mirror materialization after solve

**Workaround**: For models with boundary elements on symmetry planes,
use explicit element duplication instead of IMA.

## IMA + Nonlinear Solve

IMA works with nonlinear solvers (Picard, Newton, Hantila):

```python
rad.Solve(container, 0.0001, 1000, 1, image='+x-z')
# method=1 (BiCGSTAB) + IMA for quarter model
```

Verified: IMA sign is correct for both linear and nonlinear problems.
"""


# ============================================================
# Topic: Symmetry Reductions (1/2, 1/4, 1/8)
# ============================================================
SYMMETRY_REDUCTIONS = """
# C-type Magnet Symmetry Reductions: 1/2, 1/4, 1/8

## Two distinct Kelvin panel paths (do NOT conflate)

The Radia EM stack ships **two unrelated panel modes** that both use
Kelvin transformations.  Requests for "1/8 support" almost always
mean the EM panel path, not the benchmark path.

| Panel mode | Sample location | Purpose | Reductions |
|------------|----------------|---------|------------|
| **EM panel FEM/HDiv-VIM** | `panels/samples/em/em_1-{1,2,4,8}.jou` | Production C-yoke dipole / quad analysis | 1/1, 1/2, 1/4, **1/8** |
| **EM panel "Kelvin Benchmark"** | `panels/samples/kelvin_benchmark_sphere_1_{2,4}.vol` | Verify Kelvin pipeline against analytical sphere-in-uniform-Hz | 1/2, 1/4 only |

The Kelvin Benchmark mode does **NOT** ship a 1/8 reduction because
the magnetic-sphere-in-uniform-Hz BVP is fundamentally not
1/8-symmetric: the source `H0 z_hat` reverses sign under the z=0
mirror.  This is a property of the BVP, not a Cubit / NGSolve issue,
and no formulation tweak can recover it.  See
`mcp-server-radia-ngsolve` topic `kelvin_transformation` ->
"Why 1/8 is unsupported for the sphere benchmark".

## C-yoke 1/8 reduction (shipped, golden-tested)

The C-type dipole geometry has full octant symmetry (yoke and coil
both invariant under x->-x, y->-y, z->-z modulo current direction),
so the 1/8 reduction is geometric and stable.  Sample:

    panels/samples/em/em_1-8_eighth.jou
    tests/panels/golden/em_eighth_mu1000.json

Validated 2026-04-26: 774/774 Kelvin pairs identified at machine
precision (post deterministic-anchor fix), 56,369 elements, 11,708
DOFs at p=1, B_origin_mag = 4.892e-04 T at I=2000 A (regression
guard, NOT physics accuracy -- the reduced coil BCs do not match
the full racetrack).

## ELF CEFC 2020 BC convention (C-yoke)

The 1/8 C-yoke uses these symmetry boundary conditions:

| Plane | BC | Meaning |
|-------|-----|---------|
| x=0 | `ht=0_x` | tangential H zero (mirror, B_x continuous) |
| y=0 | `ht=0_y` | tangential H zero (mirror, B_y continuous) |
| z=0 | `bn=0_z` | normal B zero (anti-mirror) |

In the launcher config:

    "kelvin_reduction": {"x": "ht=0", "y": "ht=0", "z": "bn=0"}

The Kelvin offset goes along **x** (the first axis listed), and the
Cubit auto-Kelvin script webcuts at x=offset_dist (`kelvin_far`,
Dirichlet anchor), y=0 (`sym_ht=0_y`), z=0 (`sym_bn=0_z`).  GND
nodeset sits at the Kelvin centre `(offset_dist, 0, 0)`.

This is opposite in sign pattern to the sphere benchmark (which
uses `bn=0` on both x and y for the lower-order reductions, then
would need `ht=0_z` on z for 1/8 -- which fails as discussed
above).  **The two BVPs need different BC patterns; copying the
sphere convention onto the C-yoke (or vice-versa) gives wrong
results.**

## Reduction ladder for the EM panel

| Reduction | jou | sym BCs | Kelvin path |
|-----------|-----|---------|-------------|
| 1/1 (full) | `em_1-1_full.jou` | none | full Kelvin shell or octant offset |
| 1/2 (split y) | `em_1-2_half_z.jou` | `ht=0_y` | Kelvin offset along x |
| 1/4 (split x,y) | `em_1-4_quarter_xz.jou` | `ht=0_x, ht=0_y` | Kelvin offset along x |
| 1/8 (split x,y,z) | `em_1-8_eighth.jou` | `ht=0_x, ht=0_y, bn=0_z` | Kelvin offset along x |

All four are golden-tested in `tests/panels/test_em_*_golden.py`.

## Don't add a 1/8 sphere benchmark

If a future task asks to "add 1/8 to the Kelvin benchmark panel",
the answer is no -- the BVP doesn't admit it.  See
`memory/feedback_kelvin_1_8_not_a_real_use_case.md` for the
multi-hour debug trail (rho_min sweep + h-refinement both
disproved numerical-issue hypotheses) before the BVP-symmetry
realisation landed.

The build script `panels/samples/kelvin_benchmark_sphere_1_8_build.py`
exists for parameterised-core completeness (`--frac 1_2|1_4|1_8`)
but the resulting `.vol` is intentionally NOT in the wheel and
NOT in `panels/samples/`.
"""


# ============================================================
# Topic: Harmonics
# ============================================================
HARMONICS = """
# Field Harmonics / Multipole Analysis

## Overview

Accelerator magnets are characterized by their field harmonics (multipole
content). The field in the aperture is expanded in a Fourier series:

```
B_y + i*B_x = B_ref * sum_n (b_n + i*a_n) * (x + iy)^(n-1) / R_ref^(n-1)
```

Where:
- B_ref: reference field (typically B at R_ref)
- b_n: normal multipole coefficients (units: 10^-4 relative)
- a_n: skew multipole coefficients (units: 10^-4 relative)
- R_ref: reference radius (typically 2/3 of aperture radius)
- n=1: dipole, n=2: quadrupole, n=3: sextupole, ...

## Extracting Harmonics from Radia/NGSolve Field

```python
import numpy as np
import radia as rad

# Sample field on a circle at z=0 (mid-plane)
R_ref = 0.010  # 10 mm reference radius
N_pts = 64     # sampling points
theta = np.linspace(0, 2*np.pi, N_pts, endpoint=False)
pts = np.column_stack([R_ref * np.cos(theta),
                       R_ref * np.sin(theta),
                       np.zeros(N_pts)])

# Get field at sample points
B = np.array([rad.Fld(container, 'b', p.tolist()) for p in pts])
# or use batch: B = rad.Fld(container, 'b', pts.tolist())

# Complex field: By + i*Bx
B_complex = B[:, 1] + 1j * B[:, 0]

# FFT to get harmonics
coeffs = np.fft.fft(B_complex) / N_pts
B_ref = abs(coeffs[1])  # dipole component magnitude

# Normal and skew coefficients in units of 10^-4
for n in range(1, 11):
    b_n = coeffs[n].real / B_ref * 1e4
    a_n = coeffs[n].imag / B_ref * 1e4
    if n == 1:
        print(f"n={n} (dipole):     b_{n} = {b_n:.2f}, a_{n} = {a_n:.2f}")
    else:
        print(f"n={n}:              b_{n} = {b_n:.2f}, a_{n} = {a_n:.2f}")
```

## From NGSolve GridFunction

```python
from ngsolve import *

# Assuming B_total is a CoefficientFunction from FEM solve
B_vals = np.array([B_total(mesh(*p)) for p in pts])
B_complex = B_vals[:, 1] + 1j * B_vals[:, 0]
coeffs = np.fft.fft(B_complex) / N_pts
# ... same harmonic extraction as above
```

## Quality Criteria (Typical)

| Harmonic | Dipole spec | Quadrupole spec |
|----------|-------------|-----------------|
| b_1 (dipole) | 10000 (reference) | < 1 |
| b_2 (quad) | < 1 | 10000 (reference) |
| b_3 (sextupole) | < 5 | < 5 |
| b_5 (decapole) | < 1 | < 1 |
| a_n (skew) | < 1 | < 1 |

Allowed harmonics depend on magnet symmetry:
- Dipole (2-fold): b_1, b_3, b_5, b_7, ... (odd normal)
- Quadrupole (4-fold): b_2, b_6, b_10, ... (n = 2 + 4k)

## Important Notes

- Sample at z=0 (magnet center) for 2D harmonics
- For 3D field quality, sample at multiple z positions
- R_ref should be 2/3 of the good field region radius
- Use enough sample points: N_pts >= 2 * max_harmonic_order + 1
- B_ref convention: real part of n=1 coefficient for dipole
"""


# ============================================================
# Topic: Clebsch hodograph panel mode
# ============================================================
CLEBSCH_HODOGRAPH = """
# Clebsch Hodograph Mode in `radia-em`

`radia-em` includes a **Clebsch hodograph** Method backed by
`src/radia/panels/calc_clebsch_hodograph.py`.  It is intentionally a
verified forward-analysis panel mode, not the entire inverse pole-design
research program:

- The mode self-meshes its demonstration geometry, so it does not require a
  user-supplied `.vol` file.
- It solves the reduced scalar potential `Omega`, recovers the field, and
  post-processes conjugate Clebsch potentials (`psi`, `Phi`) as an
  a-posteriori field-quality diagnostic.
- It emits JSON plus a GMSH `.msh` visualization through the same Layer-3 /
  Layer-4 panel contract as the other `radia-em` methods.

For pole-face inverse design, end-pack construction, FFAG sectors, and
saturation studies, use the accelerator knowledge server:

```python
accelerator("two_plane_design")
accelerator("endpack_two_plane")
accelerator("endpack_cobake")
accelerator("sector_saturation")
```

## Field-plane (Chaplygin) design direction: VERIFIED end-to-end (2026-07-23)

The field-plane hodograph `(q, theta)` linearizes a saturable `mu(q)` exactly
(the material law becomes a known coefficient of the independent variable).
Three facts make it practical, all locked by
`validation_test/clebsch_legendre/`:

- **Legendre potential** `chi = H.r - Psi`: one scalar unknown whose FIRST
  derivatives are the physical coordinates (`e_H.r = chi_q`,
  `e_perp.r = chi_theta/q`, `Psi = q chi_q - chi`). Self-adjoint equation
  `d/dq(mu q chi_q) + ((mu q)'/q) chi_thth = 0`.
- **B-radial form** feeds a measured secant curve `mu_s(B)` directly
  (`mu_d = mu_s^2/(mu_s - B dmu_s/dB)`), no table inversion.
- **End-to-end check** (`verify_chaplygin_bend_design.py`): a 90-deg
  saturable flux-guide bend designed by ONE linear solve, then re-solved by
  an independent nonlinear FEM on the designed shape. Wall |B| matches the
  spec to 0.29-0.43 % mean (5-85 deg core), MMF to 0.12 % (0.03 % flat
  case), constant-mu sanity reproduces an exact annulus to 1.8e-9,
  mesh-converged, J single-signed (no folding). Golden-banded; committed
  sidecars `results_*.json`.

Pitfall: the forward nonlinear FEM check REQUIRES damped Picard
(omega ~ 0.35). The undamped iteration stalls and fakes a ~140 % "design
error" -- see bug pattern `reference-secant-picard-oscillation`
(`bug_patterns_lookup(topic="nonlinear")`).

Scope: the DESIGN direction on a known hodograph domain (spec posed in
`(|B|, theta)`: saturable throats, flux regulators, bend guides). The
ANALYSIS direction (fixed pole shape, hodograph image unknown) remains open
research; do not present it as solved.

Repository anchors:

- `docs/clebsch_hodograph/public_demo.ipynb` (result-saved public demo entry point)
- `docs/clebsch_hodograph/examples_catalog.ipynb` (full source/result catalog)
- `docs/clebsch_hodograph/DESIGN_METHODOLOGY.md`
- `docs/clebsch_hodograph/HODOGRAPH_BACKBONE.md`
- `docs/clebsch_hodograph/HDIV_VIM_CLEBSCH_BRIDGE.md`
- `docs/clebsch_hodograph/demos/` (runnable companion set)
- `tests/feec/test_clebsch_hodograph_research.py`
- `validation_test/clebsch_legendre/verify_chaplygin_bend_design.py`
  (field-plane design -> shape -> independent nonlinear FEM, golden-banded)
- `validation_test/clebsch_legendre/verify_chi_modesum_solver.py`
  (Legendre-potential mode-sum machinery vs closed forms)
"""


# ============================================================
# Dispatcher
# ============================================================
_TOPICS = {
    "overview": OVERVIEW,
    "coilbuilder": COILBUILDER,
    "kelvin_workflow": KELVIN_WORKFLOW,
    "hantila": HANTILA,
    "hysteresis": HYSTERESIS,
    "ima": IMA,
    "harmonics": HARMONICS,
    "clebsch_hodograph": CLEBSCH_HODOGRAPH,
    "clebsch": CLEBSCH_HODOGRAPH,
    "hodograph": CLEBSCH_HODOGRAPH,
    "symmetry_reductions": SYMMETRY_REDUCTIONS,
    "1_8": SYMMETRY_REDUCTIONS,  # alias
    "eighth": SYMMETRY_REDUCTIONS,  # alias
    "c_yoke_symmetry": SYMMETRY_REDUCTIONS,  # alias
    "kelvin_benchmark_vs_em": SYMMETRY_REDUCTIONS,  # alias
}


def get_electromagnet_documentation(topic: str = "overview") -> str:
    """Return electromagnet documentation for the given topic.

    Args:
        topic: One of: overview, coilbuilder, kelvin_workflow, hantila,
               hysteresis, ima, harmonics, clebsch_hodograph,
               symmetry_reductions, all

    Returns:
        Documentation string for the requested topic.
    """
    topic = topic.strip().lower()

    if topic == "all":
        parts = []
        for name, content in _TOPICS.items():
            parts.append(f"{'=' * 60}\n# Topic: {name}\n{'=' * 60}")
            parts.append(content)
        return "\n\n".join(parts)

    if topic in _TOPICS:
        return _TOPICS[topic]

    available = ", ".join(sorted(_TOPICS.keys()))
    return (
        f"Unknown topic: '{topic}'. "
        f"Available topics: {available}, all"
    )
