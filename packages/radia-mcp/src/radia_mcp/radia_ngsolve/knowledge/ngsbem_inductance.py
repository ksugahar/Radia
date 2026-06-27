"""
NGSolve ngsolve.bem boundary element method for inductance extraction.

Covers: LaplaceSL BEM operators, HDivSurface basis, self/mutual inductance,
Cubit -> SetGeomInfo -> Curve -> BEM pipeline, and practical examples.

ngsolve.bem is NGSolve's native boundary element module. It provides:
- Laplace/Helmholtz single/double layer potentials
- Natural open boundary treatment (no PML, no truncation)
- Low-frequency stable (Laplace kernel for MQS regime)
- Direct integration with NGSolve FE spaces (HDivSurface, SurfaceL2)

Key references:
  - ngsolve.bem documentation: https://docu.ngsolve.org/latest/how_to/ngsbem.html
  - Lucy Weggler's stabilized BEM: https://github.com/Weggler/docu-ngsbem
  - SetGeomInfo for externally imported meshes
"""

NGBEM_OVERVIEW = """
# ngsolve.bem: Boundary Element Method for Inductance Extraction

## What Is ngsolve.bem

NGSolve's native BEM module providing integral equation operators on surfaces.
Key operators:

| Operator | Kernel | Use |
|----------|--------|-----|
| `LaplaceSL` | 1/(4*pi*r) | MQS inductance (L matrix) |
| `LaplaceDL` | d/dn[1/(4*pi*r)] | Scalar DL for scalar BIE + SIBC; **not** the MFIE `n x curl(SL)` operator |
| `HelmholtzSL` | exp(-jkr)/(4*pi*r) | Full-wave BEM |
| `HelmholtzDL` | d/dn[exp(-jkr)/(4*pi*r)] | Full-wave BEM |
| `MaxwellSL` | Full Maxwell kernel | High-frequency |

For inductance extraction at DC to ~1 MHz, **LaplaceSL is sufficient**
(MQS/Darwin regime). No Helmholtz kernel needed.

### CAVEAT: closed-surface LaplaceSL inductance is rank-deficient on dense meshes (verified)

The `L = 1/(e^T L^{-1} e)` total-inductance extraction on a **closed** conductor
surface (e.g. a torus) is **numerically unreliable as the surface mesh is refined**:
the `LaplaceSL` Gram matrix becomes rank-deficient on dense closed-surface meshes
(ngbem surface-integration on closed surfaces), and the inductance error blows up.
Measured on circular loops (Neumann reference `L = mu_0 R (ln(8R/a) - 2)`):
coarse mesh (`curvaturesafety=0.5`, ~89 tris) gave ~+15%, but refined meshes hit
"Rank-deficient matrix (714/715)" and errors of **-9% to -66%**. **Do NOT fix this
by refining the triangulation** -- it makes it worse. Mitigations: keep the mesh
coarse, use **p-refinement** (`order>0`) or **quad elements** (from Cubit), or --
preferred for inductance -- use the **Radia PEEC filament/panel extractor**
(`radia.peec_*`, Neumann-formula based) which does not have this closed-surface
rank-deficiency. (Verified 2026-06-27; rendered in
`docs/bem_extractor/bem_inductance_limitations.ipynb`.
This is a negative/limitation result, kept as knowledge.)

## When to Use ngsolve.bem

| Task | Use ngsolve.bem? | Alternative |
|------|-----------|-------------|
| Self-inductance of a coil | Yes | Neumann formula (analytical) |
| Mutual inductance (multi-conductor) | Yes | Neumann + filament model |
| Inductance with magnetic core | Partial | FEM-BEM coupling |
| Capacitance extraction | Yes | - |
| Permanent magnet field | No | Use Radia (BEM/MSC) |
| Nonlinear iron core | No | Use NGSolve FEM |

## Architecture

```
Cubit / Netgen OCC
  |
  | export mesh (tet/hex volume or surface-only)
  v
NGSolve Mesh (1st order from Cubit, curved by Netgen via mesh.Curve(p))
  |
  | HDivSurface(mesh, order=0)  -- RT0 surface current basis
  v
ngsolve.bem operators (LaplaceSL, use_fmm=False)
  |
  | SL = L_op.mat.ToDense().NumPy()  -- dense matrix extraction
  v
Energy method: L = mu_0 * J^T * SL * J  (for I=1 toroidal current)
  |
  v
Self-inductance [H]
```

### Port Detection

| Cubit Blocks | Mode | Description |
|-------------|------|-------------|
| `source` + `sink` | **1-port** | Gap conductor (current in at source, out at sink) |
| No `sink` block | **loop** | Closed loop (e.g., complete torus) |

## Key Concept: Laplace Single Layer for Inductance

The Laplace single layer potential integral:

```
(LaplaceSL J)(x) = integral_S  1/(4*pi*|x-y|) * J(y) dS_y
```

When J is a surface current on a conductor, the bilinear form
`<LaplaceSL(J), J>` gives the vector potential energy, which is
proportional to inductance:

```
L_matrix = mu_0 * <LaplaceSL(J_trial), J_test>
```

This is the **MQS (Magneto-Quasi-Static)** approximation:
- No wave propagation (Laplace kernel, not Helmholtz)
- Valid for conductor dimensions << wavelength
- Accurate from DC to ~1 MHz for typical PCB/inductor geometries

## Comparison with Radia PEEC

| Feature | ngsolve.bem BEM | Radia PEEC |
|---------|-----------|-----------|
| Kernel | LaplaceSL (continuous) | Filament mutual inductance |
| Mesh | Surface triangles/quads | Segments (1D filaments) |
| Proximity effect | Natural (surface current) | Multi-filament (nwinc/nhinc) |
| Skin effect | Needs SIBC | Built-in (Bessel/Dowell) |
| SPICE output | Manual | Built-in netlist |
| Speed (small) | Slower (dense BEM) | Fast (analytical) |
| Speed (large) | FMM acceleration (use_fmm=True) | H-matrix (HACApK) |
| Accuracy | High (surface integral) | Moderate (filament approx) |
"""

NGBEM_API = """
# ngsolve.bem Python API Reference

## Import

```python
from ngsolve import Mesh, HDivSurface, TaskManager, ds, Integrate, CF, BND
from ngsolve.bem import LaplaceSL  # Core BEM operator
import numpy as np
```

## BEM Operator: LaplaceSL

```python
fes = HDivSurface(mesh, order=0)  # RT0 basis on surface
j_trial = fes.TrialFunction()
j_test = fes.TestFunction()

# CRITICAL: Use .Trace() for boundary integrals
with TaskManager():
    L_op = LaplaceSL(
        j_trial.Trace() * ds("conductor")
    ) * j_test.Trace() * ds("conductor")
```

**Common mistake**: Forgetting `.Trace()` gives silently wrong results.
Always use `j.Trace() * ds(label)`, NOT `j * ds(label)`.

**CRITICAL**: Do NOT use `BilinearForm(...).Assemble()` for BEM operators.
This HANGS indefinitely. Use the pattern above: `LaplaceSL(u*ds) * v*ds`
returns an IntegralOperator whose `.mat` gives the dense matrix directly.

```python
# CORRECT: Direct IntegralOperator (returns immediately)
L_op = LaplaceSL(u.Trace() * ds) * v.Trace() * ds
L_dense = extract_dense(L_op.mat, ndof)

# WRONG: BilinearForm.Assemble() HANGS forever
slp = BilinearForm(LaplaceSL(u.Trace()*ds) * v.Trace() * ds).Assemble()  # NEVER returns
```

## Cubit Surface Mesh for BEM

BEM only needs surface elements (no volume mesh):
```
cubit.cmd("surface all scheme trimesh")
cubit.cmd(f"surface all size {element_size}")
cubit.cmd("mesh surface all")
cubit.cmd("block 1 add tri all")
```

Cubit's automatic curvature refinement makes coarse meshes difficult.
Use `curve all interval N` for direct control.

**DOF guideline** (dense BEM, single thread extract):

| Surface DOFs | Extract Time | Status |
|-------------|-------------|--------|
| 300 | ~5s | Fast |
| 1000 | ~30s | OK |
| 3000 | ~5min | Limit |
| 10000+ | hours | Need H-matrix |

## CRITICAL: Mesh Dimension for BEM

BEM inductance requires **dim=2 surface mesh**, NOT dim=3 volume mesh.

```python
# CORRECT: dim=2 surface mesh (all HDivSurface DOFs are boundary edges)
# Created by: OCCGeometry(surface_shape).GenerateMesh()
# Or: export netgen "mesh.vol" (surface-only via sideset)
mesh = Mesh(ngmesh)  # mesh.dim == 2
fes = HDivSurface(mesh, order=0)
# fes.ndof = number of surface edges only (all active in BEM)

# WRONG: dim=3 volume mesh (HDivSurface includes interior edges)
# Created by: export netgen "mesh.vol" with tet volume mesh
mesh = Mesh(ngmesh)  # mesh.dim == 3
fes = HDivSurface(mesh, order=0)
# fes.ndof = ALL edges (interior + boundary)
# FreeDofs = boundary edges only, but L matrix has zeros for interior
# -> rank-deficient L matrix -> wrong inductance (3000%+ error)
```

**For Cubit meshes**: Use `export netgen "mesh.vol" (surface-only via sideset)` which
creates a dim=2 surface mesh from Cubit surface blocks (triangles/quads only).

**For OCC (netgen.occ) meshes**: Use `Glue(shape.faces)` to extract surface:
```python
from netgen.occ import Box, Pnt, Glue, OCCGeometry
wire = Box(Pnt(0, -0.5e-3, -0.5e-3), Pnt(0.1, 0.5e-3, 0.5e-3))
# WRONG: OCCGeometry(wire) -> volume mesh -> BEM cond = 1e17
# CORRECT: surface-only
geo = OCCGeometry(Glue(wire.faces))
mesh = Mesh(geo.GenerateMesh(maxh=0.5e-3))  # maxh <= min_cross_section / 2
```

**maxh rule**: Set `maxh <= smallest_cross_section / 2` to avoid elongated
triangles. For 1mm wire: `maxh=0.5e-3`. Larger gives cond ~1e17 (singular).

## PITFALL: ds(label) Boundary Name Mismatch

When using `export netgen "mesh.vol" order N`, boundary labels come from
Cubit block names. Using `ds('conductor')` when the boundary block has
a different name causes LaplaceSL to **hang indefinitely** (no error, no timeout).

Always verify boundary names before using labeled ds:
```python
mesh = Mesh(ngmesh)
print(mesh.GetBoundaries())  # Check actual names!
# Then use: ds('conductor')  -- must match a name in GetBoundaries()
```

If labels don't match, use `ds` (no label) for all boundaries.

## CRITICAL: CalcSurfacesOfNode() After Manual Element2D Addition

When building a Netgen mesh manually, you MUST call `CalcSurfacesOfNode()`
after adding all Element2D elements. The `export netgen` command handles this internally.

Without this call, internal topology tables are not built, causing HDivSurface
edge orientation to be inconsistent.

```python
# export netgen handles this automatically:
import tempfile
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)

# For manual mesh construction (rare), call explicitly:
# ngmesh.CalcSurfacesOfNode()
# ngmesh.RebuildSurfaceElementLists()
```

## CRITICAL: TaskManager + use_fmm=False (Joachim Schoeberl, 2026-03-22)

**use_fmm=False** is required for:
1. Reproducible results (FMM causes non-deterministic floating-point summation)
2. Faster dense matrix extraction (FMM recomputes farfield on every MatVec)

**TaskManager** must be used for BOTH setup and extraction, or neither:

```python
# CORRECT: TaskManager wraps both setup and extraction, use_fmm=False
with TaskManager():
    L_op = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds
    SL = L_op.mat.ToDense().NumPy()  # Optimized dense extraction

# WRONG: TaskManager only for setup, or use_fmm=True (non-deterministic)
with TaskManager():
    L_op = LaplaceSL(u.Trace() * ds) * v.Trace() * ds  # FMM default -> fluctuates
```

## Dense Matrix Extraction: COO (NOT ToDense)

ngsolve.bem stores BEM matrices as `SparseMatrixdouble` with 100% fill
(every entry is nonzero).  `ToDense()` is ~2500x slower than necessary
because it internally performs N column-by-column MatVecs instead of a
direct memory copy.  Use `COO()` + scipy instead:

```python
from scipy.sparse import coo_matrix

# CORRECT: COO extraction (~0.06s at N=5085)
with TaskManager():
    L_op = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds
rows, cols, vals = L_op.mat.COO()
SL = coo_matrix((vals, (rows, cols)),
                shape=(L_op.mat.height, L_op.mat.width)).toarray()

# SLOW: ToDense() (~144s at N=5085, internally does N MatVecs)
SL = L_op.mat.ToDense().NumPy()  # ~2500x slower than COO

# SLOWEST: Manual column-by-column
for i in range(n):
    ei[:] = 0; ei[i] = 1.0
    L_op.mat.Mult(ei, col)  # Each call = O(N^2) kernel evaluation
```

**Benchmark** (N=5085 DOFs, LaplaceSL, use_fmm=False):

| Method | Time | Relative |
|--------|------|----------|
| Operator creation | 22s | (BEM integral assembly) |
| COO -> dense | 0.06s | 1x |
| ToDense() | 144s | 2500x slower |
| N x MatVec (manual) | 143s | ~= ToDense |

**Why SparseMatrix?** ngsolve.bem reuses NGSolve's FEM sparse matrix
infrastructure.  BEM matrices are dense by nature, but stored in CSR
sparse format with 100% fill.  This is a known design limitation.
Forum report: https://forum.ngsolve.org/t/...

## Self-Inductance: Source/Sink Saddle Point EFIE (Recommended)

For conductors with a gap (source/sink ports), use the constrained EFIE:

```
[SL  D^T] [J] = [0]
[D   0  ] [p] = [g]
```

where:
- SL = LaplaceSL matrix (HDivSurface)
- D = divergence matrix (HDivSurface -> SurfaceL2)
- g = source/sink current injection (+1/A_src at source, -1/A_snk at sink)
- J = surface current (unknowns), p = Lagrange multiplier

Inductance: `L = mu_0 * J^T @ SL @ J`

```python
from scipy.linalg import solve as scipy_solve

# Divergence matrix
bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
bf_D += div(u_J.Trace()) * q * ds
bf_D.Assemble()
D = bf_D.mat.ToDense().NumPy()

# Source/sink RHS (unit current injection)
f_src = LinearForm(fes_L2); f_src += q * ds("source"); f_src.Assemble()
f_snk = LinearForm(fes_L2); f_snk += q * ds("sink"); f_snk.Assemble()
g = f_src.vec.FV().NumPy() / sum(f_src.vec.FV().NumPy()) \
  - f_snk.vec.FV().NumPy() / sum(f_snk.vec.FV().NumPy())

# Saddle point solve (remove last constraint for regularity)
D_red, g_red = D[:-1, :], g[:-1]
K = np.block([[SL, D_red.T], [D_red, np.zeros((len(g_red), len(g_red)))]])
rhs = np.concatenate([np.zeros(n_J), g_red])
x = scipy_solve(K, rhs)
J = x[:n_J]
L = MU_0 * J @ SL @ J
```

**SurfaceL2 constraint**: Always use `order=0` (element-wise). Higher order
causes rank deficiency in D matrix.

## B-Distribution: Direct Biot-Savart

After BEM solve, compute B field in air volume via direct Biot-Savart:

```
B(x) = mu_0/(4*pi) * sum_e J_e x (x - c_e) / |x - c_e|^3 * A_e
```

**Do NOT use curl(A)**: Monopole A approximation + numerical curl = noise.
Direct B computation gives +0.2% accuracy at center of circular loop.

```python
# Per-element J extraction
elem_A = Integrate(CF(1), mesh, VOL_or_BND=BND, element_wise=True)
elem_Jx = Integrate(gf_J[0], mesh, VOL_or_BND=BND, element_wise=True)
# ... (extract centroids, areas, J vectors)

# Direct B at observation point
dx = obs - centroids  # (n_elem, 3)
r = np.sqrt(np.sum(dx**2, axis=1))
cross = np.cross(J_vecs, dx)  # J x (obs - centroid)
B = MU_0 * INV_4PI * np.sum(cross * (areas / r**3)[:, None], axis=0)
```

Volume mesh: OCC Box around conductor, Biot-Savart at each vertex.
Export via GmshPostExport(mesh_vol, boundary=False).

## PotentialCF: Cross-Mesh Limitation

`LaplaceSL(jt.Trace()*ds)` returns a PotentialOperator.
`pot_op(gf_J)` returns PotentialCF (dim=3 CoefficientFunction).

**PotentialCF does NOT work across meshes**: `gf.Set(A_cf)` on a different
volume mesh gives all zeros. Combined mesh (Box-Torus) causes DOF explosion
(354k DOFs) because `definedon` does not reduce DOF count.

Use direct Biot-Savart instead for volume B field evaluation.

## GMSH Visualization: Combined .geo

All results in one GMSH window via .geo that merges multiple .msh files:

```geo
// inductance.geo
Merge "inductance_B.msh";    // volume |B|, B
Merge "inductance_J.msh";    // surface |J|, J
Merge "inductance_coil.msh"; // coil wireframe (1D lines)
Mesh.NumSubEdges = 4;        // curved element display
Mesh.VolumeEdges = 0;        // clean volume rendering
```

GMSH Post-processing tree shows all views independently toggleable.
Coil wireframe (1D elements) is always visible on top of 3D volume.

## .geo Companion File

GmshPostExport automatically writes a companion .geo file alongside
each .msh when the mesh has Curve(2+). This sets `Mesh.NumSubEdges = 4`
for correct high-order element display (Tri6 appears curved, not flat).

## Legacy: Energy Method (Closed Loops Only)

For closed conductors (no gap), the energy method with toroidal current works:

```python
# Toroidal current for I=1: J = e_phi / (2*pi*a)
r_cf = sqrt(x*x + y*y)
J_toroidal = CF((-y/r_cf, x/r_cf, 0)) / (2 * math.pi * a)
gf_J = GridFunction(fes)
gf_J.Set(J_toroidal, definedon=mesh.Boundaries(".*"), dual=True)
J_vec = gf_J.vec.FV().NumPy().copy()
L_total = MU_0 * float(J_vec @ SL @ J_vec)
```

**Verified convergence** (R=50mm, a=5mm torus):

| Mesh | ndof | L (BEM) | L (Neumann) | Error |
|------|------|---------|-------------|-------|
| OCC cs=0.5 | 269 | 106.5 nH | 149.7 nH | -28.9% |
| OCC cs=1.0 | 1,790 | 142.1 nH | 149.7 nH | -5.0% |
| OCC cs=2.0 | 10,977 | 148.7 nH | 149.7 nH | -0.6% |
| Cubit (4611 DOF) | 4,611 | 148.7 nH | 149.7 nH | -0.7% |
| Cubit (8082 DOF) | 8,082 | 149.5 nH | 149.7 nH | -0.1% |

## Boundary Label Selection

```python
# Get available boundary labels
labels = mesh.GetBoundaries()
unique_labels = list(set(labels))
print(f"Boundaries: {unique_labels}")

# Use specific label
L_op = LaplaceSL(j.Trace() * ds("coil")) * j_test.Trace() * ds("coil")
```

## Surface Area Verification

Always verify mesh quality by checking surface area:

```python
area = Integrate(CF(1), mesh, VOL_or_BND=BND)
area_analytical = 4 * pi**2 * R * a  # Torus
error = abs(area - area_analytical) / area_analytical * 100
print(f"Area error: {error:.4f}%")
```

If area error > 1%, the inductance will also be inaccurate.
Use `mesh.Curve(order)` with SetGeomInfo to reduce geometric error.
"""

NGBEM_CUBIT_WORKFLOW = """
# Cubit -> ngsolve.bem BEM Inductance Pipeline

## Complete Workflow

```
1. Cubit: Create geometry (torus, cylinder, helix, ...)
2. Cubit: Tet or hex mesh
3. Cubit: Define blocks (domain, conductor surface)
4. Python: mesh = export netgen "mesh.vol" order 3
5. Python: ngsolve.bem LaplaceSL -> inductance extraction
```

No STEP files, no OCC geometry, no SetGeomInfo needed. export netgen
handles everything via Cubit's ACIS kernel + CallbackGeometry.

## Step-by-Step Code

```python
import sys, os, math, tempfile, numpy as np

# Import NGSolve BEFORE Cubit (DLL conflict avoidance)
from ngsolve import (Mesh, Integrate, CF, BND, HDivSurface, TaskManager,
                     GridFunction, ds, sqrt, x, y, z)
from ngsolve.bem import LaplaceSL

cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
import cubit

MU_0 = 4.0 * math.pi * 1e-7
R = 0.05   # Major radius [m]
a = 0.005  # Minor radius [m]

# --- Step 1: Create geometry ---
cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create torus major radius {R} minor radius {a}")

# --- Step 2-3: Mesh and define blocks ---
cubit.cmd("volume all scheme tetmesh")
cubit.cmd(f"volume all size {a/2}")
cubit.cmd("mesh volume all")
cubit.cmd('block 1 add volume 1')

# --- Step 4: Export with curving ---
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 2 overwrite')
mesh = Mesh(vol_path)

# --- Step 5: BEM inductance (energy method) ---
fes = HDivSurface(mesh, order=0)
u, v = fes.TnT()

# Toroidal current for I=1
r_cf = sqrt(x*x + y*y)
J_tor = CF((-y/r_cf, x/r_cf, 0)) / (2 * math.pi * a)
gf_J = GridFunction(fes)
gf_J.Set(J_tor, definedon=mesh.Boundaries(".*"), dual=True)
J_vec = gf_J.vec.FV().NumPy().copy()

# LaplaceSL + ToDense (use_fmm=False for speed and reproducibility)
with TaskManager():
    L_op = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds
    SL = L_op.mat.ToDense().NumPy()

# Energy method: L = mu_0 * J^T * SL * J
L_total = MU_0 * float(J_vec @ SL @ J_vec)
print(f"L = {L_total*1e9:.4f} nH")
```

## High-Order Curving: Area vs Inductance

p=2 curving dramatically improves **geometric accuracy** (surface area),
but BEM inductance accuracy is dominated by **DOF density** (mesh refinement):

| Curve Order | Surface Type | Area Improvement | L Improvement |
|-------------|-------------|-----------------|--------------|
| 1 (linear) | Flat facets | baseline | baseline |
| 2 (quadratic) | Parabolic patches | **4-56x better** | marginal |

**Key finding**: For fixed DOF count, p=2 improves area but not inductance
significantly. HDivSurface order=0 (RT0) current resolution is the bottleneck.
Mesh refinement (more DOFs) is more effective for L accuracy.

`export netgen "mesh.vol" order 2` provides quadratic-accurate surfaces
via ACIS CallbackGeometry. order=2 is the recommended default for BEM.


## Tri and Quad Surface Meshes

BEM supports both tri (from tet mesh) and quad (from hex mesh) surfaces.
`export netgen` exports volume mesh; BEM solver auto-extracts surface elements.

```python
# Tet mesh
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

# Export volume mesh (BEM solver extracts BND surface automatically)
import tempfile
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 2 overwrite')
mesh = Mesh(vol_path)
```

## Deleted APIs

The old workflow using `export_NetgenMesh()`, `set_*_geominfo()`,
STEP reimport, and OCCGeometry has been completely removed.
Use `export netgen` for all Cubit-to-NGSolve mesh transfers.
"""

NGBEM_CURVE_ORDER_STUDY = """
# Curve Order Study: Area vs Inductance

## Purpose

Demonstrates the effect of mesh curving (p=1 vs p=2) on BEM inductance.
Both tri (tet mesh) and quad (hex sweep mesh) surfaces are tested.

## Test Case: Torus (Cubit)

| Parameter | Value |
|-----------|-------|
| Major radius R | 50 mm |
| Minor radius a | 5 mm |
| R/a ratio | 10 |
| Analytical L | Neumann: L = mu_0*R*(ln(8R/a) - 2) = 149.7 nH |

## Measured Results (Cubit mesh -> export netgen -> LaplaceSL)

### Tri surface (tet mesh, gap torus)

| interval | p | nse | area_err | L_err | area improvement |
|----------|---|-----|----------|-------|------------------|
| 3 | 1 | 2498 | -2.55% | -5.95% | |
| 3 | 2 | 2498 | -0.36% | -6.37% | 7x |
| 8 | 1 | 2790 | -1.14% | -4.11% | |
| 8 | 2 | 2790 | -0.02% | -4.73% | 56x |

### Quad surface (hex sweep, 350-degree torus)

| config | p | nse | area_err | L_err | area improvement |
|--------|---|-----|----------|-------|------------------|
| 2x12 | 1 | 368 | -4.92% | -8.30% | |
| 2x12 | 2 | 368 | -1.16% | -7.79% | 4x |
| 3x18 | 1 | 784 | -2.86% | -6.20% | |
| 3x18 | 2 | 784 | -1.18% | -5.71% | 2x |

## Key Observations

1. **p=2 dramatically improves area** (4-56x), confirming curving works.

2. **L accuracy is dominated by DOF density**, not curve order.
   For fixed DOF count, p=2 doesn't improve L significantly.
   HDivSurface order=0 (RT0) current resolution is the bottleneck.

3. **Mesh refinement is more effective** for L than curving:
   OCC cs=2.0 (10977 DOF): L_err = -0.6%
   Cubit 8082 DOF: L_err = -0.1%

## Analytical References

### Neumann Formula (Circular Loop Self-Inductance)

```
L = mu_0 * R * (ln(8*R/a) - 2)
```

Valid for thin wire (a << R). For R=50mm, a=5mm: L ~ 135 nH.

### Torus Surface Area

```
A = 4 * pi^2 * R * a
```

### Torus Volume

```
V = 2 * pi^2 * R * a^2
```

## Full Example Script

See: `examples/cubit_panels/inductance/inductance_torus.py`

This script:
1. Creates a torus in Cubit
2. Tet meshes, defines blocks
3. Uses export netgen "mesh.vol" order N for each order
4. Compares order=1, order=2, order=3:
   - Surface area vs analytical
   - Volume vs analytical
   - BEM inductance vs Neumann formula
7. Prints summary convergence table
"""

NGBEM_STABILIZED = """
# Stabilized BEM for Low-Frequency Robustness

## The Low-Frequency Problem

Classical BEM with Helmholtz kernel:
```
V = V_1 - (1/kappa^2) * V_2
```

At low frequency (kappa -> 0), `1/kappa^2` diverges -> catastrophic
cancellation -> O(kappa^{-2}) condition number blow-up.

## Weggler's Stabilized Formulation

Lucy Weggler (ngsolve.bem developer) proposed a block stabilization:

```
[A_kappa,    Q_kappa   ] [J]     [b]
[Q_kappa^T,  kappa^2*V ] [rho] = [0]
```

Where:
- `A_kappa` = vector potential operator (loop part)
- `Q_kappa` = coupling operator (divergence constraint)
- `kappa^2 * V` = regularized scalar potential (NO 1/kappa^2!)
- `J` = surface current (HDivSurface)
- `rho` = surface charge (SurfaceL2)

At kappa=0 (MQS limit):
```
[A_0,    Q_0  ] [J]     [b]
[Q_0^T,  0    ] [rho] = [0]
```

This is a saddle-point system with O(1) condition number for all kappa.

## MQS Simplification

For MQS regime (DC to ~1 MHz), the system simplifies to:
- `A_0 = LaplaceSL` (Laplace kernel)
- `kappa^2 * V -> 0` (capacitive effects vanish)
- System reduces to **loop-only**: `A_0 * J = b`

This is exactly what our inductance extraction does:
```python
L_op = LaplaceSL(j.Trace() * ds("cond")) * j_test.Trace() * ds("cond")
```

## When You Need Stabilized BEM

- **MQS (inductance only)**: LaplaceSL is sufficient, no stabilization needed
- **Capacitive coupling**: Need finite kappa, use stabilized block system
- **Full-wave**: Need HelmholtzSL with stabilized formulation

## Implementation

```python
from ngsolve import *
from ngsolve.bem import LaplaceSL

# Product space: HDivSurface x SurfaceL2
fes_J = HDivSurface(mesh, order=0)
fes_rho = SurfaceL2(mesh, order=0)
fes = fes_J * fes_rho

(j, rho), (jt, rhot) = fes.TnT()

# Block system
with TaskManager():
    A = LaplaceSL(j.Trace() * ds) * jt.Trace() * ds      # [1,1] block
    Q = LaplaceSL(j.Trace() * ds) * rhot * ds             # [1,2] block
    V = LaplaceSL(rho * ds) * rhot * ds                   # [2,2] block (x kappa^2)
```

Reference: https://github.com/Weggler/docu-ngsbem/blob/main/demos/Maxwell_DtN_Stabilized.ipynb
"""

NGBEM_EXAMPLES = """
# ngsolve.bem Inductance Extraction Examples

## Example 1: Circular Loop (Netgen OCC, No Cubit)

Standalone example using Netgen OCC to create a torus (no Cubit needed):

```python
import math
import numpy as np
from netgen.occ import WorkPlane, Axes, Axis, Pnt, Dir, OCCGeometry
from netgen.meshing import MeshingParameters
from ngsolve import (Mesh, HDivSurface, TaskManager, ds, Integrate,
                     CF, BND, GridFunction, sqrt, x, y, z)
from ngsolve.bem import LaplaceSL

MU_0 = 4.0 * math.pi * 1e-7
R, a = 0.05, 0.005  # Major/minor radius [m]

# Analytical reference
L_ref = MU_0 * R * (math.log(8.0 * R / a) - 2.0)

# Create torus mesh via OCC
wp = WorkPlane(Axes(p=Pnt(R, 0, 0), n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
circle = wp.Circle(a).Face()
torus = circle.Revolve(Axis(p=Pnt(0, 0, 0), d=Dir(0, 0, 1)), 360)
geo = OCCGeometry(torus)
ngmesh = geo.GenerateMesh(mp=MeshingParameters(maxh=1.0, curvaturesafety=1.0))
mesh = Mesh(ngmesh)

# BEM inductance (energy method)
fes = HDivSurface(mesh, order=0)
u, v = fes.TnT()

r_cf = sqrt(x*x + y*y)
J_tor = CF((-y/r_cf, x/r_cf, 0)) / (2 * math.pi * a)
gf_J = GridFunction(fes)
gf_J.Set(J_tor, definedon=mesh.Boundaries(".*"), dual=True)
J_vec = gf_J.vec.FV().NumPy().copy()

with TaskManager():
    L_op = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds
    SL = L_op.mat.ToDense().NumPy()

L_total = MU_0 * float(J_vec @ SL @ J_vec)
error = abs(L_total - L_ref) / abs(L_ref) * 100
print(f"BEM L = {L_total*1e9:.2f} nH (ref: {L_ref*1e9:.2f} nH, error: {error:.1f}%)")
```

## Example 2: Cubit Torus with Source/Sink Port Detection

See `examples/induction_heating/bem_reference/bem_inductance.py`

Key workflow:
```python
# Cubit: create torus with gap, mesh
cubit.cmd(f"create torus major radius {R} minor radius {a}")
cubit.cmd(f"brick x {3*a} y {gap} z {3*a}")
cubit.cmd(f"move Volume 2 x {R} y 0 z 0 include_merged")
cubit.cmd("subtract volume 2 from volume 1")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")

# Register blocks and sidesets with source/sink
cubit.cmd('block 1 add volume 1; block 1 name "conductor"')
cubit.cmd('sideset 1 add surface {source_sid}; sideset 1 name "source"')
cubit.cmd('sideset 2 add surface {sink_sid};   sideset 2 name "sink"')

# Export with curving
import tempfile
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 2 overwrite')
mesh = Mesh(vol_path)

# Energy method
L_total = MU_0 * float(J_vec @ SL @ J_vec)
```
"""

NGBEM_BEST_PRACTICES = """
# ngsolve.bem Inductance Extraction Best Practices

## Checklist

1. **Always use .Trace()** on HDivSurface trial/test functions in BEM forms
   - Without .Trace(), boundary DOFs get corrupted silently
   - This is the #1 cause of wrong BEM results

2. **Verify surface area** before computing inductance
   - If area error > 1%, inductance will be inaccurate
   - Use `Integrate(CF(1), mesh, VOL_or_BND=BND)` vs analytical

3. **Set curvaturesafety for OCC meshes** (prevents degenerate elements)
   - `geo.GenerateMesh(maxh=maxh, curvaturesafety=1)`
   - Without it, curved surfaces may produce zero-eigenvalue BEM operators

4. **Use Curve(3) for curved surfaces**
   - Curve(1) = polygon approximation -> 5-15% inductance error
   - Curve(2) = quadratic -> ~0.5% error
   - Curve(3) = cubic -> ~0.05% error (usually sufficient)

5. **Use export netgen** for Cubit meshes
   - Handles curving automatically via ACIS CallbackGeometry
   - No SetGeomInfo, no STEP files, no OCC geometry needed

6. **Laplace kernel for MQS** (DC to ~1 MHz)
   - Do NOT use HelmholtzSL for inductance extraction
   - LaplaceSL is exact in the MQS regime

7. **HDivSurface order=0** (RT0) is sufficient for inductance
   - Higher order (order=1, 2) improves current distribution accuracy
   - But increases DOF count significantly (quadratic growth)
   - For total inductance (scalar), order=0 is usually adequate

8. **Check matrix rank** before solving
   ```python
   rank = np.linalg.matrix_rank(L_matrix)
   if rank < n_dof:
       print(f"WARNING: rank-deficient ({rank}/{n_dof})")
   ```

## Common Pitfalls

### Missing .Trace()
```python
# WRONG - silently gives wrong results
L_op = LaplaceSL(j_trial * ds("cond")) * j_test * ds("cond")

# CORRECT
L_op = LaplaceSL(j_trial.Trace() * ds("cond")) * j_test.Trace() * ds("cond")
```

### Forgetting mu_0 Factor
```python
# LaplaceSL gives the 1/(4*pi*r) kernel integral
# Must multiply by mu_0 for inductance
L_matrix = MU_0 * L_dense  # NOT just L_dense!
```

### Wrong Boundary Label
```python
# Check what labels exist
print(mesh.GetBoundaries())  # e.g., ('conductor', 'conductor', ...)

# Use the correct label
L_op = LaplaceSL(j.Trace() * ds("conductor")) * jt.Trace() * ds("conductor")
```

### Using Old API Instead of export netgen
```python
# OLD (DELETED): Do not use export_NetgenMesh or set_*_geominfo
# cubit_mesh_export.set_torus_geominfo(ngmesh, ...)  # DELETED

# NEW: Single APREPRO command handles everything
import tempfile
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)
```

## use_fmm: FMM vs Dense (No H-matrix in ngsolve.bem)

**ngsolve.bem has FMM only. H-matrix (ACA) is NOT implemented.**

| | `use_fmm=False` | `use_fmm=True` |
|--|-----------------|----------------|
| Matrix type | `SparseMatrixdouble` (100% fill) | `SumMatrix` (FMM + near-field sparse) |
| Structure | Dense matrix in CSR format | `FMM_Operator` (far) + `SparseMatrix` (near, 0.6% fill) |
| MatVec accuracy | Exact | ~2.6% relative error (multipole truncation) |
| Dense extraction | COO -> scipy (fast) | **Not possible** (SumMatrix has no COO) |
| LU / direct solve | Yes | **No** (matrix-free only) |
| `L = mu_0 * J^T SL J` | Yes (dense matrix) | **Not directly** |
| Best for | N < ~10,000 | N > ~10,000 (iterative solver) |

**Verified** (N=5067, LaplaceSL):
```python
# use_fmm=True returns SumMatrix:
#   SumMatrix = T^T @ FMM_Operator @ T  +  SparseMatrix(nze=146799)
# GetOperatorInfo():
#   SumMatrix
#     ProductMatrix (far-field: FMM_Operator LaplaceSL)
#     SparseMatrixdouble (near-field: nze=146799, 0.57% fill)
```

**FMM internals** (from `mptools.hpp`):
- Multipole expansion: `SingularMLExpansion` / `RegularMLExpansion`
- Octree: 8-child spatial decomposition
- Spherical harmonics basis
- `FMM_Parameters`: `maxdirect=100`, `minorder=20`

**POLICY: Use `use_fmm=False` for N < 10,000.**
For small-to-medium BEM (our inductance extraction target: N = 1,000-10,000):
- Dense matrix is needed for LU solve and energy inner product
- COO extraction is fast (0.18s at N=5085)
- FMM accuracy loss (2.6%) is unacceptable for engineering inductance
- FMM overhead exceeds dense assembly for small N

## Performance Tips

- BEM matrix assembly is O(N^2) where N = number of surface DOFs
- **COO extraction (NOT ToDense)**: `mat.COO()` + scipy is ~2500x faster than
  `mat.ToDense().NumPy()`.  ToDense() internally does N MatVecs.
  See "Dense Matrix Extraction" section above.
- **TaskManager**: Wrap operator setup for parallel BEM assembly.
  TaskManager gives ~5x speedup on 8 cores for BEM integral assembly.

## Performance Reference: BEM Assembly Times

**ngsolve.bem is competitive with state-of-the-art BEM libraries.**

Benchmark: Laplace single layer (P0/RT0), dense assembly, double precision:

| Library | N (DOF) | Cores | Time | Notes |
|---------|---------|-------|------|-------|
| **ngsolve.bem** | 5,085 | 8 (Xeon) | **21s** | TaskManager parallel |
| ngsolve.bem | 5,085 | 1 | 111s | Sequential |
| bempp-cl (PoCL) | ~5,000 | 8 (i9-9980HK, AVX-512) | ~15-20s | Hand-tuned OpenCL SIMD |
| bempp-cl (Numba) | 2,048 | 8 | 0.25s | (smaller problem) |
| bempp-cl (PoCL) | 2,048 | 8 | 0.08s | (smaller problem) |
| bempp-cl (PoCL) | 32,768 | 8 | ~20s | Figure 6, P0 basis |

Source: Betcke & Scroggs, "Designing a High-Performance Boundary Element
Library With OpenCL and Numba", IEEE CiSE 23(4), 2021.

**Key findings**:
- 5,000 DOF dense BEM in ~20s on 8 cores is typical/fast for current software
- ngsolve.bem (TaskManager) matches bempp-cl (hand-tuned AVX-512 OpenCL)
- Assembly is O(N^2) and dominated by Gauss quadrature (singular integrals
  use >1000 quadrature points per triangle pair)
- TaskManager scaling: 1t=111s, 2t=58s, 4t=31s, 8t=21s (5.2x on 8 cores)

**Design intent**: ngsolve.bem is designed for iterative solvers (FMM + CG).
For small-to-medium problems (N < 10,000) where the dense matrix IS needed
(LU, energy inner product), use `mat.COO()` extraction instead of `ToDense()`.

**Total pipeline time** (N=5085, inductance extraction):

| Step | Time | Notes |
|------|------|-------|
| BEM assembly | 21s | O(N^2) Gauss quadrature, TaskManager parallel |
| COO extraction | 0.18s | Was 144s with ToDense() |
| LU solve | 5s | scipy LAPACK (MKL) |
| **Total** | **~26s** | |

## Standard Output Format: GMSH .msh v4.1

**GMSH .msh v4.1 is the standard field output format** for this project.
All BEM/FEM field visualization uses `GmshPostExport`, not VTK/VTS.

**Why GMSH, not VTK**: GMSH natively supports arbitrary-order curved elements
(Tri6, Tri10, Tri15, ...). VTK approximates high-order elements as linear facets.

### Format Version Policy (2026-04: v4.1 only)

| Direction | Format | Purpose | Tool |
|-----------|--------|---------|------|
| **Input** (-> NGSolve) | **.vol** | Mesh import | `Mesh("model.vol")` |
| **Output** (NGSolve ->) | **.msh v4.1** | Field visualization | `GmshPostExport.write()` / `vol2msh()` -> GMSH GUI |

Lab-wide standard (2026-04): GMSH .msh v4.1 is the ONLY supported format.
`GmshPostExport.write_v22()` has been removed.  netgen I/O is always via
.vol (never via .msh).  v4.1 supports the same arbitrary-order elements
(Tri6, Tri10, Tri15, ...) for GMSH visualization.

### Supported Orders (Lagrange Triangles)

| Order | GMSH Type | Nodes | Use Case |
|-------|-----------|-------|----------|
| 1 | 2 (Tri3) | 3 | Flat mesh |
| 2 | 9 (Tri6) | 6 | Standard curved |
| 3 | 21 (Tri10) | 10 | High accuracy |
| 4 | 23 (Tri15) | 15 | Research |
| 5 | 25 (Tri21) | 21 | Research |

Export curved BEM surface mesh with field data to GMSH .msh v4.1 for visualization.
Use `GmshPostExport(mesh, boundary=True)` to export BND elements from a volume mesh.

### GMSH Display Setting for High-Order Elements

GMSH defaults to straight-line rendering between nodes. To see curved surfaces:

```
Mesh.NumSubEdges = 4;
```

Enter in GMSH console (Tools -> Command Line) or set in GUI:
Tools -> Options -> Mesh -> Visibility tab or General settings.

This subdivides each element edge into 4 segments for display.
- p=1 (Tri3): stays flat regardless of NumSubEdges
- p=2+ (Tri6, Tri10, ...): curved surface becomes visible
- Higher values (6-8) give smoother rendering at the cost of display speed

**Comparison tip**: Open p=1 and p=4 .msh files side by side on the
same coarse mesh (e.g. interval=4). With NumSubEdges=4, p=1 shows
polygonal facets while p=4 shows smooth curved surfaces.

```python
from radia.gmsh_post_export import GmshPostExport
post = GmshPostExport(mesh, boundary=True)  # boundary=True for BND from volume mesh
post.add_field("|J|", node_J, ncomp=1)
post.write("results.msh")
```

### Extracting High-Order Nodes: GetTrafo + GMSH Reference Coordinates

To output Tri6/Tri10/Tri15/Tri21 elements with correct curved positions,
evaluate `mesh.GetTrafo(el)` at GMSH's equidistant reference points:

```python
import gmsh
gmsh.initialize()
p = mesh.GetCurveOrder()  # e.g. 2, 3, 4, 5
etype = gmsh.model.mesh.getElementType('Triangle', p)
props = gmsh.model.mesh.getElementProperties(etype)
ref_coords = np.array(props[4]).reshape(props[3], -1)  # (n_nodes, 2)
gmsh.finalize()

# For each BND element, evaluate transformation at GMSH reference points
for el in mesh.Elements(BND):
    trafo = mesh.GetTrafo(el)
    for i in range(3, len(ref_coords)):  # skip 3 corners
        u, v = ref_coords[i]
        ir = IntegrationRule([(u, v)], [1.0])
        for ip in ir:
            mip = trafo(ip)
            x, y, z = mip.point[0], mip.point[1], mip.point[2]
```

**H1 DOF structure for order p on TRIG**:
- Vertex DOFs: 3 (same as corners)
- Edge DOFs: (p-1) per edge, 3 edges → 3*(p-1) total
- Interior DOFs: (p-1)*(p-2)/2
- Total: (p+1)*(p+2)/2 = number of GMSH Lagrange triangle nodes

| Order p | Vertex | Edge (3 edges) | Interior | Total | GMSH Type |
|---------|--------|----------------|----------|-------|-----------|
| 1 | 3 | 0 | 0 | 3 | 2 (Tri3) |
| 2 | 3 | 3 | 0 | 6 | 9 (Tri6) |
| 3 | 3 | 6 | 1 | 10 | 21 (Tri10) |
| 4 | 3 | 9 | 3 | 15 | 23 (Tri15) |
| 5 | 3 | 12 | 6 | 21 | 25 (Tri21) |

**Why GetTrafo, not H1 GridFunction?**
- `GridFunction.Set()` does L2 projection + averaging, which is exact for p=2
  but **breaks for p>=4** (coordinates drift far from the surface).
- `GetTrafo` evaluates the exact curved mapping at any reference point.
- Edge nodes are cached across shared elements (keyed by vertex pair).
  Direction is corrected by matching ref corners to physical vertices.
- GMSH reference coordinates obtained via `gmsh.model.mesh.getElementProperties()`,
  guaranteeing exact match with GMSH's node ordering.

**NGSolve TRIG edge ordering**: `el.edges[i]` is **opposite** to `el.vertices[i]`.
- edge[0] connects (verts[1], verts[2])
- edge[1] connects (verts[0], verts[2])
- edge[2] connects (verts[0], verts[1])

**Edge DOF direction**: For p >= 3, each edge has multiple DOFs. They must be
sorted by distance from the starting vertex to match GMSH's edge node ordering
(equidistant from start to end vertex). `GmshPostExport` handles this internally.

**Interior DOF ordering**: For p >= 3, interior DOFs exist. GMSH has a specific
recursive ordering pattern. For p=3 (1 interior node) this is trivial. For p >= 4,
the NGSolve DOF order may not match GMSH exactly — visually acceptable but not
formally verified for p > 3.

**Vertex matching for GetTrafo**: NGSolve BND elements have a vertex-to-reference
mapping that varies per element. GetTrafo is evaluated at ref corners (0,0), (1,0),
(0,1) and matched to physical vertex positions via nearest-neighbor. This
determines the correct ref-to-physical permutation for each element.

**Do NOT use H1 GridFunction.Set() for p>=4**: L2 projection + averaging corrupts
high-order node coordinates. GetTrafo is exact for all orders.

## Validation: Neumann Formula

For circular loops, validate against the analytical Neumann formula:

```
L = mu_0 * R * (ln(8*R/a) - 2)
```

| R [mm] | a [mm] | R/a | L [nH] |
|--------|--------|-----|--------|
| 50 | 5 | 10 | 149.7 |
| 50 | 2 | 25 | 193.2 |
| 50 | 1 | 50 | 237.5 |
| 100 | 10 | 10 | 299.3 |

Higher R/a ratios mean thinner wires -> more accurate Neumann formula.
For R/a < 5, mutual inductance corrections become significant.
"""


NGBEM_HODGE_DECOMPOSITION = """
# Hodge Decomposition for Closed Conductor Inductance

## Motivation

For **closed conductors** (no gap, no source/sink ports), the source/sink
saddle point EFIE cannot be applied. The Hodge decomposition approach extracts
toroidal and poloidal current modes from the topology of the surface mesh.

This is useful for:
- Closed loops (torus, ring-shaped conductors)
- Multi-turn coils with no explicit port
- Topological mode analysis (genus >= 1 surfaces)

## Mathematical Background

Surface currents decompose into 3 orthogonal subspaces (Hodge-Helmholtz):

```
J = J_grad + J_curl + J_harm
```

- `J_grad = grad(phi)`: gradient (irrotational) -- carries no net flux
- `J_curl = curl(psi)`: co-exact -- localized vortices
- `J_harm`: harmonic -- topological currents (toroidal, poloidal)

For a genus-g surface, dim(harmonic) = 2g. A torus (g=1) has exactly
2 harmonic modes: toroidal and poloidal.

## Discrete Hodge Decomposition

Using HDivSurface (RT0) basis on a surface mesh:

```
D = divergence matrix (n_faces x n_edges)  -- SurfaceL2 test x HDivSurface trial
C = incidence matrix (n_edges x n_vertices) -- edge-vertex connectivity
M_J = mass matrix (n_edges x n_edges)       -- HDivSurface inner product
```

The harmonic subspace is:

```
V_harm = null_space([D; C^T @ M_J])
```

- `D @ J = 0`: divergence-free (no sources/sinks)
- `C^T @ M_J @ J = 0`: curl-free in the L2 sense

## Eigenvalue Problem for Inductance

Project LaplaceSL onto the harmonic subspace:

```python
SL_harm = V_harm.T @ SL @ V_harm   # projected single-layer
M_harm  = V_harm.T @ M_J @ V_harm  # projected mass

eigvals, eigvecs = scipy.linalg.eigh(SL_harm, M_harm)
```

For a torus (2 harmonic modes):
- **Mode 0** (smaller eigenvalue): poloidal current
- **Mode 1** (larger eigenvalue): toroidal current -> **this gives inductance**

```
L = mu_0 * eigvals[1] * R / a
```

The R/a factor normalizes by the current density (J ~ 1/(2*pi*a) for unit current).

## Implementation Sketch

```python
from scipy.linalg import null_space, eigh

# Spaces
fes_J  = HDivSurface(mesh, order=0)
fes_L2 = SurfaceL2(mesh, order=0)
n_J, n_f, nv = fes_J.ndof, fes_L2.ndof, mesh.nv

# Divergence matrix D (n_f x n_J)
bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
bf_D += div(u_J.Trace()) * q * ds
bf_D.Assemble()
D = bf_D.mat.ToDense().NumPy()

# Incidence matrix C (n_J x nv): edge -> vertex connectivity
C = np.zeros((n_J, nv))
for i, e in enumerate(mesh.edges):
    vv = list(e.vertices)
    if i < n_J:
        C[i, vv[0].nr] = -1
        C[i, vv[1].nr] = +1

# Mass matrix M_J (n_J x n_J)
bf_M = BilinearForm(fes_J)
bf_M += InnerProduct(u.Trace(), v.Trace()) * ds
bf_M.Assemble()
M_J = bf_M.mat.ToDense().NumPy()

# Harmonic subspace
V_harm = null_space(np.vstack([D, C.T @ M_J]), rcond=1e-10)
# For torus: V_harm.shape[1] == 2

# LaplaceSL
with TaskManager():
    V_op = LaplaceSL(jt.Trace() * ds, use_fmm=False) * jv.Trace() * ds
    SL = V_op.mat.ToDense().NumPy()

# Eigenvalue problem on harmonic subspace
SL_harm = V_harm.T @ SL @ V_harm
M_harm  = V_harm.T @ M_J @ V_harm
eigvals, eigvecs = eigh(SL_harm, M_harm)

# Toroidal mode -> inductance
L = MU_0 * eigvals[1] * R / a
```

## Practical Notes

1. **Closed surface required**: The mesh must be a closed surface (no boundary
   edges). Use `surface_only=True` in `export netgen` or OCC `Glue()`.

2. **Euler characteristic**: For genus-g: V - E + F = 2 - 2g. Torus: Euler = 0.

3. **Mode identification**: The poloidal mode has smaller eigenvalue than the
   toroidal mode. Visualize |J| in GMSH to confirm which is which.

4. **Comparison with source/sink**: For conductors with a gap, the source/sink
   saddle point EFIE is simpler and more robust. Hodge decomposition is mainly
   useful when there is no natural port location.

5. **Normalization**: The eigenvalue-to-inductance conversion `L = mu_0 * lambda * R/a`
   depends on the geometry. For non-circular cross-sections, a different normalization
   is needed (e.g., integrate the toroidal mode to get total current).
"""


NGBEM_ESIM_WORKPIECE = """
# EFIE-SIBC: Two-Way Eddy Current BEM Solver

## Overview

EFIE-SIBC provides **two-way coupled** eddy current analysis for induction heating.
The workpiece surface currents are solved self-consistently with the coil excitation,
capturing the screening effect that one-way models miss.

**File**: `examples/cubit_panels/inductance/efie_sibc.py`
**Verification**: `examples/cubit_panels/inductance/verify_laplace_bem.py`

## Formulation (Saddle-Point EFIE + SIBC)

```
[Z_s*M + jw*mu0*SL,  D^T] [J] = [-jw * A_inc]
[D,                   0  ] [p]   [0          ]
```

- **SL**: `LaplaceSL(u.Trace()*ds, use_fmm=False) * v.Trace()*ds` (self-inductance of surface currents)
- **M**: HDivSurface mass matrix (surface impedance coupling)
- **D**: divergence matrix HDivSurface -> SurfaceL2 (enforces div J = 0)
  - Closed surface: rank = n_elem - 1 (one redundant row removed)
- **Z_s**: per-element complex surface impedance from ESIM cell problem
- **A_inc**: incident vector potential from coil via Biot-Savart:
  `A(r) = (mu0/4pi) * sum(J_coil / |r-r'| * dA)`
- **Karl iteration**: solve -> H_t = |J| per element -> update Z_s from ESIM -> repeat
  - Relaxation factor 0.5, convergence in 4-5 iterations (dZ/Z < 1e-3)

## Critical: LaplaceDL is NOT MFIE

`LaplaceDL` in ngsolve.bem for HDivSurface computes the **scalar double layer**
(dG/dn_y applied to vector basis), with eigenvalue **-1/6** on unit sphere l=1 mode.

This is **NOT** the MFIE K operator `n x curl(SL)` which would have eigenvalue +1/6.
The MFIE K operator is not available in current ngsolve.bem.

**Consequence**: The original PMCHWT formulation `(1/2 M + K)*J = n x H_inc` was
incorrect because `LaplaceDL != MFIE K`. The correct formulation is the EFIE above.

**Validation snapshot (2026-06-25)**:
`examples/cubit_panels/inductance/verify_laplace_bem.py` on a unit sphere
(`maxh=0.25`, 456 surface elements, 684 HDivSurface DOFs) confirms:

- scalar `LaplaceDL` l=1 eigenvalue: `-0.166384` vs exact `-1/6`
- EFIE PEC current ratio: `J/J_pec = 1.002233`
- EFIE PEC shape error: `1.7465e-02`
- SIBC screening sweep: `|J/J_pec| = 1.0022` at `Z_s=1e-7`, `2.63e-4` at `Z_s=10`
- Overall: PASS

## KNOWN LIMITATION: BEM EFIE-SIBC Incorrect for Finite Z_s (2026-03-28)

**The EFIE `Z_s*J + jw*mu0*SL(J) = -jw*A_inc` is fundamentally wrong for finite Z_s
in MQS.** The representation `A_scat = mu0*SL(J_s)` uses the Laplace SL operator,
whose eigenvalue for l=1 on a sphere of radius R is R/3 (not R). This gives an
effective denominator of `(3*Z_s + jw*mu0*R)` instead of `(Z_s + jw*mu0*R)`.

**BEM/Analytical ratio** on sphere for various Z_s/(jw*mu0*R):

| Z_s/(jw*mu0*R) | BEM/Analytical | Notes |
|-----------------|----------------|-------|
| 0.01 (PEC-like) | ~1.00 | Z_s negligible, EFIE correct |
| 0.1 | ~0.97 | Still acceptable |
| 1.0 | ~0.75 | Factor-of-3 effect visible |
| 5.0 | ~0.50 | EFIE significantly wrong |
| 10.0 | ~0.35 | Steel at 7 kHz regime |

**Only correct for PEC (Z_s -> 0)** where the Z_s*M term vanishes.

**Fix requires MFIE** (not available in ngsolve.bem) or a Stratton-Chu formulation
with both SL and DL operators. Until MFIE is available, use **FEM-SIBC** instead
for finite Z_s problems.

**Diagnostic script**: `examples/cubit_panels/inductance/efie_sibc.py`
**Verification script**: `examples/cubit_panels/inductance/verify_sphere_sibc.py`

## FIX: Scalar Potential BIE + SIBC (2026-03-29)

**The scalar potential BIE correctly handles all Z_s values using EXISTING
ngsolve.bem operators.** No new C++ code needed.

### System
```
(1/2*M - DL + gamma * SL * M^{-1} * K) phi = rhs
gamma = Z_s / (jw * mu0)
```

- `M` = H1 surface mass matrix
- `K` = H1 surface stiffness (Laplace-Beltrami)
- `DL` = LaplaceDL (scalar, H1 trial/test)
- `SL` = LaplaceSL (scalar, H1 trial/test)
- Gauge: Lagrange multiplier for `int(phi) dS = 0`
- Unknown `phi` = exterior scalar magnetic potential (H = -grad phi)

### Physics
SIBC enters via surface Laplacian:
```
E_t = Z_s * J_s = -Z_s * (n x grad_s phi)
Faraday: dphi/dn = -(Z_s / (jw*mu0)) * Delta_s(phi)
Weak form: M*g = gamma * K * phi  (integration by parts)
```

### Why EFIE (HDivSurface) fails but Scalar BIE (H1) works
RT0 (HDivSurface order=0) has **zero surface curl**: curl_s(J) = 0 for RT0
basis functions. The SIBC requires curl_s(J)·n = Delta_s(phi) which vanishes
for RT0. The H1 scalar potential avoids this because grad_s and Delta_s
are well-defined for H1 (order >= 1).

### Validated accuracy
- Sphere: <0.1% error for ALL Z_s (PEC to transparent)
- Cylinder workpiece + coil: +7% vs FEM-SIBC (mesh density dependent)

### Usage (reusable solver module)
```python
from radia.bem_sibc_solver import ScalarBIESIBCSolver

solver = ScalarBIESIBCSolver(mesh_wp, order=1)   # assemble once
result = solver.solve(phi_inc, Z_s=Z_s, omega=omega)  # fast per Z_s
H_t_rms = result['H_t_rms']
P_density = result['P_density']   # [W/m^2]
```

### phi_inc from coil surface current
For the incident scalar potential at workpiece from EFIE-solved coil current:
```python
from radia.bem_sibc_solver import compute_phi_inc_from_surface_J

# Extract per-element J from EFIE solution
coil_c, coil_a, coil_J = extract_element_J(mesh_coil, gf_J)
# Compute phi_inc via Biot-Savart + two-stage path integration
phi_inc = compute_phi_inc_from_surface_J(wp_nodes, coil_c, coil_a, coil_J)
```

**NOTE**: ngsolve.bem does not have a grad(G) kernel for cross-mesh evaluation.
LaplaceDL gives dG/dn' (normal component only), which is insufficient for
cross-mesh Biot-Savart (n_src != n_obs). Feature request submitted to Joachim.
Direct numerical Biot-Savart (vectorized NumPy) is used as workaround.

### Coupled two-body BEM (coil + workpiece) — sign-correct (2026-04-12)

`bem_coupled_solver.CoupledBEMSolver` implements iterative coil EFIE
+ workpiece scalar BIE+SIBC with **per-DOF back-reaction RHS**.

```python
from radia.bem_coupled_solver import CoupledBEMSolver

solver = CoupledBEMSolver(mesh_coil, mesh_wp)
result = solver.solve(Z_s=Z_s, omega=omega, max_iter=10, tol=1e-3,
                      relax=0.5)
# result keys:
#   L_air      coil-only inductance (no workpiece)
#   L_total    coupled coil terminal inductance (= L_air + Delta_L)
#   Delta_L    coil terminal inductance change due to workpiece
#   P_total    workpiece eddy power loss [W]
#   H_t_rms    surface tangential H rms [A/m]
#   iterations number of Picard iterations to convergence
#   J_coil_re  back-reacted coil current (real part)
#   J_coil_im  back-reacted coil current (imaginary part)
```

#### Picard iteration

```
1. J_coil(0) = LU^{-1}(g_red)              # uncoupled (air-only) solve
2. for k = 0..max_iter:
     phi_inc = Biot-Savart(J_coil) at workpiece nodes
     phi_wp  = ScalarBIESIBCSolver.solve(phi_inc, Z_s, omega)
     J_wp    = n x H_scat = n x (-grad_s(phi_wp - phi_inc))
     f_back[i] = int v_i.Trace() . A_wp dS_coil          # PER-DOF
     J_coil(k+1) = LU^{-1}(g_red - f_back)               # re-solve
     under-relax with relax=0.5
   until |L_total - L_prev| / |L_prev| < tol
```

#### Back-reaction RHS (the 2026-04-12 fix)

Built as a NGSolve LinearForm that integrates the workpiece-induced
vector potential against the HDivSurface test space on the coil:

```python
A_components = [CoefficientFunction(0.0) for _ in range(3)]
for j in range(M):                          # M ~= 280 wp panels
    cx, cy, cz = wp_c[j]
    r_inv = 1.0 / sqrt((x-cx)**2 + (y-cy)**2 + (z-cz)**2)
    weight = (mu0/4pi) * wp_a[j] * r_inv
    for k in range(3):
        A_components[k] += weight * wp_J[j, k]
A_cf = CoefficientFunction(tuple(A_components))

f_form = LinearForm(fes_J)
f_form += InnerProduct(u_J.Trace(), A_cf) * ds
f_form.Assemble()
f_back = f_form.vec.FV().NumPy().copy()
```

This is a per-HDivSurface-DOF vector. **Do NOT use a scalar rescale
of `SL @ J_coil`** — the previous (v1) implementation did exactly that
and produced wrong-signed Delta_L (saved as
`bem_coupled_solver_v1_buggy.py.bak` for reference).

#### L_total formula

```python
L_self = mu0 * (J_re^T SL J_re + J_im^T SL J_im)
mutual = f_back_re . J_re + f_back_im . J_im
L_total = L_self + mutual
Delta_L = L_total - L_air
```

The sum of L_self drop + negative mutual gives the Lenz reduction for
non-magnetic conductors, while a positive mutual (from large Im(Z_s)
when mu_r >> 1) gives the flux concentration limit for ferromagnetic
workpieces.

#### Verified sign behavior (2026-04-12)

Frequency sweep (copper, mu_r=1, R_coil=30mm, R_wp=10mm, H_wp=20mm):

| freq    | delta    | L_air     | L_total   | Delta_L      |
|---------|----------|-----------|-----------|--------------|
| 100 Hz  | 6.61 mm  | 86.671 nH | 86.356 nH | -0.316 nH    |
| 1 kHz   | 2.09 mm  | 86.671 nH | 85.899 nH | -0.772 nH    |
| 10 kHz  | 0.66 mm  | 86.671 nH | 85.729 nH | -0.942 nH    |
| 100 kHz | 0.21 mm  | 86.671 nH | 85.673 nH | -0.998 nH    |
| 1 MHz   | 0.066 mm | 86.671 nH | 85.655 nH | -1.016 nH    |

Negative for all frequencies (Lenz screening), monotonically growing
with frequency, asymptoting to the PEC limit at high frequency.

mu_r sweep (steel, sigma=2e6, half=5mm, f=50 kHz):

| mu_r | |Z_s|       | L_air     | L_total   | Delta_L      |
|------|-------------|-----------|-----------|--------------|
| 1    | 4.43e-04 Ω  | 86.671 nH | 85.841 nH | -0.831 nH    |
| 10   | 1.41e-03 Ω  | 86.671 nH | 86.209 nH | -0.462 nH    |
| 100  | 4.44e-03 Ω  | 86.671 nH | 87.013 nH | **+0.342 nH** |
| 1000 | 1.41e-02 Ω  | 86.671 nH | 87.976 nH | +1.304 nH    |

Sign change between mu_r=10 and mu_r=100, marking the cross-over from
Lenz-dominated screening to flux-concentration-dominated storage in the
ferromagnetic skin layer.

### Panel integration

`calc_inductance.py` calls `_run_coupled_bem` when the user picks
`Method=BEM, Workpiece=SIBC` in the IH panel. The headline
`inductance_H` becomes `L_total`, and `L_air_H` keeps the coil-only
value so the panel can display both. `radia_gui_base.py::_on_finished`
shows L (air), delta L, L (eff), R (added), iters, skin depth.

For ESIM (nonlinear cell problem) the coupled solver is NOT used —
falls back to one-way uncoupled estimator with R reported only.

### Files
- `src/radia/bem_sibc_solver.py`: ScalarBIESIBCSolver + phi_inc computation
- `src/radia/bem_coupled_solver.py`: CoupledBEMSolver (per-DOF f_back, sign correct)
- `src/radia/bem_coupled_solver_v1_buggy.py.bak`: archived buggy v1 (scalar rescale)
- `src/radia/panels/calc_inductance.py::_run_coupled_bem`: wrapper for IH panel
- `examples/cubit_panels/inductance/scalar_bie_sibc.py`: sphere validation
- `examples/cubit_panels/inductance/bem_sibc_workpiece.py`: coil + workpiece demo
- `examples/cubit_panels/inductance/experiment_coupled_bem.py`: FEM-ESIM cross-check (still uses old API)

## Screening Physics

The key dimensionless parameter is `Z_s / (jw * mu0 * a)` where `a` is the workpiece
characteristic size (radius for cylinder).

| Z_s / (jw*mu0*a) | Behavior | One-way accuracy |
|-------------------|----------|-----------------|
| < 0.3 | Weak screening | One-way OK (-11%) |
| 0.3 - 3 | Transition | One-way unreliable |
| > 3 | Strong screening | **One-way fails (100x+ error)** |

One-way models (BEM-ESIM, FEM-ESIM) use `H_t = H_inc` (PEC approximation).
The correct surface current is:
```
H_t = omega * mu0 * H_inc * R_wp / (2 * Z_s)   (EMF / loop impedance)
```
For steel at 7kHz: H_t = 0.77 A/m, not 18 A/m. One-way overestimates P by 300x.

## Validated Results (EFIE-SIBC BEM)

| Condition | Z_s/(jw*mu0*a) | H_t [A/m] | P [W] | vs one-way |
|-----------|----------------|-----------|-------|------------|
| Sphere PEC (Z_s=0) | 0 | J/J_pec = 1.002 | - | 0.2% error |
| Copper 1kHz (xi=4.8) | 0.1 | 13.93 | 1.35e-6 | -11% |
| Steel 7kHz (xi=2.4) | 8.9 | 0.77 | 1.96e-6 | -99.7% |

## Cross-Validation: EFIE-SIBC (BEM) vs FEM-SIBC

| Material | Freq | EFIE-SIBC (BEM) | FEM-SIBC | Diff |
|----------|------|-----------------|----------|------|
| Steel | 7 kHz | 1.96e-6 W | 1.76e-6 W | -9.9% |
| Copper | 1 kHz | 1.35e-6 W | 1.26e-6 W | -6.8% |

Both methods agree within ~10%. FEM-SIBC is slightly lower due to the thin-shell
approximation (zero-thickness interface vs BEM's SL self-inductance).

**CAUTION**: The BEM EFIE-SIBC values above are affected by the SL eigenvalue
mismatch (see "KNOWN LIMITATION" section above). For finite Z_s, BEM EFIE-SIBC
underestimates the true surface current. The cross-validation agreement is partly
coincidental -- both methods have different sources of error that partially cancel.
For authoritative results, use the analytical solution on spheres or FEM-SIBC with
the total-field formulation on general geometries.

## Usage

```python
# Standalone
python efie_sibc.py --material steel --freq 7000
python efie_sibc.py --material copper --freq 1000

# As module
from efie_sibc import run
result = run(material='steel', frequency=7000)
# result keys: P_total, Q_total, H_t_rms, Z_s_elem, J_sol, ndof, ne, area
```

## Pipeline

```
Coil BEM (source/sink EFIE)
  -> gf_J (coil surface current for 1A)
  -> Biot-Savart -> A_inc, H_inc at workpiece

Workpiece surface mesh (OCC Cylinder -> Glue -> GenerateMesh)
  -> HDivSurface order=0, SurfaceL2 order=0

EFIE-SIBC saddle point + Karl iteration
  -> J (surface current), Z_s (per-element impedance)
  -> P = sum(P'_i * area_i), Q = sum(Q'_i * area_i)
```

## FEM-SIBC: Two-Way FEM with Surface Impedance

**File**: `examples/cubit_panels/inductance/fem_esim_3d.py`

FEM-SIBC solves `curl(nu*curl(A)) = J_source` with SIBC penalty on an internal
interface.

### Critical: Scattered-Field RHS Must Include Both Terms (2026-03-28)

In a **scattered-field** formulation `A = A_inc + A_scat`, the FEM RHS must include
**both** surface terms:

```
f(v) = -(jw/Z_s) * <A_inc, v>_sibc   (surface impedance term)
     + (-1)      * <n x H_inc, v>_sibc (incident field boundary term)
```

The second term `<n x H_inc, v>` was missing and caused a **factor-of-3 error** on
the sphere benchmark. This term arises from the curl-curl integration by parts:
`-int nu*curl(A_inc)*curl(v) dx = -<n x H_inc, v>_boundary`.

For the Kelvin formulation (mapped domain), the boundary form is preferred over the
volume form to avoid incorrect contributions from the Kelvin domain.

**The total-field formulation** (used in `fem_esim_3d.py` with coil current source)
does NOT have this issue because the source is `J_source` in the volume, not a
scattered-field decomposition.

**Verification**: `verify_sphere_sibc.py` -- FEM vs BEM vs analytical on sphere.

**Fix applied in**: `verify_sphere_sibc.py` line 162-163.

### SIBC = Robin BC on Conductor Surface (2026-04-14 Update)

**SIBC is a Robin BC.  Conductor interior is NOT solved.**

Use hole approach: subtract workpiece from mesh, Robin BC on hole boundary.
Validated 2026-04-14 (2D axisym Kelvin): L < 1%, P < 2% for Cu/Steel/Al.

```python
# Hole approach (correct): workpiece removed from mesh
air = air_sphere - torus - wp_cyl    # wp_cyl subtracted = hole
shape = Glue([air, torus])
# Robin BC on hole boundary:
a += (1j * omega / Z_s) * u.Trace() * v.Trace() * ds("sibc")
# Z_s for solid cylinder: rho*gamma*I1(ga)/I0(ga)
```

### Critical: H_t from SIBC Relation, NOT curl(A)

**WRONG** (gives incident field ~18 A/m, not surface current ~0.7 A/m):
```python
H_cf = nu0 * curl(gfu)
H_mag_sq = sum(H_cf[i].real**2 + H_cf[i].imag**2 for i in range(3))
H_t_rms = sqrt(Integrate(H_mag_sq, mesh, BND, definedon=wp_region) / A_wp)
# This gives the TOTAL tangential H dominated by incident coil field.
```

**CORRECT** (gives physical surface current for ESIM input):
```python
# J_s = -(jw/Z_s) * A_t  =>  H_t = |jw/Z_s| * |A_t|
At_sq = sum(gfu[i].real**2 + gfu[i].imag**2 for i in range(3))
At_rms = sqrt(Integrate(At_sq, mesh, BND, definedon=wp_region) / A_wp)
H_t_rms = abs(1j * omega / Z_s) * At_rms
# H_t is the thin-shell surface current density = ESIM input H_0
```

**Why curl(A) is wrong on internal interface**: On an internal SIBC interface,
`curl(A)/mu0` from adjacent elements gives the local H field in each element.
This includes the incident field from the coil (~18 A/m) plus a small scattered
contribution. The ESIM cell problem needs the physical surface current
`J_s = H_t_surface`, not the total H. The SIBC relation `J_s = -(jw/Z_s)*A_t`
extracts the correct surface current from the continuous tangential A.

### Formulation

```
int nu0 * curl(A) . curl(v) dx + (jw/Z_s) * int A_t . v_t ds("wp_surface")
  = int J_source . v dx("coil")

dx: all domains (air + workpiece_interior + coil)
ds("wp_surface"): internal interface between air and workpiece
```

**Thin conducting shell model**: The penalty term represents a jump condition
`[n x H] = J_s = -(jw/Z_s) * A_t` on the internal interface.

### ESIM Geometry: Slab vs Cylinder

The `--geometry` flag selects the 1D cell problem ODE:
- `cylinder` (default): Bessel I0/I1, for cylindrical workpieces
- `slab`: cosh/sinh, for flat plates

When delta/R < 0.1 (thin skin): slab ~ cylinder (<2% difference).
When delta/R > 0.1: curvature matters (up to 11% P' difference for delta/R=0.2).

### Usage

```python
python fem_esim_3d.py --material steel --freq 7000
python fem_esim_3d.py --material copper --freq 1000 --geometry slab

from fem_esim_3d import run
result = run(material='steel', frequency=7000, sigma=2e6, esim_geometry='cylinder')
# result keys: P_total, Q_total, L, H_t_rms, Z_s, ndof, ne
```

### When to Use FEM-SIBC vs EFIE-SIBC (BEM)

| Criterion | FEM-SIBC | EFIE-SIBC (BEM) |
|-----------|----------|-----------------|
| Accuracy (finite Z_s) | Reference (total-field) | **Wrong** (SL eigenvalue, see KNOWN LIMITATION) |
| Accuracy (PEC, Z_s->0) | Good | Good (0.2% on sphere) |
| Speed (this geometry) | ~300s | ~5s |
| Magnetic core coupling | Easy (FEM volume) | Needs FEM-BEM coupling |
| Complex coil geometry | Easy (volume J source) | Needs Biot-Savart |
| Self-inductance of eddy current | Approximate (thin shell) | SL operator (but eigenvalue wrong) |
| Recommended for | **All SIBC problems** | PEC only, or quick estimates |

**Note**: FEM-SIBC with total-field formulation (`fem_esim_3d.py`) is now the
recommended method for two-way SIBC analysis. BEM EFIE-SIBC (`efie_sibc.py`) is
valid only for PEC (Z_s -> 0) or as a quick estimate with known SL bias.
For scattered-field FEM, ensure both RHS terms are included (see above).
"""


NGBEM_MQS_FORMULATION_LIMITS = """
# BEM Formulation Selection for MQS-SIBC (2026-03-28)

## Overview

BEM formulations for eddy currents with Surface Impedance Boundary Condition (SIBC)
behave very differently in MQS (Magneto-Quasi-Static) vs full-wave regimes. The key
parameter is `Z_s / (jw * mu0 * R)` where R is the conductor characteristic size.

## SOLVED: Scalar Potential BIE + SIBC (<0.1% error for ALL Z_s)

The scalar potential BIE uses phi (H1 on surface) instead of J (HDivSurface).
SIBC enters via the surface Laplacian (stiffness matrix K).

**System**: `(1/2*M - DL + gamma * SL * M_inv * K) phi = rhs`
where `gamma = Z_s / (jw * mu0)`.

All operators are existing ngsolve.bem: LaplaceSL (scalar H1) + LaplaceDL (scalar H1).
No new C++ code needed.

**Module**: `src/radia/scalar_bie_sibc.py` (ScalarBIE_SIBC class)
**Verification**: `examples/cubit_panels/inductance/scalar_bie_sibc.py`

### Sphere Benchmark (R=10mm, maxh=R/5, 363 DOFs)

| Z_s/(jw*mu0*R) | Analytical | Scalar BIE | Error |
|-----------------|------------|------------|-------|
| 0 (PEC) | 974.62 | 975.47 | +0.09% |
| 0.1 | 908.28 | 908.95 | +0.07% |
| 1.0 | 527.46 | 527.46 | +0.00% |
| 10.0 | 90.83 | 90.76 | -0.07% |

### Usage

```python
from radia.scalar_bie_sibc import ScalarBIE_SIBC

solver = ScalarBIE_SIBC(mesh, order=1)  # assemble once
result = solver.solve(H_inc_cf, Zs, omega)  # solve per frequency
# result['H_rms'], result['P_loss'], result['Q_reactive']

# Frequency sweep (Z_s auto-computed from sigma)
results = solver.frequency_sweep(H_inc_cf, freqs=[1e3, 10e3], sigma=5.8e7)
```

## Formulation Comparison Table

| Formulation | Unknown | MQS Validity | Error | ngsolve.bem? |
|-------------|---------|-------------|-------|-------------|
| **Scalar BIE + SIBC** | **phi (H1)** | **All Z_s** | **<0.1%** | **Yes** |
| EFIE-SIBC | J (HDivSurface) | Z_s/(jw*mu0*R) < 0.1 | 65% at ratio 10 | Yes |
| MFIE tangential | J (HDivSurface) | PEC only | 0% (PEC) | Yes |
| PMCHWT-SIBC | J + M | Impossible in MQS | N/A | N/A |
| FEM-SIBC | H (volume) | All Z_s | ~1% | N/A (pure FEM) |

## Why Scalar BIE Works (and Others Failed)

1. **EFIE-SIBC fails**: SL eigenvalue R/3 for l=1 gives 3x Z_s error
2. **MFIE fails**: RT0 (order=0 HDivSurface) has zero surface curl -> no Z_s dependence
3. **PMCHWT fails**: jw*eps*SL(M) ~ 10^{-14} in MQS -> M doesn't contribute to H
4. **Scalar BIE works**: H1 stiffness K provides nonzero surface Laplacian for SIBC coupling.
   No SL*J term (avoids R/3 eigenvalue). Single-equation formulation (no M coupling needed).

## Decision Tree

```
Workpiece eddy current problem:
  Is Z_s field-dependent (nonlinear mu)?
    Yes -> FEM-SIBC + ESIM Karl iteration (fem_esim_3d.py)
    No  -> Scalar BIE + SIBC (scalar_bie_sibc.py)
           - Surface-only (no volume mesh)
           - Fast frequency sweep (BEM operators cached)
           - All Z_s values (<0.1% error)
```

## Verification Scripts

- `examples/cubit_panels/inductance/scalar_bie_sibc.py` -- Scalar BIE verification (RECOMMENDED)
- `examples/cubit_panels/inductance/verify_sphere_sibc.py` -- Analytical vs BEM vs FEM on sphere
- `examples/cubit_panels/inductance/efie_sibc.py` -- EFIE eigenvalue analysis
- `examples/cubit_panels/inductance/fem_esim_3d.py` -- FEM-SIBC reference implementation
"""


NGBEM_KNOWN_LIMITATIONS = """
# ngsolve.bem Known Limitations and Workarounds (2026-03-29)

## 1. curvaturesafety Required for OCC Surface Meshes

Netgen's default OCC meshing can produce **degenerate elements** on curved surfaces
when `curvaturesafety` is not set. This leads to:
- Zero eigenvalues in BEM operators (LaplaceSL, LaplaceDL)
- Non-invertible matrices / random results from LU solve

**Fix**: Always set `curvaturesafety=1` (or higher) for BEM surface meshes:

```python
# WRONG - may produce degenerate mesh on curved surfaces
mesh = Mesh(geo.GenerateMesh(maxh=maxh))

# CORRECT - curvaturesafety prevents degenerate elements
mesh = Mesh(geo.GenerateMesh(maxh=maxh, curvaturesafety=1))
mesh.Curve(3)
```

**Origin**: Reported by Joachim Schöberl (2026-03-29). Coarse torus mesh produced
topologically incorrect mesh with many zero eigenvalues in SL operator.

## 2. TaskManager Non-Determinism with BEM Operators

`TaskManager()` causes **non-deterministic results** with BEM operator assembly
(LaplaceSL, LaplaceDL) due to floating-point summation order in parallel integration.

| Threads | Fluctuation | Spread (370 nH case) |
|---------|-------------|---------------------|
| 1 | 0 | 0 |
| 2-4 | ~0.3% | ~1 nH |
| 8 | ~1.4% | ~5 nH |

**Workaround**: For reproducible results, do NOT use TaskManager for BEM assembly:

```python
# For reproducibility - no TaskManager
DL_bf = LaplaceDL(u.Trace() * ds) * v.Trace() * ds
SL_bf = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds

# If TaskManager is needed for speed, pin to 1 thread:
from ngsolve import SetNumThreads
with TaskManager():
    SetNumThreads(1)
    SL_bf = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds
```

TaskManager does improve assembly speed for large problems (~5x at 8 cores for
N=5000), but the non-determinism may be unacceptable for validation notebooks.

## 3. ET_QUAD Elements Hang ngsolve.bem (Triangles Only)

ngsolve.bem currently supports **only triangular surface elements** (ET_TRIG).
Passing quad elements (ET_QUAD) to LaplaceSL or LaplaceDL causes the computation
to **hang indefinitely** (no error, no exception).

**Root cause** (ngbem.cpp):
- Line ~152: `ET_QUAD` sets `classnr = -1` (skipped from element grouping)
- Line ~558: `IntegrationRule irtrig(ET_TRIG, intorder)` -- hardcoded to TRIG
- Line ~797: `MappedIntegrationRule<2,3>` -- TRIG reference element only

**Workaround**: Split quads into triangles before BEM assembly.

Minimal reproducer (no Cubit needed):
```python
from netgen.meshing import Mesh as NetgenMesh, MeshPoint, Pnt, Element2D, FaceDescriptor
from ngsolve import Mesh, H1, ds
from ngsolve.bem import LaplaceSL

ngm = NetgenMesh(dim=3)
ngm.Add(FaceDescriptor(surfnr=1, domin=0, domout=0, bc=1))
fd = ngm.Add(FaceDescriptor(surfnr=1, domin=0, domout=0, bc=1))
pts = [ngm.Add(MeshPoint(Pnt(x, y, 0))) for x, y in [(0,0), (1,0), (1,1), (0,1)]]
ngm.Add(Element2D(fd, pts))  # quad element -> will hang
mesh = Mesh(ngm)
fes = H1(mesh, order=1)
u, v = fes.TnT()
sl = LaplaceSL(u.Trace() * ds) * v.Trace() * ds  # HANGS HERE
```

## 4. grad(G) Kernel Not Available (Cross-Mesh Biot-Savart)

ngsolve.bem provides:
- `LaplaceSL`: kernel `G = 1/(4*pi*r)` ✓
- `LaplaceDL`: kernel `dG/dn' = n' . grad' G` ✓ (normal component only)
- `grad G = -(x-x') / (4*pi*|x-x'|^3)` ✗ NOT AVAILABLE

This limits cross-mesh field evaluation. On a **single mesh**, `n_src = n_obs`
so the DL operator provides the needed normal derivative. For **cross-mesh**
evaluation (e.g., coil surface → workpiece surface), the full gradient is needed
because `n_src ≠ n_obs`.

**Use case**: Biot-Savart law for surface currents:
```
H_inc(x_wp) = (1/4pi) * int_{S_coil} J(x') x grad_x G(x, x') dS'
```

**Current workaround**: Direct numerical Biot-Savart with vectorized NumPy,
using per-element J solved by EFIE. This works but lacks ngsolve.bem's
optimized quadrature and singularity handling.

**Status**: Feature request submitted to Joachim Schöberl (2026-03-29).
Cross-mesh LaplaceSL/DL assembly already works (verified). Only the full
gradient kernel is missing.
"""


def get_ngsbem_inductance_documentation(topic: str = "all") -> str:
    """Return ngsolve.bem inductance extraction documentation by topic."""
    topics = {
        "overview": NGBEM_OVERVIEW,
        "api": NGBEM_API,
        "cubit_workflow": NGBEM_CUBIT_WORKFLOW,
        "curve_order": NGBEM_CURVE_ORDER_STUDY,
        "stabilized": NGBEM_STABILIZED,
        "examples": NGBEM_EXAMPLES,
        "best_practices": NGBEM_BEST_PRACTICES,
        "hodge": NGBEM_HODGE_DECOMPOSITION,
        "esim_workpiece": NGBEM_ESIM_WORKPIECE,
        "mqs_formulation_limits": NGBEM_MQS_FORMULATION_LIMITS,
        "known_limitations": NGBEM_KNOWN_LIMITATIONS,
        # Aliases
        "cubit": NGBEM_CUBIT_WORKFLOW,
        "setgeominfo": NGBEM_CUBIT_WORKFLOW,
        "inductance": NGBEM_OVERVIEW,
        "laplace": NGBEM_API,
        "weggler": NGBEM_STABILIZED,
        "harmonic": NGBEM_HODGE_DECOMPOSITION,
        "topology": NGBEM_HODGE_DECOMPOSITION,
        "workpiece": NGBEM_ESIM_WORKPIECE,
        "induction_heating": NGBEM_ESIM_WORKPIECE,
        "efie_sibc": NGBEM_ESIM_WORKPIECE,
        "fem_sibc": NGBEM_ESIM_WORKPIECE,
        "eddy_current": NGBEM_ESIM_WORKPIECE,
        "screening": NGBEM_ESIM_WORKPIECE,
        "mfie": NGBEM_MQS_FORMULATION_LIMITS,
        "pmchwt": NGBEM_MQS_FORMULATION_LIMITS,
        "formulation": NGBEM_MQS_FORMULATION_LIMITS,
        "mqs_sibc": NGBEM_MQS_FORMULATION_LIMITS,
        "scalar_bie": NGBEM_MQS_FORMULATION_LIMITS,
        "scalar_bie_sibc": NGBEM_MQS_FORMULATION_LIMITS,
        "scalar_potential": NGBEM_MQS_FORMULATION_LIMITS,
        "curvaturesafety": NGBEM_KNOWN_LIMITATIONS,
        "taskmanager": NGBEM_KNOWN_LIMITATIONS,
        "quad": NGBEM_KNOWN_LIMITATIONS,
        "grad_g": NGBEM_KNOWN_LIMITATIONS,
        "biot_savart": NGBEM_KNOWN_LIMITATIONS,
        "limitations": NGBEM_KNOWN_LIMITATIONS,
    }

    topic = topic.lower().strip()
    if topic == "all":
        # Return main topics (not aliases)
        main = [
            NGBEM_OVERVIEW, NGBEM_API, NGBEM_CUBIT_WORKFLOW,
            NGBEM_CURVE_ORDER_STUDY, NGBEM_STABILIZED,
            NGBEM_EXAMPLES, NGBEM_BEST_PRACTICES,
            NGBEM_HODGE_DECOMPOSITION, NGBEM_ESIM_WORKPIECE,
            NGBEM_MQS_FORMULATION_LIMITS, NGBEM_KNOWN_LIMITATIONS,
        ]
        return "\n\n".join(main)
    elif topic in topics:
        return topics[topic]
    else:
        available = [k for k in topics if k not in (
            "cubit", "setgeominfo", "inductance", "laplace", "weggler",
            "workpiece", "induction_heating",
        )]
        return f"Unknown topic: '{topic}'. Available: {', '.join(available)}"
