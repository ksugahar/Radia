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
  - ngsolve.bem documentation: https://docu.ngsolve.org/latest/how_to/ngbem.html
  - Lucy Weggler's stabilized BEM: https://github.com/Weggler/docu-ngsbem
  - Netgen PR#232: SetGeomInfo for externally imported meshes
"""

NGBEM_OVERVIEW = """
# ngsolve.bem: Boundary Element Method for Inductance Extraction

## What Is ngsolve.bem

NGSolve's native BEM module providing integral equation operators on surfaces.
Key operators:

| Operator | Kernel | Use |
|----------|--------|-----|
| `LaplaceSL` | 1/(4*pi*r) | MQS inductance (L matrix) |
| `LaplaceDL` | d/dn[1/(4*pi*r)] | Potential problems |
| `HelmholtzSL` | exp(-jkr)/(4*pi*r) | Full-wave BEM |
| `HelmholtzDL` | d/dn[exp(-jkr)/(4*pi*r)] | Full-wave BEM |
| `MaxwellSL` | Full Maxwell kernel | High-frequency |

For inductance extraction at DC to ~1 MHz, **LaplaceSL is sufficient**
(MQS/Darwin regime). No Helmholtz kernel needed.

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
| Speed (large) | H-matrix acceleration | H-matrix (HACApK) |
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
# Or: export_NGSolveCurvedMesh(cubit, surface_only=True)
mesh = Mesh(ngmesh)  # mesh.dim == 2
fes = HDivSurface(mesh, order=0)
# fes.ndof = number of surface edges only (all active in BEM)

# WRONG: dim=3 volume mesh (HDivSurface includes interior edges)
# Created by: export_NGSolveCurvedMesh(cubit) with tet volume mesh
mesh = Mesh(ngmesh)  # mesh.dim == 3
fes = HDivSurface(mesh, order=0)
# fes.ndof = ALL edges (interior + boundary)
# FreeDofs = boundary edges only, but L matrix has zeros for interior
# -> rank-deficient L matrix -> wrong inductance (3000%+ error)
```

**For Cubit meshes**: Use `export_NGSolveCurvedMesh(cubit, surface_only=True)` which
creates a dim=2 surface mesh from Cubit surface blocks (triangles/quads only).

## PITFALL: ds(label) Boundary Name Mismatch

When using `export_NGSolveCurvedMesh(cubit, order=N)`, boundary labels come from
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
after adding all Element2D elements. `export_NGSolveCurvedMesh()` handles this internally.

Without this call, internal topology tables are not built, causing HDivSurface
edge orientation to be inconsistent.

```python
# export_NGSolveCurvedMesh handles this automatically:
mesh = cubit_mesh_export.export_NGSolveCurvedMesh(cubit, order=3, surface_only=True)

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
4. Python: mesh = export_NGSolveCurvedMesh(cubit, order=3)
5. Python: ngsolve.bem LaplaceSL -> inductance extraction
```

No STEP files, no OCC geometry, no SetGeomInfo needed. export_NGSolveCurvedMesh()
handles everything via Cubit's ACIS kernel + CallbackGeometry.

## Step-by-Step Code

```python
import sys, os, math, numpy as np

# Import NGSolve BEFORE Cubit (DLL conflict avoidance)
from ngsolve import (Mesh, Integrate, CF, BND, HDivSurface, TaskManager,
                     GridFunction, ds, sqrt, x, y, z)
from ngsolve.bem import LaplaceSL

cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
import cubit
import cubit_mesh_export

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
cubit.cmd('block 2 add tri all')
cubit.cmd('block 2 name "conductor"')

# --- Step 4: Export with curving (1st order mesh, curved by Netgen) ---
mesh = cubit_mesh_export.export_NGSolveCurvedMesh(
    cubit, order=2, surface_only=True, split_quads=True)

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

`export_NGSolveCurvedMesh(cubit, order=2)` provides quadratic-accurate surfaces
via ACIS CallbackGeometry. order=2 is the recommended default for BEM.

See: `examples/cubit_panels/inductance/demo_curving_effect.py`

## Tri and Quad Surface Meshes

BEM supports both tri (from tet mesh) and quad (from hex mesh) surfaces.
`export_NGSolveCurvedMesh(split_quads=True)` splits quads into triangles
for ngsolve.bem compatibility (HDivSurface requires triangles).

```python
# Tet mesh -> tri surface
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")
cubit.cmd("block 2 add tri all")

# Hex mesh -> quad surface (split to tri for BEM)
cubit.cmd("volume all scheme sweep")
cubit.cmd("mesh volume all")
cubit.cmd("block 2 add face all")  # quad elements

# Export handles both:
mesh = cubit_mesh_export.export_NGSolveCurvedMesh(
    cubit, order=2, surface_only=True, split_quads=True)
```

Register **both** tri and face blocks to support mixed meshes:
```python
cubit.cmd("block 2 add tri all")   # from tet volumes
cubit.cmd("block 2 add face all")  # from hex volumes (Cubit silently skips if none)
```

## Deleted APIs

The old workflow using `export_NetgenMesh()`, `set_*_geominfo()`,
STEP reimport, and OCCGeometry has been completely removed.
Use `export_NGSolveCurvedMesh()` for all Cubit-to-NGSolve mesh transfers.
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

## Measured Results (Cubit mesh -> export_NGSolveCurvedMesh -> LaplaceSL)

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

See: `Radia/examples/cubit/netgen_torus_bem_inductance.py`

This script:
1. Creates a torus in Cubit
2. Tet meshes, defines blocks
3. Uses export_NGSolveCurvedMesh(cubit, order=N) for each order
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

See `Radia/examples/cubit_panels/inductance/test_bem_inductance.py`

Key workflow:
```python
# Cubit: create torus with gap, mesh
cubit.cmd(f"create torus major radius {R} minor radius {a}")
cubit.cmd(f"brick x {3*a} y {gap} z {3*a}")
cubit.cmd(f"move Volume 2 x {R} y 0 z 0 include_merged")
cubit.cmd("subtract volume 2 from volume 1")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")

# Register blocks with source/sink
cubit.cmd('block 1 add volume 1; block 1 name "conductor"')
cubit.cmd('block 2 add tri all;  block 2 name "boundary"')
cubit.cmd('block 3 add tri in surface {source_sid}; block 3 name "source"')
cubit.cmd('block 4 add tri in surface {sink_sid};   block 4 name "sink"')

# Export with curving (1st order mesh, curved by Netgen)
mesh = cubit_mesh_export.export_NGSolveCurvedMesh(
    cubit, order=2, surface_only=True, split_quads=True)

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

3. **Use Curve(3) for curved surfaces**
   - Curve(1) = polygon approximation -> 5-15% inductance error
   - Curve(2) = quadratic -> ~0.5% error
   - Curve(3) = cubic -> ~0.05% error (usually sufficient)

4. **Use export_NGSolveCurvedMesh()** for Cubit meshes
   - Handles curving automatically via ACIS CallbackGeometry
   - No SetGeomInfo, no STEP files, no OCC geometry needed

5. **Laplace kernel for MQS** (DC to ~1 MHz)
   - Do NOT use HelmholtzSL for inductance extraction
   - LaplaceSL is exact in the MQS regime

6. **HDivSurface order=0** (RT0) is sufficient for inductance
   - Higher order (order=1, 2) improves current distribution accuracy
   - But increases DOF count significantly (quadratic growth)
   - For total inductance (scalar), order=0 is usually adequate

7. **Check matrix rank** before solving
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

### Using Old API Instead of export_NGSolveCurvedMesh
```python
# OLD (DELETED): Do not use export_NetgenMesh or set_*_geominfo
# cubit_mesh_export.set_torus_geominfo(ngmesh, ...)

# NEW: Single function call handles everything
mesh = cubit_mesh_export.export_NGSolveCurvedMesh(cubit, order=3)
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

### Format Version Policy

| Direction | Format | Purpose | Tool |
|-----------|--------|---------|------|
| **Input** (-> NGSolve) | **v2.2** | Mesh import | `GmshPostExport.write_v22()` -> `ReadGmsh()` |
| **Output** (NGSolve ->) | **v4.1** | Field visualization | `GmshPostExport.write()` -> GMSH GUI |
| **Output** (NGSolve ->) | **v2.2** | High-order mesh exchange | `GmshPostExport.write_v22()` |

`ReadGmsh()` supports v2.2 only. `GmshPostExport.write_v22()` outputs v2.2 with
arbitrary-order high-order elements (Tri6, Tri10, Tri15, ...).
Element type codes are identical in both versions.

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
   edges). Use `surface_only=True` in `export_NGSolveCurvedMesh()` or OCC `Glue()`.

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


def get_ngbem_inductance_documentation(topic: str = "all") -> str:
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
        # Aliases
        "cubit": NGBEM_CUBIT_WORKFLOW,
        "setgeominfo": NGBEM_CUBIT_WORKFLOW,
        "inductance": NGBEM_OVERVIEW,
        "laplace": NGBEM_API,
        "weggler": NGBEM_STABILIZED,
        "harmonic": NGBEM_HODGE_DECOMPOSITION,
        "topology": NGBEM_HODGE_DECOMPOSITION,
    }

    topic = topic.lower().strip()
    if topic == "all":
        # Return main topics (not aliases)
        main = [
            NGBEM_OVERVIEW, NGBEM_API, NGBEM_CUBIT_WORKFLOW,
            NGBEM_CURVE_ORDER_STUDY, NGBEM_STABILIZED,
            NGBEM_EXAMPLES, NGBEM_BEST_PRACTICES,
            NGBEM_HODGE_DECOMPOSITION,
        ]
        return "\n\n".join(main)
    elif topic in topics:
        return topics[topic]
    else:
        available = [k for k in topics if k not in (
            "cubit", "setgeominfo", "inductance", "laplace", "weggler"
        )]
        return f"Unknown topic: '{topic}'. Available: {', '.join(available)}"
