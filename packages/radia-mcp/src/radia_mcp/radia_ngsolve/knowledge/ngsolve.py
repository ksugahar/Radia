"""
NGSolve usage knowledge base for the Radia MCP server.

Covers: finite element spaces, Maxwell/magnetostatics, solvers,
preconditioners, BEM, mesh generation, nonlinear solving,
electromagnetic formulations (A, Omega, A-Phi, T-Omega),
Kelvin transformation, and common pitfalls gathered from
official docs, community forums, and the EMPY_Analysis project.

Sources:
  - https://docu.ngsolve.org/latest/i-tutorials/
  - https://ngsolve.org/documentation.html (NGS24, Treasure Trove, iFEM)
  - https://forum.ngsolve.org/
  - https://weggler.github.io/ngsbem/intro.html
  - public-safe curated corpus (K. Sugahara's EM formulations)
"""

NGSOLVE_OVERVIEW = """
# NGSolve Overview

NGSolve is a high-performance finite element library with Python bindings.
Key features relevant to electromagnetic simulation:

## Installation

```bash
pip install ngsolve==6.2.2606 netgen-mesher==6.2.2606
# Radia pins both packages together because its native extensions use their C++ ABI.
# NGSolve 6.2.2606 uses ngsolve-openblas; Radia's own HACApK/PARDISO kernels use
# the separately installed MKL 2026 runtime. Do not reuse Radia binaries built
# against NGSolve 6.2.2604.
# The Periodic BC regression in 6.2.2406--6.2.2501 is fixed in this release.
# See: https://forum.ngsolve.org/t/3805
# Compact AMS / COCR ships inside Radia: import radia.sparsesolv_ngsolve as ssn
```

## Basic Workflow

```python
from ngsolve import *
from netgen.occ import *

# 1. Geometry (OCC or CSG)
box = Box(Pnt(0,0,0), Pnt(1,1,1))
box.faces.name = "outer"

# 2. Mesh
mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.3))

# 3. Finite element space
fes = H1(mesh, order=2, dirichlet="outer")
u, v = fes.TnT()

# 4. Bilinear + linear forms
a = BilinearForm(fes)
a += grad(u)*grad(v)*dx

f = LinearForm(fes)
f += 1*v*dx

# 5. Assemble & solve
a.Assemble()
f.Assemble()
gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
```

## Available Direct Solvers

| Solver | Flag | Notes |
|--------|------|-------|
| UMFPACK | `inverse="umfpack"` | Default, single-threaded, robust |
| PARDISO | `inverse="pardiso"` | Multi-threaded (MKL), 4-8x faster |
| Sparse Cholesky | `inverse="sparsecholesky"` | SPD matrices only |
| MUMPS | `inverse="mumps"` | Distributed memory (MPI) |

**Tip**: For production runs, always use PARDISO if available:
```python
gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f.vec
```
"""

NGSOLVE_FE_SPACES = """
# Finite Element Spaces in NGSolve

## De Rham Complex

NGSolve implements the exact de Rham sequence:

```
H1 --grad--> HCurl --curl--> HDiv --div--> L2
```

Each space has continuous and discontinuous variants.

## Space Summary

| Space | Continuity | DOF Location | Canonical Derivative | EM Use |
|-------|-----------|-------------|---------------------|--------|
| `H1` | Point values | Vertices + edges + faces | `grad` | Scalar potential |
| `HCurl` | Tangential | Edges + faces + cells | `curl` | Vector potential A |
| `HDiv` | Normal | Faces + cells | `div` | Magnetic flux B |
| `L2` | None | Cells | - | DG, source terms |
| `HDivSurface` | Normal (surface) | Surface edges | `div` | BEM currents (RT0/RWG) |
| `SurfaceL2` | None (surface) | Surface faces | - | BEM charges |

## HCurl Space (Nedelec Elements)

Used for vector potential A, electric field E.

```python
fes = HCurl(mesh, order=2, dirichlet="outer")
u, v = fes.TnT()

# Curl-curl bilinear form
a = BilinearForm(fes)
a += curl(u)*curl(v)*dx
```

### CRITICAL: nograds Parameter

For magnetostatics (curl-curl only), use `nograds=True` to remove gradient
basis functions and eliminate gauge freedom:

```python
# Magnetostatics: MUST use nograds=True
fes = HCurl(mesh, order=3, dirichlet="outer", nograds=True)

# Time-harmonic Maxwell: do NOT use nograds
fes = HCurl(mesh, order=3, dirichlet="outer", complex=True)
```

**Why**: Without `nograds=True`, the curl-curl matrix is singular (null space =
gradient fields). The regularization term `1e-8*u*v*dx` alone may not suffice.

## HDiv Space (Raviart-Thomas / BDM)

Used for magnetic flux density B.

```python
# BDM (default)
fes = HDiv(mesh, order=2)

# Raviart-Thomas
fes = HDiv(mesh, order=2, RT=True)
```

## HDivSurface (Surface H(div))

Surface version of HDiv for BEM. Corresponds to RWG (Rao-Wilton-Glisson) or
RT0 basis functions on surface triangles.

```python
fes_J = HDivSurface(mesh, order=0)  # RWG edge-based

# CRITICAL: Use .Trace() for BEM operators
j_trial, j_test = fes_J.TnT()
V_op = LaplaceSL(j_trial.Trace() * ds("conductor")) * j_test.Trace() * ds("conductor")
```

## SurfaceL2

Face-based L2 space on surfaces. Used for charge density in PEEC.

```python
fes_L2 = SurfaceL2(mesh, order=0, dual_mapping=True,
                    definedon=mesh.Boundaries("conductor"))
```

**dual_mapping=True**: Correct L2 discretization for dual space in BEM.

## Compound (Product) Spaces

```python
fes = fesH1 * fesL2    # Product space
u, dudn = fes.TrialFunction()
v, dvdn = fes.TestFunction()
```

## DOF Access

```python
# Get DOFs for a mesh entity
edge_dofs = fes.GetDofNrs(NodeId(EDGE, 10))
face_dofs = fes.GetDofNrs(NodeId(FACE, 10))

# Available operators
gf.Operators()  # List available operators for GridFunction
```
"""

NGSOLVE_MAXWELL = """
# Maxwell Equations & Magnetostatics in NGSolve

## Magnetostatics (A-Formulation)

Solves: curl(nu * curl(A)) = J, where nu = 1/(mu_0 * mu_r).

```python
from ngsolve import *
from netgen.occ import *

MU_0 = 4e-7 * 3.14159265358979

# Geometry
core = Box(Pnt(-0.02,-0.02,-0.02), Pnt(0.02,0.02,0.02))
core.solids.name = "core"
core.faces.maxh = 2e-3   # Local mesh refinement on core surface

air = Sphere(Pnt(0,0,0), 0.12) - core
air.solids.name = "air"
geo = Glue([core, air])
geo.faces.name = "outer"

mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=8e-3))

# Material properties
mu_r = mesh.MaterialCF({"core": 1000}, default=1)
nu = 1 / (MU_0 * mu_r)

# HCurl space with nograds=True (magnetostatics!)
fes = HCurl(mesh, order=2, dirichlet="outer", nograds=True)
u, v = fes.TnT()

# Bilinear form with regularization
a = BilinearForm(fes)
a += nu * curl(u) * curl(v) * dx
a += 1e-8 * nu * u * v * dx   # Regularization for well-posedness

# BDDC preconditioner (register BEFORE assembly)
c = Preconditioner(a, "bddc")

# Source term (magnetization or current density)
J_src = mesh.MaterialCF({"core": (0, 0, 1e6)}, default=(0,0,0))
f = LinearForm(fes)
f += J_src * v * dx

with TaskManager():
    a.Assemble()
    f.Assemble()

# Solve with preconditioned CG
gfu = GridFunction(fes)
with TaskManager():
    solvers.CG(sol=gfu.vec, rhs=f.vec, mat=a.mat,
               pre=c.mat, maxsteps=500, printrates=True, tol=1e-10)

# B = curl(A)
B = curl(gfu)
```

## Permanent Magnet Source

```python
# Magnetization as source term: RHS = integral(M * curl(v))
mag = mesh.MaterialCF({"magnet": (1, 0, 0)}, default=(0,0,0))
M_A_per_m = 954930  # Br=1.2T -> M = Br/mu_0

f = LinearForm(fes)
f += M_A_per_m * mag * curl(v) * dx("magnet")
```

## Time-Harmonic Maxwell

```python
# Complex HCurl space (no nograds!)
fes = HCurl(mesh, order=2, dirichlet="outer", complex=True)
u, v = fes.TnT()

omega = 2 * 3.14159 * 1e6  # 1 MHz
sigma = 5.8e7  # Cu conductivity

# curl-curl + j*omega*sigma mass
a = BilinearForm(fes)
a += 1/MU_0 * curl(u) * curl(v) * dx
a += 1j * omega * sigma * u * v * dx("conductor")
a.Assemble()

# Solve with PARDISO (complex direct)
gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f.vec
```

## Nonlinear B-H Curve (Newton Method)

```python
# Define nonlinear nu(|B|^2)
from ngsolve import Norm, sqrt

# Example: simple tanh model
B_sat = 2.0  # Tesla
nu_0 = 1 / MU_0
nu_r_max = 1000

def nu_nonlinear(B2):
    # nu = 1/(mu_0 * mu_r(B))
    # Simple model: mu_r decreases with B
    return nu_0 * (1 + (nu_r_max - 1) / (1 + B2 / B_sat**2))

# In bilinear form:
a = BilinearForm(fes)
a += nu_nonlinear(InnerProduct(curl(u), curl(u))) * curl(u) * curl(v) * dx

# Newton iteration
gfu = GridFunction(fes)
for it in range(20):
    a.Apply(gfu.vec, res)
    a.AssembleLinearization(gfu.vec)
    inv = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
    du = inv * res
    gfu.vec.data -= du
    err = sqrt(abs(InnerProduct(du, res)))
    print(f"Newton iter {it}: err = {err:.2e}")
    if err < 1e-10:
        break
```

## Local Mesh Refinement (OCC)

```python
# Refine only on specific faces (much cheaper than global refinement)
core.faces.maxh = 2e-3        # Fine mesh on core surface
air.faces.maxh = 8e-3         # Coarser air boundary

# This adds only ~6% more elements vs global refinement
```
"""

NGSOLVE_SOLVERS = """
# Solver Selection in NGSolve

## Direct Solvers

```python
# UMFPACK (default, single-threaded)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * f.vec

# PARDISO (multi-threaded, MKL, 4-8x faster)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f.vec

# Sparse Cholesky (SPD only, efficient)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
```

### PARDISO vs UMFPACK

| Feature | PARDISO | UMFPACK |
|---------|---------|---------|
| Threading | Multi-threaded (MKL) | Single-threaded |
| Speed | 4-8x faster | Baseline |
| Complex | Yes | Yes |
| Availability | Requires MKL | Always available |
| Stability | Excellent | Excellent |

**Rule**: Use PARDISO for production, UMFPACK for debugging.

## Iterative Solvers

### CG (Conjugate Gradient) - SPD systems

```python
pre = Preconditioner(a, "bddc")
a.Assemble()

solvers.CG(mat=a.mat, rhs=f.vec, sol=gfu.vec,
           pre=pre.mat, maxsteps=500, tol=1e-10, printrates=True)
```

### GMRes - Non-symmetric / indefinite systems

```python
from ngsolve.solvers import GMRes

result = GMRes(A=system_mat, b=rhs_vec, pre=pre_mat,
               tol=1e-6, maxsteps=200, printrates=True)
```

### MinRes - Symmetric indefinite

```python
solvers.MinRes(mat=a.mat, rhs=f.vec, sol=gfu.vec,
               pre=pre.mat, maxsteps=500, tol=1e-10)
```

## Solver Selection Guide

| System Type | Direct | Iterative |
|-------------|--------|-----------|
| SPD (Poisson, curl-curl) | sparsecholesky/pardiso | CG + BDDC |
| Symmetric indefinite (saddle point) | pardiso/umfpack | MinRes + block precond |
| Non-symmetric (Navier-Stokes) | pardiso/umfpack | GMRes |
| Complex symmetric (eddy current) | pardiso | CG + BDDC (complex) |
| BEM dense | N/A | GMRes + Calderon precond |
| Small (<10K DOF) | pardiso | Not needed |
| Large (>100K DOF) | May run out of memory | BDDC + AMG |

## Eigenvalue Problems

```python
# Generalized eigenvalue: a*u = lambda*m*u
a = BilinearForm(fes)
a += curl(u)*curl(v)*dx
a.Assemble()

m = BilinearForm(fes)
m += u*v*dx
m.Assemble()

# PINVIT (preconditioned inverse iteration)
pre = Preconditioner(a, "bddc")
evals, evecs = solvers.PINVIT(a.mat, m.mat, pre.mat,
                               num=10, maxit=50, printrates=True)
```
"""

NGSOLVE_PRECONDITIONERS = """
# Preconditioners in NGSolve

## Available Preconditioners

| Type | Name | Best For |
|------|------|----------|
| Jacobi | `"local"` | Quick testing, simple problems |
| Block Jacobi | `"local"` + condense | Higher-order H1 |
| Multigrid | `"multigrid"` | h-refined H1, optimal complexity |
| AMG | `"h1amg"` | Large H1 problems, automatic hierarchy |
| BDDC | `"bddc"` | General purpose, parallel, curl-curl |
| Direct | `"direct"` | Coarse level, small problems |

## BDDC Preconditioner (Recommended for EM)

```python
# CRITICAL: Register BEFORE assembly!
c = Preconditioner(a, "bddc")
a.Assemble()

# Use with CG
solvers.CG(mat=a.mat, rhs=f.vec, sol=gfu.vec,
           pre=c.mat, maxsteps=500, tol=1e-10)
```

### BDDC Parameters

```python
# Static condensation for better conditioning
a = BilinearForm(fes, eliminate_internal=True)

# AMG coarse solver for large problems
c = Preconditioner(a, "bddc", coarsetype="h1amg")

# Control wirebasket DOFs
fes = H1(mesh, order=5, wb_withedges=False)  # Only vertex wirebasket
```

### BDDC DOF Classification

| Type | Role |
|------|------|
| LOCAL_DOF | Condensable (interior) |
| WIREBASKET_DOF | Conforming (vertices, first edge DOF) |
| INTERFACE_DOF | Non-conforming (remaining edge/face DOFs) |

### Custom DOF Classification

```python
for ed in mesh.edges:
    dofs = fes.GetDofNrs(ed)
    for d in dofs:
        fes.SetCouplingType(d, COUPLING_TYPE.INTERFACE_DOF)
    # First edge DOF = wirebasket
    fes.SetCouplingType(dofs[0], COUPLING_TYPE.WIREBASKET_DOF)
```

## HCurl AMG Preconditioner (for Large EM Problems)

For large HCurl systems (>100K DOFs), use algebraic multigrid with Hiptmair smoother:

```python
fes = HCurl(mesh, order=0, nograds=True)
u, v = fes.TnT()

a = BilinearForm(fes)
a += 1/mu * curl(u) * curl(v) * dx + 1e-6/mu * u * v * dx

# HCurl AMG with Hiptmair smoother (Reitzinger-Schoeberl)
pre = preconditioners.HCurlAMG(a,
    smoothingsteps=3,           # Gauss-Seidel sweeps per level
    smoothedprolongation=True,  # Improve prolongation operator
    maxcoarse=10,               # Max coarse-level DOFs
    maxlevel=10,                # Max AMG levels
    potentialsmoother='amg',    # AMG on potential (H1) space
    verbose=2)                  # Print level details

with TaskManager():
    a.Assemble()
    inv = CGSolver(a.mat, pre.mat, tol=1e-8, maxiter=100)
    gfu.vec.data = inv * f.vec
```

**Hiptmair smoother**: Decomposes HCurl = grad(H1) + complement;
applies separate AMG smoothers to each part. Standard multigrid fails
on HCurl because the curl kernel (gradient fields) is not handled.

### h1amg for Large H1 Problems

```python
# Use as standalone or as BDDC coarse solver
c = Preconditioner(a, "bddc", coarsetype="h1amg")

# Or standalone:
pre = preconditioners.h1amg(a)
```

**h1amg** uses unsmoothed agglomeration with edge-weight Schur complements.
Works well for scalar H1 problems (Poisson, scalar potential).

## FEM-BEM Preconditioner

For coupled FEM-BEM systems (indefinite):

```python
# Block diagonal preconditioner
bfpre = BilinearForm(grad(u)*grad(v)*dx + 1e-10*u*v*dx +
                     dudn*dvdn*ds("outer")).Assemble()
pre = bfpre.mat.Inverse(freedofs=fes.FreeDofs(), inverse="sparsecholesky")

# Solve with GMRes (indefinite system)
GMRes(A=system, b=rhs, pre=pre, tol=1e-6, maxsteps=200)
```
"""

NGSOLVE_BEM = """
# Boundary Element Method in NGSolve (ngsolve.bem)

## Overview

**ngsolve.bem is bundled with NGSolve** (compiled into ngslib.pyd). It is NOT a separate
package. Do NOT use `import ngbem` -- use `from ngsolve.bem import ...`.

Developer: Lucy Weggler (https://weggler.github.io/ngsbem/intro.html)

```python
# CORRECT: ngsolve.bem is accessed via ngsolve.bem
from ngsolve.bem import LaplaceSL, LaplaceDL, HelmholtzSL, HelmholtzDL

# WRONG: ngsolve.bem is not a standalone package
# import ngbem  # ImportError!
```

## Complete API (ngsolve.bem exports, 31 symbols)

### Laplace Kernel (G = 1/(4*pi*r)) -- MQS/Darwin regime
- `LaplaceSL` -- Single Layer potential (variational form)
- `LaplaceDL` -- Double Layer potential (variational form)
- `SingleLayerPotentialOperator` -- Legacy API, returns IntegralOperator
- `DoubleLayerPotentialOperator` -- Legacy API
- `HypersingularOperator` -- Legacy API

### Helmholtz Kernel (G_k = e^{-jkr}/(4*pi*r)) -- Full wave
- `HelmholtzSL` -- Single Layer (kappa parameter)
- `HelmholtzDL` -- Double Layer (kappa parameter)

### Maxwell Operators -- EFIE/MFIE
- `MaxwellSingleLayerPotentialOperator`
- `MaxwellDoubleLayerPotentialOperator`

### Elasticity
- `LameSL` -- Lame single layer for elasticity BEM

### Evaluation & Utility
- `BiotSavartCF` -- Multipole-accelerated H field from wire currents (returns H in A/m)
- `PotentialOperator` -- Modern variational form wrapper
- `IntegralOperator` -- Legacy operator wrapper (has `.mat`, `.GetPotential()`)

### Two API Styles

```python
# Modern variational API (recommended):
V = LaplaceSL(u*ds("surf")) * v*ds("surf")  # Returns IntegralOperator
V.mat  # Access NGSolve BaseMatrix

# Legacy API:
slp = SingleLayerPotentialOperator(fes, intorder=3)
slp.mat  # Dense matrix
```

## BiotSavartCF -- Wire Coil H-Field

Multipole-accelerated CoefficientFunction for Biot-Savart from wire segments.
Returns H field (A/m). Multiply by mu_0 for B (Tesla).

```python
from ngsolve.bem import BiotSavartCF
from ngsolve.bla import Vec3D

# Create CF with multipole expansion
bs = BiotSavartCF(order=15, kappa=1e-10, center=Vec3D(0,0,0), rad=0.05)
# IMPORTANT: kappa=0 gives NaN. Use kappa=1e-10 for static fields.

# Add wire segments (p1, p2, current, n_quad)
bs.AddCurrent(Vec3D(x1,y1,z1), Vec3D(x2,y2,z2), I, 10)

# Evaluate as CoefficientFunction on any mesh
H = bs(mesh(x, y, z))       # H field [A/m]
B = mu0 * bs(mesh(x, y, z))  # B field [T]
```

Accuracy: +0.24% at r=1.7*R_coil. Multipole expansion requires r_obs > r_source.

Radia wrapper: `from radia.biot_savart import biot_savart_wire, biot_savart_loop`

## MaxwellDL -- Biot-Savart BEM Operator

BEM integral operator for curl(G) kernel. Creates IntegralOperator for EFIE assembly.

```python
from ngsolve.bem import MaxwellDL

fes = HCurl(mesh, order=0)
u, v = fes.TnT()

# Surface source (BND)
io = MaxwellDL(u*ds, kappa) * v * ds

# Curve source (BBND -- edges on boundary)
io = MaxwellDL(u*ds(element_vb=BBND), kappa) * v * ds

# Volume source (dx): NOT SUPPORTED
# MaxwellDL(u*dx, kappa) -> error "need boundary integral"

io.mat           # BEM matrix
io.GetPotential(gf)  # PotentialCF (evaluate on BND via Integrate only)
```

## Available BEM Operators (Summary)

| Operator | Function | Kernel | Source |
|----------|----------|--------|--------|
| Single Layer (V) | `LaplaceSL` / `HelmholtzSL` | G(r) | ds |
| Double Layer (K) | `LaplaceDL` / `HelmholtzDL` | dG/dn | ds |
| Hypersingular (D) | Via curl of SL | d^2G/dn^2 | ds |
| Biot-Savart (BEM) | `MaxwellDL` | curl G(r) | ds, ds(BBND) |
| Biot-Savart (CF) | `BiotSavartCF` | curl G(r) | AddCurrent() |

## Basic BEM Assembly

```python
from ngsolve import *
from ngsolve.bem import LaplaceSL

mesh = Mesh(...)
fes = SurfaceL2(mesh, order=0, definedon=mesh.Boundaries("surf"))
u, v = fes.TnT()

with TaskManager():
    V = LaplaceSL(u*ds("surf")) * v*ds("surf")

# Extract dense matrix
V_mat = V.mat
```

## CRITICAL: .Trace() for HDivSurface BEM

When using LaplaceSL/HelmholtzSL with HDivSurface, you MUST use .Trace():

```python
fes_J = HDivSurface(mesh, order=0)
j_trial, j_test = fes_J.TnT()

# CORRECT - with .Trace()
V_op = LaplaceSL(j_trial.Trace() * ds("conductor")) * j_test.Trace() * ds("conductor")

# WRONG - without .Trace() -> corrupted boundary-edge DOFs
V_op = LaplaceSL(j_trial * ds("conductor")) * j_test * ds("conductor")
```

**Bug symptom**: Without `.Trace()`, boundary-edge DOFs get corrupted diagonal
entries (e.g., -1.1e+11 instead of ~1e-10), causing wildly wrong results.

**Diagnostic**: Check diagonal of assembled matrix:
```python
M = np.array(V_op.mat.IsComplex() and V_op.mat or V_op.mat)
diag = np.real(np.diag(M))
print(f"min diag = {diag.min():.3e}")  # Must be positive
```

## bonus_intorder Parameter

Controls extra quadrature points for near-singular integrals:

```python
# Standard: no bonus
V = LaplaceSL(u*ds("surf")) * v*ds("surf")

# With bonus integration order (more accurate, slower)
V = LaplaceSL(u*ds("surf", bonus_intorder=4)) * v*ds("surf", bonus_intorder=4)
```

| bonus_intorder | Use Case |
|---------------|----------|
| 0 (default) | Far interactions, quick tests |
| 2-3 | Preconditioner assembly |
| 4 | EFIE operator assembly |
| 6+ | High accuracy validation |

## FEM-BEM Coupling Pattern

```python
from ngsolve.bem import LaplaceSL, LaplaceDL

# Interior FEM space
fesH1 = H1(mesh, order=4, dirichlet="top|bottom")

# Boundary BEM space
fesL2 = SurfaceL2(mesh, order=3, dual_mapping=True,
                  definedon=mesh.Boundaries("coupling"))

# Product space
fes = fesH1 * fesL2
u, dudn = fes.TrialFunction()
v, dvdn = fes.TestFunction()

# FEM part
a_fem = BilinearForm(grad(u)*grad(v)*dx).Assemble()

# BEM operators
n = specialcf.normal(3)
with TaskManager():
    V = LaplaceSL(dudn*ds("coupling")) * dvdn*ds("coupling")
    K = LaplaceDL(u*ds("coupling")) * dvdn*ds("coupling")
    D = LaplaceSL(Cross(grad(u).Trace(),n)*ds("coupling")) * \\
        Cross(grad(v).Trace(),n)*ds("coupling")

# Mass matrix for coupling
M = BilinearForm(u*dvdn*ds("coupling")).Assemble()

# Combined system (symmetric indefinite)
system = a_fem.mat + D.mat - (0.5*M.mat + K.mat).T - \\
         (0.5*M.mat + K.mat) - V.mat
```

## Extracting Dense BEM Matrix

```python
import numpy as np

def extract_dense(mat, n):
    \"\"\"Extract dense numpy matrix from NGSolve operator.\"\"\"
    M = np.zeros((n, n), dtype=complex)
    for i in range(n):
        ei = mat.CreateColVector()
        ei[:] = 0
        ei[i] = 1
        col = mat * ei
        for j in range(n):
            M[j, i] = col[j]
    return M
```

## EFIE Calderon Preconditioner (Schoeberl)

For EFIE with HelmholtzSL on HDivSurface, Calderon preconditioning reduces
GMRES iterations from O(N) to O(1):

```python
# EFIE operators
V1 = HelmholtzSL(u.Trace()*ds(bonus_intorder=4), kappa) * v.Trace()*ds(bonus_intorder=4)
V2 = HelmholtzSL(div(u.Trace())*ds(bonus_intorder=4), kappa) * div(v.Trace())*ds(bonus_intorder=4)

lhs = 1j*kappa * V1.mat + 1/(1j*kappa) * V2.mat

# Preconditioner: rotated SL + potential correction
invM = BilinearForm(u.Trace()*v.Trace()*ds).Assemble().mat.Inverse(inverse="sparsecholesky")
n = specialcf.normal(3)
Vrot = HelmholtzSL(Cross(u.Trace(),n)*ds, kappa) * Cross(v.Trace(),n)*ds
Vpot = HelmholtzSL(upot*ds, kappa) * vpot*ds
surfcurl = BilinearForm(Cross(grad(upot).Trace(),n) * v.Trace()*ds).Assemble().mat
invMH1 = BilinearForm(upot*vpot*ds).Assemble().mat.Inverse(inverse="sparsecholesky")

pre = invM @ (kappa*Vrot.mat - 1/kappa*surfcurl@invMH1@Vpot.mat@invMH1@surfcurl.T) @ invM

# Solve with GMRES
gfj.vec[:] = solvers.GMRes(A=lhs, b=rhs.vec, pre=pre, maxsteps=500, tol=1e-8)
```

## BEM Inductance Extraction

Self-inductance of a conductor loop via LaplaceSL on HDivSurface:

```python
import numpy as np
from ngsolve import Mesh, HDivSurface, TaskManager, ds
from ngsolve.bem import LaplaceSL

MU_0 = 4.0 * np.pi * 1e-7

# Surface mesh of conductor (torus, coil, etc.)
mesh = ...  # NGSolve surface mesh
fes = HDivSurface(mesh, order=0)  # RT0 edge DOFs
j_trial, j_test = fes.TnT()

label = list(set(mesh.GetBoundaries()))[0]  # or specific name

# Assemble Laplace SL (MQS kernel: G = 1/(4*pi*r))
with TaskManager():
    L_op = LaplaceSL(j_trial.Trace()*ds(label)) * j_test.Trace()*ds(label)

# Extract dense matrix and scale by mu_0
n = fes.ndof
L_mat = MU_0 * extract_dense(L_op.mat, n)

# Total inductance with uniform excitation (order=0 only!)
e = np.ones(n) / n
L_total = 1.0 / (e @ np.linalg.solve(L_mat, e))
```

**IMPORTANT**: The `e = ones/n` uniform excitation is only valid for order=0 (RT0).
For higher-order HDivSurface (order >= 1), project the current excitation properly.

### Verified accuracy (circular loop, Neumann formula L = mu_0*R*(ln(8R/a) - 2)):
- 89 elements, order=0, Curve(1): +16% error (very coarse mesh)
- For < 5% error: use finer mesh (curvaturesafety >= 1.0) or higher order

## Stabilized BEM for Low-Frequency (Weggler)

Classical EFIE has O(kappa^{-2}) condition number blow-up. Lucy Weggler's stabilized
formulation uses product space HDivSurface x SurfaceL2:

```python
# Stabilized block system: [A_k, Q_k; Q_k^T, k^2*V_k]
# At k=0 (MQS): [A_0, Q_0; Q_0^T, 0] -> saddle point (pure L+R)

from ngsolve import HDivSurface, SurfaceL2, ds
from ngsolve.bem import LaplaceSL

fes_J = HDivSurface(mesh, order=0)  # Current (loop)
fes_rho = SurfaceL2(mesh, order=0)  # Charge (star)

# Block (1,1): A_0 = LaplaceSL on HDivSurface
A_0 = LaplaceSL(j_trial.Trace()*ds) * j_test.Trace()*ds

# Block (1,2) / (2,1): Q_0 = coupling (div J -> rho)
# Built via: BilinearForm(div(u.Trace()) * v_l2 * ds)

# Block (2,2): k^2 * V_0 = scalar single layer (vanishes at k=0)
V_0 = LaplaceSL(rho_trial*ds) * rho_test*ds
```

Reference: https://github.com/Weggler/docu-ngsbem/blob/main/demos/Maxwell_DtN_Stabilized.ipynb
"""

NGSOLVE_BEM_VOL_50 = """
# .vol contract, NGSolve.BEM comparison, and vibro-acoustic teaching lane

## Solver mesh and post-processing boundary

`.vol` is Radia's solver-facing Netgen/NGSolve mesh format, not the public
post-processing format. Do not install a Radia desktop viewer or a Windows
double-click association for `.vol` or `.sol`. Before a solve, run `check-vol`
and retain its `cubit-mesh-export.vol-check.v1` JSON report with the run. The
headless preflight owns element counts, labels, bounding box, orientation, and
source digest.

An NGSolve GridFunction `.sol` contains coefficients but no mesh or finite-
element-space definition. Keep it only as an internal restart/interchange
artifact with its exact companion `.vol` and declared space/order. A public run
that computes a spatial field must also write checked GMSH `.msh v4.1` output
inside the run directory and index it from `result.json`. Human field
inspection uses that GMSH artifact; the solver continues to consume `.vol`.

This split prevents a visualization path from becoming an accidental solver
input and keeps Radia free of a standalone Qt, Tk, PyVista, or VTK desktop
surface. Coreform Cubit's embedded PySide toolbar remains the sole GUI-runtime
exception and owns mesh export, not field analysis.

For the MATLAB acoustic FEM/BEM teaching lane, the same first-order Netgen
`.vol` is read twice: volume tetrahedra form the FEM view and boundary triangles
form the BEM view.  radia-ngsolve/NGSolve.BEM should read the same `.vol` so
node order, boundary labels, and surface orientation become part of the
validation contract.

## 50-case BEM-first milestone and 100-case catalog

Fifty examples are a sensible first milestone before promoting the full
100-case catalog.  Use the BEM-heavy half:

- GYP-031..040: Laplace P1 BEM, dense assembly.
- GYP-041..050: Laplace P1 BEM, H-matrix/compression behavior.
- GYP-051..060: acoustic low-frequency stability checks.
- GYP-061..070: acoustic Helmholtz P1 BEM.
- GYP-091..100: NGSolve.BEM reference and convention checks.

The full 100-case MATLAB acoustic FEM/BEM catalog is:

- GYP-001..010: `.vol` mesh/topology, including negative quad/hex rejection.
- GYP-011..020: H1 P1 tetra FEM.
- GYP-021..030: HCurl/Nedelec edge FEM.
- GYP-031..040: Laplace P1 BEM dense operators.
- GYP-041..050: readable Laplace H-matrix blocks.
- GYP-051..060: acoustic low-frequency kernels.
- GYP-061..070: Helmholtz/acoustic BEM.
- GYP-071..080: scalar FEM/BEM coupling.
- GYP-081..090: RWG/HCurl trace maps.
- GYP-091..100: NGSolve/NGSolve.BEM reference smoke gates.

This is a cross-code validation ladder, not a promise that every case has a new
closed-form solution.  Keep the analytic anchors that remain useful -- sphere
and ball references, Laplace capacity, point-source reproduction, reciprocity,
symmetry, and convergence identities -- but promote the rest as measured
reference artifacts only after versions, dates, timing breakdowns, schema ids,
and convention ids are stored.

## Promotion rule

A result may teach radia-mcp only after four gates pass:

1. `.vol` mesh summary agrees with the intended tri/tet first-order contract.
2. The NGSolve.BEM operator convention gate passes (`.Trace()` where required,
   sign/kappa convention measured rather than assumed).
3. MATLAB/Gypsilab-style result and radia-ngsolve/NGSolve.BEM reference agree
   inside a documented tolerance or show a useful convergence trend.
4. The result manifest is self-contained enough to rerun: tool versions, run
   date, timing breakdown, input digest, output schema, and convention id.

Do not call a cross-code artifact "analytic"; call it a stored reference unless
it is truly backed by a closed-form or manufactured exact identity.

## Vibro-acoustic drum example

The next high-value FEM/BEM example should feel like a struck drum: a baffled
circular membrane or thin elastic plate whose structural FEM mode supplies
normal velocity on the radiating face, while acoustic P1 BEM supplies the
exterior Helmholtz radiation condition.  Gate it by radiation impedance,
acoustic power balance, and far-field directivity, then compare the same
`.vol`/surface labels against NGSolve.BEM.  This is a better teaching example
than another sphere-only smoke test because students can see the chain:
"hit membrane -> structural mode -> radiated pressure -> sound power".
Standing split for drums: the drum structure is FEM, and the air radiation is
acoustic BEM.  Do not turn the drum body into an acoustic volume-FEM problem
unless the lesson is specifically an interior-air/cavity transmission problem.

For the first time-domain visualization, keep it in MATLAB rather than Gmsh:
use a step-force structural modal response plus the causal Rayleigh
retarded-potential integral, then draw r-z pressure snapshots with a MATLAB
figure or an animated GIF written directly from the pressure array.  The GIF
path is a good default artifact for reports because it is browser-viewable,
  headless-friendly, and still backed by readable `.m` code. GMSH can remain
  an export path for later 3D fields, but the teaching visualization
should be readable `.m` code first.

For the bounded-domain FEM visualization, show the outer spherical truncation
as Radia's high-order impedance / absorbing-boundary lane, never as a Kelvin
boundary.  For acoustic waves, Kelvin boundary language is not allowed:
use BEM, Sommerfeld, spherical DtN, or high-order impedance boundary terms.
The teaching picture is: cylindrical drum body + struck top membrane at z=0
+ outgoing pressure wave + full spherical high-order impedance boundary.  The
report GIF should be axis-equal, and the boundary should appear as a full
sphere/circle in the cut plane, not a hemisphere.  The first rung is a
one-sided baffled top-head radiation model.  The real-drum FEM/BEM rung then
adds bottom membrane, cylindrical shell side leakage, explicit damping ratios,
and lower-half radiation; this is the intended demo for explaining FEM/BEM
coupling before moving to a full structural-FEM/acoustic-BEM drum model.  Do not model
this teaching drum as an internal cavity-pressure oscillator unless a true
acoustic cavity mesh is present.  Be precise about the transient: the current
MATLAB drum demo is a reduced FEM ODE integrated by ode45, with damped
membrane/shell oscillators and exterior air pressure reconstructed from causal
retarded boundary potentials.  The BEM layer must not split the observation
field by source direction: top, bottom, and side boundary sources are all
evaluated at every exterior air observation point with the same retarded
Green kernel, then superposed.  A top source may wrap to the side/lower field
and a side source may radiate upward/downward; direction-only painting is a
visualization bug, not FEM/BEM coupling.  Its GIF should color-map only propagating air
pressure outside the drum; the interior air volume is geometric only, not a
cavity pressure DOF.  The same modeling split can be implemented in NGSolve:
membrane/shell FEM or modal damped oscillators provide Neumann data,
ngsolve.bem supplies the exterior radiation operator, and time response comes
from frequency sweeps with inverse FFT or from an externally assembled CQ loop
over Laplace-domain BEM operators.  Be precise about dimensionality: the drum GIF is a
3D axisymmetric physical picture rendered on an r-z slice.  Disk and side
sources are integrated over azimuth, so it is not a 2D Cartesian wave cartoon,
but it is not yet a full 3D structural-FEM/acoustic-BEM drum mesh.

The parallel acoustic-volume teaching lane is not a periodic sine-wave animation:
read the Netgen `.vol` tri/tet mesh, assemble H1/P1 volume FEM plus boundary
P1 BEM, solve the frequency-domain Helmholtz FEM/BEM system over many
frequency bins, multiply by a real pulse spectrum, enforce Hermitian symmetry,
and inverse FFT to obtain a real pressure time history.  In the MATLAB
teaching repo this is the role of `volFemBemIfftResponse`.  It is still not convolution-quadrature TD-BEM, but it is the readable frequency-sweep/iFFT
route for acoustic transmission and interior-fluid examples, not the preferred
drum structural model.

The first convolution-quadrature TD-BEM rung is the `.vol` boundary P1 BEM
path: sample the BDF generating function `delta(zeta)`, evaluate
Laplace-domain retarded single-layer matrices `V(s)` at
`s = delta(zeta)/dt` with `Re(s)>0`, solve `V(s) qhat = ghat`, and
FFT-recover the boundary density and exterior pressure.  In the MATLAB
teaching repo this is `volTdBemConvolutionQuadrature`.  It is real Lubich CQ
TD-BEM for the exterior Dirichlet single-layer problem; coupled volume-FEM CQ
is the next step.

The production-form volume-FEM/interior coupled CQ rung is
`volFemBemCoupledConvolutionQuadrature`: use `.vol` tetrahedra as H1/P1
interior wave FEM, `.vol` boundary triangles as P1 exterior BEM, then solve
`[A+(s/c1)^2 M, -T'Mb; (1/2 Mb-K(s))*T, V(s)] [u_hat; q_hat] = [F_hat; 0]`
at every CQ Laplace point.  A volume source pulse drives the interior and
the BEM flux density radiates through the exterior representation
`-S(s)q + D(s)Tu`.  This is the Calderon/Johnson-Nedelec coupled CQ system
with the retarded double-layer `K(s)`.  Keep the old single-layer row only as
a `SingleLayerTeaching` regression contrast.  `writeVolFemBemCqGmsh3dArtifact`
exports that acoustic-volume CQ result as Gmsh v4.1 NodeData and is useful for
radia-acoustic artifact contracts.  For a drum, reuse the artifact contract but
replace the interior acoustic volume FEM with structural membrane/shell FEM and
use the BEM unknowns for exterior sound.

## Curve-only high-order .vol geometry

Curve-only high-order geometry is useful, but it should be optional and
isolated.  Keep solution unknowns first order (H1 P1 FEM, P1 BEM surface
unknowns) and use high-order `.vol` data only for curved surface geometry,
normals, Jacobians, and quadrature.  A small adapter such as
`CurvedBoundaryView` keeps readability: the default `.vol` parser can continue
to fail loudly on high-order records unless `EnableCurvedGeometry=true` (or the
Python equivalent) is explicitly requested.
"""

NGSOLVE_MESH = """
# Mesh Generation in NGSolve (Netgen + OCC)

## OCC Geometry (Recommended)

```python
from netgen.occ import *

# Basic shapes
box = Box(Pnt(0,0,0), Pnt(1,1,1))
sphere = Sphere(Pnt(0,0,0), 1.0)
cylinder = Cylinder(Pnt(0,0,0), Z, r=0.5, h=2.0)

# Boolean operations
result = box - sphere          # Difference
result = box * sphere          # Intersection
result = box + sphere          # Fusion (single solid, NO internal interface)
result = Glue([box, sphere])   # Glue (keeps internal interface between solids)
```

**CRITICAL: Glue vs Fusion (+)**

| Operation | Internal Interface | Material Regions | Use Case |
|-----------|-------------------|-----------------|----------|
| `Glue([a, b])` | Preserved | Separate materials | Multi-material EM |
| `a + b` (fusion) | Removed | Single material | Uniform body |
| `Compound([a, b])` | Independent | Separate meshes | Components |

For electromagnetic problems with different mu/sigma, ALWAYS use `Glue()`:
```python
# CORRECT: Glue preserves interface for material properties
iron = Box(...); iron.solids.name = "iron"
air = Sphere(...) - iron; air.solids.name = "air"
geo = Glue([iron, air])  # dx("iron") and dx("air") work

# WRONG: Fusion removes the interface
geo = iron + air  # Only one material region!

# Face naming (for boundary conditions)
box.faces.name = "outer"
box.faces.Min(Z).name = "bottom"
box.faces.Max(Z).name = "top"

# Local mesh size
box.faces.maxh = 0.1           # All faces
box.faces.Min(Z).maxh = 0.02   # Only bottom face
box.solids.maxh = 0.05         # Volume constraint

# Generate mesh
mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.3))
```

## 2D Mesh Generation

```python
from netgen.occ import *

# CRITICAL: specify dim=2 for 2D!
rect = Rectangle(3, 5).Face()
mesh = Mesh(OCCGeometry(rect, dim=2).GenerateMesh(maxh=1))

# Without dim=2, you get a 3D mesh of a surface -> wrong!
```

## STEP Import

```python
from netgen.occ import OCCGeometry

geo = OCCGeometry("model.step")
mesh = Mesh(geo.GenerateMesh(maxh=0.5))
```

## Structured Meshes

```python
from ngsolve import *

# 2D structured
mesh = MakeStructured2DMesh(nx=10, ny=10,
                             mapping=lambda x,y: (x, y))

# 3D structured
mesh = MakeStructured3DMesh(nx=5, ny=5, nz=5,
                             mapping=lambda x,y,z: (x, y, z))
```

## Surface Mesh for BEM

```python
from netgen.occ import *

# Surface-only mesh (no volume)
face = Glue(Sphere(Pnt(0,0,0), 1).faces)
mesh = Mesh(OCCGeometry(face).GenerateMesh(maxh=0.2))
mesh.Curve(3)  # Curved elements for accuracy

# Ring mesh for BEM inductance
torus = Torus(R=0.01, r=0.001)   # Not always available
# Alternative: Revolution of circle
from netgen.occ import *
circle = Circle(Pnt(R, 0, 0), r).Face()
ring = circle.Revolve(Axis(Pnt(0,0,0), Z), 360)
ring_mesh = Mesh(OCCGeometry(Glue(ring.faces)).GenerateMesh(maxh=maxh))
```

## Mesh Refinement

```python
# Global h-refinement
mesh.Refine()

# HP refinement (geometric grading toward edges/corners)
mesh.RefineHP(levels=2)

# Curved elements (higher-order geometry approximation)
mesh.Curve(order=3)

# Adaptive refinement based on error estimator
for l in range(max_levels):
    # ... solve ...
    # ... compute error estimator ...
    mesh.Refine()
```

## Mesh Import (.vol)

Load Netgen .vol meshes (from Cubit export or Netgen):

```python
from ngsolve import Mesh

# Load .vol mesh (exported from Cubit or Netgen)
mesh = Mesh("model.vol")

# Block/sideset names become material names and boundary names
print("Materials:", mesh.GetMaterials())
print("Boundaries:", mesh.GetBoundaries())
```

**Key points:**
- .vol files preserve material labels and boundary labels
- Supports high-order curved elements (order 1-5)
- Use `dx("region_name")` and `ds("boundary_name")` in weak forms
- Material properties via `CoefficientFunction([... for mat in mesh.GetMaterials()])`

For Cubit/Coreform exports, run the lightweight parser before solving when
you need a readable FEM/BEM intake gate:

```python
from radia_mcp.radia_ngsolve.netgen_vol import read_netgen_tri_tet_vol

v = read_netgen_tri_tet_vol("model.vol")
print(v.tetrahedron_quality_summary())
print(v.surface_triangle_quality_summary())
```

`surface_triangle_quality_summary()` checks the boundary trace triangles
used by scalar-BEM/RWG operators: area, edge ratio, angle range, and
`2 * inradius / circumradius` quality. This complements Cubit 2026.6's
stronger triangle/tet quality tooling with an independent `.vol`-side
validation after export.

```python
# Material properties from Gmsh physical groups
mu_d = {"iron": 4e-7*3.14*1000, "air": 4e-7*3.14, "coil": 4e-7*3.14}
mu = CoefficientFunction([mu_d.get(mat, 4e-7*3.14) for mat in mesh.GetMaterials()])
```

## Upstream NGSolve `VTKOutput` Reference

`VTKOutput` is a valid NGSolve feature and is shown here as upstream API
documentation. Radia-owned applications, validation artifacts, and examples use
`radia.gmsh_post_export.GmshPostExport` and checked GMSH `.msh v4.1` output.

```python
from ngsolve import VTKOutput

# Basic output (single snapshot)
vtk = VTKOutput(mesh,
                coefs=[B, J, Q],
                names=["B", "J", "Q"],
                filename="result",
                subdivision=1,     # Subdivide for smooth higher-order plots
                legacy=False)      # .vtu (XML) format, smaller files
vtk.Do()
```

### Time Series with .pvd File

Use `Do(time=t)` to create a .pvd collection file that ParaView reads as animation:

```python
# Create VTKOutput ONCE, call Do() each time step
vtk = VTKOutput(mesh,
                coefs=[gfT],
                names=["Temperature"],
                filename="thermal_series",
                subdivision=1, legacy=False)

for step in range(n_steps):
    # ... solve time step ...
    t += dt
    vtk.Do(time=t)  # Writes thermal_series_N.vtu + updates .pvd

# Open thermal_series.pvd in ParaView to play animation
```

### VTKOutput Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `subdivision` | Bisections per element (0=linear, 2=fine) | 0 |
| `legacy` | True=.vtk (legacy), False=.vtu (XML) | True |
| `floatsize` | "single" or "double" precision | "double" |
| `only_element` | Export specific element index | None |

### Do() Parameters

| Parameter | Description |
|-----------|-------------|
| `time` | Associate timestamp for .pvd animation |
| `vb` | VOL (volume, default) or BND (boundary) |
| `drawelems` | BitArray to select specific elements |

**Note**: With `subdivision>0`, the visualized mesh is finer than
the computational mesh. This is cosmetic only.

## Mesh Information

```python
print(f"Elements: {mesh.ne}")
print(f"Vertices: {mesh.nv}")
print(f"Edges: {len(list(mesh.edges))}")
print(f"Faces: {len(list(mesh.faces))}")
print(f"Boundary regions: {mesh.GetBoundaries()}")
print(f"Material regions: {mesh.GetMaterials()}")
```
"""

NGSOLVE_NONLINEAR = """
# Nonlinear Problem Solving in NGSolve

## Newton's Method (Built-in)

```python
from ngsolve.solvers import Newton

# Define nonlinear form
a = BilinearForm(fes)
a += nu(InnerProduct(curl(u),curl(u))) * curl(u) * curl(v) * dx

# Built-in Newton solver
Newton(a, gfu, printing=True, maxerr=1e-10, maxit=20, dampfactor=1.0)
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `maxerr` | 1e-10 | Convergence tolerance |
| `maxit` | 20 | Maximum iterations |
| `dampfactor` | 1.0 | Damping (0<d<=1, use <1 for difficult problems) |
| `printing` | False | Print convergence info |

## Manual Newton Implementation

```python
def newton_solve(gfu, a, fes, tol=1e-10, maxiter=20, damping=1.0):
    res = gfu.vec.CreateVector()
    du = gfu.vec.CreateVector()

    for it in range(maxiter):
        a.Apply(gfu.vec, res)
        a.AssembleLinearization(gfu.vec)

        inv = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
        du.data = inv * res

        gfu.vec.data -= damping * du

        err = sqrt(abs(InnerProduct(du, res)))
        print(f"Newton {it}: err = {err:.2e}")
        if err < tol:
            return it + 1

    raise RuntimeError(f"Newton did not converge after {maxiter} iterations")
```

## Key Technical Points

1. **`a.Apply(vec, res)`**: Evaluates the nonlinear form at `vec`, stores result in `res`
2. **`a.AssembleLinearization(vec)`**: Assembles the Jacobian at `vec` via symbolic differentiation
3. **Symbolic differentiation**: NGSolve automatically differentiates nonlinear CoefficientFunctions
4. **Damping**: Use `dampfactor < 1.0` for strong nonlinearities (e.g., 0.5 for first few iterations)
5. **Convergence**: Well-posed problems show quadratic convergence (error squares each iteration)

## Nonlinear Material Example (B-H Curve)

### Method 1: Analytical Tanh Model (simple)

```python
from ngsolve import sqrt

def nu_BH(B_squared):
    \"\"\"Reluctivity nu = 1/(mu_0 * mu_r(|B|)).\"\"\"
    mu_0 = 4e-7 * 3.14159
    B_mag = sqrt(B_squared + 1e-20)  # Avoid division by zero
    M_sat = 1.6e6  # A/m, saturation magnetization
    chi = M_sat / (B_mag / mu_0 + 1e-10)  # Susceptibility
    mu_r = 1 + chi
    return 1 / (mu_0 * mu_r)
```

### Method 2: BSpline Interpolation (measured data, TEAM benchmark)

```python
from netgen.occ import BSplineCurve1D

# B-H data points (measured)
B_data = [0, 0.1, 0.3, 0.5, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
H_data = [0, 80, 160, 300, 500, 800, 1500, 3000, 8000, 20000, 40000]

# Create BSpline interpolation
nu_spline = BSplineCurve1D(B_data, [h/b if b > 0 else H_data[1]/B_data[1]
                                     for b, h in zip(B_data, H_data)])

# Use in bilinear form (symbolic differentiation works!)
B2 = InnerProduct(curl(u), curl(u))
B_mag = sqrt(B2 + 1e-20)
nu = nu_spline(B_mag)

a = BilinearForm(fes)
a += nu * curl(u) * curl(v) * dx
```

**IMPORTANT**: Do NOT use scipy.interp1d or Python callbacks.
NGSolve CoefficientFunctions require symbolic math (sqrt, **, atan, exp, etc.)
or BSpline for automatic differentiation in Newton's method.

### Method 3: SymbolicEnergy + BSpline.Integrate() (RECOMMENDED)

The most robust approach for measured B-H data. Uses the magnetic co-energy
density w*(B) = integral H(B') dB' as the variational energy functional.
NGSolve automatically differentiates the energy for Newton's Jacobian.

**Key idea**: Instead of assembling `nu(|B|)*curl(u)*curl(v)` (Method 1/2),
define the ENERGY density and let NGSolve compute both residual and Jacobian.

```python
import pandas as pd
from ngsolve import BSpline

mu0 = 4e-7 * 3.14159265

# Load B-H data (column 0 = H in A/m, column 1 = B in Tesla)
df = pd.read_csv('BH.txt', sep='\\t')
H_data = list(df.iloc[:, 0])
B_data = list(df.iloc[:, 1])

# BSpline: H as function of B (= reluctivity * B)
bh_curve = BSpline(2, [0] + B_data, H_data)

# Co-energy density: w*(B) = integral_0^B H(B') dB'
energy_dens = bh_curve.Integrate()

# HCurl FE space
fes = HCurl(mesh, order=2, nograds=True, dirichlet="symmetry|outer")
u, v = fes.TnT()
sol = GridFunction(fes, "A")
sol.vec[:] = 0

# Bilinear form (symmetric for energy-based Newton)
a = BilinearForm(fes, symmetric=True)
# Linear regions (air, coil): standard curl-curl with constant nu
a += SymbolicBFI(1/mu0 * curl(u) * curl(v),
                 definedon=~mesh.Materials("iron"))
# Nonlinear iron: SymbolicEnergy from co-energy density
a += SymbolicEnergy(
    energy_dens(sqrt(1e-12 + curl(u) * curl(u))),
    definedon=mesh.Materials("iron"))
# Coulomb gauging (stabilization)
a += SymbolicBFI(1e-1 * u * v)

c = Preconditioner(a, type="bddc", inverse="sparsecholesky")

# Source: current density in coil
f = LinearForm(-J0 * J_normalized * v * dx("coil"))
```

**BSpline vs BSplineCurve1D**: `BSpline` (from ngsolve) has `.Integrate()` for
exact antiderivative. `BSplineCurve1D` (from netgen.occ) is for geometry only.

**SymbolicEnergy vs SymbolicBFI**: SymbolicEnergy defines a scalar energy density;
NGSolve derives the weak form (1st derivative) and Jacobian (2nd derivative)
automatically via symbolic differentiation. This is more robust than manually
writing `nu(|B|)*curl(u)*curl(v)` because the Jacobian is exact.

## Newton with Energy-Based Line Search (RECOMMENDED)

For nonlinear problems, Newton with energy-based line search guarantees
monotone energy decrease at every iteration. This is much more robust than
fixed damping (`dampfactor`).

**Reference**: `2024_09_21_C_type_gmsh_nonlinear.py`

```python
au = sol.vec.CreateVector()
r = sol.vec.CreateVector()
w = sol.vec.CreateVector()
sol_new = sol.vec.CreateVector()

with TaskManager():
    f.Assemble()
    err = 1
    it = 0

    while err > 1e-4:
        it += 1

        # Current total energy: E = W(A) - <f, A>
        E0 = a.Energy(sol.vec) - InnerProduct(f.vec, sol.vec)

        # Jacobian and residual
        a.AssembleLinearization(sol.vec)
        a.Apply(sol.vec, au)
        r.data = f.vec - au

        # Newton step (CG + BDDC preconditioner)
        inv = CGSolver(mat=a.mat, pre=c.mat)
        w.data = inv * r

        # Convergence check (energy norm of residual)
        err = InnerProduct(w, r)
        print(f"Newton iter {it}: err = {err:.2e}")

        # Energy-based line search (halving)
        sol_new.data = sol.vec + w
        E = a.Energy(sol_new) - InnerProduct(f.vec, sol_new)
        tau = 1
        while E > E0:
            tau = 0.5 * tau
            sol_new.data = sol.vec + tau * w
            E = a.Energy(sol_new) - InnerProduct(f.vec, sol_new)
            print(f"  tau = {tau}, E = {E:.6e}")

        sol.vec.data = sol_new
```

**Why energy line search is superior to fixed damping:**

| Feature | Fixed damping | Energy line search |
|---------|--------------|-------------------|
| Convergence | Sensitive to dampfactor | Guaranteed monotone decrease |
| Tuning | Manual (0.5-1.0) | Automatic |
| Deep saturation | May diverge | Robust |
| Near convergence | Quadratic (full step) | Quadratic (tau=1 accepted) |

**Convergence criterion**: `InnerProduct(w, r)` is the energy norm of the
residual. This is the natural measure for energy-based problems. Typical
threshold: 1e-4 for engineering, 1e-8 for high accuracy.

## Newton with BDDC Preconditioner (Large-Scale)

For large nonlinear problems (>100K DOFs), use BDDC + CG per Newton step:

```python
a = BilinearForm(fes)
a += nu_BH(InnerProduct(curl(u), curl(u))) * curl(u) * curl(v) * dx
c = Preconditioner(a, "bddc")

gfu = GridFunction(fes)
for it in range(20):
    a.Apply(gfu.vec, res)
    a.AssembleLinearization(gfu.vec)
    c.Update()  # Re-update preconditioner with new Jacobian

    inv = CGSolver(a.mat, c.mat, tol=1e-8, maxiter=200)
    du.data = inv * res
    gfu.vec.data -= du

    err = sqrt(abs(InnerProduct(du, res)))
    print(f"Newton {it}: err = {err:.2e}")
    if err < 1e-10:
        break
```

## Implicit Time-Dependent Newton (Backward Euler + Newton)

For time-dependent nonlinear problems (e.g., TEAM 24 motor benchmark):

```python
tau = 0.005  # time step [s]

# Time-dependent eddy current: sigma/tau * (A - A_old) + curl(nu*curl(A)) = J
a = BilinearForm(fes)
a += sigma/tau * u * v * dx("conductor")       # Mass/time term
a += nu_BH(curl(u)**2) * curl(u) * curl(v) * dx  # Nonlinear stiffness
a += -sigma/tau * gfA_old * v * dx("conductor")   # Old time step
a += -J_source * v * dx("coil")                   # Source

# Newton at each time step
for t_step in range(n_steps):
    Newton(a, gfu, maxerr=1e-8, maxit=15, dampfactor=1.0)
    gfA_old.vec.data = gfu.vec
```
"""

NGSOLVE_PITFALLS = """
# Common Pitfalls in NGSolve

## 1. Assignment Operator (CRITICAL)

```python
# WRONG - creates symbolic expression, NOT evaluated!
y = 2 * x_vec

# CORRECT - evaluates and stores result
y.data = 2 * x_vec

# WRONG - nested operations need temporaries
result.data = mat * (u1 + u2)  # FAILS

# CORRECT - separate steps
temp.data = u1 + u2
result.data = mat * temp
```

**Rule**: Always use `.data =` to write to existing vectors.

## 2. Overwriting x, y, z Coordinate Variables

```python
# WRONG - destroys x as coordinate!
for x in range(10):
    pass
# Now x is 9, not the coordinate CoefficientFunction

# CORRECT - use different variable name
for xi in range(10):
    pass
# x, y, z still work as coordinates
```

## 3. Missing dim=2 for 2D OCC Geometry

```python
# WRONG - produces 3D surface mesh
mesh = Mesh(OCCGeometry(rect).GenerateMesh(maxh=1))

# CORRECT - produces 2D mesh
mesh = Mesh(OCCGeometry(rect, dim=2).GenerateMesh(maxh=1))
```

## 4. Preconditioner Registration Order

```python
# WRONG - preconditioner doesn't see element matrices
a.Assemble()
c = Preconditioner(a, "bddc")  # Too late!

# CORRECT - register before assembly
c = Preconditioner(a, "bddc")
a.Assemble()
```

## 5. FreeDofs vs AllDofs

```python
# CORRECT - respects Dirichlet BC
inv = a.mat.Inverse(fes.FreeDofs())

# WRONG - includes constrained DOFs -> singular matrix
inv = a.mat.Inverse()  # Dangerous!
```

## 6. TaskManager Scope

```python
# CORRECT - wrap assembly and solve
with TaskManager():
    a.Assemble()
    f.Assemble()
    solvers.CG(...)

# Also valid - global TaskManager
SetNumThreads(8)
with TaskManager():
    # All parallel work here
    pass
```

## 7. BEM .Trace() Requirement

```python
# WRONG - corrupts boundary DOFs
V = LaplaceSL(j_trial * ds) * j_test * ds

# CORRECT
V = LaplaceSL(j_trial.Trace() * ds) * j_test.Trace() * ds
```

## 8. nograds for Magnetostatics

```python
# WRONG - singular system for curl-curl
fes = HCurl(mesh, order=2, dirichlet="outer")

# CORRECT - removes gradient null space
fes = HCurl(mesh, order=2, dirichlet="outer", nograds=True)
```

## 9. GridFunction.Set() vs Direct Assignment

```python
# CORRECT - projects function onto FE space
gfu.Set(sin(x)*cos(y), definedon=mesh.Materials("domain"))

# CORRECT - set with BND flag for Dirichlet
gfu.Set(x, BND)

# WRONG - BND only sets on Dirichlet-marked boundaries
gfu.Set(x, BND)  # Only sets on dirichlet="..." boundaries!

# To set on ALL boundaries:
gfu.Set(x, mesh.Boundaries(".*"))
```

## 10. Complex Spaces

```python
# For eddy current / time-harmonic: MUST use complex=True
fes = HCurl(mesh, order=2, complex=True, dirichlet="outer")

# WRONG - real space for complex problem
fes = HCurl(mesh, order=2, dirichlet="outer")
a += 1j * omega * sigma * u * v * dx  # Will fail or give wrong results!
```

## 11. Mesh.MaterialCF Spelling

```python
# Material names must match EXACTLY
mu_r = mesh.MaterialCF({"core": 1000}, default=1)

# Check available materials:
print(mesh.GetMaterials())  # ['air', 'core', 'coil']
```

## 12. InnerProduct vs Multiply

```python
# For vector-valued CoefficientFunctions, use InnerProduct:
a += InnerProduct(curl(u), curl(v)) * dx

# Or shorthand (works for simple cases):
a += curl(u) * curl(v) * dx

# For matrix-valued, always use InnerProduct
```

## 13. Upstream NGSolve VTK Export

This is an NGSolve API example. Use `GmshPostExport` for Radia-owned output.

```python
# Export to VTK for ParaView
vtk = VTKOutput(mesh, coefs=[gfu, curl(gfu)],
                names=["A", "B"],
                filename="output",
                subdivision=2)
vtk.Do()
```

## 14. Regularization Term for curl-curl (Forum: Thread 3237)

Pure curl-curl without L2 regularization causes singular matrix errors:

```python
# WRONG - singular matrix with iterative solvers
a += (1/mu_0) * curl(u)*curl(v)*dx

# CORRECT - add small L2 regularization
a += (1/mu_0) * (curl(u)*curl(v)*dx + 1e-6*u*v*dx)
#                                      ^^^^^^^^ L2 term
```

**Note**: `1e-6*curl(u)*curl(v)*dx` is NOT the same as `1e-6*u*v*dx`.
The L2 mass term regularizes the gradient null space; a second curl-curl
term does nothing.

## 15. curl(E) Has No Trace on Edges in HCurl (Forum: Thread 3757)

```python
# WRONG - "don't know how I shall evaluate" error
line_integral = Integrate(curl(E), mesh, definedon=edge)

# CORRECT - interpolate H field first, then integrate
H_gf = GridFunction(HCurl(mesh, order=order))
H_gf.Set(1/(1j*omega*mu) * curl(E))
current = Integrate(H_gf, mesh, definedon=path)
```

**Rule**: curl of an HCurl field is element-wise constant. You cannot
evaluate it on lower-dimensional entities (edges, faces). Always
interpolate into an appropriate space first.

## 16. bonus_intorder for Curved HCurl Elements (Forum: Thread 2068)

```python
# Default integration order is insufficient for curved tets
# Causes O(1e-6) errors

# CORRECT - add bonus_intorder=2 for curved meshes
a += curl(u)*curl(v)*dx(bonus_intorder=2)
f += source*v*dx(bonus_intorder=2)
```

**When needed**: Only for curved tetrahedral meshes (`mesh.Curve(order>=2)`).
Linear geometry (order=1) or boundary layer meshes work fine with defaults.

## 17. BH Curve: Use Analytical Fit, NOT scipy interp (Forum: Thread 2821)

```python
# WRONG - scipy interp doesn't work with NGSolve CoefficientFunctions
from scipy.interpolate import interp1d
f = interp1d(H_data, B_data)
mu = f(H_gf)  # ValueError: object arrays not supported

# CORRECT - use NGSolve symbolic math
B_gf = GridFunction(fes)
mu_gf = GridFunction(fes)
mu_gf.Set(160*sqrt(1-((B_gf/2.5) + 0.17*(B_gf/2.5)**4)**2))
```

**Rule**: NGSolve CoefficientFunctions support `sqrt`, `**`, `+`, `*`,
`atan`, `exp`, `log`, but NOT arbitrary Python callbacks. Fit your
B-H data to an analytical expression.

## 18. GridFunction is NOT a NumPy Array (Forum: Thread 1859)

For p >= 2, `gfu.vec` contains polynomial coefficients, NOT point values.

```python
# WRONG - gives coefficients, not field values
values = np.array(gfu.vec.FV().NumPy())  # Coefficients!

# CORRECT - use IntegrationRuleSpace for point values
fesir = IntegrationRuleSpace(mesh, order=order)
gfuir = GridFunction(fesir)
gfuir.Interpolate(gfu)  # Now contains actual point values
point_values = np.array(gfuir.vec.data[:])
```

## 18b. GridFunction Point Evaluation: Unpack Tuple (Radia)

`mesh.vertices[v].point` returns a tuple `(x, y, z)`. GridFunction `__call__`
does NOT accept a single tuple — it expects separate `x, y, z` arguments
or a `MeshPoint`.

```python
# WRONG - tuple passed as single argument
val = gf(mesh.vertices[v.nr].point)       # TypeError!

# CORRECT - unpack tuple with *
val = gf(*mesh.vertices[v.nr].point)      # gf(x, y, z)

# ALSO CORRECT - via MeshPoint
val = gf(mesh(*mesh.vertices[v.nr].point))
```

This bug broke GMSH |B| export in both `calc_fem_kelvin.py` and
`calc_accel_magnet.py` (fixed 2026-04-09).

## 18c. Surface-Restricted H1 GridFunction (qsurf-Style Save/Load)

Use case: an EM solver computes a surface quantity (q_surf for IH-SIBC,
B_n on a defect surface, ...) and wants to ship it as a stand-alone
.sol file to a downstream solver (thermal, post-processing) without
shipping the mesh inside the .sol.  NGSolve .sol files are raw
coefficient vectors -- no mesh, no fes-order header.

Save side (Producer; e.g. `calc_fem_kelvin.py` IH):

```python
fes_q = H1(mesh, order=fes_order)
gf_q = GridFunction(fes_q)
gf_q.vec[:] = 0                          # interior DOFs stay 0
gf_q.Set(q_surf_cf, definedon=wp_region) # only boundary DOFs touched
gf_q.Save("q.sol")                       # raw coefficient vector
```

`wp_region = mesh.Boundaries("sibc")` is a Region; passing it to
`.Set(definedon=...)` triggers a boundary projection -- only the H1
DOFs that live on that boundary (vertex DOFs on boundary nodes, edge
DOFs on boundary edges, face DOFs on boundary faces; order >= 2) get
populated.  Interior bubbles stay 0.

Load side (Consumer; e.g. `calc_heat.py` Phase B):

```python
em_mesh = Mesh("em.vol")                  # MUST be passed explicitly --
                                          # .sol is mesh-free.
fes_q_em = H1(em_mesh, order=fes_order)   # MUST match save-side order.
gf_q_em = GridFunction(fes_q_em)
gf_q_em.Load("q.sol")
```

THREE contracts that must hold:

1. The EM .vol companion MUST be passed alongside the .sol; the .sol
   has no mesh.  Auto-locating siblings by filename convention is a
   silent-fallback footgun (`<stem>_fem.vol` next to `<stem>_qsurf.sol`
   was the convention; tightened to required 2026-05-20 per
   CLAUDE.md "No Fallbacks").
2. The FES order MUST match the producer's.  No header in .sol means
   a mismatch loads garbage silently (no NGSolve error).  Default
   for radia panel chain is 1 on both sides; for `--fes-order 2` EM
   solves, pass matching `--qsurf-order 2` to the thermal consumer.
3. The producer must `gf.vec[:] = 0` BEFORE `Set(definedon=...)`.
   Otherwise the interior DOFs carry whatever the previous state was;
   downstream consumers expect interior=0.

Point-evaluating the loaded GF at boundary face vertices:

```python
# wp_surf_point is a body-frame surface vertex coord.  em_mesh(...)
# returns a volume MeshPoint inside the element touching the boundary.
em_mip = em_mesh(*wp_surf_point)
val = gf_q_em(em_mip)
```

WHY THIS WORKS WITHOUT `.Trace()`: at a face point on a boundary
element, the H1 volume basis decomposition is
``v_total = sum(boundary_DOFs * boundary_shape_fns)
         + sum(interior_DOFs * interior_bubbles)``.
The interior bubbles are zero at face points by construction (bubble
support is the element interior); the interior DOFs themselves are
also zero (we never wrote them).  So the volume-basis evaluation
EQUALS the trace evaluation -- no `.Trace()` call needed at point
sample.

WHEN YOU **DO** NEED `.Trace()`:

* Building a boundary-only CoefficientFunction for symbolic algebra
  outside `ds(...)`:
  ``q_trace_cf = gf_q_em.Trace()``  -- valid only on boundaries.
* Wiring into HDivSurface BEM operators (LaplaceSL / HelmholtzSL --
  see Section 18 above): ``j_trial.Trace() * ds("conductor")`` etc.
* Coulomb-gauge coupling: ``Cross(grad(u).Trace(), n) * ds("coupling")``.

For LinearForm assembly with surface integrals, `.Trace()` is implicit
in `ds(label)`:

```python
f_form = LinearForm(fes_T)
f_form += q_cf * v * ds(surface_label)   # NGSolve auto-traces both
                                         # q_cf and v's volume bases.
f_form.Assemble()
```

Diagnostic value of POINT-EVALUATION over `.Trace()`: when wp surface
vertices fall slightly OUTSIDE the EM mesh (mesh-mismatch, geometry
drift), `em_mesh(x,y,z)` raises and the consumer can count failures
and report `Q_SURF:projected N/M surface vertices (M-N outside EM
mesh, set to 0)`.  Using `gf.Trace()` would not give this diagnostic
without extra plumbing.

## 19. Python Krylov Solvers for Nested BlockMatrix (Forum: Thread 3399)

```python
# C++ GMRes cannot handle nested BlockMatrix
# Use Python implementation instead:
from ngsolve.krylovspace import CG, GMRes

lhs = BlockMatrix([[A_mat, -B_mat, None],
                   [C_mat, None, -D_mat],
                   [None, M_mat, V_mat]])

sol = GMRes(A=lhs, b=rhs, pre=pre, tol=1e-8, maxsteps=400)
```

## 20. SurfaceL2 for HDiv Normal Matching (Forum: Thread 3379)

Use SurfaceL2 as Lagrange multiplier space for matching HDiv normals:

```python
# SurfaceL2 is the proper space for HDiv normal-trace continuity
fes_lambda = SurfaceL2(mesh, order=order,
                       definedon=mesh.Boundaries("interface"))
```

## 21. type1=True for Nedelec Elements (EMPY)

For magnetostatic/eddy current A-formulations, use first-kind Nedelec elements:

```python
# RECOMMENDED for electromagnetic problems
fesA = HCurl(mesh, order=feOrder, type1=True, nograds=True)

# type1=True: tangential continuity, lower DOF count
# Default (type1=False): second-kind Nedelec, more DOFs
```

**When**: All A-based formulations (A-1, A-2, A-Phi, A-Omega). The T-Omega
formulation also uses type1=True for the electric vector potential T.

## 22. Dirichlet DOF Zeroing After Partial Assembly (EMPY)

When assembling RHS with both Dirichlet and Neumann contributions, the
Dirichlet part must be isolated before adding Neumann terms:

```python
# Step 1: Assemble Dirichlet contribution
gfA.Set(Av, BND, mesh.Boundaries(total_boundary))
f = LinearForm(fes)
f += 1/Mu * curl(gfA) * curl(N) * dx(reduced_region)
f.Assemble()

# Step 2: CRITICAL - Zero out boundary DOF contributions
fcut = np.array(f.vec.FV())[fes.FreeDofs()]
np.array(f.vec.FV(), copy=False)[:] = 0
np.array(f.vec.FV(), copy=False)[fes.FreeDofs()] = fcut

# Step 3: Now add Neumann terms
f += Cross(N.Trace(), Hv) * normal * ds(total_boundary)
f.Assemble()
```

**Why**: Without zeroing, Dirichlet DOF entries leak into the RHS,
corrupting the Neumann contribution.

## 23. CloseSurfaces for Quasi-2D Problems (EMPY)

To solve a 2D problem in 3D NGSolve, identify top/bottom faces
and use ZRefine to create a single element layer:

```python
# Geometry: identify top and bottom
iron.faces.Min(Z).Identify(
    iron.faces.Max(Z), "bot-top", type=IdentificationType.CLOSESURFACES)

# Mesh: single z-layer
ngmesh = OCCGeometry(geo).GenerateMesh(maxh=0.1)
ngmesh.ZRefine("bot-top", [])  # Empty list = no intermediate layers
mesh = Mesh(ngmesh).Curve(1)
```

**Pitfall**: Forgetting `ZRefine` creates multiple z-layers, wasting DOFs.
Always verify `mesh.ne` is the expected number.

## 24. A-Omega Interface Coupling is INDEFINITE (EMPY)

The A-Omega mixed formulation produces a saddle-point system:

```python
# The interface coupling term:
a += -(Omega * (curl(N).Trace() * normal) +
       o * (curl(A).Trace() * normal)) * ds(total_boundary)
```

**Consequence**: The system matrix is INDEFINITE. CG will fail or diverge.
Use MinRes, GMRES, or a custom solver. Extract CSR and use scipy:

```python
A_csr = sp.csr_matrix(a.mat.CSR())
x = sp.linalg.spsolve(A_csr[np.ix_(free, free)], f_free)
```

## 25. HtoOmega Surface Projection (EMPY)

Convert a vector H field to scalar potential Omega on a boundary surface:

```python
fes_surf = H1(mesh, order=feOrder, definedon=mesh.Boundaries(boundary))
omega, psi = fes_surf.TnT()

# Surface Laplacian: find Omega such that grad(Omega) ~= H on surface
a = BilinearForm(fes_surf)
a += grad(omega).Trace() * grad(psi).Trace() * ds(boundary)

f = LinearForm(fes_surf)
f += grad(psi).Trace() * H * ds(boundary)
```

**Use case**: Setting Dirichlet BC for Omega from a known H field
in the total/reduced domain decomposition.

## 26. EMPY_Field Version Incompatibility (EMPY)

C++ pybind11 modules compiled against one NGSolve version may break with
newer versions due to internal API changes (CoefficientFunction vtable, etc.).

```
# Symptom: import crashes or method not found errors
ImportError: DLL load failed while importing EMPY_Field

# Fix: rebuild .pyd against the target NGSolve version
# Or: pin NGSolve version in requirements.txt
ngsolve==6.2.2405
```

**Rule**: Always pin NGSolve version for C++ extension compatibility.

## 27. Static Condensation Pattern (i-tutorials unit-1.4)

Static condensation eliminates internal DOFs, solving only interface DOFs:

```python
a = BilinearForm(grad(u)*grad(v)*dx, condense=True).Assemble()

# Solution recovery: 3-step process
invS = a.mat.Inverse(freedofs=fes.FreeDofs(coupling=True))
ext = IdentityMatrix() + a.harmonic_extension
extT = IdentityMatrix() + a.harmonic_extension_trans
invA = ext @ invS @ extT + a.inner_solve
u.vec.data = invA * f.vec
```

| Attribute | Purpose |
|-----------|---------|
| `a.harmonic_extension` | Maps interface solution -> local DOFs |
| `a.harmonic_extension_trans` | Maps local RHS -> interface RHS |
| `a.inner_solve` | Direct solve of local DOFs |

**Key**: Use `fes.FreeDofs(coupling=True)` for interface-only free DOFs.

## 28. VectorH1 is WRONG for EM Fields (iFEM theory)

```python
# WRONG: 3 independent C^0 scalar components
fes = VectorH1(mesh, order=2)  # Over-constrains tangential continuity!

# CORRECT for E, A fields:
fes = HCurl(mesh, order=2)  # Tangential continuity, normal can jump

# CORRECT for B, J fields:
fes = HDiv(mesh, order=2)   # Normal continuity, tangential can jump
```

**Why**: VectorH1 enforces full C^0 continuity on ALL components. Physical
electromagnetic fields only have tangential (HCurl) or normal (HDiv)
continuity across material interfaces. Using VectorH1 introduces spurious
modes and incorrect interface conditions.

## 29. Higher Order Only Helps for Smooth Solutions (iFEM theory)

Convergence rates depend on solution regularity:
- **P1**: H1 error = O(h), L2 error = O(h^2)
- **P2**: H1 error = O(h^2), L2 error = O(h^3)
- **Pk**: H1 error = O(h^k), L2 error = O(h^{k+1})

**BUT**: At material interfaces (mu jumps), regularity is limited to
H^{1+alpha} where alpha < 1. Higher p gives no benefit there.

```python
# For problems with mu jumps: h-refine near interface, not p-refine
core.faces.maxh = 2e-3    # Fine mesh at mu-jump surface
air.solids.maxh = 8e-3    # Coarser elsewhere
```

## 30. HCurl AMG Preconditioner (i-tutorials unit-2.1.5)

For large HCurl systems, use the HCurl AMG preconditioner:

```python
from ngsolve import preconditioners

fes = HCurl(mesh, order=0)
a = BilinearForm(1/mu*curl(u)*curl(v)*dx + 1e-6/mu*u*v*dx)

# HCurl AMG with Hiptmair smoother
pre = preconditioners.HCurlAMG(a, verbose=2, smoothingsteps=3,
    smoothedprolongation=True, maxcoarse=10, maxlevel=10,
    potentialsmoother='amg')

with TaskManager():
    a.Assemble()
    inv = solvers.CGSolver(a.mat, pre.mat, tol=1e-8, maxiter=100)
    gfu.vec.data = inv * f.vec
```

**Hiptmair smoother**: Decomposes HCurl = grad(H1) + complement;
applies separate smoothers to each. Standard multigrid fails on HCurl
because the curl kernel (gradient fields) is not handled properly.

## 31. Maxwell Eigenvalue: Gradient Projection Required (i-tutorials unit-2.4)

The curl-curl eigenvalue problem has spurious zero eigenvalues from gradient
fields. Must project them out before PINVIT/LOBPCG:

```python
# Build gradient projection
gradmat, fesh1 = fes.CreateGradient()
gradmattrans = gradmat.CreateTranspose()
math1 = gradmattrans @ m.mat @ gradmat
math1[0,0] += 1  # Pin one DOF
invh1 = math1.Inverse(inverse="sparsecholesky")

# Projected preconditioner
proj = IdentityMatrix() - gradmat @ invh1 @ gradmattrans @ m.mat
projpre = proj @ pre.mat

# Solve
evals, evecs = solvers.PINVIT(a.mat, m.mat, pre=projpre, num=12, maxit=20)
```

**Why**: Without projection, PINVIT returns many spurious near-zero
eigenvalues (gradient modes) mixed with physical modes.

## 32. Symmetric COO Export Gives Half Matrix (i-tutorials unit-1.6)

```python
# WRONG - only upper triangle when symmetric=True
a = BilinearForm(fes, symmetric=True)
a.Assemble()
rows, cols, vals = a.mat.COO()  # Half the entries!

# CORRECT - reconstruct full matrix
import scipy.sparse as sp
A = sp.csr_matrix((vals, (rows, cols)), shape=(ndof, ndof))
A = A + A.T - sp.diags(A.diagonal())  # Symmetrize
```

**Rule**: When using `a.mat.COO()` with `symmetric=True`, the result is
only the upper triangle. You must explicitly symmetrize.

## 33. Periodic BC: Need >1 Element in Periodic Direction (i-tutorials unit-4.4)

```python
# WRONG - only 1 element in periodic direction -> singular system
mesh = Mesh(geo.GenerateMesh(maxh=0.6))  # maxh too large for thin domain

# CORRECT - ensure at least 2 elements
mesh = Mesh(geo.GenerateMesh(maxh=0.3))  # Small enough for 2+ layers
```

**Rule**: Periodic boundary conditions require at least 2 elements in the
periodic direction. If maxh is larger than half the domain width, the mesh
may have only 1 element and the periodic identification fails silently.

## 34. VectorH1 for Electromagnetic Fields (i-tutorials unit-2.3)

```python
# WRONG for electromagnetic fields
fes = VectorH1(mesh, order=2)  # All components C^0 continuous

# CORRECT
fes = HCurl(mesh, order=2)  # Tangential continuity only (E, A fields)
fes = HDiv(mesh, order=2)   # Normal continuity only (B, J fields)
```

**Why**: VectorH1 enforces full vector continuity across elements. EM fields
only have tangential (HCurl) or normal (HDiv) continuity at material
interfaces. VectorH1 introduces spurious constraint and wrong interface
conditions.

## 35. Surface-Only Mesh Generation (i-tutorials unit-11.1)

For BEM, you need a surface mesh without volume elements:

```python
from netgen.meshing import MeshingStep

# Surface mesh only (no tetrahedra)
mesh = Mesh(geo.GenerateMesh(maxh=0.5,
            perfstepsend=MeshingStep.MESHSURFACE))
mesh.Curve(2)  # Curved for accuracy
```

**Pitfall**: Using `geo.GenerateMesh(maxh=0.5)` without `perfstepsend`
creates a full volume mesh. BEM operators then see both volume and surface
elements, leading to incorrect results or excessive memory usage.

## 36. Compress() for Subdomain FE Spaces (i-tutorials unit-1.5)

```python
# Space restricted to subdomain has many UNUSED_DOF entries
fes1 = H1(mesh, definedon="conductor")
print(fes1.ndof)  # Same as full mesh ndof (many unused)

# CORRECT - use Compress() to remove unused DOFs
fescomp = Compress(fes1)
print(fescomp.ndof)  # Much smaller, only active DOFs
```

**Why**: Without `Compress()`, the stiffness matrix is much larger than
necessary, wasting memory and solver time. Always compress subdomain spaces.

## 37. mstar Sparsity Pattern Must Match (i-tutorials unit-3.1)

For implicit time stepping `(M + dt*K) * u^{n+1} = ...`, both matrices
must have identical sparsity patterns:

```python
# CORRECT: Use symmetric=False on both to ensure matching patterns
a = BilinearForm(fes, symmetric=False)
a += grad(u)*grad(v)*dx
a.Assemble()

m = BilinearForm(fes, symmetric=False)
m += u*v*dx
m.Assemble()

# Combine via AsVector (element-wise addition of nonzero entries)
mstar = m.mat.CreateMatrix()
mstar.AsVector().data = m.mat.AsVector() + dt * a.mat.AsVector()
```

**Pitfall**: If the sparsity patterns differ (e.g., one is symmetric=True),
`mstar.AsVector().data = ...` silently gives wrong results or crashes.

## 38. NumberSpace for Circuit Coupling (NGS24 TEAM benchmark)

For coupling FEM with external circuits (voltage/current source):

```python
# NumberSpace = single DOF representing a scalar (e.g., current I)
N = NumberSpace(mesh)
fes = HCurl(mesh, order=3, complex=True) * N
(A, I), (v, w) = fes.TnT()

# Bilinear form with circuit equation
a = BilinearForm(fes)
a += nu * curl(A) * curl(v) * dx           # EM equation
a += I * J_coil * v * dx("coil")           # Current source
a += A * J_coil * w * dx("coil")           # Flux linkage
a += R * I * w * dx                         # Resistance

# I (current) is solved simultaneously with A (vector potential)
```

**Use case**: Motor/transformer simulation with prescribed voltage.
The coil current is unknown and coupled to the field via flux linkage.

## 39. Glue vs Fusion for Multi-Material Geometry (i-tutorials unit-4.4)

```python
# WRONG: Fusion (+) removes internal interface
geo = iron + air  # Single material, no interface!

# CORRECT: Glue preserves internal interface
geo = Glue([iron, air])  # dx("iron") and dx("air") work

# CORRECT: Use - for subtraction (different from fusion)
air = Sphere(...) - iron  # Air = Sphere minus iron region
```

## 40. Pickling of NGSolve Objects (appendix)

```python
import pickle

# Save mesh + solution
with open("solution.pkl", "wb") as f:
    pickle.dump([mesh, gfu], f)

# Load (shared references preserved: gfu.space.mesh == mesh)
with open("solution.pkl", "rb") as f:
    mesh2, gfu2 = pickle.load(f)

# CoefficientFunction expression trees also supported
func = x * gfu + y
pickle.dump([mesh, func], outfile)
```

**Supported**: Mesh, FESpace, GridFunction, CoefficientFunction expressions.
Shared references (e.g., two GridFunctions on same space) are preserved.

## 41. GridFunction Mutable Reference Trap in CoefficientFunctions (CRITICAL)

`GridFunction` objects are **mutable**. When embedded in a
`CoefficientFunction` expression, the CF stores a **reference**, not a
copy, of the GridFunction's data. Re-assigning `gf.vec.data` between
uses silently changes the value of every CF that references it.

```python
# WRONG — gf_acc is updated in the loop, every prior expression breaks
gf_acc = GridFunction(fes)
gf_acc.vec[:] = 0.0
J_history = []                        # list of CFs we want to keep
for n in range(N):
    ...
    gf_acc.vec.data += R_n * Aphi_n.vec        # MUTATES the SHARED gf_acc
    proj_cf = gf_acc - grad(phi_acc_n)         # CF holds a REFERENCE
    J_n_cf = J_n_cf - sigma * proj_cf / L_n    # captures gf_acc by ref
    J_history.append(J_n_cf)                    # later evaluations see
                                                # the LATEST gf_acc value!

# CORRECT — snapshot to a fresh GridFunction each iteration
for n in range(N):
    ...
    gf_acc.vec.data += R_n * Aphi_n.vec
    snap = GridFunction(fes, name=f"snap_{n}")  # fresh object
    snap.vec.data = gf_acc.vec                  # explicit copy
    proj_cf = snap - grad(phi_acc_n)            # frozen reference
    J_n_cf = J_n_cf - sigma * proj_cf / L_n
    J_history.append(J_n_cf)                    # safe: snap never mutates
```

**Symptom**: in iterative schemes (Kameari accumulation CLN, Picard
non-linear loops, time stepping that reuses prior states), `stage 0`
matches the analytical answer to machine precision, but `stage 1+`
returns wrong-sign or wildly wrong values that look like algorithmic
divergence — when in fact the prior `J_n_cf` expressions are silently
re-evaluating against the *current* (updated) `gf_acc`, not the
`gf_acc` value at the time the CF was built.

**Verified failure mode** (2026-05-10, Cu sphere a=10mm, B0=1T, 3D
Kameari + Kelvin):
- Without snapshot: τ_0 = 693.95 μs (matches Stoll Cauer-I 694.14
  to 0.027 %), τ_1 = −796 μs (sign flip; analytical 154.6 μs),
  stages 2+ diverge by orders of magnitude.
- With snapshot: τ_0 = 693.95 (−0.027 %), τ_1 = 154.46 (−0.094 %),
  τ_2 = 63.45 (−0.97 %), τ_3 = 33.08 (−4.06 %); stage 4 hits the
  FP64 Hankel-Padé wall.

**Diagnostic strategy**: if a 3D iteration matches axisym at stage 0
but breaks at stage 1+, suspect this trap before suspecting H-H
projection, gauge fixing, ORDER, Schmidt orthogonalisation, or
solver tolerance. None of those help if the GridFunction the CFs
refer to keeps mutating beneath them.

**Cross-references**:
- `cln_3d.py` knowledge — Kameari accumulation snapshot pattern
- public-safe curated corpus
  cln_team28_kelvin.py — reference fix in `Apot_acc -> snap_acc`
"""

NGSOLVE_LINALG = """
# Linear Algebra in NGSolve

## Vector Operations

```python
# Create from GridFunction
vu = gfu.vec
help_vec = vu.CreateVector()

# CRITICAL: Use .data for assignment
help_vec.data = 3 * vu              # Scalar multiply
help_vec.data += mat * vu            # Matrix-vector multiply
help_vec.data = 3 * u1 + 4 * u2     # Linear combination

# Inner product
ip = InnerProduct(help_vec, vu)

# Norm
from ngsolve import Norm
n = Norm(vu)
```

## NumPy Interop

```python
import numpy as np

# GridFunction -> NumPy (shared memory!)
arr = gfu.vec.FV().NumPy()

# Modify in-place
arr[:] = np.sin(np.linspace(0, 1, len(arr)))

# Sparse matrix -> NumPy
rows, cols, vals = a.mat.COO()
import scipy.sparse as sp
A_sp = sp.csr_matrix((vals, (rows, cols)))
```

## Block Matrices

```python
# Create block system from product space
fes = fesA * fesB
a = BilinearForm(fes)
# ... assemble ...

# Access blocks (0-indexed)
A00 = a.mat[0,0]   # Block (0,0)
A01 = a.mat[0,1]   # Block (0,1)

# Block preconditioner
pre = BlockMatrix([[inv_A, None],
                   [None, inv_S]])  # Schur complement
```

## GridFunction Operations

```python
gfu = GridFunction(fes)

# Set from CoefficientFunction
gfu.Set(sin(x)*cos(y))

# Evaluate at point
val = gfu(mesh(0.5, 0.5, 0.5))

# Integrate
integral = Integrate(gfu, mesh)
L2_err = sqrt(Integrate((gfu - exact)**2, mesh))

# Differentiate
grad_gfu = grad(gfu)       # H1
curl_gfu = curl(gfu)       # HCurl
div_gfu = div(gfu)         # HDiv
```
"""


NGSOLVE_EM_FORMULATIONS = """
# Electromagnetic Formulations in NGSolve

Reference implementations from EMPY_Analysis (K. Sugahara).

## Formulation Summary

| Formulation | Primary Unknown | FE Space | Domain | Application |
|-------------|-----------------|----------|--------|-------------|
| **A-1** | A (vector potential) | HCurl | Single | Magnetostatic |
| **A-2** | A | HCurl | Total + Reduced | Magnetostatic |
| **Omega-1** | Omega (scalar potential) | H1 | Single | Magnetostatic |
| **Omega-2** | Omega | H1 | Total + Reduced | Magnetostatic |
| **A-Omega** | A (total) + Omega (reduced) | HCurl * H1 | Hybrid | Magnetostatic |
| **A-Phi** | A + Phi (conductor) | HCurl * H1 | Total + Reduced | Eddy current |
| **T-Omega-1** | T (conductor) + Omega | HCurl * H1 | Single | Eddy current |
| **T-Omega-2** | T + Omega + Loop fields | HCurl * H1 + augmented | Total + Reduced | Eddy current (multiply connected) |

## Total/Reduced Domain Decomposition

Core pattern: split the domain into near-field (total) and far-field (reduced) regions.
The analytical source field (Biot-Savart) is subtracted in the reduced region, so FEM
solves only the material perturbation.

```python
# Geometry: iron + near-field air + far-field air
iron = Box(...).mat("iron")
total_air = Sphere(r=r_total) - iron
total_air.mat("A_domain")
reduced_air = Sphere(r=r_outer) - Sphere(r=r_total)
reduced_air.mat("Omega_domain")
geo = Glue([iron, total_air, reduced_air])

# Boundary naming
iron.faces.Min(X).name = "Bn0"      # B_n = 0 symmetry plane
iron.faces.Min(Z).name = "Ht0"      # H_t = 0 symmetry plane
```

## A-2 Formulation (Vector Potential, Two Regions)

```python
fesA = HCurl(mesh, order=feOrder, type1=True, nograds=True,
             dirichlet=Bn0_boundary)

# Bilinear form
a = BilinearForm(fesA)
a += 1/Mu * curl(A) * curl(N) * dx(total_region)
a += 1/Mu * curl(A) * curl(N) * dx(reduced_region)

# Dirichlet: set known source field on total-reduced boundary
gfA.Set(Av, BND, mesh.Boundaries(total_boundary))

# RHS: Dirichlet contribution + Neumann surface integral
f = LinearForm(fesA)
f += 1/Mu * curl(gfA) * curl(N) * dx(reduced_region)
f.Assemble()

# Isolate free DOFs (critical: zero out boundary DOFs)
fcut = np.array(f.vec.FV())[fesA.FreeDofs()]
np.array(f.vec.FV(), copy=False)[fesA.FreeDofs()] = fcut

# Add Neumann BC: tangential H continuity
f += Cross(N.Trace(), Hv) * normal * ds(total_boundary)
f.Assemble()
```

**Post-processing**: reconstruct full field from solved + source:
```python
BField = Bt + Br + Bs  # solved_total + solved_reduced + source(Biot-Savart)
```

## A-Omega Mixed Formulation (HCurl + H1 Hybrid)

Couples vector potential A in the total region with scalar potential Omega
in the reduced region via interface surface terms:

```python
fesA = HCurl(mesh, order=feOrder, type1=True, nograds=True,
             definedon=total_region, dirichlet=Bn0_boundary)
fesOmega = H1(mesh, order=feOrder, definedon=reduced_region,
              dirichlet=Ht0_boundary)
fes = fesA * fesOmega
(A, Omega), (N, o) = fes.TnT()

a = BilinearForm(fes)
a += 1/Mu * curl(A) * curl(N) * dx(total_region)
a += -Mu * grad(Omega) * grad(o) * dx(reduced_region)

# Interface coupling (saddle-point structure)
normal = specialcf.normal(mesh.dim)
a += -(Omega * (curl(N).Trace() * normal) +
       o * (curl(A).Trace() * normal)) * ds(total_boundary)

# RHS
f = LinearForm(fes)
f += -o * (Bv * normal) * ds(total_boundary)         # Normal B
f += Cross(N.Trace(), Hs) * normal * ds(total_boundary)  # Tangential H
```

**Note**: The A-Omega system is INDEFINITE (saddle-point). Use GMRES or MinRes,
not CG. Or use a custom solver (e.g., ICCG with augmented system).

## A-Phi Eddy Current Formulation

Product space HCurl * H1 where Phi is defined only in the conductor:

```python
fesA = HCurl(mesh, order=feOrder, type1=True, nograds=True,
             dirichlet=Bn0_boundary, complex=True)
fesPhi = H1(mesh, order=feOrder, definedon=conductive_region, complex=True)
fes = fesA * fesPhi
(A, phi), (N, psi) = fes.TnT()

s = 2j * pi * freq  # Complex frequency

a = BilinearForm(fes)
a += 1/Mu * curl(A) * curl(N) * dx
a += s * Sigma * (A + grad(phi)) * (N + grad(psi)) * dx(conductive_region)
```

**Current density**: `J = -s * Sigma * (A + grad(phi))`

**Joule loss** (complex):
```python
WJ = Integrate((J.real*J.real + J.imag*J.imag) / Sigma * dx(cond), mesh) / 2
```

## T-Omega Eddy Current Formulation

Electric vector potential T (conductor only) + magnetic scalar potential Omega:

```python
fesT = HCurl(mesh, order=feOrder, nograds=True,
             definedon="conductor", dirichlet="conductorBND", complex=True)
fesOmega = H1(mesh, order=feOrder, complex=True)
fes = fesT * fesOmega
(T, omega), (W, psi) = fes.TnT()

a = BilinearForm(fes)
a += Mu * grad(omega) * grad(psi) * dx("air")                        # air region
a += Mu * (T + grad(omega)) * (W + grad(psi)) * dx("conductor")      # conductor
a += 1/(s*Sigma) * curl(T) * curl(W) * dx("conductor")               # Faraday's law
```

**Physical fields**: `B = mu * (T + grad(Omega))`, `J = curl(T)`

**Impedance extraction**:
```python
R = Integrate(1/sigma * curl(gfT)**2 * dx("conductor"), mesh)
L = Integrate(mu * (gfT + grad(gfOmega) + loopField)**2 * dx, mesh)
Z = R + s * L
```

## Loop Fields for Multiply-Connected Conductors

For conductors with holes (T-Omega-2), topological loop fields must be computed:

```python
# 1. Compute surface genus from Euler characteristic
g = p - (nv - ne + nf) / 2  # genus = number of independent loops

# 2. For each loop: seed an edge DOF, solve curl-curl, project out gradients
for k in range(g):
    edge_dofs = fes.GetDofNrs(random_boundary_edge)
    gfu.vec[edge_dofs[0]] = 1
    fes.FreeDofs().__setitem__(edge_dofs[0], False)  # Fix DOF

    # Solve curl-curl in air region
    a = BilinearForm(curl(u) * curl(v) * dx("air"))

    # Helmholtz decomposition: subtract gradient part
    a_grad = BilinearForm(grad(phi) * grad(psi) * dx)
    f_grad = LinearForm(grad(psi) * gfu * dx)
    loop_field = gfu - grad(gfPhi)

    # Gram-Schmidt orthogonalization
    for prev_loop in loops:
        prod = Integrate(loop_field * prev_loop * dx, mesh)
        loop_field -= prod * prev_loop
```

**Circuit coupling**: Loop fields are coupled to the main system via
matrix augmentation (additional rows/columns for loop current amplitudes).

## Kelvin Transformation for Open Boundary

### Type 1: Inversion Kelvin Transform

Maps exterior domain via inversion through a sphere:

```python
# Geometry: identify outer sphere with Kelvin region
external_domain.faces[0].Identify(
    Omega_domain.faces[0], "kelvin", IdentificationType.PERIODIC)
fes = Periodic(fes)

# Modified bilinear form in Kelvin domain
rs = model.rKelvin
xs = 2 * rs  # center of inversion sphere
r = sqrt((x-xs)*(x-xs) + y*y + z*z)
fac = rs*rs / (r*r)

a += 1/(Mu*fac) * curl(A) * curl(N) * dx("Kelvin")       # A formulation
a += Mu * fac * grad(omega) * grad(psi) * dx("Kelvin")    # Omega formulation
```

### Type 2: Radial Scaling Kelvin Transform

Anisotropic permeability in a spherical shell (ra to rb):

```python
r = sqrt(x*x + y*y + z*z)
rvec = CF((x, y, z)) / r

# Radial vs tangential permeability
mut = Mu * (rb-ra) * ra / ((rb-r)*(rb-r))   # tangential
mun = Mu * (rb-ra) * ra / (r*r)             # normal

B = curl(A)
Bn = (B * rvec) * rvec          # normal component
Bt = B - Bn                      # tangential component

# CRITICAL: bonus_intorder=4 for accurate integration of varying coefficients
a += (1/mun * Bn*Nn + 1/mut * Bt*Nt) * dx("Kelvin", bonus_intorder=4)
```

## Ht Surface Regularization

Projects tangential H onto the surface to improve Neumann BC quality:

```python
normal = specialcf.normal(mesh.dim)
fes_reg = H1(mesh, order=feOrder, definedon=mesh.Boundaries(boundary))
u, v = fes_reg.TnT()

a = BilinearForm(fes_reg)
a += Cross(normal, grad(u).Trace()) * Cross(normal, grad(v).Trace()) * ds

f = LinearForm(fes_reg)
f += -H * Cross(normal, grad(v).Trace()) * ds

# Corrected tangential H
Ht_corrected = H + Cross(normal, grad(gfu).Trace())
```

This is a surface Hodge decomposition: finds scalar potential u such that
H + n x grad(u) is exactly tangential on the boundary.

## Flux Measurement Through Internal Faces

Uses gradient trick (divergence theorem) to measure flux without surface normals:

```python
class MeasureFace:
    def CalcFlux(self, field):
        # Create H1 function = 1 on face, 0 elsewhere (from one side)
        fes = H1(mesh, order=1, definedon=self.from_side, dirichlet=self.face)
        gf = GridFunction(fes)
        gf.Set(1, definedon=mesh.Boundaries(self.face))
        int1 = Integrate(field * grad(gf) * dx(self.from_side), mesh)

        # Same from other side
        fes = H1(mesh, order=1, definedon=self.to_side, dirichlet=self.face)
        gf = GridFunction(fes)
        gf.Set(1, definedon=mesh.Boundaries(self.face))
        int2 = Integrate(field * grad(gf) * dx(self.to_side), mesh)

        return (int1 - int2) / 2  # Average from both sides
```

Same technique works for current through cross-sections (impedance calculation).

## A Posteriori Error Estimation (ZZ-type)

Projects computed field onto a lower-order space; the difference is the error:

```python
# H-field error (for curl-curl formulations)
fesH = HCurl(mesh, order=feOrder-1, type1=True, nograds=True, definedon=region)
Hd = GridFunction(fesH)
Hd.Set(1/Mu * BField)
dH = 1/Mu * BField - Hd
err = Mu * dH * dH / 2
eta = Integrate(err, mesh, VOL, element_wise=True)

# B-field error (for scalar potential formulations)
fesBd = HDiv(mesh, order=feOrder-1, definedon=region)
Bd = GridFunction(fesBd)
Bd.Set(BField)
err = 1/Mu * (BField - Bd) * (BField - Bd)
eta = Integrate(err, mesh, VOL, element_wise=True)

# Adaptive refinement
maxerr = max(eta)
for el in mesh.Elements():
    mesh.SetRefinementFlag(el, eta[el.nr] > 0.25 * maxerr)
mesh.Refine()
mesh.Curve(curveOrder)
```

## CSR Matrix Extraction to External Solvers

```python
import scipy.sparse as sp

# Extract sparse matrix
asci = sp.csr_matrix(A.mat.CSR())

# Restrict to free DOFs (remove Dirichlet rows/columns)
free = list(fes.FreeDofs())
Acut = asci[np.ix_(free, free)]

# Extract RHS
fcut = np.array(f.vec.FV())[free]

# Solve with external solver (e.g., scipy, custom ICCG)
from scipy.sparse.linalg import spsolve
ucut = spsolve(Acut, fcut)

# Write solution back
np.array(gfu.vec.FV(), copy=False)[free] = ucut
```

## Quasi-2D with CloseSurfaces Identification

Reduce a 3D problem to 2D by identifying top/bottom faces:

```python
iron.faces.Min(Z).Identify(
    iron.faces.Max(Z), "bot-top", type=IdentificationType.CLOSESURFACES)

# Mesh with single element in z-direction
ngmesh = OCCGeometry(geo).GenerateMesh(maxh=0.1)
ngmesh.ZRefine("bot-top", [])  # No intermediate z-layers
mesh = Mesh(ngmesh).Curve(1)
```

## Key HCurl Flags for EM

```python
# type1=True: First-kind Nedelec elements (tangential continuity)
# nograds=True: Removes gradient DOFs (tree-cotree gauge for curl-curl)
# complex=True: Complex-valued space for eddy current (jw analysis)
# definedon=...: Region restriction (e.g., T only in conductor)

fesA = HCurl(mesh, order=2, type1=True, nograds=True,
             dirichlet="Bn0", complex=True)
```

## Material Properties as CoefficientFunction

```python
# Standard pattern: dictionary -> list ordered by mesh.GetMaterials()
mu_d = {"iron": mu0*mur, "A_domain": mu0, "Omega_domain": mu0,
        "Kelvin": mu0, 'default': mu0}
Mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])

sigma_d = {"conductor": 5.8e7, "air": 0, 'default': 0}
Sigma = CoefficientFunction([sigma_d[mat] for mat in mesh.GetMaterials()])
```

## Boundary Condition Summary

| BC Name | Physical Meaning | Space | Method |
|---------|-----------------|-------|--------|
| **Bn=0** | Normal B = 0 (symmetry) | `HCurl(dirichlet="Bn0")` | Tangential A = 0 |
| **Ht=0** | Tangential H = 0 (symmetry) | `H1(dirichlet="Ht0")` | Omega = 0 |
| **n x H** | Tangential H continuity | Neumann | `Cross(N.Trace(), H) * n * ds` |
| **B . n** | Normal B continuity | Neumann | `(n * B) * psi * ds` |
| **Periodic** | Kelvin transform | `Periodic(fes)` | `Identify(PERIODIC)` |
| **CloseSurfaces** | Quasi-2D | `Identify(CLOSESURFACES)` | `ZRefine([])` |

## Time-Domain Eddy Current (Backward Euler + Newton)

For time-transient eddy current with nonlinear materials (e.g., TEAM 24):

```python
tau = 0.005  # time step [s]
sigma = CoefficientFunction(...)  # conductivity

# FE space: HCurl * NumberSpace (A + circuit current I)
fes = HCurl(mesh, order=3, dirichlet="outer") * NumberSpace(mesh)
(A, I), (v, w) = fes.TnT()

# Nonlinear bilinear form (time-dependent)
a = BilinearForm(fes)
# Eddy current: sigma/tau * (A - A_old)
a += sigma/tau * A * v * dx("conductor")
a += -sigma/tau * gfA_old * v * dx("conductor")
# Nonlinear reluctivity
a += nu_BH(curl(A)**2) * curl(A) * curl(v) * dx
# Circuit coupling: U = R*I + N/S * d(flux)/dt
a += I * J_coil * v * dx("coil")
a += A * J_coil * w * dx("coil")
a += R * I * w * dx
a += -U_applied * w * dx

# Newton solve at each time step
for step in range(n_steps):
    Newton(a, gfAPhi, maxerr=1e-8, maxit=15)
    gfA_old.vec.data = gfAPhi.components[0].vec
```

## Maxwell Stress Tensor for Torque Computation

```python
# Torque on rotor via Maxwell stress tensor integration
# T = integral_S (1/mu_0) * [(B.n)*B - 0.5*|B|^2 * n] x r dS

mu0 = 4e-7 * 3.14159265358979
n = specialcf.normal(mesh.dim)
B = curl(gfA)

# For 2D (z-axis torque):
Bn = B * n  # Normal component
Bt = B - Bn * n  # Tangential component

# Torque density: (1/mu0) * Bn * Bt * r (simplified for axisym)
# In 3D: cross product of force and position vector
r_vec = CF((x, y, 0))  # Position from rotation axis
stress = (1/mu0) * (OuterProduct(B, B) - 0.5 * InnerProduct(B, B) * Id(3))
force = stress * n  # Surface force density
torque_density = Cross(r_vec, force)[2]  # Z-component of torque

torque = Integrate(torque_density * ds("rotor_surface"), mesh)
# Multiply by 2 if using symmetry reduction
```
"""

NGSOLVE_ADAPTIVE = """
# Adaptive Mesh Refinement in NGSolve for EM

## Error Estimation

### ZZ-type Estimator (Equilibrium Residual)

Project the computed field onto a lower-order conforming space.
The difference indicates inter-element discontinuity (error):

```python
# For A formulation: project H = (1/mu)*curl(A) into HCurl(order-1)
fesH = HCurl(mesh, order=feOrder-1, type1=True, nograds=True,
             definedon=total_region)
Hd = GridFunction(fesH, "Hprojected")
H = 1/Mu * curl(gfA)
Hd.Set(H)
dH = H - Hd

# Energy norm error indicator
err = Mu * dH * dH / 2
eta = Integrate(err, mesh, VOL, element_wise=True)

# For complex fields (eddy current):
err = Mu * (dH.real*dH.real + dH.imag*dH.imag) / 4
```

### For Scalar Potential: Project B into HDiv(order-1)

```python
fesBd = HDiv(mesh, order=feOrder-1, definedon=total_region)
Bd = GridFunction(fesBd, "Bprojected")
Bd.Set(BField)
err = 1/Mu * (BField - Bd) * (BField - Bd)
eta = Integrate(err, mesh, VOL, element_wise=True)
```

### Combined Error for Eddy Current (Flux + Current)

```python
# Flux error
eta_B = Integrate(Mu*dH*dH/2, mesh, VOL, element_wise=True)

# Current density error
fesJ = HDiv(mesh, order=feOrder-1, definedon=cond_region, complex=True)
Jd = GridFunction(fesJ)
Jd.Set(JField)
dJ = JField - Jd
eta_J = Integrate(1/(Sigma*omega) * (dJ.real*dJ.real + dJ.imag*dJ.imag)/2,
                  mesh, VOL, element_wise=True)

eta = eta_B + eta_J  # Combined error indicator
```

## Marking Strategy

Fixed-fraction marking (Doerfler-style):

```python
maxerr = max(eta)
errorRatio = 0.25  # Refine elements with error > 25% of max

for el in mesh.Elements():
    mesh.SetRefinementFlag(el, eta[el.nr] > errorRatio * maxerr)

mesh.Refine()
mesh.Curve(curveOrder)  # Re-curve after refinement!
```

## Adaptive Loop Pattern

```python
for level in range(max_levels):
    # 1. Create solver and solve
    solver = A_ReducedOmega_Method(model, coil, feOrder=feOrder)
    solver.Solve()

    # 2. Estimate error
    maxerr, eta = solver.CalcError()
    print(f"Level {level}: ndof={solver.fes.ndof}, maxerr={maxerr:.2e}")

    # 3. Check DOF limit
    if solver.fes.ndof > 1000000:
        break

    # 4. Refine and re-create solver
    solver.Refine(maxerr, eta, errorRatio=0.25)
```

**Key**: The solver object must be re-created after refinement because the
FE space changes. This is different from a simple `mesh.Refine()` loop.

## Kelvin Domain Error Scaling

Error indicators in the Kelvin-transformed region need Jacobian scaling:

```python
# Inside Kelvin region: scale by Jacobian factor
H_kelvin = 1/(Mu*fac) * curl(gfA)  # fac = rs^2/r^2
dH_kelvin = H_kelvin - Hd_kelvin
err_kelvin = Mu * fac * dH_kelvin * dH_kelvin / 2
eta_kelvin = Integrate(err_kelvin, mesh, VOL,
                       definedon=mesh.Materials("Kelvin"), element_wise=True)

eta_total = eta_inner + eta_kelvin
```
"""


NGSOLVE_DARWIN_SIBC = """
# Darwin Approximation and Surface Impedance Method in FEM

Reference: A. Kameari, H. Kaimori (SSIL), S. Hiruma (Hokkaido Univ.),
H. Nagamine (Gifu Univ.), IEEJ SA-26-002 / RM-26-002 (2026).

## Darwin Approximation

Maxwell equations without electromagnetic radiation (displacement current dD/dt removed).

**Full Maxwell** A-Phi energy functional:
```
F_Maxwell = (1/2)*integral[ (1/mu)*|curl(A)|^2 + sigma*|s*A + s*grad(phi)|^2
            + (1/c^2)*|s*A|^2 ] dV
```

**Darwin approximation** removes the `(1/c^2)*|s*A|^2` term:
```
F_Darwin = (1/2)*integral[ (1/mu)*|curl(A)|^2 + sigma*|s*A + s*grad(phi)|^2 ] dV
```

**Key properties**:
- Coulomb gauge (div(A) = 0) is imposed **implicitly** by the Darwin functional
- Electromagnetic radiation is completely removed
- Condition number is very poor (Hiruma, IEEE Trans. Magn., 2024)
- Valid up to frequencies where radiation/resonance effects are negligible

### NGSolve Implementation

```python
from ngsolve import *
from netgen.occ import *

# A-Phi formulation for Darwin / Full Maxwell eddy current
fesA = HCurl(mesh, order=2, type1=True, nograds=True,
             dirichlet="outer", complex=True)
fesPhi = H1(mesh, order=2, definedon="conductor",
            dirichlet="...", complex=True)
fes = fesA * fesPhi
(A, phi), (N, psi) = fes.TnT()

s = 2j * pi * freq
mu0 = 4e-7 * pi

# Bilinear form (common to Darwin and Full Maxwell)
a = BilinearForm(fes)
a += 1/Mu * curl(A) * curl(N) * dx                                    # curl-curl
a += Sigma * s * (A + grad(phi)) * (N + grad(psi)) * dx("conductor")  # eddy current

# For Full Maxwell, add displacement current term:
# a += (1/c**2) * s**2 * A * N * dx   # (omit for Darwin)
```

## ICCG Convergence: Two-Stage Behavior

Darwin approximation solved by ICCG exhibits **two distinct convergence stages**:

1. **Fast initial convergence** to the Full Maxwell solution (residual drops rapidly)
2. **Very slow subsequent convergence** toward the Coulomb gauge solution

**Practical implication**: The Full Maxwell solution (stage 1) is already accurate
for low frequencies. The slow stage 2 convergence to Coulomb gauge is unnecessary
for engineering accuracy. This explains why Darwin ICCG appears to "converge" quickly
but the residual then stalls -- it is the Coulomb gauge convergence that is slow.

**Frequency dependence**:
- Higher frequency -> initial minimum is higher (further from Full Maxwell)
- Higher frequency -> oscillation period in stage 2 decreases
  (Coulomb gauge convergence actually slightly improves)
- At very high frequencies, the Darwin ICCG may fail to converge at all

**CAUTION**: Adding a gauge condition (e.g., penalty on div(A)) and using a direct
solver can make the system **unstable** (Hiptmair et al., IEEE Access, 2022).
Use ICCG without explicit gauge for robustness.

## Surface Impedance Method (SIBC in FEM)

Replaces volumetric conductor meshing with a surface boundary condition.
Avoids the need to resolve skin depth in the mesh.

### When to Use

| Approach | Frequency Range | Limitation |
|----------|----------------|------------|
| **Bulk conductor mesh** | DC - ~10^5 Hz | Mesh cannot resolve skin depth above ~10^5 Hz; ICCG diverges |
| **Surface impedance BC** | ~10^2 Hz and above | Invalid below ~10^2 Hz (skin depth > conductor size) |
| **Overlap region** | 10^2 - 10^5 Hz | Both methods valid; choose based on accuracy needs |

**Key advantage**: ICCG converges normally with SIBC even at high frequencies,
whereas bulk conductor analysis fails (mesh under-resolution + ill-conditioning).

### Surface Impedance Formulation

The conductor surface is replaced by an impedance boundary condition:
```
n x H = Zs * (n x (n x E))   on conductor surface
```

where `Zs` is the surface impedance:
```python
import numpy as np
mu0 = 4e-7 * np.pi
omega = 2 * np.pi * freq
sigma = 5.8e7       # conductivity [S/m]
mu_r = 1.0           # relative permeability

delta = np.sqrt(2 / (omega * mu0 * mu_r * sigma))  # skin depth
Zs = (1 + 1j) / (sigma * delta)                     # surface impedance
```

### NGSolve Implementation (SIBC as Robin BC)

```python
# Instead of volumetric conductor elements, apply surface BC:
# The conductor surface replaces the conductor volume

fesA = HCurl(mesh, order=2, type1=True, nograds=True,
             dirichlet="outer", complex=True)

a = BilinearForm(fesA)
a += 1/mu0 * curl(A) * curl(N) * dx                           # curl-curl in air
a += (1/Zs) * A.Trace() * N.Trace() * ds("conductor_surface") # SIBC Robin BC

# RHS: source current (e.g., coil)
f = LinearForm(fesA)
f += J_source * N * dx("coil")
```

**Notes**:
- The conductor volume is excluded from the mesh (or meshed coarsely as air)
- SIBC appears as a Robin-type boundary condition on the conductor surface
- Zs is complex and frequency-dependent
- For thin conductors (thickness t ~ delta), use slab impedance:
  `Zs_slab = Zs * coth(gamma * t)` where `gamma = (1+j)/delta`

## Extended Darwin Approximation

For higher frequencies where standard Darwin deviates from Full Maxwell:

**Method**: Add auxiliary variable chi to correct the scalar potential Phi.
The corrected Phi distribution matches Full Maxwell exactly.

```
F_extended = F_Darwin + correction_term(chi)
```

**Properties**:
- Phi distribution identical to Full Maxwell solution
- E computed as E = -s*(A + grad(Phi)) (same as standard)
- ICCG convergence identical to Full Maxwell (no Coulomb gauge issue)
- Valid up to ~10^9 Hz (no resonance/radiation in Darwin)
- Above 10^9 Hz, deviates from Full Maxwell (electromagnetic radiation dominates)

## Boundary Condition Dependence in Darwin

**CRITICAL**: Darwin approximation treats longitudinal and transverse E-fields
differently, leading to BC-dependent results.

### Two Types of Uniform E-field Excitation

| BC Type | Physical Meaning | Darwin Behavior |
|---------|-----------------|-----------------|
| **Voltage BC** (Phi on faces) | Charge-origin E-field (longitudinal) | Captures correctly |
| **Current BC** (Az on sides) | Current-origin E-field (transverse) | Captures correctly |

A uniform E-field (harmonic field) is **both** longitudinal AND transverse.
Darwin approximation separates these, so the two BC types give **different results**.
This is physically correct behavior of the Darwin model, not a bug.

**At high frequencies** (>10^9 Hz with extended Darwin): BC dependence becomes
visible as the model diverges from Full Maxwell, which does not have this separation.

## Practical Recommendations

1. **Low frequency (DC - 10^5 Hz)**: Use standard A-Phi with bulk conductor mesh.
   Darwin ICCG converges fast enough (stage 1 = Full Maxwell accuracy).

2. **Medium frequency (10^2 - 10^5 Hz)**: Either bulk mesh or SIBC.
   SIBC preferred if skin depth is much smaller than conductor dimensions.

3. **High frequency (10^5 Hz+)**: Use SIBC (bulk mesh fails).
   Standard Darwin + SIBC with ICCG converges normally.

4. **Very high frequency (10^5 - 10^9 Hz)**: Extended Darwin + SIBC.
   Avoids Coulomb gauge convergence issues entirely.

5. **Above 10^9 Hz**: Darwin approximation invalid.
   Use Full Maxwell (radiation/resonance effects cannot be ignored).

## Nonlinear Surface Impedance (ESIM)

For nonlinear magnetic conductors (steel, iron, ferrite), the linear SIBC
formula Zs = (1+j)/(sigma*delta) is inaccurate because mu_r depends on H.

Radia provides ESIM (Effective Surface Impedance Method) which computes
H-dependent Zs(H, omega) from a 1D cell problem. This Zs can be used as
a Robin BC in **any** NGSolve formulation (A-Phi, T-Omega, Darwin, etc.).

See topic `"esim"` for full ESIM API, cell problem details, and NGSolve
integration code examples.

## IEC/SDT Method for Convergence

IEC (Iterative Edge Correction) / SDT (Source Distribution Technique) can
improve ICCG convergence for Darwin systems. See Kaimori et al. (2024) for
the low-frequency stabilized formulation of Darwin model.

## References

1. H. Kaimori, T. Mifune, A. Kameari, S. Wakao: "Low-frequency stabilized
   formulations of Darwin model in time-domain electromagnetic FEM",
   IEEE Trans. Magn., Vol.60, No.3 (2024)
2. R. Hiptmair, F. Kraemer, J. Ostrowski: "A robust Maxwell formulation
   for all frequencies", IEEE Access, Vol.10 (2022)
3. S. Hiruma, T. Mifune, T. Matsuo: "Estimation of condition number of
   quasi-static Darwin model", IEEE Trans. Magn., Vol.60, No.12 (2024)
4. K. Hollaus, O. Biro: "Derivation of a complex permeability from
   the Preisach model", IEEE Trans. Magn. (2002)
"""


NGSOLVE_ESIM = """
# ESIM: Nonlinear Surface Impedance for FEM

ESIM (Effective Surface Impedance Method) extends SIBC to **nonlinear magnetic
conductors** (B-H curve dependent mu). The 1D cell problem computes Zs(H, omega)
which replaces the linear Zs = (1+j)/(sigma*delta) in the Robin BC.

ESIM is formulation-independent: it works with A-Phi, T-Omega, Darwin,
reduced vector potential, or any other FEM formulation that supports
surface impedance as a Robin BC.

## Physics: 1D Cell Problem

The cell problem solves H penetration into a semi-infinite conductor:

```
rho * d^2H/dz^2 + j*omega*mu(|H|)*H = 0    z in [0, inf)
H(0) = H0    (surface tangential field)
H(inf) = 0   (decay into bulk)
```

- Linear mu: analytical Zs = (1+j)/(sigma*delta)
- Nonlinear mu(H): ESIM solves numerically, returns Zs(H0)
- Complex mu = mu' - j*mu": includes magnetic losses (hysteresis, grain eddy current)

## Radia ESIM API

### ESI Table Generation

```python
from radia import generate_esi_table_from_bh_curve, ESITable

# BH curve for steel (H in A/m, B in Tesla)
bh_curve = [[0, 0], [100, 0.5], [500, 1.2], [2000, 1.5], [50000, 2.0]]
sigma = 5e6        # conductivity [S/m]
freq = 10000       # frequency [Hz]

# Generate ESI table: Zs(H0) for H0 = 1..100000 A/m (log-spaced)
esi_table = generate_esi_table_from_bh_curve(bh_curve, sigma, freq, n_points=50)

# Query impedance at specific surface H
Zs = esi_table.get_impedance(H0=1000)  # complex Zs [Ohm]
P, Q = esi_table.get_power_loss(H0=1000)  # W/m^2, var/m^2

# Save/load for reuse across simulations
esi_table.save('steel_10kHz.esi')
esi_table = ESITable.load('steel_10kHz.esi')
```

### Finite Slab (Thin Conductors)

For conductors where skin depth ~ thickness (transition region):

```python
from radia.esim_cell_problem import ESIMFiniteSlabSolver

solver = ESIMFiniteSlabSolver(
    half_thickness=0.005,  # 5mm half-thickness [m]
    bh_curve=bh_curve, sigma=sigma, frequency=freq
)
result = solver.solve(H0=1000)
Zs = result['Z']        # includes coth(gamma*t) correction
rac_rdc = result['R_ac_over_R_dc']
H_profile = result['H_profile']  # H(z) depth distribution
```

### Complex Permeability

For materials with magnetic losses (hysteresis, grain-boundary eddy currents):

```python
from radia.esim_cell_problem import ESIMCellProblemSolver

# Constant complex mu
solver = ESIMCellProblemSolver(
    complex_mu=(500, 50),  # mu'_r=500, mu"_r=50
    sigma=5e6, frequency=10000
)

# H-dependent complex mu
complex_mu_data = [
    [100, 800, 30],   # [H(A/m), mu'_r, mu"_r]
    [1000, 500, 50],
    [5000, 200, 20],
]
solver = ESIMCellProblemSolver(
    complex_mu=complex_mu_data,
    sigma=5e6, frequency=10000
)
```

## NGSolve Integration: Nonlinear Robin BC

### Applicable Formulations

ESIM provides Zs(H) as a surface impedance. In NGSolve, this appears as a
Robin-type boundary condition on the conductor surface:

| Formulation | Robin BC Term | Unknown |
|-------------|--------------|---------|
| **A-Phi** (eddy current) | `(1/Zs) * A.Trace() * N.Trace() * ds` | A in HCurl |
| **T-Omega** | `Zs * curl(T).Trace() * curl(W).Trace() * ds` | T in HCurl |
| **Darwin A-Phi** | Same as A-Phi (no displacement current) | A in HCurl |
| **Reduced A** (magnetostatic+AC) | `(1/Zs) * A.Trace() * N.Trace() * ds` | A_r in HCurl |

### A-Phi Implementation with ESIM

```python
from ngsolve import *
from radia import generate_esi_table_from_bh_curve

# 1. Pre-compute ESI table
bh_curve = [[0, 0], [100, 0.5], [500, 1.2], [2000, 1.5], [50000, 2.0]]
esi_table = generate_esi_table_from_bh_curve(bh_curve, sigma=5e6, frequency=10e3)

# 2. FEM setup (air mesh only; conductor excluded, surface is boundary)
fesA = HCurl(mesh, order=2, nograds=True, dirichlet="outer", complex=True)
A, N = fesA.TnT()
mu0 = 4e-7 * pi
s = 2j * pi * freq

# 3. Per-surface-element Zs storage
fes_bnd = SurfaceL2(mesh, order=0, definedon=mesh.Boundaries("workpiece"))
Zs_real_gf = GridFunction(fes_bnd)
Zs_imag_gf = GridFunction(fes_bnd)
Zs_init = esi_table.get_impedance(100)  # initial guess
Zs_real_gf.Set(Zs_init.real)
Zs_imag_gf.Set(Zs_init.imag)

# 4. Picard iteration (fixed-point on Zs)
gfA = GridFunction(fesA)

for iteration in range(50):
    Zs_cf = Zs_real_gf + 1j * Zs_imag_gf

    a = BilinearForm(fesA)
    a += 1/mu0 * curl(A) * curl(N) * dx                        # curl-curl in air
    a += (1/Zs_cf) * A.Trace() * N.Trace() * ds("workpiece")   # ESIM Robin BC
    a.Assemble()

    f = LinearForm(fesA)
    f += J_source * N * dx("coil")
    f.Assemble()

    gfA.vec.data = a.mat.Inverse(fesA.FreeDofs()) * f.vec

    # Update Zs from surface |H_t|
    H_surf = 1/mu0 * curl(gfA)
    max_change = 0
    for el in mesh.Elements(BND):
        if mesh.GetBCName(el.index) != "workpiece":
            continue
        mip = mesh(*el_centroid(el))
        H_vec = H_surf(mip)
        H_mag = abs(sqrt(sum(h**2 for h in H_vec)))
        Zs_new = esi_table.get_impedance(max(H_mag, 1.0))

        # Under-relaxation for stability
        alpha = 0.3
        Zs_old_r = Zs_real_gf.vec[el.nr]
        Zs_old_i = Zs_imag_gf.vec[el.nr]
        Zs_real_gf.vec[el.nr] = alpha * Zs_new.real + (1-alpha) * Zs_old_r
        Zs_imag_gf.vec[el.nr] = alpha * Zs_new.imag + (1-alpha) * Zs_old_i
        max_change = max(max_change, abs(Zs_new - complex(Zs_old_r, Zs_old_i)))

    if max_change < 1e-6 * abs(Zs_init):
        print(f"ESIM converged at iteration {iteration}")
        break

# 5. Post-processing: Joule loss from ESIM
total_loss = 0
for el in mesh.Elements(BND):
    if mesh.GetBCName(el.index) != "workpiece":
        continue
    mip = mesh(*el_centroid(el))
    H_mag = abs(sqrt(sum(h**2 for h in H_surf(mip))))
    P_loss, _ = esi_table.get_power_loss(max(H_mag, 1.0))
    total_loss += P_loss * el_area(el)  # integrate P' over surface
```

### Induction Heating Workflow with ESIM

```
1. Gmsh/Cubit mesh (air + coil only, no conductor volume)
2. Generate ESI table from workpiece B-H curve + sigma + freq
3. NGSolve A-Phi solve with ESIM Robin BC (Picard iteration)
4. Post-process: Joule loss from ESI table P'(H) * surface_area
5. Optional: thermal analysis with Joule loss as surface heat source
```

Compared to bulk conductor meshing:
- No boundary layer mesh for skin depth resolution
- No frequency-dependent re-meshing
- Handles nonlinear mu(H) automatically
- Joule loss directly from ESI table (no J field needed)

## ESIM vs Linear SIBC

| Property | Linear SIBC | ESIM |
|----------|-------------|------|
| Material | Constant mu_r | B-H curve (nonlinear) |
| Zs formula | `(1+j)/(sigma*delta)` | Numerical 1D cell problem |
| FEM solver | Single linear solve | Picard iteration |
| Complex mu | No | Yes (mu' - j*mu" for losses) |
| Application | Cu, Al conductors | Steel, iron, ferrite |
| Accuracy at high H | Good (linear) | Captures saturation |

## When to Use ESIM vs Bulk Mesh

| Regime | Method | When |
|--------|--------|------|
| delta >> conductor size | Bulk mesh + nonlinear Newton | Low freq, thin conductor |
| delta ~ element size | Fine boundary layer mesh | Moderate freq, needs J profile |
| delta << element size | **ESIM + Robin BC** | High freq, thick conductor |

**Rule of thumb**: If you need > 3 boundary layers to resolve skin depth
and the material is nonlinear, use ESIM.

## Key Radia Files

| File | Class/Function | Purpose |
|------|---------------|---------|
| `esim_cell_problem.py` | `ESIMCellProblemSolver` | Semi-infinite 1D solver |
| `esim_cell_problem.py` | `ESIMFiniteSlabSolver` | Finite thickness solver |
| `esim_cell_problem.py` | `ESITable` | Zs(H) lookup with interpolation |
| `esim_cell_problem.py` | `generate_esi_table_from_bh_curve()` | Main entry point |
| `esim_coupled_solver.py` | `ESIMCoupledSolver` | PEEC + ESIM coupled solver |
| `esim_workpiece.py` | `ESIMWorkpiece` | 3D workpiece with ESI tables |
"""


NGSOLVE_TREE_COTREE = """
# Tree-Cotree Splitting and Low-Frequency Stabilization

Reference: J.-M. Jin, "The Finite Element Method in Electromagnetics,"
3rd ed., Ch.12.9.3 (Wiley, 2014).

## The Low-Frequency Breakdown Problem

The A formulation curl-curl + mass system:

```
[S + jw*T]{A} = {f}
where S_ij = integral(1/mu * curl(Ni) . curl(Nj) dV)   (curl-curl)
      T_ij = integral(sigma * Ni . Nj dV)               (mass/conductivity)
```

At low frequency (w -> 0) or large time step (dt -> infinity), the mass term
vanishes, and [S] alone is singular because curl(grad(phi)) = 0. The
null space of S consists of all gradient fields. This causes:
- Direct solvers: loss of accuracy due to near-singular matrix
- Iterative solvers: very slow convergence or divergence

## Tree-Cotree Splitting (TCS)

TCS decomposes edge DOFs into **gradient** (tree) and **solenoidal** (cotree) parts:

1. Build a spanning tree connecting all mesh nodes without loops
2. Tree edges -> gradient basis functions (replace with grad(Ni))
3. Cotree edges -> solenoidal basis functions (linearly independent curls)
4. Total DOF count unchanged: N_tree = N_free_nodes

```
E = sum_{free nodes} U_i * grad(Ni) + sum_{cotree edges} E_j * Nj
    |--- irrotational part ---|   |--- solenoidal part ---|
```

After TCS + diagonal scaling, the iteration count is INDEPENDENT of
frequency or time step size (verified from dt=50ps to dt=5000ns).

## NGSolve Implementation

NGSolve's hierarchical basis (order >= 2) already separates gradient
and rotational DOFs, so explicit TCS is only needed for lowest-order (order=0).

```python
# For order >= 1: nograds=True removes gradient DOFs automatically
fes = HCurl(mesh, order=1, nograds=True)
# This is equivalent to TCS for higher-order elements

# For lowest-order (order=0): use type1=True + nograds=True
fes = HCurl(mesh, order=0, type1=True, nograds=True)
```

### The Schoberl-Zaglmayr Hierarchical Basis

NGSolve uses the Schoberl-Zaglmayr high-order HCurl basis which
decomposes into gradient and non-gradient (rotational-like) DOFs at
each polynomial order. This is the generalization of TCS to arbitrary order.

```
HCurl DOFs = {edge DOFs (order 0)} + {face DOFs} + {volume DOFs}
           = {gradient (tree)} + {non-gradient (cotree)}
           = {low-order grad} + {low-order rot} + {high-order grad} + {high-order rot}
```

The `nograds=True` flag removes ALL gradient DOFs from the space,
leaving only solenoidal components. This:
- Eliminates the null space of curl-curl
- Improves conditioning at low frequency
- Is ESSENTIAL for magnetostatic A formulation

### When NOT to use nograds=True

- A-Phi eddy current: gradient DOFs couple with Phi -> keep them
- T-Omega: T space needs nograds=True (curl(T) = J must be solenoidal)
- Time-domain: gradient DOFs needed for dA/dt -> use regularization instead

## Regularization Alternative (When Gradients Must Be Kept)

Add a small mass term to regularize the curl-curl operator:

```python
# Regularization parameter: epsilon * mass term
eps = 1e-6 / mu_0  # Small enough not to affect physics
a += eps * u * v * dx

# Or use the Coulomb gauge penalty:
a += alpha * div(u) * div(v) * dx  # Penalty on div(A) != 0
```

## Field-Circuit Coupling (Hybrid Analysis)

Embed lumped circuit elements (R, L, C) into the FEM system by modifying
edge DOFs along the circuit path C:

```python
# For a resistor R on edge k:
# K_kk += Y * L_k^2   where Y = 1/R and L_k = edge length
# For a capacitor C: modify mass matrix T_kk += C * C_kk
# For an inductor L: modify stiffness S_kk += (1/L) * C_kk
```

### MNA Coupling with External Circuit

For complex circuits, use Modified Nodal Analysis:

```
[K0    0    delta*C ] {E  }   {b_prev}
[0     Y    B^T     ] {V  } = {J_CKT }
[delta*C^T  B    0  ] {J_CP}   {0     }
```

where:
- K0: FEM matrix
- Y: circuit admittance matrix (MNA)
- C: edge-to-port coupling matrix (C_ip = integral Ni . L_hat dl on port p)
- B: port selection matrix
- delta = 1/(2*dt) for Newmark-beta time stepping

The system is symmetric if Y is symmetric, and always stable.

## Key References

- Jin Ch.12.9.3: Tree-cotree splitting for low-frequency stabilization
- Jin Ch.12.9.1-2: Hybrid field-circuit analysis (frequency/time domain)
- Schoberl & Zaglmayr (2005): High-order Nedelec elements with
  explicit gradient-rotational decomposition
- Wang & Jin (2010): TCS for time-domain FEM, 5 orders of magnitude in dt
"""

NGSOLVE_PML = """
# Perfectly Matched Layers (PML) in NGSolve

Reference: J.-M. Jin, "The Finite Element Method in Electromagnetics,"
3rd ed., Ch.9.6 (Wiley, 2014).

## Concept

PML is a coordinate-stretching technique that creates a reflectionless
absorbing layer. Complex-valued stretching factors attenuate waves:

```
s_x = s' - j*s"    (s' >= 1 for evanescent, s" >= 0 for propagating)
```

The wave impedance is NOT affected by stretching -> zero reflection
at the PML interface (regardless of angle, frequency, polarization).

## Anisotropic Absorber Interpretation

PML is equivalent to an anisotropic material with:

```
epsilon_PML = eps * Lambda
mu_PML = mu * Lambda

where Lambda = diag(s_y*s_z/s_x, s_x*s_z/s_y, s_x*s_y/s_z)
```

This interpretation makes FEM implementation straightforward: just
modify the material tensors in the PML region.

## NGSolve PML Implementation

NGSolve has built-in PML support via `mesh.SetPML()`:

```python
from ngsolve import *

# Create geometry with PML region
geo = SplineGeometry()
# ... define inner domain + PML shell ...
mesh = Mesh(geo.GenerateMesh())

# Set PML with polynomial grading
mesh.SetPML(pml.Cartesian(mins=(-1,-1), maxs=(1,1),
                           alpha=2j))  # alpha = absorption strength

# Or radial PML for spherical domains:
mesh.SetPML(pml.Radial(rad=1.0, alpha=1j))

# FEM formulation is unchanged - PML modifies the Jacobian
fes = HCurl(mesh, order=2, complex=True)
u, v = fes.TnT()

a = BilinearForm(fes)
a += 1/mu * curl(u) * curl(v) * dx - omega**2 * eps * u * v * dx
a.Assemble()
```

## PML Parameter Selection

### Polynomial Grading (Essential for Discretization Error)

```
sigma(x) = sigma_max * (x/L)^m    m = 2 (parabolic, recommended)
```

A constant sigma causes large discretization error at the PML interface.
Polynomial grading ensures smooth transition.

### Key Parameters

| Parameter | Typical Value | Effect |
|-----------|---------------|--------|
| PML thickness L | 0.25*lambda | 10-20 elements in PML |
| sigma_max | Optimized (5-6 for lambda/20 mesh) | Absorption strength |
| Polynomial order m | 2 (parabolic) | Smoothness of grading |
| s' (real part) | 1.0 | Set > 1 to attenuate evanescent waves |

### Metal-backed PML Reflection

```
|R| = exp(-2 * k * s" * L * cos(theta))
```

Key implications:
- Reflection is angle-dependent (worse at grazing incidence)
- Increase s"*L product to reduce reflection
- But too large s" increases discretization error

## Advanced Techniques

### PML + ABC (Combined)

Replace metal backing with an ABC for better absorption at grazing angles:

```python
# 1st-order ABC on PML outer boundary
a += -1j*k0 * u.Trace() * v.Trace() * ds("pml_outer")
```

Combined PML+ABC is significantly better than either alone.

### Complementary PML (COM-PML)

Solve the problem TWICE:
1. PEC-backed PML (reflection = +R)
2. PMC-backed PML (reflection = -R)
3. Average the solutions -> error reduces from O(R) to O(R^2)

```python
# PEC-backed: standard Dirichlet on outer boundary
# PMC-backed: Neumann on outer boundary
# Average: u = (u_PEC + u_PMC) / 2
```

COM-PML achieves ~50 dB error reduction with only 2x computational cost.

## When to Use PML vs Other Open Boundary Methods

| Method | Best For | Limitation |
|--------|----------|------------|
| **PML** | Full-wave, broadband | Adds mesh volume, needs optimization |
| **ABC (1st/2nd order)** | Simple, narrowband | Low accuracy at wide angles |
| **FEM-BEM coupling** | Exact open boundary | Dense BEM matrix |
| **Radia (BEM)** | Magnetostatic/MQS | No wave propagation |

**For Radia users**: PML is NOT needed for magnetostatic or MQS problems.
Radia's integral method provides exact open boundary naturally.
PML is relevant for high-frequency (wave) problems in NGSolve.

## Practical Tips

1. **Start with radial PML** for 3D scattering/radiation problems
2. **Use polynomial grading** (m=2) - never constant sigma
3. **Optimize sigma_max** for your mesh density (start with 5-6 for lambda/20)
4. **Place PML at least 0.25*lambda from scatterer** (closer with COM-PML)
5. **Verify by varying PML thickness** - solution should be insensitive
6. **For waveguide problems**: use PML at waveguide ends, not sides

## Key References

- Berenger (1994): Original PML concept
- Chew & Weedon (1994): Coordinate stretching interpretation
- Sacks et al. (1995): Anisotropic absorber interpretation
- Jin & Chew (1996): PML optimization, combined PML+ABC
- Jin et al. (1997): Complementary PML
"""

NGSOLVE_DOMAIN_DECOMP = """
# Domain Decomposition for Large-Scale EM Problems

Reference: J.-M. Jin, "The Finite Element Method in Electromagnetics,"
3rd ed., Ch.14.3 (Wiley, 2014).

## FETI-DP Method (Finite Element Tearing and Interconnecting - Dual-Primal)

FETI-DP decomposes the computational domain into subdomains and solves
a global interface system iteratively. It is the most scalable domain
decomposition method for electromagnetics.

### Algorithm Overview

1. **Tearing**: Partition domain into Ns subdomains (METIS)
2. **Local factorization**: Factor each subdomain matrix [K_rr] independently
3. **Interface system**: Solve for Lagrange multipliers {lambda} on interfaces
4. **Recovery**: Compute subdomain solutions from factored matrices

```
Three types of unknowns:
- Interior: solved locally per subdomain
- Interface (dual): connected by Lagrange multipliers
- Corner (primal): assembled globally for coarse problem
```

### Interface System

```
[F_rr - F_rc * K_cc^{-1} * F_rc^T] {lambda} = {d_r - F_rc * K_cc^{-1} * f_c}
```

where F_rr and F_rc involve subdomain-local forward/backward substitutions.
Solved by CG (SPD) or GMRES (non-SPD).

### Low-Frequency Stability

For magnetostatic and eddy current problems, tree-cotree gauge must be
applied WITHIN each subdomain to avoid low-frequency breakdown:

```python
# Each subdomain independently applies tree-cotree splitting
# Corner edges are always cotree edges (globally connected)
```

### Numerical Scalability

| Test | Result |
|------|--------|
| Fixed size, varying #subdomains | Iterations constant (3-5) |
| Fixed #subdomains, varying size | Iterations constant (3-5) |
| Unit-load-per-processor | >80% parallel efficiency up to 128 cores |
| Nonlinear (Newton-FETI-DP) | Works with B-H curves |

Validated on:
- Eddy current: 1M DOFs, 128 subdomains, 4-8 iterations
- SRM motor (nonlinear): 3.3M DOFs, 128 cores, 400s total
- Phased arrays: up to 3.6 billion DOFs (300x300 array)

### Limitations

- Loses scalability when subdomains are electrically large (support resonant modes)
- Solution: FETI-DP2L (two Lagrange multipliers) for high-frequency
- For low-frequency EM (magnetostatic, eddy current): always scalable

## NGSolve Domain Decomposition

### BDDC (Balancing Domain Decomposition by Constraints)

BDDC is NGSolve's primary domain decomposition preconditioner, closely
related to FETI-DP (dual of each other):

```python
fes = HCurl(mesh, order=2)
a = BilinearForm(fes)
a += 1/mu * curl(u) * curl(v) * dx

# BDDC preconditioner
pre = Preconditioner(a, "bddc")
a.Assemble()

with TaskManager():
    solvers.CG(mat=a.mat, rhs=f.vec, sol=gfu.vec,
               pre=pre.mat, maxsteps=500, tol=1e-10)
```

### BDDC + AMG for Large Problems

```python
# AMG coarse solver for scalability beyond single-node
pre = Preconditioner(a, "bddc", coarsetype="h1amg")
```

### Parallel Execution with MPI

```python
# Run with: mpirun -np 16 python script.py
from ngsolve import *

# Mesh partitioned automatically
mesh = Mesh("large_problem.vol", comm=MPI_COMM_WORLD)

# Same code, parallelized automatically with BDDC
fes = HCurl(mesh, order=2)
pre = Preconditioner(a, "bddc")
```

## Dual-Field Domain Decomposition (DFDD)

For time-domain problems, DFDD computes both E and H fields:

1. E on integer time steps, H on half-integer steps (leapfrog)
2. Subdomains coupled via surface equivalence principle
3. No global system to solve -> naturally parallel
4. Each subdomain uses direct sparse solver (LU factored once)

Key advantage: memory and factorization time scale as O(N/M) where
M = number of subdomains (super-linear reduction for direct solvers).

## Solver Selection for Large Problems

| Problem Size | Recommendation |
|-------------|----------------|
| < 50K DOFs | Direct solver (PARDISO/UMFPACK) |
| 50K - 500K DOFs | BDDC + CG/GMRes |
| 500K - 5M DOFs | BDDC + AMG coarse + TaskManager |
| > 5M DOFs | MPI parallel BDDC + METIS partitioning |
| > 100M DOFs | FETI-DP (external, dedicated solver) |

## Fast Frequency Sweep (AWE/SSP)

For broadband analysis, avoid solving at each frequency independently.

### Asymptotic Waveform Evaluation (AWE)

```
1. Solve at expansion point k0:    m0 = A^{-1}(k0) * y(k0)
2. Compute moment vectors:         mn = A^{-1}(k0) * [y'(k0) - sum(A'*m_{n-i})]
3. Taylor series:                  x(k) = sum(mn * (k-k0)^n)
```

Matrix A factored ONCE at k0, reused for all moments.
Use Complex Frequency Hopping (CFH) for wide bands.

### Solution Space Projection (SSP)

More robust than AWE for wide frequency bands:
- Compute solutions at a few sampled frequencies
- Project onto reduced basis (Galerkin)
- Solve small reduced system at each frequency
- Adaptive binary search for sampling point selection

## Key References

- Farhat et al. (2001): FETI-DP for electromagnetics
- Li & Jin (2007): FETI-DP for low-frequency and high-frequency EM
- Yao et al. (2012): Nonlinear FETI-DP for SRM motors
- Lou & Jin (2006): Dual-field domain decomposition
- Slone et al. (2003): AWE for FEM frequency sweep
"""


NGSOLVE_MATERIAL_MODELING = """
# Material Modeling for Magnetic FEM

Based on: Takahashi Norio, "Magnetics FEM" (2013), "3D FEM" (2006).

## B-H Curve Handling in FEM

### Input Data Format
B-H curve data: list of [H(A/m), B(T)] pairs, sorted by H ascending.

```python
# NGSolve BH curve via CoefficientFunction
from ngsolve import *

# Piecewise linear BH data (H, B pairs)
BH_DATA = [(0, 0), (100, 0.5), (500, 1.0), (2000, 1.4),
           (10000, 1.7), (50000, 2.0)]

# Create nu(B) = H/B as function of |B|^2
def create_reluctivity(bh_data):
    H_vals, B_vals = zip(*bh_data)
    # Build piecewise linear interpolation of nu(B^2)
    # nu = H/B, with B^2 as independent variable for Newton
    b2_vals = [b**2 for b in B_vals]
    nu_vals = [h/b if b > 0 else H_vals[1]/B_vals[1]
               for h, b in zip(H_vals, B_vals)]
    # Use NGSolve's CF interpolation...
    return nu_vals, b2_vals
```

### Key Pitfall: Saturation Extrapolation
When B exceeds measured data range, the BH curve slope MUST approach mu_0:
```
For B > B_max_measured:
    H = H_last + (B - B_last) / mu_0
```
Saturation magnetization Ms values (Takahashi Table 6.2):
- Pure iron: Ms = 2.158 T
- 3% Si steel (35A250): Ms = 2.03 T
- Si formula: Ms = 2.158 - 0.048 * (Si content %)

## Anisotropy Modeling

### Method 1: Two B-H Curves (Simple but Inaccurate)
Uses rolling direction (RD) and transverse direction (TD) curves independently:
```
Hx = nu_RD(Bx) * Bx   (RD direction)
Hy = nu_TD(By) * By   (TD direction)
```

**PITFALL**: At high flux density, this gives physically wrong results -
permeability in hard direction exceeds easy direction. Use only for
rough estimates.

### Method 2: Multiple B-H Curves (Accurate)
Measure H(B, theta_B) and theta_H(B, theta_B) for multiple angles.
Requires 2D magnetic property measurement device.

Newton-Raphson needs 4 partial derivatives:
dH/dB, dH/d(theta_B), d(theta_H)/dB, d(theta_H)/d(theta_B)

**Key issue**: Coefficient matrix becomes NON-SYMMETRIC for anisotropy.
Must use BiCGSTAB (not CG/ICCG) as linear solver.

### Method 3: Fixed-Point Method (Recommended for Anisotropy)
Fixed-Point iteration avoids non-symmetric matrices:
```
H = nu_FP * B + H_FP(B)   (linear part + residual)
```

**Advantages over Newton-Raphson for anisotropic materials**:
1. Coefficient matrix is SYMMETRIC -> use ICCG (faster per iteration)
2. No Jacobian computation needed
3. No relaxation factor needed
4. Convergence: ~5x more iterations but ~2x faster total (ICCG vs BiCGSTAB)

Example: IPM motor with 35A300 steel, Fixed-Point 236 iterations (3915 ms)
vs Newton-Raphson 42 iterations (8900 ms) -> FP is 2.3x faster total.

**Implementation in NGSolve**:
```python
# Fixed-Point iteration for anisotropic materials
nu_FP = 1 / (mu_0 * mu_r_avg)  # constant linear part
H_FP = GridFunction(fes_H)      # residual field

for iteration in range(max_iter):
    # Solve linear system (symmetric!)
    a_FP = BilinearForm(nu_FP * curl(u) * curl(v) * dx)
    f_FP = LinearForm(J * v * dx + H_FP * curl(v) * dx)
    # ...solve...
    # Update H_FP from B-H curve residual
    B_new = curl(gf_A)
    H_exact = evaluate_anisotropic_BH(B_new)
    H_FP.Set(H_exact - nu_FP * B_new)
```

## E&S Model (Enokizono & Soda)

For rotating flux conditions (motors), the E&S model represents:
```
Hx = nu_xr * Bx + nu_xi * integral(Bx dt)
Hy = nu_yr * By + nu_yi * integral(By dt)
```
where nu_xr, nu_xi, nu_yr, nu_yi are measured coefficients that
vary with B_max, theta_B, and axis ratio alpha.

This captures both anisotropy AND rotational hysteresis from measured data.

## Key References

- Takahashi (2013): "Magnetics FEM", Ch.6 (Asakura Publishing)
- Enokizono & Soda (1990): E&S model for rotating flux
- Muramatsu et al. (2015): Fixed-Point method for anisotropic IPM motor
"""

NGSOLVE_IRON_LOSS = """
# Iron Loss Estimation in FEM

Based on: Takahashi Norio, "Magnetics FEM" (2013), Ch.6.4.

## Iron Loss Decomposition (Classical)

Total iron loss = Hysteresis loss + Eddy current loss + Anomalous loss

```
W_total = W_h + W_e + W_a

W_h = k_h * f * B_max^alpha   [W/kg]  (Steinmetz, alpha ~ 1.6-2.0)
W_e = k_e * f^2 * B_max^2     [W/kg]  (classical eddy current)
W_a = k_a * f^1.5 * B_max^1.5 [W/kg]  (anomalous/excess loss)
```

### Classical Eddy Current Loss (from Maxwell's equations)
For sinusoidal flux in a laminated steel sheet of thickness d:
```
W_e = (pi * d * B_max * f)^2 / (6 * rho * gamma)  [W/m^3]
    = sigma * d^2 * (dB/dt)^2 / 12                  [W/m^3, instantaneous]
```
where rho = resistivity, gamma = density, sigma = conductivity.

### Anomalous Eddy Current Loss
Domain wall motion creates localized eddy currents beyond classical prediction.
For non-oriented steel at 50 Hz, anomalous loss can be 40-60% of total eddy
current loss. Classical eddy current alone underestimates significantly.

## FEM-Based Iron Loss Computation

### Post-Processing Method (Standard)
```python
# After time-stepping FEM solve:
# 1. Extract B(t) waveform at each element
# 2. Compute FFT -> B_n for each harmonic n
# 3. Apply loss formula per harmonic

import numpy as np

def compute_iron_loss_element(B_waveform, dt, bh_coeffs):
    \"\"\"Compute iron loss for one element from B(t) waveform.\"\"\"
    N = len(B_waveform)
    f0 = 1.0 / (N * dt)

    # FFT of |B| waveform
    B_fft = np.fft.rfft(B_waveform)
    B_n = 2.0 * np.abs(B_fft) / N  # harmonic amplitudes

    k_h, alpha, k_e, k_a = bh_coeffs
    W = 0.0
    for n in range(1, len(B_n)):
        fn = n * f0
        Bn = B_n[n]
        W += k_h * fn * Bn**alpha      # hysteresis
        W += k_e * fn**2 * Bn**2        # eddy current
        W += k_a * fn**1.5 * Bn**1.5    # anomalous
    return W  # W/kg
```

### Locus-Based Method (for Rotating Flux)
For rotating machines, B traces an elliptical locus, not 1D alternation.
Decompose into orthogonal components:
```
W_total = W_alternating(B_major) + W_alternating(B_minor) * k_rot
```
where k_rot accounts for rotational hysteresis loss.

### Direct Integration Method (Most Accurate)
```
W_h = integral(H * dB)  per cycle   [J/m^3, from hysteresis loop area]
W_e = integral(sigma * d^2/12 * (dB/dt)^2 * dt)  per cycle
```

## Electrical Steel Classification (JIS)

### Non-Oriented (Motor Cores)
Format: [thickness*100]A[loss*100]
- 35A250: 0.35mm, W_15/50 <= 2.50 W/kg, B50 = 1.62 T
- 35A300: 0.35mm, W_15/50 <= 3.00 W/kg
- 50A350: 0.50mm, W_15/50 <= 3.50 W/kg

### Grain-Oriented (Transformer Cores)
Format: [thickness*100]P/G[loss*100]
- 35P135: 0.35mm high-oriented, W_17/50 <= 1.35 W/kg, B8 = 1.88 T
- 30G130: 0.30mm standard, W_17/50 <= 1.30 W/kg

### High-Frequency Steel
- 6.5% Si steel: 0.10mm, for motors > 400 Hz
- Ms = 1.8 T, conductivity = 0.122e7 S/m

## Stress Effects on Magnetic Properties

### Cutting Stress
Punching/laser cutting creates residual stress zone (~1-3 mm from edge).
Effect: increased coercivity, decreased permeability near cut edges.
For narrow teeth (< 5mm), degradation can be 20-40%.

**FEM modeling**: Apply degraded BH curve to elements near cut edges.

### Compressive Stress
External compressive stress significantly degrades BH curve.
At 10 MPa compression: permeability drops ~30% for non-oriented steel.
At 40 MPa: permeability drops ~60%.

Stress-dependent BH modeling:
```
B(H, sigma) = B_0(H) * k_sigma(sigma)
```
where k_sigma is stress degradation factor (measured).

## PWM Inverter Excitation

PWM carrier frequency harmonics cause additional iron loss:
- Carrier harmonics (e.g., 5 kHz): small B amplitude but high frequency
- Additional loss: 10-30% above sinusoidal excitation
- Must include harmonics up to 2x carrier frequency in loss calculation

## Key References

- Takahashi (2013): "Magnetics FEM", Ch.6.4 (comprehensive iron loss)
- Bertotti (1988): Three-component loss model
- Steinmetz (1892): Original hysteresis loss formula
"""

NGSOLVE_PRACTICAL_TECHNIQUES = """
# Practical FEM Techniques for Electrical Machines

Based on: Takahashi Norio, "Magnetics FEM" (2013) Ch.7,
          "3D FEM" (2006) Ch.6.

## Voltage Source Coupling (Field-Circuit)

When coil current is unknown (voltage-driven), couple FEM with circuit:

### Basic Equation
```
d(phi)/dt + R*I + L_end*dI/dt = V_source
```
where phi = flux linkage from FEM, R = coil resistance,
L_end = end-winding inductance, V_source = applied voltage.

### NGSolve Implementation
```python
# Coupled field-circuit in NGSolve
fes_A = HCurl(mesh, order=2, complex=True)
# Circuit variable: current as additional unknown
fes_I = NumberSpace(mesh)
fes = fes_A * fes_I

(u, dI), (v, dv) = fes.TnT()

# FEM part: curl-curl
a = BilinearForm(fes)
a += nu * curl(u) * curl(v) * dx
a += sigma * u * v * dx          # eddy current term

# Coupling: current source term
n_turns = 100
S_coil = coil_area
J0 = n_turns / S_coil
a += -J0 * dI * v * dx_coil      # current -> field
a += J0 * u * dv * dx_coil       # field -> circuit (flux linkage)

# Circuit equation
a += (R + 1j*omega*L_end) * dI * dv
f = LinearForm(fes)
f += V_source * dv
```

### Current Direction Vector T
For 3D coils, define current direction via current vector potential T:
```
J0 = curl(T),  integral(T . dl) = I0  (along coil path)
```
This ensures div(J0) = 0 automatically (required for edge elements).

## Electromagnetic Force and Torque

### Method 1: Maxwell Stress Tensor (MST)
```
F_i = oint_S T_ij * n_j dS

T_ij = (1/mu_0) * (B_i*B_j - 0.5*delta_ij*|B|^2)
```
Integration surface S must be entirely in air.

```python
# Maxwell stress tensor in NGSolve
B = curl(gf_A)
Bx, By, Bz = B[0], B[1], B[2]
B2 = Bx**2 + By**2 + Bz**2

# Force in x-direction
Txx = (1/mu_0) * (Bx*Bx - 0.5*B2)
Txy = (1/mu_0) * (Bx*By)

Fx = Integrate(Txx*n[0] + Txy*n[1], mesh, BND,
               definedon=mesh.Boundaries("air_gap"))
```

**PITFALL**: Integration surface must NOT cross material boundaries.
Place surface in air gap for best accuracy.

### Method 2: Virtual Work Principle
```
F = -dW/dx  (W = magnetic energy or coenergy)

Torque = -dW/d(theta)
```
Compute by two solutions with small displacement:
```
F_x = -(W(x+dx) - W(x-dx)) / (2*dx)
```

For coenergy (more accurate with nonlinear materials):
```
W_co = integral(B * dH)  (coenergy, integrate B with respect to H)
F = +dW_co/dx            (note: positive sign for coenergy)
```

### Method 3: Nodal Force (Eggshell + SetDeformation)
Uses Coulomb virtual work: F_i = -dW/dx_i with frozen A (no re-solve).

**MST vs Nodal Force -- accuracy is equivalent for linear problems**
(Henrotte & Hameyer 2004). Both extract force from the same FEM solution,
so the accuracy bottleneck is the field quality, not the force method.
Nodal force advantages are practical, not accuracy:
- 3D complex geometry: no need to define a clean integration surface in air
- Structural coupling: per-node forces -> direct FEM structural load
- Torque: virtual rotation, no lever-arm calculation
- Nonlinear: coenergy formulation is rigorous for B-H curves

**Key idea**: Surround the body with a thin "eggshell" air layer.
Create a blending function (1 inside body, 0 outside eggshell).
Use mesh.SetDeformation() + central difference for the energy derivative.

**CRITICAL**: Use `mesh.Curve(3)` for curved geometry. Without curving,
polygon approximation of circles loses ~2% area -> ~9% force error.

```python
# --- Eggshell nodal force method (verified implementation) ---
# Geometry: create eggshell ring around iron body via OCC boolean ops
# iron_shape = wp.MoveTo(d, 0).Circle(R_cyl).Face()
# iron_shape.name = "iron"
# egg_outer = wp.MoveTo(d, 0).Circle(R_cyl + t_egg).Face()
# eggshell = egg_outer - iron_shape
# eggshell.name = "eggshell"
# ... then Glue([wire, iron, eggshell, air])
# mesh.Curve(3)  # ESSENTIAL for curved boundaries

# After solving A-formulation (gf_A = solution):
# 1. Blending function: 1 on iron, 0 outside eggshell
fes_blend = H1(mesh, order=1)
gf_blend = GridFunction(fes_blend)
gf_blend.Set(1.0, definedon=mesh.Materials("iron"))

# 2. Virtual displacement V = (blend * delta, 0)
fes_def = VectorH1(mesh, order=1)
gf_def = GridFunction(fes_def)
delta = 1e-7  # perturbation [m]

nu_cf = 1.0 / mu_cf  # reluctivity
energy_density = 0.5 * nu_cf * InnerProduct(grad(gf_A), grad(gf_A))

# 3. Central difference
gf_def.components[0].vec.data = delta * gf_blend.vec
gf_def.components[1].vec[:] = 0.0
mesh.SetDeformation(gf_def)
W_plus = Integrate(energy_density, mesh, order=2*order+4)
mesh.UnsetDeformation()

gf_def.components[0].vec.data = -delta * gf_blend.vec
mesh.SetDeformation(gf_def)
W_minus = Integrate(energy_density, mesh, order=2*order+4)
mesh.UnsetDeformation()

Fx = -(W_plus - W_minus) / (2 * delta)
# For torque: use tangential virtual displacement instead
```

**Verified**: Wire + iron cylinder benchmark (2D, linear, mu_r=1000;
I=1000 A, R_cyl=0.05 m, d=0.3 m; analytic F_x = -1.90096e-02 N/m):
- nodal force matches analytical to ~0.01% even on the COARSEST mesh
  (685 elems, +0.01%), where MST still has +0.28% error
- at the finest mesh MST and nodal force agree within ~6e-4%
- both match analytical to ~0.01% (with mesh.Curve(3))
- without mesh.Curve(): ~9% systematic error for BOTH methods equally
See the self-contained notebook
`docs/nodal_force/nodal_force.ipynb` (executed code and embedded plots) and
`validation_test/nodal_force/nodal_force_results.json` (checked numerical evidence).

**When to choose which method:**
| Situation | Recommended | Reason |
|-----------|-------------|--------|
| Quick global force check | MST | Simple, one-line Integrate |
| 3D iron body (complex shape) | Nodal force | No integration surface needed |
| Structural FEM coupling | Nodal force | Per-node force distribution |
| Torque in rotating machines | Nodal force | Virtual rotation is natural |
| Nonlinear B-H material | Nodal force | Coenergy formulation rigorous |
| Linear, simple geometry | Either | Equivalent accuracy |

## Rotation Handling for Motors

### Method 1: Remeshing (simple but slow)
Regenerate air gap mesh at each angular position.

### Method 2: Sliding Surface with Unknown Equipotential
```python
# Floating-point nodes on stator-rotor interface
# Connect via constraint equations
# A_stator(theta) = A_rotor(theta + theta_mech)
```

### Method 3: Locked-Step (Moving Band)
Rotate mesh by one element width per time step.
Requires air gap mesh with uniform angular divisions.

## Gap Elements (Thin Air Gap Modeling)

For thin gaps (gap << surrounding dimensions), use line elements:
```
W_gap = (1/2) * nu_0 * (A1 - A2)^2 / L_gap * D_gap * L_element
```
No area needed - attach to existing mesh without remeshing.

**Advantages**:
1. Add/remove gaps without mesh modification
2. Easy gap length parametric studies
3. More accurate than flat triangles for very thin gaps

## Magneto-Thermal Coupled Analysis

### Coupling Strategy
```
1. Solve magnetic field -> compute losses (W_e + W_h)
2. Use losses as heat source -> solve thermal
3. Update material properties at new temperature
4. Repeat until temperature converges
```

### Temperature-Dependent Properties
- Conductivity: sigma(T) = sigma_20 / (1 + alpha*(T - 20))
  (alpha = 0.004 for copper, 0.0039 for aluminum)
- Permanent magnet: Br(T) = Br_20 * (1 + alpha_Br*(T - 20))
  (alpha_Br ~ -0.001/K for NdFeB, -0.0004/K for ferrite)
- BH curve: shifts toward lower B at higher temperature

### Thermal FEM in NGSolve
```python
# Thermal analysis
fes_T = H1(mesh, order=2, dirichlet="outer")
T, s = fes_T.TnT()

# Heat conduction
k_thermal = mesh.MaterialCF({"iron": 50, "copper": 400,
                              "air": 0.026})
a_th = BilinearForm(k_thermal * grad(T) * grad(s) * dx)

# Convection BC: q = h*(T - T_ambient)
h_conv = 10  # W/(m^2*K)
a_th += h_conv * T * s * ds("surface")

# Heat source from EM solve
f_th = LinearForm(Q_density * s * dx)
f_th += h_conv * T_ambient * s * ds("surface")
```

## Inductance and EMF Computation

### Flux Linkage from Vector Potential
```
phi = (n_turns / S_coil) * integral(A . n_coil dV)
     = (n_turns / S_coil) * integral(A_z dS)   (2D)
```

### Inductance
```
L = phi / I  = (n^2 / S^2) * integral(A . J0 dV) / I^2
  = 2*W_magnetic / I^2  (from energy)
```

### Back-EMF
```
EMF = -d(phi)/dt
    = -(n/S) * integral(dA/dt . n_coil dV)
```
In time-stepping: use backward difference (A^n - A^{n-1})/dt.

## Key References

- Takahashi (2013): "Magnetics FEM", Ch.7 (Asakura Publishing)
- Takahashi (2006): "3D FEM", Ch.6 (IEEJ)
- Nakata & Takahashi (1995): "FEM for EE", Ch.6-7 (Morikita Publishing)
"""


NGSOLVE_ELECTROMAGNETIC_FORCE = """
# Electromagnetic force/torque extraction map in radia-ngsolve

Use this as the front door for force work.  The detailed validation tables live
in `force_validation("method_map")` and `force_validation("cross_validation")`;
this topic connects the practical NGSolve helper functions to the physics.

## Fast method selection

| Problem | Use | Function/topic |
|---|---|---|
| Global force on a body surrounded by air | weighted Maxwell stress | `force.eggshell_force`, `force.eggshell_force_2d` |
| Global torque from the same field | weighted Maxwell stress torque | `force.eggshell_torque`, `force.eggshell_torque_2d` |
| Clean closed surface in air | Maxwell surface stress | `force.maxwell_surface_force` |
| Complex harmonic field | time-averaged Maxwell stress | `force.maxwell_surface_force_harmonic` |
| Current conductor / busbar | Lorentz body force | `force.lorentz_force_2d`, `ngsolve_usage("lorentz_force")` |
| Axisymmetric actuator or coil pull | `2*pi*r` weighted stress | `force.eggshell_force_axi` |
| Air-gap holding force | magnetic pressure | `force.air_gap_force_summary`, `ngsolve_usage("air_gap_force")` |
| Motor dq operating torque | lumped dq torque | `solve.dq_torque`, `ngsolve_usage("dq_torque")` |
| Electrostatic/MEMS force | electric weighted stress | `force.electrostatic_eggshell_force_2d` |

## The key convention

Maxwell stress is a stress in the surrounding air.  Put the integration surface
or eggshell band fully in air:

```text
T = (1/mu0) (B tensor B - 0.5 |B|^2 I)
F = integral_S T n dS
```

The eggshell form is the same physics written as a volume integral with a
smooth weight `g` that goes from 1 near the body to 0 outside the band:

```text
F_k = - integral_band [(1/mu0) B_k (B.grad g)
                       - (1/(2 mu0)) |B|^2 d_k g] dV
```

This is usually the robust FEM answer: it avoids noisy field traces on material
interfaces and averages over a band of elements.

## Checks before trusting a force number

- Curve circular/curved geometry; force is quadratic in `B`, so geometry error is
  amplified.
- Compare against at least one analytic anchor when possible: two-wire force,
  air-gap pressure, magnetized sphere, dq torque identity, or energy/inductance.
- In 2D planar solves, report force per unit length.  Axisymmetric helpers return
  full 3D quantities because they include the `2*pi*r` weight.
- For nonlinear B-H, keep the body directly in air and use an analytic eggshell
  band; do not introduce a nested artificial shell that isolates the HCurl field.
- For harmonic fields, use time-average factors (`1/(2 mu0)` and `1/(4 mu0)`),
  not the static stress tensor.

## Good companion topics

- `force_validation("method_map")`: fuller method map and validation anchors.
- `force_validation("eggshell")`: derivation and implementation details.
- `ngsolve_usage("lorentz_force")`: busbar and image-force recipes.
- `ngsolve_usage("air_gap_force")`: pole-face pressure and magnetic circuit force.
- `ngsolve_usage("electrostatic_force")`: MEMS electric-force chain.
"""


NGSOLVE_AXISYMMETRIC = """
# Axisymmetric Magnetostatics in NGSolve

NGSolve does NOT have a built-in axisymmetric mode. Axisymmetric problems must be
formulated manually in the 2D (r,z) half-plane with appropriate weighted weak forms.

## The Problem: 1/r Singularity

In axisymmetric coordinates with A = A_theta(r,z) e_theta:
- B_r = -dA/dz
- B_z = (1/r) d(rA)/dr = A/r + dA/dr

The energy functional contains an A^2/r term (singular at r=0).
FEMM handles this with **analytical integration** (closed-form log-mean radius).
NGSolve uses **numerical quadrature**, so a different approach is needed.

## Recommended: Mixed Formulation (Joachim Schoeberl)

Uses phi = rA (flux function) as H1 unknown and B as VectorL2 unknown.
The 1/r singularity is eliminated by the formulation structure.

Reference: TEAM Problem 28 benchmark (2024, K. Sugahara / J. Schoeberl).

```python
from ngsolve import *
from netgen.occ import *

# Geometry: 2D half-plane (r >= 0), x = r, y = z
# Build geometry with OCC (dim=2)
mesh = Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=0.005))

# Mixed FE spaces
p = 3  # polynomial order (p >= 2 recommended)
fesPhi = H1(mesh, order=p, dirichlet="outer")    # phi = rA (flux function)
fesB = VectorL2(mesh, order=p-1)                  # B = (Br, Bz) directly
fes = fesPhi * fesB
(u, B), (v, dB) = fes.TnT()

# B operator: maps phi -> B via curl of (phi/r * e_theta)
# Br = -(1/r) * dphi/dz,  Bz = (1/r) * dphi/dr
def Bop(phi):
    return CF((-1/x * grad(phi)[1], 1/x * grad(phi)[0]))

# Material: reluctivity nu = 1/mu
mu0 = 4 * pi * 1e-7
nu = mesh.MaterialCF({"iron": 1/(1000*mu0)}, default=1/mu0)

# Bilinear form: saddle-point system
# nu * r * (B . dB - Bop(u) . dB + Bop(v) . B)
bfa = BilinearForm(nu * x * (B*dB - Bop(u)*dB + Bop(v)*B) * dx).Assemble()

# Right-hand side: current density J_theta (A/m^2)
Jz = mesh.MaterialCF({"coil": 1e6}, default=0)  # 1 MA/m^2
lff = LinearForm(v * Jz * dx).Assemble()

# Solve
gfu = GridFunction(fes)
gfu.vec.data = bfa.mat.Inverse(fes.FreeDofs()) * lff.vec

phi_sol = gfu.components[0]  # phi = rA (flux function)
B_sol = gfu.components[1]    # B = (Br, Bz) directly available
```

## Why Mixed Formulation Works

1. **No 1/r singularity**: Bop has 1/x but multiplied by grad(phi);
   the bilinear form weight `x` (= r) cancels the singularity.
2. **B is computed directly**: No need to differentiate A to get B.
   VectorL2 gives B per element without post-processing accuracy loss.
3. **VectorL2 is discontinuous**: B can be eliminated element-wise
   via static condensation (Schur complement), reducing to a phi-only system.
4. **Arbitrary polynomial order**: p=2,3,4,... unlike FEMM (fixed p=2).

## Saddle-Point Structure

The system has the structure:
```
[ M    -G^T ] [B  ]   [0]
[ G     0   ] [phi] = [f]
```
where M is the VectorL2 mass matrix (block-diagonal, element-local)
and G is the gradient operator from Bop. Since M is block-diagonal,
B can be eliminated cheaply: G M^{-1} G^T phi = f.

## Coordinate Convention

- `x` = r (radial direction, r >= 0)
- `y` = z (axial direction)
- Geometry is the right half-plane (r >= 0)
- Use OCC 2D geometry: `OCCGeometry(shape, dim=2)`
- Right half only: subtract left half from full shapes
  ```python
  circle = MoveTo(0, 0).Circle(R).Face()
  left_half = MoveTo(-R, -R).Rectangle(R, 2*R).Face()
  half_circle = circle - left_half
  ```

## Boundary Conditions

- **Axis (r=0)**: phi = rA = 0 is NATURAL (no Dirichlet needed on axis).
  The flux function phi = rA vanishes at r=0 by definition.
- **Outer boundary**: Dirichlet phi = 0 (far-field decay)
- **Symmetry planes**: Neumann (natural BC, nothing to specify)

## Recovering A and Flux

```python
# Vector potential: A_theta = phi / r
r_coord = IfPos(x - 1e-10, x, 1e-10)
A_theta = phi_sol / r_coord

# Magnetic flux through circle at (r, z):
# Phi_flux = 2*pi*r*A = 2*pi*phi
Phi_flux = 2 * pi * phi_sol

# B-field is directly available from B_sol = (Br, Bz)
# For plotting: Draw(B_sol, mesh)
```

## Comparison with FEMM

| Feature | FEMM | NGSolve Mixed |
|---------|------|---------------|
| Dimension | 2D only | 2D + 3D |
| Element order | p=2 fixed | p=2,3,4,... arbitrary |
| Integration | Analytical (exact for p=2) | Numerical quadrature |
| 1/r handling | Log-mean radius (R_hat) | Mixed formulation eliminates |
| B computation | Post-processed from A | Direct (VectorL2 variable) |
| Nonlinear | Newton-Raphson | Newton (via NonlinearSolve) |
| Open boundary | Asymptotic BC (Kelvin-like) | Kelvin transformation |

**Performance note**: At p=2 with the same mesh, FEMM produces smoother results
due to analytical integration. NGSolve catches up at p=3 or with finer mesh.
For p >= 3, NGSolve has no FEMM equivalent (FEMM is fixed at p=2).

## Nonlinear Axisymmetric (BH Curve)

```python
# Newton's method for nonlinear mu(B)
def solve_nonlinear_axi(mesh, bh_curve, Jz, p=3, tol=1e-6, maxiter=50):
    fesPhi = H1(mesh, order=p, dirichlet="outer")
    fesB = VectorL2(mesh, order=p-1)
    fes = fesPhi * fesB
    (u, B), (v, dB) = fes.TnT()

    gfu = GridFunction(fes)

    for it in range(maxiter):
        B_val = gfu.components[1]
        Bmag = sqrt(B_val[0]**2 + B_val[1]**2 + 1e-20)
        nu_nl = bh_curve.get_nu(Bmag)  # user-defined nu(|B|)

        bfa = BilinearForm(nu_nl * x * (B*dB - Bop(u)*dB + Bop(v)*B) * dx)
        bfa.Assemble()
        lff = LinearForm(v * Jz * dx).Assemble()

        res = lff.vec - bfa.mat * gfu.vec
        if Norm(res) < tol:
            break
        gfu.vec.data += bfa.mat.Inverse(fes.FreeDofs()) * res

    return gfu
```

## Common Mistakes

1. **Using H1 scalar with weighted form instead of mixed**:
   `a += nu/r * grad(u)*grad(v) * dx` works but is less accurate than
   mixed formulation at the same polynomial order, and has 1/r issues near axis.

2. **Forgetting x = r convention**: In NGSolve 2D, `x` is the first coordinate.
   For axisymmetric, `x` must be the radial direction (r >= 0).

3. **Using VectorH1 instead of VectorL2 for B**: VectorL2 allows element-local
   elimination. VectorH1 would over-constrain B continuity.

4. **Applying Dirichlet on axis for phi**: phi = rA = 0 at r=0 is natural;
   DO NOT set Dirichlet on the axis boundary (it's redundant and can cause issues
   with the mixed system).
"""


NGSOLVE_ACCELERATOR_MAGNET = """
# Accelerator Magnet Solver (Radia + NGSolve + Cubit)

Complete pipeline for accelerator electromagnet analysis.

## Architecture

```
CoilBuilder → Biot-Savart Hs (no coil mesh)
Cubit → hex mesh + export netgen "mesh.vol" order N
NGSolve → Omega-reduced Omega + Kelvin + Energy-based B-input Play
```

## Key Components

### 1. CoilBuilder (radia.coil_builder)
```python
from radia.coil_builder import CoilBuilder

coil = (CoilBuilder(current=10000)
    .set_cross_section(width=15e-3, height=25e-3)
    .add_straight(80e-3)
    .add_arc(radius=20e-3, arc_angle=180)
    .add_straight(80e-3)
    .add_arc(radius=20e-3, arc_angle=180)
    .close())                    # verify/optimize loop closure

lower = coil.mirror('xy')       # dipole symmetry
coil.write_step("coil.step")    # GMSH visualization
radia_objects = coil.to_radia() # Biot-Savart source (NO mesh)
```

### 2. Cubit Hex Mesh + Curving
```python
import tempfile
from ngsolve import Mesh
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)
# ACIS CallbackGeometry: no STEP, no OCC, no seam problems
```

### 3. Omega-reduced Omega + Kelvin
```python
# 2-scalar method (fastest formulation for magnetostatics)
# Source: Hs from Radia Biot-Savart (CoilBuilder.to_radia())
# Iron: H = nu_rev * B + H_irr(B)  (Energy-based Play model)
# Kelvin: open boundary transformation (no PML)
# LU factored ONCE, iteration = back-substitution only
```

### 4. Energy-based B-input Play Hysteresis
```
H(B) = H_rev(B) + H_irr(B, history)
H_rev = nu_rev * B              (reversible, constant)
H_irr = sum_k g_k(p_k(B))      (irreversible, Play operators)

Separation ensures:
- g_k >= 0 for all k (convex energy)
- Fast inverse: Picard 2-3 iter (vs Newton ~100)
- Hantila polarization: LU once, back-sub iteration
```

## Unique Capabilities

| Feature | OPERA | ANSYS | COMSOL | Radia+NGSolve |
|---------|-------|-------|--------|---------------|
| Hex + arbitrary order | - | - | - | Yes |
| No coil mesh | - | - | - | Yes (Biot-Savart) |
| No STEP/OCC needed | - | - | - | Yes (ACIS direct) |
| Kelvin transform | - | - | - | Yes |
| 2-scalar method | - | - | - | Yes |
| B-input hysteresis | - | - | - | Yes (Play model) |
| Energy-based | - | - | - | Yes (convex) |
| Open source | - | - | - | Yes (pip install) |

## Reference
"""


NGSOLVE_TEAM7 = """
# TEAM Problem 7 - Asymmetrical Conductor with a Hole (A-Formulation)

Reference: https://www.compumag.org/wp/wp-content/uploads/2018/06/problem7.pdf
NGSolve implementation: https://ngsolve.github.io/TEAM-problems/TEAM-7/team7.html

## Problem Description

3D eddy current benchmark: a rectangular aluminium plate with a rectangular hole,
excited by a racetrack coil above it. Solve the frequency-domain A-formulation:

    curl(nu * curl(A)) + j*omega*sigma*A = J_src

Boundary condition: A x n = 0 on the far boundary (airbox).

### Physical Parameters

| Parameter | Value |
|-----------|-------|
| mu_0 | 4*pi*1e-7 H/m |
| sigma_Al | 35.26e6 S/m |
| Frequency | 50 Hz and 200 Hz |
| Coil turns | 2742 |

## Geometry (NetGen OCC)

Three regions: aluminium plate, racetrack coil, air box.

```python
from netgen.occ import *
from ngsolve import *
import numpy as np

# ---- Aluminium plate (294 x 294 x 19 mm) with rectangular hole ----
wpplate = WorkPlane(Axes((0,0,0), Z, X))
f = wpplate.Rectangle(0.294, 0.294).MoveTo(0.018, 0.018).Rectangle(0.108, 0.108).Reverse().Face()
plate = f.Extrude(0.019)
plate.faces.maxh = 0.04
plate.solids.name = "plate"

# ---- Racetrack coil (25 mm cross-section, 50 mm offset, 100 mm extent) ----
wp = WorkPlane(Axes(p=(0.294, 0.000, 0.049), n=Z, h=Y))
coil2di = wp.MoveTo(0.100, 0.100).RectangleC(0.100, 0.100).Offset(0.025).Face()
coil2di.edges.name = "coili"
coil2do = wp.MoveTo(0.100, 0.100).RectangleC(0.100, 0.100).Offset(0.050).Face()
coil2d = coil2do - coil2di

# Split coil at corners for current direction assignment
cutout = wp.MoveTo(0.100, 0.100).RectangleC(0.300, 0.100).Face() + \\
         wp.MoveTo(0.100, 0.100).RectangleC(0.100, 0.300).Face()
coil2d = Glue([coil2d * cutout, coil2d - cutout])
coil = coil2d.Extrude(0.100)
coil.solids.name = "coil"

# ---- Air box ----
airbox = Box(Pnt(-0.2, -0.2, -0.2), Pnt(0.5, 0.5, 0.5))
airbox.faces.name = "outer"
airbox.solids.name = "air"
airbox = airbox - coil - plate

shape = Glue([coil, plate, airbox])
geo = OCCGeometry(shape)
```

## Mesh Generation

All-tet mesh with high-order curving:

```python
mesh = Mesh(geo.GenerateMesh(maxh=0.2))
mesh.Curve(5)
print("materials:", set(mesh.GetMaterials()))
print("number of elements:", mesh.ne)
```

## Coil Current Model

The coil current density uses a projection-based tangential vector field and
a divergence-free correction via a scalar potential:

```python
coilrect_xmin, coilrect_xmax = 0.294 - 0.150, 0.294 - 0.050
coilrect_ymin, coilrect_ymax = 0.050, 0.150

def Project(val, minval, maxval):
    return IfPos(val - minval, IfPos(val - maxval, maxval, val), minval)

projx = Project(x, coilrect_xmin, coilrect_xmax)
projy = Project(y, coilrect_ymin, coilrect_ymax)

# Tangential direction (follows racetrack path)
tau_coil = CF((projy - y, x - projx, 0))
tau_coil /= Norm(tau_coil)

# Scalar potential for div-free correction
pot_coil = CF((0, 0, sqrt((projy - y)**2 + (projx - x)**2) - 0.05))
```

## FEM Solve (A-Formulation)

Key features of this implementation:
- HCurl order=5 for high accuracy
- `gradientdomains="plate"` for tree-cotree gauge in the conducting region
- Boundary current + div-free correction (not volume current density)
- BDDC preconditioner with sparse Cholesky

```python
mu = 4 * pi * 1e-7
sigma = mesh.MaterialCF({"plate": 35.26e6}, default=1)
turns = 2742

def Solve(freq):
    omega = 2 * pi * freq
    fes = HCurl(mesh, order=5, complex=True, dirichlet="outer",
                gradientdomains="plate")
    print("free dofs:", sum(fes.FreeDofs()))
    u, v = fes.TnT()

    a = BilinearForm(fes, symmetric=True, condense=True)
    a += 1/mu * curl(u) * curl(v) * dx
    a += 1j * omega * sigma * u * v * dx

    # Coil excitation: boundary current + div-free correction
    f = LinearForm(
        -turns / 0.100 * tau_coil * v.Trace() * ds("coili.*", bonus_intorder=4)
        + turns / (0.025 * 0.100) * pot_coil * curl(v) * dx("coil.*", bonus_intorder=4)
    )

    A = GridFunction(fes)
    pre = Preconditioner(a, type="bddc", inverse="sparsecholesky")
    solvers.BVP(bf=a, lf=f, gf=A, pre=pre,
                solver=solvers.CGSolver,
                solver_flags={"plotrates": True, "tol": 1e-12})

    B = curl(A)
    J = -1j * omega * sigma * A
    return {"A": A, "B": B, "J": J}

with TaskManager():
    ret = Solve(freq=50)
B, J = ret["B"], ret["J"]
```

### Notes on the LinearForm (Coil Excitation)

Two equivalent approaches are available:

1. **Volume current density** (simpler but requires coil mesh refinement):
   ```python
   f = LinearForm(-turns / (0.025 * 0.100) * tau_coil * v * dx("coil.*", bonus_intorder=4))
   ```

2. **Boundary current + div-free correction** (used in TEAM-7, more accurate):
   ```python
   f = LinearForm(
       -turns / 0.100 * tau_coil * v.Trace() * ds("coili.*", bonus_intorder=4)
       + turns / (0.025 * 0.100) * pot_coil * curl(v) * dx("coil.*", bonus_intorder=4)
   )
   ```
   The second form uses a surface integral on the inner coil boundary ("coili")
   plus a curl(v) volume integral with a scalar potential for divergence-free correction.
   This is more accurate because it avoids discretization error from approximating
   a uniform current density in the coil cross-section.

### Key Implementation Details

| Feature | Choice | Reason |
|---------|--------|--------|
| `symmetric=True` | Complex symmetric (not Hermitian) | A-formulation with jωσ gives complex symmetric matrix |
| `condense=True` | Static condensation | Eliminates internal DOFs, reduces system size |
| `gradientdomains="plate"` | Tree-cotree gauge | Removes gradient null space in conducting region |
| `dirichlet="outer"` | Essential BC on airbox | A x n = 0 (zero tangential A) |
| `bonus_intorder=4` | Extra quadrature | Accurate integration of coil source term |
| `inverse="sparsecholesky"` | BDDC subdomain solver | Exact local solves for robustness |
| `order=5` | High polynomial order | Benchmark accuracy with coarse mesh |

## Post-Processing: Comparison with Benchmark Data

Measurement lines for validation:
- A1-B1: Bz at y=72mm, z=34mm (above plate)
- A2-B2: Bz at y=144mm, z=34mm (above plate, offset)
- A3-B3: Jy at y=72mm, z=19mm (plate top surface)
- A4-B4: Jy at y=72mm, z=0mm (plate bottom surface)

```python
xi_msm = np.array([0, 18, 36, 54, 72, 90, 108, 126, 144, 162,
                    180, 198, 216, 234, 252, 270, 288]) * 1e-3
xi_sim = np.linspace(0, 0.288, 500)

# Evaluate Bz along line A1-B1
Bz_A1B1_sim = np.array([B[2](mesh(xi, 0.072, 0.034)) for xi in xi_sim])

# Evaluate Jy along line A3-B3 (plate top)
Jy_A3B3_sim = np.array([1e-6 * 1j * J[1](mesh(xi, 0.072, 0.019 - 1e-5))
                         for xi in xi_sim])
```

### Reference Values (50 Hz)

```python
# Bz on A1-B1 (Gauss = 1e-4 T)
Bz_A1B1_50Hz_ref = np.array([
    -4.9-1.16j, -17.88+2.48j, -22.13+4.15j, -20.19+4.j, -15.67+3.07j,
    0.36+2.31j, 43.64+1.89j, 78.11+4.97j, 71.55+12.61j, 60.44+14.15j,
    53.91+13.04j, 52.62+12.4j, 53.81+12.05j, 56.91+12.27j, 59.24+12.66j,
    52.78+9.96j, 27.61+2.26j])

# Bz on A2-B2 (Gauss)
Bz_A2B2_50Hz_ref = np.array([
    -1.83-1.63j, -8.5-0.6j, -13.6-0.43j, -15.21+0.11j, -14.48+1.26j,
    -5.62+3.4j, 28.77+6.53j, 60.34+10.25j, 61.84+11.83j, 56.64+11.83j,
    53.4+11.01j, 52.36+10.58j, 53.93+10.8j, 56.82+10.54j, 59.48+10.62j,
    52.08+9.03j, 26.56+1.79j])

# Jy on A3-B3 (1e6 A/m^2)
Jy_A3B3_50Hz_ref = np.array([
    0.249-0.629j, 0.685-0.873j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j,
    -0.015-0.593j, -0.103-0.249j, -0.061-0.101j, -0.004-0.001j,
    0.051+0.087j, 0.095+0.182j, 0.135+0.322j, 0.104+0.555j,
    -0.321+0.822j, -0.687+0.855j])

# Jy on A4-B4 (1e6 A/m^2)
Jy_A4B4_50Hz_ref = np.array([
    0.461-0.662j, 0.621-0.664j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j,
    1.573-1.027j, 0.556-0.757j, 0.237-0.364j, 0.097-0.149j,
    -0.034+0.015j, -0.157+0.154j, -0.305+0.311j, -0.478+0.508j,
    -0.66+0.747j, -1.217+1.034j])
```

## Applicability to Shifted AMS + COCR

This TEAM-7 problem is an ideal test case for the shifted AMS preconditioner:
- Complex symmetric system (A-formulation with jωσ)
- Air + conductor domains (σ=0 in air makes system singular)
- Standard approach uses `gradientdomains` for gauge fixing; shifted AMS
  provides an alternative that regularizes only the preconditioner

To use shifted AMS + COCR instead of BDDC:

```python
from radia.sparsesolv_ngsolve import CompactAMSPreconditioner, COCRSolver

# Original system (no regularization)
a = BilinearForm(fes, symmetric=True)
a += 1/mu * curl(u) * curl(v) * dx
a += 1j * omega * sigma * u * v * dx

# Shifted system (preconditioner only)
a_shift = BilinearForm(fes, symmetric=True)
a_shift += 1/mu * curl(u) * curl(v) * dx
a_shift += 1j * omega * sigma * u * v * dx
a_shift += eps_shift * u * v * dx  # shift in all domains

a.Assemble()
a_shift.Assemble()

ams = CompactAMSPreconditioner(a_shift.mat, fes, coarsetype="sparsecholesky")
cocr = COCRSolver(a.mat, ams, tol=1e-12, maxiter=200)
A.vec.data = cocr * f.vec
```

## Reference

- COMPUMAG TEAM Problem 7: https://www.compumag.org/wp/wp-content/uploads/2018/06/problem7.pdf
- NGSolve TEAM Problems: https://ngsolve.github.io/TEAM-problems/TEAM-7/team7.html
"""


NGSOLVE_BOOLEAN_POLICY = """
# Boolean ops in netgen.occ: Glue vs Fuse vs Compound vs ACIS-unite

Three OCC / netgen.occ operations combine multiple shapes; each behaves
differently under "trampolining" (passing geometry through pipelines).
Choosing wrong silently breaks FEM meshing: interface nodes duplicate,
UV charts degenerate, or material labels collapse.

## The three OCC operations

| Op | Boolean runs? | Solids in output | UV preserved? | Interface shared? |
|----|---------------|------------------|----------------|-------------------|
| `Glue([s1, s2, ...])` | NO (annotation only) | N (all inputs kept) | YES | YES, if inputs already share faces |
| `Fuse(s1, s2)` / `s1 + s2` | YES | 1 | Re-edited | One solid, no interface |
| `Fuse(..., SetGlue(BOPAlgo_GlueFull))` (OCCT C++) | YES, with hint | 1 | Preserved | One solid, no interface |
| `Compound([s1, s2])` | NO | N (independent) | YES | NO -- meshes duplicate at touching faces |

## Decision tree

```
Do the input shapes already share faces (BRep-wise, same Face handle)?
|-- YES: use Glue(...)
|         Fastest. No Boolean runs. Multi-solid preserved.
|         Best for in-memory OCC pipelines.
|
|-- NO, but they touch spatially:
|   Do you need multi-solid output (different materials, etc.)?
|   |-- YES: You cannot get a glued sharing from disjoint inputs in
|   |         one step.  Either:
|   |           (a) rebuild shapes as shared (e.g. extract the shared face
|   |               once and reference it from both sides), or
|   |           (b) mesh each solid alone and stitch the .vol text yourself.
|   |
|   |-- NO (one solid OK): use Fuse(...) with SetGlue(GlueFull) when
|             available, or plain `s1 + s2` for simple shapes.
|
|-- NO, don't want any interaction: Compound(...)
```

## Glue and WriteStep do not commute

NGSolve's `Glue()` hint is an **in-memory** BRep annotation.  Roundtrip
through STEP loses it:

```python
geo = Glue([iron, air])
# ... mesh here: interface shared (22 nodes on a circular interface)
geo.WriteStep("model.step")
geo2 = OCCGeometry("model.step")
# ... mesh here: interface NOT shared (42 nodes, duplicated on both sides)
```

This is an OCCT STEP writer limitation (the STEP format can carry shared
topology, but OCCT's writer emits each solid's faces independently).

**Workarounds** (confirmed with NGSolve dev team, Christopher):

1. Rebuild Glue after load:
   ```python
   geo = Glue(OCCGeometry("model.step").solids)
   ```
   Works only if original solids coincide at matching faces (they usually
   do after OCCT read).

2. Use BREP (OCC native):
   ```python
   # write
   iron.WriteBrep("iron.brep")   # preserves BRep topology
   # read
   from OCC.Core.BRepTools import breptools
   # ... shared faces preserved exactly
   ```
   Caveat: Cubit (ACIS) cannot write BREP.  Only useful when both sides
   are OCC.

3. Stay in memory; export the final mesh directly to `.vol` from
   `OCCGeometry(glue).GenerateMesh(...)`.

## Glue is NOT a solution for ACIS-derived STEP

```python
# Cubit STEP (from ACIS) has duplicated faces at every interface.
# Glue on OCC-loaded-STEP is a no-op: no existing sharing to "keep".
geo = Glue(OCCGeometry("cubit_export.step").solids)
# ^^ interfaces are still duplicated; meshing will have 42 nodes.
```

Cubit's `imprint+merge` creates ACIS-kernel sharing, but OCCT's STEP
reader does not propagate that sharing into BRep identity.  For ACIS
origins, the **correct pipeline is Cubit -> export netgen -> .vol**
(no STEP round-trip, no OCC load, `.vol` text preserves sharing natively).

## Correspondence with Cubit

| OCC / netgen.occ | Cubit (ACIS) | Semantic |
|------------------|--------------|----------|
| `Glue([a, b])` | `imprint all; merge all; compress` (no unite) | multi-solid with shared interfaces |
| `Fuse(a, b)` / `a + b` | `unite volume all` | single solid |
| `Compound([a, b])` | (no direct equivalent; disjoint bodies) | independent meshes |

But note: Cubit `unite` on complex lofted geometry (3turnCoil class)
produces UV-degenerate merged faces, which break downstream high-order
curving.  This is an ACIS-specific failure mode; OCC `Fuse` with
`SetGlue(GlueFull)` handles the same cases cleanly, because OCCT is
designed around BRep topology first.

## Correct pattern for multi-material FEM (e.g. coil + air + iron)

```python
from netgen.occ import *
coil_core = Cylinder(...).Face().Revolve(...)  # or from build
iron_yoke = Box(...) - coil_core                 # carve hole
air_domain = Sphere(...) - iron_yoke - coil_core  # exterior

# All three must share the interfaces they were carved along.
# Since we built them by subtraction, the faces ARE shared BRep-wise.
coil_core.solids.name = "coil"
iron_yoke.solids.name = "iron"
air_domain.solids.name = "air"

geo = Glue([coil_core, iron_yoke, air_domain])
mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=0.01))

# Result: dx("coil"), dx("iron"), dx("air") all work; Kelvin interfaces
# are conformal; mu * curl(u) * curl(v) integrates cleanly over each.
```

## Pitfall: "I used Fuse by accident"

```python
# WRONG (seen often in tutorials that don't need materials):
geo = iron + air
# This is a Fuse, not a Glue.  The interface between iron and air is
# consumed.  dx("iron") will cover the whole domain or raise.

# CORRECT:
geo = Glue([iron, air])
```

The `+` operator is Fuse, not Glue.  No way around it.  Always type the
`Glue(...)` call explicitly.

## Pitfall: Glue on INDEPENDENT solids is a no-op

`Glue(...)` annotates BRep topology that already exists.  It cannot
create sharing from scratch.  Two solids that were built independently
(e.g. `Cylinder(...)` + `Box(...)` without carve) have disjoint Face
handles -- calling `Glue` on them does nothing.

```python
# WRONG: disjoint solids, Glue is a no-op
c = Cylinder(Pnt(0,0,0), Z, r=0.5, h=1.0)
b = Box(Pnt(-1,-1,0), Pnt(1,1,1))
# c and b are constructed independently; their touching face is
# represented by TWO separate Face handles.
geo = Glue([c, b])
# The interface is NOT shared. Mesh will have duplicated nodes.

# CORRECT: carve to establish BRep sharing
air = b - c                        # now "air" has a carved face
geo = Glue([c, air])
# c.faces contains one cylinder lateral face;
# air.faces contains the carved hole referencing the SAME face.
# Glue sees the shared handle and keeps them merged.
```

**Rule**: the solids passed to `Glue(...)` must have been constructed
so that their shared faces are actually the same OCC `Face` object
(usually via Boolean subtract / intersection).  If you built two solids
and they happen to touch in space, you need to Fuse (with `SetGlue`)
to *create* the sharing; `Glue` alone will not.

## Pitfall: STEP/BREP round-trip and Glue identity

Even when the original solids were correctly carved, serialising and
reloading breaks OCC Face identity:

| Format | Reloaded `Glue([...])` works? |
|--------|-------------------------------|
| In-memory (no serialisation) | YES |
| `.brep` (OCC native) | YES -- Face handles preserved |
| `.step` (OCCT writer) | NO -- each solid gets fresh Face handles |
| `.step` (Cubit-exported) | NO -- same reason + ACIS sharing is already lost upstream |
| `.vol` (netgen) | N/A -- `.vol` already IS a mesh, not a shape |

So: the trampoline via STEP always loses Glue.  BREP is the only
serialisation format that carries it (Christopher's recommendation).

**Practical implication**: if you must serialise AND want shared
interfaces, either:
- Use `.brep` (Python side only; Cubit cannot write BREP)
- Rebuild the carve logic after `OCCGeometry(step).solids`
- Use the Cubit -> `export netgen` -> `.vol` pipeline entirely
  (skips OCC on the Python side)

## See also

- topic `mesh` for general OCC meshing recipes
- topic `pitfalls` for Glue-vs-Fusion basic cases (section 39)
- mcp-server-cubit topic `boolean_policy` for the ACIS side (unite vs
  imprint+merge)
"""


NGSOLVE_CLN_CAUER = r"""
# Cauer Ladder Network (CLN) Validation in NGSolve

CLN is a model-order-reduction technique for eddy current problems
(Kameari et al. 2018, Köster et al. 2021). The full-order eddy current
PDE is approximated by an infinite RL ladder, with each rung representing
an eigenmode of the diffusion operator.

## Two flavors of CLN

| Flavor | Excitation | Method | NGSolve fit |
|---|---|---|---|
| **Köster A-T recursion** | voltage-driven coil | static A-T alternating recursion | YES (Köster 2021 §III.B) |
| **Direct frequency sweep** | any (incl. external B) | solve harmonic eddy current at each ω | YES (standard FEM) |
| **Eigenvalue + Foster sum** | any | analytical eigenfunction expansion | analytical only |

**Key insight**: Köster's A-T recursion **does not directly apply to
"isolated conductor in external B" (Case B)** because there is no
terminal voltage. For Case B validation, use direct frequency sweep
and compare with analytical Foster sum.

## Köster A-T recursion (voltage-driven)

For voltage-driven coil + conductor problems, recursive static problems:

```python
# Stage 0: T_0 with H_s × n BC from Biot-Savart of impressed coil
fes_T = HCurl(mesh, order=2,
              definedon=mesh.Materials("conductor"),
              dirichlet="conductor_surface",
              nograds=True)
# Solve: ⟨σ⁻¹ curl T̃_0, curl v⟩_C = 0
#   with essential BC T̃_0 × n = H_s × n on Γ_NC

# Stage n: T̃_(n+1) with BC T̃_(n+1) × n = -T_n × n  (KEY for orthogonality)
# RHS: ⟨σ⁻¹ curl T̃_(n+1), curl v⟩ = -⟨L_n⁻¹ A_n, curl v⟩

# Stage n: A_n on full domain
fes_A = HCurl(mesh, order=2, dirichlet="outer_box", nograds=True)
# Solve: ⟨μ⁻¹ curl Ã_n, curl w⟩_Ω = R_n ⟨T_n, curl w⟩_Ω

# Each stage: R_n = ⟨σ⁻¹ curl T_n, curl T_n⟩_C
#             L_n = ⟨μ⁻¹ curl A_n, curl A_n⟩_Ω
```

## Element-by-element preconditioning

For each static sub-problem, NGSolve `Preconditioner(a, "local")` gives
a block-Jacobi (element-by-element) preconditioner. Use with CG:

```python
a = BilinearForm(fes_T)
a += sigma_inv * curl(u) * curl(v) * dx("conductor")
c = Preconditioner(a, type="local")
# CG with local preconditioner = element-by-element iteration
solvers.CG(sol=gf.vec, rhs=f.vec, mat=a.mat, pre=c.mat,
           tol=1e-8, maxsteps=10000)
```

Avoids global LU/sparse direct solver. Suitable for large 3D problems.

## Direct frequency sweep validation (Case B: external B field)

For uniform external B_z on isolated cuboid, validate the analytical
Foster sum P(ω)/B_0² = (ω²/2) Re[Y_eq(jω)] by direct FEM:

```python
fes = HCurl(mesh, order=2, dirichlet="conductor_surface",
            complex=True, nograds=True)
u, v = fes.TnT()

A_ext = CoefficientFunction((-B0/2 * y, B0/2 * x, 0))  # uniform B_z gauge

# Harmonic eddy current PDE for A_ind (induced)
a = BilinearForm(fes, symmetric=False)
a += (1/mu0) * curl(u) * curl(v) * dx
a += 1j * omega * sigma * u * v * dx
c = Preconditioner(a, type="local")

f = LinearForm(fes)
f += -1j * omega * sigma * A_ext * v * dx

with TaskManager():
    a.Assemble(); f.Assemble()
gf_Aind = GridFunction(fes)
GMRes(A=a.mat, b=f.vec, pre=c.mat, x=gf_Aind.vec, tol=1e-8)

# Loss
A_tot = gf_Aind + A_ext
loss_density = sigma * (omega**2 / 2.0) * (
    A_tot[0]*Conj(A_tot[0])
    + A_tot[1]*Conj(A_tot[1])
    + A_tot[2]*Conj(A_tot[2]))
P = Integrate(loss_density, mesh).real
```

## Critical pitfalls

### 1. Coordinate origin / gauge choice for A_ext

For uniform B_z field, A_ext can be in any gauge:
- A_ext = (B_0/2)(-y, x, 0) — rotates about origin (0,0,0)
- A_ext = (B_0/2)(-(y-b/2), x-a/2, 0) — rotates about (a/2, b/2)

Both give same B_z = B_0 (gauge invariance of B). However, with
**Dirichlet BC A_ind = 0 on conductor surface**, the loss DEPENDS on
the gauge choice (because A_total = A_ind + A_ext at boundary).

→ **Match coordinate convention between analytical and FEM**:
   if Mathematica derived for cuboid at [0,a]×[0,b]×[0,c], place
   NGSolve cuboid at the same location (NOT centered).

Wrong placement gives factor-8 disagreement in loss.

### 2. Factor of 1/2 in P(ω) ↔ Y_eq(s)

The relation is:
  **P(ω) = (ω²/2) Re[Y_eq(jω)] · |B_0|²**

NOT P = ω² Re[Y_eq] · B_0². The 1/2 comes from time-averaging
|E|² for harmonic excitation: ⟨|E|²⟩_t = |E_phasor|²/2.

### 3. Mesh resolution for skin depth

At high frequency, skin depth δ = √(2/(ωμσ)) becomes small.
For Cu at 1 MHz: δ = 66 µm. Mesh size must satisfy h < δ/4 for
accurate FEM. If h > δ, FEM over-predicts loss (currents can't
concentrate at surface as physics requires).

→ For wide-band validation, use adaptive mesh refinement or
   stop frequency sweep at ω where δ = 4·h.

### 4. complex=True for harmonic problems

NGSolve must use complex FE space for jω terms:
```python
fes = HCurl(mesh, order=2, complex=True, ...)
```
Use GMRes (not CG) for non-Hermitian complex systems.

### 5. nograds=True for HCurl

The curl-curl operator has gradient-field nullspace. Use
`nograds=True` to remove gradient DOFs from the FE space:
```python
fes = HCurl(mesh, order=2, dirichlet="...", nograds=True)
```
Otherwise CG/GMRes can't converge.

## Reference implementations

- **Tanimoto's penalty CLN**: `public-safe curated corpus
  20231211_A_(Penalty)_CLN.ipynb` — A-formulation cylinder example
- **Tanimoto's gauge CLN**: `..._A_gauge_CLN.ipynb` — adds H1 gauge
  potential for div-free A
- **Cuboid Case B validation**: `public-safe curated corpus
  2026_04_01_長方形CLN/ngsolve_validation/cuboid_CaseB_freq_sweep.py`
  — direct frequency sweep, agrees with Mathematica Foster sum
  (1 Hz - 100 kHz, <10% error with N_modes=21 truncation)

## References

- A. Kameari, H. Ebrahimi, K. Sugahara, Y. Shindo, T. Matsuo,
  "Cauer ladder network representation of eddy-current fields...",
  IEEE Trans. Magn. 54(3), 7201804 (2018).
- N. Köster, O. König, O. Bíró, "Proper Generalized Decomposition
  with Cauer Ladder Network Applied to Eddy Current Problems",
  IEEE Trans. Magn. 57(6), 6300904 (2021).
- H. Ebrahimi, K. Sugahara, T. Matsuo, H. Kaimori, A. Kameari,
  "Modal decomposition of 3-D quasi-static Maxwell equations by
  Cauer ladder network representation", IEEE Trans. Magn. 56(3),
  7513004 (2020).

## Validation findings (cuboid + cube, 2026-04-27)

### EBE-only solve does NOT converge for general 3D problems

Bíró (verbal): "EBE preconditioning alone is sufficient."
Niels (verbal): "Sometimes works, sometimes doesn't."

**Empirical**: For a 5x2x20 mm Cu cuboid + air with Köster A-T recursion,
applying `Preconditioner(a, "local")` ONCE (no CG iteration) gives R_0
off by 22 orders of magnitude. The HCurl local block on tetrahedra
shares edges with neighbors, so block-diagonal inverse doesn't capture
global coupling.

**Recommendation**: Always wrap with CG or GMRes when using "local"
preconditioner. EBE-only is suitable for very specific symmetric
problems (1D axisymmetric, quasi-1D) but not general 3D.

### Köster A-T is voltage-driven-coil-only

Köster's Eq. (8) BC `T_tilde_n × n = -T_(n-1) × n on Gamma_NC` plus
T_0 BC `T_0 × n = H_s × n` (Biot-Savart from impressed coil) assume
a coil source. For "isolated conductor in external uniform B field"
(Case B), there's no terminal voltage, so Köster's recursion does
not directly apply. Use:
  - **Direct frequency sweep**: solve harmonic eddy current at each ω
  - **Kameari A-only iterative basis** (original, Tanimoto pattern)

### Stage 0 exact validation (1mm Cu cube, Case A all-Dirichlet)

NGSolve Kameari recursion vs Mathematica analytical:

| Quantity | Mathematica | NGSolve | Match |
|---|---|---|---|
| `R_0` | `1/(σV) = 17.24 Ω` | 17.24 Ω | 4-digit ✓ |
| `L_0` | (need conversion) | 2.56 µH | (different topology) |

R_0 matches exactly because both compute the same Y_eq(0) = σ V_C
(Parseval identity for constant-1 expansion in eigenbasis).

### Higher stages numerically unstable

Kameari's J update `J_(n+1) = J_n - σ A_n / L_n` accumulates roundoff.
For cube at h=100µm, order=1: L_1 came out NEGATIVE (catastrophic
cancellation). Remedies:
  - Higher polynomial order (3+) with refined mesh
  - Explicit Coulomb gauge (Tanimoto's H1 potential pattern)
  - WorkingPrecision = 50 in symbolic computation (Mathematica), then
    transfer to NGSolve as exact rational coefficients

### Verified arithmetic with Mathematica Interval[]

Mathematica `Interval[]` propagation through truncated Foster sum
gives provable bounds on Y_eq(0):
```
Y(0) Parseval exact = 0.05800 S·m²
Y(0) Interval sum (N=21^3) = [0.05486, 0.05486]
Truncation deficit: 5.42% with N=21 truncation
```
Combine with NGSolve discretization for end-to-end verified Cauer
ladder values.

## EBE compatibility for div A = 0 / div T = 0 enforcement

Köster orthogonality theorem requires div(curl T) = 0 i.e. T must be
div-free. Standard methods to enforce:

| Method | EBE compatible? | Cost | Exactness |
|---|---|---|---|
| Tokumasu penalty (γ‖div u‖²) | YES (same matrix) | +1 term | Approx (γ tuning) |
| Helmholtz-Hodge projection | NO (extra H1 Poisson) | +1 global solve | Exact |
| Tree-cotree gauge | YES | O(N) preprocess | Exact |
| Gram-Schmidt re-orthogonalization | YES (inner products only) | O(N²) inner products | Recovers orthogonality but residual noise |
| `nograds=True` (NGSolve) | YES | None | Removes high-order grads only |
| A-V mixed with Lagrange multiplier | NO (mixed system) | 2x DoF | Exact |
| PINVIT direct eigensolve | YES (Krylov + EBE pre) | similar to CG | Exact (with order ≥ 3) |

### Empirical findings (1mm Cu cube validation)

- **Stage 0**: R_0 = 17.24 Ω with EBE+CG matches Mathematica
  Y_eq(0) = σV exactly to 4 digits.
- **Stage 1+**: Kameari iteration `J_(n+1) = J_n - σA_n/L_n` suffers
  from cancellation error. Result: L_1 < 0 (unphysical), L grows
  by 10^20 per stage.
- **Gram-Schmidt fix attempt**: Restores orthogonality (GS coef ~ 1e-9)
  but residual J vector itself becomes noise-dominated, breaking
  subsequent stages.
- **PINVIT eigensolve at order=1**: Smallest eigenvalues are gradient
  kernel pollution (1e-4 instead of physical 4e5). Needs order≥3 +
  nograds=True to eliminate.

### Recommendation

For **practical robust validation** of CLN circuit constants:
1. Use NGSolve `HCurl(order=3, nograds=True, dirichlet="all_boundary")`
2. PINVIT for direct eigenvalue/eigenvector extraction
3. Compute R_n = 1/(σ|β_n|²), L_n = μ/(λ_n²|β_n|²) per mode
4. This is fully EBE compatible (PINVIT + EBE preconditioner)

For **strict-EBE Kameari iteration**: tree-cotree gauge required
for higher-stage stability. Implementation: build spanning tree of
mesh edges, mask tree edge DoFs to zero. Non-trivial in 3D.

For **simple robust** (accepting non-EBE): Tanimoto's gauge_CLN
pattern (Helmholtz-Hodge auxiliary H1 Poisson per stage) works well.

## Tree-cotree gauge: WINNING APPROACH (verified 2026-04-27)

**Problem**: `nograds=True` removes only HIGHER-ORDER (p≥2) gradient
bubbles. Lowest-order vertex gradients ∇(linear hat) remain in the
HCurl basis → curl-curl operator has gradient kernel → Kameari
iteration explodes at Stage 1+ (L_n flips negative, grows to 10^21).

**Solution**: Tree-cotree gauge via mesh edge spanning tree.

| Method | EBE pure? | Stable stages | Cost/stage | Comments |
|---|---|---|---|---|
| Plain Kameari | YES | 0 | 1× | Fails immediately |
| Gram-Schmidt re-orthog | YES | 0 | 1× + ⟨,⟩ | Restores orthogonality but residual = noise |
| Helmholtz-Hodge | NO | 3 | **2×** | Tanimoto's gauge_CLN.ipynb pattern |
| **Tree-cotree gauge** | **YES** | **5+** | **1×** + 1-time BFS | **WINNER** |

### Implementation

```python
from ngsolve import HCurl, BitArray
from collections import deque

def build_spanning_tree(mesh):
    # BFS spanning tree of mesh edge graph
    nv = mesh.nv
    visited = [False] * nv
    tree_edges = []
    adj = [[] for _ in range(nv)]
    for ed in mesh.edges:
        v0, v1 = ed.vertices[0].nr, ed.vertices[1].nr
        adj[v0].append((v1, ed.nr))
        adj[v1].append((v0, ed.nr))
    visited[0] = True
    queue = deque([0])
    while queue:
        v = queue.popleft()
        for vn, edn in adj[v]:
            if not visited[vn]:
                visited[vn] = True
                tree_edges.append(edn)
                queue.append(vn)
    return tree_edges  # ≈ nv-1 edges

# ONE-TIME setup
fes = HCurl(mesh, order=2, dirichlet="conductor_surface", nograds=True)
free = BitArray(fes.FreeDofs())
for edge_nr in build_spanning_tree(mesh):
    edge = mesh.edges[edge_nr]
    dofs = fes.GetDofNrs(edge)
    if dofs and free[dofs[0]]:  # lowest-order DoF on this edge
        free[dofs[0]] = False    # tree edge = essential 0

# Each Kameari stage: curl-curl in cotree subspace + ACCUMULATED Apotential
# (Tanimoto pattern, validated against analytic Cauer-I for both circular and
#  rectangular cross sections — see public-safe curated corpus)
Apot = None
for n in range(N_STAGES):
    a = BilinearForm(fes); a += (1/mu) * curl(u) * curl(v) * dx
    pre = Preconditioner(a, "local")  # EBE
    f = LinearForm(fes); f += J * v * dx
    a.Assemble(); f.Assemble()
    inv = a.mat.Inverse(freedofs=free, inverse="sparsecholesky")
    gfA.vec.data = inv * f.vec
    R_n = 1/Integrate(J*J/sigma * dx, mesh)
    Apot = R_n*gfA if Apot is None else Apot + R_n*gfA  # accumulated, weighted by R
    L_n = R_n * Integrate(J * Apot * dx, mesh)          # = R_n * <J_n, Σ R_k A_k>
    J = J - sigma * Apot / L_n                           # subtract accumulated, not bare A_n
```

### CRITICAL: L formula — accumulated A, NOT bare A

A common mistake (which I made initially!) is to compute
`L_n = R_n * <J_n, A_n>` and update `J -= σ A_n/L_n` using only the **current**
stage's potential A_n. This **does not** reproduce the analytic Cauer-I ladder
of Z(s). The correct Tanimoto formula uses the **R-weighted accumulated**
potential `Apot_n = Σ_(k≤n) R_k A_k`:

```
L_n = R_n · ⟨J_n, Apot_n⟩      # Tanimoto, validated
J_(n+1) = J_n - σ · Apot_n / L_n   # subtract accumulated potential
```

Empirical verification (5×2×1 mm Cu, Stage 0):
- Wrong formula `L = R·⟨J,A⟩`: L_0 = 4.67 µH (3D) / 18.16 nH/m (2D) — does
  NOT match Mathematica Cauer-I analytic
- Correct formula `L = R²·⟨J,A⟩` (= Tanimoto for stage 0): L_0 = 8.05 µH /
  31.34 nH/m — **matches Mathematica analytic limit μ·Σ(β²/λ²)/V² to 0.05%**

For stage 1+ the accumulation matters because J_(n+1) needs to be subtracted
against the cumulative response, not just the latest mode.

### Why it works

`nograds=True` + tree-cotree mask = **fully div-free HCurl basis at all orders**:
- High-order grad bubbles: removed by `nograds=True`
- Lowest-order vertex grads: removed by tree edge masking
- Spanning tree edges: each represents one ∇(vertex hat) direction
- Cotree edges: span the orthogonal complement = div-free subspace
- curl-curl on cotree is **strictly positive definite** → CG/direct converges fast

### Geometry caveat: a ≠ b ≠ c required

For cube (a=b=c), eigenvalues λ²(mx,my,mz) = (mx²+my²+mz²)π²/a² have
heavy degeneracy: (1,1,3)=(1,3,1)=(3,1,1) all give same λ². Kameari/PINVIT
cannot distinguish degenerate modes → mode mixing → bad convergence.

**Use a ≠ b ≠ c** (e.g., 5×2×1 mm cuboid) for clean validation.

### Verified result — 2D rectangular bar, 5×2 mm Cu (Case A, per unit length)

With the **corrected Tanimoto formula** (accumulated Apotential):

- **Stage 0**: NGSolve R_0 = 1.7241×10⁻³ Ω/m = 1/(σab) analytic — exact match
- **Stage 0**: NGSolve L_0 = 31.34 nH/m = Mathematica analytic μ·Σ(β²/λ²)/(ab)²
  = 31.32 nH/m — **0.05% match**
- **Stages 0–11**: all positive R_n, L_n with Tanimoto-pattern J update

Use this as the canonical CLN validation case. The 3D 5×2×1 mm cuboid case
in `cuboid_521_treecotree_extended.py` was using the **incorrect** L formula
(no accumulation) — needs to be re-run with the Tanimoto pattern for proper
coefficient-by-coefficient comparison against analytic Cauer-I.

### NGSolve `CreateGradient`: building block

```python
G_matrix, fes_H1 = fes_HCurl.CreateGradient()
# G is sparse (HCurl_ndof × H1_ndof) discrete gradient operator
# Image of G in HCurl = gradient subspace = curl-curl kernel
# tree-cotree picks one HCurl edge per H1 vertex (column of G)
```

### Reference files

- `public-safe curated corpus`
- `public-safe curated corpus`
"""


NGSOLVE_MULTIPHYSICS = r"""
# Multiphysics couplings in NGSolve (COMSOL-class problems)

Reusable building blocks live in ``radia_mcp.radia_ngsolve.multiphysics``;
worked, validated demos.

## Induction heating (magneto-thermal, one-way EM -> heat)

The canonical COMSOL induction-heating tutorial, reproduced and validated.
Pipeline:

    eddy solve (A-Phi harmonic)  ->  joule_loss_density()  ->  solve_heat_steady()

```python
from radia_mcp.radia_ngsolve.solve import solve_eddy_current_harmonic_APhi
from radia_mcp.radia_ngsolve.multiphysics import joule_loss_density, solve_heat_steady

# 1) EM: scattered harmonic eddy in a conductor, axial background B0
A0 = CoefficientFunction((-B0*y/2, B0*x/2, 0))          # curl A0 = B0 z_hat
gfAPhi = solve_eddy_current_harmonic_APhi(mesh, nu, sigma_cf, omega, add_source, ...)
gfA, gfPhi = gfAPhi.components

# 2) Joule loss density q = 1/2 sigma |E|^2  (W/m^3)
q = joule_loss_density(gfA, gfPhi, sigma_cf, omega, A0=A0)   # <-- A0 is essential, see below

# 3) Steady heat -div(k grad T) = q on the conductor, T=0 on its surface
gT = solve_heat_steady(mesh, q, k_thermal, conductor="cyl", dirichlet="cyl_surf")
```

### THE gotcha: scattered-field E must add back A0 (cost ~10x error)

In the SCATTERED eddy formulation the source is a background potential A0
(curl A0 = applied B0) and the solved ``gfA`` is the SCATTERED potential --
that is why field probes write ``curl(gfA) + B0``. By the same token the TOTAL
electric field is

    E = -j omega (A0 + A_scattered + grad Phi)

NOT ``-j omega (gfA + grad Phi)``. Omitting A0 overestimates the Joule loss by
about an order of magnitude. (The TEAM-21 eddy-loss tests use a COIL source, so
their ``gfA`` is already the TOTAL field and they correctly use A0=None -- do not
copy that pattern for background-field problems.) ``joule_loss_density`` takes
``A0=`` precisely so this is not forgotten.

### Validation

Non-magnetic cylinder (R=20mm, sigma=1.6e7, f=500Hz, B0=0.1T axial, R/delta=3.6)
vs closed forms:
  * P/L (central section) vs exact cylinder eddy loss (complex I0/I1 Bessel,
    power-series -- scipy.special.iv is unreliable for complex args): 8.9 %
  * T_centre vs the 1-D radial heat equation driven by the same q(r): 9.6 %
Compare the *central* section's P/L (not the finite-cylinder total) against the
infinite-cylinder analytic -- the caps carry end-effect power.

## Rotating-machine / motor torque (air-gap Maxwell stress)

Electromagnetic torque on a rotor via the weighted Maxwell-stress (eggshell)
band, ``radia_ngsolve.force.eggshell_torque`` (3D) / ``eggshell_torque_2d``:

```python
from radia_mcp.radia_ngsolve.force import eggshell_torque
Tx, Ty, Tz = eggshell_torque(B, mesh, center, r_inner, r_outer, pivot=(0,0,0))
# band r_inner<|r-center|<r_outer must lie in the AIR GAP enclosing the rotor.
```

Core physics validated: a HARD PM
rotor (Br along x) in a uniform applied field B0 at angle theta gives the
torque-angle characteristic tau_z = m B0 sin(theta), m=(Br/mu0)*V, reproduced
to 1-3 % at theta=15..90 deg.

KEY efficiency: a HARD magnet (mu_r=1) does not respond to the applied field,
so fields SUPERPOSE -- ONE magnetostatic solve of the PM gives B_PM, then the
whole angle sweep is B_total = B_PM + B0(theta) with no re-solve. (For a soft
or saturating rotor this breaks; re-solve per angle.) Slotted stators / windings
/ cogging extend the same air-gap-torque method; the validated piece is the
torque operator itself.

## Magnetic shielding (high-mu shell, scattered-field uniform applied field)

Shielding factor S = B_applied/B_cavity of a permeable shell in a uniform field --
the COMSOL AC/DC "magnetic shielding" benchmark. Driven by the PERMEABILITY JUMP
(no coil): ``solve_scattered_uniform_field`` solves curl(nu curl A~) =
curl(nu_pert B0) with nu_pert = nu0(1-1/mu_r) in the shell; total B = curl(A~)+B0.

```python
from radia_mcp.radia_ngsolve.solve import (solve_scattered_uniform_field,
                                           shell_shielding_factor)
# 3 regions: cavity (r<a) / shell (a<r<b, the only permeable mat) / air
gfA  = solve_scattered_uniform_field(mesh, {"shell": mu_r}, (0, 0, B0))
B_in = Integrate((curl(gfA)[2] + B0) * dx("cavity"), mesh) / vol_cavity  # VOLUME AVERAGE
S_fem   = B0 / B_in
S_exact = shell_shielding_factor(mu_r, a, b, "sphere")                   # closed form
```

Exact (sphere): S = [(2 mu_r+1)(mu_r+2) - 2(a/b)^3(mu_r-1)^2]/(9 mu_r);
(cylinder, transverse): S = [(mu_r+1)^2 - (a/b)^2(mu_r-1)^2]/(4 mu_r). Both -> 1 at
mu_r=1. Validated: mu_r=50/100/200,
a/b=2/3 -> S=8.5/16.3/32.0, FEM matches to **0.9 %**. KEY: read the cavity field as
a VOLUME AVERAGE, never a point eval -- at high mu_r the interior is a small residual
of B0 + curl(A~) (near cancellation), so point values are noisy. Same scattered-field
machinery as the permeable-sphere demag check (N=1/3).

## Halbach cylinder (PM dipole, 2D planar)

Ideal dipole Halbach PM cylinder -- the accelerator/undulator workhorse. Remanence Br
with the easy axis rotating as 2*phi over the annulus Ri<r<Ro gives a UNIFORM bore
field B = Br ln(Ro/Ri) (and ZERO field outside, so the box BC is nearly exact).

```python
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic,
    halbach_dipole_magnetization, halbach_bore_field, NU0)
M_cf = halbach_dipole_magnetization(mesh, Br, region="magnet")   # rotating in-plane M
gfu  = solve_planar_magnetostatic(mesh, NU0, magnet_cf=M_cf, order=3)   # 2D A_z
B     = CF((grad(gfu)[1], -grad(gfu)[0]))
B_ref = halbach_bore_field(Br, Ri, Ro)        # = Br ln(Ro/Ri)
```

The new ``magnet_cf`` arg of solve_planar_magnetostatic accepts ANY spatially-varying
in-plane magnetization s = nu0*Br (source += (Mx*dv/dy - My*dv/dx)). cos2phi/sin2phi
are built as (x^2-y^2)/r^2, 2xy/r^2 -- no atan2 branch cut. Validated:
Br=1.2, Ri/Ro=20/40 mm -> 0.832 T uniform
to 0.3 % across the bore. A finite-LENGTH 3D Halbach has end effects (B < the 2D ln
formula); this exact result is the 2D / infinite-length limit.

## Roadmap (same pattern, more couplings)
  * magneto-mechanical (Maxwell-stress body force -> linear elasticity)
  * temperature-dependent sigma(T) -> two-way EM<->thermal Picard loop
  * Joule + convective/radiative surface BC (Robin) instead of fixed T.
  * slotted PM machine: cogging torque vs angle (re-solve per step).
"""


NGSOLVE_ELECTROSTATICS_3D = r"""
# 3D electrostatics & capacitance (COMSOL AC/DC "Computing Capacitance")

The 2D/axisymmetric electrostatics is in ``scalar_fem2d`` (FEMM csolv analog). The
genuine 3D capability is ``radia_ngsolve.electrostatic3d``: conductors as Dirichlet
boundaries, capacitance from the field ENERGY (the route COMSOL uses for the
capacitance matrix).

```python
from radia_mcp.radia_ngsolve.electrostatic3d import (
    solve_electrostatic_3d, capacitance_from_energy, spherical_capacitor_C, EPS0)
from radia_mcp.radia_ngsolve.electrostatics import spherical_capacitor_energy_summary
# conductors = named boundaries; eps_cf = EPS0 (or mesh.MaterialCF for dielectrics)
gfV = solve_electrostatic_3d(mesh, EPS0, {"inner": V0, "outer": 0.0}, order=2)
C   = capacitance_from_energy(gfV, EPS0, V0, mesh)      # C = 2 W / V0^2
row = spherical_capacitor_energy_summary(er, a, b, V0, sample_radius=(a*b)**0.5)
```

-div(eps grad V) = rho ;  E = -grad V ;  W = 1/2 integral eps |grad V|^2 ;  C = 2W/V^2.

GOTCHA (baked into solve_electrostatic_3d): ``gf.Set(cf, definedon=region)`` ZEROS the
whole vector first, so setting conductor potentials in a per-boundary LOOP wipes each
previous boundary (energy -> 0). Set ALL Dirichlet potentials in ONE call via
``mesh.BoundaryCF({"inner":V0,"outer":0}, default=0)``.

Validated: spherical capacitor
C = 4 pi eps0 a b/(b-a), FEM matches to <0.2 % over a gap-ratio sweep. DIELECTRICS
(eps_cf = mesh.MaterialCF({mat: eps0*eps_r})): layered spherical capacitor
C = 4 pi eps0 / [(1/er1)(1/a-1/c) + (1/er2)(1/c-1/b)] matches to 0.1 %
(capacitance_dielectric.py). ISOLATED
conductor self-capacitance: put the ground boundary far away (C_sphere -> 4 pi eps0 R
as the box -> infinity).

Spherical-capacitor JSON gate: ``spherical_capacitor_energy_summary`` records
``C=4*pi*eps*a*b/(b-a)``, ``W=0.5*C*V^2``, ``E(r)=V*a*b/((b-a)r^2)``, and
``p_inner/p_outer=(b/a)^4``.  Use it as a compact 3D field/energy check before
trusting mesh-derived capacitance rows or notebook results.

Capacitance MATRIX: ``capacitance_matrix(mesh, eps_cf, [names], order)`` energises each
conductor to 1 V (others grounded) and reads EVERY conductor's charge from the FEM
REACTION  Q_i = <K V, chi_i>  (chi_i = the FE function 1 on conductor i) -- more
accurate/robust than surface-integrating eps dV/dn, and it comes out exactly symmetric.
C[i,j]=Q_i when V_j=1; C_ii>0, C_ij<0, each row sums to the conductor's capacitance to
ground (0 for a closed system). Validated:
closed spherical capacitor -> C = [[C0,-C0],[-C0,C0]] to 0.03 % with 0 % asymmetry and
0 % row-sum (charge conservation).

Electrostatic FORCE (MEMS actuator): ``force.electrostatic_eggshell_force(E, mesh,
gradg)`` -- the electric twin of the magnetic eggshell, F_k = -int_air[eps0 E_k(E.gradg)
- eps0/2 |E|^2 d_k g]. ``E = -grad(gfV)``; ``gradg`` is grad of a smooth weight band in
air (axis-aligned ramp across a gap, or spherical for a compact body). Validated:
a parallel-plate capacitor with NEUMANN
side walls (exact 1-D field) gives |F| = eps0 A V^2/(2 d^2), matched to 0.4 %, attractive.
"""


NGSOLVE_FIELD_QUALITY = r"""
# Multipole / FIELD-QUALITY analysis of a 2D magnet (accelerator post-processing)

What an accelerator-magnet designer runs after every solve: decompose B on a
reference circle into NORMAL/SKEW multipoles (b_n, a_n), then read the main
harmonic and the field errors (allowed + forbidden harmonics, in 1e-4 units).
Module ``radia_ngsolve.fieldquality``.

Convention (CERN / "European", used in beam dynamics):

    B_y(z) + i B_x(z) = sum_{n>=1} (b_n + i a_n) (z/R_ref)^{n-1},  z = x + i y

so n=1 dipole, n=2 quadrupole, n=3 sextupole; b_n NORMAL, a_n SKEW, each [T]
(the field that pole produces at R_ref). (The US convention shifts the index by
one -- their n-1 is this n.)

```python
from radia_mcp.radia_ngsolve.fieldquality import (
    multipole_coefficients, line_current_multipoles, superpose_multipoles, field_errors)
B = CF((grad(gfu)[1], -grad(gfu)[0]))                       # planar A_z flux density
mp = multipole_coefficients(B, R_ref, n_max=10, mesh=mesh)  # samples |z|=R_ref, DFT
fe = field_errors(mp)        # fe['main_n']==2 for a quad; fe['rel_C'][n] in 1e-4 units
```

The extractor is a DFT: sample f_k = By + i Bx around the circle, then
(b_n + i a_n) is the (n-1)-th bin. ``multipole_coefficients`` takes either an
NGSolve CF (+ ``mesh``, sampled via ``B(mesh(px,py))``) OR a python callable
``phi->(Bx,By)`` (so it works with no FEM, e.g. on measured field).

EXACT reference (no FEM, no convention ambiguity): a line current ``I`` at
``z0 = r0 e^{i phi0}`` gives ``By + i Bx = (mu0 I)/(2 pi)/(z - z0)``, whose
expansion about the origin (R_ref < r0) is
``b_n + i a_n = -(mu0 I)/(2 pi) R_ref^{n-1}/z0^n`` -- ``line_current_multipoles``.
``superpose_multipoles`` sums several wires -> the exact free-space content of any
current arrangement (the analytic the FEM is held to). A current ON the +x axis is
purely NORMAL (a_n=0); rotating the magnet rotates normal<->skew.

ENGINEERING GOTCHAS:
  * Current SOURCE normalization. A wire region driven by ``Jz = I/(pi rw^2)``
    delivers ``J * (discrete material area)``, and the meshed area != pi*rw^2 on a
    coarse mesh -> the MAIN harmonic comes out several % low. Fix: normalize on
    the ACTUAL area, ``J_k = I / Integrate(MaterialCF({wire:1}) * dx)`` so
    ``int Jz = I`` exactly. (By the exterior theorem the bore field depends only on
    each wire's TOTAL current, not its shape, so finite wires == line currents for
    R_ref < r0 -- only the net current must be right.)
  * R_ref must lie in a SOURCE-FREE annulus (no currents/iron between the centre
    and R_ref); harmonics are only defined there.
  * Use n_samples > 2*n_max (guarded) and a high element order (order 4 here) --
    the harmonics are derivatives of A, so they need the field shape resolved.

ALLOWED harmonics from symmetry: a 2N-pole with the matching 2N-fold sign pattern
keeps only n = N, 3N, 5N, ...  A normal quadrupole (4 wires +,-,+,- at
0/90/180/270 deg) -> n = 2, 6, 10; the 12-pole (n=6) sits at exactly (R_ref/r0)^4
of the main; all "forbidden" n (1,3,4,5,7,...) must vanish. This symmetry table is
the first sanity check on any magnet design.

Validated (tests/test_fieldquality.py,
validation/force/validate_force_xval.py::validate_multipole_quadrupole_fem): the DFT reproduces
line_current_multipoles to ~1e-15; an air-cored FEM quadrupole gives main b2 to
**0.02 %**, the allowed 12-pole at 256 units (== (R_ref/r0)^4) vs 256.00 exact, and
forbidden leakage below **0.3 units** (1e-4) -- numerical noise.
"""


NGSOLVE_PM_AXIAL = r"""
# Cylinder permanent magnet: on-axis field (axisymmetric) + the axi extraction recipe

The PM-design workhorse: on-axis B_z of a uniformly axially-magnetized cylinder
(radius R, length L, remanence Br), the start of every pull/holding-force and gap
estimate. Closed form (surface-charge model, recoil mu_r=1):

    B_z(z) = (Br/2)[ (z+L/2)/sqrt(R^2+(z+L/2)^2) - (z-L/2)/sqrt(R^2+(z-L/2)^2) ]

``solve.cylinder_magnet_axial_field(Br,R,L,z)`` is the closed form; centre value
Br L/(2 sqrt(R^2+(L/2)^2)), far field -> the dipole Br R^2 L/(2 z^3).

```python
from radia_mcp.radia_ngsolve.solve import solve_axi_magnetostatic, cylinder_magnet_axial_field, NU0, MU0
# axisymmetric (r,z): magnet rectangle in air; rigid axial M = Br/mu0 via theta=90 deg
gfu = solve_axi_magnetostatic(mesh, CoefficientFunction(NU0),
                              magnets={"magnet": (Br/MU0, 90.0)}, order=2)
```

AXI EXTRACTION RECIPE (hard-won; reuse for ANY axisymmetric post-processing):
  * The H1Henrotte axi space gives B_z = dA/dr + A/r, RELIABLE INSIDE AN INTEGRAL
    but its GridFunction GRADIENT is NOT point-evaluable (it falls back to a base
    diffop -> wrong/noisy point values; you'll see "called base class apply ...
    AxiHenrotteDiffOpGradient"). So PROJECT the derived field onto a standard space
    first -- ``bzgf = GridFunction(L2(mesh,order=2)); bzgf.Set(grad(gfu)[0]+gfu/x)``
    (the .Set uses the working integral path) -- then point-sample ``bzgf``.
  * On-axis fields: sample JUST OFF axis (r>0), never exactly r=0 (there
    dA/dr and A/r each -> B_z0/2; an exact-axis eval double-counts/zeros). Average a
    few near-axis radii (e.g. R/25..R/16) to wash out the L2 inter-element jump that
    can spike one unlucky sample ~10 %. (Equivalent flux view: 2 A_phi(rho)/rho =
    <B_z>_disk(rho) -> B_z(0) as rho->0.)
  * Mesh: refine the r=0 AXIS edge (``face.edges.Min(X).maxh = ...``) so the first
    element row is finer than the sampling radius, and keep all sample points inside
    a uniformly-fine NEAR region (away from the coarse-far interface) -- the EXTERNAL
    on-axis field is the sensitive part (a too-coarse near gives a uniform ~10 % low).

DEAD ENDS (tried, FAILED -- don't repeat; "learn from failures too"):
  * Point-evaluating ``grad(gfu)`` directly -> wrong/noisy (the diffop fallback).
  * Richardson-extrapolating the disk average ``(4*davg(rho)-davg(2*rho))/3`` ->
    AMPLIFIES the L2 sampling noise, made every point 8-15 % worse.
  * H1 (continuous) global projection instead of L2 -> SPREADS the pole-corner
    singularity error over a wide z-band (z/L=0.75..1.5 got worse). L2 keeps it local.
  The winner is the plain L2 projection + near-axis radial MEAN above.

Validated (validation/force/validate_force_xval.py::
test_cylinder_magnet_axial_field): centre field to **0.01 %**, on-axis profile to a
few % out to z=3L. THE ONE caveat: right at the pole exit (z ~ L) a sharp-cornered
magnet has a charge SINGULARITY at the (R,L/2) corner and the field bulges off-axis,
so the on-axis value there is ~8-10 % sensitive (refine the pole corner, or read the
centre/far field which are robust). Same magnets={mat:(Hc,theta)} source as the
validated magnetized-sphere (B_in = 2 mu0 mu_r Hc/(mu_r+2)); here Hc = Br/mu0, mu_r=1.
"""


NGSOLVE_SOLENOID = r"""
# Finite air-core SOLENOID: on-axis field + self-inductance (axisymmetric)

Air-cored coil design. A solenoid (radius a, length L, n turns/m, current I) is an
azimuthal current sheet K = nI -> model the winding as a thin annulus carrying
J_phi = nI/t, solved with ``solve_axi_magnetostatic(mesh, NU0, Jr=J_phi_cf)``.

```python
from radia_mcp.radia_ngsolve.solve import (solve_axi_magnetostatic,
    solenoid_axial_field, solenoid_inductance_long, NU0)
from radia_mcp.radia_ngsolve.force import inductance_axi
Jr  = mesh.MaterialCF({"coil": n*I/t}, default=0.0)
gfu = solve_axi_magnetostatic(mesh, CF(NU0), Jr=Jr, order=2)
B   = CF((grad(gfu)[0] + gfu/x, -grad(gfu)[1]))          # (B_z, B_r)
L_ind = inductance_axi(B, mesh, NU0, I)                  # 2 W / I^2  (per-turn I !)
```

Two validated outputs (validation/force/validate_force_xval.py::validate_solenoid_field_and_inductance):
  * ON-AXIS B_z(z) = (mu0 nI/2)[(L/2-z)/sqrt((L/2-z)^2+a^2)+(L/2+z)/sqrt((L/2+z)^2+a^2)]
    (``solenoid_axial_field``). Smooth on axis (no pole corner like the bar magnet
    #14), so the L2-project + near-axis radial-mean extraction is clean: centre
    **0.7 %**, whole profile **<1.7 %**. Centre = mu0 nI L/sqrt(L^2+D^2); long limit
    mu0 nI; coil end (z=L/2) ~ half centre.
  * SELF-INDUCTANCE by the field ENERGY (``inductance_axi``, 2W/I^2 -- robust, NO
    axis sampling). Pass the PER-TURN current I (not nI): the source already carries
    the n turns, so 2W/I^2 ~ N^2. Ratio L_ind/L_long (L_long = mu0 n^2 pi a^2 L =
    ``solenoid_inductance_long``) is the Nagaoka coefficient k_N(2a/L) < 1, matched
    to Wheeler 1/(1+0.9 a/L) to **0.5 %**.

KEYS: (1) capture the EXTERNAL flux energy -- a too-small box truncates W and reads L
~3 % low (use a box several coil-lengths out). (2) a thin winding (t/a small) sits on
the thin-sheet Nagaoka value; a thick multilayer winding is a few % below it (real
effect, not error). (3) per-turn current in inductance_axi, as above.
"""


NGSOLVE_BUSBAR_LORENTZ = r"""
# Parallel busbar / conductor Lorentz force (short-circuit / fault force)

The force between two long parallel conductors d apart carrying I1, I2:
``F/L = mu0 I1 I2/(2 pi d)`` (attractive if parallel, repulsive if anti-parallel) --
the force that mechanically wrecks busbars during a short circuit.

```python
from radia_mcp.radia_ngsolve.force import lorentz_force_2d
B = CF((grad(gfu)[1], -grad(gfu)[0]))               # planar A_z flux density
Fx, Fy = lorentz_force_2d(Jz, B, mesh, "w2")        # N/m on conductor w2
```

``lorentz_force_2d(Jz, B, mesh, region)`` = the J x B volume integral
``Fx=-int Jz By, Fy=int Jz Bx`` over the conductor -- the direct current-source twin
of the Maxwell-stress ``eggshell_force_2d``. Use the TOTAL field B: a conductor's own
self-field exerts ZERO net force (symmetry), so it drops out automatically.

KEY: area-normalize the current source ``Jz = I/Integrate(MaterialCF({wire:1})*dx)``
so int Jz = I exactly on the discrete wire (same gotcha as the multipole quad). Use a
fine-enough wire and d >> rw so B1 is ~uniform across wire 2. Validated
(validation/force/validate_force_xval.py::validate_busbar_lorentz_force):
both attractive (parallel) and repulsive (anti-parallel) match mu0 I1 I2/(2 pi d) to
~1 %, correct signs, transverse force ~0.

METHOD OF IMAGES (wire above an iron plane): a wire (current I) a height h above an IRON
(mu->inf) half-plane is attracted with ``F = mu0 I^2/(4 pi h)`` (``image_force_wire_iron``)
-- the iron mirrors a SAME-sign image current at -h (perpendicular field lines), so it's the
two-parallel-currents force at separation 2h. (A perfect CONDUCTOR mirrors an OPPOSITE-sign
image -> equal repulsion.) Same recipe: ``lorentz_force_2d(Jz, B, mesh, "wire")`` (the wire's
self-field nets to zero). With a large-but-FINITE iron block the FE force is ~4 % below the
infinite-half-plane ideal (-> exact as the block grows); validated
tests/test_image_force.py.
"""


NGSOLVE_HELMHOLTZ = r"""
# Helmholtz coil pair: uniform-field generation + uniformity (axisymmetric)

Two coaxial N-turn loops (radius a) at z = +/- a/2 in the SAME sense -> a maximally
UNIFORM field at the centre (d^2B/dz^2 = 0, the Helmholtz condition). The design tool
for uniform fields (NMR shim, sensor/magnetometer calibration, beam optics).

``coil_pair_axial_field(N, I, a, separation, z)`` is the exact on-axis field
``B_z = (mu0 N I a^2/2)[((z-s/2)^2+a^2)^-3/2 + ((z+s/2)^2+a^2)^-3/2]``; Helmholtz =
``separation == a``, centre ``B0 = (4/5)^(3/2) mu0 N I/a = 0.7155 mu0 NI/a``. (Pass
``-I`` to one loop for ANTI-Helmholtz: a linear gradient / magnetic bottle / axial
quadrupole.) Solve with ``solve_axi_magnetostatic(Jr=...)`` (two loop regions), extract
on-axis B_z with the L2-project + near-axis radial-mean recipe.

KEY (uniformity): the figure of merit is the field FLATNESS over the central working
volume -- ``max|B(z)/B0 - 1|`` over ``|z| <= a/4`` is ~0.42 % for an ideal Helmholtz pair.
To measure it from FEM, the centre value must be clean: average the near-axis radii AND a
small +/-dz pair, because the EXACT symmetry plane (z=0) is a sampling-artifact hot spot
(a lone z=0 sample can read ~1.5 % high). Validated
(validation/force/validate_force_xval.py::validate_helmholtz_uniformity):
on-axis B_z <0.5 %, centre 0.38 %, FEM uniformity 0.84 % (= the 0.42 % ideal + sub-%
numerical scatter) -- the field is flat to <1 % across the bore.
"""


NGSOLVE_C_MAGNET = r"""
# Iron-yoke window-frame DIPOLE: air-gap field (magnetic circuit / reluctance)

The accelerator-dipole / actuator / relay workhorse: an iron yoke (mu_r) closes the
coil's flux and forces it across a small air gap. Reluctance model (Ampere loop, B
continuous): ``B_gap = mu0 N I/(gap + iron_path/mu_r)`` -> ``magnetic_circuit_gap_field``;
``mu_r -> inf`` gives ``mu0 NI/gap``.

```python
from radia_mcp.radia_ngsolve.solve import solve_planar_magnetostatic, magnetic_circuit_gap_field
nu = mesh.MaterialCF({"iron": NU0/mu_r}, default=NU0)      # high-mu yoke
gfu = solve_planar_magnetostatic(mesh, nu, Jz=Jz, order=3, dirichlet="outer")
B = CF((grad(gfu)[1], -grad(gfu)[0]))
Bx_gap = Integrate(MaterialCF({"gap":1})*B[0]*dx)/gap_area  # average over the gap
```

BUILD: 2D planar A_z, iron = ``big - window - gap`` (OCC boolean), a coil threading one
leg as +NI inside the window and -NI outside the frame (area-normalize ``Jz = I/Integrate
(MaterialCF({coil:1})*dx)`` so int Jz = NI exactly). Read the gap field as the gap-AVERAGE
of the across-gap B component (not a point). Validated
(validation/force/validate_force_xval.py::validate_c_magnet_gap_field):
g=6 mm, mu_r=2000 -> B_gap matches the reluctance model to **0.5 %**, sitting just below
the mu_r->inf ideal (the deficit = gap fringing/leakage -- the lumped model's blind spot,
and exactly what a COMSOL cross-check resolves). Wider gaps -> more fringing -> bigger
deficit; that trend IS the engineering content (effective gap > geometric gap).

HOLDING FORCE (the force companion -- relay / solenoid / electromagnet pull): the pole faces
attract with the Maxwell-stress holding force ``F = B_gap^2 * A_face/(2 mu0)`` (A_face =
pole-face area = face height * stack depth). ``magnetic_circuit_gap_force(N,I,gap,iron_path,
mu_r,area)`` = mu0(NI)^2 A/[2(gap+iron_path/mu_r)^2] ~ 1/gap^2. The FE force (B_gap_FE^2
A/2mu0 from the solved gap field) matches the reluctance model to **<1.4 %**, and
F*(gap+iron_path/mu_r)^2 is constant to **<1 %** across a 2.25x gap range -- the 1/g^2
reluctance force law. Validated tests/test_reluctance_actuator.py.

If ``B_gap`` is already known from a nonlinear magnetic circuit or an FE solve,
use the dependency-free force helpers directly:
``air_gap_maxwell_pressure(B_T)``, ``air_gap_holding_force(B_T, area_m2, faces=1)``,
and ``air_gap_force_summary``.  `validation_air_gap_force_sweep.py` connects the
nonlinear B-H circuit to this pressure law: for a 96 A-turn example with
`A=2.5e-4 m^2` and two active faces, force falls from `562.442916 N` at closed
gap to `0.723099826 N` at 2 mm; `p(1 T)=397887.35772973835 Pa`; pressure/force
identities are exact to machine precision.

GAPPED-CORE INDUCTANCE (#39) -- the inductance twin of the actuator force, SAME window
frame: ``magnetic_circuit_inductance(N,A,gap,iron_path,mu_r)`` = N^2/R =
``N^2 mu0 A/(gap + iron_path/mu_r)`` (the gap-force effective gap, again; mu_r->inf =>
gap-limited L -> N^2 mu0 A/gap). radia reads the true L with ``inductance_2d`` (lambda/I on
the coil bundles). KEY CONTRAST WITH THE FORCE: the gap FORCE reads only B_gap, so the
leakage-free reluctance model is good to <1.4 %; but the coil INDUCTANCE also LINKS THE
WINDOW LEAKAGE FLUX, so the FE L sits WELL ABOVE the reluctance model -- **~1.27x at
g=4 mm growing to ~1.52x at g=9 mm** (more flux short-cuts the window as the gap opens).
That excess IS the physical content the lumped model misses. Tool-independent gate: the FE
L bounded BELOW by L0 (and < ~1.8 L0), L falling with g, and L/L0 RISING with g; a radia
self-regression pins the values and an independent 2D FE solver matches them to ~0.4 %
(capability parity, recorded in the internal cross-val). Validated
validation_test/radia_mcp/test_gapped_core_inductance.py.

PM LOAD LINE / operating point (#42) -- the PERMANENT-MAGNET-driven member, completing the
magnetic-circuit trilogy (coil force #27, coil inductance #39, PM operating point now). A
magnet (Br, recoil mu_r=1) in one leg drives flux across the gap; Ampere's law + the recoil
line B=Br+mu0 mu_rec H give ``pm_circuit_loadline_gap_field`` =
``Br l_m/(l_m + mu_rec g + l_fe/mu_r)`` (the NI of #27 replaced by the magnet's Br l_m/mu_rec
MMF). A larger gap is a STEEPER load line -> lower B_gap (more demagnetisation). radia drives
it with ``solve_planar_magnetostatic(magnets={"magnet":(Br/mu0, phi_deg)})``, the magnet
magnetised ALONG its leg (phi 90 deg = +y), and reads the gap Bx. LEAKAGE (opposite sign to
the inductance #39): magnet flux short-cuts the window instead of reaching the gap, so the FE
B_gap is BELOW the leakage-free load line, ~0.86x (g=2 mm) -> ~0.68x (g=8 mm), the deficit
GROWING with the gap. Tool-independent gate: FE below the load line (upper bound), B_gap and
B_gap/loadline both FALLING with the gap, + a radia self-regression; an independent 2D FE
solver matches radia to ~0.3 % (capability parity, internal cross-val). Validated
validation_test/radia_mcp/test_pm_loadline.py.

For hot-magnet demagnetisation screening, keep the compact student-facing summary:
``pm_loadline_demag_risk_summary(rows)`` checks that safe load-line rows form a prefix,
records ``safe_prefix_count`` / ``first_unsafe_gap_m`` / ``largest_safe_gap_m``, and labels
the minimum margin as green (>100 kA/m), amber (0..100 kA/m), or red (<0).  This is a
pre-run margin contract only; local per-element field checks are still required for a final
irreversible-demag claim.

Continuous-loop slot 142 adds the RunResult handoff rule for ELF/MAGIC-style notebooks: before
turning local product result rows into a load-line or demag-margin panel, normalize the row into
``case_id``, ``temperature_C``, ``gap_m``, ``B_gap_T``, ``H_pm_A_per_m``,
``H_knee_A_per_m``, ``recoil_mu_r``, and ``safe_against_knee``.  Run
``pm_loadline_metadata_gate`` first, then ``pm_temperature_demag_sweep_summary`` or
``pm_recoil_demag_step_summary``.  A parsed RunResult without H units, demag sign convention,
magnetization axis, knee reference, and recoil permeability is not solver-ready evidence.
"""


NGSOLVE_JOULE_ELECTROTHERMAL = r"""
# Electro-thermal (Joule heating) multiphysics -- chaining the two solvers

The canonical COMSOL "Joule Heating" coupling: current dissipates ohmic heat, which sets
the temperature. radia-ngsolve does it by CHAINING the existing solvers (not hand-imposing
the heat):

    solve_current_flow  ->  joule_heat_source(q = sigma|grad V|^2)  ->  solve_heat_steady

```python
from radia_mcp.radia_ngsolve.scalar_fem2d import solve_current_flow
from radia_mcp.radia_ngsolve.multiphysics import (solve_heat_steady, joule_heat_source,
                                                  joule_bar_temperature_rise)
V  = solve_current_flow(mesh, sigma, {"left": V0, "right": 0.0})   # -div(sigma grad V)=0
q  = joule_heat_source(V, sigma)                                   # sigma |grad V|^2 [W/m^3]
dT = solve_heat_steady(mesh, q, k, conductor="bar", dirichlet="left|right")  # T=0 -> RISE
```

``joule_heat_source`` is the DC twin of ``joule_loss_density`` (AC induction heating, #1).
``solve_heat_steady`` returns the temperature RISE (T=0 on the cooled ``dirichlet``).

Validated (validation/force/validate_force_xval.py::
test_joule_heating_electrothermal): a uniform bar (voltage V over length L, both ends cold,
sides insulated) -> uniform q = sigma(V/L)^2 and the EXACT parabolic rise
``dT(x) = (q/2k) x (L-x)``, peak ``sigma V^2/(8k)`` (length-INDEPENDENT). FEM matches to
**0.00 %** -- order-2 H1 represents the quadratic profile and uniform source exactly. The
clean two-physics chain; for AC/eddy use ``joule_loss_density`` instead of the DC source,
and for temperature-dependent sigma(T)/k(T) wrap the chain in a Picard loop (two-way).

CYLINDRICAL sibling (#41) -- a round conductor (radius a) carrying a uniform DC current density
J: uniform Joule source q=J^2/sigma, surface held at the coolant temperature, steady RADIAL
conduction -div(k grad T)=q gives the parabolic rise dT(r)=q(a^2-r^2)/(4k), centre
``dT_peak = q a^2/(4 k)`` (``multiphysics.joule_heated_cylinder_peak_dT``) -- the classic
current-carrying-wire / ampacity temperature rise, the cylindrical analog of the slab's
sigma V^2/(8k). solve_heat_steady on a disk (face "wire", edge "surf") reproduces dT(r) to
**0.00 %**. VALIDATED LIVE 3-way (the corrected COMSOL policy: `comsol_preflight` -> gate1
LiveLink up + gate2 module licensed [core ht, no add-on] -> live): COMSOL 6.4 LiveLink
**25.0000 K == NGSolve == analytic q a^2/4k, 0.000 %** (COMSOL number internal/unattributed;
public gate is the closed form). tests/test_cylinder_joule_heating.py.

CONVECTIVE cooling (#44) -- the realistic BC. Instead of a fixed surface temperature, the slab
is cooled by CONVECTION (Newton cooling, film coefficient h): a ROBIN boundary
``-k dT/dn = h(T - T_inf)`` so the surface FLOATS to balance the convected heat.
``solve_heat_steady_robin(mesh, q, k, conductor, robin, h, T_inf)`` adds ``h*T*s`` to the
stiffness and ``h*T_inf*s`` to the load on the ``robin`` boundary (other boundaries natural =
insulated). For a slab (thickness L, both faces convective) the centre rise is
``convective_slab_peak_dT`` = ``qL^2/(8k) + qL/(2h)`` = conduction parabola PLUS the convective
FILM rise qL/(2h); the surface sits qL/(2h) above the coolant, and h->inf recovers the fixed-T
parabola qL^2/8k. VALIDATED LIVE 3-way (`comsol_preflight` gates core ht -> live): COMSOL 6.4
LiveLink **22.5000 K == NGSolve == analytic, 0.000 %** (the convective BC = ht `HeatFluxBoundary`
with `HeatFluxType='ConvectiveHeatFlux'`, `h`, `Text`; COMSOL number internal, public gate is the
closed form). tests/test_convective_electrothermal.py.

COOLING FIN / extended surface (#47) -- conduction ALONG a thin fin + distributed surface
convection: the 1-D fin equation, dT(x) = dT_base cosh(m(L-x))/cosh(mL), m = sqrt(2h/(k t))
(``fin_parameter``), adiabatic tip dT_tip = dT_base/cosh(mL) (``fin_tip_temperature_rise``), and
the FIN EFFICIENCY eta = tanh(mL)/(mL) (``fin_efficiency``) = actual heat / heat-if-all-at-base
(eta->1 short/fat/conductive, ->0 long/thin -- the heat-sink design trade-off). radia uses
``solve_heat_steady_robin`` with the fin BASE as an inhomogeneous Dirichlet (``dirichlet_value``)
+ convective long faces; the fin efficiency comes from the convected heat
``Integrate(h*T*ds(conv))/(h*P*L*dT_base)``. VALIDATED LIVE 3-way (`comsol_preflight` gates core
ht -> live): COMSOL 6.4 LiveLink tip **16.3307 K == NGSolve 16.3307 == 1-D analytic 16.307**
(0.14 %, the 2-D-vs-1-D fin approximation), eta FE 0.3945 vs tanh(mL)/(mL) 0.3946 (0.04 %).
COMSOL recipe: base `TemperatureBoundary` + long-face `HeatFluxBoundary`/`ConvectiveHeatFlux`.
Generalises the convective slab (#44, no along-fin gradient). tests/test_cooling_fin.py.

MULTILAYER thermal resistance / DIE-STACK junction temperature (#50) -- the electronics-cooling
junction-to-ambient metric. A heat-generating top layer conducts DOWN through a passive layer
stack to a cooled base; the heat flux Q=q*L1 flows through the SERIES thermal resistances
R_i = L_i/k_i (``thermal_resistance_series`` = sum L_i/k_i, the thermal-circuit analog of series
electrical R). The temperature DROP across passive layer i is Q*L_i/k_i, and within the active
top layer the self-heating drop is q*L1^2/(2 k1), so

    T_junction = q L1^2/(2 k1) + Q*(L2/k2 + L3/k3),   Q = q L1.

Low-k layers (solder, TIM) dominate the rise; high-k spreaders (Cu) drop little. radia uses
``solve_heat_steady`` with a region-wise k CF (``mesh.MaterialCF({lay1:k1,...})``) + the source
only in the active layer. VALIDATED LIVE 3-way (`comsol_preflight` gates core ht -> live): COMSOL
6.4 LiveLink == NGSolve == analytic, **0.000 %** at every interface (the conduction profile is
piecewise-linear/parabolic, represented exactly). COMSOL recipe: per-layer k via
`SolidHeatTransferModel` features (one per domain, `k_mat='userdef'`) + `HeatSource` in the
active layer + base `TemperatureBoundary`. tests/test_die_stack_thermal.py.
"""


NGSOLVE_ELASTICITY = r"""
# Linear ELASTICITY + thermo/electro-mechanical coupling (structural multiphysics)

The structural half of the couplings (thermal stress, magneto/electro-mechanical force ->
deflection). ``radia_ngsolve.elasticity``: 2D small-strain isotropic, plane stress/strain,
with body force, boundary tractions, and an optional thermal pre-strain.

```python
from radia_mcp.radia_ngsolve.elasticity import solve_linear_elasticity, stress_2d
u = solve_linear_elasticity(mesh, E, nu, dirichletx="left", dirichlety="bottom",
                            traction={"right": (s0, 0.0)}, plane="stress")   # uniaxial
sig = stress_2d(u, E, nu, "stress")          # sig[0]=sigma_xx, sig[3]=sigma_yy
```

sigma = 2 mu eps(u) + lam tr(eps) I - 2(mu+lam) alpha dT I; mu=E/2(1+nu), plane-stress
lam=E nu/(1-nu^2), plane-strain lam=E nu/((1+nu)(1-2nu)). Per-component clamps:
``dirichletx``/``dirichlety`` fix u_x/u_y=0 (use dirichletx="left", dirichlety="bottom"
for the statically-determinate uniaxial setup; clamp u_y on a full edge only when the
thermal expansion there is genuinely zero). ``thermal=(alpha, dT)`` adds the thermal load.

Validated EXACT (order-2 captures the linear/quadratic fields):
  * uniaxial tension: u_x(L)=s0 L/E, lateral u_y=-nu s0/E y, sigma_xx=s0 (0.00 %).
  * thermal stress, bar clamped in x both ends, uniform dT: sigma_xx=-E alpha dT (0.00 %).
  * ELECTRO-THERMO-MECHANICAL chain (validation/force/validate_force_xval.py::validate_electro_thermo_mechanical_chain): solve_current_flow ->
    joule_heat_source -> solve_heat_steady -> solve_linear_elasticity. A Joule-heated bar
    bolted at both ends -> peak rise sigma V^2/(8k), then axial stress
    ``constrained_bar_thermal_stress = -E alpha <dT>`` (<dT>=2/3 dT_max for the parabola),
    both to 0.00 %. The same boundaries are electrode = heat sink = clamp.

MAGNETO-MECHANICAL (validation/force/validate_force_xval.py::validate_magneto_mechanical_beam): a current-carrying cantilever in
a transverse field B0 -> Lorentz body force ``f_y = -J_x B0`` -> deflection. Tip matches
Euler-Bernoulli ``cantilever_tip_deflection = w L^4/(8EI)`` (w = I B0, I = h^3/12) to
**0.02 %** for a slender beam (L/h=25). The magnetic-actuator / loudspeaker / galvanometer
principle. For a distributed Maxwell-stress load use the ``force.eggshell_*`` body force;
for MEMS, the electrostatic force. Two-way (large deflection / moving mesh) needs a
Picard/ALE loop on top.

THERMAL BENDING (bimorph / bimetallic actuator): a LINEAR through-thickness gradient
T(y)=dT*(y/h-1/2) on a free cantilever bends it STRESS-FREE into constant curvature
kappa=alpha*dT/h; tip delta=kappa L^2/2 = alpha*dT*L^2/(2h) (``thermal_bending_tip_deflection``).
Pass the y-linear field as ``thermal=(alpha, dT*(y/h-1/2))`` and clamp one end; FE matches beam
theory to **0.24 %** (tests/test_thermal_bending.py).
Distinct from the AXIAL blocked-bar stress above: a through-thickness GRADIENT -> curvature
(stress-free), whereas a uniform/parabolic dT in a CLAMPED bar -> stress (no bending).

BIMETALLIC STRIP (thermostat / thermal-switch): TWO bonded layers of equal thickness/modulus
but different expansion (alpha1, alpha2) under a UNIFORM rise dT curve with Timoshenko (1925)
``kappa = 3(alpha2-alpha1)dT/(2h)`` (h=total thickness); cantilever tip delta=kappa L^2/2
(``bimetal_tip_deflection``). Pass alpha as a region-wise CF ``mesh.MaterialCF({lay1:a1,
lay2:a2})`` to ``solve_linear_elasticity(thermal=(alpha, dT))``. FE vs Timoshenko **1.5 %**
(tests/test_bimetallic_strip.py). The TWO-MATERIAL
Delta-alpha bending -- vs thermal_bending's ONE-material through-thickness gradient.

BIAXIAL thermal stress (fully in-plane-constrained plate, uniform dT): u=0 is exact, so
sigma_x=sigma_y = -E alpha dT/(1-nu) (PLANE STRESS, ``biaxial_thermal_stress``) -- the 2D
constraint factor 1/(1-nu) vs the 1D uniaxial -E alpha dT (constrained_bar_thermal_stress).
Validated to **0.000 %** (tests/test_biaxial_thermal_stress.py). KNOWN BUG (found here): radia's plane='strain' thermal
load is (1+nu) TOO SMALL -- a fully-constrained plane-strain plate should give -E alpha dT/(1-2nu)
but yields -E alpha dT/((1+nu)(1-2nu)). Plane-strain thermal problems need the effective
alpha* = (1+nu) alpha (or the thermal coefficient (3 lam + 2 mu) instead of 2(mu+lam)); flagged
for a separate fix. Plane STRESS is correct; use plane='stress' for thermal stress until fixed.

LAMINATE / COMPOSITE RESIDUAL stress (FREE bonded bar, uniform dT, NO external constraint):
the layers are forced to a common strain by the bond, so a CTE mismatch alone builds internal
stress. Stiffness-weighted rule of mixtures (``rule_of_mixtures_cte``):
``alpha_eff = sum(E_i A_i alpha_i)/sum(E_i A_i)``; per-layer ``sigma_i = E_i (alpha_eff - alpha_i)
dT`` (``laminate_residual_stress``) -- high-alpha layers COMPRESSED, low-alpha in tension, the
section self-equilibrated sum(A_i sigma_i)=0. The FREE-bar counterpart of the externally blocked
``constrained_bar_thermal_stress``; the coating / die-attach / composite-laminate workhorse.
Build a SYMMETRIC stack (e.g. face|core|face) so it does NOT bend (asymmetric -> bimetal curl,
above), make it long (L>>h) and read sigma_xx in the uniform CENTRAL window (Saint-Venant away
from the ends), passing region-wise E and alpha CFs to solve_linear_elasticity(thermal=(alpha,dT)).
Validated to **0.000 %** (sigma per layer AND the common central strain = alpha_eff*dT)
(tests/test_laminate_thermal_stress.py).
"""


NGSOLVE_RELUCTANCE_TORQUE = r"""
# Cogging / reluctance TORQUE by rotor-sweep -- validated against a regression reference

The classic rotating-machine workflow: ROTATE the rotor and read the air-gap
Maxwell-stress torque at each angle. radia-ngsolve does it with EXISTING capability
and NO new solver code -- rotate the magnet's ``phi_deg`` (no remesh) and integrate
``eggshell_torque_2d`` (weighted Maxwell stress, Henrotte) over an air-gap band:

```python
from radia_mcp.radia_ngsolve.solve import reluctivity, solve_planar_magnetostatic
from radia_mcp.radia_ngsolve.force import eggshell_torque_2d
nu = reluctivity(mesh, {"iron": 1000.0})
for th in thetas:                                   # "rotate" = re-aim the PM, no remesh
    A = solve_planar_magnetostatic(mesh, nu, magnets={"rotor": (Br/MU0, th)}, order=3)
    B = CF((grad(A)[1], -grad(A)[0]))
    tau_mm = eggshell_torque_2d(B, mesh, (0,0), r_in, r_out)   # band INSIDE the air gap
```

UNITS: meshed in mm, ``eggshell_torque_2d`` returns a 2D stress-torque in mm-units. To
the N*m of a depth-``d`` stack, the 2D torque scales as length^2 (lever-arm * area):
``tau[N*m] = tau_mm * 1e-6 * d_m`` (mm^2->m^2 = 1e-6, then * depth). The PM field is
SCALE-INVARIANT, so the mm mesh gives the right B; only the torque needs this geometric
factor.

VALIDATION: a diametric PM rotor (R=10 mm, Br=1 T) in a clean two-teeth salient iron
stator (mu_r=1000), 1 mm gap; the radia eggshell band (10.2->10.8) over theta=0..180 vs a
stored regression reference (an independent weighted-stress-tensor solve):
  * amplitude T0 = 7.6296e-3 (radia) vs 7.5919e-3 N*m (reference) -> **0.50 %**
  * normalised tau(theta) shape -> **1.3 % RMS**, both clean ``T0 sin(2 theta)``
    (reluctance torque: the diametric magnet is pulled to align with the low-reluctance
    iron axis; zeros at 0/90/180, antisymmetric extrema at +-45/135).
Two FULLY INDEPENDENT torque methods + independent meshes agreeing sub-percent also
validates the unit conversion. Locked in:
``validation_test/radia_mcp/test_motor_cogging_torque.py::test_reluctance_torque_amplitude`` (the reference
number is stored -> CI-portable).

EXTEND to true slotted cogging (multi-tooth stator x PM poles): same method, period
360/LCM(slots,poles), >=4 air-gap mesh layers, band fully inside the gap. For a soft/
saturating rotor you must RE-SOLVE per angle (superposition breaks); a hard PM (mu_r=1)
in a LINEAR stator could superpose, but the salient-iron case here is re-solved per angle
anyway (the iron responds). GOTCHA avoided: do not let the salient iron pinch to a zero-
thickness tangency (rect-minus-bore with height = bore diameter) -- it meshes badly and
will not cross-validate; use clearly separated teeth.
"""


NGSOLVE_MEMS_ELECTROMECH = r"""
# MEMS electro-mechanical chain (electrostatic force -> deflection) -- the electric twin of magneto-mechanical

The electrostatic actuator: a parallel-plate pull deflects a movable electrode. ONE-WAY
chain from existing solvers (no new physics):

    solve_electrostatic  ->  electrostatic_eggshell_force_2d  ->  solve_linear_elasticity

```python
from radia_mcp.radia_ngsolve.scalar_fem2d import solve_electrostatic
from radia_mcp.radia_ngsolve.force import electrostatic_eggshell_force_2d, EPS0
from radia_mcp.radia_ngsolve.elasticity import solve_linear_elasticity, cantilever_tip_deflection
V  = solve_electrostatic(mesh, EPS0, {"top": V0, "bottom": 0.0}, order=2)  # gap, sides Neumann
E  = CF((-grad(V)[0], -grad(V)[1]))
Fx, Fy = electrostatic_eggshell_force_2d(E, mesh, CF((0, 1.0/d)), air_region="air")  # gradg=(0,1/d)
P  = abs(Fy)/W                                   # electrostatic pressure [Pa]
u  = solve_linear_elasticity(beam, E_mod, nu, dirichlet="clamp",
                             traction={"loaded": (0, -P)}, plane="stress")  # -> deflection
```

Validated EXACT (validation/force/validate_force_xval.py::validate_mems_electro_mechanical): parallel plate (sides = natural-
Neumann symmetry -> 1-D field, no fringing) gives P = 1/2 eps0 (V0/d)^2 to **0.00 %**, and
the clamped-cantilever tip under that uniform pressure matches Euler-Bernoulli
``w L^4/(8EI)`` (w = P, unit depth) to **0.05 %**.

DEAD END (recorded so it is not repeated): computing the electrostatic force as a BOUNDARY
trace ``oint 1/2 eps0 |E|^2 n dS`` over the electrode FAILS -- on a conductor V = const, so the
boundary trace of ``grad(V)`` in a ``ds`` integral keeps only the TANGENTIAL derivative (=0)
and DROPS the normal field; ``avg|E|^2`` on the face reads ~1e-17 while the interior is ~1e11,
so the surface stress integrates to ~0. FIX: the VOLUME weighted-stress ("eggshell") band --
the interior ``grad(V)`` is the true field. (Same reason the magnetic eggshell beats a B.n
surface trace.)

QUADRATURE GOTCHA: a discontinuous band weight (IfPos: g=0/1/(y1-y0) over y0<y<y1) whose edges
STRADDLE coarse elements integrates inexactly -- 4.9 % error at maxh=d/4, 1.3 % at d/12. FIX:
span the ramp over the WHOLE source-free gap, ``g = y/d`` so ``gradg = (0, 1/d)`` is CONSTANT
and the ramp edges sit on the (mesh-aligned) conductor surfaces -> exact at any mesh. Two-way
pull-in (deformation changes the gap changes the force) is NONLINEAR -> Picard/moving-mesh loop.
"""


NGSOLVE_BACKEMF_FLUX = r"""
# No-load FLUX LINKAGE lambda(theta) / back-EMF -- validated against a regression reference

The other half of PMSM characterisation (the torque pair is ngsolve_usage("cogging")).
Rotate the PM (phi_deg, no remesh), read a stator/search coil's flux linkage at each
angle; for a 2-pole machine lambda(theta)=lambda0 cos(theta) and the no-load back-EMF is
e = -d(lambda)/dt = -omega d(lambda)/d(theta) = omega lambda0 sin(theta).

```python
from radia_mcp.radia_ngsolve.solve import (reluctivity, solve_planar_magnetostatic,
                                           coil_flux_linkage_2d)
nu = reluctivity(mesh, {"iron": 1000.0})
for th in thetas:
    A = solve_planar_magnetostatic(mesh, nu, magnets={"rotor": (Br/MU0, th)}, order=3)
    lam = coil_flux_linkage_2d(A, mesh, "coilpos", "coilneg", n_turns, depth=1.0) * 1e-6
```

``coil_flux_linkage_2d`` = ``N*depth*(<A_z>_pos - <A_z>_neg)``, the AREA-AVERAGE of A_z over
the two coil sides -- the standard circuit "flux linkage" of an N-turn coil (the magnetic
vector potential integrated over the go/return coil sides).

UNITS (same mm->SI factor as the torque): meshed in mm, A_z(mm-solve) = 1000*A_phys (B is
scale-invariant so the mm B-field is right, but A_z is a factor 1e3 too big); so
``lambda[Wb] = coil_flux_linkage_2d(depth=1) * 1e-6`` (depth 1 mm = 1e-3 m; A_z /1e3).

VALIDATION: 2-pole PM rotor + iron stator + air-gap search coil (N=100), theta=0..360 vs a
stored regression reference: lambda0 **0.10 %**, normalised lambda(theta) shape **0.00 % RMS**
(perfect cos(theta) at every angle). Locked in:
``tests/test_motor_backemf.py::test_flux_linkage_amplitude`` (reference number stored -> CI-portable).
For a real winding sum the per-phase coil sides (distributed/short-pitch -> winding factor);
the back-EMF constant ke = omega lambda0.

Kt = Ke (TORQUE CONSTANT == BACK-EMF CONSTANT): by energy conservation a current I in the coil
feels the co-energy torque T = I dlambda_m/dtheta, so Kt = (T(I)-T(0))/I = dlambda_m/dtheta = Ke
[Nm/A == Wb/rad == V s]. ``pm_machine_constant(lam_plus, lam_minus, dtheta)`` gives Ke from the
no-load flux-linkage FD; the FE torque-per-amp (eggshell torque WITH coil current, minus the
PM-only cogging, /I) matches it to **1.5 %**. This MIXES a PM with a CURRENT source -> mesh in
SI METRES (current-driven B is not scale-invariant; depth cancels in Kt=Ke). Validated
tests/test_kt_ke.py.
"""


NGSOLVE_LDLQ = r"""
# Ld / Lq synchronous inductance + frozen-permeability dq -- validated vs a regression reference

The third PMSM quantity (with torque = ngsolve_usage("cogging") and back-EMF =
ngsolve_usage("back_emf")). Drive a stator coil with current, read its inductance
``L = lambda/I`` (``solve.inductance_2d``); the rotor SALIENCY makes Ld (rotor d-axis / iron
along the coil flux) != Lq (rotor q-axis). Ld/Lq is the basis of reluctance torque + MTPA.

```python
from radia_mcp.radia_ngsolve.solve import reluctivity, solve_planar_magnetostatic, inductance_2d
nu = reluctivity(mesh, {"rotor": 1000.0, "stator": 1000.0})
Jz = CF([{"coilpos": j0, "coilneg": -j0}.get(m, 0.0) for m in mesh.GetMaterials()])  # j0=N*I/area
A  = solve_planar_magnetostatic(mesh, nu, Jz=Jz, order=3)
L  = inductance_2d(A, mesh, "coilpos", "coilneg", n_turns=N, current=I, depth=d)   # = lambda/I
# Ld from a mesh with the rotor d-axis along the coil-flux axis; Lq from the rotor rotated 90 deg.
```

MESH IN METRES: a current-driven B is NOT scale-invariant (unlike the PM cases), so solve in
SI metres and L is in henries directly (no mm fudge factor).

FROZEN-PERMEABILITY (standard FE dq method, for a SATURATED machine):
in saturation lambda is NOT a linear PM+armature sum (superposition broken), so you cannot
difference. Instead: (1) full NONLINEAR solve at the operating point (magnet + dq currents) ->
per-element converged mu; (2) FREEZE that mu to a fixed CF; (3) linear dq-perturbation solves
give lambda_pm, Ld=dlambda_d/di_d, Lq=dlambda_q/di_q, cross-coupling Ldq (lambda Park->dq).
These feed the control layer: T=(3/2)p[lambda_pm i_q + (Ld-Lq) i_d i_q], MTPA / field-weakening.

VALIDATION: salient iron rotor + stator coil, radia inductance lambda/I vs a stored regression
reference: **Ld 0.20 %, Lq 0.05 %, saliency Ld/Lq 0.26 %**. Locked in:
tests/test_motor_ldlq.py::test_ldlq_saliency. LESSON: inductance is sensitive to the IN-COIL
field (self/internal inductance) -- a coarse coil mesh undershoots L by ~2 %; refine the coil
region (~7 elements across the wire) to converge. (For a linear machine like this no freeze is
needed; the freeze matters once the iron saturates.)
"""


NGSOLVE_MTPA = r"""
# MTPA (Maximum-Torque-Per-Ampere) -- the dq control capstone of a salient PM machine

Given the FE-extracted machine parameters -- PM flux linkage ``lambda_m`` (from a PM-only
``coil_flux_linkage_2d``) and the d-/q-axis inductances ``Ld``/``Lq`` (from ``inductance_2d``
along the rotor d-/q-axis; use the FROZEN-PERMEABILITY values once the iron saturates,
ngsolve_usage("ld_lq")) -- MTPA picks the stator-current ANGLE that maximises torque per amp.

```python
from radia_mcp.radia_ngsolve.solve import dq_torque, mtpa_operating_point
# dq torque (motor convention), current angle g from the q-axis (id=-I sin g, iq=I cos g):
#   T = (3/2) p [ lambda_m iq + (Ld - Lq) id iq ]   (magnet/alignment + reluctance torque)
g, id_, iq, T = mtpa_operating_point(lambda_m, Ld, Lq, current=I, pole_pairs=p)
```

CLOSED FORM: dT/dg=0 -> lambda_m sin g + (Ld-Lq) I cos 2g = 0, a quadratic in s=sin g:
``2(Ld-Lq) I s^2 - lambda_m s - (Ld-Lq) I = 0``; the T-maximising real root is taken, so it is
ROBUST TO BOTH SALIENCY SIGNS:
  * IPM (Ld < Lq): g > 0 (id < 0, demagnetising) exploits the reluctance torque.
  * salient-iron / inverse saliency (Ld > Lq): g < 0 (id > 0) -- e.g. the radia ld_lq rotor.
  * synchronous reluctance (lambda_m = 0): g = -/+45 deg (pure reluctance).
  * non-salient (Ld = Lq): g = 0 (pure q-axis; all torque is magnet torque).
Behaviour: at low current the magnet term dominates (g ~ 0); as I rises the reluctance term is
increasingly worth exploiting (|g| grows) and MTPA beats pure-q-axis by a widening margin.

VALIDATION:
  * validation_test/radia_mcp/test_motor_mtpa.py -- closed form == independent NUMERICAL argmax over the current
    angle (0.000 %), the stationarity residual ~ 0, current magnitude honoured, MTPA >= pure-q,
    and the limits above. PURE dq theory -> tool-independent (no FE solve, no reference values).
  * An MTPA-locus demo FE-extracts (lambda_m, Ld, Lq) with radia then sweeps
    I; for the small salient rotor the MTPA torque exceeds pure-q-axis by up to ~37 % at high
    current (Ld/Lq=1.23). PITFALL (found building it): the (0,+/-rc) coil links X-directed flux,
    so its d-axis is x -- magnetise the PM along the d-axis (+x), not +y, or lambda_m reads ~0
    (a +y magnet gives A_z ~ x/r^2 = 0 on the x=0 coil axis).
"""


NGSOLVE_FIELD_WEAKENING = r"""
# Field weakening -- the dq operating-region map above base speed (MTPA -> FW -> MTPV)

MTPA (ngsolve_usage("mtpa")) is the LOW-speed answer. Above base speed the inverter VOLTAGE
limit caps the reachable dq flux, so torque-max must trade flux for speed. Two constraints
bound the (id,iq) plane (stator R neglected, electrical speed omega_e):

    current:  id^2 + iq^2 <= Imax^2                                   (a circle)
    voltage:  (Ld id + lambda_m)^2 + (Lq iq)^2 <= (Vmax/omega_e)^2    (an ellipse, centre
              id = -lambda_m/Ld = -Ich, shrinking as omega_e rises)

```python
from radia_mcp.radia_ngsolve.solve import (field_weakening_operating_point,
                                           field_weakening_speed_capability,
                                           base_speed_electrical)
wb = base_speed_electrical(lambda_m, Ld, Lq, Imax, Vmax, p)        # MTPA@Imax hits the V limit
cap = field_weakening_speed_capability(lambda_m, Ld, Imax)         # Ich / CPSR / MTPV verdict
op = field_weakening_operating_point(lambda_m, Ld, Lq, Imax, Vmax, omega_e, p)
# op = (id, iq, torque, region) or None (None = omega_e beyond the machine's max speed)
```

THREE REGIONS (the optimum is always on the feasible boundary, so it is the best of three
closed-form candidates):
  * "MTPA" -- omega_e <= base speed: the MTPA-at-Imax point is still voltage-feasible
    (constant torque). base speed = Vmax / |(Ld id_mtpa + lambda_m, Lq iq_mtpa)|.
  * "FW"   -- on the current circle INTERSECT voltage ellipse (constant ~power; torque falls
    as id is driven more negative). Solve iq^2=Imax^2-id^2 into the ellipse -> a quadratic in id.
  * "MTPV" -- max-torque-per-voltage tangent point INSIDE the current circle. Appears IFF the
    characteristic current Ich = lambda_m/Ld < Imax (wide CPSR). Found by the change of variables
    X = Ld id + lambda_m, Y = Lq iq, which turns the voltage ellipse into a circle X^2+Y^2=Vlim^2
    and the torque into Y*(lambda_m Lq + (Ld-Lq) X); maximising on the circle -> a quadratic in
    cos(theta). If Ich >= Imax there is NO MTPV region and a FINITE max speed (feasible set goes
    empty -> the helper returns None) -- the rule of thumb Ich = Imax is the "infinite max speed"
    design target.

``field_weakening_speed_capability`` exposes the geometry in one readable dict:
``characteristic_current`` (Ich), ``infinite_speed_possible`` (Ich <= Imax),
``finite_max_speed`` (Ich > Imax), ``mtpv_possible`` (Ich < Imax), and the current margin/ratio.

VALIDATION (validation_test/radia_mcp/test_field_weakening.py): the
closed-form operating point == an INDEPENDENT brute-force numeric constrained argmax of T over
the feasible (id,iq) to **0.000 %**, across a wide-CPSR IPM (Ich<Imax, exercises MTPV), a
current-limited IPM (Ich>Imax, finite max speed -- helper-None == numeric-None), and a
non-salient PMSM. Torque is monotone-decreasing through FW; MTPV presence matches the Ich<Imax
criterion. PURE dq theory -> tool-independent. The high-speed sequel to MTPA (#26); feeds the
T-N envelope / efficiency-map layer.
"""


NGSOLVE_INDUCTION_MACHINE = r"""
# Induction machine -- torque-slip curve + breakdown (single-cage equivalent circuit)

A NEW motor class beside the PM dq blocks (MTPA, field weakening): the 3-phase induction
machine. The single-cage equivalent circuit, reduced by Thevenin (shunt Xm -> source):

```python
from radia_mcp.radia_ngsolve.solve import (induction_machine_thevenin,
    induction_machine_torque, induction_machine_breakdown,
    induction_machine_noload_lockedrotor_summary)
tests = induction_machine_noload_lockedrotor_summary(
    V0_LL, I0_line, P0_total, Vlr_LL, Ilr_line, Plr_total, f_line, pole_pairs)
Vth, Rth, Xth = induction_machine_thevenin(V1, R1, X1, Xm)
T  = induction_machine_torque(V1, R1, X1, R2, X2, Xm, omega_s, slip)   # air-gap torque [N.m]
s_max, T_max = induction_machine_breakdown(V1, R1, X1, R2, X2, Xm, omega_s)
```

    T(s) = (3/omega_s) Vth^2 (R2/s) / [(Rth+R2/s)^2 + (Xth+X2)^2]   (omega_s = SYNC mech. speed)

with the BREAKDOWN (pull-out) maximum at
    s_max = R2 / sqrt(Rth^2 + (Xth+X2)^2),
    T_max = 3 Vth^2 / (2 omega_s [Rth + sqrt(Rth^2 + (Xth+X2)^2)]).

KEY INVARIANT: T_max is INDEPENDENT of the rotor resistance R2 -- only the breakdown SLIP
scales with R2 (s_max ∝ R2). That is exactly how wound-rotor / deep-bar designs push peak
torque toward standstill for high starting torque WITHOUT changing its height. Curve: T -> 0
at synchronism (s->0), rises to T_max at s_max, the s=1 value is the starting torque.

VALIDATION (tests/test_induction_machine.py):
the closed-form (s_max, T_max) == an INDEPENDENT numeric slip sweep argmax to <1e-3 / <1e-4,
and the R2-invariance of T_max (with s_max ∝ R2) holds exactly. PURE circuit theory ->
tool-independent. Continuous-loop slot 36 adds the no-load / locked-rotor front gate:
400 V/5 A/900 W no-load gives Rc=177.7777777777778 ohm and Xm=47.830502044476084 ohm;
90 V/20 A/1200 W locked-rotor gives Req=1.0 ohm, Xeq=2.3979157616563596 ohm, and
R2'=0.6 ohm if R1=0.4 ohm. The analytic basis for the JMAG-roadmap IM 1/4-model FE study
(anti-periodic + slip + cage), the induction-machine companion to dq_torque for PM machines.
"""


NGSOLVE_POWER_ANGLE = r"""
# Synchronous machine power-angle (delta) curve + pull-out -- the grid/voltage frame

The voltage-frame characteristic of a synchronous machine, COMPLEMENT to the current-frame
MTPA (#26) / field weakening (#37): torque vs the LOAD ANGLE delta (between terminal voltage V
and excitation EMF E).

```python
from radia_mcp.radia_ngsolve.solve import synchronous_power_angle_torque, synchronous_pullout
T = synchronous_power_angle_torque(V, E, Xd, Xq, omega_e, p, delta)   # [N.m]
d_max, T_max = synchronous_pullout(V, E, Xd, Xq, omega_e, p)          # stability limit
```

    T(d) = (3 p/omega_e) [ V E/Xd sin d  +  (V^2/2)(1/Xq - 1/Xd) sin 2d ]

= FIELD/excitation torque (~ sin d, from the EMF E = omega_e*lambda_m for a PM rotor) + the
RELUCTANCE torque (~ sin 2d, from the saliency Xd != Xq). The steady-state stability limit is
the PULL-OUT (the peak); beyond it the rotor slips poles.

PULL-OUT closed form (dT/dd=0 -> (V E/Xd) cos d + V^2(1/Xq-1/Xd) cos 2d = 0, a quadratic in
cos d). LIMITS that anchor it: NON-salient (Xd=Xq) -> pull-out at exactly 90 deg (the classic
round-rotor result, all field torque); pure RELUCTANCE (E=0, synchronous reluctance machine) ->
45 deg; a salient field machine sits BETWEEN (the reluctance term advances the peak below 90 deg,
e.g. 67 deg for Xd/Xq=1.67). Stronger saliency advances it further.

VALIDATION (tests/test_synchronous_power_angle.py):
closed-form (d_max, T_max) == an INDEPENDENT numeric angle sweep argmax to <1e-3/<1e-4; the
non-salient -> 90 deg and reluctance-only -> 45 deg limits hold exactly; field+reluctance
decomposition is exact. PURE circuit theory -> tool-independent. Same torque as dq_torque (#26)
re-expressed in the (V, E, delta) grid frame; the synchronous-machine companion to the
induction-machine torque-slip (#40) -- together they span the AC-machine torque characteristics.
"""


NGSOLVE_DQ_OPERATING_POINT = r"""
# dq steady-state operating point -- terminal voltage, power, power factor, efficiency

Where MTPA (#26) / field weakening (#37) CHOOSE the dq currents, this EVALUATES the terminal
quantities the inverter sees at (id, iq, omega_e) -- the voltage/power layer of the dq model:

```python
from radia_mcp.radia_ngsolve.solve import dq_voltages, dq_operating_point
vd, vq = dq_voltages(R, Ld, Lq, lambda_m, id, iq, omega_e)
op = dq_operating_point(R, Ld, Lq, lambda_m, id, iq, omega_e, p)   # dict: Vmag, P_in, PF, ...
```

PM flux-linkage constants (motor report sanity checks): with peak phase dq
variables, no-load ``e_phase_peak = omega_e lambda_m = omega_mech p lambda_m`` and
surface-PM ``T/iq_peak = (3/2) p lambda_m``.  Use
``pm_flux_linkage_constants(lambda_m, p)`` for phase-peak, phase-RMS,
line-line-RMS Ke and peak/RMS Kt, and ``pm_no_load_back_emf(lambda_m, omega_mech, p)``
for speed rows.  `validation_pm_emf_constant_table.py` checks dq voltage and power
identities: baseline `lambda_m=0.1`, `p=4` gives `Ke_ll_rms=0.4898979485566356`,
`Kt_peak=0.6000000000000001`, and `Kt_rms/Ke_ll_rms=sqrt(3)`; no-load vq equals
phase-peak EMF exactly and `T*omega_mech=(3/2)e_phase_peak iq_peak` to roundoff.

    lambda_d = Ld id + lambda_m,  lambda_q = Lq iq
    v_d = R id - omega_e lambda_q,   v_q = R iq + omega_e lambda_d
    P_in = (3/2)(v_d i_d + v_q i_q),  P_em = (3/2) omega_e (lambda_d i_q - lambda_q i_d) = T omega_mech,
    P_cu = (3/2) R |I|^2,  power_factor = P_in/((3/2)|V||I|),  efficiency = P_em/P_in.

Identities (the self-consistency, all to machine precision): **P_in = P_cu + P_em** (energy
balance), **P_em = torque*omega_mech** (air-gap power), and the torque equals dq_torque (#26).
At R=0 (lossless) P_in = P_em and efficiency = 1. |V| rises ~linearly with speed (back-EMF) --
when it hits the inverter ceiling you are at the field-weakening limit (#37); |I| vs the device
rating. The terminal (V, I, P, PF) layer for inverter/drive SIZING; closes the dq loop
currents (MTPA/FW) -> terminals. Validated tests/test_dq_operating_point.py.

For map exports, run ``pm_drive_efficiency_map_health(sweep)`` before comparing rows from a FE
tool.  It checks the table contract: efficiencies bounded in [0,1], ``P_in = P_em + P_cu``,
``omega_e = p omega_mech``, nondecreasing speed rows, and required regions such as MTPA/FW/MTPV.
Continuous-loop slot 37 fixes the motor speed-map gate: regions [MTPA, FW, FW, MTPV],
max efficiency 0.994425732150756 in MTPV, max output power 2597.750669164122 W in FW, and
power/speed contract errors below roundoff.

Continuous-loop slot 141 extends the notebook handoff rule: a useful FE operating-point row is
not only an efficiency number.  Preserve an ``operating_point_id`` plus region, speed, torque,
``id/iq``, terminal-voltage and power-factor columns, loss buckets, and DC-bus voltage-margin
metadata.  Run ``pm_drive_efficiency_map_health`` together with
``pm_drive_terminal_table_health`` and ``pm_drive_loss_bucket_efficiency_gate`` so the notebook
can show why the row is feasible, where the losses went, and which voltage-utilization value is
report-only versus selector feasibility.
"""


NGSOLVE_SHORT_CIRCUIT = r"""
# PM machine 3-phase short circuit -- fault current, braking torque, demagnetisation

The FAULT / protection regime (vs the operating point above). With the terminals SHORTED
(v_d = v_q = 0) at electrical speed omega_e, the dq voltage equations give

```python
from radia_mcp.radia_ngsolve.solve import (short_circuit_dq_currents, characteristic_current,
                                           dq_torque)
id, iq = short_circuit_dq_currents(R, Ld, Lq, lambda_m, omega_e)
T_brake = dq_torque(lambda_m, Ld, Lq, id, iq, p)
Ich = characteristic_current(lambda_m, Ld)          # = lambda_m/Ld
```

    id = -we^2 Lq lm /(R^2 + we^2 Ld Lq),   iq = -we R lm /(R^2 + we^2 Ld Lq).

HIGH-SPEED limit (we -> inf, or R -> 0): id -> -lambda_m/Ld = -Ich (the CHARACTERISTIC current,
PURE d-axis DEMAGNETISING -- it cancels the magnet flux, lambda_d = 0), iq -> 0, so the steady
short-circuit current |Isc| -> Ich. Two design checks fall out: (1) the DEMAG current Ich vs the
magnet knee (a machine with small Ich = lambda_m/Ld survives a short), and (2) the BRAKING TORQUE
= dq_torque(id, iq), which rises from 0, PEAKS near the critical speed (exactly we = R/Ls for a
non-salient machine, Ld=Lq=Ls), then decays to 0 at high speed (iq -> 0). Note Ich = Imax is also
the wide-CPSR / 'infinite max speed' field-weakening target (#37).

VALIDATION (tests/test_short_circuit.py): the
shorted dq equations hold exactly, |Isc| -> Ich to <1e-3, the braking torque peaks then decays,
and the non-salient critical speed we=R/Ls is reproduced. PURE dq theory -> tool-independent.
"""


NGSOLVE_WINDING_FACTOR = r"""
# AC winding factor k_w = k_d * k_p -- what scales the back-EMF and shapes its harmonics

The winding-geometry factor behind the back-EMF / torque constant (Kt=Ke, #34, F3): the n-th
harmonic EMF is k_w,n times the concentrated-full-pitch value.

```python
from radia_mcp.radia_ngsolve.solve import (winding_distribution_factor, winding_pitch_factor,
                                           winding_factor)
kw1 = winding_factor(1, q=slots_per_pole_per_phase, slot_angle_elec_deg=gamma, pitch_fraction=beta)
```

    DISTRIBUTION  k_d,n = sin(n q gamma/2)/(q sin(n gamma/2)) = |sum_{i<q} exp(j n i gamma)|/q
    PITCH         k_p,n = sin(n beta pi/2)              (beta = coil span / pole pitch)
    k_w,n = k_d,n * k_p,n

q = slots per pole per phase, gamma = electrical slot angle (= p*360/Q deg). LIMITS: q=1
(concentrated) -> k_d=1; beta=1 (full pitch) -> k_p,1=1. The fundamental k_w,1 ~ 0.92-0.96 for a
good 3-phase winding (it scales the useful EMF/torque a bit below 1). The POINT is harmonic
suppression: distributing (q>1) and SHORT-PITCHING (beta<1) kill the low harmonics -- the classic
5/6 pitch makes |k_p,5| = |k_p,7| = sin(pi/12) = 0.259, so the 5th & 7th EMF harmonics (and their
parasitic torques) nearly vanish while k_w,1 stays ~0.93.

VALIDATION (tests/test_winding_factor.py): k_d == the
direct phasor sum of the q slot EMFs (exact), the q=1 / beta=1 limits, and the 5/6-pitch 5th/7th
suppression. PURE winding geometry -> tool-independent; feeds the back-EMF amplitude (F3) and its
harmonic content.
"""


NGSOLVE_FROZEN_PERM = r"""
# Frozen-permeability superposition (saturated-machine Ld/Lq engine)

In a SATURATED magnetic circuit the flux linkage is NOT a linear sum of its sources --
nu(|B|) couples them -- so you cannot separate PM / d-axis / q-axis flux by differencing.
The frozen-permeability method restores superposition: take a CONVERGED nonlinear solve,
FREEZE nu at its operating-point field, and the linear frozen-nu problem superposes exactly.

```python
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic_nonlinear,
    solve_planar_magnetostatic, coil_flux_linkage_2d, frozen_reluctivity)
A0 = solve_planar_magnetostatic_nonlinear(mesh, nu_of_B, Jz=JzA+JzB, order=2)  # operating point
nu_fz = frozen_reluctivity(A0, nu_of_B)                       # snapshot nu(|B0|), DIRECT CF
A_a = solve_planar_magnetostatic(mesh, nu_fz, Jz=JzA)         # linear, source A only
A_b = solve_planar_magnetostatic(mesh, nu_fz, Jz=JzB)         # linear, source B only
# lambda_A(nonlinear) == lambda_A(frozen, A) + lambda_A(frozen, B)
```

Self-validating, NO external data (tests/test_frozen_permeability.py): a 2-winding saturable window-frame core at a strongly
nonlinear operating point (iron mu_r 1000 -> ~45) -- the nonlinear coil flux linkage
decomposes into superposable frozen-nu parts to **0.51 %** (and the frozen-nu solve
reproduces the nonlinear operating point to the same 0.51 %, which IS the convergence check).
Extends to dq: freeze at (id,iq,PM), perturb each axis -> Ld(id,iq), Lq, Ldq, lambda_pm maps
-> the MTPA / T-N control layer (ngsolve_usage("ld_lq")).

DEAD ENDS / lessons:
  * L2-projecting a clamped/discontinuous nu SMEARS it -> the freeze stops reproducing the
    operating point (saw ~20 % error). FIX: use the DIRECT CF nu_of_B(curl A0) (frozen_reluctivity).
  * A HARD clamp (IfPos at Bsat) makes Picard non-convergent (the fixed point oscillates;
    freeze-reproduces blew to >100 %). FIX: a SMOOTH BOUNDED nu(B) (mu_r0 -> 1 via
    B^n/(B^n+Bk^n)) -- Picard converges and the freeze is exact. The freeze-reproduces-nonlinear
    check is the built-in correctness gate: if it is off, the nonlinear solve has not converged.

SATURATING INDUCTANCE CURVE L(I): the SECANT inductance L_sec(I)=lambda_nl(I)/I (one nonlinear
solve per current -- ``solve.saturated_secant_inductance``) falls below the unsaturated linear
L0 as the iron saturates -- the Ld(i) knee underlying saturated dq / MTPA maps. Self-consistent:
low-I L_sec -> L0 (nu=nu0); L_sec monotone-decreasing; deep saturation pulls it well below L0
(mu_r_op << mu_r0). Validated tests/test_saturating_inductance.py. NB the frozen-perm APPARENT inductance (freeze the SECANT
nu(B)) reproduces L_sec at the operating point; the INCREMENTAL dlambda/di needs the TANGENT
reluctivity dH/dB (not the secant nu) -- freezing the secant nu gives the apparent/secant L,
NOT the differential -- a noted refinement for the incremental dq map.

INCREMENTAL inductance + CO-ENERGY (the differential refinement): on the saturating lambda(I),
L_inc = dlambda/dI (``incremental_inductance``, central FD across nonlinear solves at I0+/-dI)
is the small-signal/control inductance; the co-energy W'(I)=int_0^I lambda dI'
(``magnetic_coenergy``) and energy W=lambda*I-W'. Robust self-consistent physics of a CONCAVE
lambda(I): **L_inc < L_sec everywhere**, both fall monotonically (saturation), **co-energy W' >
energy W** once saturating, and the exact identity **W+W'=lambda*I**. Validated
tests/test_incremental_inductance.py. (FD is
the ground-truth L_inc; the frozen-perm L_inc would freeze the tangent differential-reluctivity
TENSOR dH/dB -- isotropic ALONG B, secant PERPENDICULAR -- still a future tensor refinement.)

INCREMENTAL inductance MATRIX + cross-saturation + RECIPROCITY (#52, the multi-port version):
``incremental_inductance_matrix(flux_linkages, currents, di)`` builds the n x n
L_jk = dlambda_j/di_k by central FD of the nonlinear flux-linkage map (one nonlinear solve per
+-perturbation). Two exact properties on a SHARED SATURABLE core (two windings / dq axes):
  * **RECIPROCITY L_jk = L_kj** -- the matrix is SYMMETRIC (the magnetic co-energy Hessian
    d^2 W'/di_j di_k); the off-diagonals are the cross-(saturation) inductances (L_dq, ZERO
    without saturation, growing as shared iron saturates). FE gives |L_AB-L_BA|/L_AB < 0.01 %.
  * **INCREMENTAL << APPARENT(secant)**: the frozen-perm secant inductance (freeze nu at the
    operating |B|) is ~4x the tangent incremental in deep saturation (mu_r 1000->45), and the
    incremental cross inductance COLLAPSES (~0.04x the unsaturated linear mutual). Use the
    tangent matrix for small-signal / dq-current-loop / control design, the secant for the
    operating flux. Resolves the #28/#31 secant-vs-incremental note for the MATRIX. Validated
    tests/test_cross_saturation_tensor.py.
"""


NGSOLVE_ELECTROTHERMAL_2WAY = r"""
# Two-way temperature-dependent sigma(T) electro-thermal

The genuinely 2-way electro-thermal coupling: electrical conductivity falls with temperature,
sigma(T)=sigma0/(1+alpha T), so the Joule density q = J^2/sigma(T) = (J^2/sigma0)(1+alpha T)
RISES with T and feeds back into the heat (thermal <-> electrical).  Solved by a Picard fixed
point with the temperature-dependent source:

```python
from radia_mcp.radia_ngsolve.multiphysics import solve_heat_steady_nonlinear
T = solve_heat_steady_nonlinear(mesh, lambda Tcf: (J*J/sigma0)*(1.0 + alpha*Tcf),
                                k, conductor="bar", dirichlet="ends", order=2)
```

``q_of_T(T_cf) -> q_cf`` builds the source from the current temperature; the loop re-solves
``solve_heat_steady`` until ||T|| converges (a few iters for mild feedback).

EXACT closed form (1-D bar, uniform J, cooled ends): the ODE -k T'' = (J^2/sigma0)(1+alpha T)
is LINEAR -> T(x)=(1/alpha)[cos(s x)+((1-cos sL)/sin sL) sin(s x)-1], T_max=(1/alpha)(sec(sL/2)-1),
s=sqrt(alpha J^2/(sigma0 k)); as alpha->0 this recovers the constant-sigma parabola peak
J^2 L^2/(8 sigma0 k).  Validated **0.00 %** (validation/force/validate_force_xval.py::validate_electrothermal_sigmaT_2way) -- the feedback raises the peak ~7 %
over the no-feedback parabola (alpha*T_max ~ 0.09), so the 2-way effect is real and captured.

``solve_heat_steady_nonlinear`` is GENERAL (any T-dependent volumetric source: sigma(T) Joule,
volumetric reaction/radiative terms).  The 2-way twin of the one-way ``solve_heat_steady``.  For
full 2-D 2-way with current redistribution, also put sigma(T) in ``solve_current_flow`` inside the
loop; in 1-D the current density is uniform (conservation) so q=J^2/sigma(T) is exact.
"""


NGSOLVE_COAX_LINE = r"""
# Coaxial transmission line: external inductance + characteristic impedance

The magnetic companion to the coax capacitance.  Inner conductor carries +I, shield -I; in
the dielectric a<r<b the field is exactly the Ampere field B = mu0 I/(2 pi r), so the external
inductance per unit length is L_ext = 2W/I^2 = (mu0/2pi) ln(b/a),  W = int_diel |B|^2/(2mu0) dA.

```python
from radia_mcp.radia_ngsolve.force import magnetic_energy_2d
A = solve_planar_magnetostatic(mesh, reluctivity(mesh, {}), Jz=Jz_coax, order=3)  # inner +I, shield -I
B = CF((grad(A)[1], -grad(A)[0]))
L = 2*magnetic_energy_2d(B, mesh, "diel")/I**2          # = (mu0/2pi) ln(b/a)
```

With the coax capacitance C = 2 pi eps0/ln(b/a): Z0 = sqrt(L/C) = (1/2pi) sqrt(mu0/eps0) ln(b/a)
= **60 ln(b/a) ohm**, and the wave speed v = 1/sqrt(LC) = **c**.  Validated
(validation/force/validate_force_xval.py::validate_coax_line_inductance):
L_ext **0.46 %**, Z0 0.30 %, v=c to 0.2 %.

LESSON: the |B|^2 energy of a 1/r field peaks at small r -> refine the mesh near the inner
radius (energy converged 1.8 % -> 0.46 % under refinement).  ``magnetic_energy_2d`` is the
general energy->inductance helper (L = 2W/I^2), the energy twin of the lambda/I ``inductance_2d``.

TWO-WIRE LINE (the parallel-wire dual): two wires radius a, centre separation D carrying +-I
form a 2D LINE-DIPOLE whose field energy CONVERGES (B ~ 1/r^2 far field), so the external
inductance L_ext = 2W/I^2 = (mu0/pi) acosh(D/2a) = ``two_wire_external_inductance`` -- the
energy is taken over the WHOLE air with NO outer-boundary dependence (unlike a single wire's
log-divergent energy). With C = pi eps0/acosh(D/2a): v = 1/sqrt(LC) = c, Z0 = 120 acosh(D/2a) ohm.
Validated tests/test_two_wire_inductance.py:
L_ext **1.8 %**, v=c 0.9 % -- the magnetic twin of the two-wire capacitance.

WIRE OVER GROUND (#63) -- the single-ended / microstrip-style line: a wire (radius a) at height h
above a PERFECT GROUND. The IMAGE METHOD: the ground (A=0 plane) reflects a -I image at -h, so the
field above is the +I/-I pair at separation 2h but only HALF the energy is above the plane ->
``wire_over_ground_inductance(h, a) = (mu0/2 pi) acosh(h/a) = 0.5*two_wire_external_inductance(2h, a)``.
radia reads it from ``magnetic_energy_2d`` over the air ABOVE the Dirichlet-A=0 ground (L=2W/I^2);
Z0 = c L_ext is the single-conductor line. FE vs closed form to a few % (open-domain truncation,
FE BELOW the h->inf-domain value, same as the two-wire) -- tests/test_wire_over_ground.py. The MAGNETIC image method (-I image) on the SAME footing as the
current-wire-above-iron force (#33) and the coax/two-wire family.

TOTAL LOOP inductance (#48) -- what a meter reads: the external field inductance PLUS BOTH
wires' internal: ``two_wire_loop_inductance`` = ``(mu0/pi) acosh(D/2a) + 2*(mu0/8pi)`` =
L_ext + mu0/(4 pi). radia PARTITIONS the FE field energy: air -> external, the two wire
interiors -> internal. The internal term mu0/4pi (radius-independent, the DC value -- rolls off
via the skin effect, ``skin_effect_internal_inductance_ratio``) is got to **0.18 %** (exact,
Ampere fixes the interior B), the external acosh to ~1.8 % (energy truncation), the total loop
to **1.6 %**. tests/test_two_wire_loop_inductance.py.

INTERNAL inductance (the energy stored INSIDE the conductor): a round wire with a UNIFORM (DC)
current has B(r)=mu0 I r/(2 pi R^2) for r<R, so the interior energy gives
L_int = 2 W_in/I^2 = **mu0/(8 pi)** [H/m] = ``internal_inductance_round_wire()`` -- INDEPENDENT
of the radius. radia reads it straight off the wire-region field energy, which is
boundary-independent (Ampere fixes the interior B from the enclosed current):

```python
A = solve_planar_magnetostatic(mesh, reluctivity(mesh, {}), Jz=I/(pi*R**2) on "wire", order=3)
B = CF((grad(A)[1], -grad(A)[0]))
L_int = 2*magnetic_energy_2d(B, mesh, "wire")/I**2     # = mu0/(8 pi), any R
```

The total low-frequency line inductance is L_ext + L_int; L_int is what ROLLS OFF as the skin
effect expels the current to the surface at high frequency (the AC-resistance companion).
Validated tests/test_internal_inductance.py:
**0.24 %** at two radii (radius-independence confirmed). Completes the inductance set:
external coax (mu0/2pi)ln(b/a), external two-wire (mu0/pi)acosh(D/2a), internal mu0/(8pi).

AC ROLL-OFF (#45) -- mu0/(8 pi) is the DC limit; at frequency the skin effect expels current
from the core, so the INTERNAL inductance ROLLS OFF while the AC RESISTANCE rises. For a round
wire at q = a sqrt(omega mu0 sigma) = sqrt(2) a/delta, the Kelvin (ber/bei) ratios are

    skin_effect_resistance_ratio(q)          = Rac/Rdc        = (q/2)(ber bei'-bei ber')/(ber'^2+bei'^2)
    skin_effect_internal_inductance_ratio(q) = L_int/L_int_dc = (4/q)(ber ber'+bei bei')/(ber'^2+bei'^2)

so Z_internal(omega) = Rdc*Rac_ratio + j*omega*(mu0/8pi)*Lint_ratio. q->0: both ->1 (DC,
mu0/8pi); q>>1: Rac/Rdc -> q/(2 sqrt2)+1/4, L_int/L_dc -> 2 sqrt2/q (->0). radia gets BOTH from
ONE current-driven `solve_planar_eddy` on the WIRE ALONE with A_z=0 on its surface (so Vc/I is
the pure INTERNAL impedance, no external term to subtract): Rac=Re(Vc/I), L_int=Im(Vc/I)/omega.
FE vs Kelvin **<0.06 %** L_int and **<1 %** Rac over q=0.5..8 (tests/test_ac_internal_inductance.py). The high-frequency complement
to the DC internal inductance (#36); the AC resistance is the existing skin-effect test
(validation_test/radia_mcp/test_planar_eddy.py), also vs the Kelvin closed form.
"""


NGSOLVE_READABLE_P1_FEM = r"""
# Readable P1 tetrahedron FEM blocks -- MATLAB/Gypsilab-style teaching code

The small clean-room blocks in ``radia_mcp.radia_ngsolve.scalar_fem3d`` are meant to be mirrored
line-for-line in MATLAB teaching scripts: no solver magic, only local P1 tetrahedron formulas.

```python
from radia_mcp.radia_ngsolve.scalar_fem3d import (
    p1_tetrahedron_geometry, p1_tetrahedron_stiffness,
    p1_tetrahedron_gradient, p1_tetrahedron_flux,
    p1_tetrahedron_boundary_fluxes)

grad_u = p1_tetrahedron_gradient(vertices, nodal_values)
q = p1_tetrahedron_flux(vertices, nodal_values, coeff=k)              # q = -k grad(u)
rows = p1_tetrahedron_boundary_fluxes(vertices, nodal_values, coeff=k) # outward Neumann rows
```

For an affine field, ``grad_u`` is exact on each tetrahedron and the stiffness energy satisfies
``u^T K u = k volume |grad_u|^2``. The face rows provide the FEM->BEM trace primitive:
``integrated_flux = q . outward_area_vector`` and the closed tetrahedron sum is zero for a constant
flux. `validation_p1_tet_flux_trace.py` checks these identities on a Netgen `.vol` unit tetrahedron:
gradient/flux errors < roundoff, energy error < roundoff, four face rows matching the surface
triangles, and total integrated outward flux ~0. This is the readable bridge from volume P1 FEM to
boundary Neumann data.
"""


NGSOLVE_ACOUSTIC_BEM = r"""
# Acoustic FEM/BEM readable gates -- Helmholtz DtN, radiation impedance, Robin rows

Scalar acoustics is the cleanest FEM/BEM coupling classroom: the FEM side owns a pressure trace,
and the exterior radiation/BEM side returns a normal derivative or impedance.  The public helpers
in ``radia_mcp.radia_ngsolve.acoustics`` keep the sign conventions explicit:

```python
from radia_mcp.radia_ngsolve.acoustics import (
    acoustic_dtn_from_impedance, acoustic_impedance_from_dtn,
    planar_helmholtz_dtn_symbol, spherical_helmholtz_dtn_eigenvalue,
    planar_mode_radiation_impedance, spherical_mode_radiation_impedance)

lam = acoustic_dtn_from_impedance(f, specific_impedance=z, rho=rho)["dtn_eigenvalue"]
z = acoustic_impedance_from_dtn(f, lam, rho=rho)["specific_impedance"]
```

With ``exp(+i omega t)``, Euler's equation gives ``v_n = i (partial_n p)/(omega rho)``.  Therefore
an impedance/admittance boundary is the Robin/DtN coefficient
``partial_n p = lambda p`` with ``lambda = -i omega rho / z = -i omega rho Y``.  The validation
example `validation_acoustic_impedance_dtn_bridge.py` checks this bridge against planar and
spherical radiation modes and a baffled piston average impedance round trip.  Companion examples
cover low-frequency Helmholtz kernel splitting, spherical/planar DtN modes, pulsating spheres, and
baffled-piston radiation.
"""


NGSOLVE_CROSS_VALIDATION_REGISTRY = r"""
# Cross-validation registry -- what was recovered and how MCP uses it

The loop output is useful only when it lands in reusable artifacts.  For
radia-ngsolve, a completed validation item should leave three public-safe
breadcrumbs:

1. a runnable `validation_test/.../validation_*.py` script,
2. a machine-readable `validation_test/.../validation_*_summary.json`, and
3. a knowledge/API hook that explains how the number guards future FEM/BEM,
   force, mesh, CAD, wave, or readable-MATLAB work.

Private lab campaign notes may also exist, but public MCP knowledge should not
depend on those paths or on commercial-tool provenance.  In public text, keep
the reusable math, API shape, tolerances, and failure modes; leave private
solver names, internal paths, and licensed benchmark attribution out of the
repo.

## Reusable validation families already harvested

| Family | Reusable artifacts | MCP topic/API that learned from it |
|---|---|---|
| Acoustic FEM/BEM | `validation_test/acoustic_bem/validation_*_summary.json` | `ngsolve_usage("acoustic")`, `radia_ngsolve.acoustics` |
| Low-frequency Helmholtz kernel | `validation_low_frequency_helmholtz_kernel_summary.json` | stable singular-plus-regular BEM kernel split |
| Spherical/planar DtN | `validation_spherical_dtn_modes_summary.json`, `validation_planar_dtn_symbol_summary.json` | exact radiation Robin/DtN sign gates |
| Baffled piston and boundary power | `validation_baffled_piston_radiation_summary.json`, `validation_acoustic_boundary_power_summary.json` | acoustic power and impedance conventions |
| Netgen `.vol` tri/tet mesh | `validation_test/cubit_mesh_export/validation_vol_*_summary.json` | `.vol` tri/tet parser, boundary inventory, incidence, quality |
| build123d to mesh/CAD checks | `validation_test/build123d_netgen_gmsh_flow/validation_*_summary.json` | CAD face pressure/traction rows, area/volume consistency |
| Electromagnetic force and torque | `validation_test/electric_machine/validation_*_summary.json` | `force_validation("method_map")`, `force_validation("cross_validation")` |
| RF/waveguide momentum | `validation_test/rf_waveguide/validation_*_summary.json` | radiation pressure, S-parameter momentum, TE/TM cutoff gates |
| MCP server fleet quality | `packages/radia-mcp/tests/test_meta_health.py` plus per-server `--selftest` | `radia_mcp_golden_gate()`, catalog and publish-boundary regressions |
| Readable first-order FEM | `test_p1_*`, `test_surface_triangle_*`, `test_tetrahedron_lorentz_force.py` | P1 triangle/tet energy, flux, HCurl/RWG trace teaching primitives |
| FEM/BEM open boundary | `test_fem_bem_coupling.py`, `fem_bem_schur("api")` | dense reference Schur coupling and Kelvin/BEM DtN comparisons |

## What "MCP became smarter" means here

The target is not to grow a private solver wrapper inside radia-mcp.  The target
is that future agents can ask the public MCP for the right gate before coding:

```text
ngsolve_usage("cross_validation_registry")  -> where reusable artifacts live
ngsolve_usage("acoustic")                   -> acoustic FEM/BEM sign gates
force_validation("method_map")              -> force/torque method selection
force_validation("cross_validation")        -> stored neutral regression cases
fem_bem_schur("api")                        -> exact-open-boundary Schur recipe
ngsolve_usage("readable_fem")               -> P1 tri/tet educational primitives
panel_review("cubit_boundary")              -> checked .vol and Simulink boundary
radia_mcp_golden_gate()                     -> public MCP fleet quality gate
```

When a loop slot produces only a one-off note, it is not finished.  Promote the
lesson into one of the rows above, or add a new row with a summary JSON and a
knowledge topic.

## Optuna boundary

The old in-package Optuna MCP path is retired.  There must be no
`radia_mcp.optuna` package, no `mcp-server-optuna` console script, and no Optuna
runtime dependency in radia-mcp. Every shared operation present in the
official external `optuna/optuna-mcp` live `tools/list` belongs to that server.
`mcp-server-radia-matlab` may support only the differences of the standalone
MATLAB distribution: table/MAT persistence, Simulink monitoring and failure
telemetry, MATLAB parallel execution, `optuna_mex`, and Radia CAE artifacts.
radia-mcp may also keep CAE objective helpers, analytic gates, and example
scripts that are optional/plain-Python, but it must not proxy, rename, or
reimplement an upstream Optuna MCP tool.

The policy guard in `packages/radia-mcp/tools/policy_lint.py` enforces this
boundary with synthetic tests for server entry points, dependencies, source
imports, and a returned `radia_mcp.optuna` package.
"""


NGSOLVE_WAVEGUIDE = r"""
# Waveguide / cavity cutoff modes -- the 2D Helmholtz EIGENVALUE problem

The wave-physics entry: a hollow metallic waveguide's modes solve the TRANSVERSE Helmholtz
eigenproblem on the cross-section,

    -nabla_t^2 psi = k_c^2 psi ,

an EIGENVALUE problem (not a source-driven BVP). The wall BC picks the mode family (PEC walls):

    * TM modes (E_z != 0): E_z = 0 on the wall      ->  DIRICHLET Laplacian
    * TE modes (H_z != 0): dH_z/dn = 0 on the wall  ->  NEUMANN  Laplacian

The eigenvalues are the squared cutoff wavenumbers k_c^2; below the cutoff frequency
f_c = c k_c/(2 pi) a mode is evanescent. ``radia_ngsolve.waveguide``:

```python
from radia_mcp.radia_ngsolve.waveguide import (
    rectangular_waveguide_cutoff, cutoff_frequency, helmholtz_cutoff_wavenumbers_2d)

# closed form (rectangle a x b): f_c,mn = (c/2) sqrt((m/a)^2 + (n/b)^2)
f10 = rectangular_waveguide_cutoff(a, b, 1, 0)      # dominant TE10 = c/(2a)

# FE eigensolve on ANY cross-section (rectangular, ridged, circular, L-shaped):
te = helmholtz_cutoff_wavenumbers_2d(mesh, 3, bc="neumann")    # TE modes (k_c, ascending)
tm = helmholtz_cutoff_wavenumbers_2d(mesh, 2, bc="dirichlet")  # TM modes
f = cutoff_frequency(te[0])                                    # c k_c/(2 pi)
```

Rectangular guide spectrum (EXACT): k_c^2 = (m pi/a)^2 + (n pi/b)^2. TE needs m,n>=0 (not 0,0),
TM needs m,n>=1. With a>b the dominant mode is TE10 (f_c=c/2a); the gap up to the next mode
(TE20 at c/a, or TE01 at c/2b) is the SINGLE-MODE bandwidth. WR-90 X-band (22.86x10.16 mm):
TE10 6.557, TE20 13.11, TE01 14.75, TE11/TM11 16.14 GHz.

EIGENSOLVE recipe (``helmholtz_cutoff_wavenumbers_2d``): assemble stiffness ``grad(u).grad(v)``
and mass ``u v``, solve the generalised eigenproblem with NGSolve's ``ArnoldiSolver`` using a
small NEGATIVE shift -- then A - shift*M = A + |shift|*M is SPD and never singular, and the
eigenvalues NEAREST the shift are the algebraically smallest (the lowest cutoffs). The Neumann
(TE) space always carries a trivial CONSTANT null mode (k_c ~ 0, no field) -- it is discarded;
pad ``nev`` to recover the wanted count after the drop. Order-3 H1 on maxh=min(a,b)/16 gives the
cutoffs to <0.1 %.

The SAME eigenproblem is the 2D cavity TM resonance and the vibrating-membrane (drum) spectrum;
only the reading of k_c changes (cutoff wavenumber / resonance / mode shape). Validated
tests/test_waveguide_cutoff.py: TE10/TE20/TE01/TM11
all match the exact spectrum to FE accuracy, and TM11 lies above every TE mode (the BC is the
physics). Stored regression reference cross-checked against an independent eigenvalue solver.

3-D CAVITY (full-wave Maxwell, #59): a closed PEC box a x b x d resonates at
``f_mnp = (c/2) sqrt((m/a)^2+(n/b)^2+(p/d)^2)`` -- the genuine VECTOR-EM cavity, ``curl curl E =
(omega/c)^2 E`` with ``n x E = 0`` on the walls, on an HCurl (edge) space:

```python
from radia_mcp.radia_ngsolve.waveguide import rectangular_cavity_frequency, maxwell_cavity_modes_3d
f101 = rectangular_cavity_frequency(a, b, d, 1, 0, 1)              # dominant TE101 = (c/2)sqrt(1/a^2+1/d^2)
shift = (pi/a)**2 + (pi/d)**2                                       # ~ fundamental k^2
freqs = maxwell_cavity_modes_3d(mesh, 4, shift)                     # HCurl curl-curl eigenvalue
```

CRITICAL (the vector trap): the curl-curl operator has a HUGE GRADIENT KERNEL (curl grad = 0 -> a
cloud of spurious ZERO modes). Do NOT search near 0 -- give ``ArnoldiSolver`` a ``shift`` near the
target k^2 so shift-invert amplifies the PHYSICAL resonances and starves the null space (then drop
any residual eigenvalue << shift). This is why an EM cavity needs HCurl + a shift, unlike the 2-D
scalar H1 cutoff above. Validated: TE101/TM110/TE011/
TE111 of a 40x20x30 mm box to < 0.1 % (order-2 HCurl), gradient modes cleanly suppressed.

CIRCULAR waveguide (#61): a ROUND guide's modes solve the SAME 2D Helmholtz eigenproblem on the
DISK cross-section, so ``helmholtz_cutoff_wavenumbers_2d`` runs on a disk mesh unchanged -- only
the closed form differs (Bessel zeros, not sin):

    TM_mn: k_c = j_mn/a  (n-th zero of J_m),   TE_mn: k_c = j'_mn/a  (n-th zero of J'_m).

``circular_waveguide_cutoff(a, 'TE'|'TM', m, n)`` (scipy ``jn_zeros``/``jnp_zeros``). Dominant TE11
(j'_11 = 1.8412 -> f_c = 0.293 c/a); next TM01 (j_01 = 2.4048). Circular guides have the famous
DEGENERACY TE_0n / TM_1n (because j'_0n = j_1n), and TE11 is a degenerate POLARISATION pair -- both
appear in the FE spectrum. Validated tests/test_circular_waveguide.py: TE11/TM01/TE21/TM11 == the Bessel-zero cutoffs to 0.017 %. The
circular sibling of the rectangular cutoff (#53) -- one eigensolver, two cross-sections.

DISPERSION above cutoff (#65): a cutoff EIGENVALUE only says WHETHER a mode propagates; the
propagation sequel says HOW. Above f_c the axial constant is beta = (2 pi/c) sqrt(f^2 - f_c^2), so
the GUIDE WAVELENGTH lambda_g = lambda0/sqrt(1-(f_c/f)^2) > lambda0 (obeying 1/lambda_g^2 =
1/lambda0^2 - 1/lambda_c^2), the PHASE velocity v_p = c/sqrt(1-(f_c/f)^2) > c (carries no energy),
the GROUP velocity v_g = c sqrt(1-(f_c/f)^2) < c (energy/signal), with the reciprocal identity
v_p * v_g = c^2. At cutoff v_g -> 0, v_p -> inf, lambda_g -> inf (a standing wave, no axial
transport); for f >> f_c both -> c (TEM-like). BELOW cutoff the mode is evanescent, exp(-alpha z),
alpha = (2 pi/c) sqrt(f_c^2 - f^2) -- the guide as a HIGH-PASS FILTER (alpha -> k_c as f -> 0).
``waveguide_dispersion(f, fc)`` (dict: beta, lambda_g, v_phase, v_group), ``guide_wavelength(f, fc)``,
``waveguide_evanescent_attenuation(f, fc)`` -- pure closed forms, while the cutoff fc itself is read
from the FE eigensolver above (the analytic dispersion curve rides on the FE-computed cutoff).
Validated tests/test_waveguide_dispersion.py:
v_p*v_g = c^2 and the lambda_g relation exact, v_g == d omega/d beta numerically, FE-anchored
fc to +0.000 %. Sets lambda_g/4 transformers, iris/slot spacings, and pulse-spreading (v_g) in
waveguide components.

S-PARAMETER post-processing (#221): the same propagation constants feed simple RF traces:
``waveguide_dielectric_slab_sparams`` for a TE10 dielectric section/cascade,
``waveguide_offset_short_s11`` for a shorted offset line, and ``reflection_metrics`` for return
loss / delivered power / VSWR.  For phase traces, ``sparameter_group_delay(frequencies, s_values)``
unwraps ``arg(S)`` and returns ``tau_g = -d arg(S)/d omega``; the offset-short check gives the
round-trip delay ``2d/v_group``.

Continuous-loop slot 143 adds the FAR-FIELD notebook handoff rule: after
``farfield_pattern_metadata_gate`` confirms frequency, degree/radian convention, theta/phi cuts,
polarization basis, field components, row count, and accepted-power normalization, a notebook lobe
row must still carry ``lobe_id``, ``theta_deg``, ``phi_deg``, ``gain_unit=dBi``,
``directivity_unit=dBi``, radiated/accepted powers, and the identity ``G = eta_rad * D``.  Use
``farfield_lobe_notebook_handoff_gate`` before plotting or ranking antenna/EMC lobes; negative
controls should include missing lobe identity, a phi/theta location outside the exported cut grid,
and gain above directivity.

Continuous-loop slot 144 mirrors the same rule in the MATLAB/Gypsilab teaching lane as
``educationalFarfieldLobeNotebookHandoff``.  The purpose is not performance; it is readability:
students see metadata, lobe identity, angular cut membership, polarization basis, and
``G = eta_rad * D`` in one MATLAB struct before a notebook panel plots the row.

Continuous-loop slot 145 adds the source-table variant of the same rule: a lobe row
should carry source dataset identity and gain/directivity expression identity together
with the public row checks.  This keeps docs/panel notebooks replayable without relying
on column position or a naked gain scalar.

TE10 conductor loss: finite-conductivity rectangular-guide walls are a matched lossy line section.
``rectangular_waveguide_te10_conductor_loss(f, a, b, sigma, length=None)`` integrates
``(R_s/2)|H_t|^2`` over the four walls and returns ``alpha_np_per_m``, ``alpha_db_per_m``,
``surface_resistance_ohm``, ``skin_depth_m``, and optional ``S21_mag`` / insertion loss for a
length. The checks to remember are: loss diverges near cutoff because transmitted power collapses,
insertion loss is linear in length in dB, and alpha scales as ``1/sqrt(sigma)`` at fixed frequency.

TE10 port NORMALIZATION: ``rectangular_waveguide_te10_port_normalization(f, a, b, power_w=1)``
returns the peak ``E_y`` / ``H_x`` / wall ``H_z`` amplitudes that carry the requested time-average
power. The closed-form anchor is ``P=(1/2) integral E_y H_x dA = a b E0^2/(4 Z_TE)`` with
``E_y=E0 sin(pi x/a)`` and ``H_x=E_y/Z_TE``.  `validation_waveguide_te10_port_normalization.py`
checks WR-90 1 W rows, power scaling, ``H_z/H_x=k_c/beta``, and frequency trends to roundoff.

QUALITY FACTOR / wall loss (#68): the lossless resonance (#59) only gives WHERE a cavity rings; the
finite-conductivity walls give HOW SHARPLY. The unloaded Q = omega U / P_wall (stored energy over
per-cycle wall dissipation). For the TE101 rectangular cavity

    Q = (k a d)^3 b eta / (2 pi^2 R_s (2 a^3 b + 2 b d^3 + a^3 d + a d^3)),

R_s = sqrt(omega mu/(2 sigma)) = 1/(sigma*delta) the wall surface resistance, eta = mu0 c. So
Q ~ (volume/surface)/delta ~ sqrt(sigma) -- a better conductor or a colder/SC wall (smaller R_s)
sharpens the resonance; copper X-band cavities sit at Q ~ 1e4, superconducting RF cavities reach 1e9+.
``rectangular_cavity_q(a, b, d, sigma)`` with ``surface_resistance(f, sigma)`` / ``skin_depth(f, sigma)``
-- the BRIDGE from the wave family to the skin depth (#45/#55: same delta = sqrt(2/(omega mu sigma))).
Validated: the closed form == the numerically-integrated TE101 field energy U and wall loss P to
machine precision, and Q(4 sigma) = 2 Q(sigma) (tests/test_cavity_quality_factor.py). The figure of merit for resonators, filters, and accelerator /
NMR RF cavities.
"""


NGSOLVE_SLOT_LEAKAGE = r"""
# Slot-leakage inductance of a machine slot (the leakage reactance)

The flux that crosses a stator/rotor slot from tooth to tooth and links the conductor WITHOUT
reaching the air gap -- the leakage reactance X_sigma=omega L every winding carries on top of
the magnetising reactance. Set by the slot GEOMETRY (not the iron), it limits the short-circuit
current and shapes the transient response. For a rectangular slot the specific (per-axial-metre)
permeance is

    lambda_s = h_c/(3 w_s)  +  h_o/w_o ,     L = mu0 N^2 l_stk lambda_s ,

h_c=conductor height, w_s=slot width, (h_o,w_o)=tooth-tip opening above the conductor.

```python
from radia_mcp.radia_ngsolve.solve import slot_leakage_permeance, slot_leakage_inductance
lam = slot_leakage_permeance(h_c, w_s, opening_height=h_o, opening_width=w_o)
L   = slot_leakage_inductance(h_c, w_s, axial_length=l_stk, turns=N, opening_height=h_o, opening_width=w_o)
```

The **1/3 factor is exact** and is the whole point: a uniform-current conductor builds the
cross-slot MMF LINEARLY (0 at the slot bottom to NI at the conductor top), so
B_x(y)=mu0(NI/w_s)(y/h_c); integrating B^2 over that TRIANGULAR field gives 1/3 (vs 1 for a
current-free region of the same shape -> the opening's h_o/w_o series term).

FE recipe (the well-posed 1-D-equivalent BVP): solve
``solve_planar_magnetostatic`` on the slot rectangle with uniform Jz; the high-permeability tooth
walls become NEUMANN edges (B enters the iron perpendicular -> natural BC), and the yoke flux
line at the slot TOP is DIRICHLET A=0. Then L' = 2*magnetic_energy_2d(B)/I^2. The 1/3 conductor
term and a SAME-WIDTH current-free extension reproduce the closed form to **0.000 %**; a NARROWED
opening (w_o<w_s) makes the FE sit a few % ABOVE the simple series form (tooth-shoulder fringing
the lumped permeance omits), so the closed form is then a documented LOWER bound -- the same
signed-leakage story as the gapped-core inductance (coil links extra window flux) and the PM
load-line (magnet flux short-cuts the gap). validation_test/radia_mcp/test_slot_leakage.py.

KEY GOTCHA: do NOT model the slot as an isolated current in a closed iron block -- a net slot
current with no return path produces a global field (and a meaningless energy). The leakage flux
must CLOSE through the teeth/yoke; the Neumann-walls + Dirichlet-top BVP encodes exactly that.
The motor-design companion to the winding factor (#51, the EMF) and the magnetic-circuit
inductance (#39, the air-gap term).
"""


NGSOLVE_DEEP_BAR = r"""
# Deep-bar (slot-conductor) AC resistance & reactance -- the deep-bar / double-cage effect

The AC sequel to the slot-leakage inductance (#54): a SOLID conductor filling an iron slot, at
frequency, has its current pushed toward the slot OPENING by the slot-leakage field -- a 1-D skin
effect over the slot DEPTH. Field's closed form (xi = h/delta, delta the skin depth):

    k_R = Rac/Rdc = xi (sinh2xi + sin2xi)/(cosh2xi - cos2xi)
    k_X = Lac/Ldc = (3/2xi)(sinh2xi - sin2xi)/(cosh2xi - cos2xi)

```python
from radia_mcp.radia_ngsolve.solve import (deep_bar_skin_ratio,
    deep_bar_resistance_factor, deep_bar_reactance_factor, slot_leakage_permeance, MU0)
xi  = deep_bar_skin_ratio(h, f, sigma)            # = h sqrt(pi f mu0 mu_r sigma)
Rac = deep_bar_resistance_factor(xi) / (sigma*w*h)             # per metre
Lac = deep_bar_reactance_factor(xi) * MU0*slot_leakage_permeance(h, w)   # mu0 h/3w * k_X
```

LIMITS: xi->0 -> k_R,k_X -> 1 (DC: Rdc=1/(sigma w h), the #54 slot leakage mu0 h/3w). xi>>1 ->
k_R -> xi (Rac ~ sqrt(f), current in a one-skin-depth surface layer) and k_X -> 3/(2 xi) -> 0.

This is the DEEP-BAR / double-cage rotor effect: at high slip the rotor frequency is high, so the
bar Rac RISES (high STARTING resistance/torque, #40 deep-bar IM) while the slot-leakage reactance
FALLS; near synchronous speed the frequency -> 0 and the bar reverts to its DC resistance (high
running efficiency). Also the AC copper loss of stator slot conductors at high electrical
frequency (high-speed machines).

FE recipe (current-driven eddy on the iron-slot 1-D BVP):
``solve_planar_eddy(mesh, NU0, sigma, omega, driven_region="bar", total_current=I,
dirichlet="top")`` -- the SAME slot BVP as #54 (top edge Dirichlet yoke / Neumann tooth walls)
but time-harmonic; read Z = Vc/I from the NumberSpace, Rac=Re(Z), Lac=Im(Z)/omega. radia FE ==
Field's k_R/k_X to **0.00 %** over xi = 0.7..4.3 (tests/test_deep_bar.py). The rectangular-bar /
slot sibling of the round-wire Kelvin skin effect (#45, ber/bei) -- different geometry, the
sinh/sin (not ber/bei) functions; both are the AC frequency axis of the DC inductance set.

DOWELL (#60, the m-LAYER winding): stack m foils in series in the window and the leakage field
builds up layer by layer -> the outer layers sit in a large field and dissipate heavily (PROXIMITY
loss on top of each foil's skin loss):

    F_R = Rac/Rdc = Delta[ (sinh2D+sin2D)/(cosh2D-cos2D) + (2/3)(m^2-1)(sinhD-sinD)/(coshD+cosD) ]

(Delta = h/delta, m = layers). ``dowell_resistance_factor(skin_ratio, layers)``; the first term is
exactly ``deep_bar_resistance_factor`` (m=1), the second is proximity, growing ~m^2 -- why thick
multi-layer windings are AC-loss disasters (m=6 at Delta=1.5 -> F_R=17.7) and why litz wire / thin
foils exist. FE = an m-conductor SERIES eddy circuit (``solve_planar_eddy_multi(..., connection=
'series', total_current=I)``, Rac=2*sum_k P_k/I^2) on the SAME slot BVP (top Dirichlet yoke / Neumann
walls); matches Dowell to **0.00 %** for m=1..6, Delta=1.0/1.5 (tests/test_dowell.py). The transformer/inductor winding-loss workhorse; the multi-layer generalisation
of the deep-bar (#55).
"""


NGSOLVE_PHI_THETA = r"""
# phi-theta (voltage-temperature) relation -- constriction supertemperature

A current-carrying conductor / electrical CONSTRICTION whose two terminals (at potentials U and 0)
are both held at T0 develops a hot spot whose temperature is GEOMETRY-INDEPENDENT:

    T(V)^2 = T0^2 + V (U - V)/L ,      T_max^2 - T0^2 = U^2/(4 L)   (at the V = U/2 surface),

L = Wiedemann-Franz Lorenz number = 2.44e-8 V^2/K^2. The temperature is a function of the local
POTENTIAL ALONE -- independent of the conductor's shape, size, conductivity, or current.

```python
from radia_mcp.radia_ngsolve.multiphysics import (
    phi_theta_max_temperature, phi_theta_local_temperature, melting_voltage)
Tmax = phi_theta_max_temperature(U, T0)          # sqrt(T0^2 + U^2/4L)
U_soft = melting_voltage(463.0, T0)              # ~0.11 V softens copper; ~0.43 V melts it (1356 K)
```

WHY it works (the elegant bit): under Wiedemann-Franz (k = L sigma T) the Kirchhoff transform
psi = L sigma T^2/2 turns the coupled electro-thermal problem into a LINEAR Poisson problem,
-lap(psi) = sigma|grad V|^2, whose solution is QUADRATIC in V -- so T^2 is quadratic in V with the
boundary values pinned to T0, giving T^2 = T0^2 + V(U-V)/L exactly.

FE recipe: solve V with
``scalar_fem2d.solve_poisson_2d(mesh, sigma, {"hi":U, "lo":0})`` (insulated sides), then the
Kirchhoff potential ``solve_poisson_2d(mesh, 1.0, {"hi":psi0,"lo":psi0}, source=sigma*|grad V|^2)``
with psi0 = L sigma T0^2/2, then T = sqrt(2 psi/(L sigma)). radia FE matches the closed form to
**0.000 %** on BOTH a uniform bar AND a dog-bone constriction -- the SAME U gives the SAME T_max
(the geometry independence) -- validation_test/radia_mcp/test_phi_theta.py. The Wiedemann-Franz-exact generalisation
of the constant-property Joule bar peak sigma V^2/8k (#19); the basis of contact-voltage limits
(connectors, relays, fuses, motor terminals -- a metal softens/melts at a voltage set by the
material alone).
"""


NGSOLVE_CARTER = r"""
# Carter coefficient -- slotting makes the air gap magnetically larger

The slot openings of a machine make the air gap behave LARGER: as the flux crosses a slot opening
it spreads (adding reluctance), so the average gap flux for a given MMF drops as if the gap were
k_C * g. The standard machine-design correction (Carter 1901):

    gamma = (b_o/g)^2/(5 + b_o/g),   k_C = tau_s/(tau_s - gamma g),   g_eff = k_C g,

tau_s = slot pitch, g = gap, b_o = slot opening.

```python
from radia_mcp.radia_ngsolve.solve import (
    carter_coefficient, effective_air_gap, slotted_air_gap_permeance_factor)
kc   = carter_coefficient(tau_s, g, b_o)             # >= 1
geff = effective_air_gap(tau_s, g, b_o)              # kc * g -> use in L_m calc
kp   = slotted_air_gap_permeance_factor(tau_s, g, b_o)  # P_slot/P_smooth = 1/kc
```

k_C feeds the magnetising inductance, the no-load flux, and the slot-ripple permeance ratio;
typical openings lengthen the gap a few % to ~30-40 %.

FE validation (the gap PERMEANCE): a SCALAR magnetic
potential phi_m solves Laplace on the air region (gap + a deep slot) with phi_m = F on the stator
iron surface (tooth face + slot walls/top), 0 on the rotor, and NEUMANN on the tooth-centre and
slot-centre symmetry planes (a half slot pitch). The permeance P = integral|grad phi_m|^2 (F=1,
the magnetic analogue of capacitance C=2W/V^2), and k_C = P_smooth/P_slot with P_smooth = (tau_s/2)/g.
``scalar_fem2d.solve_poisson_2d`` matches Carter's k_C to **< 0.3 %** over b_o/g = 1..6
(tests/test_carter.py). The static slotting companion to the slot-leakage inductance (#54, the
conductor-in-slot) and the deep-bar slot AC resistance (#55); the reluctance source of cogging.
"""


NGSOLVE_SKEW = r"""
# Skew factor -- axial skew suppresses EMF harmonics and cancels cogging

Skewing the stator slots (or rotor magnets) continuously over an electrical angle theta_sk
averages the n-th air-gap harmonic along the STACK, scaling the n-th EMF / torque harmonic by

    k_sk,n = sin(n theta_sk/2)/(n theta_sk/2)        (the sinc / averaging factor).

The AXIAL twin of the distribution factor (#51, which spreads a phase over slots in the slot
plane); here the conductor is spread along the stack length.

```python
from radia_mcp.radia_ngsolve.solve import (skew_factor, slot_pitch_skew_angle,
                                           skewed_winding_factor)
th = slot_pitch_skew_angle(Q, p)           # one slot pitch elec = 2 pi p/Q (the usual skew)
k_sk1 = skew_factor(1, th)                  # small fundamental cost (<1)
k_w_sk = skewed_winding_factor(n, q, slot_angle_deg, th)   # total = k_w,n * k_sk,n
```

KEY: a harmonic whose period equals the skew (n theta_sk = 2 pi) is NULLED. Skewing by ONE SLOT
PITCH (theta_sk = 2 pi p/Q) makes the SLOT-PASSING / cogging harmonic (order Q/p) see k_sk = 0 ->
the cogging torque is cancelled, at a small fundamental cost (e.g. 36-slot/4-pole: k_sk,1 = 0.995,
only 0.5 %, while the 5th/7th/11th/13th are further suppressed on top of the distribution factor).

Validated self-consistently against the phasor average |(1/theta_sk) integral exp(j n phi) dphi|
to 1e-3 and the exact null at n theta_sk = 2 pi (tests/test_skew_factor.py). Completes the
harmonic-mitigation set: distribution + pitch (#51, slot plane), Carter slotting (#57, gap),
skew (#58, stack). Skew is a 1-D axial AVERAGING, so it has no 2-D FE -- the phasor integral is
the ground truth (as for the winding factor #51).

Result-saved electric-machine demo: docs/electric_machine/cogging_skew_demo.ipynb runs the
finite-element reluctance-torque curve and applies MachineScaling. Its public values are
embedded in the notebook; checked order-2 ripple and skew_factor evidence lives in
validation_test/electric_machine/cogging_skew_demo_results.json with the executable
machine force/torque validation corpus.
"""


NGSOLVE_CORE_LOSS = r"""
LAMINATION EDDY-CURRENT LOSS / CORE LOSS (#66) -- why iron cores are laminated, and the f^2 'eddy'
term of the Steinmetz core-loss separation P_core = P_hyst(k_h f B^a) + P_eddy(k_e f^2 B^2).

A steel sheet of thickness d carrying a sinusoidal flux of peak B at frequency f dissipates the
CLASSICAL eddy loss (uniform flux, thin sheet d << skin depth):

    P_eddy = pi^2 sigma d^2 f^2 B_peak^2 / 6   [W/m^3]   (= sigma omega^2 d^2 B^2 / 24).

The d^2 law is the whole reason for lamination: halve the sheet thickness and you QUARTER the eddy
loss. ``classical_eddy_loss_density(sigma, d, f, B_peak)``.

SKIN-EFFECT roll-off: once d is no longer thin vs the skin depth delta = sqrt(2/(omega mu sigma)),
the flux cannot fully penetrate and the loss is reduced by the exact 1-D magnetic-diffusion factor

    F(xi) = (3/xi) (sinh xi - sin xi)/(cosh xi - cos xi),   xi = d/delta,

with F -> 1 as xi -> 0 (the classical d^2 law) and F -> 3/xi as xi -> inf (heavy skin effect: flux
excluded, loss grows only as d, so the classical formula OVER-predicts). Helpers
``lamination_eddy_skin_factor(d, delta)`` and ``lamination_eddy_loss_density(sigma, d, f, B_peak,
mu_r)`` (= classical * F, with delta from mu = mu0 mu_r).

This is the FIELD-DRIVEN counterpart of the CURRENT-DRIVEN deep-bar / Dowell AC-resistance family
(#55/#60): the SAME magnetic-diffusion equation, opposite drive (imposed flux vs imposed current), so
the same sinh/sin special functions appear -- here as a LOSS reduction rather than an R_ac rise.
Validated tool-independently: the closed-form F == the numerically-integrated exact slab loss to 1e-6,
and the thin limit == the classical formula
(validation_test/eddy_current_analytical_validation/test_lamination_eddy_loss.py).
Typical 0.35 mm electrical steel at 50/60 Hz sits in
the thin regime (xi ~ 0.5, F ~ 1); thick sheets or PWM-harmonic frequencies push xi up and F down.
For full motor core loss, add the hysteresis term (a material B-H model, not a field solve).
"""


NGSOLVE_NONLINEAR_CIRCUIT = r"""
NONLINEAR MAGNETIC CIRCUIT / SATURATING IRON-FRAME ELECTROMAGNET (#67) -- the saturating counterpart
of the constant-mu reluctance circuits (#27 gap force, #39 inductance, #42 PM load line).

A closed iron loop wound with N turns at current I reaches an operating flux density B set by
Ampere's law around the loop together with the iron B->H curve:

    NI = H_iron(B) * l_fe + (B/mu0) * gap .

With gap=0 every ampere-turn drops in the iron, so B sits on the BH KNEE (no air-gap term to hold it
back); a gap adds the large linear term B*gap/mu0 that DE-SATURATES the iron fast (a few tenths of a
millimetre of gap collapses B -- the leakage-free #27/#39/#42 story, now with a real saturating BH).
``magnetic_circuit_bh_operating_point(mmf, iron_path, h_of_b, gap)`` bisects for B (the right side is
monotone in B); ``h_of_b`` is the iron curve, e.g. H = 51 B + 2.5 B^15 (mu_r ~ 15600 unsaturated,
hard saturation above ~1.5 T).

Validated: self-consistent (Ampere residual -> 0), the constant-mu limit reproduces the linear
reluctance NI/(l_fe/mu0/mu_r + gap/mu0) to machine precision, and -- as a CROSS-TOOL anchor -- a live
nonlinear FE of a closed iron-frame electromagnet with this BH curve drives the iron to B ~ 1.56 T at
NI = 32 A-turns over a ~15.6 mm mean path (a stored regression reference), exactly the operating point
the lumped solver returns. The lumped B is the loop average; a 2-D FE adds corner/leakage detail (the
two legs carry +-B in opposite directions around the frame). The nonlinear-FE companion to the
constant-mu circuits and to the saturating-inductance knee (#28). tests/test_nonlinear_magnetic_circuit.py.
"""


NGSOLVE_MAGNETIZING_INDUCTANCE = r"""
MAGNETIZING (AIR-GAP) INDUCTANCE (#69) -- the dominant part of the dq L_md / L_mq: the per-phase
main-flux inductance an AC winding builds across the air gap.

From the fundamental MMF -> air-gap flux density -> flux per pole -> flux-linkage chain, the per-phase
SELF magnetizing inductance is

    L_mu = (2 mu0/pi) (kw1 N_ph)^2 D L_stk / (p^2 g_eff),

D = air-gap diameter, L_stk = stack length, g_eff = effective (Carter) gap, p = pole PAIRS, kw1 = the
fundamental winding factor, N_ph = series turns per phase. The SYNCHRONOUS magnetizing inductance
L_m = (m/2) L_mu (3/2 for three phases) adds the mutual coupling of the m phases sharing the common
rotating air-gap field; L_m is the dq L_md (= L_mq for a non-salient machine), and the terminal
Ld = L_m + L_leakage (slot leakage #54 + end-winding + harmonic leakage).

``magnetizing_inductance_per_phase(D, L_stk, g_eff, p, kw1, N_ph)`` /
``synchronous_magnetizing_inductance(..., phases=3)``. This block SYNTHESISES the campaign: kw1 from
the winding factor (#51), g_eff from Carter (#57); it scales as 1/g_eff (slotting lowers it), 1/p^2,
and (kw1 N_ph)^2 -- the machine-design levers. Validated against the explicit derivation chain to
machine precision (tests/test_magnetizing_inductance.py).
The air-gap (main) counterpart of the leakage inductances (slot #54, gapped-core #39).
"""


NGSOLVE_DEMAG_FACTOR = r"""
DEMAGNETIZING FACTOR & INTERNAL FIELD OF A MAGNETIZED BODY (#71) -- a permanent magnet sets up a
self-demagnetizing field H = -N M inside itself, so it only holds

    B_in = Br (1 - N),   H_in = -N Br/mu0,   N = demagnetizing factor (Nx+Ny+Nz = 1 for any ellipsoid).

Exact limiting shapes: sphere (N = 1/3 -> 2 Br/3); infinite cylinder TRANSVERSE (N = 1/2 -> Br/2) vs
AXIAL (N = 0 -> full Br); thin slab through-thickness (N = 1 -> ~0). So an axially-magnetized ROD
keeps its remanence while a transversely-magnetized DISC loses half -- the basis of shape anisotropy
and of the PM operating point / load line (#42, where the permeance coefficient plays N's role; the
larger N, the deeper into the second quadrant, the closer to the demag knee).
``demagnetizing_factor(shape)`` (1st component = the magnetized axis), ``magnetized_body_internal_field(
Br, N)``, ``demagnetizing_field(Br, N)``. Validated: a transverse-magnetized infinite cylinder solved
with radia ``magnets=`` (mu_r = 1) holds a UNIFORM internal B = Br/2 to ~0.2 % (open-domain
truncation, FE below); tests/test_demagnetizing_factor.py.
The uniform-magnetization counterpart of the iron magnetic circuits (#42 PM load line, #67 nonlinear).
"""


NGSOLVE_MMF_HARMONICS = r"""
WINDING MMF SPACE HARMONICS & ROTATING-FIELD SYNTHESIS (#72) -- the air-gap MMF a polyphase winding
lays down, and how its harmonics rotate. The n-th SPACE harmonic of one phase's MMF is

    F_phi,n = (4/pi) (k_w,n / n) (N_ph I / 2p),

the (4/pi)/n square-wave envelope of a full-pitch coil scaled by the winding factor k_w,n (#51). In a
balanced m-phase winding the harmonics combine (the rotating-field theorem) into ROTATING waves of
amplitude (m/2) F_phi,n: the fundamental rotates FORWARD, the triplen (n a multiple of m: 3rd, 9th)
CANCEL, and the rest split into forward (when n-1 is a multiple of m) and backward (when n+1 is). For
the 3-phase winding the ODD harmonics it carries are n = 6k +/- 1: forward 6k+1 (1, 7, 13), backward
6k-1 (5, 11). The rotating fundamental (3/2) F_phi,1 is the per-phase MMF behind the magnetizing
inductance L_m = (3/2) L_mu (#69); the backward 5th / forward 7th drive the parasitic asynchronous
torques and the rotor-surface / magnet-eddy losses. ``single_phase_mmf_harmonic(n, kw_n, N_ph, I, p)``,
``rotating_mmf_amplitude(n, m, F)``, ``mmf_harmonic_direction(n, m)``. Validated against the direct
rotating-field trig sum to machine precision (tests/test_mmf_harmonics.py). The MMF (excitation) companion of the winding factor (#51, EMF) and skew
(#58); the source term of the air-gap field that drives the magnetizing inductance and torque.
"""


NGSOLVE_ACDC_CROSS_LEARNING = r"""
AC/DC CROSS-LEARNING CATALOG -- public-safe 100-case learning queue distilled from private
source-native metadata.  Query the code artifact via
``radia_mcp.radia_ngsolve.acdc_cross_learning.public_acdc_problem_catalog()`` and gate it with
``acdc_problem_catalog_manifest_gate()`` before promoting cases to solver-ready validations.

The catalog deliberately stores only scrubbed engineering families, solver lanes, observables, and
validation gates.  Source-native titles, URLs, and benchmark numbers stay in the private
cross-validation lane.  Important families include coil inductance, electric-current conduction,
induction heating, transformer/magnetic circuits, rotating machines, capacitance matrices,
open-boundary BEM/high-order impedance boundaries, material models, transmission-line/cable fields,
eddy-current braking, and EM multiphysics coupling.

Representative next promotions:
  - coil cases -> inductance energy/flux-linkage reciprocity
  - eddy/heating cases -> skin-depth/loss monotonicity and energy balance
  - transformer cases -> open/short equivalent parameters and flux conservation
  - motor cases -> torque/coenergy identity plus AGE/HDiv lane agreement
  - electrostatic cases -> capacitance-matrix symmetry and positive semidefinite checks
  - open-boundary cases -> FEM/BEM reciprocity and the lab high-order impedance boundary policy
"""


def get_ngsolve_documentation(topic: str = "all") -> str:
    """Return NGSolve usage documentation by topic."""
    topics = {
        "coax": NGSOLVE_COAX_LINE,
        "transmission_line": NGSOLVE_COAX_LINE,
        "two_wire": NGSOLVE_COAX_LINE,
        "two_wire_inductance": NGSOLVE_COAX_LINE,
        "internal_inductance": NGSOLVE_COAX_LINE,
        "wire_inductance": NGSOLVE_COAX_LINE,
        "skin_effect": NGSOLVE_COAX_LINE,
        "ac_internal_inductance": NGSOLVE_COAX_LINE,
        "ac_resistance": NGSOLVE_COAX_LINE,
        "loop_inductance": NGSOLVE_COAX_LINE,
        "two_wire_loop_inductance": NGSOLVE_COAX_LINE,
        "wire_over_ground": NGSOLVE_COAX_LINE,
        "microstrip": NGSOLVE_COAX_LINE,
        "single_ended_line": NGSOLVE_COAX_LINE,
        "ground_plane": NGSOLVE_COAX_LINE,
        "transmission_line_inductance": NGSOLVE_COAX_LINE,
        "skin_effect_inductance": NGSOLVE_COAX_LINE,
        "kelvin_functions": NGSOLVE_COAX_LINE,
        "proximity_skin": NGSOLVE_COAX_LINE,
        "characteristic_impedance": NGSOLVE_COAX_LINE,
        "magnetic_energy": NGSOLVE_COAX_LINE,
        "inductance_energy": NGSOLVE_COAX_LINE,
        "waveguide": NGSOLVE_WAVEGUIDE,
        "cutoff_frequency": NGSOLVE_WAVEGUIDE,
        "cutoff": NGSOLVE_WAVEGUIDE,
        "rectangular_waveguide": NGSOLVE_WAVEGUIDE,
        "circular_waveguide": NGSOLVE_WAVEGUIDE,
        "bessel_cutoff": NGSOLVE_WAVEGUIDE,
        "te11": NGSOLVE_WAVEGUIDE,
        "round_waveguide": NGSOLVE_WAVEGUIDE,
        "te_mode": NGSOLVE_WAVEGUIDE,
        "tm_mode": NGSOLVE_WAVEGUIDE,
        "te10": NGSOLVE_WAVEGUIDE,
        "helmholtz_eigenvalue": NGSOLVE_WAVEGUIDE,
        "helmholtz_equation": NGSOLVE_WAVEGUIDE,
        "waveguide_mode": NGSOLVE_WAVEGUIDE,
        "eigenmode": NGSOLVE_WAVEGUIDE,
        "eigenmodes": NGSOLVE_WAVEGUIDE,
        "cavity_resonance": NGSOLVE_WAVEGUIDE,
        "cavity_mode": NGSOLVE_WAVEGUIDE,
        "rectangular_cavity": NGSOLVE_WAVEGUIDE,
        "resonant_cavity": NGSOLVE_WAVEGUIDE,
        "maxwell_eigenvalue": NGSOLVE_WAVEGUIDE,
        "curl_curl_eigenvalue": NGSOLVE_WAVEGUIDE,
        "hcurl_eigenvalue": NGSOLVE_WAVEGUIDE,
        "te101": NGSOLVE_WAVEGUIDE,
        "microwave_cavity": NGSOLVE_WAVEGUIDE,
        "membrane_mode": NGSOLVE_WAVEGUIDE,
        "single_mode_bandwidth": NGSOLVE_WAVEGUIDE,
        "waveguide_dispersion": NGSOLVE_WAVEGUIDE,
        "dispersion": NGSOLVE_WAVEGUIDE,
        "guide_wavelength": NGSOLVE_WAVEGUIDE,
        "group_velocity": NGSOLVE_WAVEGUIDE,
        "phase_velocity": NGSOLVE_WAVEGUIDE,
        "evanescent": NGSOLVE_WAVEGUIDE,
        "propagation_constant": NGSOLVE_WAVEGUIDE,
        "cutoff_attenuation": NGSOLVE_WAVEGUIDE,
        "quality_factor": NGSOLVE_WAVEGUIDE,
        "cavity_q": NGSOLVE_WAVEGUIDE,
        "q_factor": NGSOLVE_WAVEGUIDE,
        "surface_resistance": NGSOLVE_WAVEGUIDE,
        "skin_depth": NGSOLVE_WAVEGUIDE,
        "cavity_loss": NGSOLVE_WAVEGUIDE,
        "sparameter": NGSOLVE_WAVEGUIDE,
        "s_parameter": NGSOLVE_WAVEGUIDE,
        "s11": NGSOLVE_WAVEGUIDE,
        "s21": NGSOLVE_WAVEGUIDE,
        "return_loss": NGSOLVE_WAVEGUIDE,
        "group_delay": NGSOLVE_WAVEGUIDE,
        "sparameter_group_delay": NGSOLVE_WAVEGUIDE,
        "offset_short": NGSOLVE_WAVEGUIDE,
        "readable_fem": NGSOLVE_READABLE_P1_FEM,
        "p1_tetrahedron": NGSOLVE_READABLE_P1_FEM,
        "p1_tet": NGSOLVE_READABLE_P1_FEM,
        "p1_flux_trace": NGSOLVE_READABLE_P1_FEM,
        "fem_bem_trace": NGSOLVE_READABLE_P1_FEM,
        "acoustic": NGSOLVE_ACOUSTIC_BEM,
        "acoustics": NGSOLVE_ACOUSTIC_BEM,
        "acoustic_bem": NGSOLVE_ACOUSTIC_BEM,
        "helmholtz_dtn": NGSOLVE_ACOUSTIC_BEM,
        "acoustic_impedance": NGSOLVE_ACOUSTIC_BEM,
        "impedance_dtn": NGSOLVE_ACOUSTIC_BEM,
        "radiation_impedance": NGSOLVE_ACOUSTIC_BEM,
        "cross_validation": NGSOLVE_CROSS_VALIDATION_REGISTRY,
        "cross_validation_registry": NGSOLVE_CROSS_VALIDATION_REGISTRY,
        "validation_registry": NGSOLVE_CROSS_VALIDATION_REGISTRY,
        "xval_registry": NGSOLVE_CROSS_VALIDATION_REGISTRY,
        "loop_learning": NGSOLVE_CROSS_VALIDATION_REGISTRY,
        "acdc_cross_learning": NGSOLVE_ACDC_CROSS_LEARNING,
        "acdc_problem_catalog": NGSOLVE_ACDC_CROSS_LEARNING,
        "acdc_learning_catalog": NGSOLVE_ACDC_CROSS_LEARNING,
        "optuna_policy": NGSOLVE_CROSS_VALIDATION_REGISTRY,
        "force": NGSOLVE_ELECTROMAGNETIC_FORCE,
        "forces": NGSOLVE_ELECTROMAGNETIC_FORCE,
        "electromagnetic_force": NGSOLVE_ELECTROMAGNETIC_FORCE,
        "em_force": NGSOLVE_ELECTROMAGNETIC_FORCE,
        "maxwell_stress": NGSOLVE_ELECTROMAGNETIC_FORCE,
        "maxwell_stress_tensor": NGSOLVE_ELECTROMAGNETIC_FORCE,
        "mst": NGSOLVE_ELECTROMAGNETIC_FORCE,
        "eggshell_force": NGSOLVE_ELECTROMAGNETIC_FORCE,
        "weighted_stress": NGSOLVE_ELECTROMAGNETIC_FORCE,
        "nodal_force": NGSOLVE_ELECTROMAGNETIC_FORCE,
        "force_method": NGSOLVE_ELECTROMAGNETIC_FORCE,
        "force_extraction": NGSOLVE_ELECTROMAGNETIC_FORCE,
        "torque": NGSOLVE_ELECTROMAGNETIC_FORCE,
        "core_loss": NGSOLVE_CORE_LOSS,
        "eddy_loss": NGSOLVE_CORE_LOSS,
        "eddy_current_loss": NGSOLVE_CORE_LOSS,
        "lamination": NGSOLVE_CORE_LOSS,
        "lamination_loss": NGSOLVE_CORE_LOSS,
        "iron_loss": NGSOLVE_CORE_LOSS,
        "steinmetz": NGSOLVE_CORE_LOSS,
        "classical_eddy_loss": NGSOLVE_CORE_LOSS,
        "nonlinear_magnetic_circuit": NGSOLVE_NONLINEAR_CIRCUIT,
        "saturating_circuit": NGSOLVE_NONLINEAR_CIRCUIT,
        "saturating_iron_frame": NGSOLVE_NONLINEAR_CIRCUIT,
        "bh_operating_point": NGSOLVE_NONLINEAR_CIRCUIT,
        "magnetic_circuit_operating_point": NGSOLVE_NONLINEAR_CIRCUIT,
        "magnetizing_inductance": NGSOLVE_MAGNETIZING_INDUCTANCE,
        "magnetising_inductance": NGSOLVE_MAGNETIZING_INDUCTANCE,
        "air_gap_inductance": NGSOLVE_MAGNETIZING_INDUCTANCE,
        "main_inductance": NGSOLVE_MAGNETIZING_INDUCTANCE,
        "l_md": NGSOLVE_MAGNETIZING_INDUCTANCE,
        "lmd": NGSOLVE_MAGNETIZING_INDUCTANCE,
        "demagnetizing_factor": NGSOLVE_DEMAG_FACTOR,
        "demag_factor": NGSOLVE_DEMAG_FACTOR,
        "demagnetizing_field": NGSOLVE_DEMAG_FACTOR,
        "magnetized_body": NGSOLVE_DEMAG_FACTOR,
        "shape_anisotropy": NGSOLVE_DEMAG_FACTOR,
        "mmf_harmonics": NGSOLVE_MMF_HARMONICS,
        "mmf": NGSOLVE_MMF_HARMONICS,
        "rotating_field": NGSOLVE_MMF_HARMONICS,
        "space_harmonics": NGSOLVE_MMF_HARMONICS,
        "magnetomotive_force": NGSOLVE_MMF_HARMONICS,
        "rotating_mmf": NGSOLVE_MMF_HARMONICS,
        "electrothermal_2way": NGSOLVE_ELECTROTHERMAL_2WAY,
        "sigma_of_t": NGSOLVE_ELECTROTHERMAL_2WAY,
        "two_way_coupling": NGSOLVE_ELECTROTHERMAL_2WAY,
        "temperature_dependent": NGSOLVE_ELECTROTHERMAL_2WAY,
        "multiphysics": NGSOLVE_MULTIPHYSICS,
        "induction_heating": NGSOLVE_MULTIPHYSICS,
        "phi_theta": NGSOLVE_PHI_THETA,
        "voltage_temperature": NGSOLVE_PHI_THETA,
        "supertemperature": NGSOLVE_PHI_THETA,
        "constriction": NGSOLVE_PHI_THETA,
        "contact_resistance_heating": NGSOLVE_PHI_THETA,
        "melting_voltage": NGSOLVE_PHI_THETA,
        "wiedemann_franz": NGSOLVE_PHI_THETA,
        "kohlrausch": NGSOLVE_PHI_THETA,
        "holm_contact": NGSOLVE_PHI_THETA,
        "hot_spot": NGSOLVE_PHI_THETA,
        "kirchhoff_transform": NGSOLVE_PHI_THETA,
        "electrostatics": NGSOLVE_ELECTROSTATICS_3D,
        "capacitance": NGSOLVE_ELECTROSTATICS_3D,
        "field_quality": NGSOLVE_FIELD_QUALITY,
        "multipole": NGSOLVE_FIELD_QUALITY,
        "harmonics": NGSOLVE_FIELD_QUALITY,
        "field_errors": NGSOLVE_FIELD_QUALITY,
        "permanent_magnet": NGSOLVE_PM_AXIAL,
        "cylinder_magnet": NGSOLVE_PM_AXIAL,
        "axi_extraction": NGSOLVE_PM_AXIAL,
        "solenoid": NGSOLVE_SOLENOID,
        "coil": NGSOLVE_SOLENOID,
        "inductance": NGSOLVE_SOLENOID,
        "busbar": NGSOLVE_BUSBAR_LORENTZ,
        "lorentz_force": NGSOLVE_BUSBAR_LORENTZ,
        "method_of_images": NGSOLVE_BUSBAR_LORENTZ,
        "image_force": NGSOLVE_BUSBAR_LORENTZ,
        "image_method": NGSOLVE_BUSBAR_LORENTZ,
        "helmholtz": NGSOLVE_HELMHOLTZ,
        "uniform_field": NGSOLVE_HELMHOLTZ,
        "coil_pair": NGSOLVE_HELMHOLTZ,
        "c_magnet": NGSOLVE_C_MAGNET,
        "magnetic_circuit": NGSOLVE_C_MAGNET,
        "reluctance": NGSOLVE_C_MAGNET,
        "reluctance_actuator": NGSOLVE_C_MAGNET,
        "holding_force": NGSOLVE_C_MAGNET,
        "air_gap_force": NGSOLVE_C_MAGNET,
        "air_gap_holding_force": NGSOLVE_C_MAGNET,
        "maxwell_pressure": NGSOLVE_C_MAGNET,
        "gapped_core_inductance": NGSOLVE_C_MAGNET,
        "magnetic_circuit_inductance": NGSOLVE_C_MAGNET,
        "core_inductor": NGSOLVE_C_MAGNET,
        "magnetic_circuit_leakage": NGSOLVE_C_MAGNET,
        "pm_loadline": NGSOLVE_C_MAGNET,
        "load_line": NGSOLVE_C_MAGNET,
        "permeance_coefficient": NGSOLVE_C_MAGNET,
        "operating_point": NGSOLVE_C_MAGNET,
        "demagnetization": NGSOLVE_C_MAGNET,
        "electromagnet": NGSOLVE_C_MAGNET,
        "solenoid_force": NGSOLVE_C_MAGNET,
        "dipole_magnet": NGSOLVE_C_MAGNET,
        "joule_heating": NGSOLVE_JOULE_ELECTROTHERMAL,
        "electrothermal": NGSOLVE_JOULE_ELECTROTHERMAL,
        "cylinder_joule_heating": NGSOLVE_JOULE_ELECTROTHERMAL,
        "wire_heating": NGSOLVE_JOULE_ELECTROTHERMAL,
        "ampacity": NGSOLVE_JOULE_ELECTROTHERMAL,
        "convective_cooling": NGSOLVE_JOULE_ELECTROTHERMAL,
        "robin_bc": NGSOLVE_JOULE_ELECTROTHERMAL,
        "newton_cooling": NGSOLVE_JOULE_ELECTROTHERMAL,
        "convective_electrothermal": NGSOLVE_JOULE_ELECTROTHERMAL,
        "film_coefficient": NGSOLVE_JOULE_ELECTROTHERMAL,
        "cooling_fin": NGSOLVE_JOULE_ELECTROTHERMAL,
        "fin_efficiency": NGSOLVE_JOULE_ELECTROTHERMAL,
        "extended_surface": NGSOLVE_JOULE_ELECTROTHERMAL,
        "heat_sink": NGSOLVE_JOULE_ELECTROTHERMAL,
        "thermal_resistance": NGSOLVE_JOULE_ELECTROTHERMAL,
        "die_stack": NGSOLVE_JOULE_ELECTROTHERMAL,
        "junction_temperature": NGSOLVE_JOULE_ELECTROTHERMAL,
        "multilayer_thermal": NGSOLVE_JOULE_ELECTROTHERMAL,
        "electro_thermal": NGSOLVE_JOULE_ELECTROTHERMAL,
        "elasticity": NGSOLVE_ELASTICITY,
        "thermal_stress": NGSOLVE_ELASTICITY,
        "thermo_mechanical": NGSOLVE_ELASTICITY,
        "structural": NGSOLVE_ELASTICITY,
        "thermal_bending": NGSOLVE_ELASTICITY,
        "bimorph": NGSOLVE_ELASTICITY,
        "bimetal": NGSOLVE_ELASTICITY,
        "bimetallic": NGSOLVE_ELASTICITY,
        "thermostat": NGSOLVE_ELASTICITY,
        "biaxial_thermal_stress": NGSOLVE_ELASTICITY,
        "biaxial": NGSOLVE_ELASTICITY,
        "laminate": NGSOLVE_ELASTICITY,
        "laminate_residual_stress": NGSOLVE_ELASTICITY,
        "rule_of_mixtures": NGSOLVE_ELASTICITY,
        "composite_thermal_stress": NGSOLVE_ELASTICITY,
        "residual_stress": NGSOLVE_ELASTICITY,
        "magneto_mechanical": NGSOLVE_ELASTICITY,
        "beam_deflection": NGSOLVE_ELASTICITY,
        "actuator": NGSOLVE_ELASTICITY,
        "reluctance_torque": NGSOLVE_RELUCTANCE_TORQUE,
        "cogging": NGSOLVE_RELUCTANCE_TORQUE,
        "cogging_torque": NGSOLVE_RELUCTANCE_TORQUE,
        "motor_torque": NGSOLVE_RELUCTANCE_TORQUE,
        "rotor_sweep": NGSOLVE_RELUCTANCE_TORQUE,
        "stress_tensor_torque": NGSOLVE_RELUCTANCE_TORQUE,
        "electro_mechanical": NGSOLVE_MEMS_ELECTROMECH,
        "mems": NGSOLVE_MEMS_ELECTROMECH,
        "electrostatic_actuator": NGSOLVE_MEMS_ELECTROMECH,
        "pull_in": NGSOLVE_MEMS_ELECTROMECH,
        "electrostatic_force": NGSOLVE_MEMS_ELECTROMECH,
        "back_emf": NGSOLVE_BACKEMF_FLUX,
        "kt_ke": NGSOLVE_BACKEMF_FLUX,
        "torque_constant": NGSOLVE_BACKEMF_FLUX,
        "machine_constant": NGSOLVE_BACKEMF_FLUX,
        "backemf": NGSOLVE_BACKEMF_FLUX,
        "flux_linkage": NGSOLVE_BACKEMF_FLUX,
        "emf": NGSOLVE_BACKEMF_FLUX,
        "ld_lq": NGSOLVE_LDLQ,
        "ldlq": NGSOLVE_LDLQ,
        "inductance_saliency": NGSOLVE_LDLQ,
        "saliency": NGSOLVE_LDLQ,
        "mtpa": NGSOLVE_MTPA,
        "maximum_torque_per_ampere": NGSOLVE_MTPA,
        "current_angle": NGSOLVE_MTPA,
        "field_weakening": NGSOLVE_FIELD_WEAKENING,
        "mtpv": NGSOLVE_FIELD_WEAKENING,
        "maximum_torque_per_voltage": NGSOLVE_FIELD_WEAKENING,
        "base_speed": NGSOLVE_FIELD_WEAKENING,
        "constant_power": NGSOLVE_FIELD_WEAKENING,
        "operating_region": NGSOLVE_FIELD_WEAKENING,
        "cpsr": NGSOLVE_FIELD_WEAKENING,
        "characteristic_current_cpsr": NGSOLVE_FIELD_WEAKENING,
        "field_weakening_speed_capability": NGSOLVE_FIELD_WEAKENING,
        "induction_machine": NGSOLVE_INDUCTION_MACHINE,
        "induction_motor": NGSOLVE_INDUCTION_MACHINE,
        "torque_slip": NGSOLVE_INDUCTION_MACHINE,
        "breakdown_torque": NGSOLVE_INDUCTION_MACHINE,
        "slip": NGSOLVE_INDUCTION_MACHINE,
        "equivalent_circuit": NGSOLVE_INDUCTION_MACHINE,
        "power_angle": NGSOLVE_POWER_ANGLE,
        "load_angle": NGSOLVE_POWER_ANGLE,
        "pull_out": NGSOLVE_POWER_ANGLE,
        "synchronous_machine": NGSOLVE_POWER_ANGLE,
        "delta_curve": NGSOLVE_POWER_ANGLE,
        "dq_operating_point": NGSOLVE_DQ_OPERATING_POINT,
        "operating_point_dq": NGSOLVE_DQ_OPERATING_POINT,
        "power_factor": NGSOLVE_DQ_OPERATING_POINT,
        "dq_voltage": NGSOLVE_DQ_OPERATING_POINT,
        "terminal_voltage": NGSOLVE_DQ_OPERATING_POINT,
        "back_emf_constant": NGSOLVE_DQ_OPERATING_POINT,
        "pm_back_emf": NGSOLVE_DQ_OPERATING_POINT,
        "pm_emf_constant": NGSOLVE_DQ_OPERATING_POINT,
        "motor_efficiency": NGSOLVE_DQ_OPERATING_POINT,
        "short_circuit": NGSOLVE_SHORT_CIRCUIT,
        "characteristic_current": NGSOLVE_SHORT_CIRCUIT,
        "demagnetization_current": NGSOLVE_SHORT_CIRCUIT,
        "braking_torque": NGSOLVE_SHORT_CIRCUIT,
        "fault_current": NGSOLVE_SHORT_CIRCUIT,
        "winding_factor": NGSOLVE_WINDING_FACTOR,
        "distribution_factor": NGSOLVE_WINDING_FACTOR,
        "pitch_factor": NGSOLVE_WINDING_FACTOR,
        "chording": NGSOLVE_WINDING_FACTOR,
        "harmonic_winding": NGSOLVE_WINDING_FACTOR,
        "slot_leakage": NGSOLVE_SLOT_LEAKAGE,
        "leakage_inductance": NGSOLVE_SLOT_LEAKAGE,
        "leakage_reactance": NGSOLVE_SLOT_LEAKAGE,
        "slot_permeance": NGSOLVE_SLOT_LEAKAGE,
        "slot_inductance": NGSOLVE_SLOT_LEAKAGE,
        "slot_reactance": NGSOLVE_SLOT_LEAKAGE,
        "tooth_leakage": NGSOLVE_SLOT_LEAKAGE,
        "deep_bar": NGSOLVE_DEEP_BAR,
        "deep_bar_effect": NGSOLVE_DEEP_BAR,
        "double_cage": NGSOLVE_DEEP_BAR,
        "slot_ac_resistance": NGSOLVE_DEEP_BAR,
        "ac_resistance_factor": NGSOLVE_DEEP_BAR,
        "skin_effect_slot": NGSOLVE_DEEP_BAR,
        "field_formula": NGSOLVE_DEEP_BAR,
        "starting_torque": NGSOLVE_DEEP_BAR,
        "rotor_bar": NGSOLVE_DEEP_BAR,
        "dowell": NGSOLVE_DEEP_BAR,
        "dowell_equation": NGSOLVE_DEEP_BAR,
        "proximity_effect": NGSOLVE_DEEP_BAR,
        "winding_ac_resistance": NGSOLVE_DEEP_BAR,
        "winding_loss": NGSOLVE_DEEP_BAR,
        "foil_winding": NGSOLVE_DEEP_BAR,
        "litz": NGSOLVE_DEEP_BAR,
        "transformer_winding": NGSOLVE_DEEP_BAR,
        "carter": NGSOLVE_CARTER,
        "carter_coefficient": NGSOLVE_CARTER,
        "carter_factor": NGSOLVE_CARTER,
        "effective_air_gap": NGSOLVE_CARTER,
        "slotting": NGSOLVE_CARTER,
        "slot_opening": NGSOLVE_CARTER,
        "air_gap_permeance": NGSOLVE_CARTER,
        "skew": NGSOLVE_SKEW,
        "skew_factor": NGSOLVE_SKEW,
        "skewing": NGSOLVE_SKEW,
        "cogging_cancellation": NGSOLVE_SKEW,
        "skewed_winding": NGSOLVE_SKEW,
        "stack_skew": NGSOLVE_SKEW,
        "harmonic_suppression": NGSOLVE_SKEW,
        "frozen_permeability": NGSOLVE_FROZEN_PERM,
        "frozen_perm": NGSOLVE_FROZEN_PERM,
        "saturated_inductance": NGSOLVE_FROZEN_PERM,
        "saturating_inductance": NGSOLVE_FROZEN_PERM,
        "saturation_knee": NGSOLVE_FROZEN_PERM,
        "incremental_inductance": NGSOLVE_FROZEN_PERM,
        "incremental_inductance_matrix": NGSOLVE_FROZEN_PERM,
        "cross_saturation": NGSOLVE_FROZEN_PERM,
        "ldq": NGSOLVE_FROZEN_PERM,
        "reciprocity": NGSOLVE_FROZEN_PERM,
        "inductance_matrix": NGSOLVE_FROZEN_PERM,
        "coenergy": NGSOLVE_FROZEN_PERM,
        "co_energy": NGSOLVE_FROZEN_PERM,
        "superposition": NGSOLVE_FROZEN_PERM,
        "overview": NGSOLVE_OVERVIEW,
        "spaces": NGSOLVE_FE_SPACES,
        "maxwell": NGSOLVE_MAXWELL,
        "solvers": NGSOLVE_SOLVERS,
        "preconditioners": NGSOLVE_PRECONDITIONERS,
        "bem": NGSOLVE_BEM,
        "ngsolve_bem_50": NGSOLVE_BEM_VOL_50,
        "ngsolve_bem_100": NGSOLVE_BEM_VOL_50,
        "bem_50": NGSOLVE_BEM_VOL_50,
        "bem_100": NGSOLVE_BEM_VOL_50,
        "vol_bem_50": NGSOLVE_BEM_VOL_50,
        "vol_visualization": NGSOLVE_BEM_VOL_50,
        "vol_double_click": NGSOLVE_BEM_VOL_50,
        "sol_double_click": NGSOLVE_BEM_VOL_50,
        "vibroacoustic_drum": NGSOLVE_BEM_VOL_50,
        "drum": NGSOLVE_BEM_VOL_50,
        "curved_vol_geometry": NGSOLVE_BEM_VOL_50,
        "curved_vol": NGSOLVE_BEM_VOL_50,
        "mesh": NGSOLVE_MESH,
        "nonlinear": NGSOLVE_NONLINEAR,
        "pitfalls": NGSOLVE_PITFALLS,
        "linalg": NGSOLVE_LINALG,
        "formulations": NGSOLVE_EM_FORMULATIONS,
        "adaptive": NGSOLVE_ADAPTIVE,
        "darwin": NGSOLVE_DARWIN_SIBC,
        "esim": NGSOLVE_ESIM,
        "treecotree": NGSOLVE_TREE_COTREE,
        "pml": NGSOLVE_PML,
        "decomposition": NGSOLVE_DOMAIN_DECOMP,
        "material": NGSOLVE_MATERIAL_MODELING,
        "ironloss": NGSOLVE_IRON_LOSS,
        "practical": NGSOLVE_PRACTICAL_TECHNIQUES,
        "axisymmetric": NGSOLVE_AXISYMMETRIC,
        "accelerator": NGSOLVE_ACCELERATOR_MAGNET,
        "team7": NGSOLVE_TEAM7,
        "boolean_policy": NGSOLVE_BOOLEAN_POLICY,
        "boolean": NGSOLVE_BOOLEAN_POLICY,
        "glue": NGSOLVE_BOOLEAN_POLICY,
        "fuse": NGSOLVE_BOOLEAN_POLICY,
        "compound": NGSOLVE_BOOLEAN_POLICY,
        "trampoline": NGSOLVE_BOOLEAN_POLICY,
        "cln": NGSOLVE_CLN_CAUER,
        "cauer": NGSOLVE_CLN_CAUER,
        "ladder": NGSOLVE_CLN_CAUER,
        "mor": NGSOLVE_CLN_CAUER,
        "eddy_current_mor": NGSOLVE_CLN_CAUER,
        "tree_cotree": NGSOLVE_CLN_CAUER,
        "treecotree_gauge": NGSOLVE_CLN_CAUER,
        "kameari": NGSOLVE_CLN_CAUER,
    }

    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join(topics.values())
    elif topic in topics:
        return topics[topic]
    else:
        return (
            f"Unknown topic: '{topic}'. "
            f"Available: all, {', '.join(topics.keys())}"
        )
