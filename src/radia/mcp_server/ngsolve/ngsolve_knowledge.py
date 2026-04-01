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
  - S:\\NGSolve\\EMPY\\EMPY_Analysis (K. Sugahara's EM formulations)
"""

NGSOLVE_OVERVIEW = """
# NGSolve Overview

NGSolve is a high-performance finite element library with Python bindings.
Key features relevant to electromagnetic simulation:

## Installation

```bash
pip install ngsolve   # Use official PyPI ngsolve >= 6.2.2601
# Version 6.2.2601 includes: MKL + PARDISO + Periodic BC fix
# Periodic BC regression (6.2.2406~6.2.2501) is FIXED in 6.2.2601+
# See: https://forum.ngsolve.org/t/3805
# Optional: pip install ngsolve-sparsesolv  # Compact AMS/COCR solvers
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

## Gmsh Mesh Import (ReadGmsh)

Load external Gmsh v2 format meshes (from Cubit, Gmsh, or other generators):

```python
from netgen.read_gmsh import ReadGmsh
from ngsolve import Mesh

# Load Gmsh v2 mesh (.msh)
ngmesh = ReadGmsh("model.msh")
mesh = Mesh(ngmesh)

# Physical groups become material names and boundary names
print("Materials:", mesh.GetMaterials())
print("Boundaries:", mesh.GetBoundaries())
```

**Key points:**
- Physical groups in Gmsh map to `mesh.GetMaterials()` (volumes) and
  `mesh.GetBoundaries()` (surfaces)
- Supports 1st and 2nd order elements (TET4/TET10, HEX8/HEX20, etc.)
- Use `dx("region_name")` and `ds("boundary_name")` in weak forms
- Material properties via `CoefficientFunction([... for mat in mesh.GetMaterials()])`

```python
# Material properties from Gmsh physical groups
mu_d = {"iron": 4e-7*3.14*1000, "air": 4e-7*3.14, "coil": 4e-7*3.14}
mu = CoefficientFunction([mu_d.get(mat, 4e-7*3.14) for mat in mesh.GetMaterials()])
```

## VTK Output for ParaView

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

## 13. VTK Export

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

**Verified**: Wire + iron cylinder benchmark (2D, linear, mu_r=1000):
- MST and nodal force agree within 0.001% (both use same FEM field)
- Both match analytical to 0.01% (with mesh.Curve(3))
- Without mesh.Curve(): ~9% systematic error for BOTH methods equally
- Accuracy difference between methods: none (linear problem)
See `Radia/examples/nodal_force/demo_nodal_force_2d.py`.

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
Cubit → hex mesh + export_NGSolveCurvedMesh(order=N)
NGSolve → Omega-reduced Omega + Kelvin + Energy-based B-input Play
```

## Key Components

### 1. CoilBuilder (radia.radia_coil_builder)
```python
from radia.radia_coil_builder import CoilBuilder

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
from cubit_netgen_bridge import extract_curved_mesh
from ngsolve import Mesh
mesh = Mesh(extract_curved_mesh(cubit, order=3))
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
- docs/research/ENERGY_BASED_B_INPUT_PLAY.md
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
from ngsolve.la import CompactAMSPreconditioner, COCRSolver

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


def get_ngsolve_documentation(topic: str = "all") -> str:
    """Return NGSolve usage documentation by topic."""
    topics = {
        "overview": NGSOLVE_OVERVIEW,
        "spaces": NGSOLVE_FE_SPACES,
        "maxwell": NGSOLVE_MAXWELL,
        "solvers": NGSOLVE_SOLVERS,
        "preconditioners": NGSOLVE_PRECONDITIONERS,
        "bem": NGSOLVE_BEM,
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
