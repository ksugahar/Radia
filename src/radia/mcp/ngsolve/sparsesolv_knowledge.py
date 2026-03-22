"""
ngsolve-sparsesolv knowledge base for MCP server.

Repository: https://github.com/ksugahar/ngsolve-sparsesolv
Version: 3.1.1
License: MPL 2.0
Based on: JP-MARs/SparseSolv

Import: from sparsesolv_ngsolve import CompactAMSPreconditioner, COCRSolver, etc.
Install: pip install sparsesolv-ngsolve (PyPI standalone add-on package)
Requires: ngsolve >= 6.2.2601
"""

SPARSESOLV_OVERVIEW = """
# ngsolve-sparsesolv

## What It Is

A standalone PyPI add-on package for NGSolve providing iterative solvers
and preconditioners. Install: `pip install sparsesolv-ngsolve`

```python
import sparsesolv_ngsolve as ssn
# or:
from sparsesolv_ngsolve import (
    SparseSolvSolver, ICPreconditioner, SGSPreconditioner,
    CompactAMSPreconditioner, ComplexCompactAMSPreconditioner,
    CompactAMGPreconditioner, COCRSolver, GMRESSolver,
)
```

**Key capabilities:**

| Capability           | Description                                       |
|----------------------|---------------------------------------------------|
| IC preconditioner    | IC with auto-shift, diagonal scaling              |
| SGS preconditioner   | Standalone SGS and SGS-MRTR                       |
| Compact AMS          | HCurl auxiliary-space Maxwell (header-only, no HYPRE) |
| Compact AMG          | Algebraic multigrid (PMIS + classical interp)     |
| COCR solver          | Complex-symmetric short recurrence (3.5x vs GMRES)|
| ABMC ordering        | Parallel triangular solve with level scheduling   |
| Divergence detection | Stagnation-based early termination                |

Repository: https://github.com/ksugahar/ngsolve-sparsesolv
Based on: JP-MARs/SparseSolv (https://github.com/JP-MARs/SparseSolv)

## Architecture

- **Header-only C++17 core**: `include/sparsesolv/` (no separate .cpp compilation)
- **Standalone PyPI package**: `pip install sparsesolv-ngsolve`, lightweight add-on
- **NGSolve integration**: Links against NGSolve's `SparseMatrix`, `BaseVector`,
  `BitArray`, `BilinearForm`, `FESpace` for seamless interop
- **Parallel backend**: Compile-time dispatch to NGSolve TaskManager, OpenMP, or serial

## Solver Methods

| Method    | Description                                      | Best For                        |
|-----------|--------------------------------------------------|---------------------------------|
| ICCG      | CG + Incomplete Cholesky preconditioning         | SPD and semi-definite systems   |
| SGSMRTR   | Min. Residual Truncated + Symmetric Gauss-Seidel | When IC factorization fails     |
| CG        | Conjugate Gradient (no preconditioning)          | Well-conditioned SPD systems    |
| COCR      | Conjugate Orth. Conjugate Residual (C++ native)  | **Complex-symmetric (A^T=A)**   |
| GMRES     | GMRES(40) with MGS orthogonalization             | Non-symmetric systems           |

## Preconditioners

| Preconditioner    | Description                                                       |
|-------------------|-------------------------------------------------------------------|
| IC                | Incomplete Cholesky with auto-shift for semi-definite systems     |
| SGS               | Symmetric Gauss-Seidel, factorization-free                        |
| **Compact AMS**   | Auxiliary-space Maxwell Solver for HCurl (header-only, no HYPRE)  |
| **Compact AMG**   | Algebraic multigrid (PMIS + classical interp + l1-Jacobi)        |

## Key Features

- **Add-on to NGSolve**: Independent module, no NGSolve source modification needed
- Auto-shift IC decomposition for semi-definite matrices (curl-curl problems)
- Complex-valued system support (auto-dispatch real/complex)
- Conjugate inner product switch for Hermitian vs complex-symmetric systems
- Best-iterate tracking (returns best solution even if not fully converged)
- Divergence/stagnation detection with early termination
- FreeDofs (Dirichlet boundary condition) support
- Parallel triangular solves via level scheduling or ABMC ordering
- ABMC (Algebraic Block Multi-Color) ordering for enhanced parallel triangular solves
- RCM (Reverse Cuthill-McKee) bandwidth reduction preprocessing for ABMC
- **Compact AMS** for HCurl: header-only, TaskManager-native, no HYPRE dependency
- **Compact AMG** for H1 subspace problems (PMIS coarsening + classical interpolation)
- **COCR** solver: 3.5x faster than GMRES for complex-symmetric eddy current
- Newton-optimized IC setup: skip reordering when sparsity pattern unchanged
- Fused SpMV+dot and AXPY+norm kernels: 5 kernels/CG-iteration (was 9 naive)
"""

SPARSESOLV_API = """
# Python API Reference

## Import

```python
# sparsesolv-ngsolve is a standalone PyPI package
from sparsesolv_ngsolve import (
    SparseSolvSolver,      # All-in-one solver factory
    ICPreconditioner,      # IC preconditioner for use with NGSolve CGSolver
    SGSPreconditioner,     # SGS preconditioner for use with NGSolve CGSolver
    SparseSolvResult,      # Result container

    # Compact AMS/AMG (HCurl/H1, header-only, no HYPRE)
    CompactAMSPreconditioner,            # Real AMS for HCurl
    ComplexCompactAMSPreconditioner,     # Complex Re/Im TaskManager parallel
    CompactAMGPreconditioner,            # H1 AMG (PMIS + classical interp)

    COCRSolver,            # COCR for complex-symmetric (C++ native)
    GMRESSolver,           # GMRES(40) for non-symmetric systems
)
```

Note: `import ngsolve` must be done before `import sparsesolv_ngsolve`.

## SparseSolvSolver

Factory function that auto-detects real/complex from mat.IsComplex().

### Constructor

```python
solver = SparseSolvSolver(
    mat,                          # NGSolve SparseMatrix (required)
    method="ICCG",                # "ICCG", "SGSMRTR", "CG", "COCR"
    freedofs=None,                # BitArray for Dirichlet BCs (optional)
    tol=1e-10,                    # Relative convergence tolerance
    maxiter=1000,                 # Maximum iterations
    shift=1.05,                   # IC preconditioner shift parameter
    save_best_result=True,        # Track best solution during iteration
    save_residual_history=False,  # Record residual at each iteration
    printrates=False,             # Print convergence info
    conjugate=False,              # True for Hermitian, False for complex-symmetric
    use_abmc=False,               # Enable ABMC ordering for parallel triangular solves
    abmc_block_size=4,            # Rows per block in ABMC aggregation
    abmc_num_colors=4,            # Target number of colors for ABMC coloring
    abmc_use_rcm=False,           # RCM bandwidth reduction before ABMC
)
```

### Properties (read/write)

| Property              | Type  | Default  | Description                              |
|-----------------------|-------|----------|------------------------------------------|
| method                | str   | "ICCG"   | Solver method                            |
| tol                   | float | 1e-10    | Relative tolerance                       |
| maxiter               | int   | 1000     | Maximum iterations                       |
| shift                 | float | 1.05     | IC shift parameter                       |
| save_best_result      | bool  | True     | Track best solution                      |
| save_residual_history | bool  | False    | Record residual history                  |
| auto_shift            | bool  | False    | Auto-adjust for semi-definite systems    |
| diagonal_scaling      | bool  | False    | Enable 1/sqrt(A[i,i]) scaling            |
| divergence_check      | bool  | False    | Stagnation-based early termination       |
| divergence_threshold  | float | 1000.0   | Stagnation detection multiplier          |
| conjugate             | bool  | False    | Conjugated inner product (Hermitian)     |
| divergence_count      | int   | 100      | Bad iterations before stopping           |
| printrates            | bool  | False    | Print convergence info                   |
| use_abmc              | bool  | False    | ABMC ordering for parallel tri. solves   |
| abmc_block_size       | int   | 4        | Rows per block in ABMC aggregation       |
| abmc_num_colors       | int   | 4        | Target colors for ABMC graph coloring    |
| abmc_use_rcm          | bool  | False    | RCM bandwidth reduction before ABMC      |

### Methods

```python
# Solve with initial guess
result = solver.Solve(f.vec, gfu.vec)

# Operator interface (zero initial guess)
gfu.vec.data = solver * f.vec
```

### Usage

```python
# Operator interface (zero initial guess)
gfu.vec.data = solver * f.vec

# Solve with initial guess
result = solver.Solve(f.vec, gfu.vec)

# Check result
result = solver.last_result  # SparseSolvResult
print(result.converged, result.iterations, result.final_residual)
```

## ICPreconditioner

```python
pre = ICPreconditioner(mat, freedofs=None, shift=1.05)
pre.Update()  # Must call to compute factorization

# Use with NGSolve CGSolver
from ngsolve.krylovspace import CGSolver
inv = CGSolver(a.mat, pre, tol=1e-10, maxiter=2000)
gfu.vec.data = inv * f.vec
```

## SGSPreconditioner

```python
pre = SGSPreconditioner(mat, freedofs=None)
pre.Update()  # Must call to set up preconditioner

# Use with NGSolve CGSolver
inv = CGSolver(a.mat, pre, tol=1e-10)
gfu.vec.data = inv * f.vec
```

## SparseSolvResult

```python
result.converged          # bool
result.iterations         # int
result.final_residual     # float
result.residual_history   # list[float] (if save_residual_history=True)
```
"""

SPARSESOLV_EXAMPLES = """
# Usage Examples

## 1. Basic Poisson Problem with ICCG

```python
from ngsolve import *
from sparsesolv_ngsolve import SparseSolvSolver

mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
fes = H1(mesh, order=2, dirichlet="bottom|right|top|left")
u, v = fes.TnT()

a = BilinearForm(fes)
a += grad(u) * grad(v) * dx
a.Assemble()

f = LinearForm(fes)
f += 2 * pi * pi * sin(pi * x) * sin(pi * y) * v * dx
f.Assemble()

solver = SparseSolvSolver(a.mat, method="ICCG",
                           freedofs=fes.FreeDofs(),
                           tol=1e-10, maxiter=2000)

gfu = GridFunction(fes)
gfu.vec.data = solver * f.vec
print(f"ICCG: {solver.last_result.iterations} iterations")
```

## 2. Auto-Shift for Curl-Curl (Electromagnetic) Problems

```python
from netgen.occ import Box, Pnt

box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
for face in box.faces:
    face.name = "outer"
mesh = box.GenerateMesh(maxh=0.4)

fes = HCurl(mesh, order=1, dirichlet="outer", nograds=True)
u, v = fes.TnT()

a = BilinearForm(fes)
a += curl(u) * curl(v) * dx
a.Assemble()

f = LinearForm(fes)
f += CF((0, 0, 1)) * v * dx
f.Assemble()

solver = SparseSolvSolver(a.mat, method="ICCG",
                           freedofs=fes.FreeDofs(),
                           tol=1e-8, maxiter=2000, shift=1.0)
solver.auto_shift = True         # Critical for semi-definite curl-curl
solver.diagonal_scaling = True

gfu = GridFunction(fes)
gfu.vec.data = solver * f.vec
assert solver.last_result.converged
```

## 3. Complex Eddy Current Problem (Complex-Symmetric)

```python
fes = HCurl(mesh, order=1, complex=True, dirichlet="outer", nograds=True)
u, v = fes.TnT()

a = BilinearForm(fes)
a += InnerProduct(curl(u), curl(v)) * dx
a += 1j * InnerProduct(u, v) * dx    # complex-symmetric (A^T = A)
a.Assemble()

f = LinearForm(fes)
f += InnerProduct(CF((1, 0, 0)), v) * dx
f.Assemble()

# Complex-symmetric: conjugate=False (default)
solver = SparseSolvSolver(a.mat, method="ICCG",
                           freedofs=fes.FreeDofs(),
                           tol=1e-8, maxiter=5000,
                           save_best_result=True)
gfu = GridFunction(fes)
gfu.vec.data = solver * f.vec
```

## 3b. Hermitian System (conjugate=True)

```python
# Real-coefficient problem in complex FE space => Hermitian (A^H = A)
fes = H1(mesh, order=2, complex=True, dirichlet="bottom|right|top|left")
u, v = fes.TnT()

a = BilinearForm(fes)
a += grad(u) * grad(v) * dx + u * v * dx  # real coefficients => Hermitian
a.Assemble()

# Hermitian: use conjugate=True (conjugated inner product a^H * b)
solver = SparseSolvSolver(a.mat, method="ICCG",
                           freedofs=fes.FreeDofs(),
                           tol=1e-10, conjugate=True)
gfu = GridFunction(fes)
gfu.vec.data = solver * f.vec
```

## 4. Divergence Detection

```python
solver = SparseSolvSolver(a.mat, method="CG",
                           freedofs=fes.FreeDofs(),
                           tol=1e-12, maxiter=5000)
solver.divergence_check = True
solver.divergence_threshold = 10.0
solver.divergence_count = 20

result = solver.Solve(f.vec, gfu.vec)
# Early termination if stagnating
```

## 5. IC Preconditioner with NGSolve CGSolver

```python
from ngsolve.krylovspace import CGSolver
from ngsolve.la import ICPreconditioner

pre = ICPreconditioner(a.mat, freedofs=fes.FreeDofs(), shift=1.05)
pre.Update()
inv = CGSolver(a.mat, pre, printrates=True, tol=1e-10, maxiter=2000)
gfu.vec.data = inv * f.vec
```

## 6. Residual History Recording

```python
solver = SparseSolvSolver(a.mat, method="ICCG",
                           freedofs=fes.FreeDofs(),
                           tol=1e-10, maxiter=2000,
                           save_residual_history=True)
result = solver.Solve(f.vec, gfu.vec)
print(f"History: {len(result.residual_history)} entries")
# Plot: plt.semilogy(result.residual_history)
```

## 7. ABMC Ordering for Parallel Triangular Solves

```python
# Enable ABMC ordering in IC preconditioner
# Beneficial for problems with >25,000 DOFs and multi-threaded execution
solver = SparseSolvSolver(a.mat, method="ICCG",
                           freedofs=fes.FreeDofs(),
                           tol=1e-10, maxiter=2000,
                           use_abmc=True,
                           abmc_block_size=4,
                           abmc_num_colors=4)

# Also works with auto-shift and diagonal scaling
solver.auto_shift = True
solver.diagonal_scaling = True

gfu = GridFunction(fes)
gfu.vec.data = solver * f.vec
```

"""

SPARSESOLV_ABMC = """
# ABMC Ordering: Parallel Triangular Solve Optimization

## Overview

ABMC (Algebraic Block Multi-Color) ordering enables parallel forward/backward
substitution in IC and SGS preconditioners. Without ABMC, triangular solves
are inherently sequential -- the main bottleneck for multi-threaded IC/SGS.

Algorithm (Iwashita, Nakashima, Takahashi, IPDPS 2012):
1. **BFS Aggregation**: Groups nearby rows into blocks of ~block_size rows
2. **Graph Coloring**: Colors the block adjacency graph so blocks of the same
   color have no data dependencies
3. **Parallel Execution**: Colors processed sequentially, blocks within each
   color processed in parallel

## When to Use ABMC

| DOFs     | Threads | ABMC Benefit    |
|----------|---------|-----------------|
| <10,000  | Any     | No benefit      |
| 10K-25K  | 4-8     | Marginal        |
| >25,000  | 4+      | 1.2-1.8x faster |
| >100,000 | 8+      | Significant     |

**Crossover point**: ~25,000 DOFs with 8 threads. Below this, the overhead of
reordering exceeds the parallelism benefit.

## Performance Benchmarks (8 threads)

| Problem                      | DOFs    | Without ABMC | With ABMC | Speedup |
|------------------------------|---------|--------------|-----------|---------|
| Poisson (unit cube, H1 p=2) | 186,219 | 111 ms/iter  | 87 ms     | 1.28x   |
| Toroidal coil (HCurl p=1)   | 148,439 | 145 ms/iter  | 81 ms     | 1.80x   |

The toroidal coil shows higher speedup because HCurl mesh topology creates
longer dependency chains that benefit more from multi-color parallelism.

## Configuration Parameters

| Parameter            | Default | Description                                  |
|----------------------|---------|----------------------------------------------|
| use_abmc             | False   | Enable ABMC ordering                         |
| abmc_block_size      | 4       | Rows per block (2-16 typical)                |
| abmc_num_colors      | 4       | Target colors (actual may be higher)         |
| abmc_use_rcm         | False   | RCM bandwidth reduction before ABMC          |

## RCM Preprocessing

RCM (Reverse Cuthill-McKee) bandwidth reduction is available via `abmc_use_rcm=True`.
It reorders the matrix to reduce bandwidth before ABMC coloring, improving cache
locality for both SpMV and triangular solves.

## Known Limitations

- `abmc_num_colors` is a lower bound; the actual number of colors depends on
  the graph structure and may be higher for complex meshes
- `abmc_reorder_spmv=False` (default) is usually faster because FEM mesh
  ordering has better cache locality for SpMV than ABMC ordering

## Implementation Details

The ABMC pipeline in the IC preconditioner:
1. **Setup phase** (once): Build ABMCSchedule, reorder L factor, pre-allocate work vectors
2. **Apply phase** (each iteration): Permute input -> parallel triangular solve -> unpermute output

Parallel execution uses a **persistent parallel region** with SpinBarrier
synchronization (one `parallel_for(nthreads, ...)` call with barriers between
color levels), avoiding per-level thread dispatch overhead.

Fused-dot optimization: The preconditioner apply can fuse dot(r, z) computation
into the backward substitution phase, eliminating one kernel launch per CG iteration
(5 -> 4 kernels/iteration).
"""

SPARSESOLV_BEST_PRACTICES = """
# Best Practices and Recommendations

## Preconditioner Selection Guide

| Problem Type                    | Recommended                    | Why                              |
|---------------------------------|--------------------------------|----------------------------------|
| H1 Poisson/Elasticity (small)  | IC (shift=1.05)                | Fast, reliable                   |
| H1 Poisson/Elasticity (large)  | Compact AMG                    | Header-only AMG, no HYPRE needed |
| HCurl curl-curl (semi-definite)| IC + auto_shift                | Handles kernel automatically     |
| HCurl magnetostatics (real)    | **Compact AMS + CG**           | Mesh-independent iterations      |
| HCurl eddy current (complex)   | **Compact AMS + COCR**         | Best: 52 iters stable at any scale |
| HCurl eddy current (complex)   | ICCG (fallback)                | Simple but O(h^-1) iteration growth |
| Large parallel (H1)            | ABMC + IC                      | Parallel triangular solve        |
| Any problem, IC fails          | SGS-MRTR                       | No factorization needed          |

## Complex-Valued Problems: conjugate Setting

This is the most common source of confusion:

- **Complex-symmetric** (A^T = A): `conjugate=False` (default)
  - Examples: eddy current with `1j * sigma * u * v * dx`
  - Inner product: a^T * b (unconjugated)

- **Hermitian** (A^H = A): `conjugate=True`
  - Examples: real-coefficient problem assembled in complex FE space
  - Inner product: a^H * b (conjugated)

Using the wrong setting causes divergence or extremely slow convergence
(e.g., 5000 iterations instead of 58).

## Curl-Curl (HCurl) Problems

1. **Always use `nograds=True`**: Fewer DOFs, better conditioned
2. **Source construction matters**: `int_coil J.v dx` (domain integral) may not
   be discretely div-free, causing IC convergence failure
3. **Safe source approaches**:
   - curl-based: `int T.curl(v) dx`
   - Helmholtz correction: Remove gradient component via Laplacian solve
   - HDiv projection: Mixed FEM to enforce div=0
4. **Eddy current**: The `sigma * u * v * dx` term naturally regularizes

## Shift Parameter for IC

- Default `shift=1.05` works for most SPD problems
- For semi-definite (curl-curl): use `auto_shift=True` with `shift=1.0`
- Higher shift = more stable but slower convergence
- `diagonal_scaling=True` often helps for poorly scaled matrices
- Auto-shift uses exponential backoff (increment *= 2), so restart count is O(log n)

## ABMC Ordering

- Only enable for >25K DOFs with multi-threaded execution
- Default `abmc_block_size=4, abmc_num_colors=4` works well for most cases
- Consider `abmc_use_rcm=True` for bandwidth reduction

## Solver Selection for Non-Symmetric Systems

| Solver   | Memory       | Best For                          |
|----------|-------------|-----------------------------------|
| COCR     | 5 vectors   | Complex-symmetric (A^T = A)       |
| GMRES(40)| 40+ vectors | Non-symmetric, robust convergence |
"""

SPARSESOLV_BUILD = """
# Installation

## PyPI (Recommended)

```bash
pip install ngsolve>=6.2.2601     # Required dependency
pip install sparsesolv-ngsolve    # Standalone add-on package
```

Verify:

```bash
python -c "from sparsesolv_ngsolve import SparseSolvSolver, CompactAMSPreconditioner; print('OK')"
```

## Building from source (development only)

```bash
git clone https://github.com/ksugahar/ngsolve-sparsesolv.git
cd ngsolve-sparsesolv
pip install -e .  # editable install with scikit-build-core
```

Repository: https://github.com/ksugahar/ngsolve-sparsesolv
"""


SPARSESOLV_EXAMPLE_POISSON = '''# 2D Poisson Problem with ICCG
from ngsolve import *
from sparsesolv_ngsolve import SparseSolvSolver

mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
fes = H1(mesh, order=2, dirichlet="bottom|right|top|left")
u, v = fes.TnT()

a = BilinearForm(fes)
a += grad(u) * grad(v) * dx
a.Assemble()

f = LinearForm(fes)
f += 2 * pi * pi * sin(pi * x) * sin(pi * y) * v * dx
f.Assemble()

# Solve with ICCG
solver = SparseSolvSolver(a.mat, method="ICCG",
                           freedofs=fes.FreeDofs(),
                           tol=1e-10, maxiter=2000)

gfu = GridFunction(fes)
gfu.vec.data = solver * f.vec

result = solver.last_result
print(f"Converged: {result.converged}, Iterations: {result.iterations}")
print(f"Final residual: {result.final_residual:.2e}")
'''

SPARSESOLV_EXAMPLE_CURLCURL = '''# 3D Curl-Curl (Electromagnetic) with Auto-Shift ICCG
from ngsolve import *
from netgen.occ import Box, Pnt
from sparsesolv_ngsolve import SparseSolvSolver

box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
for face in box.faces:
    face.name = "outer"
mesh = Mesh(box.GenerateMesh(maxh=0.4))

fes = HCurl(mesh, order=1, dirichlet="outer", nograds=True)
u, v = fes.TnT()

a = BilinearForm(fes)
a += curl(u) * curl(v) * dx
a.Assemble()

f = LinearForm(fes)
f += CF((0, 0, 1)) * v * dx
f.Assemble()

# ICCG with auto-shift for semi-definite curl-curl system
solver = SparseSolvSolver(a.mat, method="ICCG",
                           freedofs=fes.FreeDofs(),
                           tol=1e-8, maxiter=2000, shift=1.0)
solver.auto_shift = True         # Critical for semi-definite systems
solver.diagonal_scaling = True   # Improve conditioning

gfu = GridFunction(fes)
gfu.vec.data = solver * f.vec

result = solver.last_result
print(f"Converged: {result.converged}, Iterations: {result.iterations}")
'''

SPARSESOLV_EXAMPLE_EDDY = '''# Complex Eddy Current Problem
from ngsolve import *
from netgen.occ import Box, Pnt, OCCGeometry
from sparsesolv_ngsolve import SparseSolvSolver

box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
for face in box.faces:
    face.name = "outer"
mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.3))

# Complex-valued HCurl space for eddy currents
fes = HCurl(mesh, order=1, complex=True, dirichlet="outer", nograds=True)
u, v = fes.TnT()

# Complex-symmetric: curl-curl + j*omega*sigma*mass
a = BilinearForm(fes)
a += InnerProduct(curl(u), curl(v)) * dx
a += 1j * InnerProduct(u, v) * dx
a.Assemble()

f = LinearForm(fes)
f += InnerProduct(CF((1, 0, 0)), v) * dx
f.Assemble()

# Factory auto-dispatches to complex solver
solver = SparseSolvSolver(a.mat, method="ICCG",
                           freedofs=fes.FreeDofs(),
                           tol=1e-8, maxiter=5000,
                           save_best_result=True)

gfu = GridFunction(fes)
gfu.vec.data = solver * f.vec

result = solver.last_result
print(f"Converged: {result.converged}, Iterations: {result.iterations}")
'''

SPARSESOLV_EXAMPLE_PRECOND = '''# Using IC/SGS Preconditioners with NGSolve CGSolver
from ngsolve import *
from ngsolve.krylovspace import CGSolver
from sparsesolv_ngsolve import ICPreconditioner, SGSPreconditioner

mesh = Mesh(unit_square.GenerateMesh(maxh=0.05))
fes = H1(mesh, order=2, dirichlet="bottom|right|top|left")
u, v = fes.TnT()

a = BilinearForm(fes)
a += grad(u) * grad(v) * dx
a.Assemble()

f = LinearForm(fes)
f += 1 * v * dx
f.Assemble()

gfu = GridFunction(fes)

# Option 1: IC preconditioner
pre_ic = ICPreconditioner(a.mat, freedofs=fes.FreeDofs(), shift=1.05)
pre_ic.Update()
inv = CGSolver(a.mat, pre_ic, printrates=True, tol=1e-10, maxiter=2000)
gfu.vec.data = inv * f.vec
print(f"IC+CG done")

# Option 2: SGS preconditioner
pre_sgs = SGSPreconditioner(a.mat, freedofs=fes.FreeDofs())
pre_sgs.Update()
inv = CGSolver(a.mat, pre_sgs, printrates=False, tol=1e-10)
gfu.vec.data = inv * f.vec
print(f"SGS+CG done")
'''

SPARSESOLV_EXAMPLE_DIVERGENCE = '''# Divergence Detection and Early Termination
from ngsolve import *
from sparsesolv_ngsolve import SparseSolvSolver

mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
fes = H1(mesh, order=2, dirichlet="bottom|right|top|left")
u, v = fes.TnT()

a = BilinearForm(fes)
a += grad(u) * grad(v) * dx
a.Assemble()

f = LinearForm(fes)
f += 1 * v * dx
f.Assemble()

solver = SparseSolvSolver(a.mat, method="CG",
                           freedofs=fes.FreeDofs(),
                           tol=1e-12, maxiter=5000)
solver.divergence_check = True
solver.divergence_threshold = 10.0   # residual > best * 10 = stagnation
solver.divergence_count = 20         # stop after 20 bad iterations

gfu = GridFunction(fes)
result = solver.Solve(f.vec, gfu.vec)
print(f"Converged: {result.converged}")
print(f"Iterations: {result.iterations} (max: 5000)")
print(f"Final residual: {result.final_residual:.2e}")
'''


SPARSESOLV_COMPACT_AMS = """
# Compact AMS Preconditioner (HCurl Eddy Current)

## Overview

Compact AMS (Auxiliary-space Maxwell Solver) is a **header-only, HYPRE-free**
preconditioner for HCurl curl-curl + mass systems. Based on Hiptmair-Xu (2007)
auxiliary space method, using CompactAMG as subspace solver and l1-Jacobi
as fine smoother. Fully TaskManager-native, no external library dependencies.

**When to use**: HCurl order=1 problems (any scale)
- Real magnetostatics: CompactAMSPreconditioner + CG
- Complex eddy current: ComplexCompactAMSPreconditioner + COCR

**Advantage over ICCG**: Mesh-size independent iteration count.
ICCG grows as O(h^-1), Compact AMS stays constant (~52 iters).

## Theory: Hiptmair-Xu Auxiliary Space

HCurl space has Helmholtz decomposition: u = grad(phi) + z
(phi in H1, z in curl-free complement)

Compact AMS exploits this with two auxiliary spaces:
1. **G correction** (gradient): Discrete gradient G maps H1 -> HCurl.
   Scalar CompactAMG solves the kernel component.
2. **Pi correction** (Nedelec interpolation): Pi_d = diag(coord_d) * G (d = x, y, z).
   Vector CompactAMG solves each component.
3. **Fine smoother**: l1-Jacobi (fully parallel, symmetric).

Combined as multiplicative V-cycle:
cycle_type=1 (default): smooth -> G correction -> Pi correction -> G correction -> smooth

## Why COCR (not GMRES)

Compact AMS uses **l1-Jacobi** (symmetric smoother), making the preconditioner
**symmetric**. This enables **COCR** (short recurrence, 5 vectors) instead of
GMRES (m vectors). Result: 3.5x faster than GMRES.

| Solver | Preconditioner symmetry | Memory    | Iterations (197K DOFs) |
|--------|------------------------|-----------|-----------------------|
| COCR   | Symmetric (l1-Jacobi)  | 5 vectors | 168                   |
| GMRES  | Any                    | m vectors | 384                   |

## Complex Re/Im Splitting

For complex system (K + jw*sigma*M) * x = b:
1. Build real SPD auxiliary matrix: A_real = K + |omega*sigma|*M (CRITICAL SCALING!)
2. ComplexCompactAMS creates two independent AMS instances
3. Apply: y_re = AMS(x_re), y_im = AMS(x_im) in TaskManager parallel

The `a_real` matrix must be real SPD, built in a **real** HCurl space
(same mesh, order, nograds, Dirichlet BC, but complex=False).

### CRITICAL: Real Auxiliary Matrix Scaling

**The #1 cause of AMS+COCR "not converging" is wrong scaling of a_real.**

```python
# WRONG: mass=1.0 (always diverges on real problems!)
a_real += curl(u)*curl(v)*dx + 1.0 * u*v*dx     # -> 500 iters, no convergence

# CORRECT: mass=|omega*sigma| (matches complex system scale)
a_real += curl(u)*curl(v)*dx + abs(omega*sigma) * u*v*dx  # -> 14 iters, 7.5e-10 error
```

The real auxiliary matrix MUST match the spectral structure of the complex system.
For eddy current `K + jw*sigma*M`, the correct real auxiliary is:
- `A_real = K + |omega*sigma|*M` (captures same scale)
- Optionally add `+ eps*M` (eps ~ 1e-6) for IC regularization

**Verified** (uniform sigma, 2,839 DOFs, 30 kHz, sigma=1e6):
- Wrong scaling (mass=1.0): 500 iterations, COCR does not converge
- Correct scaling (mass=|omega*sigma|): **13 iterations**, relative error 2.15e-10

**Verified** (large problem, 44,056 DOFs, |omega*sigma|=1.88e11, NGSolve 6.2.2601):
- Wrong scaling (mass=1.0): ~500 iterations, no convergence
- Correct scaling (mass=|omega*sigma|=1.88e11): **14 iterations**, relative error 7.5e-10 vs direct solver

The test confirmed that mesh-size and omega*sigma scale do NOT affect iteration count
when `a_real` is scaled correctly — only 1–2 extra iterations across a 15x DOF increase.

### Why a_real is needed

AMS auxiliary space corrections (G, Pi) operate on real arithmetic.
Complex system matrix A = K + jw*sigma*M cannot be used directly.
Real auxiliary matrix `A_real = K + |omega*sigma|*M` captures
the same spectral structure without imaginary part.

### Conductor-in-Air Problems: Mass Regularization Required

When sigma is only in a subregion (conductor in air), the curl-curl system
has a near-null space in the air region. This causes:
- **COCR breakdown** at ~78 iterations (algorithm stalls, cannot reduce residual further)
- **GMRES also stalls** (converges to a different null-space component)
- Both solutions satisfy A*x = f (the difference is in the null space of A)

**Fix**: Add mass regularization `eps * u * v * dx` to ALL domains:

```python
mu0 = 4e-7 * 3.14159265
nu = 1.0 / mu0
eps_reg = 0.05 * nu   # = 0.05/mu0, minimum for COCR convergence

# Complex system (add eps to entire domain)
a += nu * curl(u) * curl(v) * dx
a += 1j * omega * sigma_cf * u * v * dx
a += eps_reg * u * v * dx                   # regularization (all domains)

# Real auxiliary (also add eps)
a_real += nu * curl(u_r) * curl(v_r) * dx
a_real += (abs(omega*sigma) + eps_reg) * u_r * v_r * dx
```

**Impact on physical solution** (B = curl(E)):

| eps/nu | COCR iters | pardiso match | B field error |
|--------|-----------|--------------|--------------|
| 0 (none) | 78 (breakdown) | 96% diff (null space) | - |
| 0.01 | 78 (breakdown) | 0.94% | - |
| **0.05** | **117** | **1.0e-04** | **0.25%** |
| 0.10 | 111 | 9.9e-05 | ~0.3% |
| 1.00 | 121 | 9.7e-07 | 4.7% |

- `eps=0.05*nu`: Best balance (B field error 0.25%, COCR converges)
- The 96% "error" without regularization is gauge freedom (E differs, curl(E) is same)
- Uniform sigma (all conducting): No regularization needed (COCR: 13 iters, 2e-10 error)

### CompactAMGPreconditioner: pybind11 Registration (Fixed in v3.1.1)

**Symptom**: `NGException: Mult or MultAdd must be implemented for BaseMatrix!`
when using `CompactAMGPreconditioner` as a preconditioner in NGSolve's `CGSolver`.

**Root cause**: `CompactAMG` C++ class was not registered as a pybind11 type.
When the factory function returned `shared_ptr<BaseMatrix>`, pybind11 could not
find the runtime type and wrapped it as the base `ngsolve.la.BaseMatrix` Python
type. Python-level virtual dispatch then called `BaseMatrix::Mult` (which throws)
instead of `CompactAMG::Mult`.

**Fix (v3.1.1)**: Register `CompactAMG` as `CompactAMGPreconditionerImpl` pybind11
class, and change factory return type to `shared_ptr<CompactAMG>`. This pattern
(register the derived class + return derived shared_ptr) is required for any
C++ class that overrides virtual methods and is used as a Python preconditioner.

**Lesson**: All sparsesolv preconditioner/solver C++ classes that override
`BaseMatrix::Mult` must be registered with pybind11. Compare: `CompactAMS` was
registered as `CompactAMSPreconditionerImpl` and worked; `CompactAMG` was not
registered and failed silently until used in `CGSolver`.

### NGSolve built-in HCurlAMG: Known Crash

NGSolve 6.2.2601's built-in `HCurlAMG` preconditioner crashes with
access violation (-1073741819). This is an NGSolve bug, NOT a sparsesolv issue.
Use `CompactAMSPreconditioner` from sparsesolv as a working replacement.

### COCRSolver API

```python
# COCRSolver returns iterations via .iterations attribute (not .last_result)
solver = COCRSolver(a.mat, pre, maxiter=500, tol=1e-10)
gfu.vec.data = solver * f.vec
print(f"Iterations: {solver.iterations}")  # Direct attribute, not .last_result
```

## Benchmark: AMS vs ICCG Scaling

Iron-core eddy current (mu_r=1000, sigma=1e6, 30 kHz):

| DOFs   | ICCG iters | AMS+COCR iters | AMS speedup         |
|-------:|---------:|--------------:|---------------------|
| 2,728  | 97        | 46            | 1.5x fewer          |
| 6,382  | 147       | 52            | 2.8x fewer          |
| 19,357 | 234       | 52            | **4.5x fewer, 2.6x faster** |
| 197K   | 17,178    | 168           | **100x fewer**      |

AMS iteration count stays **constant** (46-52) under mesh refinement.
ICCG iteration count **grows** as O(h^-1).

## Solver Comparison

| Property         | ICCG           | Compact AMS+COCR    | BDDC (NGSolve)  |
|-----------------|----------------|---------------------|-----------------|
| Target          | General SPD    | HCurl (order=1)     | Any (real only) |
| Complex support | Yes            | **Yes**             | No              |
| Iteration scaling | h-dependent  | **h-independent**   | h-independent   |
| External deps   | None           | **None (header-only)** | NGSolve      |
| Krylov solver   | CG / COCR      | CG / **COCR**       | CG              |
| Memory          | Minimal (5 vec)| Medium (AMG levels) | Medium-Large    |
| Setup cost      | Low            | Medium (AMG setup)  | High            |
| High-order (p>1)| Yes            | **p=1 only**        | Yes             |

## Python API

```python
import sparsesolv_ngsolve as ssn

# Real Compact AMS (for real SPD HCurl magnetostatics)
pre_real = ssn.CompactAMSPreconditioner(
    mat=a.mat,                     # SparseMatrix<double>
    grad_mat=G_mat,                # Discrete gradient (fes.CreateGradient())
    freedofs=fes.FreeDofs(),       # Dirichlet mask
    coord_x=cx, coord_y=cy, coord_z=cz,  # Vertex coordinates
    cycle_type=1,                  # 1 = default multiplicative cycle
    print_level=0)

# Complex Compact AMS (TaskManager Re/Im parallel)
pre_complex = ssn.ComplexCompactAMSPreconditioner(
    a_real_mat=a_real.mat,         # Real SPD auxiliary matrix (NOT complex!)
    grad_mat=G_mat,                # From fes_real.CreateGradient()
    freedofs=fes_real.FreeDofs(),  # Real space FreeDofs
    coord_x=cx, coord_y=cy, coord_z=cz,
    ndof_complex=0,                # 0 = auto-derive from a_real_mat
    cycle_type=1,
    print_level=0)

# COCR solver (complex-symmetric, short recurrence)
solver = ssn.COCRSolver(a.mat, pre_complex, maxiter=500, tol=1e-10,
                        freedofs=fes.FreeDofs())
gfu = GridFunction(fes)
with TaskManager():
    gfu.vec.data = solver * f.vec

# H1 Compact AMG (for Poisson, heat conduction)
pre_h1 = ssn.CompactAMGPreconditioner(
    mat=a_h1.mat,
    freedofs=fes_h1.FreeDofs(),
    theta=0.25,                    # Strong connection threshold
    print_level=0)

# Newton iteration Update() (preserves geometric info, fast rebuild)
a_real.Assemble()  # Matrix values changed
pre_complex.Update(a_real.mat)  # Rebuild with new matrix
```

## COCR Solver (Complex-Symmetric)

For complex-symmetric systems (A^T = A) with symmetric preconditioner,
COCR (Sogabe-Zhang 2007) is optimal: short recurrence, minimizes ||A*r||.

```python
# Standalone COCR
solver = ssn.COCRSolver(a.mat, pre, maxiter=500, tol=1e-10)
gfu.vec.data = solver * f.vec

# Via SparseSolvSolver dispatch
solver = ssn.SparseSolvSolver(a.mat, method="COCR",
    freedofs=fes.FreeDofs(), tol=1e-10, maxiter=1000)
```

Note: COCR uses unconjugated inner product (x^T y, not x^H y).
For Hermitian systems, use CG with conjugate=True.

## Limitations

1. **order=1 only**: First-order Nedelec elements required. For high-order, use
   NGSolve BDDC (real) or direct solver.
2. **nograds=True required**: Gradient DOFs are handled by AMS G correction.
3. **TaskManager() required**: Must run inside `with TaskManager():` block.
"""

SPARSESOLV_EXAMPLE_COMPACT_AMS = '''# Compact AMS + COCR for Complex Eddy Current
# ComplexCompactAMSPreconditioner + COCRSolver (TaskManager parallel Re/Im)
#
# Two cases:
#   Case A: Uniform sigma (all conducting) -> no regularization needed
#   Case B: Conductor-in-air -> mass regularization required for COCR
from ngsolve import *
import sparsesolv_ngsolve as ssn
import math

# --- Case B: Conductor-in-Air (common real-world setup) ---
from netgen.occ import Box, Pnt, Glue
conductor = Box(Pnt(0.3, 0.3, 0.3), Pnt(0.7, 0.7, 0.7))
conductor.faces.name = "inner"
air = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
air.faces.name = "outer"
shape = Glue([air - conductor, conductor])
from netgen.occ import OCCGeometry
mesh = Mesh(OCCGeometry(shape).GenerateMesh(maxh=0.15))

freq = 30e3
omega = 2 * math.pi * freq
sigma_val = 1e6
mu0 = 4e-7 * math.pi
nu = 1.0 / mu0
sigma_cf = mesh.MaterialCF({"inner": sigma_val}, default=0)

# Mass regularization: required for conductor-in-air (COCR breakdown without it)
# eps = 0.05*nu: B field error ~0.25%, COCR converges in ~117 iterations
# eps = 0: COCR breaks down at ~78 iterations (null space in air region)
eps_reg = 0.05 * nu

# Complex HCurl space
fes = HCurl(mesh, order=1, nograds=True, dirichlet="outer", complex=True)
u, v = fes.TnT()
a = BilinearForm(fes)
a += nu * curl(u) * curl(v) * dx
a += 1j * omega * sigma_cf * u * v * dx
a += eps_reg * u * v * dx                    # regularization (all domains)
a.Assemble()

f_vec = LinearForm(fes)
f_vec += CF((0, 0, 1)) * v * dx
f_vec.Assemble()

# --- Real auxiliary matrix for Compact AMS ---
# CRITICAL: mass coefficient must be |omega*sigma| to match complex system scale!
# Using mass=1.0 will cause COCR to diverge (500+ iterations, no convergence).
fes_real = HCurl(mesh, order=1, nograds=True, dirichlet="outer", complex=False)
u_r, v_r = fes_real.TnT()
a_real = BilinearForm(fes_real)
a_real += nu * curl(u_r) * curl(v_r) * dx
a_real += (abs(omega * sigma_val) + eps_reg) * u_r * v_r * dx
a_real.Assemble()

# --- Discrete gradient + vertex coordinates ---
G_mat, h1_fes = fes_real.CreateGradient()
cx = [mesh[v].point[0] for v in mesh.vertices]
cy = [mesh[v].point[1] for v in mesh.vertices]
cz = [mesh[v].point[2] for v in mesh.vertices]

# --- ComplexCompactAMSPreconditioner ---
pre = ssn.ComplexCompactAMSPreconditioner(
    a_real.mat, G_mat, freedofs=fes_real.FreeDofs(),
    coord_x=cx, coord_y=cy, coord_z=cz,
    ndof_complex=fes.ndof)

# --- Solve with COCR (Compact AMS is symmetric -> COCR optimal) ---
solver = ssn.COCRSolver(a.mat, pre, maxiter=500, tol=1e-10,
                        freedofs=fes.FreeDofs())
gfu = GridFunction(fes)
with TaskManager():
    gfu.vec.data = solver * f_vec.vec

print(f"COCR converged in {solver.iterations} iterations")
# Expected: ~117 iterations for conductor-in-air with eps=0.05*nu
# For uniform sigma (Case A): ~13 iterations, no regularization needed
'''


def get_full_documentation() -> str:
    """Return complete ngsolve-sparsesolv documentation."""
    return "\n\n".join([
        SPARSESOLV_OVERVIEW,
        SPARSESOLV_API,
        SPARSESOLV_EXAMPLES,
        SPARSESOLV_ABMC,
        SPARSESOLV_COMPACT_AMS,
        SPARSESOLV_BEST_PRACTICES,
        SPARSESOLV_BUILD,
    ])


def get_sparsesolv_documentation(topic: str = "all") -> str:
    """Return sparsesolv documentation by topic, including code examples."""
    topics = {
        "overview": SPARSESOLV_OVERVIEW,
        "api": SPARSESOLV_API,
        "examples": SPARSESOLV_EXAMPLES,
        "abmc": SPARSESOLV_ABMC,
        "best_practices": SPARSESOLV_BEST_PRACTICES,
        "build": SPARSESOLV_BUILD,
        "example_poisson": SPARSESOLV_EXAMPLE_POISSON,
        "example_curlcurl": SPARSESOLV_EXAMPLE_CURLCURL,
        "example_eddy": SPARSESOLV_EXAMPLE_EDDY,
        "example_precond": SPARSESOLV_EXAMPLE_PRECOND,
        "example_divergence": SPARSESOLV_EXAMPLE_DIVERGENCE,
        "compact_ams": SPARSESOLV_COMPACT_AMS,
        "example_compact_ams": SPARSESOLV_EXAMPLE_COMPACT_AMS,
        # Legacy aliases
        "hypre_ams": SPARSESOLV_COMPACT_AMS,
        "example_hypre_ams": SPARSESOLV_EXAMPLE_COMPACT_AMS,
    }
    topic = topic.lower().strip()
    if topic in ("best-practices", "bestpractices", "tips"):
        topic = "best_practices"
    if topic in ("ams", "cocr"):
        topic = "compact_ams"
    if topic == "all":
        return get_full_documentation()
    elif topic in topics:
        return topics[topic]
    else:
        return (
            f"Unknown topic: '{topic}'. "
            f"Available: {', '.join(k for k in topics.keys() if not k.startswith('hypre'))}"
        )
