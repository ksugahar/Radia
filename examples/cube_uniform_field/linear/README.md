# Linear Material Benchmark: Cube in Uniform External Field

This benchmark validates Radia's tetrahedral and hexahedral mesh solvers for **linear magnetic materials**.

## Problem Setup

- **Geometry**: 100mm x 100mm x 100mm cube centered at origin
- **Material**: Linear isotropic, mu_r = 100 (chi = 99)
- **External Field**: B0 = 1.0 T (uniform along z-axis)
- **Validation Metric**: External field at points OUTSIDE the cube

**Important**: Radia MMM is designed for external field calculation. The proper validation metric is the field in air regions, not inside the magnetic material.

## Benchmark Results (2025-12-02)

### External Field Comparison

Reference: Hexahedral 8x8x8 mesh (512 elements)

| Mesh Type | Elements | Avg Error | Max Error | Notes |
|-----------|----------|-----------|-----------|-------|
| **Hexahedral Meshes** |
| Hex 1x1x1 | 1 | 3.33% | 9.27% | Single element |
| Hex 2x2x2 | 8 | 2.38% | 6.80% | |
| Hex 3x3x3 | 27 | 1.17% | 3.53% | |
| Hex 4x4x4 | 64 | 0.65% | 1.98% | |
| Hex 5x5x5 | 125 | 0.37% | 1.12% | |
| Hex 6x6x6 | 216 | 0.19% | 0.59% | Good convergence |
| **Tetrahedral Meshes (Manual 5-tetra cells)** |
| Tet 1^3x5 | 5 | 3.03% | 8.18% | Single cell |
| Tet 2^3x5 | 40 | 0.43% | 0.87% | Good accuracy |
| Tet 3^3x5 | 135 | 0.23% | 0.62% | Excellent |
| **Tetrahedral Meshes (Netgen/NGSolve)** |
| maxh=0.10 | 12 | 1.93% | 5.34% | |
| maxh=0.06 | 28 | 0.20% | 0.49% | Best accuracy |
| maxh=0.04 | 80 | 0.26% | 0.67% | |

### Key Findings

1. **Hexahedral mesh convergence**: Error decreases consistently with mesh refinement
2. **Tetrahedral mesh accuracy**: Comparable to hexahedral at similar element counts
3. **External field validation**: Proper metric for MMM validation

## Mesh Configurations

### Hexahedral Mesh (ObjRecMag + ObjDivMag)

Standard Radia approach using built-in mesh subdivision:

```python
import radia as rad
rad.FldUnits('m')

cube = rad.ObjRecMag([0, 0, 0], [0.1, 0.1, 0.1], [0, 0, 0])
rad.ObjDivMag(cube, [5, 5, 5])  # 5x5x5 = 125 elements
mat = rad.MatLin(99.0)  # chi = mu_r - 1
rad.MatApl(cube, mat)
```

### Tetrahedral Mesh (Manual 5-tetra decomposition)

Each cubic cell is decomposed into 5 tetrahedra:

```python
rad.SolverTetraMethod(0)  # Method 0: polygon-based
# Create 5 tetrahedra per cubic cell using rad.ObjPolyhdr()
```

### Tetrahedral Mesh (Netgen import)

Using netgen_mesh_import for unstructured meshes:

```python
from netgen_mesh_import import netgen_mesh_to_radia
from netgen.occ import Box, OCCGeometry
from ngsolve import Mesh

cube_solid = Box((-0.05, -0.05, -0.05), (0.05, 0.05, 0.05))
geo = OCCGeometry(cube_solid)
mesh = Mesh(geo.GenerateMesh(maxh=0.06))

mag_obj = netgen_mesh_to_radia(mesh, material={'magnetization': [0, 0, 0]}, units='m')
```

## Files

| File | Description |
|------|-------------|
| `benchmark_external_field.py` | Main benchmark comparing hex/tet meshes |
| `compare_external_field.py` | External field comparison with NGSolve |
| `precision_evaluation.py` | Comprehensive precision evaluation |
| `cube_benchmark_ngsolve.vtu` | NGSolve solution (VTU format) |
| `cube_benchmark_radia.vtk` | Radia solution (VTK format) |

## Running the Benchmark

```bash
cd examples/cube_uniform_field/linear
python benchmark_external_field.py
```

Output will show:
- Hexahedral mesh convergence study (1x1x1 to 6x6x6)
- Tetrahedral mesh (manual 5-tetra decomposition)
- Tetrahedral mesh (Netgen import)
- Summary table with error statistics

## Physics Background

For a linear magnetic material in uniform external field:

- **Demagnetizing factor**: N ~ 1/3 for cube (approximate)
- **Interior H-field**: H_int = H0 / (1 + N*(mu_r - 1))
- **Magnetization**: M = chi * H_int
- **Exterior field**: Dipole field from magnetized cube

## High Permeability Behavior (mu_r >= 100)

### Convergence Characteristics

| mu_r | Hex 4x4x4 | Hex 8x8x8 | Tet 5 | Notes |
|------|-----------|-----------|-------|-------|
| 100 | 404 iter | ~20k iter | 22 iter | Fine hex mesh slow to converge |
| 500 | 2095 iter | N/C | 22 iter | |
| 1000 | 3758 iter | N/C | 22 iter | |
| 2000 | 844 iter | N/C | 20 iter | |
| 4000 | 37k iter | N/C | 21 iter | Coarse mesh required |

N/C = Not converged within 20,000 iterations

### Key Observations

1. **Coarse mesh converges faster**: Hex 4x4x4 converges while 8x8x8 may not
2. **Tet single cell**: 5-element tetra always converges quickly (~20 iter)
3. **Accuracy trade-off**: Coarser meshes converge but have lower spatial resolution

### Recommendations for High mu_r

1. Use coarser meshes (4x4x4 to 6x6x6 hex subdivisions)
2. Increase iteration limit (50,000+) for fine meshes
3. Use relaxed tolerance (0.001 instead of 0.0001)
4. For tetrahedral: single 5-tetra cell provides fast convergence

## Known Limitations

1. **High permeability (mu_r > 100)**:
   - Fine meshes may not converge within reasonable iteration limits
   - Use coarser meshes for better convergence

2. **Adjacent tetrahedra**:
   - Dense tetrahedral meshes with many shared faces may show numerical instability
   - Single 5-tetra cell decomposition is most stable

3. **Internal field**:
   - `rad.Fld()` inside magnetic materials has limited accuracy
   - Always validate using external field points

## Comparison with ELF/MAGIC Solver

### Why ELF/MAGIC is Stable for High Permeability

ELF/MAGIC uses a **direct LU decomposition solver** (LAPACK DGESV) instead of Radia's **iterative relaxation method**. This fundamental difference explains why ELF is stable for high permeability materials.

### Solver Approaches

| Aspect | Radia | ELF/MAGIC |
|--------|-------|-----------|
| **Method** | Iterative (successive substitution) | Direct (LU decomposition) |
| **Implementation** | Gauss-Seidel relaxation | LAPACK DGESV |
| **High mu_r** | Convergence issues | Always stable |
| **Memory** | O(N^2) for matrix | O(N^2) for matrix |
| **Time per solve** | O(iter * N^2) | O(N^3) |

### Key Differences

1. **Radia's Iterative Method** (`rad_relaxation_methods.cpp`):
   - Uses successive substitution: `M_new = omega * M(H) + (1-omega) * M_old`
   - Spectral radius of iteration matrix increases with mu_r
   - For high permeability, |chi * N_diag| >> 1, causing slow/no convergence
   - Multiple relaxation methods (0-4) but all iterative

2. **ELF/MAGIC's Direct Method** (`m_dll_solver.f90`):
   - Assembles full coefficient matrix: `[I + chi * N] * M = chi * H_ext`
   - Solves directly using LAPACK DGESV (LU with partial pivoting)
   - Stable regardless of permeability
   - Also supports iterative (Gauss-Seidel) and H-matrix (BiCGSTAB) modes

### Matrix Equation

The MMM equation for linear materials:
```
[I + chi * N] * M = chi * H_external
```

Where:
- `I` = identity matrix
- `chi` = susceptibility (mu_r - 1)
- `N` = demagnetization tensor (from element interactions)
- `M` = magnetization vector (unknown)
- `H_external` = applied field

For high chi, the matrix `[I + chi * N]` becomes increasingly ill-conditioned for iterative methods, but LU decomposition handles it reliably.

### ELF/MAGIC Code Reference

From `S:\ELF_MAGIC\01_Github\src\dll\m_dll_solver.f90`:
```fortran
! Solve linear system using LU decomposition (LAPACK DGESV)
subroutine solve_linear_system(A, b, x, n, ierr)
    real(wp), intent(inout) :: A(:,:)  ! Will be overwritten by LU factors
    real(wp), intent(in) :: b(:)
    real(wp), intent(out) :: x(:)
    integer, intent(in) :: n
    integer, intent(out) :: ierr

    integer :: ipiv(n)
    real(wp) :: work_b(n)

    work_b = b(1:n)
    call dgesv(n, 1, A, n, ipiv, work_b, n, ierr)

    if (ierr == 0) then
        x(1:n) = work_b
    end if
end subroutine
```

### Radia Method 9: Direct LU Solver (Implemented!)

Radia now includes **Method 9**, a direct LU solver equivalent to ELF/MAGIC's approach:

```python
import radia as rad

# Use Method 9 (direct LU solver) for high permeability
result = rad.Solve(system, 0.001, 1, 9)  # precision, max_iter, method=9
```

### Method 9 Benchmark Results (2025-12-02)

| mu_r | Mesh | Method 4 (iter) | Method 9 (LU) | Status |
|------|------|-----------------|---------------|--------|
| 100 | 4x4x4 | 404 iter, 0.01s | 1 iter, 0.002s | Both OK |
| 500 | 4x4x4 | 2095 iter, 0.03s | 1 iter, 0.002s | Both OK |
| 1000 | 4x4x4 | 3758 iter, 0.06s | 1 iter, 0.002s | Both OK |
| 4000 | 4x4x4 | 37160 iter, 0.5s | 1 iter, 0.002s | Both OK |
| 100 | 8x8x8 | N/C (20000 iter) | 1 iter, 0.7s | M9 OK |
| 500 | 8x8x8 | N/C (20000 iter) | 1 iter, 0.7s | M9 OK |
| 1000 | 8x8x8 | N/C (20000 iter) | 1 iter, 1.0s | M9 OK |

N/C = Not converged within iteration limit

### When to Use Method 9

**Use Method 9** when:
- High permeability materials (mu_r > 100)
- Fine meshes where iterative methods don't converge
- Small to medium problem sizes (< 1000 elements)
- Linear materials (`rad.MatLin()`)

**Use Method 4** (default) when:
- Large problem sizes (> 1000 elements) - iterative is faster
- Nonlinear materials (B-H curves)
- Low permeability materials

### Time Complexity

| Method | Time Complexity | Memory |
|--------|-----------------|--------|
| Method 4 (iterative) | O(iter * N^2) | O(N^2) |
| Method 9 (direct LU) | O(N^3) | O(N^2) |

For N=512 (8x8x8 mesh), Method 9 takes ~0.7s regardless of mu_r.
For N=64 (4x4x4 mesh), Method 9 takes ~2ms.

## Notes

1. **Unit System**: All examples use meters (`rad.FldUnits('m')`)
2. **Material Definition**: Use isotropic `rad.MatLin(chi)` for best results
3. **External field**: Evaluate at points outside the magnetic material for validation
4. **Solver settings**: For high mu_r, use `rad.Solve(obj, 0.001, 1, 9)` (Method 9)
