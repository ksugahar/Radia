# Cube in Uniform Field Benchmark

This folder contains benchmark examples for validating Radia's magnetostatic solver using a soft iron cube in a uniform external magnetic field.

## Subfolder Structure

```
cube_uniform_field/
  linear/       - Linear material (constant mu_r)
  nonlinear/    - Nonlinear material (BH curve with saturation)
```

## Problem Description

A soft iron cube is placed in a uniform external magnetic field H0 along the z-axis. This is a classic magnetostatic problem with known analytical solutions for the limiting cases.

---

## Linear Material (`linear/`)

Tests Radia's solver with constant permeability materials (MatLin).

### Files

| File | Description |
|------|-------------|
| `benchmark_external_field.py` | Main benchmark: hex vs tet mesh comparison |
| `benchmark_mesh_convergence.py` | Mesh refinement convergence study |
| `benchmark_solver_methods.py` | Method 0 vs Method 1 comparison |
| `benchmark_tetra_vs_hex.py` | Tetrahedral vs hexahedral accuracy |
| `benchmark_tetra_vs_ngsolve.py` | Radia vs NGSolve FEM reference |
| `benchmark_high_mu.py` | High permeability (mu_r > 1000) tests |
| `compare_external_field.py` | Field comparison at external points |
| `evaluate_perturbation_field.py` | Perturbation field analysis |
| `precision_evaluation.py` | Numerical precision evaluation |
| `test_method9_high_mu.py` | Method 9 high-mu validation |
| `test_method9_fine_mesh.py` | Method 9 fine mesh validation |
| `cube_benchmark_external_field.py` | Legacy benchmark script |
| `cube_benchmark_radia.vtk` | VTK output for visualization |
| `README.md` | Detailed documentation |
| `README_CUBE_BENCHMARK.md` | Legacy benchmark documentation |

### Key Results

- **Mesh types**: Hexahedral (ObjPolyhdr) and Tetrahedral (Netgen import)
- **Solver methods**: Method 0 (LU), Method 1 (BiCGSTAB)
- **Reference**: NGSolve FEM solution for validation
- **Accuracy**: Tetrahedral achieves <1% error vs hexahedral reference

---

## Nonlinear Material (`nonlinear/`)

Tests Radia's nonlinear solver with saturable BH curves.

### Files

| File | Description |
|------|-------------|
| `benchmark_solver_methods.py` | LU vs BiCGSTAB solver comparison |
| `benchmark_bicgstab_hex.py` | BiCGSTAB with/without H-matrix |
| `benchmark_nonlinear_tetra_vs_hex.py` | Tetra vs hex for nonlinear materials |
| `compare_radia_elfmagic_field.py` | Radia vs ELF/MAGIC comparison |
| `README.md` | Detailed documentation with benchmark results |

### Result Files (JSON)

| File | Description |
|------|-------------|
| `radia_lu_N10_results.json` | LU solver, 10x10x10 mesh |
| `radia_lu_N15_results.json` | LU solver, 15x15x15 mesh |
| `radia_lu_N20_results.json` | LU solver, 20x20x20 mesh |
| `radia_bicgstab_N10_results.json` | BiCGSTAB, 10x10x10 mesh |
| `radia_bicgstab_N15_results.json` | BiCGSTAB, 15x15x15 mesh |
| `radia_bicgstab_N20_results.json` | BiCGSTAB, 20x20x20 mesh |
| `radia_bicgstab_hmatrix_N*.json` | BiCGSTAB with H-matrix |
| `radia_bicgstab_benchmark_summary.json` | All BiCGSTAB results combined |

### Key Results (2025-12-05)

| N | Elements | LU Time | BiCGSTAB Time | Speedup |
|---|----------|---------|---------------|---------|
| 10 | 1,000 | 2.01s | 0.54s | 3.7x |
| 15 | 3,375 | 24.28s | 7.56s | 3.2x |
| 20 | 8,000 | 278.68s | 47.05s | 5.9x |

**Key Findings**:
1. BiCGSTAB is 3-6x faster than LU for nonlinear problems
2. H-matrix shows NO speedup for compact cube geometry
3. Both methods produce identical results within tolerance

---

## Quick Start

```python
import radia as rad
rad.FldUnits('m')  # Required for NGSolve compatibility

# Linear material benchmark
cd examples/cube_uniform_field/linear
python benchmark_external_field.py

# Nonlinear material benchmark
cd examples/cube_uniform_field/nonlinear
python benchmark_solver_methods.py 10 --all
```

## Related Examples

- `examples/electromagnet/` - Electromagnet with iron yoke
- `examples/simple_problems/` - Basic magnetostatic problems
- `examples/solver_benchmarks/` - H-matrix and solver performance tests

## References

1. Chubar, O., Elleaume, P. "Accurate and Efficient Computation of Synchrotron Radiation in the Near Field Region". EPAC98.
2. ELF/MAGIC - Magnetostatic Analysis using Generalized Integral Computations
3. HACApK - Lattice H-matrix Approximation on Parallel Kernel (reference implementation)

---

**Last Updated**: 2025-12-05
