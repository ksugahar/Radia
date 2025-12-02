# Cube in Uniform Field Benchmark

This folder contains benchmark examples for validating Radia's magnetostatic solver using a soft iron cube in a uniform external magnetic field.

## Subfolder Structure

```
cube_uniform_field/
  linear/       - Linear material (constant mu_r)
  nonlinear/    - Nonlinear material (BH curve with saturation)
```

## Problem Description

A soft iron cube (100mm x 100mm x 100mm) is placed in a uniform external magnetic field H0 along the z-axis. This is a classic magnetostatic problem with known analytical solutions for the limiting cases.

## Benchmarks

### Linear Material (`linear/`)

Tests Radia's solver with constant permeability materials:
- Hexahedral mesh (structured)
- Tetrahedral mesh (unstructured, via Netgen)
- Method 0 (polygon-based) vs Method 1 (analytical)
- Comparison with NGSolve FEM reference solution

### Nonlinear Material (`nonlinear/`)

Tests Radia's nonlinear solver and H-matrix acceleration:
- Soft iron BH curve with saturation
- Direct solver vs H-matrix acceleration
- Convergence behavior for iterative relaxation

## Quick Start

```python
import radia as rad
rad.FldUnits('m')  # Required for NGSolve compatibility

# Linear material example
exec(open('linear/cube_benchmark_external_field.py').read())

# Nonlinear material example
exec(open('nonlinear/cube_nonlinear_hmatrix.py').read())
```

## Related Examples

- `examples/electromagnet/` - Electromagnet with iron yoke
- `examples/simple_problems/` - Basic magnetostatic problems

## References

1. Chubar, O., Elleaume, P. "Accurate and Efficient Computation of Synchrotron Radiation in the Near Field Region". EPAC98.
2. ELF/MAGIC - Electron Linear accelerator Facility, Magnetostatic Analysis using Generalized Integral Computations
