# Solver Benchmarks

Comprehensive benchmarks comparing Radia's solver methods and H-matrix acceleration.

## Benchmarks

### 1. Solver Methods Comparison (`benchmark_solver_methods.py`)

**Problem:**
- Cube subdivided into N×N×N elements (Steel37 nonlinear material)
- Uniform background field: [10000, 0, 0] A/m
- Initial magnetization: [0, 0, 1000] A/m
- Test sizes: N=3,5,7 (27, 125, 343 elements)

**Methods Compared:**
1. **Direct calculation** - Uses initial magnetization only (no solver iterations)
2. **Standard relaxation** - Iterative solver without H-matrix (reference method)
3. **H-matrix relaxation** - Accelerated relaxation using hierarchical matrices

**Metrics:**
- Computation time
- Field accuracy (relative error vs standard relaxation)
- Memory usage (for H-matrix)
- Speedup factor

**Purpose:**
Demonstrate the trade-offs between solver methods:
- Direct: Fastest but lowest accuracy (no material interaction)
- Relaxation: Most accurate but slowest for large problems
- H-matrix: Good accuracy with significant speedup for N > 100

### 2. H-matrix Field Evaluation (`benchmark_hmatrix_field.py`)

**Problem:**
- Grid: 10×10 = 100 rectangular magnets
- Observation points: 10×10 = 100 points
- Initial magnetization: [0, 0, 1] T

**Comparison:**
1. Direct field calculation (use_hmatrix=0)
2. H-matrix field calculation (use_hmatrix=1)

**Purpose:**
Demonstrate H-matrix accuracy and performance for field evaluation.
Target accuracy: < 1% relative error.

## Running

```bash
# Solver methods comparison (recommended)
python benchmark_solver_methods.py

# H-matrix field evaluation
python benchmark_hmatrix_field.py
```

## Expected Results

### Solver Methods Comparison

For N=27 (small problem):
- Direct: < 0.001s, high error (~50-100%)
- Relaxation: ~0.005s, reference (0% error)
- H-matrix: ~0.010s, low error (< 1%), 0.5x speedup (overhead dominates)

For N=343 (large problem):
- Direct: < 0.001s, high error (~50-100%)
- Relaxation: ~0.1s, reference (0% error)
- H-matrix: ~0.03s, low error (< 1%), **3-4x speedup**

### H-matrix Field Evaluation

- Accuracy: < 1% relative error (typically 0.01-0.1%)
- Memory: ~0.2-0.3 MB for N=100
- Speedup varies with problem size

## Summary

| Benchmark | Problem | Methods | Best Method |
|-----------|---------|---------|-------------|
| Solver Methods | Nonlinear material + background field | Direct / Relax / H-matrix | H-matrix (large N) |
| Field Evaluation | Permanent magnets | Direct / H-matrix | H-matrix (N > 100) |

## Key Findings

1. **Direct Method**
   - Fastest computation
   - Does not account for material interactions
   - Suitable only for permanent magnets (no solver needed)

2. **Standard Relaxation**
   - Most accurate (reference method)
   - Computational cost: O(N²) per iteration
   - Best for small-medium problems (N < 100)

3. **H-matrix Relaxation**
   - Slightly lower accuracy (< 1% error)
   - Computational cost: ~O(N log N) per iteration
   - **Recommended for large problems (N > 100)**
   - Speedup increases with problem size

## References

- CHANGELOG.md: v1.0.10 H-matrix implementation details
- CLAUDE.md: H-matrix usage guidelines
- tests/test_relaxation.py: Relaxation solver tests
