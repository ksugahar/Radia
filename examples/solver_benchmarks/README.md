# Solver Benchmarks

H-matrix field evaluation benchmark.

## Benchmark

### H-matrix Field Evaluation (`benchmark_hmatrix_field.py`)

**Problem:**
- Grid: 10×10 = 100 rectangular magnets
- Observation points: 10×10 = 100 points
- Initial magnetization: [0, 0, 1] T

**Comparison:**
1. Direct calculation (use_hmatrix=0)
2. H-matrix calculation (use_hmatrix=1)

**Purpose:**
Demonstrate H-matrix accuracy and performance for field evaluation.

## Running

```bash
python benchmark_hmatrix_field.py
```

## Expected Results

- **Accuracy:** < 1% relative error
- **Memory:** ~0.2-0.3 MB for N=100
- **Speedup:** Variable (H-matrix overhead may dominate for small problems)

## References

- CLAUDE.md: H-matrix testing guidelines
- test_hmatrix_large.py: Original test (same as benchmark_hmatrix_field.py)
