# Field Computation Benchmark

Field computation performance benchmark for Radia, designed to match ELF_MAGIC's fmm_field_benchmark.

## Overview

Benchmarks magnetic field computation at multiple evaluation points around a magnetized cube.

**Current Status**: Direct computation only. FMM acceleration planned for future implementation.

### Computation Complexity

| Method | Complexity | Description |
|--------|------------|-------------|
| Direct | O(N_elements x N_points) | Current Radia implementation |
| FMM | O(N_elements + N_points) | Planned future implementation |

## Usage

```bash
cd examples/fmm_field_benchmark

# Hexahedral mesh
python benchmark_fmm_field.py --hex 10

# Tetrahedral mesh (requires Netgen)
python benchmark_fmm_field.py --tetra 0.2 0.15

# Custom grid size
python benchmark_fmm_field.py --hex 10 --n_grid 30
```

## Command Line Options

| Option | Default | Description |
|--------|---------|-------------|
| `--hex` | 5, 10 | Hexahedral mesh divisions (N^3 elements) |
| `--tetra` | - | Tetrahedral mesh maxh (requires Netgen) |
| `--n_grid` | 20 | Field evaluation grid size per axis |
| `--output` | fmm_benchmark_results.json | Output JSON file |

## Benchmark Conditions

- **Cube size**: 1.0 m
- **Material**: Linear, mu_r = 1000
- **External field**: H_ext = 200,000 A/m (z-direction)
- **Solver**: BiCGSTAB (Method 1)
- **Grid**: 3x3x3 m^3 centered at origin
- **Exclusion zone**: 0.5 m margin from cube surface

## Expected Results

| Mesh | Elements | Eval Points | Time | Points/sec |
|------|----------|-------------|------|------------|
| Hex N=5 | 125 | ~6,000 | TBD | TBD |
| Hex N=10 | 1,000 | ~6,000 | TBD | TBD |
| Hex N=15 | 3,375 | ~6,000 | TBD | TBD |

## Output Format

Results are saved to JSON with the following fields:

```json
{
  "mesh_type": "hex",
  "mesh_param": 10,
  "mesh_desc": "Hex N=10",
  "n_elements": 1000,
  "ndof": 6000,
  "n_grid": 20,
  "n_points": 6272,
  "t_mesh": 0.5,
  "t_solve": 2.0,
  "t_direct": 5.0,
  "direct_points_per_sec": 1254,
  "B_max": 0.35,
  "B_min": 0.18,
  "H_max": 280000,
  "H_min": 145000
}
```

## Comparison with ELF_MAGIC

This benchmark is designed to match ELF_MAGIC's `examples/fmm_field_benchmark/` for direct performance comparison.

| Feature | ELF_MAGIC | Radia |
|---------|-----------|-------|
| Direct computation | Yes | Yes |
| FMM acceleration | Yes (20x speedup) | Planned |
| Inside/outside detection | Yes | Planned |
| Scalar potential phi | Yes | Planned |
| Vector potential A | Yes | Planned |

## Implementation Plan

1. [x] Direct field computation benchmark
2. [ ] Point classification (inside/near/far)
3. [ ] FMM3D library integration
4. [ ] Hybrid FMM + near-field correction
5. [ ] Scalar potential computation
6. [ ] Vector potential computation

## Files

| File | Description |
|------|-------------|
| `benchmark_fmm_field.py` | Main benchmark script |
| `README.md` | This documentation |
| `fmm_benchmark_results.json` | Benchmark output |
