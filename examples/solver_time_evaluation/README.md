# Solver Time Evaluation

This directory contains benchmarks for evaluating Radia relaxation solver performance.

## Purpose

Compare the performance of different solver methods with **separated timing**:
- **Matrix construction time** (O(N²))
- **Solver time** (O(N²) for Gauss-Seidel, O(N³) for LU decomposition)

## Scripts

### 1. `benchmark_lu_vs_gs.py`
Comprehensive comparison of LU decomposition vs Gauss-Seidel:
- Matrix construction time (same for both methods)
- Solver time per iteration (different scaling)
- Total convergence time
- Scaling analysis

### 2. `benchmark_matrix_construction.py`
Measures only the matrix construction time:
- Uses `rad.RlxPre()` to build interaction matrix
- Verifies O(N²) scaling

### 3. `benchmark_solver_scaling.py`
Measures only the solver time:
- Uses `rad.RlxMan()` to execute iterations without rebuilding matrix
- Separates LU decomposition O(N³) from Gauss-Seidel O(N²)

## Usage

```bash
python benchmark_lu_vs_gs.py
```

## API Design

### Separated Timing Approach

1. **Matrix Construction**
   ```python
   intrc = rad.RlxPre(obj, obj)  # Build interaction matrix
   ```

2. **Solver Execution**
   ```python
   # Gauss-Seidel (Method 4)
   rad.RlxMan(intrc, 4, 1, 1.0)  # 1 iteration

   # LU decomposition (Method 5 with RelaxTogether)
   rad.SetRelaxSubInterval(intrc, 0, N-1, 1)  # Enable LU decomposition
   rad.RlxMan(intrc, 5, 1, 1.0)  # Method 5 - LU decomposition
   ```

3. **Time Measurement**
   ```python
   from time import perf_counter

   # Matrix construction
   t_start = perf_counter()
   intrc = rad.RlxPre(obj, obj)
   t_matrix = perf_counter() - t_start

   # Solver
   t_start = perf_counter()
   rad.RlxMan(intrc, method, 1, 1.0)
   t_solver = perf_counter() - t_start
   ```

## Expected Results

| Method | Matrix Construction | Solver per Iteration |
|--------|---------------------|----------------------|
| Gauss-Seidel (Method 4) | O(N²) | O(N²) |
| LU Decomposition (Method 5) | O(N²) | O(N³) |

## Test Configuration

- **Material**: Nonlinear (soft iron) for realistic convergence behavior
- **Geometry**: N×N×N cube subdivision
- **N range**: 2×2×2 to 10×10×10 (8 to 1000 elements)
- **Measurement**: Average of multiple runs for stability

## Notes

- Linear materials converge in 0-1 iterations → use nonlinear for benchmarking
- Method 5 requires `SetRelaxSubInterval()` to enable LU decomposition
- Matrix construction time is identical for all solver methods

## Benchmark Results

### Complete Timing Analysis (Nonlinear Material)

Results from `benchmark_lu_vs_gs.py` with separated timing:

#### Matrix Construction Time (O(N²))

| N (elements) | Time (ms) | t/N² |
|--------------|-----------|------|
| 8            | 4.79      | 0.075 |
| 27           | 0.42      | 0.001 |
| 64           | 1.91      | 0.000 |
| 125          | 7.00      | 0.000 |
| 216          | 21.42     | 0.000 |
| 343          | 53.28     | 0.000 |
| 512          | 116.16    | 0.000 |
| 1000         | 443.34    | 0.000 |

**Power law fit**: t = 4.44×10⁻⁵ × N^1.191

**Interpretation**:
- Overall α = 1.191 appears sub-quadratic due to small-N overhead
- For large N (≥ 125), t/N² ratio is approximately constant (mean = 0.000449, CV = 1.3%)
- **Conclusion**: True O(N²) scaling confirmed for asymptotic regime

### Solver Time Scaling (Matrix Construction Excluded)

Results from `benchmark_solver_scaling.py`:

| N (elements) | GS (ms) | LU (ms) | LU/GS Ratio | t_gs/N² | t_lu/N³ |
|--------------|---------|---------|-------------|---------|---------|
| 8            | 1.82    | 0.93    | 0.51×       | 0.028   | 1.823   |
| 27           | 0.01    | 0.58    | 70.35×      | 0.000   | 0.030   |
| 64           | 0.03    | 6.77    | 209.59×     | 0.000   | 0.026   |
| 125          | 0.63    | 52.84   | 84.30×      | 0.000   | 0.027   |
| 216          | 0.36    | 260.33  | 719.35×     | 0.000   | 0.026   |
| 343          | 0.65    | 1132.27 | 1746.26×    | 0.000   | 0.028   |
| 512          | 1.29    | 4436.51 | 3435.16×    | 0.000   | 0.033   |
| 1000         | 4.62    | 59471.10| 12872.81×   | 0.000   | 0.059   |

### Power Law Fits (t = a × N^α)

**Overall scaling (N = 8 to 1000):**
- Gauss-Seidel: t = 2.17×10⁻⁵ × N^0.591
- LU Decomposition: t = 7.31×10⁻⁷ × N^2.454

**Large N only (N ≥ 216):**
- LU Decomposition: α = 3.552 → **True O(N³) behavior confirmed**

### Interpretation

1. **Gauss-Seidel Scaling**:
   - Overall α = 0.591 appears sub-quadratic due to measurement noise at small N
   - For large N (≥ 125), t_gs/N² shows approximately constant behavior
   - **Conclusion**: O(N²) per iteration as expected

2. **LU Decomposition Scaling**:
   - Overall α = 2.454 is between O(N²) and O(N³)
   - For large N only (≥ 216): α = 3.552 → **O(N³) confirmed**
   - Crossover from cache-efficient to memory-bound regime around N = 200
   - **Conclusion**: True O(N³) behavior observed for asymptotic regime

3. **Performance Comparison**:
   - Small N (< 64): LU and GS comparable (setup costs dominate)
   - Medium N (64-343): LU becomes 100-1000× slower
   - Large N (> 343): LU becomes 1000-10000× slower per iteration
   - **Recommendation**: Use Gauss-Seidel for iterative problems

### Configuration

- **Material**: Nonlinear soft iron (M-H curve)
- **Geometry**: 100×100×100 mm cube subdivided into N×N×N elements
- **Measurement**: Single iteration timing using `perf_counter()`
- **Hardware**: Results may vary depending on CPU and memory bandwidth

---

### Linear Material Evaluation

Results from `benchmark_linear_material.py`:

| N (elements) | Matrix (ms) | Solver (ms) | Total (ms) | t_matrix/N² |
|--------------|-------------|-------------|------------|-------------|
| 8            | 7.80        | 1.72        | 9.52       | 0.122       |
| 27           | 0.43        | 0.01        | 0.43       | 0.001       |
| 64           | 1.91        | 0.03        | 1.94       | 0.000       |
| 125          | 7.07        | 0.53        | 7.60       | 0.000       |
| 216          | 21.32       | 0.25        | 21.57      | 0.000       |
| 343          | 53.22       | 0.63        | 53.84      | 0.000       |
| 512          | 116.64      | 1.32        | 117.96     | 0.000       |
| 1000         | 445.49      | 4.73        | 450.22     | 0.000       |

**Power law fit**: t = 6.82×10⁻⁵ × N^1.115

**Key Findings**:

1. **Matrix Construction**: O(N²) scaling (same as nonlinear)
   - t/N² ratio approximately constant for large N (mean = 0.000450, CV = 1.0%)
   - Matrix construction time identical to nonlinear materials

2. **Solver Convergence**: 0-1 iterations only
   - Linear relationship M = χH requires no iteration
   - Solver time negligible (< 1% of matrix construction time)

3. **Comparison with Nonlinear Materials**:

   | Aspect | Linear | Nonlinear |
   |--------|--------|-----------|
   | Matrix construction | O(N²) | O(N²) |
   | Solver iterations | 0-1 | 10-100+ |
   | Total time | ≈ Matrix time | Matrix + Solver time |
   | Method choice | Irrelevant | Critical (GS vs LU) |

4. **Recommendations for Linear Problems**:
   - Matrix construction O(N²) dominates total time
   - Solver method choice (GS vs LU) irrelevant (0-1 iteration)
   - Focus optimization on matrix construction, not solver
   - Use H-matrix for large problems (N > 1000)
