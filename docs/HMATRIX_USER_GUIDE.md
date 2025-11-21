# H-Matrix Acceleration User Guide

**Version:** 1.0.10
**Date:** 2025-11-13

## Overview

This guide explains how to use Radia's H-matrix acceleration feature for faster magnetostatic simulations. H-matrix provides up to **350x speedup** for large problems and iterative workflows.

## Table of Contents

1. [Quick Start](#quick-start)
2. [When to Use H-Matrix](#when-to-use-h-matrix)
3. [Automatic vs Manual Configuration](#automatic-vs-manual-configuration)
4. [Understanding Cache Behavior](#understanding-cache-behavior)
5. [Performance Optimization Tips](#performance-optimization-tips)
6. [Troubleshooting](#troubleshooting)
7. [Advanced Topics](#advanced-topics)

---

## Quick Start

### Basic Usage (Manual Enable)

**Important**: H-matrix acceleration must be **explicitly enabled** by the user. It is NOT automatic.

```python
import radia as rad

# Create geometry (e.g., 7×7×7 = 343 elements)
mat = rad.MatSatIsoFrm([1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759])
cube = rad.ObjRecMag([0, 0, 0], [100, 100, 100], [0, 0, 0])
rad.MatApl(cube, mat)
rad.ObjDivMag(cube, [7, 7, 7])  # 343 elements

# Enable H-matrix explicitly (required!)
rad.SolverHMatrixEnable(1, 1e-4, 30)  # enable=1, eps=1e-4, max_rank=30

# Solve with H-matrix
rad.Solve(cube, 0.0001, 1000)
```

**Console Output:**
```
[Phase 2-B] Using H-matrix solver (N=343)
[Phase 3] Cache miss (new geometry, hash=e20cad57)
Building True H-Matrix...
Construction time: 1.26 s
[Phase 3] Saved to cache (./.radia_cache/hmatrix_cache.bin)
```

### Disabling H-Matrix

```python
# Disable H-matrix (use dense solver)
rad.SolverHMatrixEnable(0)

# Solve with dense solver
rad.Solve(cube, 0.0001, 1000)
```

---

## When to Use H-Matrix

### ⚠️ Important: Manual Enable Policy

**Design Decision**: H-matrix is **intentionally NOT automatic**. Users must explicitly enable it.

**Reasons**:
1. **User control**: Different applications have different performance requirements
2. **Small models**: H-matrix construction overhead can exceed benefits for N < 200
3. **Explicit behavior**: No hidden automatic switching based on problem size
4. **Predictable performance**: Users decide when to use H-matrix

**Recommendation**: Enable H-matrix for **N ≥ 200** for best performance.

### ✅ Recommended For

| Problem Type | Element Count | Expected Speedup | Notes |
|--------------|---------------|------------------|-------|
| **Large single solve** | N ≥ 200 | 10-100x | Enable explicitly |
| **Iterative workflows** | N ≥ 100 | 25-350x | H-matrix reuse |
| **Magnetization optimization** | Any N | 33x | Only M changes |
| **Geometry exploration** | N ≥ 200 | 7-10x | Multiple geometries |

### ❌ Not Recommended For

| Problem Type | Reason | Alternative |
|--------------|--------|-------------|
| **Small problems** (N < 200) | H-matrix overhead > benefit | Use dense solver (disable H-matrix) |
| **Single small solve** (N < 100) | Construction time dominates | Use dense solver |
| **Memory-constrained systems** | H-matrix uses more memory | Use dense solver |

**Important**: For small models (N < 200), enabling H-matrix may **NOT improve performance** and can even be slower due to construction overhead. Always benchmark your specific use case.

### Performance Expectations

```
Small problems (N=125):    May be slower (construction overhead)
Large single solve (N=343): 33x faster
100 solves (magnetization): 25x faster
10 geometries × 10 solves:  7.9x faster
```

---

## H-Matrix Enabled APIs

H-matrix acceleration affects the following Radia operations when explicitly enabled via `rad.SolverHMatrixEnable(1)`:

### 1. `rad.Solve()` - Relaxation Solver

**Affected**: Yes (primary use case)

```python
rad.SolverHMatrixEnable(1, 1e-4, 30)  # Enable H-matrix
rad.Solve(container, 0.0001, 1000)    # Uses H-matrix for interaction matrix
```

**Performance Impact**: 10-350x speedup for N ≥ 200

### 2. `rad.Fld()` - Field Calculation

**Affected**: No (H-matrix does NOT accelerate field evaluation)

```python
rad.SolverHMatrixEnable(1, 1e-4, 30)  # Enable H-matrix
rad.Solve(container, 0.0001, 1000)    # H-matrix used here

# Field evaluation - H-matrix NOT used
B = rad.Fld(container, 'b', [0, 0, 10])  # Direct calculation, no H-matrix
```

**Note**: Field evaluation always uses direct calculation, regardless of H-matrix setting.

### 3. `radia_ngsolve.RadiaField()` - NGSolve Integration

**Affected**: No (H-matrix does NOT accelerate field evaluation on NGSolve mesh)

```python
rad.SolverHMatrixEnable(1, 1e-4, 30)  # Enable H-matrix
rad.Solve(container, 0.0001, 1000)    # H-matrix used for solving

# NGSolve field evaluation - H-matrix NOT used
B_cf = radia_ngsolve.RadiaField(container, 'b')  # Calls rad.Fld() per point
gf.Set(B_cf)  # No H-matrix acceleration
```

**Note**: Each mesh point evaluation calls `rad.Fld()` which does not use H-matrix.

### Summary Table

| API | H-Matrix Used? | Performance Impact | Notes |
|-----|----------------|-------------------|-------|
| `rad.Solve()` | ✓ Yes | 10-350x speedup | Interaction matrix assembly |
| `rad.Fld()` | ✗ No | No speedup | Direct field calculation |
| `radia_ngsolve.RadiaField()` | ✗ No | No speedup | Calls `rad.Fld()` internally |

---

## Manual Configuration

Users must explicitly configure H-matrix:

```python
# Enable H-matrix with custom parameters
rad.SolverHMatrixEnable(1, 1e-4, 30)  # enable=1, eps=1e-4, max_rank=30

# Use tighter tolerance for critical applications
rad.SolverHMatrixEnable(1, 1e-6, 50)

# Use relaxed tolerance for faster construction
rad.SolverHMatrixEnable(1, 5e-4, 20)

# Disable H-matrix (use dense solver)
rad.SolverHMatrixEnable(0)
```

**When to adjust parameters:**
- Debugging: Compare H-matrix vs dense results (`enable=0` vs `enable=1`)
- Performance tuning: Test different `eps` and `max_rank` values
- Accuracy requirements: Trade-off between speed and precision

---

## Understanding Cache Behavior

### Geometry Cache (In-Memory)

H-matrix is cached in memory during the session. Reused when geometry is unchanged.

**Example: Magnetization-Only Changes**
```python
# First solve - builds H-matrix
rad.Solve(container, 0.0001, 1000)  # ~1000 ms

# Change magnetization only
for element in elements:
    rad.ObjM(element, [1, 0, 0])  # New magnetization

# Second solve - reuses H-matrix
rad.Solve(container, 0.0001, 1000)  # ~30 ms (33x faster!)
```

**Console Output:**
```
[Phase 2-B] Reusing H-matrix (geometry unchanged, hash=e20cad57)
```

### Geometry Change Detection

Geometry hash invalidates cache when:
- Element count changes
- Element positions change
- Geometry topology changes

**Does NOT invalidate when:**
- Magnetization changes
- Material properties change
- External field changes

**Example: Geometry Change**
```python
# First solve
rad.Solve(container, 0.0001, 1000)

# Add new element (geometry changed)
rad.ObjAddToCnt(container, [new_element])

# Second solve - rebuilds H-matrix
rad.Solve(container, 0.0001, 1000)  # Full construction time
```

**Console Output:**
```
[Phase 2-B] Geometry changed (hash: e20cad57 -> fa2f9fa3), rebuilding...
```

### Disk Cache (Persistent Metadata)

Metadata is saved to `./.radia_cache/hmatrix_cache.bin` after each H-matrix construction.

**What's Stored:**
- Geometry hash (identifier)
- Element count
- Parameters used (eps, max_rank)
- Construction time
- Memory usage
- Compression ratio
- Timestamp

**What's NOT Stored:**
- Actual H-matrix data (too large)
- Element geometries
- Material properties

**Purpose:**
- Track geometries across sessions
- Usage statistics and performance insights
- Foundation for future ML parameter tuning

**Console Output:**
```
[Phase 3] Cache hit! (hash=e20cad57)
          Previous build: 1.26s (343 elements, eps=0.0001, rank=30)
[Phase 3] Saved to cache (./.radia_cache/hmatrix_cache.bin)
```

---

## Performance Optimization Tips

### 1. Leverage Magnetization-Only Changes

**Best Practice:**
```python
# Create geometry once
geometry = create_magnet_array(7, 7, 7)  # 343 elements

# Iterate over magnetizations (FAST)
for M in magnetization_sweep:
    set_magnetization(geometry, M)  # Change M only
    rad.Solve(geometry, 0.0001, 1000)  # Reuses H-matrix (~30 ms each)
```

**Avoid:**
```python
# Creating new geometry each time (SLOW)
for M in magnetization_sweep:
    rad.UtiDelAll()
    geometry = create_magnet_array(7, 7, 7)  # Rebuilds H-matrix
    set_magnetization(geometry, M)
    rad.Solve(geometry, 0.0001, 1000)  # Full construction (~1000 ms each)
```

### 2. Use Appropriate Subdivision

**Balance accuracy vs performance:**
```python
# Coarse subdivision (fast, less accurate)
rad.ObjDivMag(cube, [5, 5, 5])  # 125 elements, ~0.3s construction

# Medium subdivision (balanced)
rad.ObjDivMag(cube, [7, 7, 7])  # 343 elements, ~1.2s construction

# Fine subdivision (accurate, slower)
rad.ObjDivMag(cube, [10, 10, 10])  # 1000 elements, ~4.5s construction
```

### 3. Batch Multiple Solves

**Efficient workflow:**
```python
# Create all geometries first
geometries = [create_geometry(i) for i in range(10)]

# Solve each once (H-matrix built per geometry)
results = [rad.Solve(g, 0.0001, 1000) for g in geometries]

# Subsequent solves reuse H-matrices (FAST)
for g in geometries:
    for M in magnetizations:
        set_magnetization(g, M)
        rad.Solve(g, 0.0001, 1000)  # Reuses H-matrix
```

### 4. Monitor Cache Hits

**Check console output:**
```
[Phase 3] Cache hit!   -> Good! Previous construction found
[Phase 3] Cache miss   -> Expected for new geometries
[Phase 2-B] Reusing    -> Excellent! H-matrix reused within session
```

---

## Troubleshooting

### Problem: H-Matrix Not Enabled

**Symptoms:**
```
[Auto] N=125 < 200 - using optimized dense solver instead
```

**Cause:** Problem too small (N < 200)

**Solutions:**
1. Increase subdivision: `rad.ObjDivMag(cube, [7, 7, 7])`
2. Force enable: `rad.SolverHMatrixEnable(1, 1e-4, 30)`

### Problem: Slow Construction Time

**Symptoms:**
```
Construction time: 10.5 s  (Expected: ~1-2s for N=343)
```

**Possible Causes:**
1. **Too many elements:** Consider coarser subdivision
2. **Tight parameters:** Relax eps or reduce max_rank
3. **Single-threaded:** Check OpenMP is enabled

**Solutions:**
```python
# Use adaptive parameters (default)
# Or relax manually
rad.SolverHMatrixEnable(1, 5e-4, 20)
```

### Problem: Cache Not Being Reused

**Symptoms:**
```
[Phase 2-B] Geometry changed (hash: xxx -> yyy), rebuilding...
```

**Cause:** Geometry actually changed (elements moved, added, or removed)

**Check:**
1. Are you calling `rad.UtiDelAll()` between solves?
2. Are you creating new objects instead of modifying magnetization?

**Fix:**
```python
# Keep same geometry, change magnetization only
for M in magnetizations:
    rad.ObjM(element, M)  # Change M only
    rad.Solve(container, 0.0001, 1000)  # Reuses H-matrix
```

### Problem: "Failed to create Interaction Matrix"

**Symptoms:**
```
RuntimeError: Failed to create Interaction Matrix.
```

**Cause:** Permanent magnets without material setup

**Fix:**
```python
# For soft magnetic materials (nonlinear)
mat = rad.MatSatIsoFrm([1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759])
rad.MatApl(obj, mat)

# For permanent magnets (fixed magnetization)
# No material needed - magnetization set in ObjRecMag
obj = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1])
```

---

## Advanced Topics

### Understanding H-Matrix Parameters

**eps (ACA tolerance):**
- Controls approximation accuracy
- Smaller = more accurate, slower, larger memory
- Larger = less accurate, faster, smaller memory
- Recommended: 1e-4 to 5e-4 (< 1% error)

**max_rank (maximum rank):**
- Limits low-rank block size
- Smaller = better compression, may lose accuracy
- Larger = higher accuracy, less compression
- Recommended: 20 to 30

**Trade-offs:**
```
eps=1e-6, max_rank=50  ->  Accuracy: 0.1%, Time: slow, Memory: large
eps=1e-4, max_rank=30  ->  Accuracy: 1%,   Time: fast, Memory: medium
eps=5e-4, max_rank=20  ->  Accuracy: 5%,   Time: very fast, Memory: small
```

### Compression Ratio

```
Compression ratio = (H-matrix memory) / (Dense matrix memory)
```

**Typical values:**
- **< 100%:** Excellent compression (H-matrix smaller than dense)
- **100-150%:** Good compression (slight overhead acceptable for speed)
- **> 150%:** Poor compression (consider dense solver or relax parameters)

**Example:**
```
Dense matrix: 4 MB
H-matrix: 9 MB
Compression ratio: 225%  -> More memory, but 33x faster!
```

### Memory Usage

**Dense solver:**
```
Memory = 9 * N^2 * 8 bytes  (9 interaction matrices)

N=343:  9 * 343^2 * 8 = 8.5 MB
N=1000: 9 * 1000^2 * 8 = 72 MB
```

**H-matrix:**
```
Memory depends on problem geometry and parameters
Typical: 1.5-2x dense matrix memory
But much faster for N > 200
```

### Geometry Hash Details

**Algorithm:** Boost-style hash_combine

**Inputs:**
- Element count
- Element center positions (x, y, z)

**Excludes:**
- Magnetization vectors (enables M-only optimization)
- Material properties
- External fields

**Collision resistance:** Similar to SHA hash quality

---

## Examples

### Example 1: Single Large Solve

```python
import radia as rad

# Create soft iron cube (7×7×7 = 343 elements)
mat = rad.MatSatIsoFrm([1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759])
cube = rad.ObjRecMag([0, 0, 0], [100, 100, 100], [0, 0, 0])
rad.MatApl(cube, mat)
rad.ObjDivMag(cube, [7, 7, 7])

# Enable H-matrix explicitly
rad.SolverHMatrixEnable(1, 1e-4, 30)

# Solve with H-matrix
result = rad.Solve(cube, 0.0001, 1000)

# Expected output:
# [Phase 2-B] Using H-matrix solver (N=343)
# Building True H-Matrix... (1.2s)
# Result: converged
```

### Example 2: Magnetization Sweep (Fastest)

```python
# Create geometry once
cube = create_soft_iron_cube([7, 7, 7])

# Enable H-matrix
rad.SolverHMatrixEnable(1, 1e-4, 30)

# First solve - builds H-matrix
rad.Solve(cube, 0.0001, 1000)  # ~1000 ms

# Magnetization sweep - reuses H-matrix
magnetizations = [[0,0,1], [0,1,0], [1,0,0]]
for M in magnetizations:
    # Change magnetization only
    # (implementation depends on your geometry)
    update_magnetization(cube, M)

    # Solve - reuses H-matrix (~30 ms each)
    rad.Solve(cube, 0.0001, 1000)

# Total time: ~1090 ms instead of ~3000 ms
# Speedup: 2.75x
```

### Example 3: Geometry Exploration

```python
# Create multiple geometries
subdivisions = [5, 6, 7, 8]  # 125, 216, 343, 512 elements

# Enable H-matrix (recommended for N >= 200)
rad.SolverHMatrixEnable(1, 1e-4, 30)

for n in subdivisions:
    rad.UtiDelAll()
    cube = create_cube_with_subdivision(n)

    # First solve - builds H-matrix (if N >= 200)
    t0 = time.time()
    rad.Solve(cube, 0.0001, 1000)
    print(f"N={n**3}: {(time.time()-t0)*1000:.1f} ms")

    # Multiple solves - reuse H-matrix
    for i in range(10):
        change_magnetization(cube)
        rad.Solve(cube, 0.0001, 1000)  # Fast!
```

---

## FAQ

**Q: Do I need to configure H-matrix manually?**
A: **Yes!** H-matrix must be explicitly enabled via `rad.SolverHMatrixEnable(1)`. It is NOT automatic.

**Q: Why isn't H-matrix automatic?**
A: By design. This gives users full control and avoids hidden performance behavior. Small models (N < 200) may not benefit from H-matrix, so users should decide when to enable it.

**Q: How do I know if H-matrix is being used?**
A: Check console for "[Phase 2-B] Using H-matrix solver" messages.

**Q: Why is my first solve slow?**
A: H-matrix construction takes ~1s for N=343. Subsequent solves are fast (~30 ms).

**Q: Can I disable H-matrix?**
A: Yes: `rad.SolverHMatrixEnable(0)` (default is disabled)

**Q: What's the cache file?**
A: `./.radia_cache/hmatrix_cache.bin` stores metadata (not H-matrix data).

**Q: Is H-matrix always faster?**
A: No. For **N ≥ 200**, H-matrix is typically faster. For **N < 200**, enabling H-matrix may NOT improve performance due to construction overhead.

**Q: Does H-matrix accelerate `rad.Fld()` or `radia_ngsolve`?**
A: **No**. H-matrix only accelerates `rad.Solve()`. Field evaluation (`rad.Fld()`) and NGSolve integration always use direct calculation.

**Q: How accurate is H-matrix?**
A: Default: < 1% error. Adjustable with eps parameter.

**Q: Does magnetization change rebuild H-matrix?**
A: No! H-matrix is reused (33x speedup).

---

## References

- **Release Notes:** `docs/RELEASE_NOTES_v1.0.10.md`
- **Cache Design:** `docs/DISK_CACHE_DESIGN.md`
- **Phase 2-B Summary:** `docs/PHASE2B_IMPLEMENTATION_SUMMARY.md`
- **Phase 3 Analysis:** `docs/PHASE3_ANALYSIS.md`

## Support

- **Issues:** https://github.com/ksugahar/Radia_NGSolve/issues
- **Repository:** https://github.com/ksugahar/Radia_NGSolve

---

**Version:** 1.0.10
**Last Updated:** 2025-11-13
