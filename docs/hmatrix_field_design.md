# H-Matrix Field Evaluation Design

## Objective

Accelerate `rad.Fld()` field evaluation using H-matrix (Hierarchical Matrix) approximation to achieve **50-200x speedup** for large-scale problems.

## Requirements

Based on user specifications:

1. **Mandatory Caching**: H-matrix must be cached for non-linear iterations
2. **Problem Size**: Only for N ≥ 100 elements (skip small problems)
3. **User Control**: Explicit `use_hmatrix=True` API parameter
4. **Batch Evaluation Only**: Single-point evaluation not supported

## Current Status

### Existing Implementation
- **rad.Solve()**: H-matrix is already used (working)
- **rad.Fld()**: Direct summation O(M×N), **H-matrix NOT used**
- **radhmat.cpp**: Structure exists but B_comp() incomplete

### Performance Baseline
- Batch evaluation: ~4-6x speedup (Python overhead reduction only)
- Target with H-matrix: **50-200x speedup**

## Design Approach

### 1. H-Matrix Field Evaluation Method

#### Kernel Function
For magnetic field evaluation:
```
H(r) = Σ [3(m·r̂)r̂ - m] / (4π|r|³)
```

H-matrix approximates the field matrix:
```
H_ij = K(r_target_i, r_source_j, m_j)
```

where:
- `r_target_i`: observation points (M points)
- `r_source_j`: source element centers (N elements)
- `m_j`: magnetic moments

#### Low-Rank Approximation
For well-separated clusters:
```
H_ij ≈ U_i * S * V_j^T
```

Complexity: O((M+N) log(M+N))

### 2. API Design

#### Python API

```python
# Option 1: Global setting (persistent cache)
rad.SetHMatrixFieldEval(True, eps=1e-6)
H = rad.Fld(obj, 'h', observation_points)  # Uses H-matrix if enabled

# Option 2: Explicit parameter (per-call)
H = rad.Fld(obj, 'h', observation_points, use_hmatrix=True, hmatrix_eps=1e-6)

# Option 3: Batch-specific function
H = rad.FldBatch(obj, 'h', observation_points, use_hmatrix=True)

# Cache management
rad.ClearHMatrixCache()
```

**Recommendation**: Option 1 (Global setting) + Option 3 (Batch function)
- Global setting allows reuse across non-linear iterations
- Batch function clarifies that single points are not supported

#### C++ API

```cpp
// In radTApplication or radTFieldManager
class radTHMatrixFieldCache {
public:
	bool enabled;
	double epsilon;
	std::map<int, std::shared_ptr<HMatrixFieldEvaluator>> cache;

	void Enable(double eps = 1e-6);
	void Disable();
	void Clear();
};

// Field evaluation with H-matrix
void EvaluateFieldBatch(
	radTHandle<radTg>& obj,
	const std::vector<TVector3d>& obs_points,
	std::vector<TVector3d>& field_out,
	radTFieldKey field_key,
	bool use_hmatrix = false,
	double hmatrix_eps = 1e-6
);
```

### 3. Caching Strategy

#### Cache Key
```
cache_key = hash(obj_id, geometry_hash)
```

#### Cache Lifecycle
1. **Construction**: First call with `use_hmatrix=True`
   - Extract geometry from all elements
   - Build cluster tree
   - Construct H-matrix (one-time cost)

2. **Reuse**: Subsequent calls
   - Check cache validity (geometry unchanged)
   - Reuse H-matrix for new observation points

3. **Invalidation**:
   - Geometry modification
   - Explicit `ClearHMatrixCache()` call
   - Memory pressure

#### Memory Management
```
Memory = O(N log N)  for cluster tree
        + O(N log N)  for H-matrix blocks
        ≈ 100 MB for N=10,000
```

### 4. Implementation Steps

#### Phase 1: Core Infrastructure
1. Implement `HMatrixFieldEvaluator` class
2. Add batch evaluation interface
3. Implement caching mechanism

#### Phase 2: Kernel Integration
1. Define magnetic field kernel for HACApK
2. Implement cluster tree construction for observation points
3. Implement H-matrix-vector multiplication for field evaluation

#### Phase 3: API Integration
1. Add Python bindings for global setting
2. Add `rad.FldBatch()` function
3. Modify existing `rad.Fld()` to support batch mode

#### Phase 4: Optimization & Testing
1. Benchmark performance (target: 50-200x)
2. Validate field accuracy (compare with direct)
3. Test non-linear iteration caching

### 5. Technical Details

#### Observation Point Clustering
For M observation points, build separate cluster tree:
```
Cluster_obs = BuildClusterTree(obs_points, max_leaf_size=50)
Cluster_src = BuildClusterTree(source_centers, max_leaf_size=50)
```

#### Admissibility Criterion
Two clusters are admissible for low-rank approximation if:
```
max(diam(C_obs), diam(C_src)) < η * dist(C_obs, C_src)
```
where η = 0.5 (adjustable)

#### Matrix-Vector Multiplication
For field evaluation:
```
H_field = H_matrix * M_vector
```
where:
- `H_matrix`: M × N matrix (H-matrix format)
- `M_vector`: N × 3 magnetic moments
- `H_field`: M × 3 field vectors

Complexity: O((M+N) log(M+N)) instead of O(M×N)

### 6. Performance Estimates

#### Problem Sizes
| N (elements) | M (obs points) | Direct O(M×N) | H-matrix O((M+N)logN) | Speedup |
|--------------|----------------|---------------|------------------------|---------|
| 100          | 100            | 10k           | 1.3k                   | 8x      |
| 1,000        | 1,000          | 1M            | 20k                    | 50x     |
| 10,000       | 10,000         | 100M          | 200k                   | 500x    |
| 10,000       | 100            | 1M            | 140k                   | 7x      |

#### Memory Usage
| N (elements) | Cluster Tree | H-matrix | Total |
|--------------|--------------|----------|-------|
| 100          | ~50 KB       | ~500 KB  | 0.5 MB |
| 1,000        | ~500 KB      | ~5 MB    | 5 MB   |
| 10,000       | ~5 MB        | ~50 MB   | 55 MB  |

### 7. Validation Strategy

#### Accuracy Tests
1. **Direct vs H-matrix**: Compare field values
2. **Error tolerance**: |H_hmat - H_direct|/|H_direct| < ε
3. **Symmetric geometry**: Check field symmetry

#### Performance Tests
1. **Construction time**: One-time cost acceptable?
2. **Evaluation time**: 50-200x speedup achieved?
3. **Memory usage**: Within acceptable limits?

#### Integration Tests
1. **Non-linear iteration**: Cache reuse working?
2. **Geometry modification**: Cache invalidation working?
3. **Mixed evaluation**: Batch + single-point coexistence?

## Implementation Priority

1. ✅ **Design document** (this document)
2. ⏳ **Core HMatrixFieldEvaluator class**
3. ⏳ **Batch evaluation kernel**
4. ⏳ **Python API bindings**
5. ⏳ **Caching mechanism**
6. ⏳ **Performance benchmarks**

## Open Questions

1. **Should we support single-point evaluation with H-matrix?**
   - User requirement: "No" (batch only)
   - Implementation: Can be done but adds complexity

2. **How to handle field components (B, H, A, M)?**
   - H-matrix works for H-field directly
   - B = μ₀H in vacuum, scaled version
   - A (vector potential) requires different kernel

3. **OpenMP parallelization?**
   - HACApK supports OpenMP
   - Can we enable for both construction and evaluation?

4. **Adaptive recompression?**
   - As observation points change, should we recompress?
   - Trade-off: accuracy vs performance

## References

- HACApK Library: https://github.com/lwk205/HACApK-LH-Cimplm
- H-matrix theory: Bebendorf, "Hierarchical Matrices" (2008)
- Radia existing H-matrix: `src/core/radhmat.cpp`
