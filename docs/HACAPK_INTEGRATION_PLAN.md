# HACApK Integration Plan for Radia

**Date:** 2025-11-07
**Version:** 1.0
**Status:** Proposal

## 1. Executive Summary

This document outlines the plan to integrate HACApK (Hierarchical Adaptive Cross Approximation on Kernel matrices) into Radia's magnetic field computation engine to accelerate field calculations for systems with many magnetic elements.

## 2. Background

### 2.1 Current Radia Architecture

Radia's field calculation structure:
- **Base class**: `radTg3d` - all magnetic field sources inherit from this
- **Key method**: `B_comp(radTField*)` - computes magnetic field at a point
- **Group container**: `radTGroup` - contains multiple field sources
- **Field sources**:
  - `radTFlmLinCur` - filament conductors
  - `radTCoefficientFunctionFieldSource` - Python callback-based fields
  - Magnetized volumes (polyhedra, cylinders, etc.)

### 2.2 HACApK Capabilities

HACApK provides:
- **H-matrix** (Hierarchical Matrix) representation
- **ACA** (Adaptive Cross Approximation) for low-rank approximation
- **OpenMP parallelization** for shared-memory systems
- **Fast matrix-vector multiplication**: O(N log N) vs O(N²)

## 3. Integration Strategy

### 3.1 Target Use Case

**Primary application**: Accelerate magnetic field calculations for systems with many magnetic elements (N > 100).

For a system with N magnetic elements:
- **Traditional approach**: O(N) per field evaluation, O(N²M) for M points
- **HACApK approach**: O(N log N) with H-matrix precomputation

### 3.2 Integration Approach

Create a new field source class: **`radTHMatrixFieldSource`**

```
radTg3d
  └── radTHMatrixFieldSource  (NEW)
```

This class will:
1. Accept a group of magnetic elements
2. Build H-matrix representation using HACApK
3. Compute fields using fast H-matrix-vector multiplication
4. Support all field types: B, H, A, M

## 4. Detailed Design

### 4.1 New Files

Create in `src/core/`:
- `radhmat.h` - Header for H-matrix field source
- `radhmat.cpp` - Implementation

### 4.2 Class Structure

```cpp
class radTHMatrixFieldSource : public radTg3d {
public:
	// Configuration
	struct HMatrixConfig {
		double eps;         // ACA tolerance (default: 1e-6)
		int max_rank;       // Maximum rank for low-rank blocks
		int min_cluster;    // Minimum cluster size
		bool use_openmp;    // Enable OpenMP (default: true)
	};

private:
	radTmhg source_elements;           // Original magnetic elements
	hacapk::HMatrix* hmatrix;          // H-matrix representation
	std::vector<hacapk::Point3D> eval_points;  // Evaluation points
	HMatrixConfig config;
	bool is_built;

public:
	radTHMatrixFieldSource(radTGroup* group, const HMatrixConfig& cfg);
	~radTHMatrixFieldSource();

	// Build H-matrix from source elements
	int BuildHMatrix();

	// Field computation (overrides radTg3d)
	void B_comp(radTField* FieldPtr) override;
	void B_intComp(radTField* FieldPtr) override;

	// Batch field computation (efficient)
	void B_comp_batch(std::vector<radTField*>& fields);

	// Utilities
	int Type_g3d() override { return 100; }  // New type for H-matrix
	void Dump(std::ostream&, int ShortSign = 0) override;
	radTg3dGraphPresent* CreateGraphPresent() override;
	int DuplicateItself(radThg& hg, radTApplication*, char) override;
	int SizeOfThis() override { return sizeof(radTHMatrixFieldSource); }
};
```

### 4.3 Algorithm Flow

#### 4.3.1 Construction Phase

```cpp
// User creates H-matrix field source from a group
radTGroup* group = ...;  // Contains N magnetic elements
radTHMatrixFieldSource::HMatrixConfig cfg;
cfg.eps = 1e-6;
cfg.max_rank = 50;
radTHMatrixFieldSource* hmat_source = new radTHMatrixFieldSource(group, cfg);
hmat_source->BuildHMatrix();
```

#### 4.3.2 H-Matrix Construction

1. **Extract geometry** from source elements:
   ```cpp
   for each element in group:
       extract center position
       extract bounding box
       store magnetic moment / current
   ```

2. **Build cluster tree** using HACApK:
   ```cpp
   hacapk::build_cluster_tree(points, bbox, min_cluster_size)
   ```

3. **Define kernel function** for Biot-Savart law:
   ```cpp
   double kernel(int i, int j) {
       // Magnetic field from element j at position of element i
       // Uses Biot-Savart law or dipole approximation
       return compute_field_influence(elem_j, pos_i);
   }
   ```

4. **Build H-matrix** with ACA:
   ```cpp
   hacapk::build_hmatrix(cluster_tree, kernel, eps, max_rank)
   ```

#### 4.3.3 Field Evaluation

```cpp
void radTHMatrixFieldSource::B_comp(radTField* FieldPtr) {
    if (!is_built) {
        // Fallback to direct computation
        radTGroup::B_comp(FieldPtr);
        return;
    }

    // Use H-matrix for fast field evaluation
    TVector3d P = FieldPtr->P;

    // Find nearest cluster
    int cluster_idx = find_cluster_for_point(P);

    // Evaluate field using H-matrix-vector multiplication
    TVector3d B = hmatrix_matvec(P, cluster_idx);

    if (FieldPtr->FieldKey.B_) FieldPtr->B += B;
    if (FieldPtr->FieldKey.H_) FieldPtr->H += B;  // Assuming vacuum
}
```

### 4.4 Batch Evaluation (High Performance)

For M evaluation points:

```cpp
void radTHMatrixFieldSource::B_comp_batch(std::vector<radTField*>& fields) {
    int M = fields.size();

    // Prepare input vector
    std::vector<TVector3d> positions(M);
    for (int i = 0; i < M; i++) {
        positions[i] = fields[i]->P;
    }

    // H-matrix matrix-vector multiplication (parallelized with OpenMP)
    std::vector<TVector3d> B_values = hmatrix_matvec_batch(positions);

    // Store results
    #pragma omp parallel for
    for (int i = 0; i < M; i++) {
        if (fields[i]->FieldKey.B_) fields[i]->B += B_values[i];
        if (fields[i]->FieldKey.H_) fields[i]->H += B_values[i];
    }
}
```

## 5. Implementation Phases

### Phase 1: Basic Integration (1-2 weeks)
- [ ] Create `radhmat.h` and `radhmat.cpp`
- [ ] Implement `radTHMatrixFieldSource` skeleton
- [ ] Link with HACApK library in CMakeLists.txt
- [ ] Basic B_comp implementation (without optimization)
- [ ] Unit tests

### Phase 2: H-Matrix Construction (2-3 weeks)
- [ ] Implement geometry extraction from radTGroup
- [ ] Implement cluster tree building
- [ ] Implement kernel function for Biot-Savart law
- [ ] H-matrix construction with ACA
- [ ] Validation against direct calculation

### Phase 3: Optimized Field Evaluation (1-2 weeks)
- [ ] Implement fast H-matrix-vector multiplication
- [ ] Implement batch evaluation
- [ ] OpenMP optimization
- [ ] Performance benchmarks

### Phase 4: Python Interface (1 week)
- [ ] Add Python bindings (radpy.cpp)
- [ ] Create example scripts
- [ ] Documentation

### Phase 5: Testing and Validation (1-2 weeks)
- [ ] Comprehensive test suite
- [ ] Accuracy validation
- [ ] Performance benchmarks
- [ ] Documentation

## 6. Performance Expectations

### 6.1 Complexity Analysis

| Operation | Traditional | H-Matrix |
|-----------|-------------|----------|
| Construction | - | O(N log N) |
| Single point | O(N) | O(log N) |
| M points | O(N·M) | O((N+M) log N) |

### 6.2 Expected Speedup

For N = 1000 elements, M = 10000 evaluation points:
- **Traditional**: 10^7 operations
- **H-Matrix**: ~10^5 operations
- **Speedup**: ~100x

## 7. CMake Integration

Update `CMakeLists.txt`:

```cmake
# HACApK library
add_subdirectory(src/ext/HACApK_LH-Cimplm)

# Radia core with H-matrix support
set(RADIA_SOURCES
    ...
    src/core/radhmat.cpp
)

# Link HACApK
target_link_libraries(radia PRIVATE hacapk)
target_include_directories(radia PRIVATE src/ext/HACApK_LH-Cimplm)
```

## 8. Python API

```python
import radia as rad

# Create magnetic elements
group = rad.ObjCnt([
    rad.ObjRecMag([x, y, z], [10, 10, 10], [0, 0, 1])
    for x, y, z in magnet_positions
])

# Create H-matrix field source
hmat = rad.ObjHMatrix(group, eps=1e-6, max_rank=50)

# Evaluate field (accelerated)
B = rad.Fld(hmat, 'b', [0, 0, 100])

# Batch evaluation (very efficient)
points = [[x, y, z] for x, y, z in evaluation_points]
B_array = rad.FldLst(hmat, 'b', points)
```

## 9. Testing Strategy

### 9.1 Accuracy Tests

Compare H-matrix results against direct calculation:
- Single magnet
- 2 magnets (near and far)
- 10 magnets
- 100 magnets
- 1000 magnets

Tolerance: |B_hmat - B_direct| < 1e-6 T

### 9.2 Performance Tests

Measure computation time:
- N = [10, 50, 100, 500, 1000, 5000]
- M = [100, 1000, 10000]
- Compare: direct vs H-matrix

### 9.3 Convergence Tests

Test ACA convergence:
- eps = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
- Measure: rank, error, time

## 10. Risks and Mitigation

| Risk | Impact | Mitigation |
|------|--------|------------|
| H-matrix construction overhead | High | Only use for N > 100 |
| ACA accuracy issues | Medium | Adjustable tolerance parameter |
| Memory consumption | Medium | Adaptive rank control |
| Complex geometry handling | Medium | Fallback to direct calculation |

## 11. Success Criteria

- ✅ Speedup > 10x for N > 1000
- ✅ Error < 1e-6 T compared to direct calculation
- ✅ Memory usage < 2× direct method
- ✅ Python API functional
- ✅ All tests pass

## 12. Future Enhancements

- GPU acceleration using CUDA
- Distributed memory parallelization (MPI)
- Adaptive mesh refinement
- Time-dependent fields
- Optimization for specific geometries (coils, arrays)

## 13. References

1. HACApK original paper: Ida & Iwashita (2015)
2. H-matrix theory: Hackbusch (1999)
3. ACA algorithm: Bebendorf (2000)
4. Radia documentation: ESRF (1997-2025)

---

**Next Steps:**
1. Review and approve this plan
2. Create development branch
3. Begin Phase 1 implementation
