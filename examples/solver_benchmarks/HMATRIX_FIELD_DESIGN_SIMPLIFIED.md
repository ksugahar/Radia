# H-Matrix Field Evaluation - Simplified Design

**Date**: 2025-11-08
**Status**: Design Review (Simplified based on user feedback)
**Goal**: Implement H-matrix acceleration for rad.Fld() batch evaluation

---

## User Requirements (Updated)

✅ **Cache H-matrix for non-linear iterations**
✅ **N ≥ 100 only** (no small problem support needed)
✅ **Explicit API** (use_hmatrix=True/False)
✅ **Batch evaluation only** (no single-point optimization)

---

## 1. Simplified Architecture

### What We're Building

```python
# Explicit H-matrix usage
magnet = rad.ObjCnt([blocks...])  # N ≥ 100 elements
rad.Solve(magnet, precision, max_iter)

# Batch field evaluation with H-matrix
points = [[x1,y1,z1], [x2,y2,z2], ...]  # M points
B_values = rad.Fld(magnet, 'b', points, use_hmatrix=True)  # Explicit flag

# H-matrix is automatically cached for reuse
B_values2 = rad.Fld(magnet, 'b', points, use_hmatrix=True)  # Reuses cache
```

### What We're NOT Building

❌ Single-point evaluation optimization
❌ N < 100 support
❌ Automatic H-matrix detection (transparent)
❌ Fallback to direct calculation

---

## 2. Non-Linear Iteration Use Case

### Why Caching is Critical

In non-linear materials, rad.Solve() performs iterations:

```python
# Non-linear relaxation iteration
for iteration in range(max_iter):
    # 1. Compute demagnetization field at element centers
    #    → N field evaluations (N observation points)

    # 2. Update magnetization based on B-H curve

    # 3. Check convergence
```

**Without H-matrix cache**:
```
Each iteration:
  - Build H-matrix: 50 ms
  - Evaluate N points: 20 ms
  Total per iteration: 70 ms
  10 iterations: 700 ms
```

**With H-matrix cache**:
```
First iteration:
  - Build H-matrix: 50 ms (one-time)
  - Evaluate N points: 20 ms
Subsequent iterations:
  - Reuse cached H-matrix: 0 ms
  - Evaluate N points: 20 ms
  Total: 50 + 10×20 = 250 ms

Speedup: 700/250 = 2.8x
```

### NGSolve Integration

GridFunction.Set() calls Evaluate() many times:

```python
# NGSolve mesh with M=5000 vertices
mesh = Mesh(...)
fes = HCurl(mesh)
gf = GridFunction(fes)

# Without cache: 1250 H-matrix builds × 50ms = 62,500 ms
# With cache: 1 H-matrix build × 50ms + 1250 evals × 0.02ms = 75 ms

# Speedup: 833x
```

---

## 3. Implementation Plan (Simplified)

### Phase 1: H-Matrix Caching Infrastructure

**File**: `src/core/radhmat.h`

```cpp
class radTHMatrixFieldSource {
private:
	// Cached H-matrix for field evaluation
	struct FieldHMatrixCache {
		std::vector<std::unique_ptr<hacapk::HMatrix>> hmatrices;  // 9 matrices (3×3 tensor)
		std::vector<TVector3d> observation_points;  // M cached points
		size_t memory_usage;
		bool is_valid;

		FieldHMatrixCache() : memory_usage(0), is_valid(false) {}
	};

	FieldHMatrixCache field_cache;

public:
	// Build and cache H-matrix for observation points
	bool BuildAndCacheFieldHMatrix(const std::vector<TVector3d>& obs_points);

	// Check if cache is valid for these observation points
	bool IsCacheValid(const std::vector<TVector3d>& obs_points) const;

	// Invalidate cache (e.g., when magnetization changes)
	void InvalidateCache();

	// Evaluate field using cached H-matrix
	std::vector<TVector3d> EvaluateBatchCached(const std::vector<TVector3d>& obs_points);
};
```

### Phase 2: Cache Validation Logic

```cpp
bool radTHMatrixFieldSource::IsCacheValid(const std::vector<TVector3d>& obs_points) const
{
	if(!field_cache.is_valid) return false;

	// Check if observation points match exactly
	if(obs_points.size() != field_cache.observation_points.size()) {
		return false;
	}

	// Point-by-point comparison
	for(size_t i = 0; i < obs_points.size(); i++) {
		const TVector3d& p1 = obs_points[i];
		const TVector3d& p2 = field_cache.observation_points[i];

		// Check if points are identical (within tolerance)
		double dx = p1.x - p2.x;
		double dy = p1.y - p2.y;
		double dz = p1.z - p2.z;
		double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

		if(dist > 1e-10) {
			return false;  // Points don't match
		}
	}

	return true;  // Cache is valid
}
```

### Phase 3: Batch Evaluation with Caching

**File**: `src/core/radhmat.cpp`

```cpp
std::vector<TVector3d> radTHMatrixFieldSource::EvaluateBatchCached(
	const std::vector<TVector3d>& obs_points)
{
	int M = obs_points.size();

	std::cout << "[HMatrix] Batch evaluation: M=" << M << " points, N=" << num_elements << " sources" << std::endl;

	// Check cache validity
	bool cache_valid = IsCacheValid(obs_points);

	if(!cache_valid) {
		std::cout << "[HMatrix] Cache miss - building new H-matrix" << std::endl;

		auto start = std::chrono::high_resolution_clock::now();

		// Build and cache H-matrix
		if(!BuildAndCacheFieldHMatrix(obs_points)) {
			throw std::runtime_error("Failed to build field H-matrix");
		}

		auto end = std::chrono::high_resolution_clock::now();
		double build_time = std::chrono::duration<double>(end - start).count();

		std::cout << "[HMatrix] H-matrix built in " << build_time * 1000 << " ms" << std::endl;
		std::cout << "[HMatrix] Memory usage: " << field_cache.memory_usage / 1024.0 / 1024.0 << " MB" << std::endl;
	} else {
		std::cout << "[HMatrix] Cache hit - reusing H-matrix" << std::endl;
	}

	// Evaluate field using cached H-matrix
	auto start = std::chrono::high_resolution_clock::now();

	std::vector<TVector3d> B_values(M);

	// Parallel evaluation using cached H-matrix
	#pragma omp parallel for schedule(static) if(config.use_openmp && M > 10)
	for(int i = 0; i < M; i++) {
		B_values[i] = EvaluatePointHMatrix(i);
	}

	auto end = std::chrono::high_resolution_clock::now();
	double eval_time = std::chrono::duration<double>(end - start).count();

	std::cout << "[HMatrix] Evaluation: " << eval_time * 1000 << " ms" << std::endl;
	std::cout << "[HMatrix] Time per point: " << eval_time / M * 1e6 << " us" << std::endl;

	return B_values;
}
```

### Phase 4: Build and Cache H-Matrix

```cpp
bool radTHMatrixFieldSource::BuildAndCacheFieldHMatrix(
	const std::vector<TVector3d>& obs_points)
{
	int M = obs_points.size();

	// Invalidate old cache
	InvalidateCache();

	// Convert to hacapk::Point3D
	std::vector<hacapk::Point3D> target_points;
	target_points.reserve(M);
	for(const auto& p : obs_points) {
		target_points.emplace_back(p.x, p.y, p.z);
	}

	std::vector<hacapk::Point3D> source_points;
	source_points.reserve(num_elements);
	for(int i = 0; i < num_elements; i++) {
		source_points.emplace_back(
			element_positions[3*i],
			element_positions[3*i+1],
			element_positions[3*i+2]
		);
	}

	// Resize cache for 9 H-matrices (3×3 tensor)
	field_cache.hmatrices.resize(9);
	field_cache.memory_usage = 0;

	// Build 9 H-matrices in parallel
	std::cout << "[HMatrix] Building 9 H-matrices (3×3 tensor components)..." << std::endl;

	#pragma omp parallel for schedule(dynamic) if(config.use_openmp)
	for(int idx = 0; idx < 9; idx++) {
		int row = idx / 3;
		int col = idx % 3;

		#pragma omp critical
		{
			std::cout << "  Component [" << row << "][" << col << "]... " << std::flush;
		}

		// Setup Biot-Savart kernel
		BiotSavartKernel kernel;
		kernel.component_row = row;
		kernel.component_col = col;
		kernel.source_positions = &source_points;
		kernel.obs_positions = &target_points;
		kernel.magnetic_moments = &element_moments;

		// HACApK parameters
		hacapk::ControlParams params;
		params.eps_aca = config.eps;
		params.leaf_size = config.min_cluster_size;
		params.eta = 2.0;
		params.print_level = 0;  // Suppress HACApK output

		// Build H-matrix
		field_cache.hmatrices[idx] = hacapk::build_hmatrix(
			target_points,  // M observation points
			source_points,  // N source points
			kernel,
			&kernel,
			params
		);

		if(!field_cache.hmatrices[idx]) {
			#pragma omp critical
			{
				std::cerr << "Failed!" << std::endl;
			}
			return false;
		}

		size_t mem = field_cache.hmatrices[idx]->memory_usage();

		#pragma omp critical
		{
			field_cache.memory_usage += mem;
			std::cout << "rank=" << field_cache.hmatrices[idx]->ktmax
			          << ", mem=" << mem / 1024 << " KB" << std::endl;
		}
	}

	// Store observation points for cache validation
	field_cache.observation_points = obs_points;
	field_cache.is_valid = true;

	return true;
}
```

### Phase 5: Evaluate Single Point Using H-Matrix

```cpp
TVector3d radTHMatrixFieldSource::EvaluatePointHMatrix(int obs_index) const
{
	// B(r_i) = Σ_j K_ij · m_j
	// K_ij is 3×3 Biot-Savart tensor
	// We have 9 H-matrices representing this tensor

	TVector3d B(0, 0, 0);

	// For each component of B (x, y, z)
	for(int row = 0; row < 3; row++) {
		double B_component = 0.0;

		// Sum contributions from each moment component
		for(int col = 0; col < 3; col++) {
			int idx = row * 3 + col;

			// Get contribution from H-matrix
			// This requires HACApK to support extracting single row
			// For now, we use matrix-vector product for entire column

			// Create unit vector for moment component
			std::vector<double> unit_moment(num_elements, 0.0);
			for(int j = 0; j < num_elements; j++) {
				unit_moment[j] = element_moments[3*j + col];
			}

			// H-matrix × moment vector
			std::vector<double> result(field_cache.observation_points.size());
			field_cache.hmatrices[idx]->matvec(unit_moment.data(), result.data());

			B_component += result[obs_index];
		}

		B[row] = B_component;
	}

	return B;
}
```

**Note**: This is inefficient (computes all M points even for single point).
Better approach: Compute all M points at once.

### Phase 6: Optimized Batch Evaluation (All Points at Once)

```cpp
std::vector<TVector3d> radTHMatrixFieldSource::EvaluateBatchCached(
	const std::vector<TVector3d>& obs_points)
{
	int M = obs_points.size();

	// ... (cache validation code) ...

	// Evaluate ALL points at once using H-matrix-vector products
	std::vector<TVector3d> B_values(M, TVector3d(0, 0, 0));

	// For each component of the 3×3 tensor
	for(int row = 0; row < 3; row++) {
		for(int col = 0; col < 3; col++) {
			int idx = row * 3 + col;

			// Extract moment component for all elements
			std::vector<double> moment_vec(num_elements);
			for(int j = 0; j < num_elements; j++) {
				moment_vec[j] = element_moments[3*j + col];
			}

			// H-matrix × moment vector → contribution to B
			std::vector<double> contrib(M);
			field_cache.hmatrices[idx]->matvec(moment_vec.data(), contrib.data());

			// Accumulate contributions
			for(int i = 0; i < M; i++) {
				B_values[i][row] += contrib[i];
			}
		}
	}

	return B_values;
}
```

---

## 4. Python API

### 4.1. New Parameter for rad.Fld()

**File**: `src/python/radpy.cpp`

```cpp
// Modified signature
py::object RadiaFld(int ElemKey, const std::string& FieldChar, py::object points, bool use_hmatrix = false)
{
	// Check if points is a list
	if(py::isinstance<py::list>(points)) {
		py::list points_list = points.cast<py::list>();

		// Single point: [[x, y, z]] or [x, y, z]
		if(points_list.size() == 3 && py::isinstance<py::float_>(points_list[0])) {
			// Single point [x, y, z]
			if(use_hmatrix) {
				std::cerr << "[Warning] use_hmatrix=True ignored for single-point evaluation" << std::endl;
			}
			return RadiaFld_SinglePoint(ElemKey, FieldChar, points_list);
		}

		// Multiple points [[x1,y1,z1], [x2,y2,z2], ...]
		if(use_hmatrix) {
			return RadiaFld_Batch_HMatrix(ElemKey, FieldChar, points_list);
		} else {
			return RadiaFld_Batch_Direct(ElemKey, FieldChar, points_list);
		}
	}

	throw std::runtime_error("Invalid points format");
}

// Export to Python
m.def("Fld", &RadiaFld,
      py::arg("obj"),
      py::arg("field_type"),
      py::arg("points"),
      py::arg("use_hmatrix") = false,
      "Compute field at point(s). Set use_hmatrix=True for batch evaluation with H-matrix acceleration.");
```

### 4.2. Python Usage Examples

```python
import radia as rad

# Create magnet (N ≥ 100 elements)
magnet = create_magnet(n=7)  # 343 elements
rad.Solve(magnet, 0.0001, 1000)

# Example 1: Direct batch evaluation (current method)
points = [[x, y, z] for x in range(-30, 30, 2) for y in [-10, 0, 10] for z in [0]]
B_direct = rad.Fld(magnet, 'b', points, use_hmatrix=False)  # Direct sum
# Time: ~150 ms for 1000 points

# Example 2: H-matrix batch evaluation (new method)
B_hmatrix = rad.Fld(magnet, 'b', points, use_hmatrix=True)  # H-matrix
# First call: Build H-matrix (50ms) + Evaluate (20ms) = 70 ms
# Subsequent calls: Reuse cache (0ms) + Evaluate (20ms) = 20 ms

# Example 3: Repeated evaluations (non-linear iteration scenario)
for iteration in range(10):
	# Magnetization changes, but observation points stay the same
	B = rad.Fld(magnet, 'b', points, use_hmatrix=True)
	# Only first iteration builds H-matrix (50ms)
	# Iterations 2-10 reuse cache (20ms each)
	# Total: 50 + 9×20 = 230 ms vs 10×150 = 1500 ms direct
	# Speedup: 6.5x
```

---

## 5. NGSolve Integration

### 5.1. Modified radia_ngsolve.cpp

**File**: `src/python/radia_ngsolve.cpp`

```cpp
class RadiaFieldCF : public CoefficientFunction
{
private:
	int radia_obj;
	std::string field_type;
	bool use_hmatrix;  // New parameter

public:
	RadiaFieldCF(int obj, const std::string& ftype = "b", bool use_hmat = true)
		: CoefficientFunction(3)
		, radia_obj(obj)
		, field_type(ftype)
		, use_hmatrix(use_hmat)  // Default to True for NGSolve
	{
	}

	virtual void Evaluate(const BaseMappedIntegrationRule& mir,
	                      BareSliceMatrix<> result) const override
	{
		py::gil_scoped_acquire acquire;
		size_t npts = mir.Size();

		// Extract evaluation points
		py::list points_list;
		for(size_t i = 0; i < npts; i++) {
			auto pnt = mir[i].GetPoint();
			py::list coords;
			coords.append(pnt[0] * 1000.0);  // m → mm
			coords.append(pnt[1] * 1000.0);
			coords.append(pnt[2] * 1000.0);
			points_list.append(coords);
		}

		// Call rad.Fld with use_hmatrix flag
		py::module_ rad = py::module_::import("radia");
		py::object field_results = rad.attr("Fld")(
			radia_obj,
			field_type,
			points_list,
			py::arg("use_hmatrix") = use_hmatrix  // Use H-matrix for NGSolve
		);

		// Extract results
		py::list results_list = field_results.cast<py::list>();
		for(size_t i = 0; i < npts; i++) {
			py::list field_list = results_list[i].cast<py::list>();
			result(0, i) = field_list[0].cast<double>() * scale;
			result(1, i) = field_list[1].cast<double>() * scale;
			result(2, i) = field_list[2].cast<double>() * scale;
		}
	}
};

// Python export
m.def("RadiaField", [](int obj, std::string field_type, bool use_hmatrix) {
	return std::make_shared<RadiaFieldCF>(obj, field_type, use_hmatrix);
}, py::arg("obj"), py::arg("field_type") = "b", py::arg("use_hmatrix") = true);
```

### 5.2. Python Usage with NGSolve

```python
import radia_ngsolve
import ngsolve as ngs

# Create Radia magnet
magnet = rad.ObjCnt([blocks...])
rad.Solve(magnet, 0.0001, 1000)

# Create NGSolve mesh
mesh = ngs.Mesh(...)
fes = ngs.HCurl(mesh)
gf = ngs.GridFunction(fes)

# Create Radia CoefficientFunction with H-matrix
B_cf = radia_ngsolve.RadiaField(magnet, 'b', use_hmatrix=True)

# Set GridFunction
gf.Set(B_cf)  # H-matrix is cached and reused for all element evaluations!

# Speedup: 833x for M=5000 vertices (1 H-matrix build vs 1250)
```

---

## 6. Expected Performance (Simplified)

### 6.1. Batch Evaluation (No Cache)

| M | N | Direct (ms) | H-Matrix (ms) | Build | Eval | Speedup |
|---|---|-------------|---------------|-------|------|---------|
| 100 | 100 | 3 | 25 | 20 | 5 | 0.12x ❌ |
| 1000 | 500 | 150 | 40 | 25 | 15 | 3.8x ✅ |
| 5000 | 1000 | 1500 | 100 | 60 | 40 | 15x ✅ |
| 10000 | 5000 | 15000 | 300 | 180 | 120 | 50x ✅ |

### 6.2. Batch Evaluation (With Cache)

**Scenario**: 10 iterations, same observation points

| M | N | Direct (ms) | H-Matrix (ms) | Speedup |
|---|---|-------------|---------------|---------|
| 1000 | 500 | 1500 | 175 | 8.6x ✅ |
| 5000 | 1000 | 15000 | 460 | 32x ✅ |
| 10000 | 5000 | 150000 | 1380 | 109x ✅ |

**Breakdown (M=5000, N=1000)**:
```
Direct: 10 iterations × 1500 ms = 15000 ms

H-Matrix:
  First iteration: Build (60ms) + Eval (40ms) = 100 ms
  Iterations 2-10: Eval only (40ms) × 9 = 360 ms
  Total: 460 ms

Speedup: 15000 / 460 = 32x
```

### 6.3. NGSolve Integration (With Cache)

**GridFunction.Set() with M=5000 vertices, N=1000 elements**

```
Current (batch evaluation, no cache):
  1250 element evaluations × 4 points each
  Time per evaluation: 0.11 ms
  Total: 1250 × 0.11 = 137.5 ms

With H-Matrix cache:
  First evaluation: Build H-matrix (60ms) + Eval (40ms) = 100 ms
  Evaluations 2-1250: Reuse cache, eval only (0.04ms) × 1249 = 50 ms
  Total: 150 ms

Wait, this is wrong! NGSolve calls with different points each time.
```

**Correction**: NGSolve calls Evaluate() with 4-7 points per element, different points each call.

**Strategy**: Collect ALL evaluation points from ALL elements first, then evaluate once.

This requires modifying GridFunction.Set() → See SetBatch.cpp from forum post.

**With SetBatch()**:
```
Collect all M=5000 vertex positions
Call rad.Fld(magnet, 'b', all_vertices, use_hmatrix=True) once
  Build H-matrix: 60 ms
  Evaluate: 40 ms
  Total: 100 ms

Project to DOFs
  Total: ~10 ms

Total time: 110 ms vs current 137.5 ms

Speedup: 1.25x (not impressive, but overhead reduced)
```

**For larger problems (M=20000, N=5000)**:
```
Current: 20000 × 0.28 ms = 5600 ms

SetBatch + H-matrix:
  Build: 180 ms
  Eval: 120 ms
  Project: 40 ms
  Total: 340 ms

Speedup: 5600 / 340 = 16.5x
```

---

## 7. Implementation Timeline (Simplified)

### Week 1: Core Implementation (5 days)

**Day 1**: Fix geometry extraction
- ✅ Get actual element centers
- ✅ Extract magnetic moments correctly
- ✅ Test with simple magnet

**Day 2**: Implement Biot-Savart kernel
- ✅ Full 3×3 tensor kernel
- ✅ Test kernel accuracy vs analytical solution
- ✅ Verify with unit tests

**Day 3**: Implement H-matrix caching
- ✅ FieldHMatrixCache structure
- ✅ IsCacheValid() logic
- ✅ BuildAndCacheFieldHMatrix()

**Day 4**: Implement batch evaluation
- ✅ EvaluateBatchCached()
- ✅ H-matrix-vector products (9 matrices)
- ✅ Test cache hit/miss

**Day 5**: Testing and debugging
- ✅ Accuracy tests (H-matrix vs direct)
- ✅ Cache validation tests
- ✅ Memory leak checks

### Week 2: API Integration (3 days)

**Day 6**: Python API
- ✅ Add use_hmatrix parameter to rad.Fld()
- ✅ RadiaFld_Batch_HMatrix() implementation
- ✅ Python binding tests

**Day 7**: NGSolve integration
- ✅ Modify RadiaFieldCF with use_hmatrix
- ✅ Test with GridFunction.Set()
- ✅ Verify cache reuse

**Day 8**: Enable in CMakeLists.txt
- ✅ Uncomment radhmat.cpp
- ✅ Build and test
- ✅ Fix any linker errors

### Week 3: Optimization & Benchmarking (2 days)

**Day 9**: Performance optimization
- ✅ OpenMP tuning
- ✅ Memory optimization
- ✅ Cache policy tuning

**Day 10**: Benchmarking and documentation
- ✅ Run benchmark suite
- ✅ Compare with design targets
- ✅ Update documentation

**Total: 10 working days**

---

## 8. Success Criteria (Simplified)

### Performance Targets

| Scenario | Speedup Target | Priority |
|----------|----------------|----------|
| Batch (M=1000, N=500, no cache) | 3x | Medium |
| Batch (M=5000, N=1000, no cache) | 15x | High |
| Batch (M=5000, 10 iter, cache) | 30x | **Critical** |
| NGSolve (M=5000, N=1000) | 10x | High |

### Quality Targets

- ✅ Accuracy: < 0.1% error vs direct sum
- ✅ Memory: Cache size < 100 MB for typical problems
- ✅ Cache hit rate: > 90% in non-linear iterations
- ✅ No memory leaks

---

## 9. Key Decisions (Updated)

### Decision 1: Cache Strategy

**Chosen**: Mandatory caching, invalidate on magnetization change

```cpp
// Cache is always enabled
// Invalidate when:
//   1. Observation points change
//   2. Magnetization changes (rad.Solve called)
//   3. Geometry changes
```

### Decision 2: API Design

**Chosen**: Explicit use_hmatrix parameter

```python
# User explicitly enables H-matrix
B = rad.Fld(magnet, 'b', points, use_hmatrix=True)

# NGSolve: H-matrix enabled by default
B_cf = radia_ngsolve.RadiaField(magnet, 'b', use_hmatrix=True)
```

### Decision 3: Minimum Problem Size

**Chosen**: N ≥ 100, M ≥ 100

```cpp
if(N < 100 || M < 100) {
	throw std::runtime_error("H-matrix requires N ≥ 100 and M ≥ 100");
}
```

### Decision 4: Single-Point Evaluation

**Chosen**: Not supported, batch only

```python
# This will raise an error with use_hmatrix=True
B = rad.Fld(magnet, 'b', [0, 0, 0], use_hmatrix=True)  # Error!

# Use batch with single point instead
B = rad.Fld(magnet, 'b', [[0, 0, 0]], use_hmatrix=True)  # OK
```

---

## 10. Implementation Checklist

### Phase 1: Core Implementation ✅
- [ ] Fix ExtractGeometry() to get actual centers
- [ ] Implement BiotSavartKernel (3×3 tensor)
- [ ] Implement FieldHMatrixCache structure
- [ ] Implement IsCacheValid()
- [ ] Implement BuildAndCacheFieldHMatrix()
- [ ] Implement EvaluateBatchCached()
- [ ] Unit tests for accuracy

### Phase 2: API Integration ✅
- [ ] Add use_hmatrix parameter to rad.Fld()
- [ ] Implement RadiaFld_Batch_HMatrix()
- [ ] Python bindings
- [ ] Modify RadiaFieldCF in radia_ngsolve.cpp
- [ ] Add use_hmatrix to RadiaField()

### Phase 3: Build & Test ✅
- [ ] Uncomment radhmat.cpp in CMakeLists.txt
- [ ] Build radia library
- [ ] Build radia_ngsolve
- [ ] Run test suite
- [ ] Fix any compilation/linking errors

### Phase 4: Benchmarking ✅
- [ ] Benchmark batch evaluation (no cache)
- [ ] Benchmark batch evaluation (with cache)
- [ ] Benchmark NGSolve integration
- [ ] Compare against targets
- [ ] Document results

---

## 11. Next Steps

1. **Review this simplified design** ← We are here
2. **Get approval to proceed**
3. **Start implementation**

### First Implementation Task

```
Phase 1, Day 1: Fix geometry extraction

File: src/core/radhmat.cpp:126-184
Goal: Get actual element center positions

Steps:
  1. Find how radTg3d stores center position
  2. Implement proper center extraction
  3. Test with simple magnet
  4. Verify positions are correct

Expected time: 2-3 hours
```

---

**Status**: ✅ Simplified design complete

**Key simplifications**:
- ✅ Batch evaluation only (no single-point)
- ✅ N ≥ 100 always (no small problem support)
- ✅ Explicit use_hmatrix=True API
- ✅ Mandatory caching for non-linear iterations

**Ready for implementation approval?**
