# radia_ngsolve FMM Acceleration Design

## Overview

This document describes the design for FMM (Fast Multipole Method) acceleration of `GridFunction.Set()` operations in radia_ngsolve.

## Design Goals

1. **Minimal Python dependency**: All heavy computation in C++
2. **NGSolve native integration**: Hook into `GridFunction.Set()` at C++ level
3. **Future GPU support**: Architecture that enables CUDA/OpenCL acceleration
4. **Particle tracking**: Same FMM infrastructure for particle beam dynamics

## Current Problem

### Current Implementation Flow

```
GridFunction.Set(RadiaField(obj, 'b'))
    |
    v
For each element in mesh:
    For each integration point:
        RadiaFieldCF::Evaluate(point)
            |
            v
        rad.Fld(obj, 'b', point)  <-- O(N_sources) per call
```

**Problem**: Each `rad.Fld()` call is O(N_sources), resulting in O(N_targets * N_sources) total complexity.

For large problems:
- N_sources = 1000 Radia elements (dipoles)
- N_targets = 100,000 mesh DOFs
- Total: 100,000,000 dipole evaluations

### Failed Approach: Point Cache

Previous approach tried Python-side caching:
```python
# Collect DOF points from mesh
points = collect_mesh_dof_points(mesh, fes)
cf.PrepareCache(points)  # Batch call rad.Fld(obj, 'b', points_list)
gf.Set(cf)  # Hits cache
```

**Why it doesn't work well**:
1. NGSolve's `Set()` evaluates at **integration points**, not DOF points
2. Integration point locations depend on element shape and quadrature rule
3. Cache misses when integration points don't match pre-cached DOF points
4. Python-side point collection is slow and duplicates work

## Proposed Solution: C++ FMM Integration

### Key Insight

The FMM acceleration should happen at the **batch evaluation** level, not per-point:

1. **Collect all target points** from the GridFunction.Set() operation
2. **Extract all source dipoles** from Radia object
3. **Call FMM once** for all sources -> all targets
4. **Return cached results** during subsequent Evaluate() calls

### Architecture

```
                         radia_ngsolve.cpp
                               |
                               v
    +--------------------------------------------------+
    |              RadiaFieldCF                        |
    |  +--------------------------------------------+  |
    |  |  Evaluate(mir) - batch evaluation          |  |
    |  |    1. Collect all points from mir          |  |
    |  |    2. Check if FMM cache valid             |  |
    |  |    3. If not, call FMM computation         |  |
    |  |    4. Return cached results                |  |
    |  +--------------------------------------------+  |
    +--------------------------------------------------+
                               |
                               v
    +--------------------------------------------------+
    |              RadFldFMM (new C++ function)        |
    |  +--------------------------------------------+  |
    |  |  1. Extract dipoles from Radia container   |  |
    |  |  2. Call RadExaFMM::ComputeDipoleField()   |  |
    |  |  3. Return B = mu0 * H from dipole field   |  |
    |  +--------------------------------------------+  |
    +--------------------------------------------------+
                               |
                               v
    +--------------------------------------------------+
    |              RadExaFMM (existing)                |
    |  +--------------------------------------------+  |
    |  |  ComputeDipoleField() or                   |  |
    |  |  ComputeDipoleFieldDirect() (O(N*M))       |  |
    |  +--------------------------------------------+  |
    +--------------------------------------------------+
```

### API Design

#### New C API: `RadFldBatchFMM`

```cpp
// In radentry.h
/**
 * Compute B-field at multiple target points using FMM acceleration.
 *
 * This is the optimized version for large-scale field computation
 * (e.g., NGSolve GridFunction interpolation).
 *
 * @param B_out      Output B-field [Bx1,By1,Bz1, Bx2,By2,Bz2, ...] (n_targets*3)
 * @param n_targets  Number of target points
 * @param targets    Target positions [x1,y1,z1, x2,y2,z2, ...] (n_targets*3)
 * @param obj_handle Radia container handle
 * @param fmm_eps    FMM tolerance (1e-6 typical), 0 = use direct
 * @return           Error code (0 = success)
 */
EXP int CALL RadFldBatchFMM(
    double* B_out,
    int n_targets,
    double* targets,
    int obj_handle,
    double fmm_eps
);
```

#### C++ Implementation in radentry.cpp

```cpp
int RadFldBatchFMM(double* B_out, int n_targets, double* targets,
                   int obj_handle, double fmm_eps)
{
    try {
        // 1. Get container from handle
        radThg& hg = rad.RetrieveHandle(obj_handle);
        radTg3d* g3d = dynamic_cast<radTg3d*>(hg.rep);

        // 2. Collect all dipoles from container (recursive)
        std::vector<double> sources;   // [x1,y1,z1, ...]
        std::vector<double> dipoles;   // [mx1,my1,mz1, ...]
        CollectDipoles(g3d, sources, dipoles);

        int64_t n_sources = sources.size() / 3;

        // 3. Call FMM
        RadExaFMM::FMMResult result = RadExaFMM::ComputeDipoleField(
            fmm_eps, sources.data(), dipoles.data(), n_sources,
            targets, n_targets
        );

        // 4. Convert H -> B (B = mu0 * H)
        const double MU0 = 4.0 * M_PI * 1e-7;
        for (int i = 0; i < n_targets; i++) {
            B_out[i*3 + 0] = MU0 * result.gradx[i];
            B_out[i*3 + 1] = MU0 * result.grady[i];
            B_out[i*3 + 2] = MU0 * result.gradz[i];
        }

        return 0;
    } catch (...) {
        return -1;
    }
}

// Helper: Recursively collect dipoles from container
void CollectDipoles(radTg3d* obj, std::vector<double>& sources,
                    std::vector<double>& dipoles)
{
    if (radTGroup* grp = dynamic_cast<radTGroup*>(obj)) {
        // Container: recurse into children
        for (auto& ch : grp->GroupMapOfHandlers) {
            radTg3d* child = dynamic_cast<radTg3d*>(ch.second.rep);
            if (child) CollectDipoles(child, sources, dipoles);
        }
    }
    else if (radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(obj)) {
        // Polyhedron element: get centroid and magnetization
        TVector3d center = poly->GetCenterPoint();
        TVector3d M = poly->GetMagnetization();
        double V = poly->GetVolume();

        // Dipole moment: m = M * V
        sources.push_back(center.x);
        sources.push_back(center.y);
        sources.push_back(center.z);

        dipoles.push_back(M.x * V);
        dipoles.push_back(M.y * V);
        dipoles.push_back(M.z * V);
    }
    // Add other element types (RecMag, etc.) as needed
}
```

### Integration with radia_ngsolve.cpp

#### Option A: Python-side Control (Simpler)

```python
# User explicitly enables FMM mode
cf = RadiaField(obj, 'b', fmm_eps=1e-6)
gf.Set(cf)
```

In C++:
```cpp
// RadiaFieldCF with FMM mode
class RadiaFieldCF : public CoefficientFunction
{
    double fmm_eps_;  // 0 = disabled, >0 = FMM tolerance
    mutable std::vector<double> fmm_cache_;  // Cached results
    mutable bool cache_valid_;

    // Called once for entire Set() operation
    void PrepareFMMCache(const std::vector<std::array<double,3>>& all_points) {
        if (fmm_eps_ <= 0) return;

        int n = all_points.size();
        std::vector<double> targets(n * 3);
        for (int i = 0; i < n; i++) {
            targets[i*3+0] = all_points[i][0] * coord_scale_;
            targets[i*3+1] = all_points[i][1] * coord_scale_;
            targets[i*3+2] = all_points[i][2] * coord_scale_;
        }

        fmm_cache_.resize(n * 3);
        RadFldBatchFMM(fmm_cache_.data(), n, targets.data(),
                       radia_obj, fmm_eps_);
        cache_valid_ = true;
    }
};
```

#### Option B: Automatic Detection (Complex but user-friendly)

NGSolve's `GridFunction::Set()` calls:
1. `Set(cf)` -> internal loop
2. For each element: `cf.Evaluate(mir)` (batch for integration rule)

We can intercept at the `Evaluate(mir)` level:

```cpp
void Evaluate(const BaseMappedIntegrationRule& mir,
              BareSliceMatrix<> result) const override
{
    // Collect points from this integration rule
    size_t npts = mir.Size();

    if (fmm_eps_ > 0 && !cache_valid_) {
        // First call: trigger FMM computation for ALL mesh points
        // This requires knowing all points ahead of time...
        // -> Difficult without modifying NGSolve
    }
}
```

**Conclusion**: Option A (user-controlled) is simpler and more robust.

### Workflow for GridFunction.Set()

#### Recommended Usage Pattern

```python
import radia as rad
import radia_ngsolve
from ngsolve import *

# Create Radia model (many elements)
obj = rad.ObjCnt([...])  # 1000+ elements
rad.Solve(obj, 0.0001, 1000, 1)

# Create NGSolve mesh
mesh = Mesh(...)

# Method 1: Standard (slow for large problems)
cf = radia_ngsolve.RadiaField(obj, 'b')
gf = GridFunction(HDiv(mesh, order=2))
gf.Set(cf)  # O(N_targets * N_sources)

# Method 2: FMM-accelerated (proposed)
# Step 1: Collect all DOF points from mesh
points = collect_interpolation_points(mesh, fes)

# Step 2: Batch FMM computation
B_values = rad.FldBatchFMM(obj, points, eps=1e-6)

# Step 3: Create GridFunction from values
gf_fmm = interpolate_to_gridfunction(mesh, fes, points, B_values)
```

### Alternative: Extend rad.Fld() with FMM

Rather than a separate API, extend existing `rad.Fld()`:

```python
# Current API (calls directly per point or batch)
B = rad.Fld(obj, 'b', [x, y, z])
B_list = rad.Fld(obj, 'b', [[x1,y1,z1], [x2,y2,z2], ...])

# Extended API with FMM option
B_list = rad.Fld(obj, 'b', points, method='fmm', eps=1e-6)
```

This integrates naturally with the existing batch evaluation in radia_ngsolve.cpp.

## Implementation Plan

### Phase 1: Core FMM Integration (C++)

1. **Add `CollectDipoles()` helper** in `rad_application.cpp`
   - Recursively traverse container
   - Extract center, magnetization, volume for each element

2. **Add `RadFldBatchFMM()`** in `radentry.cpp`
   - Call CollectDipoles
   - Call RadExaFMM::ComputeDipoleField
   - Return B = mu0 * H

3. **Python binding** in `radpy_pyapi.cpp`
   - `rad.FldBatchFMM(obj, points, eps=1e-6)` -> returns list of [Bx,By,Bz]

### Phase 2: radia_ngsolve Integration

1. **Add FMM parameters** to RadiaFieldCF constructor
   - `fmm_eps`: FMM tolerance (0 = disabled)
   - `fmm_threshold`: Minimum points to use FMM (e.g., 1000)

2. **Modify batch Evaluate()** to use FMM cache
   - If FMM enabled and cache not valid:
     - Collect all points from mir
     - Call RadFldBatchFMM
     - Store in cache
   - Return cached values

3. **Add `PrepareFMM()` method** for explicit pre-computation
   - User can call before `gf.Set()` if they know the mesh

### Phase 3: NGSolve-side Helper Functions

1. **Python helper**: `collect_interpolation_points(mesh, fes)`
   - Returns all integration points that `Set()` will evaluate

2. **Python helper**: `interpolate_to_gridfunction(mesh, fes, points, values)`
   - Creates GridFunction from point-value pairs

## Limitations and Considerations

### FMM Approximation

- FMM computes **dipole field approximation**
- For MSC (surface charge) elements, dipole approximation is less accurate near sources
- Far-field: FMM is highly accurate
- Near-field: Need direct computation or higher-order multipoles

### When to Use FMM

| Problem Size | Recommended Method |
|--------------|-------------------|
| < 1000 targets | Direct (rad.Fld batch) |
| 1000-10000 targets | FMM with eps=1e-4 |
| > 10000 targets | FMM with eps=1e-6 |

### Mixed Source Types

Current Radia has multiple element types:
- `radTPolyhedron` (MSC hexahedra/tetrahedra)
- `radTRecMag` (rectangular magnets)
- `radTArcCur` (arc currents)
- etc.

FMM implementation needs to handle all source types, extracting:
- Position (centroid)
- Equivalent dipole moment (M * V or I * A for coils)

## Files to Modify

| File | Changes |
|------|---------|
| `src/lib/radentry.h` | Add `RadFldBatchFMM` declaration |
| `src/lib/radentry.cpp` | Implement `RadFldBatchFMM` |
| `src/core/rad_application.h` | Add `CollectDipoles` declaration |
| `src/core/rad_application.cpp` | Implement `CollectDipoles` |
| `src/radia/radpy_pyapi.cpp` | Add Python binding |
| `src/radia/radia_ngsolve.cpp` | Add FMM mode to RadiaFieldCF |

---

## Revised Design: Pure C++ FMM Integration

### Design Philosophy

**Key Principle**: Minimize Python dependency - Python only serves as configuration layer.

```
Python layer (thin):
  - RadiaField(obj, 'b', fmm_eps=1e-6)  # Just pass parameters
  - gf.Set(cf)                          # Triggers C++ code path

C++ layer (all computation):
  - NGSolve's Set() calls Evaluate()
  - FMM dispatch at first batch call
  - All dipole extraction, FMM, results in C++
```

### C++ FMM Dispatcher Architecture

```cpp
class RadiaFieldCF : public CoefficientFunction
{
private:
    // FMM acceleration
    double fmm_eps_;                    // 0 = disabled, >0 = FMM tolerance
    mutable RadFMMState* fmm_state_;    // FMM computation state (lazy init)

    // Source data (extracted once from Radia)
    mutable std::vector<double> source_pos_;      // [x,y,z, ...]
    mutable std::vector<double> source_dipoles_;  // [mx,my,mz, ...]
    mutable bool sources_extracted_;

public:
    // Batch evaluation - key entry point for Set()
    void Evaluate(const BaseMappedIntegrationRule& mir,
                  BareSliceMatrix<> result) const override
    {
        size_t npts = mir.Size();

        // 1. Extract sources from Radia (first time only)
        if (!sources_extracted_) {
            ExtractSources();  // Pure C++ - no Python
        }

        // 2. Collect target points
        std::vector<double> targets(npts * 3);
        for (size_t i = 0; i < npts; i++) {
            auto pnt = mir[i].GetPoint();
            targets[i*3+0] = pnt[0] * coord_scale_;
            targets[i*3+1] = pnt[1] * coord_scale_;
            targets[i*3+2] = pnt[2] * coord_scale_;
        }

        // 3. FMM dispatch (pure C++)
        std::vector<double> B_out(npts * 3);
        if (fmm_eps_ > 0) {
            // FMM path
            RadExaFMM::FMMResult res = RadExaFMM::ComputeDipoleField(
                fmm_eps_,
                source_pos_.data(), source_dipoles_.data(),
                source_pos_.size() / 3,
                targets.data(), npts
            );
            // H -> B
            for (size_t i = 0; i < npts; i++) {
                B_out[i*3+0] = MU0 * res.gradx[i];
                B_out[i*3+1] = MU0 * res.grady[i];
                B_out[i*3+2] = MU0 * res.gradz[i];
            }
        } else {
            // Direct path (existing Fld batch)
            RadFldBatch(B_out.data(), nullptr, npts, targets.data(),
                        radia_obj, 0);
        }

        // 4. Copy to result matrix
        for (size_t i = 0; i < npts; i++) {
            result(i, 0) = B_out[i*3+0];
            result(i, 1) = B_out[i*3+1];
            result(i, 2) = B_out[i*3+2];
        }
    }

private:
    // Extract dipole sources from Radia container (pure C++)
    void ExtractSources() const {
        source_pos_.clear();
        source_dipoles_.clear();

        // Direct C++ access to Radia internals
        // (Not through Python rad.Fld)
        radThg& hg = rad.RetrieveHandle(radia_obj);
        radTg3d* g3d = dynamic_cast<radTg3d*>(hg.rep);
        if (g3d) {
            CollectDipolesRecursive(g3d, source_pos_, source_dipoles_);
        }
        sources_extracted_ = true;
    }
};
```

### FMM State Management

```cpp
// In rad_fmm_state.h
class RadFMMState {
public:
    // Source data
    std::vector<double> sources;      // Position [x,y,z,...]
    std::vector<double> dipoles;      // Dipole moment [mx,my,mz,...]
    int64_t n_sources;

    // Pre-built FMM tree (if using ExaFMM)
    void* fmm_tree;                   // Opaque pointer to tree structure

    // Configuration
    double eps;
    int expansion_order;
    int ncrit;

    // Build FMM tree from sources
    void BuildTree();

    // Evaluate at targets (reuses tree)
    void Evaluate(const double* targets, int64_t n_targets,
                  double* B_out);

    // GPU variants (future)
    void BuildTreeGPU();
    void EvaluateGPU(const double* targets_gpu, int64_t n_targets,
                     double* B_out_gpu);
};
```

### GPU Architecture (Future)

```
                    RadFMMState
                         |
          +--------------+--------------+
          |              |              |
       CPU Direct    CPU FMM       GPU FMM
          |              |              |
     OpenMP loop    ExaFMM-t      CUDA/OpenCL
          |              |              |
          +------+-------+------+-------+
                 |               |
              B_out           B_out_gpu
                 |               |
                 +-------+-------+
                         |
                    GridFunction
```

```cpp
// Future GPU interface
class RadFMMStateGPU : public RadFMMState {
public:
    // Device memory
    double* d_sources;
    double* d_dipoles;
    double* d_targets;
    double* d_B_out;

    // Particle tracking: evaluate at many particle positions
    void EvaluateParticles(int n_particles,
                           const double* particle_pos,
                           double* B_at_particles);

    // Interpolation: field on regular grid for GPU particle pusher
    void EvaluateGrid(int nx, int ny, int nz,
                      double x0, double y0, double z0,
                      double dx, double dy, double dz,
                      double* B_grid_gpu);
};
```

### Particle Tracking Use Case

```cpp
// GPU particle tracking with FMM-accelerated field
class RadiaParticleTracker {
public:
    RadFMMStateGPU fmm_state;

    void Step(ParticleArray& particles, double dt) {
        // 1. Copy particle positions to GPU
        cudaMemcpy(d_pos, particles.pos, ...);

        // 2. Evaluate B-field at all particle positions (FMM)
        fmm_state.EvaluateParticles(particles.n, d_pos, d_B);

        // 3. Boris push (on GPU)
        BorisPushKernel<<<...>>>(particles.n, d_pos, d_vel, d_B, dt);

        // 4. Copy back to CPU (or keep on GPU)
        cudaMemcpy(particles.pos, d_pos, ...);
    }
};
```

### Summary: Implementation Priorities

**Phase 1: Pure C++ FMM in radia_ngsolve** (This PR)
1. Add `ExtractSources()` - C++ dipole extraction
2. Use existing `RadExaFMM::ComputeDipoleField()` in Evaluate()
3. Remove Python rad.Fld() calls from batch path
4. Test with large NGSolve meshes

**Phase 2: FMM Tree Reuse**
1. Add `RadFMMState` class
2. Build tree once, evaluate multiple times
3. Useful for time-dependent problems

**Phase 3: GPU Support**
1. CUDA implementation of direct dipole sum
2. GPU FMM (e.g., ExaFMM-t GPU backend)
3. GPU particle pusher integration

## References

1. ExaFMM-t library: https://github.com/exafmm/exafmm-t
2. FMM3D library: https://fmm3d.readthedocs.io/
3. NGSolve CoefficientFunction: https://docu.ngsolve.org/latest/
4. Boris particle pusher: Birdsall & Langdon, "Plasma Physics via Computer Simulation"
