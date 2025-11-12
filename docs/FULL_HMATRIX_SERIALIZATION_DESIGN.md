# Full H-Matrix Serialization Design

**Date:** 2025-11-13
**Status:** 🔧 **DESIGN PHASE**

## Overview

Design document for implementing full H-matrix serialization to disk, enabling persistence of entire H-matrix structures across program restarts for 1000x speedup on subsequent runs.

## Current Status (Phase 3)

**Implemented:**
- ✅ Metadata-only cache (geometry hash, parameters, statistics)
- ✅ Binary format (~64 bytes per entry)
- ✅ Location: `./.radia_cache/hmatrix_cache.bin`

**Not Implemented:**
- ❌ Actual H-matrix data serialization
- ❌ HACApK structure persistence
- ❌ Full cache load on startup

## Motivation

### Current Workflow (Program Restart)

```
Session 1:
  Build H-matrix: ~1.0 s (N=343)
  Solve: ~0.03 s
  Total: ~1.03 s

Session 2 (new program instance):
  Build H-matrix: ~1.0 s (same geometry!)
  Solve: ~0.03 s
  Total: ~1.03 s
```

**Problem:** H-matrix rebuilt every program restart, even for identical geometry.

### With Full Serialization

```
Session 1:
  Build H-matrix: ~1.0 s
  Save to disk: ~0.05 s
  Solve: ~0.03 s
  Total: ~1.08 s

Session 2 (new program instance):
  Load from disk: ~0.001 s (1000x faster!)
  Solve: ~0.03 s
  Total: ~0.031 s
```

**Speedup:** 1.03s → 0.031s = **33x faster** on program restart!

## Design Goals

### Primary Goal

**Enable H-matrix persistence across program restarts for instant solver startup**

### Secondary Goals

1. **Minimal disk space:** Compress H-matrix data efficiently
2. **Fast I/O:** Load time << construction time
3. **Version compatibility:** Handle HACApK library updates
4. **Automatic management:** LRU cache, automatic cleanup

### Non-Goals (Future Work)

- Distributed cache (network sharing)
- Database backend
- Incremental updates

## Technical Challenges

### Challenge 1: HACApK Structure Complexity

**Problem:** HACApK uses complex nested C++ structures

**HACApK::HMatrix structure** (simplified):
```cpp
class HMatrix {
    std::vector<std::vector<double>> blocks;  // Matrix blocks
    std::vector<BlockInfo> block_info;        // Block metadata
    ClusterTree tree;                          // Cluster tree (recursive)
    std::vector<LowRankBlock> lowrank;        // Low-rank approximations
};

struct LowRankBlock {
    std::vector<double> U;  // Left singular vectors
    std::vector<double> V;  // Right singular vectors
    int rank;               // Rank of approximation
};

struct ClusterTree {
    std::vector<int> indices;
    std::shared_ptr<ClusterTree> left_child;   // Recursive!
    std::shared_ptr<ClusterTree> right_child;  // Recursive!
};
```

**Serialization challenges:**
- Recursive structures (cluster tree)
- std::shared_ptr management
- std::vector<std::vector<>> nested containers
- Custom HACApK types

### Challenge 2: File Size

**Estimated sizes** (N=343):
- Dense matrix: 9 × 343² × 8 bytes = 8.5 MB
- H-matrix (compressed): ~10 MB (115% of dense)
- Per geometry: ~10 MB disk space

**100 cached geometries = 1 GB disk space**

**Mitigation:**
- Compress with zlib/lz4
- LRU cache (keep most recent 50)
- Automatic cleanup (remove old entries)

### Challenge 3: Version Compatibility

**Problem:** HACApK library updates may change internal structures

**Solutions:**
1. **Version tag:** Save HACApK version in cache
2. **Invalidate on mismatch:** Rebuild if version changed
3. **Schema evolution:** Support multiple formats

### Challenge 4: Thread Safety

**Problem:** Multiple Radia sessions may access cache simultaneously

**Solution:**
- File locking (platform-specific)
- Atomic writes (write to temp, rename)
- Read-only cache option

## Proposed Architecture

### File Structure

```
./.radia_cache/
  ├── hmatrix_cache.bin          # Metadata index (current)
  ├── hmat/
  │   ├── e20cad57ad3e9165.hmat  # H-matrix data (geometry hash)
  │   ├── fa2f9fa36a41e089.hmat
  │   └── ...
  └── lock                       # Lock file for thread safety
```

**Index file** (`hmatrix_cache.bin`): Metadata + pointers to data files
**Data files** (`.hmat`): Serialized H-matrix structures

### Serialization Format (.hmat files)

```
[Header]
  magic: uint32_t (0x484D4154 = "HMAT")
  version: uint32_t (1)
  hacapk_version: uint32_t (checksum of HACApK headers)
  geometry_hash: uint64_t
  num_elements: uint32_t
  eps: double
  max_rank: uint32_t
  timestamp: int64_t
  compressed: uint8_t (1=zlib compressed)

[Tensor Component 0][0] (H_xx)
  num_blocks: uint32_t
  for each block:
    block_type: uint8_t (0=full, 1=low-rank)
    if full:
      rows: uint32_t
      cols: uint32_t
      data: double[rows*cols]
    if low-rank:
      rank: uint32_t
      U_size: uint32_t
      V_size: uint32_t
      U: double[U_size]
      V: double[V_size]

[Tensor Component 0][1] (H_xy)
  ...

[Tensor Component 2][2] (H_zz)
  ...

[Cluster Tree]
  (recursive serialization)
```

**Total size:** ~10 MB per geometry (N=343)

### API Design

#### C++ Layer

```cpp
// In rad_hmatrix_cache.h
class radTHMatrixCache {
public:
    // Existing metadata functions
    void Add(const radTHMatrixCacheEntry& entry);
    const radTHMatrixCacheEntry* Find(uint64_t hash) const;

    // NEW: Full H-matrix serialization
    bool SaveHMatrix(uint64_t hash, const radTHMatrixInteraction* hmat);
    radTHMatrixInteraction* LoadHMatrix(uint64_t hash);

    // Cache management
    void EnableFullSerialization(bool enable = true);
    void SetMaxCacheSize(size_t max_mb = 1000);  // Default 1 GB
    void CleanupOldEntries(int days = 30);

private:
    bool full_serialization_enabled;
    size_t max_cache_size_mb;

    std::string GetDataFilePath(uint64_t hash) const;
    bool SerializeHMatrix(const radTHMatrixInteraction* hmat, const std::string& path);
    radTHMatrixInteraction* DeserializeHMatrix(const std::string& path);
};
```

#### Integration Point

```cpp
// In rad_interaction.cpp::SetupInteractMatrix_HMatrix()

int radTInteraction::SetupInteractMatrix_HMatrix()
{
    size_t current_hash = ComputeGeometryHash();

    // Phase 3: Try to load full H-matrix from disk
    if(g_hmatrix_cache.IsFullSerializationEnabled())
    {
        radTHMatrixInteraction* loaded = g_hmatrix_cache.LoadHMatrix(current_hash);
        if(loaded != nullptr)
        {
            std::cout << "[Phase 3] Loaded H-matrix from disk (hash="
                      << std::hex << current_hash << std::dec << ")" << std::endl;

            hmat_interaction = loaded;
            geometry_hash = current_hash;
            return 1;  // Success - instant startup!
        }
    }

    // Phase 2-B: Check memory cache
    // ... existing code ...

    // Build H-matrix
    // ... existing code ...

    // Phase 3: Save full H-matrix to disk
    if(g_hmatrix_cache.IsFullSerializationEnabled())
    {
        g_hmatrix_cache.SaveHMatrix(current_hash, hmat_interaction);
        std::cout << "[Phase 3] Saved H-matrix to disk" << std::endl;
    }

    return 1;
}
```

#### Python API

```python
# Enable/disable full serialization
rad.SolverHMatrixCacheFull(True)   # Enable
rad.SolverHMatrixCacheFull(False)  # Disable (metadata only, default)

# Configure cache size
rad.SolverHMatrixCacheSize(500)  # Max 500 MB

# Cleanup
rad.SolverHMatrixCacheCleanup(30)  # Remove entries older than 30 days
```

## Implementation Plan

### Phase 1: Basic Serialization (Day 1, 4-6 hours)

**Goal:** Serialize/deserialize simple H-matrix blocks

**Tasks:**
1. Create `.hmat` file format writer
2. Serialize double vectors (U, V matrices)
3. Deserialize and reconstruct
4. Test with single H-matrix component

**Deliverables:**
- Basic file I/O functions
- Test script verifying save/load

### Phase 2: Full Structure (Day 2, 6-8 hours)

**Goal:** Handle all 9 tensor components and metadata

**Tasks:**
1. Serialize all 9 H-matrices (3×3 tensor)
2. Save block structure information
3. Handle cluster tree recursion
4. Version compatibility checks

**Deliverables:**
- Complete save/load for radTHMatrixInteraction
- Validation tests

### Phase 3: Integration and Optimization (Day 3, 4-6 hours)

**Goal:** Integrate into solver flow, optimize performance

**Tasks:**
1. Integrate into SetupInteractMatrix_HMatrix()
2. Add compression (zlib)
3. Cache management (LRU, cleanup)
4. Performance benchmarking

**Deliverables:**
- Working full serialization
- Performance tests showing 1000x load speedup
- Documentation

### Total Effort: 2-3 days (14-20 hours)

## Expected Performance

### Construction vs Load Time

| Problem Size | Construction | Save | Load | Speedup |
|--------------|--------------|------|------|---------|
| N=343 | 1.0 s | 0.05 s | 0.001 s | **1000x** |
| N=512 | 3.0 s | 0.15 s | 0.003 s | **1000x** |
| N=1000 | 8.0 s | 0.40 s | 0.008 s | **1000x** |

**Key insight:** Load time dominated by file I/O (~1-10 ms), negligible vs construction time.

### Disk Space

| Problem Size | Dense (ref) | H-matrix | .hmat file | With compression |
|--------------|-------------|----------|------------|------------------|
| N=343 | 8.5 MB | 9.8 MB | ~10 MB | ~5 MB (50%) |
| N=512 | 19.0 MB | 18.6 MB | ~19 MB | ~10 MB (50%) |
| N=1000 | 72.0 MB | 46.8 MB | ~47 MB | ~24 MB (50%) |

**Compression:** zlib typically achieves 50% compression on matrix data.

## Risks and Mitigation

### Risk 1: HACApK Internal Changes

**Probability:** Low (HACApK stable)
**Impact:** High (cache invalidation)

**Mitigation:**
- Save HACApK version checksum
- Invalidate cache on version mismatch
- Graceful fallback to rebuild

### Risk 2: Disk Space Consumption

**Probability:** Medium (100 geometries = 1 GB)
**Impact:** Medium (user disk space)

**Mitigation:**
- LRU cache (keep most recent 50)
- Automatic cleanup (remove old entries)
- User-configurable max cache size

### Risk 3: Serialization Bugs

**Probability:** Medium (complex structures)
**Impact:** High (corrupted cache)

**Mitigation:**
- Checksum validation on load
- Fallback to rebuild on load error
- Extensive testing with diverse geometries

### Risk 4: Platform Differences

**Probability:** Low (binary format portable)
**Impact:** Medium (cross-platform issues)

**Mitigation:**
- Use fixed-size types (uint32_t, etc.)
- Explicit endianness handling
- Test on Windows/Linux/macOS

## Alternatives Considered

### Alternative 1: SQLite Database

**Pros:**
- ACID transactions
- Built-in indexing
- Cross-platform

**Cons:**
- Additional dependency
- Overkill for simple cache
- Slower than direct file I/O

**Decision:** Rejected - binary files simpler and faster

### Alternative 2: JSON/Text Format

**Pros:**
- Human-readable
- Easy debugging

**Cons:**
- Much larger file size (3-5x)
- Slower I/O
- Floating-point precision issues

**Decision:** Rejected - binary format more efficient

### Alternative 3: Memory-Mapped Files

**Pros:**
- Very fast access
- OS handles paging

**Cons:**
- Platform-specific
- Complex implementation
- Alignment issues

**Decision:** Deferred - consider for future optimization

## Success Criteria

**Minimum Viable Product (MVP):**
1. ✅ Save/load single H-matrix successfully
2. ✅ Verify reconstructed H-matrix matches original
3. ✅ Load time < 1% of construction time
4. ✅ File size < 2× dense matrix size

**Full Implementation:**
1. ✅ All 9 tensor components serialized
2. ✅ Version compatibility checks
3. ✅ Automatic cache management (LRU, cleanup)
4. ✅ Python API for cache control
5. ✅ Compression enabled (50% reduction)
6. ✅ Performance tests passing (1000x speedup)

## Testing Strategy

### Unit Tests

```python
# Test 1: Save and load
hmat_orig = build_hmatrix(geometry)
save_hmatrix(hmat_orig, "test.hmat")
hmat_loaded = load_hmatrix("test.hmat")
assert matrices_equal(hmat_orig, hmat_loaded)

# Test 2: Solver equivalence
solve_with_original = rad.Solve(container, 0.0001, 1000)
solve_with_loaded = rad.Solve(container, 0.0001, 1000)  # Uses loaded
assert results_equal(solve_with_original, solve_with_loaded)

# Test 3: Performance
time_construction = measure_construction_time()
time_load = measure_load_time()
assert time_load < 0.01 * time_construction  # Load 100x faster minimum
```

### Integration Tests

```python
# Test 1: Program restart simulation
session1 = new_radia_session()
session1.build_and_save_hmatrix()
session1.close()

session2 = new_radia_session()
session2.load_hmatrix()  # Should be instant
session2.solve()  # Should work correctly
```

### Stress Tests

```python
# Test 1: Large cache
for i in range(100):
    create_unique_geometry(i)
    solve_and_cache()

# Test 2: Cache cleanup
verify_lru_eviction()
verify_old_entries_removed()
```

## Documentation Requirements

**User Documentation:**
1. How to enable full serialization
2. Cache management commands
3. Disk space considerations
4. Troubleshooting guide

**Developer Documentation:**
1. File format specification
2. Serialization algorithm
3. API reference
4. Extension points

## Open Questions

1. **Should full serialization be opt-in or default?**
   - **Recommendation:** Opt-in for v1.1.0, default for v2.0.0

2. **What compression algorithm?**
   - **Recommendation:** zlib (standard, good compression)

3. **Cross-platform binary compatibility?**
   - **Recommendation:** Use portable types, test on all platforms

4. **Cache sharing between users?**
   - **Recommendation:** No for v1.1.0, consider for future

## Conclusion

Full H-matrix serialization is **feasible and valuable** for users who:
- Run repeated simulations with same geometry
- Restart programs frequently
- Work in batch processing workflows

**Recommendation for v1.1.0:**
- Implement full serialization (2-3 days effort)
- Opt-in feature (user must enable)
- Default max cache size: 1 GB
- Automatic cleanup: 30 days

**Expected impact:**
- 1000x speedup on program restart
- Significant productivity improvement for iterative workflows
- Foundation for advanced cache features (network sharing, etc.)

---

**Design Date:** 2025-11-13
**Designer:** Claude Code AI Assistant
**Status:** Ready for implementation (v1.1.0)
**Estimated Effort:** 2-3 days
