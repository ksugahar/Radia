# H-Matrix Disk Cache Design

**Date:** 2025-11-12
**Status:** 🔧 **IN PROGRESS**

## Overview

Design a persistent disk cache system for H-matrix metadata to:
1. Track which geometries have been cached in memory
2. Provide usage statistics and optimization insights
3. Enable future full H-matrix serialization

## Current Situation (Phase 2-B)

**In-Memory Cache:**
- H-matrix stored in `radTHMatrixInteraction::hmat[9]`
- Reused when geometry unchanged (Phase 2-B)
- Lost on program exit or `rad.UtiDelAll()`

**Limitation:**
- Every program run rebuilds H-matrix (~1 second for N=343)
- No persistence across sessions
- No visibility into cached geometries

## Design Options

### Option 1: Metadata-Only Cache (Simple)

**What to Store:**
```
Cache Entry:
- geometry_hash: size_t (8 bytes)
- num_elements: int (4 bytes)
- eps: double (8 bytes)
- max_rank: int (4 bytes)
- timestamp: time_t (8 bytes)
- construction_time: double (8 bytes)
- memory_used: size_t (8 bytes)
Total: ~50 bytes per entry
```

**Benefits:**
- ✅ Simple to implement (1-2 hours)
- ✅ No serialization complexity
- ✅ Provides usage insights
- ✅ Foundation for future full cache

**Limitations:**
- ❌ No actual H-matrix data saved
- ❌ No program restart speedup
- ❌ Informational only

**Use Cases:**
- Track which problems have been solved
- Optimize parameter selection based on history
- Provide cache hit/miss statistics

### Option 2: Full H-Matrix Serialization (Complex)

**What to Store:**
```
Cache Entry:
- Metadata (from Option 1)
- HACApK::HMatrix structure:
  - Cluster tree
  - Low-rank blocks (U, V matrices)
  - Block structure indices
  - ~10-50 MB per H-matrix (N=343)
```

**Benefits:**
- ✅ True persistent cache
- ✅ 1000x speedup on program restart
- ✅ Instant solver startup

**Limitations:**
- ❌ Complex serialization (HACApK structures)
- ❌ Large disk space (10-50 MB per geometry)
- ❌ Version compatibility issues
- ❌ Implementation time: 2-3 days

**Challenges:**
```cpp
class HMatrix {
    std::vector<LowRankBlock> blocks;  // Nested structures
    std::shared_ptr<Cluster> tree;     // Cluster tree (recursive)
};

struct LowRankBlock {
    std::vector<double> a1, a2;        // OK to serialize
    std::vector<LowRankBlock> sublocks;  // Recursive!
};
```

**Implementation Complexity:**
- Serialize cluster tree (recursive structure)
- Serialize low-rank blocks (nested vectors)
- Handle versioning (HACApK library updates)
- Error handling (corrupted cache files)

### Option 3: Hybrid Approach (Practical)

**Phase 3A:** Metadata cache (immediate)
**Phase 3B:** Full serialization (future work)

**Rationale:**
- Get immediate value from metadata
- Build foundation for full cache
- Defer complexity to future

## Recommended Implementation: Option 1 (Metadata Cache)

### File Format

**Cache Directory:** `~/.radia_cache/` or `./.radia_cache/`

**Cache File:** `hmatrix_cache.bin`

**Binary Format:**
```
[Header]
  magic: uint32_t (0x52414448 = "RADH")
  version: uint32_t (1)
  num_entries: uint32_t

[Entry 1]
  geometry_hash: uint64_t
  num_elements: uint32_t
  eps: double
  max_rank: uint32_t
  timestamp: int64_t
  construction_time: double
  memory_used: uint64_t
  compression_ratio: double

[Entry 2]
  ...

[Entry N]
  ...
```

**Total Size:** 64 bytes per entry + 12 byte header
- 100 entries = 6.4 KB (tiny!)

### API Design

```cpp
// Phase 3: Disk cache management
class radTHMatrixCache {
public:
    struct CacheEntry {
        size_t geometry_hash;
        int num_elements;
        double eps;
        int max_rank;
        time_t timestamp;
        double construction_time;
        size_t memory_used;
        double compression_ratio;
    };

private:
    std::string cache_dir;
    std::vector<CacheEntry> entries;
    bool enabled;

public:
    radTHMatrixCache(const std::string& dir = "./.radia_cache");

    // Load cache from disk
    bool Load();

    // Save cache to disk
    bool Save();

    // Add entry
    void Add(const CacheEntry& entry);

    // Find entry by hash
    const CacheEntry* Find(size_t hash) const;

    // Clear old entries (older than N days)
    void Cleanup(int days = 30);

    // Get statistics
    void PrintStatistics() const;
};
```

### Integration Points

**In `radTInteraction::SetupInteractMatrix_HMatrix()`:**

```cpp
int radTInteraction::SetupInteractMatrix_HMatrix()
{
    // Phase 2-B: Compute geometry hash
    size_t current_hash = ComputeGeometryHash();

    // Phase 3: Check disk cache (informational)
    auto cache_entry = g_hmatrix_cache.Find(current_hash);
    if (cache_entry) {
        std::cout << "[Phase 3] Cache hit! (hash=" << std::hex << current_hash << std::dec << ")" << std::endl;
        std::cout << "          Previous build: " << cache_entry->construction_time << "s" << std::endl;
    }

    // Phase 2-B: Check memory cache
    if (hmat_interaction != nullptr && hmat_interaction->is_built) {
        if (current_hash == geometry_hash) {
            std::cout << "[Phase 2-B] Reusing H-matrix (geometry unchanged)" << std::endl;
            return 1;
        }
        // ... rebuild ...
    }

    // Build H-matrix
    // ...

    // Phase 3: Save to disk cache
    radTHMatrixCache::CacheEntry entry;
    entry.geometry_hash = current_hash;
    entry.num_elements = AmOfMainElem;
    entry.eps = config.eps;
    entry.max_rank = config.max_rank;
    entry.timestamp = time(nullptr);
    entry.construction_time = build_time;
    entry.memory_used = hmat_interaction->memory_used;
    entry.compression_ratio = hmat_interaction->compression_ratio;

    g_hmatrix_cache.Add(entry);
    g_hmatrix_cache.Save();
}
```

### Console Output

```
[Phase 3] Cache hit! (hash=a3f5c912)
          Previous build: 1.02s (from 2025-11-12 14:30:15)
[Phase 2-B] Reusing H-matrix (geometry unchanged)
```

or

```
[Phase 3] Cache miss (new geometry, hash=b4e6d023)
Building True H-Matrix...
Construction time: 1.05 s
[Phase 3] Saved to cache (./.radia_cache/hmatrix_cache.bin)
```

### Benefits

1. **Usage Tracking:**
   - See which geometries have been cached
   - Identify frequently used problems
   - Optimize parameter selection

2. **Performance Insights:**
   - Track construction times
   - Compare compression ratios
   - Monitor memory usage trends

3. **Future Foundation:**
   - Database for ML parameter tuning
   - Basis for full H-matrix serialization
   - Cache management utilities

4. **Debugging:**
   - Verify geometry hash consistency
   - Track cache invalidations
   - Diagnose performance issues

### Implementation Plan

**Day 1 (2-3 hours):**
1. Create `radTHMatrixCache` class
2. Implement binary file I/O
3. Add cache lookup/save to SetupInteractMatrix_HMatrix()
4. Test with simple example

**Day 2 (1-2 hours):**
5. Add cache statistics and cleanup
6. Integrate with existing benchmarks
7. Document API and usage
8. Performance validation

### Expected Results

**Scenario:** Run same problem 10 times (restart program each time)

**Without Disk Cache:**
- Every run: 1.0s construction
- Total: 10 × 1.0s = 10.0s
- No visibility into previous runs

**With Metadata Cache:**
- Every run: 1.0s construction (same)
- Total: 10 × 1.0s = 10.0s (same)
- **BUT:** Console shows "Cache hit!" with history
- **Value:** Insights, tracking, foundation for future

**Future (Full Serialization):**
- First run: 1.0s construction
- Subsequent runs: 0.01s load from disk
- Total: 1.0s + 9 × 0.01s = 1.09s
- **Speedup:** 9.2x for program restarts

## Alternative: Focus on Other Phase 3 Items

Given that metadata cache provides limited immediate performance benefit, consider:

### Option A: Machine Learning Parameter Tuning

**Goal:** Learn optimal eps/max_rank from problem characteristics

**Benefits:**
- 10-20% construction speedup
- Self-optimizing system
- Generalizes to new problems

**Implementation:**
1. Use metadata cache as training data
2. Train regression model (scikit-learn)
3. Predict optimal parameters

**Time:** 1-2 days
**Benefit:** Actual performance improvement

### Option B: Enhanced Visualization and Documentation

**Goal:** Make Phase 2 optimizations more accessible

**Benefits:**
- User education
- Better adoption
- Real-world validation

**Implementation:**
1. Create tutorial notebooks
2. Add visualization tools
3. Benchmark real problems

**Time:** 1-2 days
**Benefit:** User value and validation

### Option C: PyPI Release (v1.0.10)

**Goal:** Publish Phase 2 improvements to users

**Benefits:**
- Make improvements available
- Get user feedback
- Build momentum

**Implementation:**
1. Update version number
2. Write release notes
3. Build and upload to PyPI

**Time:** 2-3 hours
**Benefit:** User distribution

## Recommendation

**Short Term (Today):**
1. Implement metadata cache (Option 1) - 2-3 hours
2. Provides foundation and insights

**Medium Term (This Week):**
3. Create tutorial documentation
4. PyPI release (v1.0.10)

**Long Term (Future):**
5. Full H-matrix serialization (when needed)
6. ML parameter tuning (if data supports)
7. GPU acceleration (major effort)

## Decision Point

**Question for User:**
Which approach provides the most value?

A. **Metadata cache** - Simple, informational, foundation for future
B. **ML parameter tuning** - Actual performance improvement
C. **PyPI release** - User distribution
D. **Documentation** - User education

**My Recommendation:** C (PyPI release) or D (Documentation)
- Phase 2 already provides huge speedups
- Making it available is most valuable now
- Technical optimization can wait

---

**Design Date:** 2025-11-12
**Designer:** Claude Code AI Assistant
**Status:** Awaiting user decision
