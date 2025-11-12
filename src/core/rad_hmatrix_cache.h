/*-------------------------------------------------------------------------
*
* File name:      rad_hmatrix_cache.h
*
* Project:        RADIA
*
* Description:    Disk cache for H-matrix metadata (Phase 3)
*                 Tracks cached geometries for statistics and optimization
*
* Author(s):      Claude Code AI Assistant
*
* First release:  2025-11-12
*
* Copyright (C):  2025
*
*-------------------------------------------------------------------------*/

#ifndef __RAD_HMATRIX_CACHE_H
#define __RAD_HMATRIX_CACHE_H

#include <string>
#include <vector>
#include <ctime>
#include <cstdint>

//-------------------------------------------------------------------------
// H-Matrix Cache Entry (metadata only)
//-------------------------------------------------------------------------

struct radTHMatrixCacheEntry
{
	uint64_t geometry_hash;      // Geometry hash from Phase 2-B
	uint32_t num_elements;       // Number of elements
	double eps;                  // ACA tolerance used
	uint32_t max_rank;           // Maximum rank used
	int64_t timestamp;           // When cached (Unix time)
	double construction_time;    // Build time (seconds)
	uint64_t memory_used;        // Memory usage (bytes)
	double compression_ratio;    // Compression ratio

	radTHMatrixCacheEntry();
};

//-------------------------------------------------------------------------
// H-Matrix Disk Cache (Phase 3)
//
// Purpose: Track H-matrix constructions for statistics and optimization
//
// Features:
//   - Persistent metadata storage
//   - Usage tracking and statistics
//   - Foundation for ML parameter tuning
//   - Foundation for future full H-matrix serialization
//
// File Format: Binary, ~64 bytes per entry
// Location: ./.radia_cache/hmatrix_cache.bin
//-------------------------------------------------------------------------

class radTHMatrixCache
{
private:
	std::string cache_dir;
	std::string cache_file;
	std::vector<radTHMatrixCacheEntry> entries;
	bool enabled;
	bool dirty;  // Need to save

public:
	radTHMatrixCache(const std::string& dir = "./.radia_cache");
	~radTHMatrixCache();

	// Enable/disable cache
	void Enable(bool enable = true);
	bool IsEnabled() const { return enabled; }

	// Load cache from disk
	bool Load();

	// Save cache to disk
	bool Save();

	// Add entry
	void Add(const radTHMatrixCacheEntry& entry);

	// Find entry by hash
	const radTHMatrixCacheEntry* Find(uint64_t hash) const;

	// Clear old entries (older than N days)
	void Cleanup(int days = 30);

	// Statistics
	void PrintStatistics() const;
	size_t GetNumEntries() const { return entries.size(); }

	// Get cache directory
	const std::string& GetCacheDir() const { return cache_dir; }

private:
	// Ensure cache directory exists
	bool EnsureCacheDirectory();
};

// Global cache instance
extern radTHMatrixCache g_hmatrix_cache;

#endif
