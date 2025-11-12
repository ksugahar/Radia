/*-------------------------------------------------------------------------
*
* File name:      rad_hmatrix_cache.cpp
*
* Project:        RADIA
*
* Description:    Disk cache for H-matrix metadata (Phase 3)
*
*-------------------------------------------------------------------------*/

#include "rad_hmatrix_cache.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sys/stat.h>

#ifdef _WIN32
#include <direct.h>
#define mkdir(dir, mode) _mkdir(dir)
#else
#include <sys/stat.h>
#endif

//-------------------------------------------------------------------------
// Global cache instance
//-------------------------------------------------------------------------

radTHMatrixCache g_hmatrix_cache;

//-------------------------------------------------------------------------
// Cache Entry
//-------------------------------------------------------------------------

radTHMatrixCacheEntry::radTHMatrixCacheEntry()
	: geometry_hash(0), num_elements(0), eps(0.0), max_rank(0),
	  timestamp(0), construction_time(0.0), memory_used(0), compression_ratio(0.0)
{
}

//-------------------------------------------------------------------------
// Cache Implementation
//-------------------------------------------------------------------------

radTHMatrixCache::radTHMatrixCache(const std::string& dir)
	: cache_dir(dir), enabled(true), dirty(false)
{
	cache_file = cache_dir + "/hmatrix_cache.bin";
}

radTHMatrixCache::~radTHMatrixCache()
{
	if(dirty && enabled)
	{
		Save();
	}
}

void radTHMatrixCache::Enable(bool enable)
{
	enabled = enable;
}

bool radTHMatrixCache::EnsureCacheDirectory()
{
	struct stat info;

	if(stat(cache_dir.c_str(), &info) != 0)
	{
		// Directory doesn't exist, create it
		if(mkdir(cache_dir.c_str(), 0755) != 0)
		{
			std::cerr << "[Phase 3] Warning: Failed to create cache directory: "
			          << cache_dir << std::endl;
			return false;
		}
	}
	else if(!(info.st_mode & S_IFDIR))
	{
		std::cerr << "[Phase 3] Warning: Cache path exists but is not a directory: "
		          << cache_dir << std::endl;
		return false;
	}

	return true;
}

bool radTHMatrixCache::Load()
{
	if(!enabled) return false;

	std::ifstream file(cache_file, std::ios::binary);
	if(!file.is_open())
	{
		// Cache file doesn't exist yet (first run)
		return false;
	}

	// Read header
	uint32_t magic, version, num_entries;

	file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
	file.read(reinterpret_cast<char*>(&version), sizeof(version));
	file.read(reinterpret_cast<char*>(&num_entries), sizeof(num_entries));

	// Validate header
	const uint32_t MAGIC = 0x52414448;  // "RADH"
	if(magic != MAGIC || version != 1)
	{
		std::cerr << "[Phase 3] Warning: Invalid cache file format" << std::endl;
		return false;
	}

	// Read entries
	entries.clear();
	entries.reserve(num_entries);

	for(uint32_t i = 0; i < num_entries; i++)
	{
		radTHMatrixCacheEntry entry;

		file.read(reinterpret_cast<char*>(&entry.geometry_hash), sizeof(entry.geometry_hash));
		file.read(reinterpret_cast<char*>(&entry.num_elements), sizeof(entry.num_elements));
		file.read(reinterpret_cast<char*>(&entry.eps), sizeof(entry.eps));
		file.read(reinterpret_cast<char*>(&entry.max_rank), sizeof(entry.max_rank));
		file.read(reinterpret_cast<char*>(&entry.timestamp), sizeof(entry.timestamp));
		file.read(reinterpret_cast<char*>(&entry.construction_time), sizeof(entry.construction_time));
		file.read(reinterpret_cast<char*>(&entry.memory_used), sizeof(entry.memory_used));
		file.read(reinterpret_cast<char*>(&entry.compression_ratio), sizeof(entry.compression_ratio));

		if(!file.good())
		{
			std::cerr << "[Phase 3] Warning: Error reading cache entry " << i << std::endl;
			break;
		}

		entries.push_back(entry);
	}

	file.close();

	std::cout << "[Phase 3] Loaded cache: " << entries.size() << " entries from "
	          << cache_file << std::endl;

	dirty = false;
	return true;
}

bool radTHMatrixCache::Save()
{
	if(!enabled || !dirty) return false;

	if(!EnsureCacheDirectory())
	{
		return false;
	}

	std::ofstream file(cache_file, std::ios::binary | std::ios::trunc);
	if(!file.is_open())
	{
		std::cerr << "[Phase 3] Warning: Failed to open cache file for writing: "
		          << cache_file << std::endl;
		return false;
	}

	// Write header
	const uint32_t MAGIC = 0x52414448;  // "RADH"
	const uint32_t VERSION = 1;
	uint32_t num_entries = static_cast<uint32_t>(entries.size());

	file.write(reinterpret_cast<const char*>(&MAGIC), sizeof(MAGIC));
	file.write(reinterpret_cast<const char*>(&VERSION), sizeof(VERSION));
	file.write(reinterpret_cast<const char*>(&num_entries), sizeof(num_entries));

	// Write entries
	for(const auto& entry : entries)
	{
		file.write(reinterpret_cast<const char*>(&entry.geometry_hash), sizeof(entry.geometry_hash));
		file.write(reinterpret_cast<const char*>(&entry.num_elements), sizeof(entry.num_elements));
		file.write(reinterpret_cast<const char*>(&entry.eps), sizeof(entry.eps));
		file.write(reinterpret_cast<const char*>(&entry.max_rank), sizeof(entry.max_rank));
		file.write(reinterpret_cast<const char*>(&entry.timestamp), sizeof(entry.timestamp));
		file.write(reinterpret_cast<const char*>(&entry.construction_time), sizeof(entry.construction_time));
		file.write(reinterpret_cast<const char*>(&entry.memory_used), sizeof(entry.memory_used));
		file.write(reinterpret_cast<const char*>(&entry.compression_ratio), sizeof(entry.compression_ratio));
	}

	file.close();

	dirty = false;
	return true;
}

void radTHMatrixCache::Add(const radTHMatrixCacheEntry& entry)
{
	if(!enabled) return;

	// Check if entry already exists
	auto it = std::find_if(entries.begin(), entries.end(),
		[&entry](const radTHMatrixCacheEntry& e) {
			return e.geometry_hash == entry.geometry_hash &&
			       e.eps == entry.eps &&
			       e.max_rank == entry.max_rank;
		});

	if(it != entries.end())
	{
		// Update existing entry
		*it = entry;
	}
	else
	{
		// Add new entry
		entries.push_back(entry);
	}

	dirty = true;
}

const radTHMatrixCacheEntry* radTHMatrixCache::Find(uint64_t hash) const
{
	if(!enabled) return nullptr;

	auto it = std::find_if(entries.begin(), entries.end(),
		[hash](const radTHMatrixCacheEntry& e) {
			return e.geometry_hash == hash;
		});

	return (it != entries.end()) ? &(*it) : nullptr;
}

void radTHMatrixCache::Cleanup(int days)
{
	if(!enabled) return;

	int64_t cutoff_time = std::time(nullptr) - (days * 24 * 3600);

	auto new_end = std::remove_if(entries.begin(), entries.end(),
		[cutoff_time](const radTHMatrixCacheEntry& e) {
			return e.timestamp < cutoff_time;
		});

	if(new_end != entries.end())
	{
		size_t removed = std::distance(new_end, entries.end());
		entries.erase(new_end, entries.end());

		std::cout << "[Phase 3] Cleaned up " << removed << " old cache entries" << std::endl;

		dirty = true;
	}
}

void radTHMatrixCache::PrintStatistics() const
{
	if(!enabled || entries.empty())
	{
		std::cout << "[Phase 3] Cache: No entries" << std::endl;
		return;
	}

	std::cout << "\n========================================" << std::endl;
	std::cout << "H-Matrix Cache Statistics" << std::endl;
	std::cout << "========================================" << std::endl;

	std::cout << "Total entries: " << entries.size() << std::endl;
	std::cout << "Cache file: " << cache_file << std::endl;

	// Compute statistics
	double avg_construction_time = 0.0;
	double avg_memory_used = 0.0;
	double avg_compression = 0.0;

	for(const auto& entry : entries)
	{
		avg_construction_time += entry.construction_time;
		avg_memory_used += entry.memory_used;
		avg_compression += entry.compression_ratio;
	}

	avg_construction_time /= entries.size();
	avg_memory_used /= entries.size();
	avg_compression /= entries.size();

	std::cout << "\nAverage construction time: " << avg_construction_time << " s" << std::endl;
	std::cout << "Average memory usage: " << (avg_memory_used / 1024 / 1024) << " MB" << std::endl;
	std::cout << "Average compression ratio: " << avg_compression << std::endl;

	std::cout << "========================================" << std::endl;
}

//-------------------------------------------------------------------------
