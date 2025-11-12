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
	: cache_dir(dir), enabled(true), dirty(false),
	  full_serialization_enabled(false), max_cache_size_mb(1000)
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
// Full H-Matrix Serialization (Phase 3B - v1.1.0)
//-------------------------------------------------------------------------

#include "rad_intrc_hmat.h"
#include <fstream>

bool radTHMatrixCache::EnableFullSerialization(bool enable)
{
	full_serialization_enabled = enable;
	return true;
}

void radTHMatrixCache::SetMaxCacheSize(size_t max_mb)
{
	max_cache_size_mb = max_mb;
}

std::string radTHMatrixCache::GetDataDir() const
{
	return cache_dir + "/hmat";
}

std::string radTHMatrixCache::GetDataFilePath(uint64_t hash) const
{
	char filename[256];
	snprintf(filename, sizeof(filename), "%016llx.hmat", (unsigned long long)hash);
	return GetDataDir() + "/" + std::string(filename);
}

size_t radTHMatrixCache::GetCurrentCacheSize() const
{
	// TODO: Implement directory size calculation
	return 0;
}

bool radTHMatrixCache::SaveHMatrix(uint64_t hash, const radTHMatrixInteraction* hmat)
{
	if(!full_serialization_enabled || hmat == nullptr)
		return false;

	// Ensure data directory exists
	std::string data_dir = GetDataDir();
	struct stat info;
	if(stat(data_dir.c_str(), &info) != 0)
	{
		if(mkdir(data_dir.c_str(), 0755) != 0)
		{
			std::cerr << "[Phase 3B] Failed to create data directory: " << data_dir << std::endl;
			return false;
		}
	}

	// Get data file path
	std::string filepath = GetDataFilePath(hash);

	// Open file for writing
	std::ofstream file(filepath, std::ios::binary | std::ios::trunc);
	if(!file.is_open())
	{
		std::cerr << "[Phase 3B] Failed to open file for writing: " << filepath << std::endl;
		return false;
	}

	try
	{
		// Write header
		const uint32_t MAGIC = 0x484D4154;  // "HMAT"
		const uint32_t VERSION = 1;
		const uint32_t HACAPK_VERSION = 130;  // HACApK 1.3.0

		file.write(reinterpret_cast<const char*>(&MAGIC), sizeof(MAGIC));
		file.write(reinterpret_cast<const char*>(&VERSION), sizeof(VERSION));
		file.write(reinterpret_cast<const char*>(&HACAPK_VERSION), sizeof(HACAPK_VERSION));
		file.write(reinterpret_cast<const char*>(&hash), sizeof(hash));

		uint32_t num_elements = hmat->n_elem;
		double eps = hmat->config.eps;
		uint32_t max_rank = hmat->config.max_rank;

		file.write(reinterpret_cast<const char*>(&num_elements), sizeof(num_elements));
		file.write(reinterpret_cast<const char*>(&eps), sizeof(eps));
		file.write(reinterpret_cast<const char*>(&max_rank), sizeof(max_rank));

		// Write all 9 H-matrices (3x3 tensor)
		for(int i = 0; i < 9; i++)
		{
			const hacapk::HMatrix* h = hmat->hmat[i].get();
			if(h == nullptr)
			{
				// Write null marker
				uint32_t is_null = 1;
				file.write(reinterpret_cast<const char*>(&is_null), sizeof(is_null));
				continue;
			}

			// Write non-null marker
			uint32_t is_null = 0;
			file.write(reinterpret_cast<const char*>(&is_null), sizeof(is_null));

			// Write HMatrix metadata
			file.write(reinterpret_cast<const char*>(&h->nd), sizeof(h->nd));
			file.write(reinterpret_cast<const char*>(&h->nlf), sizeof(h->nlf));
			file.write(reinterpret_cast<const char*>(&h->nlfkt), sizeof(h->nlfkt));
			file.write(reinterpret_cast<const char*>(&h->ktmax), sizeof(h->ktmax));

			// Write blocks
			uint32_t num_blocks = static_cast<uint32_t>(h->blocks.size());
			file.write(reinterpret_cast<const char*>(&num_blocks), sizeof(num_blocks));

			for(const auto& block : h->blocks)
			{
				// Write block metadata
				file.write(reinterpret_cast<const char*>(&block.ltmtx), sizeof(block.ltmtx));
				file.write(reinterpret_cast<const char*>(&block.kt), sizeof(block.kt));
				file.write(reinterpret_cast<const char*>(&block.nstrtl), sizeof(block.nstrtl));
				file.write(reinterpret_cast<const char*>(&block.ndl), sizeof(block.ndl));
				file.write(reinterpret_cast<const char*>(&block.nstrtt), sizeof(block.nstrtt));
				file.write(reinterpret_cast<const char*>(&block.ndt), sizeof(block.ndt));

				// Write a1 matrix (U)
				uint32_t a1_size = static_cast<uint32_t>(block.a1.size());
				file.write(reinterpret_cast<const char*>(&a1_size), sizeof(a1_size));
				if(a1_size > 0)
				{
					file.write(reinterpret_cast<const char*>(block.a1.data()), a1_size * sizeof(double));
				}

				// Write a2 matrix (V)
				uint32_t a2_size = static_cast<uint32_t>(block.a2.size());
				file.write(reinterpret_cast<const char*>(&a2_size), sizeof(a2_size));
				if(a2_size > 0)
				{
					file.write(reinterpret_cast<const char*>(block.a2.data()), a2_size * sizeof(double));
				}
			}

			// Write block structure
			uint32_t lbstrtl_size = static_cast<uint32_t>(h->lbstrtl.size());
			file.write(reinterpret_cast<const char*>(&lbstrtl_size), sizeof(lbstrtl_size));
			if(lbstrtl_size > 0)
			{
				file.write(reinterpret_cast<const char*>(h->lbstrtl.data()), lbstrtl_size * sizeof(int));
			}

			uint32_t lbstrtt_size = static_cast<uint32_t>(h->lbstrtt.size());
			file.write(reinterpret_cast<const char*>(&lbstrtt_size), sizeof(lbstrtt_size));
			if(lbstrtt_size > 0)
			{
				file.write(reinterpret_cast<const char*>(h->lbstrtt.data()), lbstrtt_size * sizeof(int));
			}

			uint32_t lbndl_size = static_cast<uint32_t>(h->lbndl.size());
			file.write(reinterpret_cast<const char*>(&lbndl_size), sizeof(lbndl_size));
			if(lbndl_size > 0)
			{
				file.write(reinterpret_cast<const char*>(h->lbndl.data()), lbndl_size * sizeof(int));
			}

			uint32_t lbndt_size = static_cast<uint32_t>(h->lbndt.size());
			file.write(reinterpret_cast<const char*>(&lbndt_size), sizeof(lbndt_size));
			if(lbndt_size > 0)
			{
				file.write(reinterpret_cast<const char*>(h->lbndt.data()), lbndt_size * sizeof(int));
			}
		}

		file.close();

		std::cout << "[Phase 3B] Saved H-matrix to disk: " << filepath << std::endl;
		return true;
	}
	catch(const std::exception& e)
	{
		std::cerr << "[Phase 3B] Error saving H-matrix: " << e.what() << std::endl;
		file.close();
		return false;
	}
}

radTHMatrixInteraction* radTHMatrixCache::LoadHMatrix(uint64_t hash, radTInteraction* intrct_ptr)
{
	if(!full_serialization_enabled || intrct_ptr == nullptr)
		return nullptr;

	std::string filepath = GetDataFilePath(hash);

	// Check if file exists
	std::ifstream file(filepath, std::ios::binary);
	if(!file.is_open())
	{
		// File doesn't exist - cache miss
		return nullptr;
	}

	try
	{
		// Read and validate header
		uint32_t magic, version, hacapk_version;
		uint64_t file_hash;

		file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
		file.read(reinterpret_cast<char*>(&version), sizeof(version));
		file.read(reinterpret_cast<char*>(&hacapk_version), sizeof(hacapk_version));
		file.read(reinterpret_cast<char*>(&file_hash), sizeof(file_hash));

		const uint32_t EXPECTED_MAGIC = 0x484D4154;  // "HMAT"
		const uint32_t EXPECTED_VERSION = 1;
		const uint32_t EXPECTED_HACAPK = 130;

		if(magic != EXPECTED_MAGIC)
		{
			std::cerr << "[Phase 3B] Invalid magic number in cache file" << std::endl;
			file.close();
			return nullptr;
		}

		if(version != EXPECTED_VERSION)
		{
			std::cerr << "[Phase 3B] Version mismatch (file=" << version << ", expected=" << EXPECTED_VERSION << ")" << std::endl;
			file.close();
			return nullptr;
		}

		if(hacapk_version != EXPECTED_HACAPK)
		{
			std::cerr << "[Phase 3B] HACApK version mismatch - invalidating cache" << std::endl;
			file.close();
			return nullptr;
		}

		if(file_hash != hash)
		{
			std::cerr << "[Phase 3B] Hash mismatch in cache file" << std::endl;
			file.close();
			return nullptr;
		}

		// Read metadata
		uint32_t num_elements;
		double eps;
		uint32_t max_rank;

		file.read(reinterpret_cast<char*>(&num_elements), sizeof(num_elements));
		file.read(reinterpret_cast<char*>(&eps), sizeof(eps));
		file.read(reinterpret_cast<char*>(&max_rank), sizeof(max_rank));

		// Create config
		radTHMatrixSolverConfig config;
		config.eps = eps;
		config.max_rank = max_rank;
		config.min_cluster_size = 10;
		config.use_openmp = true;
		config.num_threads = 0;

		// Create H-matrix interaction object
		radTHMatrixInteraction* hmat = new radTHMatrixInteraction(intrct_ptr, config);

		// Read all 9 H-matrices
		for(int i = 0; i < 9; i++)
		{
			uint32_t is_null;
			file.read(reinterpret_cast<char*>(&is_null), sizeof(is_null));

			if(is_null == 1)
			{
				hmat->hmat[i] = nullptr;
				continue;
			}

			// Create new HMatrix
			hmat->hmat[i] = std::make_unique<hacapk::HMatrix>();
			hacapk::HMatrix* h = hmat->hmat[i].get();

			// Read HMatrix metadata
			file.read(reinterpret_cast<char*>(&h->nd), sizeof(h->nd));
			file.read(reinterpret_cast<char*>(&h->nlf), sizeof(h->nlf));
			file.read(reinterpret_cast<char*>(&h->nlfkt), sizeof(h->nlfkt));
			file.read(reinterpret_cast<char*>(&h->ktmax), sizeof(h->ktmax));

			// Read blocks
			uint32_t num_blocks;
			file.read(reinterpret_cast<char*>(&num_blocks), sizeof(num_blocks));

			h->blocks.resize(num_blocks);
			for(uint32_t j = 0; j < num_blocks; j++)
			{
				hacapk::LowRankBlock& block = h->blocks[j];

				// Read block metadata
				file.read(reinterpret_cast<char*>(&block.ltmtx), sizeof(block.ltmtx));
				file.read(reinterpret_cast<char*>(&block.kt), sizeof(block.kt));
				file.read(reinterpret_cast<char*>(&block.nstrtl), sizeof(block.nstrtl));
				file.read(reinterpret_cast<char*>(&block.ndl), sizeof(block.ndl));
				file.read(reinterpret_cast<char*>(&block.nstrtt), sizeof(block.nstrtt));
				file.read(reinterpret_cast<char*>(&block.ndt), sizeof(block.ndt));

				// Read a1 matrix (U)
				uint32_t a1_size;
				file.read(reinterpret_cast<char*>(&a1_size), sizeof(a1_size));
				if(a1_size > 0)
				{
					block.a1.resize(a1_size);
					file.read(reinterpret_cast<char*>(block.a1.data()), a1_size * sizeof(double));
				}

				// Read a2 matrix (V)
				uint32_t a2_size;
				file.read(reinterpret_cast<char*>(&a2_size), sizeof(a2_size));
				if(a2_size > 0)
				{
					block.a2.resize(a2_size);
					file.read(reinterpret_cast<char*>(block.a2.data()), a2_size * sizeof(double));
				}
			}

			// Read block structure
			uint32_t lbstrtl_size;
			file.read(reinterpret_cast<char*>(&lbstrtl_size), sizeof(lbstrtl_size));
			if(lbstrtl_size > 0)
			{
				h->lbstrtl.resize(lbstrtl_size);
				file.read(reinterpret_cast<char*>(h->lbstrtl.data()), lbstrtl_size * sizeof(int));
			}

			uint32_t lbstrtt_size;
			file.read(reinterpret_cast<char*>(&lbstrtt_size), sizeof(lbstrtt_size));
			if(lbstrtt_size > 0)
			{
				h->lbstrtt.resize(lbstrtt_size);
				file.read(reinterpret_cast<char*>(h->lbstrtt.data()), lbstrtt_size * sizeof(int));
			}

			uint32_t lbndl_size;
			file.read(reinterpret_cast<char*>(&lbndl_size), sizeof(lbndl_size));
			if(lbndl_size > 0)
			{
				h->lbndl.resize(lbndl_size);
				file.read(reinterpret_cast<char*>(h->lbndl.data()), lbndl_size * sizeof(int));
			}

			uint32_t lbndt_size;
			file.read(reinterpret_cast<char*>(&lbndt_size), sizeof(lbndt_size));
			if(lbndt_size > 0)
			{
				h->lbndt.resize(lbndt_size);
				file.read(reinterpret_cast<char*>(h->lbndt.data()), lbndt_size * sizeof(int));
			}
		}

		// Mark as built
		hmat->is_built = true;
		hmat->n_elem = num_elements;

		// Copy statistics from cache entry
		const radTHMatrixCacheEntry* entry = Find(hash);
		if(entry)
		{
			hmat->construction_time = entry->construction_time;
			hmat->memory_used = entry->memory_used;
			hmat->compression_ratio = entry->compression_ratio;
		}

		file.close();

		std::cout << "[Phase 3B] Loaded H-matrix from disk: " << filepath << std::endl;
		return hmat;
	}
	catch(const std::exception& e)
	{
		std::cerr << "[Phase 3B] Error loading H-matrix: " << e.what() << std::endl;
		file.close();
		return nullptr;
	}
}

//-------------------------------------------------------------------------
