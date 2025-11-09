/*-------------------------------------------------------------------------
*
* File name:      radentry_hmat.cpp
*
* Project:        RADIA
*
* Description:    H-matrix field evaluation API implementation
*
* Author(s):      Radia Development Team
*
* First release:  2025
*
-------------------------------------------------------------------------*/

#include "radentry.h"
#include "radentry_hmat.h"
#include "radappl.h"
#include "radhmat.h"
#include "radgroup.h"
#include "radg.h"
#include "radcast.h"
#include "gmvect.h"
#include <cstring>
#include <map>
#include <memory>

//-------------------------------------------------------------------------

extern radTApplication rad;

//-------------------------------------------------------------------------
// Global state for H-matrix field evaluation
//-------------------------------------------------------------------------

struct HMatrixFieldGlobalState {
	bool enabled;
	double epsilon;
	std::map<int, std::shared_ptr<radTHMatrixFieldEvaluator>> cache;

	HMatrixFieldGlobalState() : enabled(false), epsilon(1e-6) {}

	void Clear() {
		cache.clear();
	}

	radTHMatrixFieldEvaluator* GetOrCreate(int obj_key, radTGroup* group) {
		auto it = cache.find(obj_key);
		if(it != cache.end()) {
			// Check if cache is still valid
			if(it->second && it->second->IsValid(group)) {
				return it->second.get();
			}
			// Invalid, remove from cache
			cache.erase(it);
		}

		// Create new evaluator
		if(!group) return nullptr;

		radTHMatrixConfig config;
		config.eps = epsilon;

		auto evaluator = std::make_shared<radTHMatrixFieldEvaluator>(config);

		// Build H-matrix (or prepare for direct calculation)
		evaluator->Build(group);

		cache[obj_key] = evaluator;
		return evaluator.get();
	}
};

static HMatrixFieldGlobalState g_hmatrix_field_state;

//-------------------------------------------------------------------------
// API Implementation
//-------------------------------------------------------------------------

EXP int CALL RadFldBatch(double* B, int* nB, int obj, char* id, double* Coords, int np, int use_hmatrix)
{
	if(!B || !nB || !id || !Coords || np <= 0) {
		return -1;  // Error
	}

	try {
		// Get object from application
		radThg hg;
		if(!rad.ValidateElemKey(obj, hg)) {
			return -1;  // Error
		}

		// Try to cast to radTGroup
		radTCast Cast;
		radTg3d* g3d = Cast.g3dCast(hg.rep);
		if(!g3d) return -1;

		radTGroup* group = Cast.GroupCast(g3d);

		// Parse field type
		char field_comp = 'h';  // Default to H-field
		if(id && strlen(id) > 0) {
			if(id[0] == 'b' || id[0] == 'B') field_comp = 'b';
			else if(id[0] == 'h' || id[0] == 'H') field_comp = 'h';
			else if(id[0] == 'a' || id[0] == 'A') field_comp = 'a';
			else if(id[0] == 'm' || id[0] == 'M') field_comp = 'm';
		}

		// Convert coordinates to vector
		std::vector<TVector3d> obs_points;
		obs_points.reserve(np);
		for(int i = 0; i < np; i++) {
			obs_points.push_back(TVector3d(
				Coords[i*3 + 0],
				Coords[i*3 + 1],
				Coords[i*3 + 2]
			));
		}

		std::vector<TVector3d> field_out;

		// Decide whether to use H-matrix
		bool should_use_hmatrix = (use_hmatrix != 0) && g_hmatrix_field_state.enabled;

		if(should_use_hmatrix && group) {
			// Try H-matrix accelerated evaluation
			radTHMatrixFieldEvaluator* evaluator = g_hmatrix_field_state.GetOrCreate(obj, group);

			if(evaluator) {
				int result = evaluator->EvaluateField(obs_points, field_out, field_comp);
				if(result == 0) {
					// H-matrix failed, fallback to direct
					should_use_hmatrix = false;
				}
			} else {
				should_use_hmatrix = false;
			}
		} else {
			should_use_hmatrix = false;
		}

		// Fallback to direct calculation if H-matrix not used
		if(!should_use_hmatrix) {
			field_out.resize(np, TVector3d(0, 0, 0));

			// Direct calculation using existing Radia functions
			TVector3d ZeroVect(0, 0, 0);
			for(int i = 0; i < np; i++) {
				radTFieldKey fieldKey;
				// Set appropriate field key based on field_comp
				if(field_comp == 'b' || field_comp == 'B') {
					fieldKey.B_ = 1;
				} else if(field_comp == 'h' || field_comp == 'H') {
					fieldKey.H_ = 1;
				} else if(field_comp == 'a' || field_comp == 'A') {
					fieldKey.A_ = 1;
				} else if(field_comp == 'm' || field_comp == 'M') {
					fieldKey.M_ = 1;
				}

				radTField field(fieldKey, obs_points[i], ZeroVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);

				// Compute field
				g3d->B_comp(&field);

				// Extract requested component
				if(field_comp == 'b' || field_comp == 'B') {
					field_out[i] = field.B;
				} else if(field_comp == 'h' || field_comp == 'H') {
					field_out[i] = field.H;
				} else if(field_comp == 'a' || field_comp == 'A') {
					field_out[i] = field.A;
				} else if(field_comp == 'm' || field_comp == 'M') {
					field_out[i] = field.M;
				}
			}
		}

		// Copy results to output array
		*nB = np * 3;
		for(int i = 0; i < np; i++) {
			B[i*3 + 0] = field_out[i].x;
			B[i*3 + 1] = field_out[i].y;
			B[i*3 + 2] = field_out[i].z;
		}

		return 0;  // Success

	} catch(...) {
		return -1;  // Error
	}
}

//-------------------------------------------------------------------------

EXP int CALL RadSetHMatrixFieldEval(int enabled, double tol)
{
	g_hmatrix_field_state.enabled = (enabled != 0);
	g_hmatrix_field_state.epsilon = (tol > 0) ? tol : 1e-6;

	if(!g_hmatrix_field_state.enabled) {
		g_hmatrix_field_state.Clear();
	}

	return 0;  // Success
}

//-------------------------------------------------------------------------

EXP int CALL RadClearHMatrixCache(void)
{
	g_hmatrix_field_state.Clear();
	return 0;  // Success
}

//-------------------------------------------------------------------------

EXP int CALL RadGetHMatrixStats(double* stats, int* nstats)
{
	if(!stats || !nstats) return -1;

	// stats[0] = is_enabled
	// stats[1] = num_cached
	// stats[2] = total_memory_MB

	stats[0] = g_hmatrix_field_state.enabled ? 1.0 : 0.0;
	stats[1] = static_cast<double>(g_hmatrix_field_state.cache.size());

	double total_memory = 0.0;
	for(const auto& pair : g_hmatrix_field_state.cache) {
		if(pair.second) {
			total_memory += pair.second->GetMemoryUsage();
		}
	}
	stats[2] = total_memory / (1024.0 * 1024.0);  // Convert to MB

	*nstats = 3;

	return 0;  // Success
}

//-------------------------------------------------------------------------

EXP int CALL RadUpdateHMatrixMagnetization(int obj)
{
	try {
		// Get object from application
		radThg hg;
		if(!rad.ValidateElemKey(obj, hg)) {
			return -1;  // Error: invalid object
		}

		// Cast to radTGroup
		radTCast Cast;
		radTg3d* g3d = Cast.g3dCast(hg.rep);
		if(!g3d) return -1;

		radTGroup* group = Cast.GroupCast(g3d);
		if(!group) {
			// H-matrix only works with groups
			return -1;
		}

		// Find cached evaluator
		auto it = g_hmatrix_field_state.cache.find(obj);
		if(it == g_hmatrix_field_state.cache.end()) {
			// No cached H-matrix, nothing to update
			return -2;  // Error: H-matrix not built yet
		}

		radTHMatrixFieldEvaluator* evaluator = it->second.get();
		if(!evaluator) return -1;

		// Update magnetization (fast, no H-matrix rebuild)
		int result = evaluator->UpdateMagnetization(group);

		return result ? 0 : -1;
	}
	catch(...) {
		return -1;
	}
}

//-------------------------------------------------------------------------
