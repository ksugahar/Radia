/*-------------------------------------------------------------------------
*
* File name:      radhmat_update.cpp
*
* Project:        RADIA
*
* Description:    Fast magnetization update for H-matrix (non-linear relaxation)
*
* Author(s):      Radia Development Team
*
* First release:  2025
*
*------------------------------------------------------------------------*/

#include "rad_hmatrix.h"
#include "rad_group.h"
#include "rad_geometry_3d.h"
#include <iostream>
#include <functional>

//-------------------------------------------------------------------------

int radTHMatrixFieldEvaluator::UpdateMagnetization(radTGroup* source_group)
{
	if(!is_built) {
		std::cerr << "[HMatrix Field] H-matrix not built, cannot update magnetization" << std::endl;
		return 0;
	}

	if(!source_group) {
		std::cerr << "[HMatrix Field] NULL source group" << std::endl;
		return 0;
	}

	std::cout << "[HMatrix Field] Updating magnetization (fast, no H-matrix rebuild)" << std::endl;

	// Clear and rebuild source moments
	source_moments.clear();
	source_moments.reserve(num_sources * 3);

	// Extract updated magnetic moments (same order as original Build())
	int extracted = 0;

	// Helper function to extract moments (same recursion as ExtractLeafElements)
	std::function<void(radTg3d*, int)> ExtractMoments = [&](radTg3d* g3d, int depth) {
		if(!g3d) return;

		radTGroup* group = dynamic_cast<radTGroup*>(g3d);

		if(group && group->GroupMapOfHandlers.size() > 0) {
			// Container - recursively extract
			for(auto& elem_pair : group->GroupMapOfHandlers) {
				radTg3d* sub_elem = (radTg3d*)(elem_pair.second.rep);
				ExtractMoments(sub_elem, depth + 1);
			}
		}
		else {
			// Leaf element
			radTg3dRelax* relaxable = dynamic_cast<radTg3dRelax*>(g3d);
			if(!relaxable) return;

			// Get UPDATED magnetization
			TVector3d M = relaxable->Magn;

			// Volume (unchanged)
			double volume = relaxable->Volume();
			double vol_m3 = volume * 1e-9;  // mm^3 to m^3

			// Magnetic moment = M * V
			TVector3d moment = M * vol_m3;

			source_moments.push_back(moment.x);
			source_moments.push_back(moment.y);
			source_moments.push_back(moment.z);

			extracted++;
		}
	};

	// Extract from all top-level elements
	for(auto& elem_pair : source_group->GroupMapOfHandlers) {
		radTg3d* g3d = (radTg3d*)(elem_pair.second.rep);
		ExtractMoments(g3d, 0);
	}

	if(extracted != num_sources) {
		std::cerr << "[HMatrix Field] ERROR: Source count mismatch!" << std::endl;
		std::cerr << "  Expected: " << num_sources << ", Got: " << extracted << std::endl;
		return 0;
	}

	std::cout << "[HMatrix Field] Updated " << num_sources << " source moments" << std::endl;
	std::cout << "[HMatrix Field] H-matrix NOT rebuilt (geometry unchanged)" << std::endl;

	return 1;
}

//-------------------------------------------------------------------------
