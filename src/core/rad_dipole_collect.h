/*-------------------------------------------------------------------------
 *
 * File name:      rad_dipole_collect.h
 *
 * Project:        RADIA
 *
 * Description:    Collect dipole data from Radia containers for FMM computation
 *
 * Author(s):      Claude Code (AI-assisted development)
 *
 * First release:  2025-12-30
 *
 * Copyright (C):  European Synchrotron Radiation Facility, France
 *
 *-------------------------------------------------------------------------
 *
 * This module provides functions to extract dipole (position, magnetization,
 * volume) data from Radia magnetic field source objects for use in FMM
 * (Fast Multipole Method) accelerated field computation.
 *
 * Each magnetic element is approximated as a point dipole:
 *   - Position: Element center point
 *   - Dipole moment: m = M * V (magnetization * volume)
 *
 *-------------------------------------------------------------------------*/

#ifndef RAD_DIPOLE_COLLECT_H
#define RAD_DIPOLE_COLLECT_H

#include <vector>
#include <cstdint>

// Forward declarations - avoid heavy includes
class radTg;
class radTg3d;
class radTGroup;
class radTg3dRelax;
class radTrans;

namespace RadDipoleCollect {

/**
 * Structure to hold dipole data for a single element
 */
struct DipoleData {
    double x, y, z;        // Position (center point) in current units
    double mx, my, mz;     // Dipole moment = M * V (magnetization * volume)
    double volume;         // Volume for reference
};

/**
 * Collection of dipoles extracted from a Radia object
 */
struct DipoleCollection {
    std::vector<DipoleData> dipoles;

    // Flattened arrays for FMM interface
    std::vector<double> positions;    // [x1,y1,z1, x2,y2,z2, ...]
    std::vector<double> moments;      // [mx1,my1,mz1, mx2,my2,mz2, ...]

    int64_t count() const { return static_cast<int64_t>(dipoles.size()); }

    // Flatten dipole data into position/moment arrays for FMM
    void flatten();

    // Clear all data
    void clear();
};

/**
 * Collect all dipoles from a Radia object (container or single element)
 *
 * Recursively traverses containers (radTGroup) and extracts dipole data
 * from relaxable elements (radTg3dRelax and derived classes).
 *
 * For elements with transformations (symmetry), all symmetry copies are included.
 *
 * @param g3dPtr      Pointer to Radia 3D object (may be container or single element)
 * @param collection  Output collection to append dipoles to
 * @param pBaseTrans  Optional base transformation (for nested symmetry handling)
 */
void CollectDipoles(radTg3d* g3dPtr, DipoleCollection& collection, radTrans* pBaseTrans = nullptr);

/**
 * Collect dipoles from a Group container
 *
 * Internal helper function - recursively processes group members.
 *
 * @param groupPtr    Pointer to radTGroup
 * @param collection  Output collection
 * @param pBaseTrans  Optional base transformation
 */
void CollectDipolesFromGroup(radTGroup* groupPtr, DipoleCollection& collection, radTrans* pBaseTrans);

/**
 * Extract dipole data from a single relaxable element
 *
 * Computes center point and dipole moment (m = M * V) for the element.
 * Handles transformation/symmetry if present.
 *
 * @param relaxPtr    Pointer to radTg3dRelax element
 * @param collection  Output collection
 * @param pBaseTrans  Optional transformation to apply
 */
void ExtractDipoleFromRelax(radTg3dRelax* relaxPtr, DipoleCollection& collection, radTrans* pBaseTrans);

/**
 * Get total dipole count without collecting full data
 *
 * Useful for pre-allocating arrays.
 *
 * @param g3dPtr      Pointer to Radia 3D object
 * @return            Total number of dipoles (including symmetry copies)
 */
int64_t CountDipoles(radTg3d* g3dPtr);

/**
 * Collect dipoles from a Radia object by its element key
 *
 * This is the main entry point for external use (e.g., from radia_ngsolve).
 * Uses the global radTApplication to look up the object.
 *
 * @param elemKey     Radia element key (object ID)
 * @param collection  Output collection to fill with dipole data
 * @return            true if successful, false if object not found
 */
bool CollectDipolesFromKey(int elemKey, DipoleCollection& collection);

} // namespace RadDipoleCollect

#endif // RAD_DIPOLE_COLLECT_H
