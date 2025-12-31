/*-------------------------------------------------------------------------
 *
 * File name:      rad_dipole_collect.cpp
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
 * FMM Dipole Collection for Tetrahedron and Hexahedron Elements
 *
 * This module extracts dipole (position, magnetization * volume) data from
 * Radia magnetic field source objects (radTg3dRelax and derived classes).
 * Currently supports:
 *   - Tetrahedral elements (4 faces, 3 DOF MMM)
 *   - Hexahedral elements (6 faces, 6 DOF MSC - converted to effective dipole)
 *
 * Each element is approximated as a point dipole at its center with
 * dipole moment m = M * V (magnetization times volume).
 *
 *-------------------------------------------------------------------------*/

#include "rad_dipole_collect.h"
#include "rad_geometry_3d.h"
#include "rad_group.h"
#include "rad_transform_def.h"
#include "rad_application.h"

// Access global Radia application instance
extern radTApplication rad;

namespace RadDipoleCollect {

//-------------------------------------------------------------------------
// Forward declarations for internal helper functions
//-------------------------------------------------------------------------

static void CollectDipolesInternal(radTg3d* g3dPtr, DipoleCollection& collection,
                                   radTrans* pTrans);

//-------------------------------------------------------------------------

void DipoleCollection::flatten()
{
    const size_t n = dipoles.size();
    positions.resize(n * 3);
    moments.resize(n * 3);

    for (size_t i = 0; i < n; ++i) {
        const DipoleData& d = dipoles[i];
        positions[i * 3 + 0] = d.x;
        positions[i * 3 + 1] = d.y;
        positions[i * 3 + 2] = d.z;
        moments[i * 3 + 0] = d.mx;
        moments[i * 3 + 1] = d.my;
        moments[i * 3 + 2] = d.mz;
    }
}

//-------------------------------------------------------------------------

void DipoleCollection::clear()
{
    dipoles.clear();
    positions.clear();
    moments.clear();
}

//-------------------------------------------------------------------------

void CollectDipoles(radTg3d* g3dPtr, DipoleCollection& collection, radTrans* pBaseTrans)
{
    if (g3dPtr == nullptr) return;

    // For elements with transformations, we need to process each copy
    if (!g3dPtr->g3dListOfTransform.empty()) {
        // Process original element first (no transformation applied)
        CollectDipolesInternal(g3dPtr, collection, pBaseTrans);

        // Process each transformation (symmetry copy)
        for (auto& trPair : g3dPtr->g3dListOfTransform) {
            int multiplicity = trPair.m;
            radTrans* pSymTrans = static_cast<radTrans*>(trPair.Handler_g.rep);
            if (pSymTrans == nullptr || multiplicity <= 0) continue;

            // Apply transformation iteratively for multiplicity
            radTrans accumTrans;
            accumTrans.SetupIdent();

            for (int i = 0; i < multiplicity; ++i) {
                // Compose: accumTrans = pSymTrans * accumTrans
                radTrans newAccum;
                TrProduct(pSymTrans, &accumTrans, newAccum);
                accumTrans = newAccum;

                // Compose with base transformation if present
                radTrans* pFinalTrans;
                radTrans composedTrans;
                if (pBaseTrans != nullptr) {
                    TrProduct(pBaseTrans, &accumTrans, composedTrans);
                    pFinalTrans = &composedTrans;
                } else {
                    pFinalTrans = &accumTrans;
                }

                // Collect from this transformed copy
                CollectDipolesInternal(g3dPtr, collection, pFinalTrans);
            }
        }
    } else {
        // No transformations - collect directly
        CollectDipolesInternal(g3dPtr, collection, pBaseTrans);
    }
}

//-------------------------------------------------------------------------

/**
 * Internal helper: Collect dipoles without processing transformation list
 *
 * This processes the element itself without iterating over its symmetry list.
 */
static void CollectDipolesInternal(radTg3d* g3dPtr, DipoleCollection& collection,
                                   radTrans* pTrans)
{
    if (g3dPtr == nullptr) return;

    // Check object type
    int type = g3dPtr->Type_g3d();

    // Type 2 = radTGroup (container)
    if (type == 2) {
        radTGroup* groupPtr = static_cast<radTGroup*>(g3dPtr);
        CollectDipolesFromGroup(groupPtr, collection, pTrans);
        return;
    }

    // Type 1 = radTg3dRelax (relaxable element with magnetization)
    if (type == 1) {
        radTg3dRelax* relaxPtr = static_cast<radTg3dRelax*>(g3dPtr);
        ExtractDipoleFromRelax(relaxPtr, collection, pTrans);
        return;
    }

    // Other types (coils, etc.) - currently not supported for FMM dipole extraction
    // Coils have their own field computation methods (elliptic integrals)
}

//-------------------------------------------------------------------------

void CollectDipolesFromGroup(radTGroup* groupPtr, DipoleCollection& collection, radTrans* pBaseTrans)
{
    if (groupPtr == nullptr) return;

    // Iterate over all members in the group
    for (auto& pair : groupPtr->GroupMapOfHandlers) {
        radTg3d* memberPtr = static_cast<radTg3d*>(pair.second.rep);
        if (memberPtr == nullptr) continue;

        // Recursively collect (handles nested groups and transformations)
        CollectDipoles(memberPtr, collection, pBaseTrans);
    }
}

//-------------------------------------------------------------------------

void ExtractDipoleFromRelax(radTg3dRelax* relaxPtr, DipoleCollection& collection, radTrans* pTrans)
{
    if (relaxPtr == nullptr) return;

    // Get center point and magnetization
    TVector3d center = relaxPtr->CentrPoint;
    TVector3d magn = relaxPtr->Magn;
    double volume = relaxPtr->Volume();

    // Skip if volume is zero or negligible
    if (volume <= 0.0) return;

    // Apply transformation if present
    if (pTrans != nullptr) {
        // Transform position: P' = M*P + V
        center = pTrans->TrPoint(center);

        // Transform magnetization vector (pseudovector/axial vector)
        // M' = s * (M * m) where s accounts for reflections
        magn = pTrans->TrVectField(magn);
    }

    // Compute dipole moment: m = M * V
    DipoleData dipole;
    dipole.x = center.x;
    dipole.y = center.y;
    dipole.z = center.z;
    dipole.mx = magn.x * volume;
    dipole.my = magn.y * volume;
    dipole.mz = magn.z * volume;
    dipole.volume = volume;

    // Skip zero-moment dipoles (but allow small non-zero values)
    double momentMagSq = dipole.mx * dipole.mx + dipole.my * dipole.my + dipole.mz * dipole.mz;
    if (momentMagSq > 0.0) {
        collection.dipoles.push_back(dipole);
    }
}

//-------------------------------------------------------------------------

int64_t CountDipoles(radTg3d* g3dPtr)
{
    if (g3dPtr == nullptr) return 0;

    int type = g3dPtr->Type_g3d();

    // Type 2 = radTGroup
    if (type == 2) {
        radTGroup* groupPtr = static_cast<radTGroup*>(g3dPtr);
        int64_t count = 0;
        for (auto& pair : groupPtr->GroupMapOfHandlers) {
            radTg3d* memberPtr = static_cast<radTg3d*>(pair.second.rep);
            count += CountDipoles(memberPtr);
        }
        return count;
    }

    // Type 1 = radTg3dRelax (tetrahedron, hexahedron, etc.)
    if (type == 1) {
        // Count symmetry multiplicity
        int totalMult = g3dPtr->TotalMultiplicity();
        return totalMult > 0 ? totalMult : 1;
    }

    return 0;
}

//-------------------------------------------------------------------------

bool CollectDipolesFromKey(int elemKey, DipoleCollection& collection)
{
    // Look up object in global Radia application
    radThg hg;
    if (!rad.UnsafeGetElemByKey(elemKey, hg)) {
        return false;  // Object not found
    }

    // Check if it's a g3d object (Type_g() == 1)
    radTg* gPtr = hg.rep;
    if (gPtr == nullptr || gPtr->Type_g() != 1) {
        return false;  // Not a 3D magnetic object
    }

    // Cast to radTg3d and collect dipoles
    radTg3d* g3dPtr = static_cast<radTg3d*>(gPtr);
    collection.clear();
    CollectDipoles(g3dPtr, collection, nullptr);
    collection.flatten();

    return true;
}

//-------------------------------------------------------------------------

} // namespace RadDipoleCollect
