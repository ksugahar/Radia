#ifndef RAD_POINT_CLASSIFY_H
#define RAD_POINT_CLASSIFY_H

/**
 * @file rad_point_classify.h
 * @brief Point classification for FMM field computation
 *
 * Classifies evaluation points as inside/near/far relative to magnetic elements.
 * Used by FMM-accelerated field computation for hybrid near/far correction.
 *
 * Classification:
 *   0 = INSIDE: Point is inside an element (use element's internal field)
 *   1 = NEAR: Point is within near_threshold * elem_size of an element
 *             (requires exact MSC field computation with near correction)
 *   2 = FAR: Point is far from all elements (FMM dipole approximation is accurate)
 */

#include "rad_polyhedron.h"
#include "rad_geom_types.h"
#include <vector>

namespace RadPointClassify {

/// Point classification enum
enum PointClass {
    INSIDE = 0,  ///< Point is inside an element
    NEAR = 1,    ///< Point is near an element (needs correction)
    FAR = 2      ///< Point is far from all elements (FMM accurate)
};

/// Classification result for a single point
struct ClassifyResult {
    PointClass classification;  ///< INSIDE, NEAR, or FAR
    int nearest_elem_id;        ///< Index of nearest element (0-based)
    double distance;            ///< Distance to nearest element center
};

/**
 * @brief Check if a point is inside a tetrahedron using signed volumes
 *
 * Uses the signed volume method: a point is inside if all four signed volumes
 * (formed with the point and three vertices) have the same sign.
 *
 * @param pt Point to test [x, y, z]
 * @param v0, v1, v2, v3 Tetrahedron vertices
 * @return true if point is inside tetrahedron
 */
bool PointInTetrahedron(const TVector3d& pt,
                        const TVector3d& v0, const TVector3d& v1,
                        const TVector3d& v2, const TVector3d& v3);

/**
 * @brief Check if a point is inside a hexahedron using 5-tetrahedra decomposition
 *
 * Decomposes the hexahedron into 5 tetrahedra and checks each.
 * More robust than Newton's method for general hexahedra.
 *
 * @param pt Point to test [x, y, z]
 * @param verts 8 hexahedron vertices (Radia ordering)
 * @return true if point is inside hexahedron
 */
bool PointInHexahedron(const TVector3d& pt, const TVector3d verts[8]);

/**
 * @brief Compute element characteristic size from AABB
 *
 * Returns the cube root of the AABB volume as the characteristic size.
 *
 * @param vertices Element vertices
 * @param n_verts Number of vertices (4 for tetra, 8 for hex)
 * @return Characteristic element size
 */
double ComputeElementSize(const TVector3d* vertices, int n_verts);

/**
 * @brief Compute AABB (Axis-Aligned Bounding Box) for element
 *
 * @param vertices Element vertices
 * @param n_verts Number of vertices
 * @param aabb_min Output: minimum corner of AABB
 * @param aabb_max Output: maximum corner of AABB
 */
void ComputeElementAABB(const TVector3d* vertices, int n_verts,
                        TVector3d& aabb_min, TVector3d& aabb_max);

/**
 * @brief Classify multiple points relative to mesh elements
 *
 * For each point, determines:
 * - INSIDE: Point is geometrically inside an element
 * - NEAR: Point is within near_threshold * elem_size of nearest element
 * - FAR: Point is outside the near zone of all elements
 *
 * Algorithm:
 * 1. Compute global AABB and average element size
 * 2. For each point:
 *    a. Quick rejection if outside global AABB + margin
 *    b. Find nearest element by center distance
 *    c. Check if inside any nearby element
 *    d. Classify as near/far based on distance threshold
 *
 * @param n_points Number of evaluation points
 * @param points Array of evaluation points [n_points][3]
 * @param elements Vector of element data (vertices, centers)
 * @param near_threshold Near zone multiplier (typically 3.0)
 * @param results Output: classification results for each point
 */
void ClassifyPoints(int n_points,
                    const double* points,
                    const std::vector<radTPolygon>& elements,
                    double near_threshold,
                    std::vector<ClassifyResult>& results);

/**
 * @brief Classify points using Radia element handles
 *
 * Wrapper that extracts element data from Radia handles and calls ClassifyPoints.
 *
 * @param n_points Number of evaluation points
 * @param points Array of evaluation points [n_points * 3]
 * @param container_handle Radia container handle
 * @param near_threshold Near zone multiplier (typically 3.0)
 * @param classification Output: classification (0=inside, 1=near, 2=far) [n_points]
 * @param nearest_elem Output: nearest element index [n_points]
 */
void ClassifyPointsFromHandle(int n_points,
                              const double* points,
                              int container_handle,
                              double near_threshold,
                              int* classification,
                              int* nearest_elem);

} // namespace RadPointClassify

#endif // RAD_POINT_CLASSIFY_H
