/**
 * @file rad_point_classify.cpp
 * @brief Point classification implementation for FMM field computation
 */

#include "rad_point_classify.h"
#include "rad_polyhedron.h"
#include "rad_group.h"
#include "rad_application.h"
#include <cmath>
#include <algorithm>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

// Access to global Radia application
extern radTApplication rad;

namespace RadPointClassify {

//-----------------------------------------------------------------------------
// Helper: Compute signed volume of tetrahedron formed by 4 points
//-----------------------------------------------------------------------------
static double SignedVolume(const TVector3d& a, const TVector3d& b,
                           const TVector3d& c, const TVector3d& d)
{
    // Volume = (1/6) * (b-a) . ((c-a) x (d-a))
    TVector3d ab = b - a;
    TVector3d ac = c - a;
    TVector3d ad = d - a;

    // Cross product ac x ad
    TVector3d cross;
    cross.x = ac.y * ad.z - ac.z * ad.y;
    cross.y = ac.z * ad.x - ac.x * ad.z;
    cross.z = ac.x * ad.y - ac.y * ad.x;

    // Dot product ab . cross
    return (ab.x * cross.x + ab.y * cross.y + ab.z * cross.z) / 6.0;
}

//-----------------------------------------------------------------------------
// Check if point is inside tetrahedron using signed volumes
//-----------------------------------------------------------------------------
bool PointInTetrahedron(const TVector3d& pt,
                        const TVector3d& v0, const TVector3d& v1,
                        const TVector3d& v2, const TVector3d& v3)
{
    // Compute signed volumes of 4 tetrahedra formed by replacing each vertex with pt
    double vol0 = SignedVolume(pt, v1, v2, v3);
    double vol1 = SignedVolume(v0, pt, v2, v3);
    double vol2 = SignedVolume(v0, v1, pt, v3);
    double vol3 = SignedVolume(v0, v1, v2, pt);

    // Point is inside if all volumes have the same sign (or are zero)
    // Use small tolerance for numerical stability
    const double eps = 1e-12;

    bool all_positive = (vol0 >= -eps && vol1 >= -eps && vol2 >= -eps && vol3 >= -eps);
    bool all_negative = (vol0 <= eps && vol1 <= eps && vol2 <= eps && vol3 <= eps);

    return all_positive || all_negative;
}

//-----------------------------------------------------------------------------
// Check if point is inside hexahedron using 5-tetrahedra decomposition
//-----------------------------------------------------------------------------
bool PointInHexahedron(const TVector3d& pt, const TVector3d verts[8])
{
    // Radia hexahedron vertex ordering (0-indexed):
    //   4----7
    //  /|   /|
    // 5----6 |
    // | 0--|-3
    // |/   |/
    // 1----2
    //
    // Decompose into 5 tetrahedra:
    // 1. (0, 1, 3, 4)
    // 2. (1, 2, 3, 6)
    // 3. (1, 4, 5, 6)
    // 4. (3, 4, 6, 7)
    // 5. (1, 3, 4, 6) - central tetrahedron

    // Tet 1: vertices 0, 1, 3, 4
    if (PointInTetrahedron(pt, verts[0], verts[1], verts[3], verts[4]))
        return true;

    // Tet 2: vertices 1, 2, 3, 6
    if (PointInTetrahedron(pt, verts[1], verts[2], verts[3], verts[6]))
        return true;

    // Tet 3: vertices 1, 4, 5, 6
    if (PointInTetrahedron(pt, verts[1], verts[4], verts[5], verts[6]))
        return true;

    // Tet 4: vertices 3, 4, 6, 7
    if (PointInTetrahedron(pt, verts[3], verts[4], verts[6], verts[7]))
        return true;

    // Tet 5: vertices 1, 3, 4, 6 (central)
    if (PointInTetrahedron(pt, verts[1], verts[3], verts[4], verts[6]))
        return true;

    return false;
}

//-----------------------------------------------------------------------------
// Compute element AABB
//-----------------------------------------------------------------------------
void ComputeElementAABB(const TVector3d* vertices, int n_verts,
                        TVector3d& aabb_min, TVector3d& aabb_max)
{
    aabb_min = vertices[0];
    aabb_max = vertices[0];

    for (int i = 1; i < n_verts; ++i) {
        aabb_min.x = std::min(aabb_min.x, vertices[i].x);
        aabb_min.y = std::min(aabb_min.y, vertices[i].y);
        aabb_min.z = std::min(aabb_min.z, vertices[i].z);
        aabb_max.x = std::max(aabb_max.x, vertices[i].x);
        aabb_max.y = std::max(aabb_max.y, vertices[i].y);
        aabb_max.z = std::max(aabb_max.z, vertices[i].z);
    }
}

//-----------------------------------------------------------------------------
// Compute element characteristic size from AABB
//-----------------------------------------------------------------------------
double ComputeElementSize(const TVector3d* vertices, int n_verts)
{
    TVector3d aabb_min, aabb_max;
    ComputeElementAABB(vertices, n_verts, aabb_min, aabb_max);

    double dx = aabb_max.x - aabb_min.x;
    double dy = aabb_max.y - aabb_min.y;
    double dz = aabb_max.z - aabb_min.z;

    // Cube root of AABB volume as characteristic size
    return std::cbrt(dx * dy * dz);
}

//-----------------------------------------------------------------------------
// Classify multiple points
//-----------------------------------------------------------------------------
void ClassifyPoints(int n_points,
                    const double* points,
                    const std::vector<ElementData>& elements,
                    double near_threshold,
                    std::vector<ClassifyResult>& results)
{
    results.resize(n_points);

    if (elements.empty()) {
        // No elements: all points are FAR
        for (int i = 0; i < n_points; ++i) {
            results[i].classification = FAR;
            results[i].nearest_elem_id = -1;
            results[i].distance = std::numeric_limits<double>::max();
        }
        return;
    }

    // Compute global AABB
    TVector3d global_min, global_max;
    global_min.x = global_min.y = global_min.z = std::numeric_limits<double>::max();
    global_max.x = global_max.y = global_max.z = std::numeric_limits<double>::lowest();

    double total_size = 0.0;
    for (size_t e = 0; e < elements.size(); ++e) {
        const ElementData& elem = elements[e];

        if (!elem.vertices.empty()) {
            TVector3d emin, emax;
            ComputeElementAABB(elem.vertices.data(), static_cast<int>(elem.vertices.size()), emin, emax);

            global_min.x = std::min(global_min.x, emin.x);
            global_min.y = std::min(global_min.y, emin.y);
            global_min.z = std::min(global_min.z, emin.z);
            global_max.x = std::max(global_max.x, emax.x);
            global_max.y = std::max(global_max.y, emax.y);
            global_max.z = std::max(global_max.z, emax.z);
        }

        total_size += elem.size;
    }

    double avg_size = total_size / elements.size();
    double margin = avg_size * near_threshold;

    // Classify each point
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int i = 0; i < n_points; ++i) {
        TVector3d pt;
        pt.x = points[i * 3 + 0];
        pt.y = points[i * 3 + 1];
        pt.z = points[i * 3 + 2];

        // Default: FAR
        results[i].classification = FAR;
        results[i].nearest_elem_id = -1;
        results[i].distance = std::numeric_limits<double>::max();

        // Step 1: Quick rejection using global AABB
        if (pt.x < global_min.x - margin || pt.x > global_max.x + margin ||
            pt.y < global_min.y - margin || pt.y > global_max.y + margin ||
            pt.z < global_min.z - margin || pt.z > global_max.z + margin) {
            continue;  // FAR
        }

        // Step 2: Find nearest element by center distance
        double min_dist = std::numeric_limits<double>::max();
        int best_eid = -1;

        for (size_t e = 0; e < elements.size(); ++e) {
            double dx = pt.x - elements[e].center.x;
            double dy = pt.y - elements[e].center.y;
            double dz = pt.z - elements[e].center.z;
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (dist < min_dist) {
                min_dist = dist;
                best_eid = static_cast<int>(e);
            }
        }

        if (best_eid < 0) continue;

        results[i].nearest_elem_id = best_eid;
        results[i].distance = min_dist;

        double elem_size = elements[best_eid].size;

        // Step 3: Check if inside any nearby element
        if (min_dist < elem_size * 0.6) {
            double search_radius = elem_size * 1.5;

            for (size_t e = 0; e < elements.size(); ++e) {
                double dx = pt.x - elements[e].center.x;
                double dy = pt.y - elements[e].center.y;
                double dz = pt.z - elements[e].center.z;
                double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

                if (dist > search_radius) continue;

                const ElementData& elem = elements[e];
                int nv = static_cast<int>(elem.vertices.size());

                bool is_inside = false;

                if (nv == 4 && elem.num_faces == 4) {
                    // Tetrahedron
                    is_inside = PointInTetrahedron(pt,
                        elem.vertices[0], elem.vertices[1],
                        elem.vertices[2], elem.vertices[3]);
                }
                else if (nv == 8 && elem.num_faces == 6) {
                    // Hexahedron
                    TVector3d verts[8];
                    for (int v = 0; v < 8; ++v) {
                        verts[v] = elem.vertices[v];
                    }
                    is_inside = PointInHexahedron(pt, verts);
                }

                if (is_inside) {
                    results[i].classification = INSIDE;
                    results[i].nearest_elem_id = static_cast<int>(e);
                    break;
                }
            }

            if (results[i].classification == INSIDE) continue;
        }

        // Step 4: Near/Far classification
        if (min_dist < elem_size * near_threshold) {
            results[i].classification = NEAR;
        }
        else {
            results[i].classification = FAR;
        }
    }
}

//-----------------------------------------------------------------------------
// Extract element data from radTPolyhedron
//-----------------------------------------------------------------------------
static void ExtractElementData(radTPolyhedron* poly, ElementData& data)
{
    // Get center point
    data.center = poly->CentrPoint;

    // Get number of faces
    data.num_faces = poly->AmOfFaces;

    // Get vertices
    radTVectorOfVector3d vertices;
    poly->VerticesInLocFrame(vertices, true);  // true = ensure unique vertices

    data.vertices.clear();
    for (size_t i = 0; i < vertices.size(); ++i) {
        data.vertices.push_back(vertices[i]);
    }

    // Compute size
    if (!data.vertices.empty()) {
        data.size = ComputeElementSize(data.vertices.data(), static_cast<int>(data.vertices.size()));
    } else {
        data.size = 0.0;
    }
}

//-----------------------------------------------------------------------------
// Helper: Recursively extract elements from group
//-----------------------------------------------------------------------------
static void ExtractElementsFromGroup(radTGroup* grp, std::vector<ElementData>& elements)
{
    for (radTmhg::iterator iter = grp->GroupMapOfHandlers.begin();
         iter != grp->GroupMapOfHandlers.end(); ++iter)
    {
        radTg* g = iter->second.rep;

        // Try as polyhedron
        radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g);
        if (poly) {
            ElementData data;
            ExtractElementData(poly, data);
            if (!data.vertices.empty()) {
                elements.push_back(data);
            }
            continue;
        }

        // Try as nested group
        radTGroup* nested_grp = dynamic_cast<radTGroup*>(g);
        if (nested_grp) {
            ExtractElementsFromGroup(nested_grp, elements);
        }
    }
}

//-----------------------------------------------------------------------------
// Classify points from Radia handle
//-----------------------------------------------------------------------------
void ClassifyPointsFromHandle(int n_points,
                              const double* points,
                              int container_handle,
                              double near_threshold,
                              int* classification,
                              int* nearest_elem)
{
    // Get element list from container
    std::vector<ElementData> elements;

    // Access Radia handle system
    radThg hg;
    int valid = rad.ValidateElemKey(container_handle, hg);
    if (!valid || hg.rep == nullptr) {
        // Invalid handle: all points are FAR
        for (int i = 0; i < n_points; ++i) {
            classification[i] = FAR;
            nearest_elem[i] = -1;
        }
        return;
    }

    radTg* g = hg.rep;

    // Try to get as group
    radTGroup* grp = dynamic_cast<radTGroup*>(g);
    if (grp) {
        ExtractElementsFromGroup(grp, elements);
    }
    else {
        // Try as single polyhedron
        radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g);
        if (poly) {
            ElementData data;
            ExtractElementData(poly, data);
            if (!data.vertices.empty()) {
                elements.push_back(data);
            }
        }
    }

    // Classify
    std::vector<ClassifyResult> results;
    ClassifyPoints(n_points, points, elements, near_threshold, results);

    // Copy results
    for (int i = 0; i < n_points; ++i) {
        classification[i] = static_cast<int>(results[i].classification);
        nearest_elem[i] = results[i].nearest_elem_id;
    }
}

} // namespace RadPointClassify
