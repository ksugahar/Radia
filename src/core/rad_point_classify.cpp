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

#include "rad_parallel.h"

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
    if (PointInTetrahedronSolidAngle(pt, verts[0], verts[1], verts[3], verts[4]))
        return true;

    // Tet 2: vertices 1, 2, 3, 6
    if (PointInTetrahedronSolidAngle(pt, verts[1], verts[2], verts[3], verts[6]))
        return true;

    // Tet 3: vertices 1, 4, 5, 6
    if (PointInTetrahedronSolidAngle(pt, verts[1], verts[4], verts[5], verts[6]))
        return true;

    // Tet 4: vertices 3, 4, 6, 7
    if (PointInTetrahedronSolidAngle(pt, verts[3], verts[4], verts[6], verts[7]))
        return true;

    // Tet 5: vertices 1, 3, 4, 6 (central)
    if (PointInTetrahedronSolidAngle(pt, verts[1], verts[3], verts[4], verts[6]))
        return true;

    return false;
}

//=============================================================================
// Solid Angle Method (Rigorous, ELF-compatible)
//=============================================================================

double ComputeTriangleSolidAngle(const TVector3d& obs,
                                  const TVector3d& v0,
                                  const TVector3d& v1,
                                  const TVector3d& v2)
{
    TVector3d r0, r1, r2;
    r0.x = v0.x - obs.x; r0.y = v0.y - obs.y; r0.z = v0.z - obs.z;
    r1.x = v1.x - obs.x; r1.y = v1.y - obs.y; r1.z = v1.z - obs.z;
    r2.x = v2.x - obs.x; r2.y = v2.y - obs.y; r2.z = v2.z - obs.z;

    double R0 = std::sqrt(r0.x*r0.x + r0.y*r0.y + r0.z*r0.z);
    double R1 = std::sqrt(r1.x*r1.x + r1.y*r1.y + r1.z*r1.z);
    double R2 = std::sqrt(r2.x*r2.x + r2.y*r2.y + r2.z*r2.z);

    const double eps = 1e-15;
    if (R0 < eps || R1 < eps || R2 < eps) {
        return 0.0;
    }

    TVector3d r1_cross_r2;
    r1_cross_r2.x = r1.y * r2.z - r1.z * r2.y;
    r1_cross_r2.y = r1.z * r2.x - r1.x * r2.z;
    r1_cross_r2.z = r1.x * r2.y - r1.y * r2.x;

    double numerator = r0.x * r1_cross_r2.x + r0.y * r1_cross_r2.y + r0.z * r1_cross_r2.z;
    double r0_dot_r1 = r0.x*r1.x + r0.y*r1.y + r0.z*r1.z;
    double r0_dot_r2 = r0.x*r2.x + r0.y*r2.y + r0.z*r2.z;
    double r1_dot_r2 = r1.x*r2.x + r1.y*r2.y + r1.z*r2.z;
    double denominator = R0*R1*R2 + R2*r0_dot_r1 + R1*r0_dot_r2 + R0*r1_dot_r2;

    return 2.0 * std::atan2(numerator, denominator);
}

bool PointInPolyhedronSolidAngle(const TVector3d& pt,
                                  const std::vector<TVector3d>& vertices,
                                  const std::vector<std::vector<int>>& faces)
{
    const double PI = 3.14159265358979323846;
    const double FOUR_PI = 4.0 * PI;
    const double tolerance = 0.1;

    double total_solid_angle = 0.0;

    for (size_t f = 0; f < faces.size(); ++f) {
        const std::vector<int>& face = faces[f];
        int n_verts = static_cast<int>(face.size());
        if (n_verts < 3) continue;

        const TVector3d& v0 = vertices[face[0]];
        for (int i = 1; i < n_verts - 1; ++i) {
            const TVector3d& v1 = vertices[face[i]];
            const TVector3d& v2 = vertices[face[i + 1]];
            total_solid_angle += ComputeTriangleSolidAngle(pt, v0, v1, v2);
        }
    }

    // Check for both winding conventions:
    // - Inward normals: sum = -4*pi for inside
    // - Outward normals: sum = +4*pi for inside
    return std::fabs(std::fabs(total_solid_angle) - FOUR_PI) < FOUR_PI * tolerance;
}

bool PointInTetrahedronSolidAngle(const TVector3d& pt,
                                   const TVector3d& v0, const TVector3d& v1,
                                   const TVector3d& v2, const TVector3d& v3)
{
    const double PI = 3.14159265358979323846;
    const double FOUR_PI = 4.0 * PI;
    const double tolerance = 0.1;

    double total_solid_angle = 0.0;
    total_solid_angle += ComputeTriangleSolidAngle(pt, v0, v2, v1);
    total_solid_angle += ComputeTriangleSolidAngle(pt, v0, v1, v3);
    total_solid_angle += ComputeTriangleSolidAngle(pt, v0, v3, v2);
    total_solid_angle += ComputeTriangleSolidAngle(pt, v1, v2, v3);

    return std::fabs(total_solid_angle + FOUR_PI) < FOUR_PI * tolerance;
}

bool PointInHexahedronSolidAngle(const TVector3d& pt, const TVector3d verts[8])
{
    const double PI = 3.14159265358979323846;
    const double FOUR_PI = 4.0 * PI;
    const double tolerance = 0.1;

    double total_solid_angle = 0.0;

    // Face 0 (bottom): 0, 1, 2, 3
    total_solid_angle += ComputeTriangleSolidAngle(pt, verts[0], verts[1], verts[2]);
    total_solid_angle += ComputeTriangleSolidAngle(pt, verts[0], verts[2], verts[3]);
    // Face 1 (top): 4, 7, 6, 5
    total_solid_angle += ComputeTriangleSolidAngle(pt, verts[4], verts[7], verts[6]);
    total_solid_angle += ComputeTriangleSolidAngle(pt, verts[4], verts[6], verts[5]);
    // Face 2 (front): 0, 4, 5, 1
    total_solid_angle += ComputeTriangleSolidAngle(pt, verts[0], verts[4], verts[5]);
    total_solid_angle += ComputeTriangleSolidAngle(pt, verts[0], verts[5], verts[1]);
    // Face 3 (back): 2, 6, 7, 3
    total_solid_angle += ComputeTriangleSolidAngle(pt, verts[2], verts[6], verts[7]);
    total_solid_angle += ComputeTriangleSolidAngle(pt, verts[2], verts[7], verts[3]);
    // Face 4 (left): 0, 3, 7, 4
    total_solid_angle += ComputeTriangleSolidAngle(pt, verts[0], verts[3], verts[7]);
    total_solid_angle += ComputeTriangleSolidAngle(pt, verts[0], verts[7], verts[4]);
    // Face 5 (right): 1, 5, 6, 2
    total_solid_angle += ComputeTriangleSolidAngle(pt, verts[1], verts[5], verts[6]);
    total_solid_angle += ComputeTriangleSolidAngle(pt, verts[1], verts[6], verts[2]);

    return std::fabs(total_solid_angle + FOUR_PI) < FOUR_PI * tolerance;
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
    ngcore::ParallelFor(ngcore::IntRange(n_points), [&](size_t i) {
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
            return;  // FAR
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

        if (best_eid < 0) return;

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

                // Use general solid angle method with actual face topology
                bool is_inside = false;

                if (!elem.faces.empty() && !elem.vertices.empty()) {
                    is_inside = PointInPolyhedronSolidAngle(pt, elem.vertices, elem.faces);
                }

                if (is_inside) {
                    results[i].classification = INSIDE;
                    results[i].nearest_elem_id = static_cast<int>(e);
                    break;
                }
            }

            if (results[i].classification == INSIDE) return;
        }

        // Step 4: Near/Far classification
        if (min_dist < elem_size * near_threshold) {
            results[i].classification = NEAR;
        }
        else {
            results[i].classification = FAR;
        }
    });
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

    // Collect all vertices and build face topology
    // We need to track vertex indices properly for the solid angle method
    data.vertices.clear();
    data.faces.clear();

    const double AbsTol = 1e-10;
    const double AbsTolE2 = AbsTol * AbsTol;

    // Helper lambda to find or add a vertex
    auto findOrAddVertex = [&data, AbsTolE2](const TVector3d& v) -> int {
        for (size_t i = 0; i < data.vertices.size(); ++i) {
            TVector3d dp;
            dp.x = data.vertices[i].x - v.x;
            dp.y = data.vertices[i].y - v.y;
            dp.z = data.vertices[i].z - v.z;
            if ((dp.x*dp.x + dp.y*dp.y + dp.z*dp.z) <= AbsTolE2) {
                return static_cast<int>(i);
            }
        }
        data.vertices.push_back(v);
        return static_cast<int>(data.vertices.size() - 1);
    };

    // Iterate through faces and extract vertices with correct topology
    for (int f = 0; f < poly->AmOfFaces; ++f) {
        radTHandlePgnAndTrans& hPgnTrans = poly->VectHandlePgnAndTrans[f];
        radTPolygon* pgn = hPgnTrans.PgnHndl.rep;
        radTrans* trans = hPgnTrans.TransHndl.rep;

        std::vector<int> face_indices;
        double locZ = pgn->CoordZ;

        for (auto& p2d : pgn->EdgePointsVector) {
            TVector3d p3d(p2d.x, p2d.y, locZ);
            p3d = trans->TrPoint(p3d);

            // Transform to global coordinates
            TVector3d global_vertex;
            global_vertex.x = p3d.x + data.center.x;
            global_vertex.y = p3d.y + data.center.y;
            global_vertex.z = p3d.z + data.center.z;

            int idx = findOrAddVertex(global_vertex);
            face_indices.push_back(idx);
        }

        data.faces.push_back(face_indices);
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
