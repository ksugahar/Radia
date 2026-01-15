/**
 * @file rad_field_unified.cpp
 * @brief Unified Field Computation Module Implementation
 *
 * Implements the unified field computation with:
 * 1. Point classification (inside/outside/near)
 * 2. FMM acceleration via ExaFMM
 * 3. Direct computation via B_genComp
 *
 * @version 1.6.0
 * @date 2026-01-14
 */

#include "rad_field_unified.h"
#include "rad_geometry_3d.h"
#include "rad_group.h"
#include "rad_polyhedron.h"
#include "rad_rectangular_block.h"
#include "rad_application.h"
#include "rad_point_classify.h"
#include "rad_dipole_collect.h"
#include "rad_exafmm.h"
#include "rad_type_cast.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cmath>
#include <unordered_map>
#include <mutex>
#include <functional>

// Access global Radia application instance
extern radTApplication rad;

namespace RadFieldUnified {

// Global cache for FMM dipole data per container
static std::unordered_map<int, RadDipoleCollect::DipoleCollection> g_fmmCache;
static std::mutex g_fmmCacheMutex;

// Global cache for element data per container
static std::unordered_map<int, std::vector<RadPointClassify::ElementData>> g_elemCache;
static std::mutex g_elemCacheMutex;

// Physical constants
static const double MU_0 = 4.0 * 3.14159265358979323846 * 1.0e-7; // T*m/A

//-----------------------------------------------------------------------------
// Helper: Extract element vertices from various Radia element types
// Uses the same pattern as rad_point_classify.cpp ExtractElementData()
//-----------------------------------------------------------------------------
static bool ExtractElementVertices(
    radTg3d* elem,
    std::vector<TVector3d>& vertices,
    std::vector<std::vector<int>>& faces,
    int& num_faces
)
{
    vertices.clear();
    faces.clear();
    num_faces = 0;

    if (elem == nullptr) return false;

    // Try rectangular block first
    radTRecMag* recMag = dynamic_cast<radTRecMag*>(elem);
    if (recMag) {
        // Rectangular block: get 8 corner vertices
        TVector3d center = recMag->CentrPoint;
        TVector3d dims = recMag->Dimensions;

        double hx = dims.x / 2.0;
        double hy = dims.y / 2.0;
        double hz = dims.z / 2.0;

        vertices.resize(8);
        vertices[0] = TVector3d(center.x - hx, center.y - hy, center.z - hz);
        vertices[1] = TVector3d(center.x + hx, center.y - hy, center.z - hz);
        vertices[2] = TVector3d(center.x + hx, center.y + hy, center.z - hz);
        vertices[3] = TVector3d(center.x - hx, center.y + hy, center.z - hz);
        vertices[4] = TVector3d(center.x - hx, center.y - hy, center.z + hz);
        vertices[5] = TVector3d(center.x + hx, center.y - hy, center.z + hz);
        vertices[6] = TVector3d(center.x + hx, center.y + hy, center.z + hz);
        vertices[7] = TVector3d(center.x - hx, center.y + hy, center.z + hz);

        // 6 faces of hexahedron (0-indexed vertices)
        faces = {
            {0, 3, 2, 1},  // -z face
            {4, 5, 6, 7},  // +z face
            {0, 1, 5, 4},  // -y face
            {2, 3, 7, 6},  // +y face
            {0, 4, 7, 3},  // -x face
            {1, 2, 6, 5}   // +x face
        };
        num_faces = 6;
        return true;
    }

    // Try polyhedron (hexahedron or tetrahedron)
    // Use the same pattern as rad_point_classify.cpp::ExtractElementData()
    radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(elem);
    if (poly) {
        const double AbsTol = 1e-10;
        const double AbsTolE2 = AbsTol * AbsTol;

        TVector3d center = poly->CentrPoint;

        // Helper lambda to find or add a vertex
        auto findOrAddVertex = [&vertices, AbsTolE2](const TVector3d& v) -> int {
            for (size_t i = 0; i < vertices.size(); ++i) {
                TVector3d dp;
                dp.x = vertices[i].x - v.x;
                dp.y = vertices[i].y - v.y;
                dp.z = vertices[i].z - v.z;
                if ((dp.x*dp.x + dp.y*dp.y + dp.z*dp.z) <= AbsTolE2) {
                    return static_cast<int>(i);
                }
            }
            vertices.push_back(v);
            return static_cast<int>(vertices.size() - 1);
        };

        // Iterate through faces and extract vertices with correct topology
        num_faces = poly->AmOfFaces;
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
                global_vertex.x = p3d.x + center.x;
                global_vertex.y = p3d.y + center.y;
                global_vertex.z = p3d.z + center.z;

                int idx = findOrAddVertex(global_vertex);
                face_indices.push_back(idx);
            }

            faces.push_back(face_indices);
        }

        return !vertices.empty();
    }

    return false;
}

//-----------------------------------------------------------------------------
// Helper: Get magnetization from element (requires radTg3dRelax or derived)
//-----------------------------------------------------------------------------
static bool GetElementMagnetization(radTg3d* elem, TVector3d& M_out)
{
    M_out = TVector3d(0, 0, 0);

    if (!elem) return false;

    // Need to cast to radTg3dRelax to access Magn member
    radTg3dRelax* relaxElem = dynamic_cast<radTg3dRelax*>(elem);
    if (relaxElem) {
        M_out = relaxElem->Magn;
        return true;
    }

    return false;
}

//-----------------------------------------------------------------------------
// BuildElementData: Extract element data from container
//-----------------------------------------------------------------------------
bool BuildElementData(
    int container_handle,
    std::vector<RadPointClassify::ElementData>& elements
)
{
    elements.clear();

    // Get application instance
    radTApplication& radApp = rad;

    // Validate handle
    radThg hg;
    if (!radApp.ValidateElemKey(container_handle, hg)) {
        return false;
    }

    radTg3d* g3dPtr = radTCast::g3dCast(hg.rep);
    if (!g3dPtr) return false;

    // Recursive function to collect elements
    std::function<void(radTg3d*)> collectElements = [&](radTg3d* elem) {
        if (!elem) return;

        // Check if it's a container
        radTGroup* group = dynamic_cast<radTGroup*>(elem);
        if (group) {
            // Iterate through children
            for (auto& child : group->GroupMapOfHandlers) {
                radTg3d* childElem = radTCast::g3dCast(child.second.rep);
                if (childElem) {
                    collectElements(childElem);
                }
            }
            return;
        }

        // Single element - extract vertices
        RadPointClassify::ElementData ed;
        if (ExtractElementVertices(elem, ed.vertices, ed.faces, ed.num_faces)) {
            // Compute center and size
            TVector3d sum(0, 0, 0);
            for (const auto& v : ed.vertices) {
                sum = sum + v;
            }
            ed.center = sum * (1.0 / ed.vertices.size());
            ed.size = RadPointClassify::ComputeElementSize(
                ed.vertices.data(), (int)ed.vertices.size());
            elements.push_back(ed);
        }
    };

    collectElements(g3dPtr);
    return !elements.empty();
}

//-----------------------------------------------------------------------------
// IsPointInsideMagnet: Check if point is inside any element
//-----------------------------------------------------------------------------
bool IsPointInsideMagnet(
    radTg3d* g3dPtr,
    const TVector3d& point,
    int& element_index
)
{
    element_index = -1;

    if (!g3dPtr) return false;

    int current_index = 0;

    // Recursive function to check all elements
    std::function<bool(radTg3d*)> checkElement = [&](radTg3d* elem) -> bool {
        if (!elem) return false;

        // Check if it's a container
        radTGroup* group = dynamic_cast<radTGroup*>(elem);
        if (group) {
            for (auto& child : group->GroupMapOfHandlers) {
                radTg3d* childElem = radTCast::g3dCast(child.second.rep);
                if (childElem && checkElement(childElem)) {
                    return true;
                }
            }
            return false;
        }

        // Single element - check inside
        std::vector<TVector3d> vertices;
        std::vector<std::vector<int>> faces;
        int num_faces;

        if (!ExtractElementVertices(elem, vertices, faces, num_faces)) {
            current_index++;
            return false;
        }

        // Use solid angle method
        if (RadPointClassify::PointInPolyhedronSolidAngle(point, vertices, faces)) {
            element_index = current_index;
            return true;
        }

        current_index++;
        return false;
    };

    return checkElement(g3dPtr);
}

//-----------------------------------------------------------------------------
// GetMagnetizationAtPoint: Get M inside an element
//-----------------------------------------------------------------------------
bool GetMagnetizationAtPoint(
    radTg3d* g3dPtr,
    const TVector3d& point,
    int element_index,
    TVector3d& M_out
)
{
    M_out = TVector3d(0, 0, 0);

    if (!g3dPtr || element_index < 0) return false;

    int current_index = 0;

    // Recursive search for element by index
    std::function<radTg3d*(radTg3d*)> findElement = [&](radTg3d* elem) -> radTg3d* {
        if (!elem) return nullptr;

        // Check if it's a container
        radTGroup* group = dynamic_cast<radTGroup*>(elem);
        if (group) {
            for (auto& child : group->GroupMapOfHandlers) {
                radTg3d* found = findElement(radTCast::g3dCast(child.second.rep));
                if (found) return found;
            }
            return nullptr;
        }

        // Single element
        if (current_index == element_index) {
            return elem;
        }
        current_index++;
        return nullptr;
    };

    radTg3d* elem = findElement(g3dPtr);
    if (!elem) return false;

    // Get magnetization from element
    return GetElementMagnetization(elem, M_out);
}

//-----------------------------------------------------------------------------
// ComputeFieldForTrajectory: B field for particle tracking
//-----------------------------------------------------------------------------
bool ComputeFieldForTrajectory(
    radTg3d* g3dPtr,
    const TVector3d& point,
    TVector3d& B_out
)
{
    B_out = TVector3d(0, 0, 0);

    if (!g3dPtr) return false;

    // Check if point is inside any magnet
    int element_index = -1;
    bool inside = IsPointInsideMagnet(g3dPtr, point, element_index);

    if (inside) {
        // Point is inside magnet - return false to signal invalid trajectory point
        // The trajectory integrator should handle this (e.g., skip or warn)
        return false;
    }

    // Point is outside - compute field normally
    TVector3d ZeroVect(0, 0, 0);
    radTFieldKey FieldKey;
    FieldKey.B_ = true;

    radTField Field(FieldKey, ZeroVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect);
    Field.P = point;

    g3dPtr->B_genComp(&Field);

    B_out = Field.B;
    return true;
}

//-----------------------------------------------------------------------------
// ComputeFieldSingle: Single point field computation
//-----------------------------------------------------------------------------
FieldResult ComputeFieldSingle(
    radTg3d* g3dPtr,
    const TVector3d& point,
    FieldType field_type,
    const ComputeConfig& config
)
{
    FieldResult result;
    result.Bx = result.By = result.Bz = 0.0;
    result.Hx = result.Hy = result.Hz = 0.0;
    result.status = STATUS_ERROR;
    result.element_id = -1;

    if (!g3dPtr) return result;

    // Check inside/outside if enabled
    if (config.check_inside) {
        int elem_idx = -1;
        if (IsPointInsideMagnet(g3dPtr, point, elem_idx)) {
            result.status = STATUS_INSIDE;
            result.element_id = elem_idx;

            if (config.return_internal_field) {
                // Return internal magnetization as B field
                TVector3d M;
                if (GetMagnetizationAtPoint(g3dPtr, point, elem_idx, M)) {
                    // Inside permanent magnet: B = mu0 * M (approximately)
                    result.Bx = MU_0 * M.x;
                    result.By = MU_0 * M.y;
                    result.Bz = MU_0 * M.z;
                    // H inside: H = B/mu0 - M = 0 (for ideal PM)
                    result.Hx = result.Hy = result.Hz = 0.0;
                }
            }
            return result;
        }
    }

    // Point is outside - compute field
    TVector3d ZeroVect(0, 0, 0);
    radTFieldKey FieldKey;

    switch (field_type) {
        case FIELD_B:
            FieldKey.B_ = true;
            break;
        case FIELD_H:
            FieldKey.H_ = true;
            break;
        case FIELD_A:
            FieldKey.A_ = true;
            break;
        case FIELD_PHI:
            FieldKey.Phi_ = true;
            break;
        default:
            FieldKey.B_ = true;
    }

    radTField Field(FieldKey, ZeroVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect);
    Field.P = point;

    g3dPtr->B_genComp(&Field);

    result.Bx = Field.B.x;
    result.By = Field.B.y;
    result.Bz = Field.B.z;
    result.Hx = Field.H.x;
    result.Hy = Field.H.y;
    result.Hz = Field.H.z;
    result.status = STATUS_OUTSIDE;

    return result;
}

//-----------------------------------------------------------------------------
// ComputeFieldBatch: Batch field computation with OpenMP
//-----------------------------------------------------------------------------
void ComputeFieldBatch(
    radTg3d* g3dPtr,
    const double* points,
    int n_points,
    FieldType field_type,
    const ComputeConfig& config,
    double* B_out,
    double* H_out,
    PointStatus* status_out
)
{
    if (!g3dPtr || n_points <= 0) return;

    // Initialize outputs
    if (B_out) std::memset(B_out, 0, n_points * 3 * sizeof(double));
    if (H_out) std::memset(H_out, 0, n_points * 3 * sizeof(double));
    if (status_out) {
        for (int i = 0; i < n_points; i++) {
            status_out[i] = STATUS_OUTSIDE;
        }
    }

    // Build element data for classification if needed
    std::vector<RadPointClassify::ElementData> elements;
    if (config.check_inside) {
        // We need to build element data for inside check
        // For now, extract elements directly (could cache in future)
        std::function<void(radTg3d*)> collectElements = [&](radTg3d* elem) {
            if (!elem) return;

            radTGroup* group = dynamic_cast<radTGroup*>(elem);
            if (group) {
                for (auto& child : group->GroupMapOfHandlers) {
                    radTg3d* childElem = radTCast::g3dCast(child.second.rep);
                    if (childElem) {
                        collectElements(childElem);
                    }
                }
                return;
            }

            RadPointClassify::ElementData ed;
            if (ExtractElementVertices(elem, ed.vertices, ed.faces, ed.num_faces)) {
                TVector3d sum(0, 0, 0);
                for (const auto& v : ed.vertices) {
                    sum = sum + v;
                }
                ed.center = sum * (1.0 / ed.vertices.size());
                ed.size = RadPointClassify::ComputeElementSize(
                    ed.vertices.data(), (int)ed.vertices.size());
                elements.push_back(ed);
            }
        };

        collectElements(g3dPtr);
    }

    // Setup field key
    radTFieldKey FieldKey;
    switch (field_type) {
        case FIELD_B: FieldKey.B_ = true; break;
        case FIELD_H: FieldKey.H_ = true; break;
        case FIELD_A: FieldKey.A_ = true; break;
        case FIELD_PHI: FieldKey.Phi_ = true; break;
        default: FieldKey.B_ = true;
    }

    TVector3d ZeroVect(0, 0, 0);

    // OpenMP parallelization
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(n_points > 100)
    #endif
    for (int i = 0; i < n_points; i++) {
        TVector3d pt(points[i*3], points[i*3+1], points[i*3+2]);

        // Check inside if enabled
        if (config.check_inside && !elements.empty()) {
            // Check if point is inside any element
            bool is_inside = false;
            for (size_t e = 0; e < elements.size(); e++) {
                if (RadPointClassify::PointInPolyhedronSolidAngle(
                        pt, elements[e].vertices, elements[e].faces)) {
                    if (status_out) status_out[i] = STATUS_INSIDE;
                    is_inside = true;

                    if (config.return_internal_field) {
                        // For now, use zero field inside
                        // TODO: Get actual magnetization
                    }
                    break;
                }
            }
            if (is_inside) continue;  // Skip to next point
        }

        // Compute field (point is outside)
        radTField Field(FieldKey, ZeroVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect);
        Field.P = pt;

        g3dPtr->B_genComp(&Field);

        if (B_out) {
            B_out[i*3 + 0] = Field.B.x;
            B_out[i*3 + 1] = Field.B.y;
            B_out[i*3 + 2] = Field.B.z;
        }
        if (H_out) {
            H_out[i*3 + 0] = Field.H.x;
            H_out[i*3 + 1] = Field.H.y;
            H_out[i*3 + 2] = Field.H.z;
        }
    }
}

//-----------------------------------------------------------------------------
// InitializeFMM: Setup FMM dipole data
//-----------------------------------------------------------------------------
bool InitializeFMM(int container_handle)
{
    std::lock_guard<std::mutex> lock(g_fmmCacheMutex);

    // Check if already initialized
    if (g_fmmCache.find(container_handle) != g_fmmCache.end()) {
        return true;
    }

    // Extract dipoles using RadDipoleCollect
    RadDipoleCollect::DipoleCollection dipoles;

    // Get application instance
    radTApplication& radApp = rad;

    radThg hg;
    if (!radApp.ValidateElemKey(container_handle, hg)) {
        return false;
    }

    radTg3d* g3dPtr = radTCast::g3dCast(hg.rep);
    if (!g3dPtr) return false;

    // Recursive extraction
    std::function<void(radTg3d*)> extractDipoles = [&](radTg3d* elem) {
        if (!elem) return;

        // Container
        radTGroup* group = dynamic_cast<radTGroup*>(elem);
        if (group) {
            for (auto& child : group->GroupMapOfHandlers) {
                extractDipoles(radTCast::g3dCast(child.second.rep));
            }
            return;
        }

        // Single element
        std::vector<TVector3d> vertices;
        std::vector<std::vector<int>> faces;
        int num_faces;

        if (!ExtractElementVertices(elem, vertices, faces, num_faces)) {
            return;
        }

        // Get magnetization
        TVector3d M;
        if (!GetElementMagnetization(elem, M)) {
            return;  // No magnetization, skip
        }

        // Compute center and volume
        TVector3d center(0, 0, 0);
        for (const auto& v : vertices) {
            center = center + v;
        }
        center = center * (1.0 / vertices.size());

        // Estimate volume (simplified using AABB)
        TVector3d aabb_min, aabb_max;
        RadPointClassify::ComputeElementAABB(vertices.data(), (int)vertices.size(),
                                             aabb_min, aabb_max);
        double vol = (aabb_max.x - aabb_min.x) *
                     (aabb_max.y - aabb_min.y) *
                     (aabb_max.z - aabb_min.z);

        if (vol > 0.0) {
            RadDipoleCollect::DipoleData d;
            d.x = center.x;
            d.y = center.y;
            d.z = center.z;
            d.mx = M.x * vol;
            d.my = M.y * vol;
            d.mz = M.z * vol;
            d.volume = vol;
            dipoles.dipoles.push_back(d);
        }
    };

    extractDipoles(g3dPtr);
    dipoles.flatten();

    if (dipoles.count() > 0) {
        g_fmmCache[container_handle] = std::move(dipoles);
        return true;
    }

    return false;
}

//-----------------------------------------------------------------------------
// ReleaseFMM: Release FMM data
//-----------------------------------------------------------------------------
void ReleaseFMM(int container_handle)
{
    std::lock_guard<std::mutex> lock(g_fmmCacheMutex);
    g_fmmCache.erase(container_handle);

    std::lock_guard<std::mutex> lock2(g_elemCacheMutex);
    g_elemCache.erase(container_handle);
}

} // namespace RadFieldUnified
