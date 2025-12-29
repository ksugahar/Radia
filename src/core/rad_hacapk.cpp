/*-------------------------------------------------------------------------
*
* File name:      rad_hacapk.cpp
*
* Project:        RADIA
*
* Description:    HACApK (H-matrix with ACA+) interface for BiCGSTAB solver
*                 Implementation of RadHACApKManager and callback functions
*
* First release:  2025
*
* Reference:      ppOpen-HPC project (MIT License)
*                 https://github.com/ppohHPC/ppOpen-HPC
*
-------------------------------------------------------------------------*/

#include "rad_hacapk.h"
#include "rad_interaction.h"
#include "rad_constants.h"
#include <cmath>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <chrono>

// Include C++ compatible HACApK wrapper header
extern "C" {
#include "../ext/HACApK/cHACApK_cpp.h"
}

//=========================================================================
// Constants for 6DOF MSC field computation
// Now using unified constants from rad_constants.h
//=========================================================================

//=========================================================================
// Global callback state (required by HACApK C interface)
//=========================================================================

namespace {
    // Thread-local storage for callback state
    RadHACApKManager* g_currentManager = nullptr;
    std::vector<double> g_invChi;
    radTInteraction* g_interaction = nullptr;
    int g_nElem = 0;
    int g_nffc = 3;  // DOF per element (default 3 for standard elements)
    // Note: lod (DOF permutation array) is now accessed via C-side accessors:
    //   HACApK_get_current_lod() and HACApK_get_current_lod_size()
    // Set during H-matrix build by HACApK_build_hmatrix_varDOF_wrapper
}

namespace RadHACApKCallback {

void SetCurrentManager(RadHACApKManager* manager) {
    g_currentManager = manager;
}

RadHACApKManager* GetCurrentManager() {
    return g_currentManager;
}

void SetInvChi(const std::vector<double>& inv_chi) {
    g_invChi = inv_chi;
}

const std::vector<double>& GetInvChi() {
    return g_invChi;
}

void SetInteraction(radTInteraction* interaction, int n_elem, int nffc) {
    g_interaction = interaction;
    g_nElem = n_elem;
    g_nffc = nffc;
}

void SetLod(int* lod, int size) {
    // Deprecated: lod is now set/cleared on the C side by HACApK_build_hmatrix_varDOF_wrapper
    // This function is kept for API compatibility but does nothing
    (void)lod;
    (void)size;
}

void ClearLod() {
    // Deprecated: lod is now cleared on the C side by HACApK_build_hmatrix_varDOF_wrapper
    // This function is kept for API compatibility but does nothing
}

double ComputeEntry(int i, int j) {
    // i, j are 1-based ORIGINAL indices from HACApK
    // The lod conversion is already done by cHACApK_fill_leafmtx_hyp:
    //   val = cHACApK_entry_ij(lodl[permuted_pos], lodt[permuted_pos], i_bemv)
    // So we receive original indices, NOT permuted indices!
    //
    // Matrix element: A(i,j) = N(i,j) + delta_ij/chi_i
    // where N already contains -K/(4*pi) from GetInteractionMatrixElement()
    // So the equation is: (-K/(4pi) + 1/chi * I) * sigma = H_ext_n

    if (g_currentManager == nullptr) {
        std::cerr << "[HACApK] Error: g_currentManager is null in ComputeEntry" << std::endl;
        return 0.0;
    }

    int ndof = g_currentManager->GetNDOF();

    // Direct 0-based conversion (indices are already original, NOT permuted)
    int i0 = i - 1;
    int j0 = j - 1;

    // Bounds check on original indices
    if (i0 < 0 || i0 >= ndof || j0 < 0 || j0 >= ndof) {
        std::cerr << "[HACApK] Error: Invalid DOF indices: i0=" << i0
                  << " j0=" << j0 << " ndof=" << ndof << std::endl;
        return 0.0;
    }

    // Get N matrix element through manager (friend class access)
    // GetInteractionMatrixElement() returns -K_ij/(4*pi) for MSC hexahedra
    double N_val = g_currentManager->GetInteractionMatrixElement(i0, j0);

    // A = N + diag(1/chi) (N already has correct sign: -K/(4pi))
    double A_val = N_val;
    if (i0 == j0 && i0 < (int)g_invChi.size()) {
        A_val += g_invChi[i0];
    }

    return A_val;
}

}  // namespace RadHACApKCallback

//=========================================================================
// C callback function for HACApK
// This is the function HACApK calls to get matrix elements
//=========================================================================

extern "C" {

double cHACApK_entry_ij(int i, int j, int i_bemv) {
    (void)i_bemv;  // Unused in Radia
    return RadHACApKCallback::ComputeEntry(i, j);
}

}  // extern "C"

//=========================================================================
// RadHACApKManager Implementation
//=========================================================================

RadHACApKManager::RadHACApKManager(radTInteraction* interaction)
    : m_interaction(interaction)
    , m_leafmtxp(nullptr)
    , m_control(nullptr)
    , m_valid(false)
    , m_ndof(0)
    , m_n_elem(0)
    , m_nffc(3)
    , m_is_6dof(false)
    , m_is_mixed_dof(false)
    , m_tri_precomputed(false)
    , m_geometry_ready(false)
    , m_geometry_3dof_ready(false)
    , m_diag_cached(false)
    , m_flat_N_ready(false)
{
    // Hash-based cache is initialized automatically
}

RadHACApKManager::~RadHACApKManager() {
    FreeResources();
}

void RadHACApKManager::FreeResources() {
    if (m_leafmtxp || m_control) {
        HACApK_free_hmatrix_wrapper(m_leafmtxp, m_control);
    }
    if (m_leafmtxp) {
        HACApK_free_leafmtxp(m_leafmtxp);
        m_leafmtxp = nullptr;
    }
    if (m_control) {
        HACApK_free_lcontrol(m_control);
        m_control = nullptr;
    }
    m_valid = false;
}

void RadHACApKManager::ExtractElementCoordinates() {
    // Extract element center coordinates for clustering
    // Supports both 3DOF tetrahedra and 6DOF MSC hexahedra
    if (!m_interaction) return;

    m_n_elem = m_interaction->AmOfMainElem;
    m_coordinates.resize(m_n_elem * 3);
    m_dof_offset.resize(m_n_elem + 1);

    // Check if using variable DOF (6DOF hexahedra) or uniform 3DOF (tetrahedra)
    bool has_variable_dof = m_interaction->HasVariableDOF();

    int total_dof = 0;
    int n_3dof = 0;
    int n_6dof = 0;

    for (int i = 0; i < m_n_elem; i++) {
        m_dof_offset[i] = total_dof;
        int elem_dof = m_interaction->GetElementDOF(i);
        total_dof += elem_dof;

        if (elem_dof == 3) {
            n_3dof++;
        } else if (elem_dof == 6) {
            n_6dof++;
        } else {
            std::cerr << "[HACApK] Error: Element " << i << " has " << elem_dof
                      << " DOF, expected 3 or 6" << std::endl;
            m_ndof = 0;
            m_nffc = 0;
            m_is_6dof = false;
            return;
        }
    }

    m_dof_offset[m_n_elem] = total_dof;
    m_ndof = total_dof;

    // Support mixed DOF elements (hex + tetra)
    if (n_3dof > 0 && n_6dof > 0) {
        // Mixed mode: variable DOF per element
        m_nffc = 0;  // Indicates variable DOF
        m_is_6dof = false;  // Not uniform 6DOF
        m_is_mixed_dof = true;
#ifdef HACAPK_RADIA_LOGGING
        std::cout << "[HACApK] Mixed DOF mode: " << n_3dof << " tetrahedra (3DOF) + "
                  << n_6dof << " hexahedra (6DOF), total " << total_dof << " DOF" << std::endl;
#endif
    } else if (n_6dof > 0) {
        m_nffc = 6;
        m_is_6dof = true;
        m_is_mixed_dof = false;
    } else {
        m_nffc = 3;
        m_is_6dof = false;
        m_is_mixed_dof = false;
    }

    // Get element centers from g3dRelaxPtrVect
    for (int i = 0; i < m_n_elem; i++) {
        radTg3dRelax* elem = m_interaction->g3dRelaxPtrVect[i];
        if (elem) {
            TVector3d center = elem->CentrPoint;
            m_coordinates[i * 3 + 0] = center.x;
            m_coordinates[i * 3 + 1] = center.y;
            m_coordinates[i * 3 + 2] = center.z;
        }
    }

    // Build O(1) DOF lookup table
    BuildDOFLookupTable();
}

//=========================================================================
// BuildDOFLookupTable: Create O(1) DOF-to-element lookup (ELF-style)
//=========================================================================

void RadHACApKManager::BuildDOFLookupTable() {
    if (m_ndof == 0) return;

    m_dof_to_elem.resize(m_ndof);
    m_dof_to_local.resize(m_ndof);

    for (int e = 0; e < m_n_elem; e++) {
        int offset = m_dof_offset[e];
        int elem_dof = m_dof_offset[e + 1] - offset;
        for (int d = 0; d < elem_dof; d++) {
            m_dof_to_elem[offset + d] = e;
            m_dof_to_local[offset + d] = d;
        }
    }
}

//=========================================================================
// PrecomputeGeometry: ELF-style pre-computed geometry for fast matrix access
// Extracts all hexahedron vertices and face geometry into contiguous arrays
//=========================================================================

void RadHACApKManager::PrecomputeGeometry() {
    if (m_geometry_ready || m_n_elem == 0 || !m_interaction) return;

    // Allocate arrays
    m_elem_centers.resize(m_n_elem * 3);
    m_face_centers.resize(m_n_elem * 6 * 3);
    m_face_normals.resize(m_n_elem * 6 * 3);
    m_face_areas.resize(m_n_elem * 6);
    m_face_vertices.resize(m_n_elem * 6 * 4 * 3);  // 6 faces, 4 vertices each, 3 coords

    for (int e = 0; e < m_n_elem; e++) {
        radTg3dRelax* elem = m_interaction->g3dRelaxPtrVect[e];
        if (!elem) continue;

        radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(elem);
        if (!poly || !poly->Use6DOF_MSC) continue;

        // Store element center
        int cIdx = e * 3;
        m_elem_centers[cIdx + 0] = poly->CentrPoint.x;
        m_elem_centers[cIdx + 1] = poly->CentrPoint.y;
        m_elem_centers[cIdx + 2] = poly->CentrPoint.z;

        // Store face data for each of the 6 faces
        for (int f = 0; f < 6; f++) {
            int fcIdx = (e * 6 + f) * 3;
            int aIdx = e * 6 + f;

            // Face center
            m_face_centers[fcIdx + 0] = poly->FaceCenter[f].x;
            m_face_centers[fcIdx + 1] = poly->FaceCenter[f].y;
            m_face_centers[fcIdx + 2] = poly->FaceCenter[f].z;

            // Face normal
            m_face_normals[fcIdx + 0] = poly->FaceNormal[f].x;
            m_face_normals[fcIdx + 1] = poly->FaceNormal[f].y;
            m_face_normals[fcIdx + 2] = poly->FaceNormal[f].z;

            // Face area
            m_face_areas[aIdx] = poly->FaceArea[f];

            // Get 4 vertices of this face from polygon
            radTHandlePgnAndTrans hpt = poly->VectHandlePgnAndTrans[f];
            radTPolygon* pgn = hpt.PgnHndl.rep;
            radTrans* tr = hpt.TransHndl.rep;

            const radTVect2dVect& verts2d = pgn->EdgePointsVector;
            int fvIdx = (e * 6 + f) * 4 * 3;  // Starting index for this face's vertices

            for (int v = 0; v < 4 && v < (int)verts2d.size(); v++) {
                TVector3d V = tr->TrPoint(TVector3d(verts2d[v].x, verts2d[v].y, pgn->CoordZ));
                m_face_vertices[fvIdx + v * 3 + 0] = V.x;
                m_face_vertices[fvIdx + v * 3 + 1] = V.y;
                m_face_vertices[fvIdx + v * 3 + 2] = V.z;
            }
        }
    }

    m_geometry_ready = true;

    // Also pre-compute triangle local coordinate systems
    PrecomputeTriangleData();
}

//=========================================================================
// PrecomputeTriangleData: Pre-compute triangle local coordinate systems
// Eliminates redundant sqrt/div operations during field computation
// Each hexahedron face is split into 2 triangles = 12 triangles per element
//=========================================================================

void RadHACApKManager::PrecomputeTriangleData() {
    if (m_tri_precomputed || m_n_elem == 0 || !m_geometry_ready) return;

    const double EPS = 1.0e-20;

    // Allocate: 12 triangles per element, 32 doubles per triangle
    m_tri_data.resize(m_n_elem * TRIS_PER_ELEM * TRI_DATA_SIZE);

    // Triangle split indices for quad face: [0,1,2], [0,2,3]
    const int tri_split[2][3] = {{0, 1, 2}, {0, 2, 3}};

    for (int e = 0; e < m_n_elem; e++) {
        const double* center = &m_elem_centers[e * 3];

        for (int f = 0; f < 6; f++) {
            // Get quad vertices
            int fvIdx = (e * 6 + f) * 4 * 3;
            const double* V[4];
            for (int v = 0; v < 4; v++) {
                V[v] = &m_face_vertices[fvIdx + v * 3];
            }

            for (int t = 0; t < 2; t++) {
                int tri_idx = e * TRIS_PER_ELEM + f * 2 + t;
                double* data = &m_tri_data[tri_idx * TRI_DATA_SIZE];

                // Get triangle vertices
                const double* v0 = V[tri_split[t][0]];
                const double* v1 = V[tri_split[t][1]];
                const double* v2 = V[tri_split[t][2]];

                // Build local coordinate system
                double e1[3] = {v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]};
                double e2[3] = {v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]};

                // basis_c = e1 x e2 (face normal)
                double* basis_c = data + 6;
                basis_c[0] = e1[1]*e2[2] - e1[2]*e2[1];
                basis_c[1] = e1[2]*e2[0] - e1[0]*e2[2];
                basis_c[2] = e1[0]*e2[1] - e1[1]*e2[0];

                double cLen = sqrt(basis_c[0]*basis_c[0] + basis_c[1]*basis_c[1] + basis_c[2]*basis_c[2]);
                if (cLen < EPS) {
                    std::memset(data, 0, TRI_DATA_SIZE * sizeof(double));
                    continue;
                }
                basis_c[0] /= cLen; basis_c[1] /= cLen; basis_c[2] /= cLen;

                // basis_a = e1 normalized
                double* basis_a = data;
                basis_a[0] = e1[0]; basis_a[1] = e1[1]; basis_a[2] = e1[2];
                double aLen = sqrt(basis_a[0]*basis_a[0] + basis_a[1]*basis_a[1] + basis_a[2]*basis_a[2]);
                if (aLen < EPS) {
                    std::memset(data, 0, TRI_DATA_SIZE * sizeof(double));
                    continue;
                }
                basis_a[0] /= aLen; basis_a[1] /= aLen; basis_a[2] /= aLen;

                // basis_b = basis_c x basis_a
                double* basis_b = data + 3;
                basis_b[0] = basis_c[1]*basis_a[2] - basis_c[2]*basis_a[1];
                basis_b[1] = basis_c[2]*basis_a[0] - basis_c[0]*basis_a[2];
                basis_b[2] = basis_c[0]*basis_a[1] - basis_c[1]*basis_a[0];

                // origin = v0
                double* origin = data + 9;
                origin[0] = v0[0]; origin[1] = v0[1]; origin[2] = v0[2];

                // 2D coordinates (v0 = origin)
                double* XY = data + 12;  // 6 doubles: XY[0..5] = {x0,y0, x1,y1, x2,y2}
                XY[0] = 0.0; XY[1] = 0.0;  // v0
                XY[2] = e1[0]*basis_a[0] + e1[1]*basis_a[1] + e1[2]*basis_a[2];
                XY[3] = e1[0]*basis_b[0] + e1[1]*basis_b[1] + e1[2]*basis_b[2];
                XY[4] = e2[0]*basis_a[0] + e2[1]*basis_a[1] + e2[2]*basis_a[2];
                XY[5] = e2[0]*basis_b[0] + e2[1]*basis_b[1] + e2[2]*basis_b[2];

                // Edge parameters
                double* DS = data + 18;  // 3 doubles
                double* AM = data + 21;  // 3 doubles
                double* XD = data + 24;  // 3 doubles
                double* YD = data + 27;  // 3 doubles
                double EPSG = 0.0;

                for (int j = 0; j < 3; j++) {
                    int l = (j + 1) % 3;
                    double dx = XY[l*2] - XY[j*2];
                    double dy = XY[l*2+1] - XY[j*2+1];
                    if (fabs(dx) < EPS) dx = (dx >= 0) ? EPS : -EPS;

                    DS[j] = sqrt(dx*dx + dy*dy);
                    AM[j] = dy / dx;
                    XD[j] = -dx / DS[j];
                    YD[j] = dy / DS[j];

                    if (DS[j] > EPSG) EPSG = DS[j];
                }

                // Store EPSG and sign
                data[30] = EPSG * 1.0e-12;  // EPSG

                // Normal orientation sign (outward from center)
                double tc[3] = {(v0[0]+v1[0]+v2[0])/3.0 - center[0],
                                (v0[1]+v1[1]+v2[1])/3.0 - center[1],
                                (v0[2]+v1[2]+v2[2])/3.0 - center[2]};
                double dot = basis_c[0]*tc[0] + basis_c[1]*tc[1] + basis_c[2]*tc[2];
                data[31] = (dot >= 0.0) ? 1.0 : -1.0;  // sign
            }
        }
    }

    m_tri_precomputed = true;
}

//=========================================================================
// PrecomputeGeometry3DOF: ELF-style pre-computed geometry for 3DOF tetrahedra
// Extracts all tetrahedron vertices and face geometry into contiguous arrays
// This avoids calling B_comp() which has significant overhead during H-matrix build
//=========================================================================

void RadHACApKManager::PrecomputeGeometry3DOF() {
    if (m_geometry_3dof_ready || m_n_elem == 0 || !m_interaction) return;

    // Allocate arrays for tetrahedra (4 triangular faces, 3 vertices each)
    m_tetra_centers.resize(m_n_elem * 3);
    m_tetra_face_vertices.resize(m_n_elem * 4 * 3 * 3);  // 4 faces, 3 vertices, 3 coords
    m_tetra_face_normals.resize(m_n_elem * 4 * 3);       // 4 faces, 3 coords (outward normals)
    m_tetra_face_areas.resize(m_n_elem * 4);              // 4 faces

    for (int e = 0; e < m_n_elem; e++) {
        radTg3dRelax* elem = m_interaction->g3dRelaxPtrVect[e];
        if (!elem) continue;

        radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(elem);
        if (!poly || poly->AmOfFaces != 4) continue;  // Skip non-tetrahedra

        // Store element center
        int cIdx = e * 3;
        m_tetra_centers[cIdx + 0] = poly->CentrPoint.x;
        m_tetra_centers[cIdx + 1] = poly->CentrPoint.y;
        m_tetra_centers[cIdx + 2] = poly->CentrPoint.z;

        // Store face data for each of the 4 triangular faces
        for (int f = 0; f < 4; f++) {
            radTHandlePgnAndTrans hpt = poly->VectHandlePgnAndTrans[f];
            radTPolygon* pgn = hpt.PgnHndl.rep;
            radTrans* tr = hpt.TransHndl.rep;

            // Get 3 vertices of this triangular face
            const radTVect2dVect& verts2d = pgn->EdgePointsVector;
            if (verts2d.size() < 3) continue;

            int fvIdx = (e * 4 + f) * 3 * 3;  // Starting index for this face's vertices
            TVector3d V[3];
            for (int v = 0; v < 3; v++) {
                V[v] = tr->TrPoint(TVector3d(verts2d[v].x, verts2d[v].y, pgn->CoordZ));
                m_tetra_face_vertices[fvIdx + v * 3 + 0] = V[v].x;
                m_tetra_face_vertices[fvIdx + v * 3 + 1] = V[v].y;
                m_tetra_face_vertices[fvIdx + v * 3 + 2] = V[v].z;
            }

            // Compute face normal (outward pointing)
            TVector3d e1 = {V[1].x - V[0].x, V[1].y - V[0].y, V[1].z - V[0].z};
            TVector3d e2 = {V[2].x - V[0].x, V[2].y - V[0].y, V[2].z - V[0].z};
            TVector3d n = {e1.y*e2.z - e1.z*e2.y, e1.z*e2.x - e1.x*e2.z, e1.x*e2.y - e1.y*e2.x};
            double nLen = sqrt(n.x*n.x + n.y*n.y + n.z*n.z);

            // Face area = 0.5 * |cross product|
            m_tetra_face_areas[e * 4 + f] = 0.5 * nLen;

            // Normalize and check orientation (outward from centroid)
            if (nLen > 1e-20) {
                n.x /= nLen; n.y /= nLen; n.z /= nLen;

                // Face center
                TVector3d fc = {(V[0].x + V[1].x + V[2].x) / 3.0,
                                (V[0].y + V[1].y + V[2].y) / 3.0,
                                (V[0].z + V[1].z + V[2].z) / 3.0};
                // Vector from centroid to face center
                TVector3d toFace = {fc.x - poly->CentrPoint.x,
                                    fc.y - poly->CentrPoint.y,
                                    fc.z - poly->CentrPoint.z};
                // If normal points inward, flip it
                if (n.x*toFace.x + n.y*toFace.y + n.z*toFace.z < 0) {
                    n.x = -n.x; n.y = -n.y; n.z = -n.z;
                }
            }

            // Store normalized outward normal
            int fnIdx = (e * 4 + f) * 3;
            m_tetra_face_normals[fnIdx + 0] = n.x;
            m_tetra_face_normals[fnIdx + 1] = n.y;
            m_tetra_face_normals[fnIdx + 2] = n.z;
        }
    }

    m_geometry_3dof_ready = true;
}

//=========================================================================
// PrecomputeFlatInteractMatrix: Flatten InteractMatrix for 3DOF tetrahedra
// Converts 2D pointer array TMatrix3df** to contiguous double array
// This eliminates pointer chasing during matrix element access
//=========================================================================

void RadHACApKManager::PrecomputeFlatInteractMatrix() {
    if (m_flat_N_ready || m_n_elem == 0 || !m_interaction) return;
    if (!m_interaction->InteractMatrix) {
        return;  // InteractMatrix not computed
    }

    // Allocate flat storage: n_elem * n_elem * 9 doubles (3x3 block per element pair)
    int64_t total_size = (int64_t)m_n_elem * m_n_elem * 9;
    m_flat_N_data.resize(total_size);

    // Copy from InteractMatrix[i][j] to flat array with sign flip (-N)
    // The system matrix is A = -N + diag(1/chi), so we store -N
    //
    // MATRIX LAYOUT FIX (2025-12-24):
    // InteractMatrix[i][j] stores TMatrix3df where:
    //   Str0 = dH/dMx = (dHx/dMx, dHy/dMx, dHz/dMx) <- COLUMN vector (response to Mx)
    //   Str1 = dH/dMy = (dHx/dMy, dHy/dMy, dHz/dMy) <- COLUMN vector (response to My)
    //   Str2 = dH/dMz = (dHx/dMz, dHy/dMz, dHz/dMz) <- COLUMN vector (response to Mz)
    //
    // For row k of the 3x3 block, we need N[i][j] element (k, l) = dH_k/dM_l
    // This requires transposed access: row k gets Str0[k], Str1[k], Str2[k]
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < m_n_elem; i++) {
        for (int j = 0; j < m_n_elem; j++) {
            const TMatrix3df& M = m_interaction->InteractMatrix[i][j];
            int64_t base_idx = ((int64_t)i * m_n_elem + j) * 9;

            // Row 0 (Hx response): -dHx/dMx, -dHx/dMy, -dHx/dMz
            m_flat_N_data[base_idx + 0] = -static_cast<double>(M.Str0.x);  // -dHx/dMx
            m_flat_N_data[base_idx + 1] = -static_cast<double>(M.Str1.x);  // -dHx/dMy
            m_flat_N_data[base_idx + 2] = -static_cast<double>(M.Str2.x);  // -dHx/dMz
            // Row 1 (Hy response): -dHy/dMx, -dHy/dMy, -dHy/dMz
            m_flat_N_data[base_idx + 3] = -static_cast<double>(M.Str0.y);  // -dHy/dMx
            m_flat_N_data[base_idx + 4] = -static_cast<double>(M.Str1.y);  // -dHy/dMy
            m_flat_N_data[base_idx + 5] = -static_cast<double>(M.Str2.y);  // -dHy/dMz
            // Row 2 (Hz response): -dHz/dMx, -dHz/dMy, -dHz/dMz
            m_flat_N_data[base_idx + 6] = -static_cast<double>(M.Str0.z);  // -dHz/dMx
            m_flat_N_data[base_idx + 7] = -static_cast<double>(M.Str1.z);  // -dHz/dMy
            m_flat_N_data[base_idx + 8] = -static_cast<double>(M.Str2.z);  // -dHz/dMz
        }
    }

    m_flat_N_ready = true;
}

//=========================================================================
// FieldFromChargedTriangleFast: Triangle field using edge-based log/atan formula
// Pure array-based computation (no object access)
// Returns field WITHOUT 4pi divisor (matches radTPolyhedron::FieldFromChargedTriangle)
//=========================================================================

void RadHACApKManager::FieldFromChargedTriangleFast(const double* obs, const double* v0, const double* v1, const double* v2, double sigma, double* H_out) const {
    // Analytic field from uniformly charged triangle using edge-based log/atan formula
    // Returns field WITHOUT 4pi divisor (4pi is applied once in matrix assembly)
    // Matches radTPolyhedron::FieldFromChargedTriangle implementation

    const double EPS = 1.0e-20;

    // Build local coordinate system
    double e1[3] = {v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]};
    double e2[3] = {v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]};

    // Face normal = e1 x e2 (basis_c)
    double basis_c[3];
    basis_c[0] = e1[1]*e2[2] - e1[2]*e2[1];
    basis_c[1] = e1[2]*e2[0] - e1[0]*e2[2];
    basis_c[2] = e1[0]*e2[1] - e1[1]*e2[0];

    double cLen = sqrt(basis_c[0]*basis_c[0] + basis_c[1]*basis_c[1] + basis_c[2]*basis_c[2]);
    if (cLen < EPS) { H_out[0] = H_out[1] = H_out[2] = 0.0; return; }
    basis_c[0] /= cLen; basis_c[1] /= cLen; basis_c[2] /= cLen;

    // basis_a = e1 normalized
    double basis_a[3] = {e1[0], e1[1], e1[2]};
    double aLen = sqrt(basis_a[0]*basis_a[0] + basis_a[1]*basis_a[1] + basis_a[2]*basis_a[2]);
    if (aLen < EPS) { H_out[0] = H_out[1] = H_out[2] = 0.0; return; }
    basis_a[0] /= aLen; basis_a[1] /= aLen; basis_a[2] /= aLen;

    // basis_b = basis_c x basis_a
    double basis_b[3];
    basis_b[0] = basis_c[1]*basis_a[2] - basis_c[2]*basis_a[1];
    basis_b[1] = basis_c[2]*basis_a[0] - basis_c[0]*basis_a[2];
    basis_b[2] = basis_c[0]*basis_a[1] - basis_c[1]*basis_a[0];
    double bLen = sqrt(basis_b[0]*basis_b[0] + basis_b[1]*basis_b[1] + basis_b[2]*basis_b[2]);
    if (bLen < EPS) { H_out[0] = H_out[1] = H_out[2] = 0.0; return; }
    basis_b[0] /= bLen; basis_b[1] /= bLen; basis_b[2] /= bLen;

    // Convert vertices to local 2D coordinates (v0 = origin)
    double xy0_x = 0.0, xy0_y = 0.0;
    double xy1_x = e1[0]*basis_a[0] + e1[1]*basis_a[1] + e1[2]*basis_a[2];
    double xy1_y = e1[0]*basis_b[0] + e1[1]*basis_b[1] + e1[2]*basis_b[2];
    double xy2_x = e2[0]*basis_a[0] + e2[1]*basis_a[1] + e2[2]*basis_a[2];
    double xy2_y = e2[0]*basis_b[0] + e2[1]*basis_b[1] + e2[2]*basis_b[2];

    double XY[3][2] = {{xy0_x, xy0_y}, {xy1_x, xy1_y}, {xy2_x, xy2_y}};
    double DS[3], AM[3], SM[3], XD[3], YD[3];
    double EPSG = 0.0;

    for (int j = 0; j < 3; j++) {
        int l = (j + 1) % 3;
        double dx = XY[l][0] - XY[j][0];
        double dy = XY[l][1] - XY[j][1];
        if (fabs(dx) < EPS) dx = (dx >= 0) ? EPS : -EPS;

        DS[j] = sqrt(dx*dx + dy*dy);
        AM[j] = dy / dx;
        SM[j] = sqrt(AM[j]*AM[j] + 1.0);
        XD[j] = -dx / DS[j];
        YD[j] =  dy / DS[j];

        if (DS[j] > EPSG) EPSG = DS[j];
    }
    EPSG *= 1.0e-12;

    // Transform observation point to local coordinates
    double d[3] = {obs[0]-v0[0], obs[1]-v0[1], obs[2]-v0[2]};
    double EE1 = d[0]*basis_a[0] + d[1]*basis_a[1] + d[2]*basis_a[2];
    double EE2 = d[0]*basis_b[0] + d[1]*basis_b[1] + d[2]*basis_b[2];
    double EE3 = d[0]*basis_c[0] + d[1]*basis_c[1] + d[2]*basis_c[2];

    double X[3], Y[3], H[3], E[3], R[3];
    for (int j = 0; j < 3; j++) {
        X[j] = EE1 - XY[j][0];
        Y[j] = EE2 - XY[j][1];
        H[j] = Y[j] * X[j];
        E[j] = EE3*EE3 + X[j]*X[j];
        R[j] = sqrt(X[j]*X[j] + Y[j]*Y[j] + EE3*EE3);
    }

    double Z = EE3;

    // Edge contributions
    double RM[3], RP[3], RR[3], AL[3];
    for (int j = 0; j < 3; j++) {
        int jp1 = (j + 1) % 3;
        RM[j] = R[j] + R[jp1] - DS[j];
        RP[j] = R[j] + R[jp1] + DS[j];
        RR[j] = (RM[j] / RP[j] > EPS) ? (RM[j] / RP[j]) : EPS;
        AL[j] = log(RR[j]);
    }

    // Field components in local frame WITHOUT 4pi divisor
    double HH1 = sigma * (-YD[0]*AL[0] - YD[1]*AL[1] - YD[2]*AL[2]);
    double HH2 = sigma * (-XD[0]*AL[0] - XD[1]*AL[1] - XD[2]*AL[2]);
    double HH3 = 0.0;

    // Normal component (atan terms) - only if not on surface
    if (fabs(Z) > EPSG) {
        double ZR[3];
        for (int j = 0; j < 3; j++) {
            ZR[j] = Z * R[j];
        }

        double AT[3], BT[3];
        for (int j = 0; j < 3; j++) {
            int jp1 = (j + 1) % 3;
            AT[j] = (AM[j]*E[j] - H[j]) / ZR[j];
            BT[j] = (AM[j]*E[jp1] - H[jp1]) / ZR[jp1];
        }

        HH3 = sigma * (-atan(AT[0]) - atan(AT[1]) - atan(AT[2])
                       +atan(BT[0]) + atan(BT[1]) + atan(BT[2]));
    }

    // Transform back to global coordinates
    H_out[0] = HH1*basis_a[0] + HH2*basis_b[0] + HH3*basis_c[0];
    H_out[1] = HH1*basis_a[1] + HH2*basis_b[1] + HH3*basis_c[1];
    H_out[2] = HH1*basis_a[2] + HH2*basis_b[2] + HH3*basis_c[2];
}

//=========================================================================
// FieldFromTrianglePrecomputed: Ultra-fast field using pre-computed data
// Uses pre-computed basis vectors and edge parameters (no sqrt/div)
// Returns field WITHOUT 4pi divisor
//=========================================================================

void RadHACApKManager::FieldFromTrianglePrecomputed(int tri_idx, const double* obs, double sigma, double* H_out) const {
    const double EPS = 1.0e-20;

    // Get pre-computed data
    const double* data = &m_tri_data[tri_idx * TRI_DATA_SIZE];

    const double* basis_a = data;       // [0..2]
    const double* basis_b = data + 3;   // [3..5]
    const double* basis_c = data + 6;   // [6..8]
    const double* origin = data + 9;    // [9..11]
    const double* XY = data + 12;       // [12..17] = {x0,y0, x1,y1, x2,y2}
    const double* DS = data + 18;       // [18..20]
    const double* AM = data + 21;       // [21..23]
    const double* XD = data + 24;       // [24..26]
    const double* YD = data + 27;       // [27..29]
    double EPSG = data[30];
    double sign = data[31];

    // Apply sign to sigma
    sigma *= sign;

    // Transform observation point to local coordinates
    double d[3] = {obs[0] - origin[0], obs[1] - origin[1], obs[2] - origin[2]};
    double EE1 = d[0]*basis_a[0] + d[1]*basis_a[1] + d[2]*basis_a[2];
    double EE2 = d[0]*basis_b[0] + d[1]*basis_b[1] + d[2]*basis_b[2];
    double EE3 = d[0]*basis_c[0] + d[1]*basis_c[1] + d[2]*basis_c[2];

    // Vertex-relative coordinates
    double X[3], Y[3], H[3], E[3], R[3];
    for (int j = 0; j < 3; j++) {
        X[j] = EE1 - XY[j*2];
        Y[j] = EE2 - XY[j*2+1];
        H[j] = Y[j] * X[j];
        E[j] = EE3*EE3 + X[j]*X[j];
        R[j] = sqrt(X[j]*X[j] + Y[j]*Y[j] + EE3*EE3);
    }

    double Z = EE3;

    // Edge contributions (log terms)
    double AL[3];
    for (int j = 0; j < 3; j++) {
        int jp1 = (j + 1) % 3;
        double RM = R[j] + R[jp1] - DS[j];
        double RP = R[j] + R[jp1] + DS[j];
        double RR = (RM / RP > EPS) ? (RM / RP) : EPS;
        AL[j] = log(RR);
    }

    // Tangential field components
    double HH1 = sigma * (-YD[0]*AL[0] - YD[1]*AL[1] - YD[2]*AL[2]);
    double HH2 = sigma * (-XD[0]*AL[0] - XD[1]*AL[1] - XD[2]*AL[2]);
    double HH3 = 0.0;

    // Normal component (atan terms) - only if not on surface
    if (fabs(Z) > EPSG) {
        double ZR[3];
        for (int j = 0; j < 3; j++) {
            ZR[j] = Z * R[j];
        }

        double AT[3], BT[3];
        for (int j = 0; j < 3; j++) {
            int jp1 = (j + 1) % 3;
            AT[j] = (AM[j]*E[j] - H[j]) / ZR[j];
            BT[j] = (AM[j]*E[jp1] - H[jp1]) / ZR[jp1];
        }

        HH3 = sigma * (-atan(AT[0]) - atan(AT[1]) - atan(AT[2])
                       +atan(BT[0]) + atan(BT[1]) + atan(BT[2]));
    }

    // Transform back to global coordinates
    H_out[0] = HH1*basis_a[0] + HH2*basis_b[0] + HH3*basis_c[0];
    H_out[1] = HH1*basis_a[1] + HH2*basis_b[1] + HH3*basis_c[1];
    H_out[2] = HH1*basis_a[2] + HH2*basis_b[2] + HH3*basis_c[2];
}

//=========================================================================
// FieldFromQuadFaceFast: Quad face field using pre-computed vertices
// Split quad into 2 triangles, sum contributions
//=========================================================================

void RadHACApKManager::FieldFromQuadFaceFast(int elem, int face, const double* obs, double sigma, double* H_out) const {
    H_out[0] = H_out[1] = H_out[2] = 0.0;

    // Use pre-computed triangle data if available (much faster)
    if (m_tri_precomputed) {
        double H_tri[3];

        // Triangle 1: indices [0,1,2] of quad
        int tri_idx1 = elem * TRIS_PER_ELEM + face * 2 + 0;
        FieldFromTrianglePrecomputed(tri_idx1, obs, sigma, H_tri);
        H_out[0] += H_tri[0]; H_out[1] += H_tri[1]; H_out[2] += H_tri[2];

        // Triangle 2: indices [0,2,3] of quad
        int tri_idx2 = elem * TRIS_PER_ELEM + face * 2 + 1;
        FieldFromTrianglePrecomputed(tri_idx2, obs, sigma, H_tri);
        H_out[0] += H_tri[0]; H_out[1] += H_tri[1]; H_out[2] += H_tri[2];

        return;
    }

    // Fallback: compute on-the-fly (for cases where precomputation wasn't done)
    int fvIdx = (elem * 6 + face) * 4 * 3;
    const double* V0 = &m_face_vertices[fvIdx + 0];
    const double* V1 = &m_face_vertices[fvIdx + 3];
    const double* V2 = &m_face_vertices[fvIdx + 6];
    const double* V3 = &m_face_vertices[fvIdx + 9];

    const double* center = &m_elem_centers[elem * 3];

    double H_tri[3];

    // Triangle 1: V0, V1, V2
    double tc1[3] = {(V0[0]+V1[0]+V2[0])/3.0, (V0[1]+V1[1]+V2[1])/3.0, (V0[2]+V1[2]+V2[2])/3.0};
    double e1_1[3] = {V1[0]-V0[0], V1[1]-V0[1], V1[2]-V0[2]};
    double e2_1[3] = {V2[0]-V0[0], V2[1]-V0[1], V2[2]-V0[2]};
    double n1[3] = {e1_1[1]*e2_1[2]-e1_1[2]*e2_1[1], e1_1[2]*e2_1[0]-e1_1[0]*e2_1[2], e1_1[0]*e2_1[1]-e1_1[1]*e2_1[0]};
    double to_c1[3] = {tc1[0]-center[0], tc1[1]-center[1], tc1[2]-center[2]};
    double dot1 = n1[0]*to_c1[0] + n1[1]*to_c1[1] + n1[2]*to_c1[2];
    double sign1 = (dot1 >= 0.0) ? 1.0 : -1.0;

    FieldFromChargedTriangleFast(obs, V0, V1, V2, sigma * sign1, H_tri);
    H_out[0] += H_tri[0]; H_out[1] += H_tri[1]; H_out[2] += H_tri[2];

    // Triangle 2: V0, V2, V3
    double tc2[3] = {(V0[0]+V2[0]+V3[0])/3.0, (V0[1]+V2[1]+V3[1])/3.0, (V0[2]+V2[2]+V3[2])/3.0};
    double e1_2[3] = {V2[0]-V0[0], V2[1]-V0[1], V2[2]-V0[2]};
    double e2_2[3] = {V3[0]-V0[0], V3[1]-V0[1], V3[2]-V0[2]};
    double n2[3] = {e1_2[1]*e2_2[2]-e1_2[2]*e2_2[1], e1_2[2]*e2_2[0]-e1_2[0]*e2_2[2], e1_2[0]*e2_2[1]-e1_2[1]*e2_2[0]};
    double to_c2[3] = {tc2[0]-center[0], tc2[1]-center[1], tc2[2]-center[2]};
    double dot2 = n2[0]*to_c2[0] + n2[1]*to_c2[1] + n2[2]*to_c2[2];
    double sign2 = (dot2 >= 0.0) ? 1.0 : -1.0;

    FieldFromChargedTriangleFast(obs, V0, V2, V3, sigma * sign2, H_tri);
    H_out[0] += H_tri[0]; H_out[1] += H_tri[1]; H_out[2] += H_tri[2];
}

//=========================================================================
// Compute6x6BlockFast: ELF-style fast 6x6 block computation
// Uses only pre-computed geometry arrays (no object access)
//=========================================================================

void RadHACApKManager::Compute6x6BlockFast(int elem_i, int elem_j, double* K_mat) const {
    if (!m_geometry_ready) {
        std::memset(K_mat, 0, 36 * sizeof(double));
        return;
    }

    // Pre-computed row element data
    const double* row_center = &m_elem_centers[elem_i * 3];

    // Pre-compute evaluation points for all 6 faces of row element
    double eval_pts[6][3];
    for (int fi = 0; fi < 6; fi++) {
        const double* fc = &m_face_centers[(elem_i * 6 + fi) * 3];
        eval_pts[fi][0] = 0.5 * (fc[0] + row_center[0]);
        eval_pts[fi][1] = 0.5 * (fc[1] + row_center[1]);
        eval_pts[fi][2] = 0.5 * (fc[2] + row_center[2]);
    }

    // Pre-load row normals
    double row_normals[6][3];
    for (int fi = 0; fi < 6; fi++) {
        const double* fn = &m_face_normals[(elem_i * 6 + fi) * 3];
        row_normals[fi][0] = fn[0];
        row_normals[fi][1] = fn[1];
        row_normals[fi][2] = fn[2];
    }

    // Get col element center for point charge
    const double* col_center = &m_elem_centers[elem_j * 3];

    // Compute all 36 elements
    for (int fj = 0; fj < 6; fj++) {
        double col_area = m_face_areas[elem_j * 6 + fj];
        double unit_point_charge = -1.0 * col_area;

        for (int fi = 0; fi < 6; fi++) {
            // Field from quad face with unit sigma
            double H_face[3];
            FieldFromQuadFaceFast(elem_j, fj, eval_pts[fi], 1.0, H_face);

            // Field from point charge at element center
            // Note: Field is computed WITHOUT 4pi divisor to match Compute6x6Block
            // The 4pi divisor is applied once at the end (K_mat[...] = -K_ij * INV_4PI)
            double r[3] = {eval_pts[fi][0] - col_center[0],
                           eval_pts[fi][1] - col_center[1],
                           eval_pts[fi][2] - col_center[2]};
            double dist = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
            double H_point[3] = {0.0, 0.0, 0.0};
            if (dist > 1e-15) {
                double scale = unit_point_charge / (dist * dist * dist);  // NO 4pi here
                H_point[0] = r[0] * scale;
                H_point[1] = r[1] * scale;
                H_point[2] = r[2] * scale;
            }

            double H_total[3] = {H_face[0] + H_point[0],
                                  H_face[1] + H_point[1],
                                  H_face[2] + H_point[2]};

            // K_ij = H_total dot normal_i
            double K_ij = H_total[0] * row_normals[fi][0] +
                          H_total[1] * row_normals[fi][1] +
                          H_total[2] * row_normals[fi][2];

            // Store -K_ij / (4*pi) in row-major order
            K_mat[fi * 6 + fj] = -K_ij * RadConst::INV_FOUR_PI;
        }
    }
}

bool RadHACApKManager::BuildHMatrix(const RadHACApKParams& params) {
    FreeResources();

    if (!m_interaction) {
        std::cerr << "[HACApK] Error: No interaction object" << std::endl;
        return false;
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    ExtractElementCoordinates();

    if (m_n_elem == 0) {
        std::cerr << "[HACApK] Error: No elements" << std::endl;
        return false;
    }

    // Set global callback state
    RadHACApKCallback::SetInteraction(m_interaction, m_n_elem, m_nffc);
    RadHACApKCallback::SetCurrentManager(this);

    // Verify DOF was configured correctly (3DOF, 6DOF, or mixed=0)
    if (m_ndof == 0 || (m_nffc != 0 && m_nffc != 3 && m_nffc != 6)) {
        std::cerr << "[HACApK] Error: Invalid DOF configuration (nffc=" << m_nffc << ")" << std::endl;
        return false;
    }

    // ELF-style pre-computation for 6DOF hexahedra
    if (m_is_6dof) {
        PrecomputeGeometry();
    }
    // ELF-style pre-computation for 3DOF tetrahedra (2025-12-26)
    // PrecomputeGeometry3DOF extracts face vertices/normals for direct field computation
    // This allows on-demand matrix element computation without O(N^2) SetupInteractMatrix()
    if (!m_is_6dof && !m_is_mixed_dof) {
        PrecomputeGeometry3DOF();
    }
    // FIX (2025-12-26): For 3DOF tetrahedra, use PrecomputeFlatInteractMatrix()
    // if InteractMatrix was already computed (fallback path)
    // This provides O(1) matrix element access during H-matrix construction
    if (!m_is_6dof && !m_geometry_3dof_ready && m_interaction->InteractMatrix != nullptr) {
        PrecomputeFlatInteractMatrix();
    }

    // Initialize inverse susceptibility with ELF-style initial chi from BH curve point 2
    // This matches ELF's initialize_chi_from_bh() which uses:
    //   chi = B2/(mu0*H2) - 1 (from 2nd point of BH curve)
    m_inv_chi.resize(m_ndof);

    for (int i = 0; i < m_n_elem; i++) {
        radTg3dRelax* g3dRelaxPtr = m_interaction->g3dRelaxPtrVect[i];
        radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

        // Use ELF-style initial chi from BH curve point 2 (matches initialize_chi_from_bh)
        double chi = 1.0;
        radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);
        if (NonlinMater != nullptr) {
            chi = NonlinMater->GetInitialChi_ELF_Style();
            if (chi <= 0) chi = 1.0;
        } else {
            // Fallback for linear materials: use DefineInstantKsiTensor with H=0
            TMatrix3d KsiTensor;
            TVector3d MrVect;
            TVector3d H_zero(0., 0., 0.);
            MaterPtr->DefineInstantKsiTensor(H_zero, KsiTensor, MrVect);
            chi = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
        }
        if (chi < 1.0e-6) chi = 1.0e-6;
        double inv_chi_val = 1.0 / chi;

        // All DOF per element use the same 1/chi
        int offset = m_dof_offset[i];
        int elem_dof = m_dof_offset[i + 1] - offset;
        for (int k = 0; k < elem_dof; k++) {
            m_inv_chi[offset + k] = inv_chi_val;
        }
    }

    RadHACApKCallback::SetInvChi(m_inv_chi);

    // Allocate opaque structures
    m_leafmtxp = HACApK_alloc_leafmtxp();
    m_control = HACApK_alloc_lcontrol();

    if (!m_leafmtxp || !m_control) {
        std::cerr << "[HACApK] Error: Failed to allocate structures" << std::endl;
        return false;
    }

    // Build H-matrix using the C wrapper
    // m_nffc: 3 for tetrahedra, 6 for hexahedra, 0 for mixed
    int ndim = 3;  // Spatial dimension

    auto t_hmatrix_start = std::chrono::high_resolution_clock::now();
    int result;

    if (m_is_mixed_dof) {
        // Variable DOF mode: use new varDOF wrapper
        result = HACApK_build_hmatrix_varDOF_wrapper(
            m_leafmtxp,
            m_control,
            m_coordinates.data(),
            m_n_elem,
            m_dof_offset.data(),
            m_ndof,
            ndim,
            params.aca_eps,
            params.leaf_size,
            params.eta,
            params.print_level
        );
    } else {
        // Uniform DOF mode (original)
        result = HACApK_build_hmatrix_wrapper(
            m_leafmtxp,
            m_control,
            m_coordinates.data(),
            m_n_elem,
            m_nffc,  // 3 for tetra, 6 for hexa
            ndim,
            params.aca_eps,
            params.leaf_size,
            params.eta,
            params.print_level
        );
    }
    auto t_hmatrix_end = std::chrono::high_resolution_clock::now();
    double t_hmatrix = std::chrono::duration<double>(t_hmatrix_end - t_hmatrix_start).count();
    if (params.print_level > 0) {
        std::cout << "[HACApK] HACApK_build_hmatrix_wrapper: " << t_hmatrix << " s" << std::endl;
    }

    if (result != 0) {
        std::cerr << "[HACApK] Error: H-matrix build failed with code " << result << std::endl;
        FreeResources();
        return false;
    }

    // Store permutation from control structure
    int* lod = HACApK_lcontrol_get_lod(m_control);
    if (lod) {
        m_permutation.resize(m_ndof);
        for (int i = 0; i < m_ndof; i++) {
            m_permutation[i] = lod[i + 1];
        }
    }

    // Get statistics from leafmtxp
    m_stats.n_dof = HACApK_leafmtxp_get_nd(m_leafmtxp);
    m_stats.n_leaves = HACApK_leafmtxp_get_nlf(m_leafmtxp);
    m_stats.n_lowrank = HACApK_leafmtxp_get_nlfkt(m_leafmtxp);
    m_stats.n_dense = m_stats.n_leaves - m_stats.n_lowrank;
    m_stats.max_rank = HACApK_leafmtxp_get_ktmax(m_leafmtxp);

    // Calculate memory usage and compression ratio (ELF-compatible)
    // Use accurate per-leaf calculation instead of rough estimation
    int64_t hmat_bytes = 0;
    int64_t dense_bytes = 0;
    HACApK_get_memory_stats(m_leafmtxp, &hmat_bytes, &dense_bytes);

    m_stats.memory_mb = (double)hmat_bytes / (1024.0 * 1024.0);
    m_stats.dense_memory_mb = (double)dense_bytes / (1024.0 * 1024.0);

    // Compression ratio = H-matrix memory / Dense matrix memory
    // Note: dense_bytes is the sum of block sizes, not full N^2 matrix
    // This matches ELF's definition where compression < 1 means memory saved
    m_stats.compression = (dense_bytes > 0) ?
        (double)hmat_bytes / (double)dense_bytes : 1.0;

    auto end_time = std::chrono::high_resolution_clock::now();
    m_stats.build_time = std::chrono::duration<double>(end_time - start_time).count();

    if (params.print_level > 0) {
        std::cout << "[HACApK] H-matrix built: "
                  << "DOF=" << m_stats.n_dof
                  << ", leaves=" << m_stats.n_leaves
                  << ", lowrank=" << m_stats.n_lowrank
                  << ", dense=" << m_stats.n_dense
                  << ", maxrank=" << m_stats.max_rank
                  << ", time=" << m_stats.build_time << "s"
                  << std::endl;
    }

    // Cache diagonal elements N_ii for Jacobi preconditioner
    // This is computed once and reused in every BiCGSTAB iteration
    m_diag_N.resize(m_ndof);
    #pragma omp parallel for
    for (int i = 0; i < m_ndof; i++) {
        m_diag_N[i] = GetInteractionMatrixElement(i, i);
    }
    m_diag_cached = true;

    m_valid = true;
    return true;
}

void RadHACApKManager::MatVec(const std::vector<double>& x, std::vector<double>& y) {
    if (!m_valid || !m_leafmtxp || !m_control) {
        std::fill(y.begin(), y.end(), 0.0);
        return;
    }

    int nd = HACApK_leafmtxp_get_nd(m_leafmtxp);
    HACApK_matvec_wrapper(m_leafmtxp, m_control, x.data(), y.data(), nd);
}

void RadHACApKManager::UpdateDiagonal(const std::vector<double>& inv_chi) {
    if (!m_valid || !m_leafmtxp || !m_control) return;

    // Update stored inverse susceptibility
    m_inv_chi = inv_chi;
    RadHACApKCallback::SetInvChi(inv_chi);

    // DEBUG: Force slow method to verify correctness
    // TODO: Re-enable fast method after debugging
    // OPTIMIZATION: Use fast diagonal update
    // Only updates true diagonal entries (i==j) using pre-computed N_ii values
    // This is O(ndof) instead of O(block_size^2 * n_diag_blocks)
    // if (m_diag_cached && m_diag_N.size() == (size_t)m_ndof) {
    //     HACApK_update_diagonal_fast_wrapper(m_leafmtxp, m_control,
    //                                          m_diag_N.data(), inv_chi.data(), m_ndof);
    // } else {
        // Fallback to slow method (recomputes all entries in diagonal blocks)
        HACApK_update_diagonal_wrapper(m_leafmtxp, m_control, cHACApK_entry_ij);
    // }
}

//=========================================================================
// GetInteractionMatrixElement: Access N matrix element
//
// This function returns the N(dof_i, dof_j) element of the interaction matrix.
// On-demand computation based on element DOF type:
// - 3DOF tetrahedra: Use B_comp() with PreRelax mode
// - 6DOF hexahedra: Use Yano-Sugahara MSC method (face-to-face interaction)
//
// On-demand computation is essential for H-matrix because:
// - HACApK uses ACA+ which only needs a subset of matrix elements
// - Pre-computing the full dense matrix would defeat the O(N log N) purpose
//=========================================================================

//=========================================================================
// Compute6x6Block: Calculate full 6x6 interaction block (ELF-style)
// K(face_i, face_j) = normal_i dot H_field(eval_pt_i, src_face_j)
//=========================================================================

void RadHACApKManager::Compute6x6Block(int elem_i, int elem_j, double* K_mat) const {
    // Bounds check for element indices
    if (elem_i < 0 || elem_i >= m_n_elem || elem_j < 0 || elem_j >= m_n_elem) {
        std::cerr << "[HACApK] Error: Invalid element indices in Compute6x6Block: "
                  << elem_i << ", " << elem_j << " (n_elem=" << m_n_elem << ")" << std::endl;
        std::memset(K_mat, 0, 36 * sizeof(double));
        return;
    }

    radTg3dRelax* elem_row = m_interaction->g3dRelaxPtrVect[elem_i];
    radTg3dRelax* elem_col = m_interaction->g3dRelaxPtrVect[elem_j];

    if (!elem_row || !elem_col) {
        std::cerr << "[HACApK] Error: Null element pointer in Compute6x6Block: "
                  << (elem_row ? "col" : "row") << " elem" << std::endl;
        std::memset(K_mat, 0, 36 * sizeof(double));
        return;
    }

    radTPolyhedron* poly_row = dynamic_cast<radTPolyhedron*>(elem_row);
    radTPolyhedron* poly_col = dynamic_cast<radTPolyhedron*>(elem_col);

    if (!poly_row || !poly_col) {
        std::memset(K_mat, 0, 36 * sizeof(double));
        return;
    }

    // Pre-compute evaluation points and normals for row element (ELF optimization)
    TVector3d eval_pts[6];
    for (int fi = 0; fi < 6; fi++) {
        eval_pts[fi].x = 0.5 * (poly_row->FaceCenter[fi].x + poly_row->CentrPoint.x);
        eval_pts[fi].y = 0.5 * (poly_row->FaceCenter[fi].y + poly_row->CentrPoint.y);
        eval_pts[fi].z = 0.5 * (poly_row->FaceCenter[fi].z + poly_row->CentrPoint.z);
    }

    // Compute all 36 elements (source face outer loop for cache locality)
    // IMPORTANT: K_mat indexing must match GetCached6x6Element access pattern
    // K_mat[fi * 6 + fj] stores K(face_i, face_j) = normal_i dot H(eval_pt_i, src_face_j)
    for (int fj = 0; fj < 6; fj++) {
        double unit_point_charge = -1.0 * poly_col->FaceArea[fj];

        for (int fi = 0; fi < 6; fi++) {
            // Field from unit sigma on face_j
            TVector3d H_face = poly_col->FieldFromQuadFace(eval_pts[fi], fj, 1.0);
            TVector3d H_point = poly_col->FieldFromPointCharge(eval_pts[fi], unit_point_charge);

            TVector3d H_total;
            H_total.x = H_face.x + H_point.x;
            H_total.y = H_face.y + H_point.y;
            H_total.z = H_face.z + H_point.z;

            // K_ij = H_total dot normal_i
            double K_ij = H_total.x * poly_row->FaceNormal[fi].x +
                          H_total.y * poly_row->FaceNormal[fi].y +
                          H_total.z * poly_row->FaceNormal[fi].z;

            // Store -K_ij / (4*pi) in row-major order: K_mat[row * 6 + col]
            K_mat[fi * 6 + fj] = -K_ij * RadConst::INV_FOUR_PI;
        }
    }
}

//=========================================================================
// GetCached6x6Element: ELF-style hash-based thread-local cache (NO locking!)
//
// ELF pattern from m_ppohBEM_user_func_unified.f90:
// - Single-entry cache checked first (most common case)
// - Hash-based cache (O(1) lookup, no LRU search)
// - NO global cache, NO critical sections
//
// This is THE key optimization: hash lookup + no locking = fast scaling
//=========================================================================

// Hash-based cache size (must be power of 2 for fast modulo)
static constexpr int TL_HASH_SIZE = 1024;
static constexpr int TL_HASH_MASK = TL_HASH_SIZE - 1;

double RadHACApKManager::GetCached6x6Element(int elem_i, int elem_j, int face_i, int face_j) const {
    // Thread-local single-entry cache (most common case: same block accessed multiple times)
    static thread_local int tl_single_elem_i = -1;
    static thread_local int tl_single_elem_j = -1;
    static thread_local double tl_single_K_mat[36];

    // Thread-local hash-based cache (O(1) lookup, no locking!)
    static thread_local int tl_cache_elem_i[TL_HASH_SIZE];
    static thread_local int tl_cache_elem_j[TL_HASH_SIZE];
    static thread_local double tl_cache_K_mat[TL_HASH_SIZE][36];
    static thread_local bool tl_initialized = false;

    // Initialize thread-local cache on first access
    if (!tl_initialized) {
        for (int i = 0; i < TL_HASH_SIZE; i++) {
            tl_cache_elem_i[i] = -1;
            tl_cache_elem_j[i] = -1;
        }
        tl_initialized = true;
    }

    // Check single-entry cache first (fastest path)
    if (tl_single_elem_i == elem_i && tl_single_elem_j == elem_j) {
        return tl_single_K_mat[face_i * 6 + face_j];
    }

    // Compute hash index (ELF-style hash function)
    int hash_idx = ((elem_i * 73856093) ^ (elem_j * 19349663)) & TL_HASH_MASK;

    // Check hash cache (O(1) lookup!)
    if (tl_cache_elem_i[hash_idx] == elem_i && tl_cache_elem_j[hash_idx] == elem_j) {
        // Cache hit - copy to single-entry cache for repeated access
        std::memcpy(tl_single_K_mat, tl_cache_K_mat[hash_idx], 36 * sizeof(double));
        tl_single_elem_i = elem_i;
        tl_single_elem_j = elem_j;
        return tl_single_K_mat[face_i * 6 + face_j];
    }

    // Cache miss - compute the block using pre-computed geometry
    if (m_geometry_ready) {
        Compute6x6BlockFast(elem_i, elem_j, tl_single_K_mat);
    } else {
        Compute6x6Block(elem_i, elem_j, tl_single_K_mat);
    }
    tl_single_elem_i = elem_i;
    tl_single_elem_j = elem_j;

    // Insert into hash cache (overwrites any existing entry at this slot)
    tl_cache_elem_i[hash_idx] = elem_i;
    tl_cache_elem_j[hash_idx] = elem_j;
    std::memcpy(tl_cache_K_mat[hash_idx], tl_single_K_mat, 36 * sizeof(double));

    return tl_single_K_mat[face_i * 6 + face_j];
}

//=========================================================================
// GetInteractionMatrixElement: Optimized with O(1) lookup and LRU cache
// Supports 3DOF tetrahedra, 6DOF hexahedra, and mixed meshes
//=========================================================================

double RadHACApKManager::GetInteractionMatrixElement(int dof_i, int dof_j) const {
    if (!m_interaction || dof_i < 0 || dof_i >= m_ndof || dof_j < 0 || dof_j >= m_ndof) {
        return 0.0;
    }

    // Safety check for lookup tables
    if (m_dof_to_elem.empty() || m_dof_to_local.empty()) {
        std::cerr << "[HACApK] Error: DOF lookup tables not initialized" << std::endl;
        return 0.0;
    }

    // O(1) DOF-to-element lookup (ELF-style optimization)
    int elem_i = m_dof_to_elem[dof_i];
    int local_i = m_dof_to_local[dof_i];
    int elem_j = m_dof_to_elem[dof_j];
    int local_j = m_dof_to_local[dof_j];

    // Safety check for element indices
    if (elem_i < 0 || elem_i >= m_n_elem || elem_j < 0 || elem_j >= m_n_elem) {
        return 0.0;
    }

    // Get element DOF counts
    int dof_elem_i = m_dof_offset[elem_i + 1] - m_dof_offset[elem_i];
    int dof_elem_j = m_dof_offset[elem_j + 1] - m_dof_offset[elem_j];

    // Dispatch based on DOF type of each element
    if (dof_elem_i == 6 && dof_elem_j == 6) {
        // 6DOF-6DOF: hex-hex interaction
        return GetCached6x6Element(elem_i, elem_j, local_i, local_j);
    } else if (dof_elem_i == 3 && dof_elem_j == 3) {
        // 3DOF-3DOF: tetra-tetra interaction
        return GetCached3x3Element(elem_i, elem_j, local_i, local_j);
    } else if (dof_elem_i == 3 && dof_elem_j == 6) {
        // 3DOF-6DOF: tetra-hex interaction (3x6 block)
        return GetMixed3x6Element(elem_i, elem_j, local_i, local_j);
    } else if (dof_elem_i == 6 && dof_elem_j == 3) {
        // 6DOF-3DOF: hex-tetra interaction (6x3 block)
        return GetMixed6x3Element(elem_i, elem_j, local_i, local_j);
    }

    return 0.0;
}

//=========================================================================
// GetCached3x3Element: On-demand 3x3 block computation with thread-local hash cache
// Similar to GetCached6x6Element for 6DOF hexahedra
// Uses O(1) hash lookup instead of O(n) LRU search
//=========================================================================

// Hash-based cache size for 3DOF (must be power of 2 for fast modulo)
static constexpr int TL_HASH_SIZE_3DOF = 1024;
static constexpr int TL_HASH_MASK_3DOF = TL_HASH_SIZE_3DOF - 1;

double RadHACApKManager::GetCached3x3Element(int elem_i, int elem_j, int comp_i, int comp_j) const {
    // If flat storage is ready (pre-computed), use O(1) direct access
    if (m_flat_N_ready) {
        int64_t base_idx = ((int64_t)elem_i * m_n_elem + elem_j) * 9;
        double val = m_flat_N_data[base_idx + comp_i * 3 + comp_j];
        return val;
    }

    // On-demand computation with thread-local hash cache
    // This path is used when SetupInteractMatrix() is NOT called (HACApK optimization)

    // Thread-local single-entry cache (most common case: same block accessed multiple times)
    static thread_local int tl_single_elem_i = -1;
    static thread_local int tl_single_elem_j = -1;
    static thread_local double tl_single_N_mat[9];

    // Thread-local hash cache (O(1) lookup, no locking needed!)
    static thread_local int tl_cache_elem_i[TL_HASH_SIZE_3DOF];
    static thread_local int tl_cache_elem_j[TL_HASH_SIZE_3DOF];
    static thread_local double tl_cache_N_mat[TL_HASH_SIZE_3DOF][9];
    static thread_local bool tl_initialized = false;

    // Initialize thread-local cache on first access
    if (!tl_initialized) {
        for (int i = 0; i < TL_HASH_SIZE_3DOF; i++) {
            tl_cache_elem_i[i] = -1;
            tl_cache_elem_j[i] = -1;
        }
        tl_initialized = true;
    }

    // Check single-entry cache first (fastest path)
    if (tl_single_elem_i == elem_i && tl_single_elem_j == elem_j) {
        return tl_single_N_mat[comp_i * 3 + comp_j];
    }

    // Compute hash index (same hash function as 6DOF cache)
    int hash_idx = ((elem_i * 73856093) ^ (elem_j * 19349663)) & TL_HASH_MASK_3DOF;

    // Check hash cache (O(1) lookup!)
    if (tl_cache_elem_i[hash_idx] == elem_i && tl_cache_elem_j[hash_idx] == elem_j) {
        // Cache hit - copy to single-entry cache for repeated access
        std::memcpy(tl_single_N_mat, tl_cache_N_mat[hash_idx], 9 * sizeof(double));
        tl_single_elem_i = elem_i;
        tl_single_elem_j = elem_j;

        return tl_single_N_mat[comp_i * 3 + comp_j];
    }

    // Cache miss - compute the 3x3 block on-demand
    // Use fast method with pre-computed geometry if available
    double N_mat[9];
    if (m_geometry_3dof_ready) {
        Compute3x3BlockFast(elem_i, elem_j, N_mat);
    } else {
        Compute3x3Block_OnDemand(elem_i, elem_j, N_mat);
    }

    // Update single-entry cache
    std::memcpy(tl_single_N_mat, N_mat, 9 * sizeof(double));
    tl_single_elem_i = elem_i;
    tl_single_elem_j = elem_j;

    // Store in hash cache (overwrites any existing entry at this hash index)
    tl_cache_elem_i[hash_idx] = elem_i;
    tl_cache_elem_j[hash_idx] = elem_j;
    std::memcpy(tl_cache_N_mat[hash_idx], N_mat, 9 * sizeof(double));

    return N_mat[comp_i * 3 + comp_j];
}

//=========================================================================
// Compute3x3Block: Compute 3x3 interaction block for tetrahedra
// N_ij = interaction matrix element between magnetization components
// Uses existing radTInteraction::InteractMatrix
//=========================================================================

void RadHACApKManager::Compute3x3Block(int elem_i, int elem_j, double* N_mat) const {
    // InteractMatrix[elem_i][elem_j] returns TMatrix3df (3x3 float matrix)
    //
    // MATRIX LAYOUT FIX (2025-12-24):
    // InteractMatrix[i][j] stores TMatrix3df where:
    //   Str0 = dH/dMx = (dHx/dMx, dHy/dMx, dHz/dMx) <- COLUMN vector (response to Mx)
    //   Str1 = dH/dMy = (dHx/dMy, dHy/dMy, dHz/dMy) <- COLUMN vector (response to My)
    //   Str2 = dH/dMz = (dHx/dMz, dHy/dMz, dHz/dMz) <- COLUMN vector (response to Mz)
    //
    // For row k of the 3x3 block, we need N[i][j] element (k, l) = dH_k/dM_l
    // This requires transposed access: row k gets Str0[k], Str1[k], Str2[k]
    //
    // IMPORTANT: Radia's InteractMatrix stores positive N, but the system matrix is:
    //   A = -N + diag(1/chi)
    // The LU solver uses: BaseMatrix = -Nij
    // For HACApK, GetInteractionMatrixElement should return -N to match.

    if (!m_interaction || !m_interaction->InteractMatrix) {
        std::memset(N_mat, 0, 9 * sizeof(double));
        return;
    }

    // Get the 3x3 block for element pair (i, j)
    // InteractMatrix is indexed by element indices, not DOF indices
    const TMatrix3df& M = m_interaction->InteractMatrix[elem_i][elem_j];

    // Extract to row-major double array with sign flip (-N) and transpose
    // Row 0 (Hx response): -dHx/dMx, -dHx/dMy, -dHx/dMz
    N_mat[0] = -static_cast<double>(M.Str0.x);  // -dHx/dMx
    N_mat[1] = -static_cast<double>(M.Str1.x);  // -dHx/dMy
    N_mat[2] = -static_cast<double>(M.Str2.x);  // -dHx/dMz
    // Row 1 (Hy response): -dHy/dMx, -dHy/dMy, -dHy/dMz
    N_mat[3] = -static_cast<double>(M.Str0.y);  // -dHy/dMx
    N_mat[4] = -static_cast<double>(M.Str1.y);  // -dHy/dMy
    N_mat[5] = -static_cast<double>(M.Str2.y);  // -dHy/dMz
    // Row 2 (Hz response): -dHz/dMx, -dHz/dMy, -dHz/dMz
    N_mat[6] = -static_cast<double>(M.Str0.z);  // -dHz/dMx
    N_mat[7] = -static_cast<double>(M.Str1.z);  // -dHz/dMy
    N_mat[8] = -static_cast<double>(M.Str2.z);  // -dHz/dMz
}

//=========================================================================
// Compute3x3Block_OnDemand: On-demand 3x3 interaction block computation
// Computes interaction directly using B_comp() without pre-computed matrix
// This is used by HACApK to avoid O(N^2) matrix pre-computation
//=========================================================================

void RadHACApKManager::Compute3x3Block_OnDemand(int elem_i, int elem_j, double* N_mat) const {
    // Compute interaction from element j to observation at element i center
    // using B_comp() directly (same approach as SetupInteractMatrix)
    //
    // IMPORTANT: Returns -N to match system matrix A = -N + diag(1/chi)

    std::memset(N_mat, 0, 9 * sizeof(double));

    if (!m_interaction || elem_i < 0 || elem_i >= m_n_elem ||
        elem_j < 0 || elem_j >= m_n_elem) {
        return;
    }

    radTg3dRelax* elem_row = m_interaction->g3dRelaxPtrVect[elem_i];
    radTg3dRelax* elem_col = m_interaction->g3dRelaxPtrVect[elem_j];
    if (!elem_row || !elem_col) return;

    // Get observation point (center of element i) in global coordinates
    TVector3d ObsPoiVect = elem_row->ReturnCentrPoint();

    // Set up field computation with interaction keys
    radTFieldKey FieldKeyInteract;
    FieldKeyInteract.B_ = FieldKeyInteract.H_ = FieldKeyInteract.PreRelax_ = 1;

    TVector3d ZeroVect(0., 0., 0.);
    radTField Field(FieldKeyInteract, m_interaction->CompCriterium, ObsPoiVect,
                   ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
    Field.AmOfIntrctElemWithSym = m_interaction->CountRelaxElemsWithSym();

    // Compute field contribution from element j
    elem_col->B_comp(&Field);

    // Field.B = response to unit Mx (dH/dMx)
    // Field.H = response to unit My (dH/dMy)
    // Field.A = response to unit Mz (dH/dMz)
    //
    // For the system matrix, we need -N:
    // N_mat[row * 3 + col] = -dH_row/dM_col

    // Row 0 (Hx response): -dHx/dMx, -dHx/dMy, -dHx/dMz
    N_mat[0] = -Field.B.x;  // -dHx/dMx
    N_mat[1] = -Field.H.x;  // -dHx/dMy
    N_mat[2] = -Field.A.x;  // -dHx/dMz
    // Row 1 (Hy response): -dHy/dMx, -dHy/dMy, -dHy/dMz
    N_mat[3] = -Field.B.y;  // -dHy/dMx
    N_mat[4] = -Field.H.y;  // -dHy/dMy
    N_mat[5] = -Field.A.y;  // -dHy/dMz
    // Row 2 (Hz response): -dHz/dMx, -dHz/dMy, -dHz/dMz
    N_mat[6] = -Field.B.z;  // -dHz/dMx
    N_mat[7] = -Field.H.z;  // -dHz/dMy
    N_mat[8] = -Field.A.z;  // -dHz/dMz
}

//=========================================================================
// Compute3x3BlockFast: ELF-style fast 3x3 block computation for tetrahedra
// Uses pre-computed geometry arrays (no object access, no B_comp overhead)
//
// The 3x3 interaction matrix N[i][j] represents how magnetization M_j creates
// demagnetizing field H at element i's center:
//   H(r_i) = sum_faces [ sigma_f * integral_triangle H_field dA ] / (4*pi)
// where sigma_f = M_j dot n_f (surface charge density)
//
// For PreRelax mode, we compute dH/dM for each unit M direction.
//=========================================================================

void RadHACApKManager::Compute3x3BlockFast(int elem_i, int elem_j, double* N_mat) const {
    std::memset(N_mat, 0, 9 * sizeof(double));

    if (!m_geometry_3dof_ready || elem_i < 0 || elem_i >= m_n_elem ||
        elem_j < 0 || elem_j >= m_n_elem) {
        return;
    }

    // Observation point: center of element i
    const double* obs = &m_tetra_centers[elem_i * 3];

    // Column element center (for point charge cancellation)
    const double* col_center = &m_tetra_centers[elem_j * 3];

    // Unit magnetization vectors
    const double M_x[3] = {1.0, 0.0, 0.0};
    const double M_y[3] = {0.0, 1.0, 0.0};
    const double M_z[3] = {0.0, 0.0, 1.0};

    // Accumulate H field for each unit M direction
    double H_from_Mx[3] = {0.0, 0.0, 0.0};
    double H_from_My[3] = {0.0, 0.0, 0.0};
    double H_from_Mz[3] = {0.0, 0.0, 0.0};

    // Track total magnetic charge for centroid cancellation
    double total_charge_Mx = 0.0;
    double total_charge_My = 0.0;
    double total_charge_Mz = 0.0;

    // Process each of the 4 triangular faces
    for (int f = 0; f < 4; f++) {
        int fnIdx = (elem_j * 4 + f) * 3;
        const double* n_f = &m_tetra_face_normals[fnIdx];

        // Surface charge density sigma = M dot n for each unit M
        double sigma_Mx = M_x[0]*n_f[0] + M_x[1]*n_f[1] + M_x[2]*n_f[2];
        double sigma_My = M_y[0]*n_f[0] + M_y[1]*n_f[1] + M_y[2]*n_f[2];
        double sigma_Mz = M_z[0]*n_f[0] + M_z[1]*n_f[1] + M_z[2]*n_f[2];

        // Accumulate total charge for each M direction
        double area = m_tetra_face_areas[elem_j * 4 + f];
        total_charge_Mx += sigma_Mx * area;
        total_charge_My += sigma_My * area;
        total_charge_Mz += sigma_Mz * area;

        // Get face vertices
        int fvIdx = (elem_j * 4 + f) * 3 * 3;
        const double* V0 = &m_tetra_face_vertices[fvIdx + 0];
        const double* V1 = &m_tetra_face_vertices[fvIdx + 3];
        const double* V2 = &m_tetra_face_vertices[fvIdx + 6];

        // Compute H field from this face (using shared triangle formula)
        // FieldFromChargedTriangleFast returns field WITHOUT 4pi divisor
        double H_f[3];

        // H from Mx contribution
        if (fabs(sigma_Mx) > 1e-20) {
            FieldFromChargedTriangleFast(obs, V0, V1, V2, sigma_Mx, H_f);
            H_from_Mx[0] += H_f[0];
            H_from_Mx[1] += H_f[1];
            H_from_Mx[2] += H_f[2];
        }

        // H from My contribution
        if (fabs(sigma_My) > 1e-20) {
            FieldFromChargedTriangleFast(obs, V0, V1, V2, sigma_My, H_f);
            H_from_My[0] += H_f[0];
            H_from_My[1] += H_f[1];
            H_from_My[2] += H_f[2];
        }

        // H from Mz contribution
        if (fabs(sigma_Mz) > 1e-20) {
            FieldFromChargedTriangleFast(obs, V0, V1, V2, sigma_Mz, H_f);
            H_from_Mz[0] += H_f[0];
            H_from_Mz[1] += H_f[1];
            H_from_Mz[2] += H_f[2];
        }
    }

    // Add point charge cancellation at centroid
    // H_point = Q * r / (4*pi * |r|^3) where Q = -total_charge (to cancel far-field)
    double r[3] = {obs[0] - col_center[0], obs[1] - col_center[1], obs[2] - col_center[2]};
    double dist_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    double dist = sqrt(dist_sq);

    if (dist > 1e-15) {
        double inv_dist3 = 1.0 / (dist * dist_sq);  // NO 4pi here (applied at end)

        // Point charge Q = -total_charge cancels far-field
        double H_pt_Mx = -total_charge_Mx * inv_dist3;
        double H_pt_My = -total_charge_My * inv_dist3;
        double H_pt_Mz = -total_charge_Mz * inv_dist3;

        H_from_Mx[0] += H_pt_Mx * r[0];
        H_from_Mx[1] += H_pt_Mx * r[1];
        H_from_Mx[2] += H_pt_Mx * r[2];

        H_from_My[0] += H_pt_My * r[0];
        H_from_My[1] += H_pt_My * r[1];
        H_from_My[2] += H_pt_My * r[2];

        H_from_Mz[0] += H_pt_Mz * r[0];
        H_from_Mz[1] += H_pt_Mz * r[1];
        H_from_Mz[2] += H_pt_Mz * r[2];
    }

    // Apply 1/(4*pi) factor and sign flip (system matrix uses -N)
    // N_mat[row * 3 + col] = -dH_row/dM_col / (4*pi)
    // Row 0 (Hx response): -dHx/dMx, -dHx/dMy, -dHx/dMz
    N_mat[0] = -H_from_Mx[0] * RadConst::INV_FOUR_PI;  // -dHx/dMx
    N_mat[1] = -H_from_My[0] * RadConst::INV_FOUR_PI;  // -dHx/dMy
    N_mat[2] = -H_from_Mz[0] * RadConst::INV_FOUR_PI;  // -dHx/dMz
    // Row 1 (Hy response): -dHy/dMx, -dHy/dMy, -dHy/dMz
    N_mat[3] = -H_from_Mx[1] * RadConst::INV_FOUR_PI;  // -dHy/dMx
    N_mat[4] = -H_from_My[1] * RadConst::INV_FOUR_PI;  // -dHy/dMy
    N_mat[5] = -H_from_Mz[1] * RadConst::INV_FOUR_PI;  // -dHy/dMz
    // Row 2 (Hz response): -dHz/dMx, -dHz/dMy, -dHz/dMz
    N_mat[6] = -H_from_Mx[2] * RadConst::INV_FOUR_PI;  // -dHz/dMx
    N_mat[7] = -H_from_My[2] * RadConst::INV_FOUR_PI;  // -dHz/dMy
    N_mat[8] = -H_from_Mz[2] * RadConst::INV_FOUR_PI;  // -dHz/dMz
}

//=========================================================================
// Mixed element methods: 3x6 and 6x3 blocks for tetra-hex interactions
// Following ELF_MAGIC convention for mixed element meshes
//=========================================================================

double RadHACApKManager::GetMixed3x6Element(int elem_tetra, int elem_hex, int comp, int face) const {
    // 3x6 block: tetra row (3DOF), hex column (6DOF)
    // K(comp, face) = H_field_comp at tetra center from unit sigma on hex face
    // Returns -K/(4*pi) to match the sign convention of the system matrix

    if (!m_interaction) return 0.0;

    // For mixed elements, we need to access the pre-computed interaction matrix
    // The variable DOF matrix m_flatInteractMatrix stores blocks sequentially
    if (m_interaction->m_flatInteractMatrix.empty()) {
        // Fall back to computing on-demand (slower but correct)
        double K_mat[18];
        Compute3x6Block(elem_tetra, elem_hex, K_mat);
        return K_mat[comp * 6 + face];
    }

    // Access from pre-computed flat matrix
    int offset_i = m_dof_offset[elem_tetra];
    int offset_j = m_dof_offset[elem_hex];
    int total_dof = m_interaction->m_totalDOF;

    // COLUMN-MAJOR access: A(row, col) at [col * total_dof + row]
    return m_interaction->m_flatInteractMatrix[(offset_j + face) * total_dof + (offset_i + comp)];
}

double RadHACApKManager::GetMixed6x3Element(int elem_hex, int elem_tetra, int face, int comp) const {
    // 6x3 block: hex row (6DOF), tetra column (3DOF)
    // K(face, comp) = normal_face dot N_mat(:, comp) where N_mat is demagnetization tensor
    // Returns -K/(4*pi) to match the sign convention

    if (!m_interaction) return 0.0;

    // Access from pre-computed flat matrix if available
    if (m_interaction->m_flatInteractMatrix.empty()) {
        double K_mat[18];
        Compute6x3Block(elem_hex, elem_tetra, K_mat);
        return K_mat[face * 3 + comp];
    }

    int offset_i = m_dof_offset[elem_hex];
    int offset_j = m_dof_offset[elem_tetra];
    int total_dof = m_interaction->m_totalDOF;

    // COLUMN-MAJOR access
    return m_interaction->m_flatInteractMatrix[(offset_j + comp) * total_dof + (offset_i + face)];
}

void RadHACApKManager::Compute3x6Block(int elem_tetra, int elem_hex, double* K_mat) const {
    // 3x6 block: H-field at tetra center from hex face charges
    // K(comp, face) = H_comp at tetra center due to unit sigma on hex face
    // Following rad_interaction.cpp 3x6 block implementation

    std::memset(K_mat, 0, 18 * sizeof(double));

    if (!m_interaction) return;

    radTg3dRelax* elem_row = m_interaction->g3dRelaxPtrVect[elem_tetra];
    radTg3dRelax* elem_col = m_interaction->g3dRelaxPtrVect[elem_hex];
    if (!elem_row || !elem_col) return;

    radTPolyhedron* poly_col = dynamic_cast<radTPolyhedron*>(elem_col);
    if (!poly_col || !poly_col->Use6DOF_MSC) return;

    // Observation point: tetra center
    TVector3d obs = elem_row->CentrPoint;

    for (int face_j = 0; face_j < 6; face_j++) {
        // Field from unit sigma on face j
        TVector3d H_face = poly_col->FieldFromQuadFace(obs, face_j, 1.0);

        // Point charge contribution (m = -sigma * area)
        double unit_charge = -1.0 * poly_col->FaceArea[face_j];
        TVector3d H_point = poly_col->FieldFromPointCharge(obs, unit_charge);

        TVector3d H_total;
        H_total.x = H_face.x + H_point.x;
        H_total.y = H_face.y + H_point.y;
        H_total.z = H_face.z + H_point.z;

        // Store with sign flip (-K/(4*pi))
        // Row-major: K_mat[comp * 6 + face]
        K_mat[0 * 6 + face_j] = -H_total.x * RadConst::INV_FOUR_PI;  // Mx
        K_mat[1 * 6 + face_j] = -H_total.y * RadConst::INV_FOUR_PI;  // My
        K_mat[2 * 6 + face_j] = -H_total.z * RadConst::INV_FOUR_PI;  // Mz
    }
}

void RadHACApKManager::Compute6x3Block(int elem_hex, int elem_tetra, double* K_mat) const {
    // 6x3 block: normal dot demagnetization tensor at hex eval points from tetra
    // K(face, comp) = normal_face dot N_mat(:, comp)
    // Following rad_interaction.cpp 6x3 block implementation

    std::memset(K_mat, 0, 18 * sizeof(double));

    if (!m_interaction) return;

    radTg3dRelax* elem_row = m_interaction->g3dRelaxPtrVect[elem_hex];
    radTg3dRelax* elem_col = m_interaction->g3dRelaxPtrVect[elem_tetra];
    if (!elem_row || !elem_col) return;

    radTPolyhedron* poly_row = dynamic_cast<radTPolyhedron*>(elem_row);
    if (!poly_row || !poly_row->Use6DOF_MSC) return;

    radTFieldKey FieldKeyInteract;
    FieldKeyInteract.B_ = FieldKeyInteract.H_ = FieldKeyInteract.PreRelax_ = 1;

    for (int face_i = 0; face_i < 6; face_i++) {
        // Yano-Sugahara evaluation point: midpoint between face center and element center
        TVector3d EvalPt;
        EvalPt.x = 0.5 * (poly_row->FaceCenter[face_i].x + poly_row->CentrPoint.x);
        EvalPt.y = 0.5 * (poly_row->FaceCenter[face_i].y + poly_row->CentrPoint.y);
        EvalPt.z = 0.5 * (poly_row->FaceCenter[face_i].z + poly_row->CentrPoint.z);

        // Compute demagnetization tensor at this point
        radTField Field(FieldKeyInteract, m_interaction->CompCriterium, EvalPt,
                       TVector3d(0., 0., 0.), TVector3d(0., 0., 0.),
                       TVector3d(0., 0., 0.), TVector3d(0., 0., 0.), 0.);
        Field.AmOfIntrctElemWithSym = m_interaction->CountRelaxElemsWithSym();

        elem_col->B_comp(&Field);

        // N_mat columns from Field (response to unit M)
        // Field.B = response to unit Mx, Field.H = My, Field.A = Mz
        TVector3d& n = poly_row->FaceNormal[face_i];

        // K(face_i, Mj) = normal dot N_mat(:, j) / (4*pi)
        double K_Mx = n.x * Field.B.x + n.y * Field.B.y + n.z * Field.B.z;
        double K_My = n.x * Field.H.x + n.y * Field.H.y + n.z * Field.H.z;
        double K_Mz = n.x * Field.A.x + n.y * Field.A.y + n.z * Field.A.z;

        // Store with sign flip (-K/(4*pi))
        // Row-major: K_mat[face * 3 + comp]
        K_mat[face_i * 3 + 0] = -K_Mx * RadConst::INV_FOUR_PI;
        K_mat[face_i * 3 + 1] = -K_My * RadConst::INV_FOUR_PI;
        K_mat[face_i * 3 + 2] = -K_Mz * RadConst::INV_FOUR_PI;
    }
}

//=========================================================================
// End of rad_hacapk.cpp
//=========================================================================
