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
#include <cmath>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <chrono>

// Include C++ compatible HACApK wrapper header
extern "C" {
#include "../ext/HACApK_LH-Cimplm/cHACApK_cpp.h"
}

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

double ComputeEntry(int i, int j) {
    // i, j are 1-based indices from HACApK (converted to 0-based)
    // Matrix element: A(i,j) = N(i,j) + delta_ij/chi_i
    // where N already contains -K/(4*pi) from GetInteractionMatrixElement()
    // So the equation is: (-K/(4pi) + 1/chi * I) * sigma = H_ext_n

    if (g_currentManager == nullptr) {
        std::cerr << "[HACApK] Error: g_currentManager is null in ComputeEntry" << std::endl;
        return 0.0;
    }

    int i0 = i - 1;  // Convert to 0-based
    int j0 = j - 1;

    // Bounds check
    int ndof = g_currentManager->GetNDOF();
    if (i0 < 0 || i0 >= ndof || j0 < 0 || j0 >= ndof) {
        std::cerr << "[HACApK] Error: Invalid DOF indices in ComputeEntry: i=" << i
                  << " j=" << j << " ndof=" << ndof << std::endl;
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
    , m_cache_access_counter(0)
{
    // Initialize block cache
    m_block_cache.resize(BLOCK_CACHE_SIZE);
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
    // HACApK is only supported for 6DOF MSC hexahedra
    if (!m_interaction) return;

    m_n_elem = m_interaction->AmOfMainElem;
    m_coordinates.resize(m_n_elem * 3);
    m_dof_offset.resize(m_n_elem + 1);

    // Verify all elements are 6DOF MSC hexahedra
    bool has_variable_dof = m_interaction->HasVariableDOF();

    if (!has_variable_dof) {
        // Standard 3DOF elements - HACApK not supported
        std::cerr << "[HACApK] Error: HACApK is only supported for 6DOF MSC hexahedra, not 3DOF elements" << std::endl;
        m_ndof = 0;
        m_nffc = 0;
        m_is_6dof = false;
        return;
    }

    // Check each element is 6DOF
    int total_dof = 0;
    int n_hex = 0;

    for (int i = 0; i < m_n_elem; i++) {
        m_dof_offset[i] = total_dof;
        int elem_dof = m_interaction->GetElementDOF(i);
        total_dof += elem_dof;

        if (elem_dof == 6) {
            n_hex++;
        } else {
            std::cerr << "[HACApK] Error: Element " << i << " has " << elem_dof
                      << " DOF, expected 6 (6DOF MSC hexahedra only)" << std::endl;
            m_ndof = 0;
            m_nffc = 0;
            m_is_6dof = false;
            return;
        }
    }
    m_dof_offset[m_n_elem] = total_dof;
    m_ndof = total_dof;
    m_nffc = 6;
    m_is_6dof = true;

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

    // Verify 6DOF was configured correctly
    if (!m_is_6dof || m_ndof == 0) {
        std::cerr << "[HACApK] Error: HACApK requires 6DOF MSC hexahedral elements" << std::endl;
        return false;
    }

    // Initialize inverse susceptibility with current values
    // 6DOF MSC: all 6 DOF per element use the same chi value
    m_inv_chi.resize(m_ndof);

    double* FlatField = m_interaction->GetFlatFieldArray();

    for (int i = 0; i < m_n_elem; i++) {
        radTg3dRelax* g3dRelaxPtr = m_interaction->g3dRelaxPtrVect[i];
        radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

        // Estimate H from field array
        int offset = m_dof_offset[i];
        TVector3d H_est(0., 0., FlatField ? FlatField[offset] : 0.);
        TMatrix3d KsiTensor;
        TVector3d MrVect;
        MaterPtr->DefineInstantKsiTensor(H_est, KsiTensor, MrVect);

        double chi = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
        if (chi < 1.0e-6) chi = 1.0e-6;
        double inv_chi_val = 1.0 / chi;

        // All 6 DOF use the same 1/chi
        for (int k = 0; k < 6; k++) {
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
    // m_nffc: 3 for tetrahedra (Mx, My, Mz), 6 for hexahedra (sigma per face)
    int ndim = 3;  // Spatial dimension

    int result = HACApK_build_hmatrix_wrapper(
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

    // Estimate compression ratio
    int64_t full_size = (int64_t)m_ndof * m_ndof;
    // Rough estimate: assume average rank and block size
    int64_t hmat_size = m_stats.n_lowrank * (m_ndof / m_stats.n_leaves) * m_stats.max_rank * 2;
    hmat_size += m_stats.n_dense * (m_ndof / m_stats.n_leaves) * (m_ndof / m_stats.n_leaves);
    m_stats.compression = (full_size > 0) ? (double)hmat_size / (double)full_size : 1.0;

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

    // Update dense diagonal blocks in H-matrix
    // This follows ELF_MAGIC's HACApK_update_diagonal_omp pattern:
    // - Only dense diagonal blocks (ltmtx==2 && nstrtl==nstrtt) are recomputed
    // - Low-rank off-diagonal blocks remain unchanged (geometry doesn't change)
    // - The callback function cHACApK_entry_ij will use the updated inv_chi values
    HACApK_update_diagonal_wrapper(m_leafmtxp, m_control, cHACApK_entry_ij);
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

// Constants for 6DOF MSC field computation
static const double PI_HACAPK = 3.14159265358979323846;
static const double INV_4PI_HACAPK = 1.0 / (4.0 * PI_HACAPK);

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

            // Store -K_ij / (4*pi) in row-major order
            K_mat[fi * 6 + fj] = -K_ij * INV_4PI_HACAPK;
        }
    }
}

//=========================================================================
// GetCached6x6Element: LRU cache lookup for 6x6 blocks (ELF-style)
// Thread-safe implementation for OpenMP parallelization in HACApK
//=========================================================================

double RadHACApKManager::GetCached6x6Element(int elem_i, int elem_j, int face_i, int face_j) const {
    // TEMPORARY: Bypass cache completely for debugging
    {
        double K_mat[36];
        Compute6x6Block(elem_i, elem_j, K_mat);
        return K_mat[face_i * 6 + face_j];
    }

#if 0  // Cache disabled for debugging
    // Ensure cache is initialized
    if (m_block_cache.empty()) {
        // Cannot resize const vector, compute directly
        double K_mat[36];
        Compute6x6Block(elem_i, elem_j, K_mat);
        return K_mat[face_i * 6 + face_j];
    }

    double result = 0.0;
    bool found = false;

    // Critical section for thread-safe cache access
    #pragma omp critical(hacapk_cache)
    {
        // Search cache
        int cache_size = (int)m_block_cache.size();
        for (int c = 0; c < cache_size; c++) {
            if (m_block_cache[c].elem_i == elem_i && m_block_cache[c].elem_j == elem_j) {
                // Cache hit
                m_cache_access_counter++;
                m_block_cache[c].access_count = m_cache_access_counter;
                result = m_block_cache[c].K_mat[face_i * 6 + face_j];
                found = true;
                break;
            }
        }

        if (!found) {
            // Cache miss - find LRU slot
            int lru_slot = 0;
            int min_access = m_block_cache[0].access_count;
            for (int c = 1; c < cache_size; c++) {
                if (m_block_cache[c].elem_i < 0) {
                    lru_slot = c;
                    break;
                }
                if (m_block_cache[c].access_count < min_access) {
                    min_access = m_block_cache[c].access_count;
                    lru_slot = c;
                }
            }

            // Compute and cache
            Compute6x6Block(elem_i, elem_j, m_block_cache[lru_slot].K_mat);
            m_block_cache[lru_slot].elem_i = elem_i;
            m_block_cache[lru_slot].elem_j = elem_j;
            m_cache_access_counter++;
            m_block_cache[lru_slot].access_count = m_cache_access_counter;

            result = m_block_cache[lru_slot].K_mat[face_i * 6 + face_j];
        }
    }

    return result;
#endif  // Cache disabled for debugging
}

//=========================================================================
// GetInteractionMatrixElement: Optimized with O(1) lookup and LRU cache
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
    int face_i = m_dof_to_local[dof_i];
    int elem_j = m_dof_to_elem[dof_j];
    int face_j = m_dof_to_local[dof_j];

    // Safety check for element indices
    if (elem_i < 0 || elem_i >= m_n_elem || elem_j < 0 || elem_j >= m_n_elem) {
        return 0.0;
    }

    // Get cached 6x6 block element
    return GetCached6x6Element(elem_i, elem_j, face_i, face_j);
}

//=========================================================================
// End of rad_hacapk.cpp
//=========================================================================
