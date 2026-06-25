/*-------------------------------------------------------------------------
*
* File name:      rad_hacapk.cpp
*
* Project:        RADIA
*
* Description:    HACApK (H-matrix with ACA+) interface for BiCGSTAB solver
*                 Implementation of RadHACApKBase / RadHACApKMMMManager and callback functions
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
#include "rad_poly_analytical.h"
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <algorithm>
#include <cstring>
#include <mkl.h>   // cblas_dgemm + LAPACKE_dsyevd for the antisym-IMA plane-slab null space (SVD-free)
#include <cstdio>
#include <cstdlib>  // std::getenv / std::atoi (RADIA_LS_BLOCK_CELLS research toggle)
#include <iostream>
#include <chrono>
#include <atomic>
#include "rad_parallel.h"

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
    // Global callback state shared across all TaskManager worker threads.
    // NOT thread_local: HACApK calls cHACApK_entry_ij from ngcore::ParallelFor
    // (TaskManager) worker threads, which must see the same manager/invChi/interaction
    // set by the main thread.
    // Thread safety for concurrent BuildHMatrix calls (multiple Python threads)
    // is ensured by Python's GIL; standalone C++ use would require a mutex.
    RadHACApKBase* g_currentManager = nullptr;
    std::vector<double> g_invChi;
    radTInteraction* g_interaction = nullptr;
    int g_nElem = 0;
    int g_nffc = 3;  // DOF per element (default 3 for standard elements)
    // Note: lod (DOF permutation array) is now accessed via C-side accessors:
    //   HACApK_get_current_lod() and HACApK_get_current_lod_size()
    // Set during H-matrix build by HACApK_build_hmatrix_varDOF_wrapper

    // FIX (2025-02-04): Generation counter for thread-local cache invalidation
    // Incremented each time a new HACApK manager is created or BuildHMatrix is called.
    // Thread-local caches check this counter to invalidate stale entries.
    // Using a counter instead of pointer comparison prevents issues with memory reuse.
    std::atomic<uint64_t> g_hacapk_generation{0};
}

namespace RadHACApKCallback {

void SetCurrentManager(RadHACApKBase* manager) {
    g_currentManager = manager;
}

RadHACApKBase* GetCurrentManager() {
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

void ClearGlobalState() {
    // Clear all global callback state to prevent interference between solves
    g_currentManager = nullptr;
    g_interaction = nullptr;
    g_nElem = 0;
    g_nffc = 3;
    g_invChi.clear();
    // Note: g_hacapk_generation is NOT reset here - it continues incrementing
}

// FIX (2025-02-04): Get current generation for cache invalidation
uint64_t GetGeneration() {
    return g_hacapk_generation.load(std::memory_order_acquire);
}

// FIX (2025-02-04): Increment generation to invalidate all thread-local caches
void IncrementGeneration() {
    g_hacapk_generation.fetch_add(1, std::memory_order_release);
}

double ComputeEntry(int i, int j) {
    // i, j are 1-based ORIGINAL indices from HACApK.
    // The lod conversion is already done by cHACApK_fill_leafmtx_hyp:
    //   val = cHACApK_entry_ij(lodl[permuted_pos], lodt[permuted_pos], i_bemv)
    // so we receive original indices, NOT permuted indices.
    //
    // The kernel-specific system-matrix convention is delegated to
    // RadHACApKBase::ComputeSystemEntry so that each subclass (MSC,
    // PEEC, future BEM) can store exactly what HACApK needs.

    if (g_currentManager == nullptr) {
        std::cerr << "[HACApK] Error: g_currentManager is null in ComputeEntry" << std::endl;
        return 0.0;
    }

    int ndof = g_currentManager->GetNDOF();

    // Direct 0-based conversion (indices are already original, NOT permuted)
    int i0 = i - 1;
    int j0 = j - 1;

    if (i0 < 0 || i0 >= ndof || j0 < 0 || j0 >= ndof) {
        std::cerr << "[HACApK] Error: Invalid DOF indices: i0=" << i0
                  << " j0=" << j0 << " ndof=" << ndof << std::endl;
        return 0.0;
    }

    return g_currentManager->ComputeSystemEntry(i0, j0);
}

}  // namespace RadHACApKCallback

//=========================================================================
// C callback function for HACApK
// This is the function HACApK calls to get matrix elements
//=========================================================================

extern "C" {

// Optional per-call entry-function override (implements the HACApK_set_entry_func
// API declared in cHACApK_cpp.h, previously unimplemented).  When non-null,
// cHACApK_acaplus / fill routines fetch matrix entries from this kernel instead
// of the default MMM/MSC system matrix.  This lets callers (e.g. the
// stream-function (ACA+)+TSVD solver) factor an arbitrary rectangular kernel
// block with HACApK's ACA+ -- keeping ACA+ a single source of truth instead of
// re-porting it.  Default null => unchanged MMM behaviour.  Set/cleared
// synchronously around one factorization (GIL-serialized; no concurrent MMM
// build), so a plain (non-thread_local) pointer is sufficient.
static HACApK_entry_func g_entry_override = NULL;

void HACApK_set_entry_func(HACApK_entry_func func) { g_entry_override = func; }
void HACApK_clear_entry_func(void) { g_entry_override = NULL; }

double cHACApK_entry_ij(int i, int j, int i_bemv) {
    if (g_entry_override != NULL) return g_entry_override(i, j, i_bemv);
    (void)i_bemv;  // Unused in Radia
    return RadHACApKCallback::ComputeEntry(i, j);
}

}  // extern "C"

//=========================================================================
// RadHACApKBase Implementation (kernel-agnostic H-matrix lifecycle)
//=========================================================================

RadHACApKBase::RadHACApKBase()
    : m_leafmtxp(nullptr)
    , m_control(nullptr)
    , m_valid(false)
    , m_ndof(0)
    , m_n_elem(0)
    , m_diag_cached(false)
    , m_defl_alpha(0.0)
    , m_defl_nplaq(0)
{
}

RadHACApKBase::~RadHACApKBase() {
    FreeResources();
}

//=========================================================================
// RadHACApKMMMManager Implementation (MMM / MSC kernel)
//=========================================================================

RadHACApKMMMManager::RadHACApKMMMManager(radTInteraction* interaction)
    : RadHACApKBase()
    , m_interaction(interaction)
    , m_nffc(3)
    , m_geometry_3dof_ready(false)
    , m_flat_N_ready(false)
{
    // Hash-based cache is initialized automatically
}

RadHACApKMMMManager::~RadHACApKMMMManager() {}

//=========================================================================
// MSC system-matrix convention: A(i, j) = -N(i, j) + delta_ij / chi_i
// where N is the physical demagnetization tensor entry returned by
// GetInteractionMatrixElement (MSC/MMM convention: +N).
//=========================================================================

double RadHACApKMMMManager::ComputeSystemEntry(int dof_i, int dof_j) const {
    double N_val = GetInteractionMatrixElement(dof_i, dof_j);
    double A_val = -N_val;
    if (dof_i == dof_j && dof_i < (int)m_inv_chi.size()) {
        A_val += m_inv_chi[dof_i];
    }
    return A_val;
}

void RadHACApKBase::FreeResources() {
    // Clear global callback state first (prevents stale state in next solve)
    RadHACApKCallback::ClearGlobalState();

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

    // Reset HACApK global state (persistent matvec buffers, lod)
    HACApK_reset_global_state();

    m_valid = false;
}

//=========================================================================
void RadHACApKMMMManager::ExtractCoordinates() {
    // Extract element center coordinates for clustering
    // Supports both 3DOF tetrahedra and 6DOF MSC hexahedra
    if (!m_interaction) return;

    m_n_elem = m_interaction->AmOfMainElem;
    m_coordinates.resize(m_n_elem * 3);
    m_dof_offset.resize(m_n_elem + 1);

    // MMM-only manager: uniform 3-DOF tetrahedra.  Surface-charge MSC (hex 6-DOF, wedge 5-DOF) is
    // solved by the multipole-moment MMM path (SolveGen forces pure-MSC to the LU/Picard moment driver) and mixed
    // MMM+MSC is rejected fail-loud in MakeAutoRelax (Error204), so ONLY 3-DOF tets reach this manager
    // (the EIEM2 surface-charge collocation kernels were retired in Phase 3b).
    int total_dof = 0;
    for (int i = 0; i < m_n_elem; i++) {
        m_dof_offset[i] = total_dof;
        int elem_dof = m_interaction->GetElementDOF(i);
        if (elem_dof != 3) {
            std::cerr << "[HACApK] Error: Element " << i << " has " << elem_dof
                      << " DOF; the MMM (HACApK) manager handles 3-DOF tetrahedra only "
                      << "(surface-charge MSC uses the multipole-moment MMM solver)" << std::endl;
            m_ndof = 0;
            m_nffc = 0;
            return;
        }
        total_dof += elem_dof;
    }

    m_dof_offset[m_n_elem] = total_dof;
    m_ndof = total_dof;
    m_nffc = 3;   // uniform 3-DOF (tetrahedron MMM)

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

void RadHACApKMMMManager::BuildDOFLookupTable() {
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
// PrecomputeGeometry3DOF: ELF-style pre-computed geometry for 3DOF tetrahedra
// Extracts all tetrahedron vertices and face geometry into contiguous arrays
// This avoids calling B_comp() which has significant overhead during H-matrix build
//=========================================================================

void RadHACApKMMMManager::PrecomputeGeometry3DOF() {
    if (m_geometry_3dof_ready || m_n_elem == 0 || !m_interaction) return;

    // Allocate arrays for tetrahedra (4 triangular faces, 3 vertices each)
    m_tetra_centers.resize(m_n_elem * 3);
    m_tetra_face_vertices.resize(m_n_elem * 4 * 3 * 3);  // 4 faces, 3 vertices, 3 coords
    m_tetra_face_normals.resize(m_n_elem * 4 * 3);       // 4 faces, 3 coords (outward normals)
    m_tetra_face_areas.resize(m_n_elem * 4);              // 4 faces

    // Parallel per-element fill: each iteration writes to disjoint slices of
    // the pre-allocated flat arrays (m_tetra_centers, m_tetra_face_vertices,
    // m_tetra_face_normals, m_tetra_face_areas), so no synchronisation needed.
    ngcore::ParallelFor(ngcore::IntRange(m_n_elem), [&](size_t e_size) {
        int e = (int)e_size;
        radTg3dRelax* elem = m_interaction->g3dRelaxPtrVect[e];
        if (!elem) return;

        radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(elem);
        if (!poly || poly->AmOfFaces != 4) return;  // Skip non-tetrahedra

        // Store element center
        int cIdx = e * 3;
        m_tetra_centers[cIdx + 0] = poly->CentrPoint.x;
        m_tetra_centers[cIdx + 1] = poly->CentrPoint.y;
        m_tetra_centers[cIdx + 2] = poly->CentrPoint.z;

        // Store face data for each of the 4 triangular faces
        for (int f = 0; f < 4; f++) {
            const radTHandlePgnAndTrans& hpt = poly->VectHandlePgnAndTrans[f];
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
    });

    m_geometry_3dof_ready = true;
}

//=========================================================================
// PrecomputeFlatInteractMatrix: Flatten InteractMatrix for 3DOF tetrahedra
// Converts 2D pointer array TMatrix3df** to contiguous double array
// This eliminates pointer chasing during matrix element access
//=========================================================================

void RadHACApKMMMManager::PrecomputeFlatInteractMatrix() {
    if (m_flat_N_ready || m_n_elem == 0 || !m_interaction) return;
    if (!m_interaction->InteractMatrix) {
        return;  // InteractMatrix not computed
    }

    // Allocate flat storage: n_elem * n_elem * 9 doubles (3x3 block per element pair)
    int64_t total_size = (int64_t)m_n_elem * m_n_elem * 9;
    m_flat_N_data.resize(total_size);

    // Copy from InteractMatrix[i][j] to flat array as +N (physical quantity).
    // ComputeEntry() handles the sign flip to -N for the system matrix.
    //
    // InteractMatrix[i][j] stores TMatrix3df where:
    //   Str0 = dH/dMx, Str1 = dH/dMy, Str2 = dH/dMz (COLUMN vectors)
    // We store row-major: N[row][col] = dH_row / dM_col
    int total_collapse = m_n_elem * m_n_elem;
    ngcore::ParallelFor(ngcore::IntRange(total_collapse), [&](size_t idx)
    {
        int i = (int)(idx / m_n_elem);
        int j = (int)(idx % m_n_elem);
            const TMatrix3df& M = m_interaction->InteractMatrix[i][j];
            int64_t base_idx = ((int64_t)i * m_n_elem + j) * 9;

            m_flat_N_data[base_idx + 0] = static_cast<double>(M.Str0.x);  // dHx/dMx
            m_flat_N_data[base_idx + 1] = static_cast<double>(M.Str1.x);  // dHx/dMy
            m_flat_N_data[base_idx + 2] = static_cast<double>(M.Str2.x);  // dHx/dMz
            m_flat_N_data[base_idx + 3] = static_cast<double>(M.Str0.y);  // dHy/dMx
            m_flat_N_data[base_idx + 4] = static_cast<double>(M.Str1.y);  // dHy/dMy
            m_flat_N_data[base_idx + 5] = static_cast<double>(M.Str2.y);  // dHy/dMz
            m_flat_N_data[base_idx + 6] = static_cast<double>(M.Str0.z);  // dHz/dMx
            m_flat_N_data[base_idx + 7] = static_cast<double>(M.Str1.z);  // dHz/dMy
            m_flat_N_data[base_idx + 8] = static_cast<double>(M.Str2.z);  // dHz/dMz
    });

    m_flat_N_ready = true;
}

bool RadHACApKBase::BuildHMatrix(const RadHACApKParams& params) {
    // TaskManager self-wrap (AGENTS.md "Parallelization: NGSolve TaskManager"): the H-matrix leaf
    // fill runs ngcore::ParallelFor, which silently falls back to single-threaded when NO
    // RegionTaskManager is active.  Stand up (or reuse the caller's) pool here so EVERY
    // HACApK build -- multipole-moment MMM, HDiv, MMM/MSC, PEEC, diagnostics -- is parallel even when
    // a non-panel caller forgot `with TaskManager()`.  Nested -> reuses the caller's (no-op).
    ngcore::RegionTaskManager rtm(radia::GetMaxThreads());

    FreeResources();

    auto start_time = std::chrono::high_resolution_clock::now();

    ExtractCoordinates();

    if (m_n_elem == 0) {
        std::cerr << "[HACApK] Error: No elements" << std::endl;
        return false;
    }
    if (m_ndof == 0) {
        std::cerr << "[HACApK] Error: Invalid DOF configuration (ndof=0)" << std::endl;
        return false;
    }

    // Set global callback state (MUST happen before OnBeforeBuild / InitializeInvChi
    // so kernel-side hooks may register additional callback state if needed)
    RadHACApKCallback::SetCurrentManager(this);

    // FIX (2025-02-04): Increment generation counter to invalidate thread-local caches
    // so that matrix-element caches from previous solves are not reused.
    RadHACApKCallback::IncrementGeneration();

    // Kernel-specific precomputation (e.g. PrecomputeHexaGeometry for MSC hex)
    OnBeforeBuild();

    // Kernel-specific initial chi (for MSC: initial susceptibility from material state)
    InitializeInvChi();

    RadHACApKCallback::SetInvChi(m_inv_chi);

    // Allocate opaque structures
    m_leafmtxp = HACApK_alloc_leafmtxp();
    m_control = HACApK_alloc_lcontrol();

    if (!m_leafmtxp || !m_control) {
        std::cerr << "[HACApK] Error: Failed to allocate structures" << std::endl;
        return false;
    }

    // Build H-matrix via HACApK wrapper. Variable-DOF mode routes through the
    // dof_offset array; uniform-DOF mode uses GetUniformNFFC().
    int ndim = 3;

    auto t_hmatrix_start = std::chrono::high_resolution_clock::now();
    int result;

    if (IsVariableDOF()) {
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
        result = HACApK_build_hmatrix_wrapper(
            m_leafmtxp,
            m_control,
            m_coordinates.data(),
            m_n_elem,
            GetUniformNFFC(),
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

    int64_t hmat_bytes = 0;
    int64_t dense_bytes = 0;
    HACApK_get_memory_stats(m_leafmtxp, &hmat_bytes, &dense_bytes);

    m_stats.memory_mb = (double)hmat_bytes / (1024.0 * 1024.0);
    m_stats.dense_memory_mb = (double)dense_bytes / (1024.0 * 1024.0);
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

    // Cache diagonal elements N_ii for Jacobi preconditioner (reused every
    // BiCGSTAB iteration). Uses virtual GetInteractionMatrixElement.
    m_diag_N.resize(m_ndof);
    ngcore::ParallelFor(ngcore::IntRange(m_ndof), [&](size_t i) {
        m_diag_N[(int)i] = GetInteractionMatrixElement((int)i, (int)i);
    });
    m_diag_cached = true;

    m_valid = true;
    return true;
}

//=========================================================================
// RadHACApKMMMManager::OnBeforeBuild
// MMM (3-DOF tet) precomputation (runs AFTER ExtractCoordinates has populated
// m_nffc=3, and BEFORE HACApK callbacks).
//=========================================================================

void RadHACApKMMMManager::OnBeforeBuild() {
    if (!m_interaction) return;

    if (m_nffc != 3) {
        std::cerr << "[HACApK] Warning: MMM manager expects 3-DOF tet (nffc=" << m_nffc << ")" << std::endl;
    }

    // Register with callback (informational; ComputeEntry uses the manager directly)
    RadHACApKCallback::SetInteraction(m_interaction, m_n_elem, m_nffc);

    // ELF-style pre-computation for 3DOF tetrahedra: extract face vertices/normals for direct
    // field computation without the O(N^2) SetupInteractMatrix().
    PrecomputeGeometry3DOF();
    // Fallback: if InteractMatrix was already computed, use flat storage for O(1) element access.
    if (!m_geometry_3dof_ready && m_interaction->InteractMatrix != nullptr) {
        PrecomputeFlatInteractMatrix();
    }
}

//=========================================================================
// RadHACApKMMMManager::InitializeInvChi
// Populate m_inv_chi from material state. Matches ELF's
// initialize_chi_from_bh() for nonlinear isotropic materials (chi from the
// 2nd BH curve point) and falls back to DefineInstantKsiTensor(H=0) for
// linear materials.
//=========================================================================

void RadHACApKMMMManager::InitializeInvChi() {
    if (!m_interaction) return;
    m_inv_chi.resize(m_ndof);

    for (int i = 0; i < m_n_elem; i++) {
        radTg3dRelax* g3dRelaxPtr = m_interaction->g3dRelaxPtrVect[i];
        radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

        double chi = 1.0;
        radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);
        if (NonlinMater != nullptr) {
            chi = NonlinMater->GetInitialChi_ELF_Style();
            if (chi <= 0) chi = 1.0;
        } else {
            TMatrix3d KsiTensor;
            TVector3d MrVect;
            TVector3d H_zero(0., 0., 0.);
            MaterPtr->DefineInstantKsiTensor(H_zero, KsiTensor, MrVect);
            chi = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
        }
        if (chi < 1.0e-6) chi = 1.0e-6;
        double inv_chi_val = 1.0 / chi;

        // All DOF on this element share the same 1/chi.
        int offset = m_dof_offset[i];
        int elem_dof = m_dof_offset[i + 1] - offset;
        for (int k = 0; k < elem_dof; k++) {
            m_inv_chi[offset + k] = inv_chi_val;
        }
    }
}

void RadHACApKBase::SetDeflation(const std::vector<int>& plaq_offsets,
                                 const std::vector<int>& dofs,
                                 const std::vector<double>& signs,
                                 double alpha) {
    m_defl_offsets = plaq_offsets;
    m_defl_dofs = dofs;
    m_defl_signs = signs;
    m_defl_alpha = alpha;
    m_defl_nplaq = plaq_offsets.empty() ? 0 : (int)plaq_offsets.size() - 1;
}

void RadHACApKBase::MatVec(const std::vector<double>& x, std::vector<double>& y) {
    if (!m_valid || !m_leafmtxp || !m_control) {
        std::fill(y.begin(), y.end(), 0.0);
        return;
    }

    int nd = HACApK_leafmtxp_get_nd(m_leafmtxp);
    HACApK_matvec_wrapper(m_leafmtxp, m_control, x.data(), y.data(), nd);

    // Matrix-free loop-mode deflation: y += alpha * L (L^T x).  L L^T acts
    // only on span(L) = the null (loop) subspace, lifting its eigenvalues
    // away from zero so the converged solution is loop-free.  O(nnz) = O(N).
    if (m_defl_nplaq > 0 && m_defl_alpha != 0.0) {
        for (int p = 0; p < m_defl_nplaq; p++) {
            double c = 0.0;
            for (int k = m_defl_offsets[p]; k < m_defl_offsets[p + 1]; k++)
                c += m_defl_signs[k] * x[m_defl_dofs[k]];
            c *= m_defl_alpha;
            for (int k = m_defl_offsets[p]; k < m_defl_offsets[p + 1]; k++)
                y[m_defl_dofs[k]] += m_defl_signs[k] * c;
        }
    }
}

void RadHACApKBase::UpdateDiagonal(const std::vector<double>& inv_chi) {
    if (!m_valid || !m_leafmtxp || !m_control) return;

    // Update stored inverse susceptibility
    m_inv_chi = inv_chi;
    RadHACApKCallback::SetInvChi(inv_chi);

    // OPTIMIZATION: Use fast diagonal update (2025-12-30)
    // Only updates true diagonal entries (i==j) using pre-computed N_ii values
    // This is O(ndof) instead of O(block_size^2 * n_diag_blocks)
    // Performance: ELF 0.39s vs old Radia 4.6s (12x slower due to slow method)
    if (m_diag_cached && m_diag_N.size() == (size_t)m_ndof) {
        HACApK_update_diagonal_fast_wrapper(m_leafmtxp, m_control,
                                             m_diag_N.data(), inv_chi.data(), m_ndof);
    } else {
        // Fallback to slow method (recomputes all entries in diagonal blocks)
        HACApK_update_diagonal_wrapper(m_leafmtxp, m_control, cHACApK_entry_ij);
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
//
// FIX (2025-02-04): Added interaction pointer tracking to invalidate cache
// when IMA settings change between solves. Previously, stale cache values
// from non-IMA solves were incorrectly used for IMA solves.
//=========================================================================

// Hash-based cache size (must be power of 2 for fast modulo)

//=========================================================================
// GetInteractionMatrixElement: Optimized with O(1) lookup and LRU cache
// Supports the MMM 3DOF tetrahedron HACApK path. Surface-charge MSC uses RadHACApKMomentSystem.
//=========================================================================

double RadHACApKMMMManager::GetInteractionMatrixElement(int dof_i, int dof_j) const {
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

    // EIEM2 retirement (Phase 3b): RadHACApKMMMManager is now MMM-only (tetrahedron, 3 DOF).  MSC
    // surface-charge models (hexahedron / wedge / pyramid) are solved by the multipole-moment MMM H-matrix
    // (RadHACApKMomentSystem) or the dense moment LU -- never this manager -- and mixed MMM+MSC is
    // rejected fail-loud in MakeAutoRelax.  So only the 3x3 (tet-tet) block can occur here.
    if (dof_elem_i == 3 && dof_elem_j == 3) {
        // 3DOF-3DOF: tetra-tetra interaction (IMA-aware via Compute3x3Block_OnDemand / B_comp)
        return GetCached3x3Element(elem_i, elem_j, local_i, local_j);
    }
    std::cerr << "[HACApK] Error: RadHACApKMMMManager received a non-MMM element pair (DOF "
              << dof_elem_i << "/" << dof_elem_j << "); surface-charge MSC is handled by the moment "
              << "solver, not this manager." << std::endl;
    return 0.0;
}

//=========================================================================
// GetCached3x3Element: On-demand 3x3 block computation with thread-local hash cache
// Similar to GetCached6x6Element for 6DOF hexahedra
// Uses O(1) hash lookup instead of O(n) LRU search
//
// FIX (2025-02-04): Added interaction pointer tracking to invalidate cache
// when IMA settings change between solves.
//=========================================================================

// Hash-based cache size for 3DOF (must be power of 2 for fast modulo)
static constexpr int TL_HASH_SIZE_3DOF = 1024;
static constexpr int TL_HASH_MASK_3DOF = TL_HASH_SIZE_3DOF - 1;

double RadHACApKMMMManager::GetCached3x3Element(int elem_i, int elem_j, int comp_i, int comp_j) const {
    // If flat storage is ready (pre-computed), use O(1) direct access
    if (m_flat_N_ready) {
        int64_t base_idx = ((int64_t)elem_i * m_n_elem + elem_j) * 9;
        double val = m_flat_N_data[base_idx + comp_i * 3 + comp_j];
        return val;
    }

    // On-demand computation with thread-local hash cache
    // This path is used when SetupInteractMatrix() is NOT called (HACApK optimization)

    // FIX (2025-02-04): Use generation counter for cache invalidation
    static thread_local uint64_t tl_cached_generation = 0;

    // Thread-local single-entry cache (most common case: same block accessed multiple times)
    static thread_local int tl_single_elem_i = -1;
    static thread_local int tl_single_elem_j = -1;
    static thread_local double tl_single_N_mat[9];

    // Thread-local hash cache (O(1) lookup, no locking needed!)
    static thread_local int tl_cache_elem_i[TL_HASH_SIZE_3DOF];
    static thread_local int tl_cache_elem_j[TL_HASH_SIZE_3DOF];
    static thread_local double tl_cache_N_mat[TL_HASH_SIZE_3DOF][9];
    static thread_local bool tl_initialized = false;

    // Check if generation changed (new solve with different settings)
    uint64_t current_gen = RadHACApKCallback::GetGeneration();
    if (tl_cached_generation != current_gen) {
        // Invalidate all caches
        tl_single_elem_i = -1;
        tl_single_elem_j = -1;
        for (int i = 0; i < TL_HASH_SIZE_3DOF; i++) {
            tl_cache_elem_i[i] = -1;
            tl_cache_elem_j[i] = -1;
        }
        tl_cached_generation = current_gen;
        tl_initialized = true;
    }

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
    // Both paths now return +N (physical quantity).
    // ComputeEntry() handles the sign flip to -N for the system matrix.

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
// GetCached5x5Element: On-demand 5x5 block for wedges with hash cache
// Same pattern as GetCached6x6Element / GetCached3x3Element
// Delegates to radTInteraction::Compute5x5BlockFast (IMA-aware)
//=========================================================================

static constexpr int TL_HASH_SIZE_5DOF = 512;
static constexpr int TL_HASH_MASK_5DOF = TL_HASH_SIZE_5DOF - 1;

//=========================================================================
// GetCachedMixedElement: On-demand mixed-DOF block with hash cache
// Delegates to radTInteraction::ComputeMixedBlockFast (IMA-aware)
// Handles all cross-DOF pairs: 3x5, 3x6, 5x3, 5x6, 6x3, 6x5
//=========================================================================

static constexpr int TL_HASH_SIZE_MIXED = 256;
static constexpr int TL_HASH_MASK_MIXED = TL_HASH_SIZE_MIXED - 1;
static constexpr int MAX_BLOCK_SIZE = 36;  // max DOF product: 6x6

//=========================================================================
// Compute3x3Block: Compute 3x3 interaction block for tetrahedra
// N_ij = interaction matrix element between magnetization components
// Uses existing radTInteraction::InteractMatrix
//=========================================================================

void RadHACApKMMMManager::Compute3x3Block(int elem_i, int elem_j, double* N_mat) const {
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
    //   A = -N - diag(1/chi) (ELF-compatible)
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

void RadHACApKMMMManager::Compute3x3Block_OnDemand(int elem_i, int elem_j, double* N_mat) const {
    // Compute interaction from element j to observation at element i center
    // using B_comp() directly (same approach as SetupInteractMatrix)
    //
    // Returns +N (physical demagnetization tensor).
    // ComputeEntry() handles the sign flip to -N for the system matrix.

    std::memset(N_mat, 0, 9 * sizeof(double));

    if (!m_interaction || elem_i < 0 || elem_i >= m_n_elem ||
        elem_j < 0 || elem_j >= m_n_elem) {
        return;
    }

    radTg3dRelax* elem_row = m_interaction->g3dRelaxPtrVect[elem_i];
    radTg3dRelax* elem_col = m_interaction->g3dRelaxPtrVect[elem_j];
    if (!elem_row || !elem_col) return;

    TVector3d ObsPoiVect = elem_row->ReturnCentrPoint();

    radTFieldKey FieldKeyInteract;
    FieldKeyInteract.B_ = FieldKeyInteract.H_ = FieldKeyInteract.PreRelax_ = 1;

    TVector3d ZeroVect(0., 0., 0.);
    radTField Field(FieldKeyInteract, m_interaction->CompCriterium, ObsPoiVect,
                   ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
    Field.AmOfIntrctElemWithSym = m_interaction->CountRelaxElemsWithSym();

    elem_col->B_comp(&Field);

    // PreRelax mode: Field.B = dH/dMx, Field.H = dH/dMy, Field.A = dH/dMz
    // Return +N (physical): N_mat[row][col] = dH_row / dM_col
    N_mat[0] = Field.B.x;  N_mat[1] = Field.H.x;  N_mat[2] = Field.A.x;
    N_mat[3] = Field.B.y;  N_mat[4] = Field.H.y;  N_mat[5] = Field.A.y;
    N_mat[6] = Field.B.z;  N_mat[7] = Field.H.z;  N_mat[8] = Field.A.z;

    // IMA mirror contributions: B_comp above computes ONLY the direct interaction.
    // For IMA, add mirror contributions using the full IMA-aware kernel, then subtract
    // the direct part (which B_comp already computed correctly).
    if (m_interaction->IsIMAEnabled() && m_interaction->IsTetraGeomReady()) {
        double N_full[9];
        m_interaction->Compute3x3BlockFast(elem_i, elem_j, N_full);
        // N_full = direct + IMA mirrors (from radTInteraction kernel)
        // N_mat = direct (from B_comp, numerically exact for LU compatibility)
        // We want: N_mat_final = N_mat(direct) + (N_full - N_direct_kernel)
        // But N_full already contains direct+mirror, and N_mat contains direct.
        // Since both compute the same physical quantity (direct N), their difference
        // is numerical noise. So N_full - N_mat ≈ pure mirror contribution.
        // Add the mirror contribution: N_final = N_mat + (N_full - N_mat) = N_full
        // This is equivalent to just using N_full, but let's use the IMA mirror
        // extraction approach for clarity and to preserve B_comp's direct accuracy.
        //
        // Actually, the simplest correct approach: use N_full directly when IMA is on.
        // Both Compute3x3BlockFast and B_comp compute the same physics; the kernel
        // version additionally includes IMA mirrors. Any numerical difference in the
        // direct part is within discretization error.
        std::memcpy(N_mat, N_full, 9 * sizeof(double));
    }
}

//=========================================================================
// Compute3x3BlockFast: ELF-style fast 3x3 block computation for tetrahedra
// Uses pre-computed geometry arrays (no object access, no B_comp overhead)
//
// OPTIMIZATION (2025-12-31): Rewritten to compute face basis ONCE per face
// and reuse it for all 3 magnetization directions. This reduces coordinate
// transformation overhead from 12x to 4x per element pair.
//
// The 3x3 interaction matrix N[i][j] represents how magnetization M_j creates
// demagnetizing field H at element i's center:
//   H(r_i) = sum_faces [ sigma_f * integral_triangle H_field dA ] / (4*pi)
// where sigma_f = M_j dot n_f (surface charge density)
//
// For PreRelax mode, we compute dH/dM for each unit M direction.
//=========================================================================

void RadHACApKMMMManager::Compute3x3BlockFast(int elem_i, int elem_j, double* N_mat) const {
    std::memset(N_mat, 0, 9 * sizeof(double));

    if (!m_geometry_3dof_ready || elem_i < 0 || elem_i >= m_n_elem ||
        elem_j < 0 || elem_j >= m_n_elem) {
        return;
    }

    // Observation point: center of element i
    const double* obs_ptr = &m_tetra_centers[elem_i * 3];
    TVector3d obsPoint(obs_ptr[0], obs_ptr[1], obs_ptr[2]);

    // Column element centroid (for normal orientation check)
    const double* col_center = &m_tetra_centers[elem_j * 3];
    TVector3d elemCentroid(col_center[0], col_center[1], col_center[2]);

    // Unit magnetization vectors
    const TVector3d unit_M[3] = {
        TVector3d(1.0, 0.0, 0.0),
        TVector3d(0.0, 1.0, 0.0),
        TVector3d(0.0, 0.0, 1.0)
    };

    // Process each of the 4 triangular faces
    // ELF-style optimization: compute face basis ONCE per face, reuse for all 3 M directions
    for (int f = 0; f < 4; f++) {
        // Get face vertices from pre-computed geometry
        int fvIdx = (elem_j * 4 + f) * 3 * 3;
        const double* V0_ptr = &m_tetra_face_vertices[fvIdx + 0];
        const double* V1_ptr = &m_tetra_face_vertices[fvIdx + 3];
        const double* V2_ptr = &m_tetra_face_vertices[fvIdx + 6];

        TVector3d V0(V0_ptr[0], V0_ptr[1], V0_ptr[2]);
        TVector3d V1(V1_ptr[0], V1_ptr[1], V1_ptr[2]);
        TVector3d V2(V2_ptr[0], V2_ptr[1], V2_ptr[2]);

        // Compute face basis ONCE (ELF optimization - saves 2 basis computations per face)
        RadTriangleFaceBasis basis;
        RadComputeTriangleFaceBasis(V0, V1, V2, elemCentroid, basis);

        if (!basis.valid) continue;

        // Get face normal for surface charge computation
        // Use pre-computed values from m_tetra_face_normals
        int fnIdx = (elem_j * 4 + f) * 3;
        TVector3d faceNormal(
            m_tetra_face_normals[fnIdx + 0],
            m_tetra_face_normals[fnIdx + 1],
            m_tetra_face_normals[fnIdx + 2]
        );

        // Compute H field for each unit magnetization direction
        // Reuse the same face basis for all 3 directions
        for (int j = 0; j < 3; j++) {
            // Surface charge: sigma = M dot n * sign_correction
            double sigma = (unit_M[j].x * faceNormal.x +
                           unit_M[j].y * faceNormal.y +
                           unit_M[j].z * faceNormal.z) * basis.charge_sign;

            if (std::abs(sigma) < 1.0e-15) continue;

            // Compute field using pre-computed basis (fast path)
            TVector3d H = RadFieldFromTriangleFaceWithBasis(basis, sigma, obsPoint);

            // Accumulate into tensor: N_mat[i*3+j] = dH_i / dM_j (+N convention)
            // ComputeEntry() handles the sign flip to -N for the system matrix.
            N_mat[0*3 + j] += H.x;
            N_mat[1*3 + j] += H.y;
            N_mat[2*3 + j] += H.z;
        }
    }
}

double RadHACApKMMMManager::GetGenericElement(int elem_i, int elem_j, int local_i, int local_j) const {
    // Generic path: access pre-computed flat interaction matrix (+N convention)
    if (!m_interaction || m_interaction->m_flatInteractMatrix.empty()) return 0.0;

    int offset_i = m_dof_offset[elem_i];
    int offset_j = m_dof_offset[elem_j];
    int total_dof = m_interaction->m_totalDOF;

    return m_interaction->m_flatInteractMatrix[(offset_i + local_i) * total_dof + (offset_j + local_j)];
}

//=========================================================================
// RadHACApKMomentSystem: the multipole-moment MMM system A_raw as a HACApK H-matrix
// (docs/multipole_moment_mmm/ACA_MOMENT_DESIGN.md, Phase 2 Increment 2).
//=========================================================================

RadHACApKMomentSystem::RadHACApKMomentSystem(radTInteraction* interaction, double chi)
    : m_interaction(interaction), m_chi(chi)
{
}

// Per-element chi (Increment 4, nonlinear Picard): each row 6h+* folds chiPerHex[h] (the row element's
// susceptibility) into A_raw via MomentSystemEntry; resolved in ExtractCoordinates once nHex is known.
RadHACApKMomentSystem::RadHACApKMomentSystem(radTInteraction* interaction, const std::vector<double>& chiPerHex)
    : m_interaction(interaction), m_chi(chiPerHex.empty() ? 1.0 : chiPerHex[0]), m_chi_in(chiPerHex)
{
}

void RadHACApKMomentSystem::ExtractCoordinates()
{
    if (!m_interaction) { m_ndof = 0; m_n_elem = 0; return; }
    const std::vector<int>& hexElem = m_interaction->GetHexaElemIndices();
    int nHex = (int)hexElem.size();
    m_n_elem = nHex;
    m_ndof = 6 * nHex;
    m_coordinates.resize((size_t)nHex * 3);
    m_dof_offset.resize(nHex + 1);
    for (int h = 0; h < nHex; h++) {
        TVector3d c = m_interaction->GetElementCenter(hexElem[h]);
        m_coordinates[(size_t)h * 3 + 0] = c.x;
        m_coordinates[(size_t)h * 3 + 1] = c.y;
        m_coordinates[(size_t)h * 3 + 2] = c.z;
        m_dof_offset[h] = 6 * h;
    }
    m_dof_offset[nHex] = 6 * nHex;
    // chi per hex for MomentSystemEntry: per-element (Increment 4, nonlinear Picard) when the vector ctor
    // supplied one of matching length, else uniform m_chi.
    if ((int)m_chi_in.size() == nHex && nHex > 0) m_chiv = m_chi_in;
    else m_chiv.assign((size_t)(nHex > 0 ? nHex : 1), m_chi);
}

void RadHACApKMomentSystem::OnBeforeBuild()
{
    if (!m_interaction) return;
    RadHACApKCallback::SetInteraction(m_interaction, m_n_elem, 6);
    m_interaction->PrecomputeMomentGeometry();
}

double RadHACApKMomentSystem::GetInteractionMatrixElement(int dof_i, int dof_j) const
{
    if (!m_interaction || m_chiv.empty()) return 0.0;
    if (dof_i < 0 || dof_j < 0 || dof_i >= m_ndof || dof_j >= m_ndof) return 0.0;
    int elem_i = dof_i / 6, elem_j = dof_j / 6;
    int local_i = dof_i - 6 * elem_i, local_j = dof_j - 6 * elem_j;

    static constexpr int TL_HASH_SIZE_MOMENT6 = 1024;
    static constexpr int TL_HASH_MASK_MOMENT6 = TL_HASH_SIZE_MOMENT6 - 1;
    static thread_local uint64_t tl_cached_generation = 0;
    static thread_local int tl_single_elem_i = -1;
    static thread_local int tl_single_elem_j = -1;
    static thread_local double tl_single_block[36];
    static thread_local int tl_cache_elem_i[TL_HASH_SIZE_MOMENT6];
    static thread_local int tl_cache_elem_j[TL_HASH_SIZE_MOMENT6];
    static thread_local double tl_cache_block[TL_HASH_SIZE_MOMENT6][36];
    static thread_local bool tl_initialized = false;

    uint64_t current_gen = RadHACApKCallback::GetGeneration();
    if (tl_cached_generation != current_gen || !tl_initialized) {
        tl_single_elem_i = -1;
        tl_single_elem_j = -1;
        for (int i = 0; i < TL_HASH_SIZE_MOMENT6; i++) {
            tl_cache_elem_i[i] = -1;
            tl_cache_elem_j[i] = -1;
        }
        tl_cached_generation = current_gen;
        tl_initialized = true;
    }

    if (tl_single_elem_i == elem_i && tl_single_elem_j == elem_j) {
        return tl_single_block[local_i * 6 + local_j];
    }

    unsigned int hash_idx = (((unsigned int)elem_i * 73856093u) ^ ((unsigned int)elem_j * 19349663u)) & TL_HASH_MASK_MOMENT6;
    if (tl_cache_elem_i[hash_idx] == elem_i && tl_cache_elem_j[hash_idx] == elem_j) {
        std::memcpy(tl_single_block, tl_cache_block[hash_idx], 36 * sizeof(double));
        tl_single_elem_i = elem_i;
        tl_single_elem_j = elem_j;
        return tl_single_block[local_i * 6 + local_j];
    }

    m_interaction->MomentSystemBlock6x6(elem_i, elem_j, m_chiv.data(), tl_single_block);
    tl_single_elem_i = elem_i;
    tl_single_elem_j = elem_j;
    tl_cache_elem_i[hash_idx] = elem_i;
    tl_cache_elem_j[hash_idx] = elem_j;
    std::memcpy(tl_cache_block[hash_idx], tl_single_block, 36 * sizeof(double));
    return tl_single_block[local_i * 6 + local_j];
}

//=========================================================================
// End of rad_hacapk.cpp
//=========================================================================
