/*-------------------------------------------------------------------------
*
* File name:      rad_hacapk.h
*
* Project:        RADIA
*
* Description:    HACApK (H-matrix with ACA+) interface for BiCGSTAB solver
*                 Provides O(N log N) matrix-vector products for large problems
*
* First release:  2025
*
* Reference:      ppOpen-HPC project (MIT License)
*                 https://github.com/ppohHPC/ppOpen-HPC
*
-------------------------------------------------------------------------*/

#ifndef __RAD_HACAPK_H
#define __RAD_HACAPK_H

#include "rad_interaction.h"
#include "rad_polyhedron.h"  // For radTPolyhedron 6DOF MSC support
#include <vector>
#include <functional>
#include <unordered_map>

// HACApK uses opaque void* pointers to avoid C/C++ struct compatibility issues
// Actual structures are defined in cHACApK_cpp.h and managed by cHACApK_cpp_impl.c

//-------------------------------------------------------------------------
// HACApK Configuration Parameters
//-------------------------------------------------------------------------

struct RadHACApKParams {
    double aca_eps;      // ACA+ compression tolerance (default: 1e-4)
    int leaf_size;       // Minimum cluster size (default: 32)
    double eta;          // Admissibility parameter (default: 2.0)
    int max_rank;        // Maximum ACA rank (default: 200)
    int print_level;     // 0=silent, 1=standard, 2=debug

    RadHACApKParams() :
        aca_eps(1.0e-4),
        leaf_size(32),   // Default: 32 (element count, not DOF)
        eta(2.0),
        max_rank(200),
        print_level(1)   // 1=standard output for debugging
    {}
};

//-------------------------------------------------------------------------
// HACApK Statistics
//-------------------------------------------------------------------------

struct RadHACApKStats {
    int n_lowrank;       // Number of low-rank (ACA compressed) blocks
    int n_dense;         // Number of dense blocks
    int max_rank;        // Maximum rank among low-rank blocks
    int n_leaves;        // Total number of leaf blocks
    int n_dof;           // Total degrees of freedom
    double compression;  // Compression ratio (H-matrix memory / dense memory)
    double build_time;   // Time to build H-matrix (seconds)
    double memory_mb;    // Actual H-matrix memory usage [MB]
    double dense_memory_mb;  // Full dense matrix memory [MB]

    RadHACApKStats() :
        n_lowrank(0), n_dense(0), max_rank(0), n_leaves(0),
        n_dof(0), compression(1.0), build_time(0.0),
        memory_mb(0.0), dense_memory_mb(0.0)
    {}
};

//-------------------------------------------------------------------------
// HACApK Manager for Radia
//-------------------------------------------------------------------------

/**
 * RadHACApKManager handles the lifecycle of H-matrix for BiCGSTAB solver
 *
 * Usage:
 *   1. Create manager with RadHACApKManager(interaction_ptr)
 *   2. Call BuildHMatrix() to construct the H-matrix
 *   3. Call MatVec(x, y) for matrix-vector products
 *   4. Call UpdateDiagonal() when chi(H) changes (nonlinear iteration)
 *
 * The H-matrix approximates A = -N + diag(1/chi) using ACA+ compression.
 * For problems where elements are spatially well-separated, this provides
 * O(N log N) complexity instead of O(N^2) for dense matrix-vector products.
 *
 * Note: For single compact objects (typical Radia use case), the admissibility
 * criterion may not be satisfied and most blocks remain dense. In such cases,
 * the dense BiCGSTAB solver (Method 1) is more efficient.
 */
class RadHACApKManager {
public:
    /**
     * Constructor
     * @param interaction Pointer to Radia interaction object (must remain valid)
     */
    RadHACApKManager(radTInteraction* interaction);

    /**
     * Destructor - frees HACApK resources
     */
    ~RadHACApKManager();

    /**
     * Build H-matrix from current interaction matrix
     * @param params Configuration parameters
     * @return true on success, false on failure
     */
    bool BuildHMatrix(const RadHACApKParams& params = RadHACApKParams());

    /**
     * Matrix-vector product: y = A * x
     * Uses H-matrix approximation for O(N log N) complexity
     * @param x Input vector (size = n_dof)
     * @param y Output vector (size = n_dof)
     */
    void MatVec(const std::vector<double>& x, std::vector<double>& y);

    /**
     * Update diagonal blocks when chi(H) changes
     * Called during nonlinear iteration when material susceptibility updates
     * @param inv_chi New inverse susceptibility values (1/chi for each DOF)
     */
    void UpdateDiagonal(const std::vector<double>& inv_chi);

    /**
     * Check if H-matrix is valid and ready for use
     */
    bool IsValid() const { return m_valid; }

    /**
     * Get the interaction object this H-matrix was built for
     */
    radTInteraction* GetInteraction() const { return m_interaction; }

    /**
     * Get H-matrix statistics
     */
    const RadHACApKStats& GetStats() const { return m_stats; }

    /**
     * Get number of degrees of freedom
     */
    int GetNDOF() const { return m_ndof; }

    /**
     * Get permutation array (for coordinate reordering)
     */
    const std::vector<int>& GetPermutation() const { return m_permutation; }

    /**
     * Get interaction matrix element N(dof_i, dof_j)
     * Uses friend access to radTInteraction::InteractMatrix
     * @param dof_i Row DOF index (0-based)
     * @param dof_j Column DOF index (0-based)
     * @return The N matrix element value
     */
    double GetInteractionMatrixElement(int dof_i, int dof_j) const;

    /**
     * Get cached diagonal elements of interaction matrix N_ii
     * These are computed once during BuildHMatrix for efficiency
     * @return Reference to cached diagonal vector (size = n_dof)
     */
    const std::vector<double>& GetDiagonalN() const { return m_diag_N; }

    /**
     * Check if diagonal elements are cached
     */
    bool IsDiagonalCached() const { return m_diag_cached; }

    /**
     * Check if flat N storage is ready (for 3DOF tetrahedra)
     */
    bool IsFlatNReady() const { return m_flat_N_ready; }

    /**
     * Flatten InteractMatrix for 3DOF tetrahedra (O(1) access)
     * Must be called AFTER InteractMatrix is computed
     * Safe to call multiple times (no-op if already ready)
     */
    void PrecomputeFlatInteractMatrix();

private:
    // Pointer to Radia interaction (not owned)
    radTInteraction* m_interaction;

    // HACApK internal structures (owned, opaque pointers)
    void* m_leafmtxp;   // st_cHACApK_leafmtxp*
    void* m_control;    // st_cHACApK_lcontrol*

    // State
    bool m_valid;
    int m_ndof;
    int m_n_elem;
    int m_nffc;  // DOF per element (3 for tetra, 6 for hexa)
    bool m_is_6dof;  // true if using 6DOF MSC hexahedra

    // DOF permutation for cluster ordering
    std::vector<int> m_permutation;

    // DOF offset per element (for variable DOF support)
    std::vector<int> m_dof_offset;

    // O(1) DOF-to-element lookup table (ELF-style)
    // m_dof_to_elem[dof] = element index
    // m_dof_to_local[dof] = local DOF index within element (0-5 for hexa)
    std::vector<int> m_dof_to_elem;
    std::vector<int> m_dof_to_local;

    // Statistics
    RadHACApKStats m_stats;

    // Current inverse susceptibility values
    std::vector<double> m_inv_chi;

    // Element coordinates for clustering
    std::vector<double> m_coordinates;  // [n_elem * 3]

    // ========================================================================
    // ELF-style pre-computed geometry (for 6DOF hexahedra)
    // ========================================================================
    // Pre-computed geometry avoids dynamic_cast and scattered memory access
    // during matrix element computation (major performance optimization)

    // Pre-computed element centers [n_elem * 3]
    std::vector<double> m_elem_centers;

    // Pre-computed vertices [n_elem * 8 * 3] (8 vertices per hexa, 3 coords each)
    std::vector<double> m_elem_vertices;

    // Pre-computed face centers [n_elem * 6 * 3] (6 faces per hexa)
    std::vector<double> m_face_centers;

    // Pre-computed face normals [n_elem * 6 * 3]
    std::vector<double> m_face_normals;

    // Pre-computed face areas [n_elem * 6]
    std::vector<double> m_face_areas;

    // Pre-computed quad face vertices [n_elem * 6 * 4 * 3] (6 faces, 4 verts, 3 coords)
    std::vector<double> m_face_vertices;

    // Flag: geometry has been pre-computed
    bool m_geometry_ready;

    // Cached diagonal elements of interaction matrix N_ii (for Jacobi preconditioner)
    // Computed once during BuildHMatrix, reused in every BiCGSTAB iteration
    std::vector<double> m_diag_N;
    bool m_diag_cached;

    // ========================================================================
    // Flat interaction matrix storage (for 3DOF tetrahedra)
    // ========================================================================
    // Pre-computed flat array for O(1) access to N_ij elements
    // Layout: m_flat_N[elem_i * n_elem + elem_j] = 3x3 block starting index
    // Actual data: m_flat_N_data[(elem_i * n_elem + elem_j) * 9 + row * 3 + col]
    std::vector<double> m_flat_N_data;  // size = n_elem * n_elem * 9
    bool m_flat_N_ready;

    // Private methods
    void FreeResources();
    void ExtractElementCoordinates();
    void BuildDOFLookupTable();
    void PrecomputeGeometry();  // ELF-style pre-computation for 6DOF hexahedra

    // 6DOF hexahedron methods
    double GetCached6x6Element(int elem_i, int elem_j, int face_i, int face_j) const;
    void Compute6x6Block(int elem_i, int elem_j, double* K_mat) const;
    void Compute6x6BlockFast(int elem_i, int elem_j, double* K_mat) const;  // Uses pre-computed geometry

    // 3DOF tetrahedron methods
    double GetCached3x3Element(int elem_i, int elem_j, int comp_i, int comp_j) const;
    void Compute3x3Block(int elem_i, int elem_j, double* N_mat) const;

    // Field computation from pre-computed face vertices
    void FieldFromQuadFaceFast(int elem, int face, const double* obs, double sigma, double* H_out) const;
    void FieldFromChargedTriangleFast(const double* obs, const double* v0, const double* v1, const double* v2, double sigma, double* H_out) const;

    // Disable copy
    RadHACApKManager(const RadHACApKManager&) = delete;
    RadHACApKManager& operator=(const RadHACApKManager&) = delete;
};

//-------------------------------------------------------------------------
// Global callback state for HACApK
// (Required because HACApK C interface uses global callback function)
//-------------------------------------------------------------------------

namespace RadHACApKCallback {
    // Set the current manager for callbacks
    void SetCurrentManager(RadHACApKManager* manager);

    // Get the current manager
    RadHACApKManager* GetCurrentManager();

    // Set inverse susceptibility for matrix element computation
    void SetInvChi(const std::vector<double>& inv_chi);

    // Get inverse susceptibility
    const std::vector<double>& GetInvChi();

    // Compute matrix element A(i,j) = -N(i,j) + delta_ij/chi_i
    // Called from cHACApK_entry_ij
    double ComputeEntry(int i, int j);
}

//-------------------------------------------------------------------------

#endif // __RAD_HACAPK_H
