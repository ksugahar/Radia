/*-------------------------------------------------------------------------
*
* File name:      rad_mmm_hacapk.h
*
* Project:        RADIA - mmm_core module
*
* Description:    HACApK H-matrix solver for standalone MMM
*                 Provides O(N log N) matrix-vector products for large problems
*
* First release:  2026
*
* Based on:       ppOpen-HPC HACApK (MIT License)
*
-------------------------------------------------------------------------*/

#ifndef RAD_MMM_HACAPK_H
#define RAD_MMM_HACAPK_H

#include "rad_mmm_matrices.h"
#include <vector>
#include <memory>

namespace radia {
namespace mmm {

//-------------------------------------------------------------------------
// HACApK Configuration Parameters
//-------------------------------------------------------------------------

struct HACApKParams {
    double aca_eps;      // ACA+ compression tolerance (default: 1e-4)
    int leaf_size;       // Minimum cluster size (default: 32)
    double eta;          // Admissibility parameter (default: 2.0)
    int max_rank;        // Maximum ACA rank (default: 200)
    int print_level;     // 0=silent, 1=standard, 2=debug

    HACApKParams() :
        aca_eps(1.0e-4),
        leaf_size(32),
        eta(2.0),
        max_rank(200),
        print_level(0)
    {}
};

//-------------------------------------------------------------------------
// HACApK Statistics
//-------------------------------------------------------------------------

struct HACApKStats {
    int n_lowrank;       // Number of low-rank (ACA compressed) blocks
    int n_dense;         // Number of dense blocks
    int max_rank;        // Maximum rank among low-rank blocks
    int n_leaves;        // Total number of leaf blocks
    int n_dof;           // Total degrees of freedom
    double compression;  // Compression ratio (H-matrix memory / dense memory)
    double build_time;   // Time to build H-matrix (seconds)
    double memory_mb;    // Actual H-matrix memory usage [MB]
    double dense_memory_mb;  // Full dense matrix memory [MB]

    HACApKStats() :
        n_lowrank(0), n_dense(0), max_rank(0), n_leaves(0),
        n_dof(0), compression(1.0), build_time(0.0),
        memory_mb(0.0), dense_memory_mb(0.0)
    {}
};

//-------------------------------------------------------------------------
// MMM HACApK Solver
//-------------------------------------------------------------------------

/**
 * MMMHACApKSolver - H-matrix accelerated BiCGSTAB solver for MMM
 *
 * This class provides H-matrix acceleration for the MMM linear system:
 *   A * M = H_ext
 * where A = (1/chi) * I - N is the system matrix.
 *
 * The H-matrix approximates A using ACA+ compression, providing
 * O(N log N) complexity for matrix-vector products instead of O(N^2).
 *
 * Usage:
 *   MMMHACApKSolver solver;
 *   solver.SetElements(elements, dof_offset);
 *   solver.SetInteractionMatrix(N);
 *   solver.BuildHMatrix(params);
 *   solver.Solve(inv_chi, H_ext, M, tol, max_iter);
 */
class MMMHACApKSolver {
public:
    MMMHACApKSolver();
    ~MMMHACApKSolver();

    /**
     * Set element data from MMMBuilder
     * @param elements Element data (centers, volumes, DOF)
     * @param dof_offset DOF offset array
     */
    void SetElements(const std::vector<MMMElement>& elements,
                     const std::vector<int>& dof_offset);

    /**
     * Set interaction matrix N (from MMMBuilder::Build())
     * @param N Interaction matrix (total_dof x total_dof)
     * @param total_dof Total degrees of freedom
     */
    void SetInteractionMatrix(const std::vector<double>& N, int total_dof);

    /**
     * Build H-matrix approximation
     * @param params H-matrix parameters
     * @param inv_chi Inverse susceptibility per element (n_elements)
     *                If empty, uses previously set value or default (0 = infinite chi)
     * @return true on success
     */
    bool BuildHMatrix(const HACApKParams& params = HACApKParams(),
                      const std::vector<double>& inv_chi = std::vector<double>());

    /**
     * Solve linear system using H-matrix accelerated BiCGSTAB
     * @param inv_chi Inverse susceptibility per element (n_elements)
     * @param H_ext External field (total_dof)
     * @param M Output magnetization (total_dof)
     * @param tol Convergence tolerance
     * @param max_iter Maximum iterations
     * @return Number of iterations (negative if not converged)
     */
    int Solve(const std::vector<double>& inv_chi,
              const std::vector<double>& H_ext,
              std::vector<double>& M,
              double tol = 1e-8,
              int max_iter = 1000);

    /**
     * Matrix-vector product y = A * x using H-matrix
     * @param x Input vector
     * @param y Output vector
     */
    void MatVec(const std::vector<double>& x, std::vector<double>& y);

    /**
     * Update diagonal when chi changes (for nonlinear iteration)
     * @param inv_chi New inverse susceptibility values
     */
    void UpdateDiagonal(const std::vector<double>& inv_chi);

    /**
     * Get H-matrix statistics
     */
    const HACApKStats& GetStats() const { return stats_; }

    /**
     * Check if H-matrix is built and valid
     */
    bool IsValid() const { return valid_; }

    /**
     * Get number of DOF
     */
    int GetNDOF() const { return total_dof_; }

    /**
     * Get interaction matrix element N(i,j)
     * Used by HACApK callback during build
     */
    double GetNElement(int i, int j) const;

    /**
     * Get current inverse chi for diagonal computation
     */
    const std::vector<double>& GetInvChi() const { return inv_chi_; }

    /**
     * Get DOF offset for element lookup
     */
    const std::vector<int>& GetDOFOffset() const { return dof_offset_; }

private:
    // Element data
    std::vector<MMMElement> elements_;
    std::vector<int> dof_offset_;
    int n_elements_;
    int total_dof_;

    // Interaction matrix (dense, for callback during H-matrix build)
    std::vector<double> N_;

    // Current inverse susceptibility
    std::vector<double> inv_chi_;

    // Element coordinates for clustering
    std::vector<double> coordinates_;

    // HACApK internal structures (opaque pointers via cHACApK_cpp.h)
    void* leafmtxp_;   // st_cHACApK_leafmtxp* - allocated via HACApK_alloc_leafmtxp()
    void* control_;    // st_cHACApK_lcontrol* - allocated via HACApK_alloc_lcontrol()

    // DOF permutation from HACApK clustering
    std::vector<int> permutation_;

    // State
    bool valid_;
    HACApKStats stats_;

    // Private methods
    void FreeHACApK();
    void ExtractCoordinates();
    void BuildDOFToElementMap();

    // DOF to element lookup
    std::vector<int> dof_to_elem_;
    std::vector<int> dof_to_local_;

    // Disable copy
    MMMHACApKSolver(const MMMHACApKSolver&) = delete;
    MMMHACApKSolver& operator=(const MMMHACApKSolver&) = delete;
};

//-------------------------------------------------------------------------
// Global callback state for HACApK
//-------------------------------------------------------------------------

namespace MMMHACApKCallback {
    void SetCurrentSolver(MMMHACApKSolver* solver);
    MMMHACApKSolver* GetCurrentSolver();
    double ComputeEntry(int i, int j);
    void SetLod(int* lod, int size);
    void ClearLod();
}

}  // namespace mmm
}  // namespace radia

#endif // RAD_MMM_HACAPK_H
