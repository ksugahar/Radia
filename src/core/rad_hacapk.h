/*-------------------------------------------------------------------------
*
* File name:      rad_hacapk.h
*
* Project:        RADIA
*
* Description:    HACApK (H-matrix with ACA+) interface for BiCGSTAB solver
*                 Provides O(N log N) matrix-vector products for large problems
*
*                 RadHACApKBase owns the kernel-agnostic H-matrix lifecycle used
*                 by the HDiv, PEEC, and BEM managers.
*
* First release:  2025
*
* Reference:      ppOpen-HPC project (MIT License)
*                 https://github.com/ppohHPC/ppOpen-HPC
*
-------------------------------------------------------------------------*/

#ifndef __RAD_HACAPK_H
#define __RAD_HACAPK_H

#include <vector>
#include <mutex>

namespace RadHACApKCallback {
std::mutex& OperationMutex();
void ClearCallbackException();
void CaptureCallbackException() noexcept;
void RethrowCallbackException();
}

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
        print_level(0)   // 0=silent (set to 1 for standard output)
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
// RadHACApKBase: kernel-agnostic H-matrix manager
//-------------------------------------------------------------------------

/**
 * RadHACApKBase holds the HACApK build / matvec / diagonal-update lifecycle
 * that is independent of the physical kernel. Subclasses supply the
 * integration kernel via four virtual hooks:
 *
 *   - ExtractCoordinates(): populate m_coordinates, m_n_elem, m_ndof,
 *     m_dof_offset (and any kernel-specific lookup tables)
 *   - OnBeforeBuild():      run kernel-specific precomputation before
 *     HACApK calls back into ComputeEntry
 *   - InitializeInvChi():   populate m_inv_chi (size == m_ndof) from the
 *     current material / constitutive state
 *   - GetInteractionMatrixElement(i, j): return the +N(i, j) physical
 *     matrix element (system matrix A = -N + diag(1/chi); the sign flip
 *     is applied once inside RadHACApKCallback::ComputeEntry)
 *
 * Concrete managers provide HDiv charge-Gram, PEEC, or BEM kernels and use
 * this class for build, matrix-vector, and diagonal-update operations.
 */
class RadHACApKBase {
public:
    virtual ~RadHACApKBase();

    /**
     * Build H-matrix. Calls ExtractCoordinates, OnBeforeBuild,
     * InitializeInvChi, then hands off to HACApK.
     */
    bool BuildHMatrix(const RadHACApKParams& params = RadHACApKParams());

    /**
     * Matrix-vector product: y = A * x (O(N log N))
     *
     * VIRTUAL (2026-08-09): subclasses may store a diagonally NORMALIZED
     * operator in the H-matrix leaves and wrap the scaling back here, so
     * base-internal applications must dispatch (see
     * RadHACApKChargeGram::MatVecSym -- the charge-basis normalization that
     * fixes the roundoff-amplified indefiniteness on extreme-size-ratio
     * meshes).  Implementations are unchanged for every other subclass.
     */
    virtual void MatVec(const std::vector<double>& x, std::vector<double>& y);
    // y = A^T x (transpose H-matvec) and y = G_sym x (EXACTLY symmetric apply built from the
    // upper-triangular leaves -- valid for a symmetric cluster tree like the charge Gram).
    virtual void MatVecTranspose(const std::vector<double>& x, std::vector<double>& y);
    virtual void MatVecSym(const std::vector<double>& x, std::vector<double>& y);
    // Row-major [nrhs][ndof] BLAS-3 symmetric H-matrix application.
    virtual void MatVecSymMany(const std::vector<double>& x, int nrhs,
                               std::vector<double>& y);
    void MatVecSymManyMasked(const std::vector<double>& x, int nrhs,
                             const std::vector<int>& active_prefix,
                             std::vector<double>& y);

    /**
     * Update diagonal blocks when 1/chi changes (nonlinear iteration)
     */
    void UpdateDiagonal(const std::vector<double>& inv_chi);

    /**
     * Set the local null-space (loop) deflation basis L for a matrix-free
     * shift: MatVec then returns (A + alpha * L L^T) x.  L is given as a
     * sparse list of plaquettes/cycles (each a few (dof, sign) entries),
     * CSR-like via plaq_offsets.  L L^T acts only on span(L) = the null
     * (loop) subspace, lifting its eigenvalues away from zero (O(nnz)=O(N)).
     * Pass empty plaq_offsets / alpha=0 to disable.  L must be in the
     * ORIGINAL DOF ordering (same ordering as MatVec's x/y).
     */
    void SetDeflation(const std::vector<int>& plaq_offsets,
                      const std::vector<int>& dofs,
                      const std::vector<double>& signs,
                      double alpha);

    /** Number of loop (plaquette + belt) cycles in the installed deflation
     *  basis L (0 if deflation is off). Diagnostic for the belted-tree count. */
    int GetDeflationCycleCount() const { return m_defl_nplaq; }

    /** The shift strength alpha actually used (after auto-scaling when the
     *  caller requested alpha <= 0). Diagnostic. */
    double GetDeflationAlpha() const { return m_defl_alpha; }

    /**
     * Return the kernel's physical +N(i, j) interaction matrix element.
     * For magnetic interaction kernels this is the demagnetization tensor contribution (system
     * matrix is -N + diag(1/chi)); for PEEC this is the +L mutual
     * inductance (system matrix is L itself, frequency-dependent
     * factors applied outside HACApK).
     *
     * The physical vs system-matrix convention is bridged by
     * ComputeSystemEntry below.
     */
    virtual double GetInteractionMatrixElement(int dof_i, int dof_j) const = 0;

    /**
     * Return the system-matrix entry A(i, j) as stored by HACApK.
     * Default implementation returns GetInteractionMatrixElement
     * unchanged (PEEC / BEM convention where the physical N is itself
     * the system matrix). Magnetostatic kernels override this to apply the sign flip
     * and the diag(1/chi) shift: A = -N + delta_ij / chi_i.
     *
     * This is the hook called by RadHACApKCallback::ComputeEntry, so
     * each kernel decides exactly what HACApK stores.
     */
    virtual double ComputeSystemEntry(int dof_i, int dof_j) const {
        return GetInteractionMatrixElement(dof_i, dof_j);
    }

    // Accessors
    bool IsValid() const { return m_valid; }
    int GetNDOF() const { return m_ndof; }
    const std::vector<int>& GetPermutation() const { return m_permutation; }
    const RadHACApKStats& GetStats() const { return m_stats; }
    const std::vector<double>& GetDiagonalN() const { return m_diag_N; }
    bool IsDiagonalCached() const { return m_diag_cached; }

    // Phase 4: opaque pointer accessors so external H-LU drivers
    // (cHACApK_hlu_run_on_hacapk) can reach the HACApK leafmtxp +
    // control without breaking encapsulation. Read-only; callers must
    // not free these.
    void* GetLeafmtxp() const { return m_leafmtxp; }
    void* GetLcontrol() const { return m_control; }

protected:
    RadHACApKBase();

    // Internal row-major symmetric apply with optional active principal mask
    // and fused S*A*S diagonal scaling.  The optional arrays contain ndof+1
    // and ndof entries and remain caller-owned for the duration of the call.
    void MatVecSymManyPrepared(
        const std::vector<double>& x, int nrhs,
        const int* active_prefix, const double* diagonal_scale,
        std::vector<double>& y);

    // Kernel hooks
    virtual void OnBuildStarting(const RadHACApKParams&) {}
    virtual void ExtractCoordinates() = 0;
    virtual void OnBeforeBuild() = 0;
    virtual void InitializeInvChi() = 0;
    virtual bool UseSymmetricFill() const { return false; }
    virtual void OnBuildFinished(bool) noexcept {}

    /**
     * Return true to use variable-DOF HACApK build (per-element DOF via
     * m_dof_offset). Return false for uniform DOF (m_nffc per element).
     */
    virtual bool IsVariableDOF() const = 0;

    /**
     * Uniform DOF count per element (used when IsVariableDOF() is false).
     * Typical: 3 (compact vector element), 1 (PEEC filament / face DOF).
     */
    virtual int GetUniformNFFC() const = 0;

    void FreeResources();

    // HACApK opaque structures (owned)
    void* m_leafmtxp;    // st_cHACApK_leafmtxp*
    void* m_control;     // st_cHACApK_lcontrol*

    // Global state
    bool m_valid;
    int m_ndof;
    int m_n_elem;

    // DOF permutation returned by HACApK clustering
    std::vector<int> m_permutation;

    // Per-element DOF offset (size = n_elem + 1); required for variable-DOF build
    std::vector<int> m_dof_offset;

    // Element center coordinates for clustering, size = n_elem * 3
    std::vector<double> m_coordinates;

    // Current inverse susceptibility (1/chi per DOF), size = n_dof
    std::vector<double> m_inv_chi;

    // Cached diagonal of N (physical, NOT system matrix), size = n_dof
    std::vector<double> m_diag_N;
    bool m_diag_cached;

    // Optional matrix-free loop-mode deflation: MatVec adds alpha * L (L^T x)
    // with L a sparse loop (plaquette/cycle) basis stored CSR-like.
    std::vector<int> m_defl_offsets;   // size n_plaq + 1 (CSR row pointers)
    std::vector<int> m_defl_dofs;      // flat DOF indices
    std::vector<double> m_defl_signs;  // flat +/-1 entries
    double m_defl_alpha;
    int m_defl_nplaq;

    RadHACApKStats m_stats;

private:
    RadHACApKBase(const RadHACApKBase&) = delete;
    RadHACApKBase& operator=(const RadHACApKBase&) = delete;
};

//-------------------------------------------------------------------------
// Global callback state for HACApK
// (Required because HACApK C interface uses global callback function)
//-------------------------------------------------------------------------

namespace RadHACApKCallback {
    // Set the current manager for callbacks (base pointer; kernel resolved via virtual dispatch)
    void SetCurrentManager(RadHACApKBase* manager);

    // Get the current manager
    RadHACApKBase* GetCurrentManager();

    // Clear all global callback state (called on manager destruction)
    void ClearGlobalState();

    // Compute the active manager's system-matrix entry. Called from
    // cHACApK_entry_ij (1-based indexing in HACApK's convention).
    double ComputeEntry(int i, int j);
}

//-------------------------------------------------------------------------

#endif // __RAD_HACAPK_H
