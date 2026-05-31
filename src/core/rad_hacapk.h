/*-------------------------------------------------------------------------
*
* File name:      rad_hacapk.h
*
* Project:        RADIA
*
* Description:    HACApK (H-matrix with ACA+) interface for BiCGSTAB solver
*                 Provides O(N log N) matrix-vector products for large problems
*
*                 Refactored 2026-04-16:
*                 - RadHACApKBase owns the kernel-agnostic H-matrix lifecycle
*                 - RadHACApKMSCManager : public RadHACApKBase implements the
*                   MMM/MSC kernel (tetra/wedge/hex magnetization moments and
*                   surface charges). A future RadHACApKPEECManager will
*                   implement Ruehli finite-filament mutual inductance.
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
 * Typical usage:
 *   RadHACApKMSCManager mgr(interaction);
 *   mgr.BuildHMatrix();
 *   mgr.MatVec(x, y);                    // y = A * x
 *   mgr.UpdateDiagonal(new_inv_chi);     // nonlinear iteration update
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
     */
    void MatVec(const std::vector<double>& x, std::vector<double>& y);

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
     * For MSC this is the demagnetization tensor contribution (system
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
     * the system matrix). MSC overrides this to apply the sign flip
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

    // Kernel hooks
    virtual void ExtractCoordinates() = 0;
    virtual void OnBeforeBuild() = 0;
    virtual void InitializeInvChi() = 0;

    /**
     * Return true to use variable-DOF HACApK build (per-element DOF via
     * m_dof_offset). Return false for uniform DOF (m_nffc per element).
     */
    virtual bool IsVariableDOF() const = 0;

    /**
     * Uniform DOF count per element (used when IsVariableDOF() is false).
     * Typical: 3 (MMM tetra), 6 (MSC hex), 1 (PEEC filament).
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
// RadHACApKMSCManager: MMM/MSC kernel (tetra 3DOF, wedge 5DOF, hex 6DOF)
//-------------------------------------------------------------------------

/**
 * RadHACApKMSCManager implements the HACApK kernel for Radia's Magnetic
 * Moment Method (MMM, tetrahedra, 3 DOF) and Magnetic Surface Charge
 * method (MSC, wedges 5 DOF / hexahedra 6 DOF) element types, including
 * mixed meshes.
 *
 * All kernel-specific precomputation (PrecomputeHexaGeometry, etc.),
 * flat matrix caching, and on-demand block computation routines live
 * here. The DOF-to-element lookup table m_dof_to_elem / m_dof_to_local
 * is built after ExtractCoordinates completes.
 *
 * For typical Radia compact geometries most H-matrix blocks are dense
 * (admissibility not satisfied). HACApK still provides a clean solver
 * pipeline, but the dense BiCGSTAB solver (Method 1) is often faster
 * on small/medium problems. Use HACApK when N > ~1000.
 */
class RadHACApKMSCManager : public RadHACApKBase {
public:
    explicit RadHACApKMSCManager(radTInteraction* interaction);
    ~RadHACApKMSCManager() override;

    radTInteraction* GetInteraction() const { return m_interaction; }

    /**
     * Build a local loop (cycle) deflation basis from the element-adjacency
     * graph's short cycles (3- and 4-cycles) and install it via SetDeflation,
     * so the HACApK MatVec applies (A + alpha L L^T). Adjacency and shared-face
     * DOF are derived from coincident face centroids (FaceCenter[f]); general
     * for any conforming hex/tet/wedge mesh. Call after BuildHMatrix.
     */
    void BuildLoopBasis(double alpha);

    /**
     * ALPHA-FREE loop-star gauge (tree-cotree split). Build a SPARSE "star" basis
     * T (columns = boundary face DOFs + symmetric internal sigma_A+sigma_B +
     * per-element divergence) spanning the complement of the loop null space,
     * solve the reduced NON-SINGULAR system T^T A T y = T^T b (K-dense LU / BiCGSTAB;
     * sigma_S = T y), then KEEP the loop content by block Gauss-Seidel iterative
     * refinement: alternately correct the star block (A_SS d_S = S^T r via the same
     * reduced solve) and the loop block (A_LL y_L = L^T r, A_LL = L^T diag(inv_chi) L,
     * by CG), each sweep recomputing r = b - A sigma. This converges to the DIRECT
     * solution -> FIELD-EXACT with the loops kept, regardless of star/loop basis
     * orthogonality. The H-matrix is UNCHANGED; each MatVec is sandwiched with the
     * sparse T / L (O(N)). Linear regime (uniform chi). Returns the reduced-solve
     * BiCGSTAB iteration count, or -1 on failure.
     */
    int BuildStarBasis();
    int SolveLoopStar(const std::vector<double>& b, std::vector<double>& sigma,
                      double tol, int max_iter,
                      const std::vector<double>& blockInverse = std::vector<double>(),
                      const std::vector<int>& blockOffsets = std::vector<int>());
    int GetStarDim() const { return m_n_star; }

    // True iff an ANTISYMMETRIC IMA plane (sign<0) is active. These add a few
    // LOCAL "plane-coupled" null modes to ker(N_ima) beyond the topological loops;
    // ComputePlaneSlabAddedModes() finds them and they are deflated in SolveLoopStar.
    bool HasAntisymmetricIMAPlane() const;

    // Antisym-IMA O(N) gauge correction: the topological star gauges the loops
    // (which stay null under IMA), but ker(N_ima) ALSO contains a few LOCAL modes
    // confined to the antisymmetric plane slab (count ~ plane area, decay from the
    // plane). This finds them: densify N_ima on the plane-touching elements ("slab")
    // via GetInteractionMatrixElement, local SVD, keep the null modes with weight on
    // the antisym-plane faces (plane-weight filter), orthonormalize. Returns the
    // count; caches them as an orthonormal CSR (m_added_*) for SolveLoopStar to
    // deflate (SetDeflation alpha L L^T -- robust because they are FEW + orthonormal).
    int ComputePlaneSlabAddedModes();

    double GetInteractionMatrixElement(int dof_i, int dof_j) const override;

    /**
     * MSC system-matrix convention: A = -N + delta_ij / chi_i.
     */
    double ComputeSystemEntry(int dof_i, int dof_j) const override;

    /**
     * Flatten InteractMatrix for 3DOF tetrahedra (O(1) element access).
     * Safe to call multiple times (no-op if already ready).
     */
    void PrecomputeFlatInteractMatrix();
    bool IsFlatNReady() const { return m_flat_N_ready; }

    // Delegates to radTInteraction::Compute6x6BlockFast (shared with LU/BiCGSTAB)
    void Compute6x6BlockFast(int elem_i, int elem_j, double* K_mat) const;

protected:
    void ExtractCoordinates() override;
    void OnBeforeBuild() override;
    void InitializeInvChi() override;
    bool IsVariableDOF() const override { return m_is_mixed_dof || m_is_5dof; }
    int GetUniformNFFC() const override { return m_nffc; }

private:
    // Pointer to Radia interaction (not owned)
    radTInteraction* m_interaction;

    // DOF per element classification
    int m_nffc;            // 3 (tetra), 5 (wedge), 6 (hex), 0 (mixed/variable)
    bool m_is_6dof;
    bool m_is_5dof;
    bool m_is_mixed_dof;

    // O(1) DOF-to-element lookup (ELF-style)
    std::vector<int> m_dof_to_elem;   // [dof] -> element index
    std::vector<int> m_dof_to_local;  // [dof] -> local DOF within element

    // Loop-star gauge: SPARSE star basis T stored column-wise (CSR). T maps
    // star coefficients -> sigma; its columns span the complement of the loop
    // null space (boundary + symmetric-internal + per-element divergence).
    std::vector<int> m_star_offsets;    // size m_n_star + 1 (column pointers)
    std::vector<int> m_star_dofs;       // flat sigma-DOF indices
    std::vector<double> m_star_coeffs;  // flat T entries
    int m_n_star = 0;                   // number of star columns

    // Loop (cycle) basis L for the KEEP-LOOPS back-substitution, stored
    // column-wise (CSR). Built once from the topological cycle space (the same
    // construction as BuildLoopBasis, copied out before the deflation shift is
    // cleared). After the reduced star solve gives sigma_S = S y_S, SolveLoopStar
    // recovers the loop part sigma_L = L y_L by the block back-substitution
    // A_LL y_L = L^T (b - A sigma_S),  A_LL = L^T diag(inv_chi) L. Keeping the
    // loops matches the direct LU/FEM solution; the earlier remove-loops gauge
    // (sigma = S y_S only) gave a ~0.5% external-field error on the C-type.
    std::vector<int> m_loop_offsets;    // size m_n_loop + 1 (column pointers)
    std::vector<int> m_loop_dofs;       // flat sigma-DOF indices
    std::vector<double> m_loop_coeffs;  // flat L entries (+/-1 per shared face)
    int m_n_loop = 0;                   // number of loop columns
    bool m_loop_built = false;          // cache guard

    // Antisym-IMA plane-slab added modes (orthonormal, CSR, sigma-DOF coords),
    // cached by ComputePlaneSlabAddedModes(). Deflated in SolveLoopStar via
    // SetDeflation (alpha L L^T). m_added_built guards one-time construction.
    std::vector<int> m_added_offsets;   // size m_n_added + 1
    std::vector<int> m_added_dofs;      // flat slab DOF indices
    std::vector<double> m_added_coeffs; // flat orthonormal mode entries
    int m_n_added = 0;                  // number of added modes
    bool m_added_built = false;         // cache guard

    // Pre-computed geometry for 3DOF tetrahedra (avoids B_comp overhead)
    std::vector<double> m_tetra_centers;         // [n_elem * 3]
    std::vector<double> m_tetra_face_vertices;   // [n_elem * 4 * 3 * 3]
    std::vector<double> m_tetra_face_normals;    // [n_elem * 4 * 3]
    std::vector<double> m_tetra_face_areas;      // [n_elem * 4]
    bool m_geometry_3dof_ready;

    // Flat interaction matrix storage (for 3DOF tetrahedra, O(1) access)
    std::vector<double> m_flat_N_data;           // [n_elem * n_elem * 9]
    bool m_flat_N_ready;

    // Private helpers
    void BuildDOFLookupTable();
    void PrecomputeGeometry3DOF();

    // 6DOF hexahedron block computation
    double GetCached6x6Element(int elem_i, int elem_j, int face_i, int face_j) const;
    void Compute6x6Block(int elem_i, int elem_j, double* K_mat) const;

    // 3DOF tetrahedron block computation
    double GetCached3x3Element(int elem_i, int elem_j, int comp_i, int comp_j) const;
    void Compute3x3Block(int elem_i, int elem_j, double* N_mat) const;
    void Compute3x3Block_OnDemand(int elem_i, int elem_j, double* N_mat) const;
    void Compute3x3BlockFast(int elem_i, int elem_j, double* N_mat) const;

    // 5DOF wedge / mixed / generic
    double GetCached5x5Element(int elem_i, int elem_j, int face_i, int face_j) const;
    double GetCachedMixedElement(int elem_i, int elem_j, int dof_i, int dof_j, int local_i, int local_j) const;
    double GetGenericElement(int elem_i, int elem_j, int local_i, int local_j) const;
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

    // Set inverse susceptibility for matrix element computation
    void SetInvChi(const std::vector<double>& inv_chi);

    // Get inverse susceptibility
    const std::vector<double>& GetInvChi();

    // Set/clear DOF permutation array for H-matrix build
    // Must be called BEFORE H-matrix build (lod is populated during build)
    void SetLod(int* lod, int size);
    void ClearLod();

    // Clear all global callback state (called on manager destruction)
    void ClearGlobalState();

    // Set interaction for callback (MSC kernel informational; PEEC adapters
    // may leave interaction null)
    void SetInteraction(radTInteraction* interaction, int n_elem, int nffc);

    // Cache generation counter (incremented on every H-matrix build so that
    // thread-local block caches in subclasses can detect stale entries).
    uint64_t GetGeneration();
    void IncrementGeneration();

    // Compute matrix element A(i,j) = -N(i,j) + delta_ij/chi_i
    // Called from cHACApK_entry_ij (1-based indexing in HACApK's convention).
    double ComputeEntry(int i, int j);
}

//-------------------------------------------------------------------------

#endif // __RAD_HACAPK_H
