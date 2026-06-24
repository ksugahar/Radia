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
*                 - RadHACApKMMMManager : public RadHACApKBase implements the
*                   MMM 3-DOF tetrahedron kernel. Surface-charge MSC is handled
*                   by the moment-yano RadHACApKMomentSystem. A future RadHACApKPEECManager will
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
 *   RadHACApKMMMManager mgr(interaction);
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
// RadHACApKMMMManager: MMM kernel (tetra 3DOF)
//-------------------------------------------------------------------------

/**
 * RadHACApKMMMManager implements the HACApK kernel for Radia's Magnetic
 * Moment Method (MMM, tetrahedra, 3 DOF).  Magnetic Surface Charge
 * elements (hex/wedge/pyramid, 5-6 DOF) use the separate moment-yano
 * RadHACApKMomentSystem or dense moment LU path.
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
class RadHACApKMMMManager : public RadHACApKBase {
public:
    explicit RadHACApKMMMManager(radTInteraction* interaction);
    ~RadHACApKMMMManager() override;

    radTInteraction* GetInteraction() const { return m_interaction; }


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

protected:
    void ExtractCoordinates() override;
    void OnBeforeBuild() override;
    void InitializeInvChi() override;
    bool IsVariableDOF() const override { return false; }   // MMM-only: uniform 3-DOF tet
    int GetUniformNFFC() const override { return m_nffc; }   // == 3

private:
    // Pointer to Radia interaction (not owned)
    radTInteraction* m_interaction;

    // DOF per element: MMM-only manager, always 3 (tetrahedron)
    int m_nffc;

    // O(1) DOF-to-element lookup (ELF-style)
    std::vector<int> m_dof_to_elem;   // [dof] -> element index
    std::vector<int> m_dof_to_local;  // [dof] -> local DOF within element


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

    // 3DOF tetrahedron block computation (MMM -- the only element type this manager solves;
    // EIEM2 surface-charge 6x6/5x5/mixed kernels were retired in Phase 3b, the moment-yano
    // H-matrix RadHACApKMomentSystem now owns hex/wedge/pyramid MSC)
    double GetCached3x3Element(int elem_i, int elem_j, int comp_i, int comp_j) const;
    void Compute3x3Block(int elem_i, int elem_j, double* N_mat) const;
    void Compute3x3Block_OnDemand(int elem_i, int elem_j, double* N_mat) const;
    void Compute3x3BlockFast(int elem_i, int elem_j, double* N_mat) const;

    double GetGenericElement(int elem_i, int elem_j, int local_i, int local_j) const;
};

//-------------------------------------------------------------------------
// RadHACApKMomentSystem: the parameter-free MOMENT-yano system A_raw as a HACApK
// H-matrix (Phase 2 of the EIEM2 full-deletion track; docs/moment_yano/ACA_MOMENT_DESIGN.md).
//-------------------------------------------------------------------------

/* The moment system A_raw = L(block-diag local moment) - chi*C(centroid field/grad coupling).
 * The off-diagonal block of well-separated element clusters is the smooth field/grad kernel
 * folded by the per-row moment functionals -> low-rank (Gate 1), so A_raw is an H-matrix with
 * the cluster tree over element (hex) centroids: one element = 6 DOF (3 dipole + 1 monopole +
 * 2 quadrupole rows on the row side, 6 face charges on the column side, co-located).  The entry
 * A_raw[i][j] is computed ON DEMAND by radTInteraction::MomentSystemEntry (no dense build, no row
 * normalization -- the row 2-norm is a diagonal scaling that leaves the direct solve invariant).
 * A_raw is NON-symmetric (rows = moment functionals, cols = charges) -- ACA+ compresses it anyway.
 * ComputeSystemEntry stores A_raw directly (no -N/+1/chi flip); the H-LU (cHACApK_hlu_*) factors it
 * (Increment 3).  HEX-ONLY; assumes m_elemDOFOffset[m_hexaElemIndices[h]] == 6*h (pure-hex moment). */
class RadHACApKMomentSystem : public RadHACApKBase {
public:
    RadHACApKMomentSystem(radTInteraction* interaction, double chi);                            // uniform chi
    RadHACApKMomentSystem(radTInteraction* interaction, const std::vector<double>& chiPerHex);  // per-element chi (Increment 4)
    ~RadHACApKMomentSystem() override {}

    radTInteraction* GetInteraction() const { return m_interaction; }

    // A_raw[i][j] on demand (the un-normalized moment system entry; rows 6*h+t, cols = face DOF).
    double GetInteractionMatrixElement(int dof_i, int dof_j) const override;
    // The H-matrix stores A_raw directly (no MSC sign flip / 1-chi shift).
    double ComputeSystemEntry(int dof_i, int dof_j) const override { return GetInteractionMatrixElement(dof_i, dof_j); }

protected:
    void ExtractCoordinates() override;   // cluster tree = hex centroids; ndof = 6*nHex
    void OnBeforeBuild() override {}
    void InitializeInvChi() override { m_inv_chi.assign(m_ndof, 0.0); }   // chi folded into A_raw
    bool IsVariableDOF() const override { return false; }
    int  GetUniformNFFC() const override { return 6; }                    // 6 DOF per hex

private:
    radTInteraction* m_interaction;   // not owned
    double m_chi;                     // uniform chi (fallback when m_chi_in is empty)
    std::vector<double> m_chi_in;     // per-element chi supplied by the ctor (Increment 4); empty -> uniform m_chi
    std::vector<double> m_chiv;       // chi per hex, resolved in ExtractCoordinates, for MomentSystemEntry
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
