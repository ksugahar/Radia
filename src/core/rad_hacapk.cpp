/*-------------------------------------------------------------------------
*
* File name:      rad_hacapk.cpp
*
* Project:        RADIA
*
* Description:    HACApK (H-matrix with ACA+) interface for BiCGSTAB solver
*                 Implementation of the kernel-agnostic RadHACApKBase and callback functions
*
* First release:  2025
*
* Reference:      ppOpen-HPC project (MIT License)
*                 https://github.com/ppohHPC/ppOpen-HPC
*
-------------------------------------------------------------------------*/

#include "rad_hacapk.h"
#include <algorithm>
#include <iostream>
#include <chrono>
#include <exception>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <utility>
#include "rad_parallel.h"

// Include C++ compatible HACApK wrapper header
extern "C" {
#include "../ext/HACApK/cHACApK_cpp.h"
}

//=========================================================================
// Global callback state (required by HACApK C interface)
//=========================================================================

namespace {
    template <class F>
    class ScopeExit {
    public:
        explicit ScopeExit(F&& callback) : m_callback(std::move(callback)) {}
        ~ScopeExit() { m_callback(); }
        ScopeExit(const ScopeExit&) = delete;
        ScopeExit& operator=(const ScopeExit&) = delete;

    private:
        F m_callback;
    };

    template <class F>
    ScopeExit(F) -> ScopeExit<F>;

    // Global callback state shared across all TaskManager worker threads.
    // NOT thread_local: HACApK calls cHACApK_entry_ij from ngcore::ParallelFor
    // (TaskManager) worker threads, which must see the manager set by the main
    // thread.
    // Python releases the GIL around every H-matrix build.  HACApK's C callback
    // ABI nevertheless requires one process-wide manager (and the symmetric
    // fill switch is process-wide too), so serialize builds while retaining
    // TaskManager parallelism inside each build.
    std::mutex g_hacapkBuildMutex;
    std::mutex g_hacapkCallbackExceptionMutex;
    std::exception_ptr g_hacapkCallbackException;
    RadHACApKBase* g_currentManager = nullptr;
}

namespace RadHACApKCallback {

std::mutex& OperationMutex() {
    return g_hacapkBuildMutex;
}

void SetCurrentManager(RadHACApKBase* manager) {
    g_currentManager = manager;
}

RadHACApKBase* GetCurrentManager() {
    return g_currentManager;
}

void ClearGlobalState() {
    // Clear callback state to prevent interference between solves.
    g_currentManager = nullptr;
    {
        std::lock_guard<std::mutex> lock(g_hacapkCallbackExceptionMutex);
        g_hacapkCallbackException = nullptr;
    }
}

void ClearCallbackException() {
    std::lock_guard<std::mutex> lock(g_hacapkCallbackExceptionMutex);
    g_hacapkCallbackException = nullptr;
}

void CaptureCallbackException() noexcept {
    try {
        std::lock_guard<std::mutex> lock(g_hacapkCallbackExceptionMutex);
        if (!g_hacapkCallbackException)
            g_hacapkCallbackException = std::current_exception();
    }
    catch (...) {
        // Never allow error bookkeeping itself to unwind through HACApK C
        // code or the TaskManager worker trampoline.
    }
}

void RethrowCallbackException() {
    std::exception_ptr failure;
    {
        std::lock_guard<std::mutex> lock(g_hacapkCallbackExceptionMutex);
        failure = g_hacapkCallbackException;
        g_hacapkCallbackException = nullptr;
    }
    if (failure) std::rethrow_exception(failure);
}

double ComputeEntry(int i, int j) {
    // i, j are 1-based ORIGINAL indices from HACApK.
    // The lod conversion is already done by cHACApK_fill_leafmtx_hyp:
    //   val = cHACApK_entry_ij(lodl[permuted_pos], lodt[permuted_pos], i_bemv)
    // so we receive original indices, NOT permuted indices.
    //
    // The kernel-specific system-matrix convention is delegated to
    // RadHACApKBase::ComputeSystemEntry so that each HDiv, PEEC, or BEM
    // subclass can store exactly what HACApK needs.

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
// of the default magnetostatic system matrix.  This lets callers (e.g. the
// stream-function (ACA+)+TSVD solver) factor an arbitrary rectangular kernel
// block with HACApK's ACA+ -- keeping ACA+ a single source of truth instead of
// re-porting it.  Default null => unchanged magnetostatic behaviour.  Set/cleared
// synchronously around one factorization (GIL-serialized; no concurrent
// build), so a plain (non-thread_local) pointer is sufficient.
static HACApK_entry_func g_entry_override = NULL;

void HACApK_set_entry_func(HACApK_entry_func func) { g_entry_override = func; }
void HACApK_clear_entry_func(void) { g_entry_override = NULL; }

double cHACApK_entry_ij(int i, int j, int i_bemv) {
    try {
        if (g_entry_override != NULL) return g_entry_override(i, j, i_bemv);
        (void)i_bemv;  // Unused in Radia
        return RadHACApKCallback::ComputeEntry(i, j);
    }
    catch (...) {
        // A C++ exception cannot safely unwind through HACApK's C fill or the
        // TaskManager worker trampoline. Preserve the first failure and
        // rethrow it on the build thread after all fill workers have joined.
        RadHACApKCallback::CaptureCallbackException();
        return std::numeric_limits<double>::quiet_NaN();
    }
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
    // Destruction mutates the same process-wide callback and HACApK C state
    // as BuildHMatrix. A manager released by another Python thread must not
    // clear that state while a different manager is filling its leaves.
    std::lock_guard<std::mutex> build_lock(RadHACApKCallback::OperationMutex());
    FreeResources();
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

bool RadHACApKBase::BuildHMatrix(const RadHACApKParams& params) {
    // The C callback ABI stores the current manager, inverse material data,
    // and symmetric-fill mode in process-wide state.  pybind deliberately
    // releases the GIL for this operation, so the C++ boundary must provide
    // the serialization.  The lock covers setup, fill, and callback-backed
    // diagonal extraction; the fill itself remains TaskManager-parallel.
    std::lock_guard<std::mutex> build_lock(RadHACApKCallback::OperationMutex());

    bool build_succeeded = false;
    cHACApK_set_sym_fill(UseSymmetricFill() ? 1 : 0);
    ScopeExit build_state_scope([this, &build_succeeded]() noexcept {
        cHACApK_set_sym_fill(0);
        OnBuildFinished(build_succeeded);
    });
    OnBuildStarting(params);

    // TaskManager self-wrap (AGENTS.md "Parallelization: NGSolve TaskManager"): the H-matrix leaf
    // fill runs ngcore::ParallelFor, which silently falls back to single-threaded when NO
    // RegionTaskManager is active.  Stand up (or reuse the caller's) pool here so EVERY
    // HACApK build -- HDiv, PEEC, compact magnetostatic kernels, diagnostics -- is parallel even when
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

    // Set global callback state before kernel-specific setup.
    RadHACApKCallback::SetCurrentManager(this);
    RadHACApKCallback::ClearCallbackException();

    // Kernel-specific precomputation.
    OnBeforeBuild();

    // Kernel-specific initial chi.
    InitializeInvChi();

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
            params.max_rank,
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
            params.max_rank,
            params.print_level
        );
    }
    auto t_hmatrix_end = std::chrono::high_resolution_clock::now();
    try {
        RadHACApKCallback::RethrowCallbackException();
    }
    catch (...) {
        FreeResources();
        throw;
    }
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
    build_succeeded = true;
    return true;
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

void RadHACApKBase::MatVecTranspose(const std::vector<double>& x, std::vector<double>& y) {
    if (!m_valid || !m_leafmtxp || !m_control) {
        std::fill(y.begin(), y.end(), 0.0);
        return;
    }
    int nd = HACApK_leafmtxp_get_nd(m_leafmtxp);
    HACApK_matvec_transpose_wrapper(m_leafmtxp, m_control, x.data(), y.data(), nd);
    // (loop-mode deflation L L^T is symmetric, so its transpose contribution is identical)
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

void RadHACApKBase::MatVecSym(const std::vector<double>& x, std::vector<double>& y) {
    if (!m_valid || !m_leafmtxp || !m_control) {
        std::fill(y.begin(), y.end(), 0.0);
        return;
    }
    int nd = HACApK_leafmtxp_get_nd(m_leafmtxp);
    HACApK_matvec_sym_wrapper(m_leafmtxp, m_control, x.data(), y.data(), nd);
    if (m_defl_nplaq > 0 && m_defl_alpha != 0.0) {   // L L^T is already symmetric
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

void RadHACApKBase::MatVecSymMany(
    const std::vector<double>& x, int nrhs, std::vector<double>& y)
{
    MatVecSymManyPrepared(x, nrhs, nullptr, nullptr, y);
}

void RadHACApKBase::MatVecSymManyPrepared(
    const std::vector<double>& x, int nrhs,
    const int* active_prefix, const double* diagonal_scale,
    std::vector<double>& y)
{
    if (!m_valid || !m_leafmtxp || !m_control || nrhs < 1 ||
        x.size() != static_cast<size_t>(nrhs)*m_ndof ||
        (active_prefix && active_prefix[0] != 0))
        throw std::runtime_error(
            "MatVecSymManyPrepared: invalid operator, batch, or active prefix");
    y.assign(x.size(), 0.0);
    HACApK_matvec_sym_many_prepared_wrapper(
        m_leafmtxp, m_control, x.data(), y.data(), m_ndof, nrhs,
        active_prefix, diagonal_scale);
    if (m_defl_nplaq > 0 && m_defl_alpha != 0.0) {
        std::vector<unsigned char> active;
        if (active_prefix) {
            const auto* control =
                static_cast<const st_cHACApK_lcontrol_t*>(m_control);
            active.assign(static_cast<size_t>(m_ndof), 0);
            for (int permuted = 0; permuted < m_ndof; ++permuted)
                if (active_prefix[permuted+1] != active_prefix[permuted])
                    active[static_cast<size_t>(
                        control->lod[permuted+1]-1)] = 1;
        }
        for (int rhs = 0; rhs < nrhs; ++rhs)
            for (int p = 0; p < m_defl_nplaq; ++p) {
                double c = 0.0;
                for (int k = m_defl_offsets[p]; k < m_defl_offsets[p+1]; ++k) {
                    const int dof = m_defl_dofs[k];
                    if (!active_prefix ||
                        active[static_cast<size_t>(dof)])
                        c += m_defl_signs[k] *
                            x[static_cast<size_t>(rhs)*m_ndof+dof] *
                            (diagonal_scale ? diagonal_scale[dof] : 1.0);
                }
                c *= m_defl_alpha;
                for (int k = m_defl_offsets[p]; k < m_defl_offsets[p+1]; ++k) {
                    const int dof = m_defl_dofs[k];
                    if (!active_prefix ||
                        active[static_cast<size_t>(dof)])
                        y[static_cast<size_t>(rhs)*m_ndof+dof] +=
                            m_defl_signs[k]*c *
                            (diagonal_scale ? diagonal_scale[dof] : 1.0);
                }
            }
    }
}

void RadHACApKBase::MatVecSymManyMasked(
    const std::vector<double>& x, int nrhs,
    const std::vector<int>& active_prefix, std::vector<double>& y)
{
    if (!m_valid || !m_leafmtxp || !m_control || nrhs < 1 ||
        x.size() != static_cast<size_t>(nrhs)*m_ndof ||
        active_prefix.size() != static_cast<size_t>(m_ndof)+1 ||
        active_prefix.front() != 0)
        throw std::runtime_error(
            "MatVecSymManyMasked: invalid operator, batch, or active prefix");
    MatVecSymManyPrepared(
        x, nrhs, active_prefix.data(), nullptr, y);
}

void RadHACApKBase::UpdateDiagonal(const std::vector<double>& inv_chi) {
    if (!m_valid || !m_leafmtxp || !m_control) return;

    // Update stored inverse susceptibility
    m_inv_chi = inv_chi;

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
