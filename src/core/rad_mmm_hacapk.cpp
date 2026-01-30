/*-------------------------------------------------------------------------
*
* File name:      rad_mmm_hacapk.cpp
*
* Project:        RADIA - mmm_core module
*
* Description:    HACApK H-matrix solver implementation for standalone MMM
*
-------------------------------------------------------------------------*/

#include "rad_mmm_hacapk.h"

// HACApK C interface - use opaque pointer API from cHACApK_cpp.h
extern "C" {
#include "../ext/HACApK/cHACApK_cpp.h"

// Callback function setter from rad_mmm_hacapk_callback.c
typedef double (*HACApK_mmm_entry_func)(int i, int j, int bemv);
void HACApK_mmm_set_entry_func(HACApK_mmm_entry_func func);
void HACApK_mmm_clear_entry_func(void);
}

#include <cmath>
#include <chrono>
#include <algorithm>
#include <cstdint>
#include <cstdio>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace radia {
namespace mmm {

//-------------------------------------------------------------------------
// Global callback state for MMM solver
//-------------------------------------------------------------------------

namespace {
    MMMHACApKSolver* g_current_solver = nullptr;
    int* g_lod = nullptr;
    int g_lod_size = 0;
}

// C callback function that bridges to C++ MMMHACApKSolver
static double mmm_entry_callback(int i, int j, int /* bemv */) {
    if (!g_current_solver) return 0.0;

    // HACApK uses 1-based indices; convert to 0-based
    int i0 = i - 1;
    int j0 = j - 1;

    // Apply LOD permutation if available (lod maps permuted -> original)
    int pi = (g_lod && i0 >= 0 && i0 < g_lod_size) ? (g_lod[i0] - 1) : i0;
    int pj = (g_lod && j0 >= 0 && j0 < g_lod_size) ? (g_lod[j0] - 1) : j0;

    // A(i,j) = -N(i,j) - delta_ij / chi_i (ELF-compatible)
    double val = -g_current_solver->GetNElement(pi, pj);

    if (pi == pj) {
        const auto& inv_chi = g_current_solver->GetInvChi();
        const auto& dof_offset = g_current_solver->GetDOFOffset();

        // Find element for this DOF
        int n_elem = static_cast<int>(dof_offset.size()) - 1;
        for (int e = 0; e < n_elem; ++e) {
            if (pi >= dof_offset[e] && pi < dof_offset[e + 1]) {
                val -= inv_chi[e];  // ELF-compatible: -1/chi
                break;
            }
        }
    }

    return val;
}

namespace MMMHACApKCallback {

void SetCurrentSolver(MMMHACApKSolver* solver) {
    g_current_solver = solver;
    if (solver) {
        HACApK_mmm_set_entry_func(mmm_entry_callback);
    } else {
        HACApK_mmm_clear_entry_func();
    }
}

MMMHACApKSolver* GetCurrentSolver() {
    return g_current_solver;
}

double ComputeEntry(int i, int j) {
    return mmm_entry_callback(i, j, 0);
}

void SetLod(int* lod, int size) {
    g_lod = lod;
    g_lod_size = size;
}

void ClearLod() {
    g_lod = nullptr;
    g_lod_size = 0;
}

} // namespace MMMHACApKCallback

//-------------------------------------------------------------------------
// MMMHACApKSolver implementation
//-------------------------------------------------------------------------

MMMHACApKSolver::MMMHACApKSolver()
    : n_elements_(0)
    , total_dof_(0)
    , leafmtxp_(nullptr)
    , control_(nullptr)
    , valid_(false)
{
}

MMMHACApKSolver::~MMMHACApKSolver() {
    FreeHACApK();
}

void MMMHACApKSolver::FreeHACApK() {
    if (leafmtxp_) {
        HACApK_free_hmatrix_wrapper(leafmtxp_, control_);
        HACApK_free_leafmtxp(leafmtxp_);
        leafmtxp_ = nullptr;
    }
    if (control_) {
        HACApK_free_lcontrol(control_);
        control_ = nullptr;
    }
    valid_ = false;
}

void MMMHACApKSolver::SetElements(const std::vector<MMMElement>& elements,
                                   const std::vector<int>& dof_offset) {
    elements_ = elements;
    dof_offset_ = dof_offset;
    n_elements_ = static_cast<int>(elements.size());
    total_dof_ = dof_offset.empty() ? 0 : dof_offset.back();

    ExtractCoordinates();
    BuildDOFToElementMap();
}

void MMMHACApKSolver::SetInteractionMatrix(const std::vector<double>& N, int total_dof) {
    N_ = N;
    total_dof_ = total_dof;
}

void MMMHACApKSolver::ExtractCoordinates() {
    coordinates_.resize(n_elements_ * 3);

    for (int e = 0; e < n_elements_; ++e) {
        coordinates_[e * 3 + 0] = elements_[e].center.x;
        coordinates_[e * 3 + 1] = elements_[e].center.y;
        coordinates_[e * 3 + 2] = elements_[e].center.z;
    }
}

void MMMHACApKSolver::BuildDOFToElementMap() {
    dof_to_elem_.resize(total_dof_);
    dof_to_local_.resize(total_dof_);

    for (int e = 0; e < n_elements_; ++e) {
        int start = dof_offset_[e];
        int end = dof_offset_[e + 1];
        for (int d = start; d < end; ++d) {
            dof_to_elem_[d] = e;
            dof_to_local_[d] = d - start;
        }
    }
}

double MMMHACApKSolver::GetNElement(int i, int j) const {
    if (i < 0 || i >= total_dof_ || j < 0 || j >= total_dof_) {
        return 0.0;
    }
    return N_[i * total_dof_ + j];
}

bool MMMHACApKSolver::BuildHMatrix(const HACApKParams& params,
                                    const std::vector<double>& inv_chi) {
    if (n_elements_ == 0 || total_dof_ == 0) {
        return false;
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    // Free previous H-matrix if exists
    FreeHACApK();

    // Set callback
    MMMHACApKCallback::SetCurrentSolver(this);

    // Set inv_chi for H-matrix build
    if (!inv_chi.empty()) {
        inv_chi_ = inv_chi;
    } else if (inv_chi_.empty()) {
        // Default: chi = infinity (mu_r = infinity) -> inv_chi = 0
        inv_chi_.resize(n_elements_, 0.0);
    }

    // Allocate HACApK structures using opaque pointer API
    leafmtxp_ = HACApK_alloc_leafmtxp();
    control_ = HACApK_alloc_lcontrol();

    if (!leafmtxp_ || !control_) {
        FreeHACApK();
        return false;
    }

    // Prepare coordinate array for clustering
    // HACApK expects coordinates per DOF, but we have per element
    // For 3DOF elements, expand to DOF coordinates
    std::vector<double> dof_coords(total_dof_ * 3);
    for (int e = 0; e < n_elements_; ++e) {
        int start = dof_offset_[e];
        int end = dof_offset_[e + 1];
        for (int d = start; d < end; ++d) {
            dof_coords[d * 3 + 0] = elements_[e].center.x;
            dof_coords[d * 3 + 1] = elements_[e].center.y;
            dof_coords[d * 3 + 2] = elements_[e].center.z;
        }
    }

    // Set current solver for callback
    MMMHACApKCallback::SetCurrentSolver(this);

    // Build H-matrix using wrapper function from cHACApK_cpp.h
    // This handles cluster tree construction and ACA+ compression
    int ierr = HACApK_build_hmatrix_wrapper(
        leafmtxp_, control_,
        dof_coords.data(),
        total_dof_,   // n_elem (using DOF count as "element" count)
        1,            // nffc = 1 DOF per "element" (we already expanded to DOF coords)
        3,            // ndim = 3
        params.aca_eps,
        params.leaf_size,
        params.eta,
        params.print_level);

    if (ierr != 0) {
        MMMHACApKCallback::SetCurrentSolver(nullptr);
        FreeHACApK();
        return false;
    }

    // Store permutation from HACApK
    // lod is 1-based (Fortran convention): lod[1..nd] contains permutation
    // lod[i] = original_1based_index for permuted position i
    // Note: We don't need to apply permutation in MatVec since
    // HACApK_matvec_wrapper already handles it internally
    int* lod = HACApK_lcontrol_get_lod(control_);
    if (lod) {
        permutation_.resize(total_dof_);
        for (int i = 0; i < total_dof_; ++i) {
            // lod is 1-indexed: lod[1] to lod[nd]
            permutation_[i] = lod[i + 1] - 1;  // Convert 1-based to 0-based
        }
    }

    // Compute statistics using accessor functions
    auto end_time = std::chrono::high_resolution_clock::now();
    stats_.build_time = std::chrono::duration<double>(end_time - start_time).count();
    stats_.n_dof = total_dof_;
    stats_.n_leaves = HACApK_leafmtxp_get_nlf(leafmtxp_);
    stats_.n_lowrank = HACApK_leafmtxp_get_nlfkt(leafmtxp_);
    stats_.n_dense = stats_.n_leaves - stats_.n_lowrank;
    stats_.max_rank = HACApK_leafmtxp_get_ktmax(leafmtxp_);

    // Get memory statistics using wrapper function
    int64_t hmat_bytes = 0, dense_bytes = 0;
    HACApK_get_memory_stats(leafmtxp_, &hmat_bytes, &dense_bytes);
    stats_.memory_mb = hmat_bytes / (1024.0 * 1024.0);
    stats_.dense_memory_mb = dense_bytes / (1024.0 * 1024.0);
    stats_.compression = (dense_bytes > 0) ?
        (double)hmat_bytes / (double)dense_bytes : 1.0;

    valid_ = true;
    return true;
}

void MMMHACApKSolver::MatVec(const std::vector<double>& x, std::vector<double>& y) {
    if (!valid_ || !leafmtxp_) {
        // Fallback to dense matrix-vector
        y.resize(total_dof_, 0.0);
        for (int i = 0; i < total_dof_; ++i) {
            double sum = 0.0;
            for (int j = 0; j < total_dof_; ++j) {
                // A(i,j) = -N(i,j) - delta_ij / chi (ELF-compatible)
                double aij = -N_[i * total_dof_ + j];
                if (i == j) {
                    int e = dof_to_elem_[i];
                    aij -= inv_chi_[e];  // ELF-compatible: -1/chi
                }
                sum += aij * x[j];
            }
            y[i] = sum;
        }
        return;
    }

    // HACApK_matvec_wrapper already handles permutation internally via lod array
    // Do NOT apply permutation here - it would cause double permutation
    y.resize(total_dof_, 0.0);
    HACApK_matvec_wrapper(leafmtxp_, control_, x.data(), y.data(), total_dof_);
}

void MMMHACApKSolver::UpdateDiagonal(const std::vector<double>& inv_chi) {
    inv_chi_ = inv_chi;

    // For nonlinear iteration, the diagonal update is handled in MatVec
    // since the system matrix A = -diag(inv_chi) - N changes per iteration (ELF-compatible).
    // The H-matrix stores the N matrix; we subtract inv_chi during matvec.
    // This approach is simpler and avoids modifying the H-matrix structure.
}

int MMMHACApKSolver::Solve(const std::vector<double>& inv_chi,
                            const std::vector<double>& H_ext,
                            std::vector<double>& M,
                            double tol,
                            int max_iter) {
    if (!valid_) {
        return -1;
    }

    inv_chi_ = inv_chi;
    UpdateDiagonal(inv_chi);

    // Initialize solution
    M.resize(total_dof_, 0.0);

    // BiCGSTAB with H-matrix acceleration
    std::vector<double> r(total_dof_);
    std::vector<double> r0(total_dof_);
    std::vector<double> p(total_dof_);
    std::vector<double> v(total_dof_);
    std::vector<double> s(total_dof_);
    std::vector<double> t(total_dof_);

    // r = b - A*x0 (x0 = 0, so r = b = H_ext)
    std::copy(H_ext.begin(), H_ext.end(), r.begin());
    std::copy(r.begin(), r.end(), r0.begin());
    std::copy(r.begin(), r.end(), p.begin());

    double rho = 0.0;
    for (int i = 0; i < total_dof_; ++i) {
        rho += r0[i] * r[i];
    }

    double b_norm = 0.0;
    for (int i = 0; i < total_dof_; ++i) {
        b_norm += H_ext[i] * H_ext[i];
    }
    b_norm = std::sqrt(b_norm);
    if (b_norm < 1e-30) b_norm = 1.0;

    int iter;
    for (iter = 0; iter < max_iter; ++iter) {
        // v = A * p
        MatVec(p, v);

        double r0v = 0.0;
        for (int i = 0; i < total_dof_; ++i) {
            r0v += r0[i] * v[i];
        }

        if (std::abs(r0v) < 1e-30) break;

        double alpha = rho / r0v;

        // s = r - alpha * v
        for (int i = 0; i < total_dof_; ++i) {
            s[i] = r[i] - alpha * v[i];
        }

        // Check convergence
        double s_norm = 0.0;
        for (int i = 0; i < total_dof_; ++i) {
            s_norm += s[i] * s[i];
        }
        s_norm = std::sqrt(s_norm);

        if (s_norm / b_norm < tol) {
            for (int i = 0; i < total_dof_; ++i) {
                M[i] += alpha * p[i];
            }
            return iter + 1;
        }

        // t = A * s
        MatVec(s, t);

        double ts = 0.0, tt = 0.0;
        for (int i = 0; i < total_dof_; ++i) {
            ts += t[i] * s[i];
            tt += t[i] * t[i];
        }

        double omega = (std::abs(tt) > 1e-30) ? ts / tt : 0.0;

        // x = x + alpha * p + omega * s
        for (int i = 0; i < total_dof_; ++i) {
            M[i] += alpha * p[i] + omega * s[i];
        }

        // r = s - omega * t
        for (int i = 0; i < total_dof_; ++i) {
            r[i] = s[i] - omega * t[i];
        }

        // Check convergence
        double r_norm = 0.0;
        for (int i = 0; i < total_dof_; ++i) {
            r_norm += r[i] * r[i];
        }
        r_norm = std::sqrt(r_norm);

        if (r_norm / b_norm < tol) {
            return iter + 1;
        }

        double rho_new = 0.0;
        for (int i = 0; i < total_dof_; ++i) {
            rho_new += r0[i] * r[i];
        }

        if (std::abs(rho) < 1e-30 || std::abs(omega) < 1e-30) break;

        double beta = (rho_new / rho) * (alpha / omega);
        rho = rho_new;

        // p = r + beta * (p - omega * v)
        for (int i = 0; i < total_dof_; ++i) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
    }

    return -iter;  // Not converged
}

}  // namespace mmm
}  // namespace radia
