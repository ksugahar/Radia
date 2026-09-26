/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/// @file compact_amg.hpp
/// @brief Compact Algebraic Multigrid (AMG) preconditioner using NGSolve TaskManager
///
/// Compact AMG implementation (~700 lines) that uses
/// NGSolve's SparseMatrix, TaskManager, and InverseMatrix infrastructure.
/// All parallelism via TaskManager (no OpenMP).
///
/// Algorithm: Classical AMG with PMIS coarsening + direct interpolation + l1-GS smoother.
/// Reference: Ruge & Stueben (1987), De Sterck et al. (2006) for PMIS.

#ifndef SPARSESOLV_COMPACT_AMG_HPP
#define SPARSESOLV_COMPACT_AMG_HPP

#include <comp.hpp>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <atomic>
#include <array>
#include <chrono>

namespace ngla {

inline bool GalerkinProductOnPattern(const SparseMatrix<double>& Pt, const SparseMatrix<double>& A,
                                     const SparseMatrix<double>& P, SparseMatrix<double>& C);
inline shared_ptr<SparseMatrix<double>> GalerkinPattern(const SparseMatrix<double>& Pt,
                                                        const SparseMatrix<double>& A,
                                                        const SparseMatrix<double>& P);

/// Pt * A * P with a native symbolic pass and the numeric pass on its pattern.
inline shared_ptr<SparseMatrix<double>> GalerkinProduct(const SparseMatrix<double>& Pt,
                                                        const SparseMatrix<double>& A,
                                                        const SparseMatrix<double>& P) {
    auto C = GalerkinPattern(Pt, A, P);
    if (!GalerkinProductOnPattern(Pt, A, P, *C))
        throw std::runtime_error("GalerkinProduct: contribution outside its own pattern");
    return C;
}

/// res = b - A x in one pass over the rows. With `values` (a float32 mirror of
/// A's values in storage order) the products read the mirror: half the value
/// traffic; vectors and accumulation stay double.
inline void ResidualWithValues(const SparseMatrix<double>& A, const float* values,
                               const BaseVector& b, const BaseVector& x, BaseVector& res) {
    auto fb = b.FVDouble();
    auto fx = x.FVDouble();
    auto fr = res.FVDouble();
    ParallelForRange(A.Height(), [&](IntRange range) {
        for (auto i : range) {
            auto cols = A.GetRowIndices(i);
            double s = fb[i];
            if (values) {
                const float* v = values + A.First(i);
                for (int j = 0; j < cols.Size(); j++) s -= double(v[j]) * fx[cols[j]];
            } else {
                auto vals = A.GetRowValues(i);
                for (int j = 0; j < cols.Size(); j++) s -= vals[j] * fx[cols[j]];
            }
            fr[i] = s;
        }
    });
}

/// float32 copy of A's values in storage order (parallel).
inline void MirrorValues(const SparseMatrix<double>& A, std::vector<float>& out) {
    auto values = A.AsVector().FVDouble();
    out.resize(values.Size());
    ParallelForRange(values.Size(), [&](IntRange range) {
        for (auto k : range) out[k] = float(values[k]);
    });
}

/// Compact CSR graph (binary adjacency, no values)
struct CSRGraph {
    int n = 0;
    std::vector<int> row_ptr;
    std::vector<int> col_idx;

    int NumNeighbors(int i) const { return row_ptr[i + 1] - row_ptr[i]; }

    const int* NeighborBegin(int i) const { return col_idx.data() + row_ptr[i]; }
    const int* NeighborEnd(int i) const { return col_idx.data() + row_ptr[i + 1]; }
};

/// True when a and b have the same sparsity pattern (row lengths and column indices).
inline bool SameSparsityPattern(const SparseMatrix<double>& a, const SparseMatrix<double>& b) {
    if (&a == &b) return true;
    if (a.Height() != b.Height() || a.Width() != b.Width() || a.NZE() != b.NZE()) return false;
    std::atomic<bool> same{true};
    ParallelForRange(a.Height(), [&](IntRange range) {
        for (auto i : range) {
            auto ca = a.GetRowIndices(i);
            auto cb = b.GetRowIndices(i);
            if (ca.Size() != cb.Size()
                || !std::equal(ca.Data(), ca.Data() + ca.Size(), cb.Data())) {
                same.store(false, std::memory_order_relaxed);
                return;
            }
        }
    });
    return same.load();
}

/// Overwrite the values of B (which has A's pattern) with A, using identity rows
/// for constrained dofs and zero couplings to them -- the values of
/// CreateBCModifiedMatrix without allocating a new matrix.
inline void FillBCModifiedValues(const SparseMatrix<double>& A, const BitArray& freedofs,
                                 SparseMatrix<double>& B) {
    ParallelFor(A.Height(), [&](size_t i) {
        auto cols = A.GetRowIndices(i);
        auto src = A.GetRowValues(i);
        auto dst = B.GetRowValues(i);
        const bool free_row = freedofs.Test(i);
        for (int j = 0; j < cols.Size(); j++) {
            if (!free_row)
                dst[j] = (cols[j] == (int)i) ? 1.0 : 0.0;
            else
                dst[j] = freedofs.Test(cols[j]) ? src[j] : 0.0;
        }
    });
}

/// Structural pattern of Pt * A * P as a zero SparseMatrix (sorted columns):
/// row u holds every j reached by u -Pt-> e -A-> k -P-> j. Two parallel passes
/// (count, fill) with a per-task column marker; no values are computed.
inline shared_ptr<SparseMatrix<double>> GalerkinPattern(const SparseMatrix<double>& Pt,
                                                        const SparseMatrix<double>& A,
                                                        const SparseMatrix<double>& P) {
    if (Pt.Width() != A.Height() || A.Width() != P.Height())
        throw std::invalid_argument("GalerkinPattern: dimension mismatch");
    const size_t height = Pt.Height();
    const int width = P.Width();
    Array<int> counts(height);
    auto visit = [&](size_t u, std::vector<int>& marker, auto&& emit) {
        auto ecols = Pt.GetRowIndices(u);
        for (int a = 0; a < ecols.Size(); a++) {
            auto kcols = A.GetRowIndices(ecols[a]);
            for (int b = 0; b < kcols.Size(); b++) {
                auto jcols = P.GetRowIndices(kcols[b]);
                for (int c = 0; c < jcols.Size(); c++) {
                    const int j = jcols[c];
                    if (marker[j] != int(u)) { marker[j] = int(u); emit(j); }
                }
            }
        }
    };
    ParallelForRange(height, [&](IntRange range) {
        std::vector<int> marker(width, -1);
        for (auto u : range) {
            int n = 0;
            visit(u, marker, [&](int) { n++; });
            counts[u] = n;
        }
    });
    auto C = make_shared<SparseMatrix<double>>(counts, width);
    ParallelForRange(height, [&](IntRange range) {
        std::vector<int> marker(width, -1);
        for (auto u : range) {
            auto cols = C->GetRowIndices(u);
            int n = 0;
            visit(u, marker, [&](int j) { cols[n++] = j; });
            std::sort(cols.Data(), cols.Data() + cols.Size());
            C->GetRowValues(u) = 0.0;
        }
    });
    return C;
}

/// Galerkin product C = Pt * A * P computed into the existing pattern of C
/// (numeric phase only: no symbolic product, no allocation). Row u of C is
/// sum_e Pt[u,e] sum_k A[e,k] P[k,:], accumulated through a per-task column
/// marker. Returns false when a nonzero contribution falls outside C's
/// pattern; the values of C are then unspecified and the caller must rebuild.
inline bool GalerkinProductOnPattern(const SparseMatrix<double>& Pt, const SparseMatrix<double>& A,
                                     const SparseMatrix<double>& P, SparseMatrix<double>& C) {
    if (Pt.Height() != C.Height() || Pt.Width() != A.Height() || A.Width() != P.Height()
        || P.Width() != C.Width())
        throw std::invalid_argument("GalerkinProductOnPattern: dimension mismatch");
    std::atomic<bool> inside{true};
    const int width = C.Width();
    ParallelForRange(C.Height(), [&](IntRange range) {
        std::vector<int> marker(width, -1);
        for (auto u : range) {
            auto ccols = C.GetRowIndices(u);
            auto cvals = C.GetRowValues(u);
            for (int p = 0; p < ccols.Size(); p++) {
                marker[ccols[p]] = p;
                cvals[p] = 0.0;
            }
            auto ecols = Pt.GetRowIndices(u);
            auto evals = Pt.GetRowValues(u);
            for (int a = 0; a < ecols.Size(); a++) {
                const double pe = evals[a];
                if (pe == 0.0) continue;
                const int e = ecols[a];
                auto kcols = A.GetRowIndices(e);
                auto kvals = A.GetRowValues(e);
                for (int b = 0; b < kcols.Size(); b++) {
                    const double w = pe * kvals[b];
                    if (w == 0.0) continue;
                    const int k = kcols[b];
                    auto jcols = P.GetRowIndices(k);
                    auto jvals = P.GetRowValues(k);
                    for (int c = 0; c < jcols.Size(); c++) {
                        const int pos = marker[jcols[c]];
                        if (pos < 0) {
                            if (jvals[c] != 0.0) inside.store(false, std::memory_order_relaxed);
                            continue;
                        }
                        cvals[pos] += w * jvals[c];
                    }
                }
            }
            for (int p = 0; p < ccols.Size(); p++) marker[ccols[p]] = -1;
        }
    });
    return inside.load();
}

/// Compact Algebraic Multigrid preconditioner for scalar H1 problems.
///
/// Uses NGSolve SparseMatrix for all matrix operations (Restrict, Mult, InverseMatrix).
/// All loops parallelized via NGSolve TaskManager.
///
/// Usage:
///   auto amg = make_shared<CompactAMG>(mat, freedofs, 0.25, 25, 50, 1);
///   amg->Setup();
///   // Use as preconditioner with CG or other Krylov solver
class CompactAMG : public BaseMatrix {
public:
    CompactAMG(shared_ptr<SparseMatrix<double>> mat,
               shared_ptr<BitArray> freedofs = nullptr,
               double theta = 0.25,
               int max_levels = 25,
               int min_coarse = 50,
               int num_smooth = 1,
               int print_level = 0)
        : mat_(mat), freedofs_(freedofs),
          theta_(theta), max_levels_(max_levels),
          min_coarse_(min_coarse), num_smooth_(num_smooth),
          print_level_(print_level)
    {
    }

    /// Build the AMG hierarchy. Must be called before Mult().
    void Setup() {
        levels_.clear();
        setup_workers_ = 0;

        {
        ngcore::RegionTaskManager tasks;

        // Level 0: finest level
        Level lev0;

        // If freedofs provided, create a modified matrix with identity rows
        // for constrained DOFs (standard AMG treatment of Dirichlet BC).
        // The original NGSolve matrix has full stiffness for all DOFs.
        if (freedofs_) {
            lev0.A = CreateBCModifiedMatrix(*mat_, *freedofs_);
        } else {
            lev0.A = mat_;
        }
        lev0.ndof = lev0.A->Height();
        ComputeL1Norms(*lev0.A, lev0.l1_norms);
        AllocWorkVectors(lev0);
        levels_.push_back(std::move(lev0));

        // Build hierarchy
        for (int l = 0; l < max_levels_ - 1; l++) {
            auto& cur = levels_[l];
            if (cur.ndof <= min_coarse_)
                break;

            // 1. Strength of connection
            auto t_phase = std::chrono::steady_clock::now();
            auto phase = [&](int k) { auto now = std::chrono::steady_clock::now(); setup_phase_s_[k] += std::chrono::duration<double>(now - t_phase).count(); t_phase = now; };
            CSRGraph S = ComputeStrength(*cur.A, theta_);
            phase(0);

            // 2. PMIS coarsening
            std::vector<int> cf_marker = CoarsenPMIS(S);
            phase(1);

            // Count C-points
            int nc = 0;
            for (int i = 0; i < cur.ndof; i++)
                if (cf_marker[i] == 1) nc++;

            if (nc == 0 || nc >= cur.ndof) {
                if (print_level_ > 0) {
                    long long total_strong = 0;
                    int empty_rows = 0;
                    for (int i = 0; i < cur.ndof; i++) {
                        int ns = S.NumNeighbors(i);
                        total_strong += ns;
                        if (ns == 0) empty_rows++;
                    }
                    std::cout << "  CompactAMG level " << l << " coarsening stopped: nc="
                              << nc << "/" << cur.ndof
                              << " strength_nnz=" << total_strong
                              << " empty=" << empty_rows << std::flush;
                }
                break;
            }

            if (print_level_ > 0) {
                std::cout << "  CompactAMG level " << l << ": " << cur.ndof
                          << " -> " << nc << " (ratio "
                          << (double)nc / cur.ndof << ")" << std::endl;
            }

            // 3. Build interpolation
            auto P = BuildClassicalInterp(*cur.A, S, cf_marker, nc);
            phase(2);
            if (!P) break;

            cur.P = P;
            cur.Pt = dynamic_pointer_cast<SparseMatrix<double>>(P->CreateTranspose(true));
            phase(3);

            // 4. Galerkin coarse matrix: A_c = P^T * A * P
            Level next_lev;
            next_lev.A = native_galerkin_ ? GalerkinProduct(*cur.Pt, *cur.A, *P)
                                          : dynamic_pointer_cast<SparseMatrix<double>>(cur.A->Restrict(*P));
            phase(4);
            if (!next_lev.A) break;

            next_lev.ndof = nc;
            ComputeL1Norms(*next_lev.A, next_lev.l1_norms);
            AllocWorkVectors(next_lev);
            phase(5);
            levels_.push_back(std::move(next_lev));
        }

        }
        // Factorization remains outside the internally owned parallel region.
        // Coarsest level: direct solver
        auto& coarsest = levels_.back();
        auto t_factor = std::chrono::steady_clock::now();
        if (coarsest.ndof <= min_coarse_ * 10) {
            // Use sparse Cholesky for coarsest level
            coarsest.A->SetInverseType("sparsecholesky");
            coarsest.inv = coarsest.A->InverseMatrix(shared_ptr<BitArray>(nullptr));
        }
        // else: just use smoother at coarsest level too
        setup_phase_s_[6] += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_factor).count();

        MirrorLevels();
        if (print_level_ > 0)
            std::cout << "\n  CompactAMG: " << levels_.size() << " levels, coarsest = "
                      << levels_.back().ndof << " DOFs" << std::endl;
    }

    /// Refresh the hierarchy for new values on the same sparsity pattern.
    /// Keeps the coarsening and the interpolation operators of the last
    /// Setup() ("frozen interpolation") and recomputes only the Galerkin
    /// coarse matrices P^T A P, the l1 smoother norms and the coarsest
    /// factorization. Falls back to Setup() when there is no hierarchy yet.
    /// When the new matrix has the sparsity pattern of the previous one, the
    /// level matrices are refreshed in place: the coarse patterns are those of
    /// the previous products, so only the numeric Galerkin phase runs.
    void Refresh(shared_ptr<SparseMatrix<double>> mat) {
        if (mat->Height() != mat_->Height() || mat->Width() != mat_->Width())
            throw std::invalid_argument("CompactAMG::Refresh: matrix dimension changed");
        if (levels_.empty()) { mat_ = mat; Setup(); return; }
        {
        ngcore::RegionTaskManager tasks;
        auto& fine = levels_[0];
        // fine.A carries the previous matrix pattern (a BC-modified copy or
        // the previous matrix itself).
        bool in_place = SameSparsityPattern(*mat, *fine.A);
        mat_ = mat;
        if (!freedofs_) fine.A = mat_;
        else if (in_place) FillBCModifiedValues(*mat_, *freedofs_, *fine.A);
        else fine.A = CreateBCModifiedMatrix(*mat_, *freedofs_);
        ComputeL1Norms(*fine.A, fine.l1_norms);
        for (size_t l = 0; l + 1 < levels_.size(); l++) {
            auto& cur = levels_[l];
            if (!cur.P || !cur.Pt)
                throw std::runtime_error("CompactAMG::Refresh: level without interpolation");
            auto& next = levels_[l + 1];
            in_place = in_place && GalerkinProductOnPattern(*cur.Pt, *cur.A, *cur.P, *next.A);
            if (!in_place) {
                auto coarse = native_galerkin_ ? GalerkinProduct(*cur.Pt, *cur.A, *cur.P)
                    : dynamic_pointer_cast<SparseMatrix<double>>(cur.A->Restrict(*cur.P));
                if (!coarse || coarse->Height() != next.ndof)
                    throw std::runtime_error("CompactAMG::Refresh: Galerkin product changed size");
                next.A = coarse;
            } else {
                in_place_products_++;
            }
            ComputeL1Norms(*next.A, next.l1_norms);
        }
        }
        auto& coarsest = levels_.back();
        if (coarsest.inv) {
            coarsest.A->SetInverseType("sparsecholesky");
            coarsest.inv = coarsest.A->InverseMatrix(shared_ptr<BitArray>(nullptr));
        }
        MirrorLevels();
        refresh_count_++;
    }

    int RefreshCount() const { return refresh_count_; }
    /// Build coarse matrices with the native symbolic + numeric Galerkin product
    /// instead of SparseMatrix::Restrict (same matrices up to summation order).
    void SetNativeGalerkin(bool on) { native_galerkin_ = on; }
    /// V-cycle residuals read float32 mirrors of the level matrices (set before Setup).
    void SetMixedPrecision(bool on) { mixed_precision_ = on; }
    /// Accumulated Setup phases: strength, coarsening, interpolation,
    /// transpose, Galerkin, l1 norms + work vectors, coarsest factorization.
    const std::array<double, 7>& SetupPhases() const { return setup_phase_s_; }
    /// Galerkin coarse matrices refreshed numerically on their existing pattern.
    int InPlaceProductCount() const { return in_place_products_; }

    // BaseMatrix interface
    int VHeight() const override { return mat_->Height(); }
    int VWidth() const override { return mat_->Width(); }
    bool IsComplex() const override { return false; }

    AutoVector CreateRowVector() const override { return mat_->CreateRowVector(); }
    AutoVector CreateColVector() const override { return mat_->CreateColVector(); }

    void Mult(const BaseVector& b, BaseVector& x) const override {
        x = 0;
        VCycle(0, b, x);

        // Zero constrained DOFs in output
        if (freedofs_) {
            auto fv = x.FVDouble();
            int n = (int)fv.Size();
            ParallelFor(n, [&](size_t i) {
                if (!freedofs_->Test(i))
                    fv[i] = 0;
            });
        }
    }

    void MultTrans(const BaseVector& b, BaseVector& x) const override {
        Mult(b, x);  // Symmetric preconditioner with l1-GS (forward pre, backward post)
    }

    /// Dual Mult: apply AMG V-cycle to two RHS simultaneously with fused SpMV.
    /// Every SpMV at every level loads matrix rows once for both RHS (halves bandwidth).
    /// Used by ComplexHypreBasedAMS for fused Re/Im processing.
    void DualMult(const BaseVector& b1, BaseVector& x1,
                  const BaseVector& b2, BaseVector& x2) const {
        x1 = 0;
        x2 = 0;
        DualVCycle(0, b1, x1, b2, x2);

        if (freedofs_) {
            auto fv1 = x1.FVDouble();
            auto fv2 = x2.FVDouble();
            int n = (int)fv1.Size();
            ParallelFor(n, [&](size_t i) {
                if (!freedofs_->Test(i)) {
                    fv1[i] = 0;
                    fv2[i] = 0;
                }
            });
        }
    }

    int NumLevels() const { return (int)levels_.size(); }
    int SetupWorkers() const { return setup_workers_; }

private:
    mutable int setup_workers_ = 0;
    int refresh_count_ = 0;
    bool native_galerkin_ = false;
    bool mixed_precision_ = false;
    std::array<double, 7> setup_phase_s_{};
    int in_place_products_ = 0;
    // =====================================================================
    // Level data
    // =====================================================================
    struct Level {
        shared_ptr<SparseMatrix<double>> A;
        shared_ptr<SparseMatrix<double>> P;   // Prolongation to this level
        shared_ptr<SparseMatrix<double>> Pt;  // P^T (restriction)
        shared_ptr<BaseMatrix> inv;           // Direct solver (coarsest only)
        int ndof = 0;
        std::vector<double> l1_norms;
        std::vector<float> values_f;  // float32 mirror of A (mixed precision only)

        // Work vectors (mutable for const Mult)
        mutable std::unique_ptr<VVector<double>> residual;
        mutable std::unique_ptr<VVector<double>> tmp;
        mutable std::unique_ptr<VVector<double>> correction;

        // Dual work vectors for fused Re/Im (DualMult)
        mutable std::unique_ptr<VVector<double>> residual2;
        mutable std::unique_ptr<VVector<double>> tmp2;
        mutable std::unique_ptr<VVector<double>> correction2;
    };

    shared_ptr<SparseMatrix<double>> mat_;
    shared_ptr<BitArray> freedofs_;
    double theta_;
    int max_levels_;
    int min_coarse_;
    int num_smooth_;
    int print_level_;
    std::vector<Level> levels_;

    void AllocWorkVectors(Level& lev) {
        lev.residual = std::make_unique<VVector<double>>(lev.ndof);
        lev.tmp = std::make_unique<VVector<double>>(lev.ndof);
        lev.correction = std::make_unique<VVector<double>>(lev.ndof);
        lev.residual2 = std::make_unique<VVector<double>>(lev.ndof);
        lev.tmp2 = std::make_unique<VVector<double>>(lev.ndof);
        lev.correction2 = std::make_unique<VVector<double>>(lev.ndof);
    }

    /// Create a copy of A with identity rows for constrained DOFs.
    /// Standard AMG treatment of Dirichlet BC: constrained DOFs are decoupled.
    shared_ptr<SparseMatrix<double>> CreateBCModifiedMatrix(
        const SparseMatrix<double>& A, const BitArray& freedofs) const
    {
        int n = A.Height();

        // Copy sparsity structure (same nnz per row)
        Array<int> cnt(n);
        for (int i = 0; i < n; i++)
            cnt[i] = A.GetRowIndices(i).Size();

        auto B = make_shared<SparseMatrix<double>>(cnt, A.Width());

        ParallelFor(n, [&](size_t i) {
            auto src_cols = A.GetRowIndices(i);
            auto src_vals = A.GetRowValues(i);
            auto dst_cols = B->GetRowIndices(i);
            auto dst_vals = B->GetRowValues(i);

            for (int j = 0; j < src_cols.Size(); j++) {
                dst_cols[j] = src_cols[j];

                if (!freedofs.Test(i)) {
                    // Constrained row: identity
                    dst_vals[j] = (src_cols[j] == (int)i) ? 1.0 : 0.0;
                } else if (!freedofs.Test(src_cols[j])) {
                    // Free row, but column is constrained: zero coupling
                    dst_vals[j] = 0.0;
                } else {
                    // Free row, free column: keep original
                    dst_vals[j] = src_vals[j];
                }
            }
        });

        return B;
    }

    // =====================================================================
    // Strength of Connection
    // =====================================================================
    CSRGraph ComputeStrength(const SparseMatrix<double>& A, double theta) const {
        int n = A.Height();
        CSRGraph S;
        S.n = n;
        S.row_ptr.resize(n + 1);

        // Pass 1: count strong connections per row
        // Standard AMG strength for M-matrices / HCurl systems:
        //   For each row i, find max negative off-diagonal: max_neg = max_{j!=i}(-a_ij)
        //   Strong connection: -a_ij >= theta * max_neg
        //   If no negative off-diagonals (rare), fallback to absolute value criterion.
        std::vector<int> row_count(n, 0);
        std::vector<int> workers(ngcore::TaskManager::GetNumThreads(), 0);

        ParallelForRange(n, [&](IntRange range) {
            if (range.Size()) workers[ngcore::TaskManager::GetThreadId()] = 1;
            for (auto i : range) {
            auto cols = A.GetRowIndices(i);
            auto vals = A.GetRowValues(i);

            // Find max negative off-diagonal (standard AMG)
            double max_neg = 0;
            double max_abs = 0;
            for (int j = 0; j < cols.Size(); j++) {
                if (cols[j] != (int)i) {
                    double v = vals[j];
                    if (-v > max_neg) max_neg = -v;  // -a_ij for negative entries
                    double av = std::abs(v);
                    if (av > max_abs) max_abs = av;
                }
            }

            // Use max_neg for M-matrices; fallback to max_abs if all positive
            double ref = (max_neg > 0) ? max_neg : max_abs;
            double threshold = theta * ref;

            int cnt = 0;
            for (int j = 0; j < cols.Size(); j++) {
                if (cols[j] != (int)i) {
                    // Strong if -a_ij >= threshold (neg entries) OR |a_ij| >= threshold (general)
                    bool strong = (max_neg > 0) ? (-vals[j] >= threshold) : (std::abs(vals[j]) >= threshold);
                    if (strong) cnt++;
                }
            }
            row_count[i] = cnt;
            }
        });
        setup_workers_ = std::max(setup_workers_,
            (int)std::count(workers.begin(), workers.end(), 1));

        // Build row_ptr
        S.row_ptr[0] = 0;
        for (int i = 0; i < n; i++)
            S.row_ptr[i + 1] = S.row_ptr[i] + row_count[i];

        // Pass 2: fill col_idx
        S.col_idx.resize(S.row_ptr[n]);

        ParallelFor(n, [&](size_t i) {
            auto cols = A.GetRowIndices(i);
            auto vals = A.GetRowValues(i);

            double max_neg = 0;
            double max_abs = 0;
            for (int j = 0; j < cols.Size(); j++) {
                if (cols[j] != (int)i) {
                    double v = vals[j];
                    if (-v > max_neg) max_neg = -v;
                    double av = std::abs(v);
                    if (av > max_abs) max_abs = av;
                }
            }
            double ref = (max_neg > 0) ? max_neg : max_abs;
            double threshold = theta * ref;

            int pos = S.row_ptr[i];
            for (int j = 0; j < cols.Size(); j++) {
                if (cols[j] != (int)i) {
                    bool strong = (max_neg > 0) ? (-vals[j] >= threshold) : (std::abs(vals[j]) >= threshold);
                    if (strong)
                        S.col_idx[pos++] = cols[j];
                }
            }
        });

        return S;
    }

    // =====================================================================
    // PMIS Coarsening (Parallel Maximal Independent Set)
    // =====================================================================
    // Returns cf_marker: 1 = C-point, 0 = F-point
    std::vector<int> CoarsenPMIS(const CSRGraph& S) const {
        int n = S.n;
        std::vector<int> cf(n, -1);  // -1 = undecided

        // Compute S^T (transpose strength graph) for measure
        CSRGraph St = TransposeGraph(S);

        // Measure = |S^T_i| (number of vertices that strongly depend on i)
        // Plus random tiebreak in [0, 1)
        std::vector<double> measure(n);
        std::mt19937 rng(42);  // Deterministic
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (int i = 0; i < n; i++) {
            measure[i] = (double)St.NumNeighbors(i) + dist(rng);
        }

        // PMIS iterations
        for (int iter = 0; iter < 100; iter++) {
            bool changed = false;

            // Phase 1: Mark new C-points (local maxima among undecided)
            std::vector<int> new_c;
            for (int i = 0; i < n; i++) {
                if (cf[i] != -1) continue;  // already decided

                bool is_max = true;
                for (const int* p = S.NeighborBegin(i); p != S.NeighborEnd(i); ++p) {
                    int j = *p;
                    if (cf[j] == -1 && measure[j] > measure[i]) {
                        is_max = false;
                        break;
                    }
                }
                // Also check transpose neighbors
                if (is_max) {
                    for (const int* p = St.NeighborBegin(i); p != St.NeighborEnd(i); ++p) {
                        int j = *p;
                        if (cf[j] == -1 && measure[j] > measure[i]) {
                            is_max = false;
                            break;
                        }
                    }
                }

                if (is_max) {
                    cf[i] = 1;  // Mark C-point immediately
                    new_c.push_back(i);
                    changed = true;
                }
            }

            // Phase 2: Mark new F-points (undecided neighbors of C-points)
            for (int i : new_c) {
                for (const int* p = S.NeighborBegin(i); p != S.NeighborEnd(i); ++p) {
                    if (cf[*p] == -1) {
                        cf[*p] = 0;  // F-point
                        changed = true;
                    }
                }
                for (const int* p = St.NeighborBegin(i); p != St.NeighborEnd(i); ++p) {
                    if (cf[*p] == -1) {
                        cf[*p] = 0;  // F-point
                        changed = true;
                    }
                }
            }

            if (!changed) break;
        }

        // Any remaining undecided → C-point
        for (int i = 0; i < n; i++)
            if (cf[i] == -1) cf[i] = 1;

        return cf;
    }

    /// Transpose a CSR graph
    CSRGraph TransposeGraph(const CSRGraph& G) const {
        int n = G.n;
        CSRGraph Gt;
        Gt.n = n;
        Gt.row_ptr.resize(n + 1, 0);

        // Count entries per row in transpose
        for (int i = 0; i < n; i++)
            for (const int* p = G.NeighborBegin(i); p != G.NeighborEnd(i); ++p)
                Gt.row_ptr[*p + 1]++;

        // Prefix sum
        for (int i = 0; i < n; i++)
            Gt.row_ptr[i + 1] += Gt.row_ptr[i];

        // Fill
        Gt.col_idx.resize(Gt.row_ptr[n]);
        std::vector<int> pos(n, 0);
        for (int i = 0; i < n; i++)
            for (const int* p = G.NeighborBegin(i); p != G.NeighborEnd(i); ++p) {
                int j = *p;
                Gt.col_idx[Gt.row_ptr[j] + pos[j]] = i;
                pos[j]++;
            }

        return Gt;
    }

    // =====================================================================
    // Classical Direct Interpolation
    // =====================================================================
    shared_ptr<SparseMatrix<double>> BuildClassicalInterp(
        const SparseMatrix<double>& A,
        const CSRGraph& S,
        const std::vector<int>& cf_marker,
        int nc) const
    {
        int nf = A.Height();

        // Build C-point -> coarse index map
        std::vector<int> coarse_idx(nf, -1);
        int cidx = 0;
        for (int i = 0; i < nf; i++)
            if (cf_marker[i] == 1)
                coarse_idx[i] = cidx++;

        // Strength flags aligned with A's value positions (flat, filled in parallel)
        std::vector<char> strong_flag(A.NZE(), 0);
        ParallelFor(nf, [&](size_t i) {
            auto cols = A.GetRowIndices(i);
            const size_t base = A.First(i);
            const int spos = S.row_ptr[i];
            const int send = S.row_ptr[i + 1];
            for (int j = 0; j < cols.Size(); j++) {
                for (int k = spos; k < send; k++) {
                    if (S.col_idx[k] == cols[j]) {
                        strong_flag[base + j] = 1;
                        break;
                    }
                }
            }
        });
        auto is_strong_at = [&](int i, int j) { return strong_flag[A.First(i) + j] != 0; };

        // Count entries per row in P
        std::vector<int> P_row_nnz(nf, 0);
        ParallelFor(nf, [&](size_t i) {
            if (cf_marker[i] == 1) {
                P_row_nnz[i] = 1;  // C-point: identity
            } else {
                // F-point: count strong C-neighbors
                auto cols = A.GetRowIndices(i);
                int count = 0;
                for (int j = 0; j < cols.Size(); j++) {
                    if (is_strong_at(int(i), j) && cf_marker[cols[j]] == 1)
                        count++;
                }
                P_row_nnz[i] = count == 0 ? 1 : count;  // Fallback: inject to nearest C
            }
        });

        // Build P as NGSolve SparseMatrix
        // First create the sparsity pattern using Table
        Array<int> cnt(nf);
        for (int i = 0; i < nf; i++)
            cnt[i] = P_row_nnz[i];

        auto P = make_shared<SparseMatrix<double>>(cnt, nc);

        // Fill P values
        ParallelFor(nf, [&](size_t i) {
            if (cf_marker[i] == 1) {
                // C-point: P(i, coarse_idx[i]) = 1
                P->GetRowIndices(i)[0] = coarse_idx[i];
                P->GetRowValues(i)[0] = 1.0;
            } else {
                // F-point: classical direct interpolation
                auto cols = A.GetRowIndices(i);
                auto vals = A.GetRowValues(i);

                // Get diagonal
                double a_ii = 0;
                for (int j = 0; j < cols.Size(); j++)
                    if (cols[j] == (int)i) { a_ii = vals[j]; break; }

                // Sum of non-interpolated connections (lumped into diagonal):
                // = weak connections + strong F-point connections
                // Only strong C-neighbors become P entries; everything else
                // is lumped into the diagonal for row-sum preservation.
                double sum_non_interp = 0;
                for (int j = 0; j < cols.Size(); j++) {
                    if (cols[j] == (int)i) continue;  // skip diagonal
                    // Skip strong C-neighbors (these become interpolation weights)
                    if (is_strong_at(int(i), j) && cf_marker[cols[j]] == 1) continue;
                    sum_non_interp += vals[j];
                }

                double denom = a_ii + sum_non_interp;
                if (std::abs(denom) < 1e-30) denom = 1.0;  // Safety

                // Fill P row: w_ij = -a_ij / denom for strong C-neighbors
                int pos = 0;
                auto p_cols = P->GetRowIndices(i);
                auto p_vals = P->GetRowValues(i);

                for (int j = 0; j < cols.Size(); j++) {
                    if (is_strong_at(int(i), j) && cf_marker[cols[j]] == 1) {
                        p_cols[pos] = coarse_idx[cols[j]];
                        p_vals[pos] = -vals[j] / denom;
                        pos++;
                    }
                }

                // Fallback: if no strong C-neighbors, inject to nearest C
                if (pos == 0) {
                    // Find nearest C-point in matrix row
                    int best_c = -1;
                    double best_val = 0;
                    for (int j = 0; j < cols.Size(); j++) {
                        if (cf_marker[cols[j]] == 1 && std::abs(vals[j]) > best_val) {
                            best_val = std::abs(vals[j]);
                            best_c = coarse_idx[cols[j]];
                        }
                    }
                    if (best_c >= 0) {
                        p_cols[0] = best_c;
                        p_vals[0] = 1.0;
                    } else {
                        p_cols[0] = 0;
                        p_vals[0] = 0.0;
                    }
                }
            }
        });

        return P;
    }

    // =====================================================================
    // l1-Jacobi Smoother (fully TaskManager parallel)
    //
    // x += D_l1^{-1} * (b - A*x) per sweep.
    // No data dependency between rows -> full ParallelFor parallelism.
    // Reference: Baker et al., SIAM J. Sci. Comput. 33(5), 2011.
    // =====================================================================
    void ComputeL1Norms(const SparseMatrix<double>& A,
                        std::vector<double>& norms) const {
        int n = A.Height();
        norms.resize(n);
        ParallelFor(n, [&](size_t i) {
            auto vals = A.GetRowValues(i);
            double sum = 0;
            for (int j = 0; j < vals.Size(); j++)
                sum += std::abs(vals[j]);
            norms[i] = (sum > 0) ? sum : 1.0;
        });
    }

    /// Copy BaseVector data via FlatVector
    static void CopyVector(const BaseVector& src, BaseVector& dst) {
        auto fv_s = src.FVDouble();
        auto fv_d = dst.FVDouble();
        ParallelFor(fv_s.Size(), [&](size_t i) { fv_d[i] = fv_s[i]; });
    }

    /// res = b - A x in one pass over the rows (no separate copy of b).
    static void ResidualInto(const SparseMatrix<double>& A, const BaseVector& b,
                             const BaseVector& x, BaseVector& res) {
        ResidualWithValues(A, nullptr, b, x, res);
    }

    /// Level residual; with mixed precision A's values come from the level's
    /// float32 mirror (vectors and accumulation stay double).
    void LevelResidual(const Level& lev, const BaseVector& b, const BaseVector& x,
                       BaseVector& res) const {
        ResidualWithValues(*lev.A, mixed_precision_ ? lev.values_f.data() : nullptr, b, x, res);
    }

    void MirrorLevels() {
        for (auto& lev : levels_) {
            if (mixed_precision_) MirrorValues(*lev.A, lev.values_f);
            else std::vector<float>().swap(lev.values_f);
        }
    }

    /// l1-Jacobi sweep: x += r / l1_norm (fully parallel, no data dependency)
    void L1JacobiSmooth(int level, const BaseVector& b, BaseVector& x,
                        bool initially_zero = false) const {
        auto& lev = levels_[level];
        if (initially_zero) {
            // Every real V-cycle starts from zero, including coarse corrections.
            auto fv_b = b.FVDouble();
            auto fv_x = x.FVDouble();
            ParallelFor(lev.ndof, [&](size_t i) {
                fv_x[i] = fv_b[i] / lev.l1_norms[i];
            });
            return;
        }
        auto& res = *lev.residual;

        // Residual: r = b - A*x in one pass
        LevelResidual(lev, b, x, res);

        // Jacobi update: x[i] += r[i] / l1_norm[i] (fully parallel)
        auto fv_x = x.FVDouble();
        auto fv_r = res.FVDouble();
        int n = lev.ndof;
        ParallelFor(n, [&](size_t i) {
            fv_x[i] += fv_r[i] / lev.l1_norms[i];
        });
    }

    // =====================================================================
    // Dual V-Cycle (fused Re/Im at every level)
    // =====================================================================
    void DualVCycle(int level, const BaseVector& b1, BaseVector& x1,
                    const BaseVector& b2, BaseVector& x2) const {
        auto& lev = levels_[level];

        if (level == (int)levels_.size() - 1) {
            // Coarsest level: direct solve (sequential, tiny problem)
            if (lev.inv) {
                lev.inv->Mult(b1, x1);
                lev.inv->Mult(b2, x2);
            } else {
                for (int s = 0; s < 10; s++)
                    DualL1JacobiSmooth(level, b1, x1, b2, x2);
            }
            return;
        }

        // Pre-smooth: fused l1-Jacobi
        for (int s = 0; s < num_smooth_; s++)
            DualL1JacobiSmooth(level, b1, x1, b2, x2);

        // Fused residual
        auto& res1 = *lev.residual;
        auto& res2 = *lev.residual2;
        DualResidual(*lev.A, b1, x1, res1, b2, x2, res2, lev.ndof);

        // Fused restrict
        auto& next = levels_[level + 1];
        auto& rc1 = *next.tmp;
        auto& rc2 = *next.tmp2;
        DualSpMV(*lev.Pt, res1, rc1, res2, rc2, lev.Pt->Height());

        // Recursive dual V-cycle
        auto& ec1 = *next.correction;
        auto& ec2 = *next.correction2;
        ec1 = 0;
        ec2 = 0;
        DualVCycle(level + 1, rc1, ec1, rc2, ec2);

        // Fused prolongate: x += P * e_c
        DualMultAdd(*lev.P, ec1, x1, ec2, x2, lev.P->Height());

        // Post-smooth: fused l1-Jacobi
        for (int s = 0; s < num_smooth_; s++)
            DualL1JacobiSmooth(level, b1, x1, b2, x2);
    }

    /// Fused l1-Jacobi smooth for two RHS: two-phase (residual then update)
    void DualL1JacobiSmooth(int level, const BaseVector& b1, BaseVector& x1,
                            const BaseVector& b2, BaseVector& x2) const {
        auto& lev = levels_[level];
        auto& res1 = *lev.residual;
        auto& res2 = *lev.residual2;

        // Phase 1: fused residual
        DualResidual(*lev.A, b1, x1, res1, b2, x2, res2, lev.ndof);

        // Phase 2: fused Jacobi update
        auto fv_x1 = x1.FVDouble();
        auto fv_x2 = x2.FVDouble();
        auto fv_r1 = res1.FVDouble();
        auto fv_r2 = res2.FVDouble();
        int n = lev.ndof;
        ParallelFor(n, [&](size_t i) {
            double inv_l1 = 1.0 / lev.l1_norms[i];
            fv_x1[i] += fv_r1[i] * inv_l1;
            fv_x2[i] += fv_r2[i] * inv_l1;
        });
    }

    /// Fused residual: res = b - A*x for two RHS in single matrix pass
    static void DualResidual(const SparseMatrix<double>& A,
                             const BaseVector& b1, const BaseVector& x1, BaseVector& res1,
                             const BaseVector& b2, const BaseVector& x2, BaseVector& res2,
                             int n) {
        auto fv_b1 = b1.FVDouble(); auto fv_x1 = x1.FVDouble(); auto fv_r1 = res1.FVDouble();
        auto fv_b2 = b2.FVDouble(); auto fv_x2 = x2.FVDouble(); auto fv_r2 = res2.FVDouble();

        ParallelFor(n, [&](size_t i) {
            auto cols = A.GetRowIndices(i);
            auto vals = A.GetRowValues(i);
            double d1 = 0, d2 = 0;
            for (int j = 0; j < cols.Size(); j++) {
                int c = cols[j];
                double v = vals[j];
                d1 += v * fv_x1[c];
                d2 += v * fv_x2[c];
            }
            fv_r1[i] = fv_b1[i] - d1;
            fv_r2[i] = fv_b2[i] - d2;
        });
    }

    /// Fused SpMV: y = A*x for two RHS in single matrix pass
    static void DualSpMV(const SparseMatrix<double>& A,
                         const BaseVector& x1, BaseVector& y1,
                         const BaseVector& x2, BaseVector& y2, int n) {
        auto fv_x1 = x1.FVDouble(); auto fv_y1 = y1.FVDouble();
        auto fv_x2 = x2.FVDouble(); auto fv_y2 = y2.FVDouble();

        ParallelFor(n, [&](size_t i) {
            auto cols = A.GetRowIndices(i);
            auto vals = A.GetRowValues(i);
            double s1 = 0, s2 = 0;
            for (int j = 0; j < cols.Size(); j++) {
                int c = cols[j];
                double v = vals[j];
                s1 += v * fv_x1[c];
                s2 += v * fv_x2[c];
            }
            fv_y1[i] = s1;
            fv_y2[i] = s2;
        });
    }

    /// Fused MultAdd: x += P*g for two RHS in single matrix pass
    static void DualMultAdd(const SparseMatrix<double>& P,
                            const BaseVector& g1, BaseVector& x1,
                            const BaseVector& g2, BaseVector& x2, int nrows) {
        auto fv_g1 = g1.FVDouble(); auto fv_x1 = x1.FVDouble();
        auto fv_g2 = g2.FVDouble(); auto fv_x2 = x2.FVDouble();

        ParallelFor(nrows, [&](size_t i) {
            auto cols = P.GetRowIndices(i);
            auto vals = P.GetRowValues(i);
            double s1 = 0, s2 = 0;
            for (int j = 0; j < cols.Size(); j++) {
                int c = cols[j];
                double v = vals[j];
                s1 += v * fv_g1[c];
                s2 += v * fv_g2[c];
            }
            fv_x1[i] += s1;
            fv_x2[i] += s2;
        });
    }

    // =====================================================================
    // V-Cycle (single RHS, original)
    // =====================================================================
    void VCycle(int level, const BaseVector& b, BaseVector& x) const {
        auto& lev = levels_[level];

        // Coarsest level: direct solve or just smooth
        if (level == (int)levels_.size() - 1) {
            if (lev.inv) {
                lev.inv->Mult(b, x);
            } else {
                // Fall back to l1-Jacobi smoothing
                for (int s = 0; s < 10; s++)
                    L1JacobiSmooth(level, b, x, s == 0);
            }
            return;
        }

        // Pre-smooth: l1-Jacobi (fully parallel)
        for (int s = 0; s < num_smooth_; s++)
            L1JacobiSmooth(level, b, x, s == 0);

        // Compute residual: r = b - A*x
        auto& res = *lev.residual;
        LevelResidual(lev, b, x, res);

        // Restrict to coarse: r_c = P^T * r
        auto& next = levels_[level + 1];
        auto& r_c = *next.tmp;
        lev.Pt->Mult(res, r_c);

        // Coarse solve: e_c = 0; VCycle(level+1, r_c, e_c)
        auto& e_c = *next.correction;
        // With a pre-smoother the coarse V-cycle assigns every entry (first
        // smoother or direct solve); zero it only when there is none.
        if (num_smooth_ < 1) e_c = 0;
        VCycle(level + 1, r_c, e_c);

        // Prolongate and add: x += P * e_c
        lev.P->MultAdd(1.0, e_c, x);

        // Post-smooth: l1-Jacobi (symmetric V-cycle)
        for (int s = 0; s < num_smooth_; s++)
            L1JacobiSmooth(level, b, x);
    }
};

}  // namespace ngla

#endif  // SPARSESOLV_COMPACT_AMG_HPP
