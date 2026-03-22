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

namespace ngla {

/// Compact CSR graph (binary adjacency, no values)
struct CSRGraph {
    int n = 0;
    std::vector<int> row_ptr;
    std::vector<int> col_idx;

    int NumNeighbors(int i) const { return row_ptr[i + 1] - row_ptr[i]; }

    const int* NeighborBegin(int i) const { return col_idx.data() + row_ptr[i]; }
    const int* NeighborEnd(int i) const { return col_idx.data() + row_ptr[i + 1]; }
};

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
            CSRGraph S = ComputeStrength(*cur.A, theta_);

            // 2. PMIS coarsening
            std::vector<int> cf_marker = CoarsenPMIS(S);

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
            if (!P) break;

            cur.P = P;
            cur.Pt = dynamic_pointer_cast<SparseMatrix<double>>(P->CreateTranspose(true));

            // 4. Galerkin coarse matrix: A_c = P^T * A * P
            Level next_lev;
            next_lev.A = dynamic_pointer_cast<SparseMatrix<double>>(cur.A->Restrict(*P));

            if (!next_lev.A) break;

            next_lev.ndof = nc;
            ComputeL1Norms(*next_lev.A, next_lev.l1_norms);
            AllocWorkVectors(next_lev);
            levels_.push_back(std::move(next_lev));
        }

        // Coarsest level: direct solver
        auto& coarsest = levels_.back();
        if (coarsest.ndof <= min_coarse_ * 10) {
            // Use sparse Cholesky for coarsest level
            coarsest.A->SetInverseType(SPARSECHOLESKY);
            coarsest.inv = coarsest.A->InverseMatrix(shared_ptr<BitArray>(nullptr));
        }
        // else: just use smoother at coarsest level too

        if (print_level_ > 0)
            std::cout << "\n  CompactAMG: " << levels_.size() << " levels, coarsest = "
                      << levels_.back().ndof << " DOFs" << std::endl;
    }

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
    /// Used by ComplexCompactAMS for fused Re/Im processing.
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

private:
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

        ParallelFor(n, [&](size_t i) {
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
        });

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

        // Build strength set for fast lookup
        std::vector<std::vector<bool>> is_strong(nf);
        for (int i = 0; i < nf; i++) {
            auto cols = A.GetRowIndices(i);
            is_strong[i].resize(cols.Size(), false);
            int spos = S.row_ptr[i];
            int send = S.row_ptr[i + 1];
            for (int j = 0; j < cols.Size(); j++) {
                for (int k = spos; k < send; k++) {
                    if (S.col_idx[k] == cols[j]) {
                        is_strong[i][j] = true;
                        break;
                    }
                }
            }
        }

        // Count entries per row in P
        std::vector<int> P_row_nnz(nf, 0);
        for (int i = 0; i < nf; i++) {
            if (cf_marker[i] == 1) {
                P_row_nnz[i] = 1;  // C-point: identity
            } else {
                // F-point: count strong C-neighbors
                auto cols = A.GetRowIndices(i);
                for (int j = 0; j < cols.Size(); j++) {
                    if (is_strong[i][j] && cf_marker[cols[j]] == 1)
                        P_row_nnz[i]++;
                }
                if (P_row_nnz[i] == 0)
                    P_row_nnz[i] = 1;  // Fallback: inject to nearest C
            }
        }

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
                    if (is_strong[i][j] && cf_marker[cols[j]] == 1) continue;
                    sum_non_interp += vals[j];
                }

                double denom = a_ii + sum_non_interp;
                if (std::abs(denom) < 1e-30) denom = 1.0;  // Safety

                // Fill P row: w_ij = -a_ij / denom for strong C-neighbors
                int pos = 0;
                auto p_cols = P->GetRowIndices(i);
                auto p_vals = P->GetRowValues(i);

                for (int j = 0; j < cols.Size(); j++) {
                    if (is_strong[i][j] && cf_marker[cols[j]] == 1) {
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

    /// l1-Jacobi sweep: x += r / l1_norm (fully parallel, no data dependency)
    void L1JacobiSmooth(int level, const BaseVector& b, BaseVector& x) const {
        auto& lev = levels_[level];
        auto& res = *lev.residual;

        // Residual: r = b - A*x (NGSolve SpMV is TaskManager-parallel)
        CopyVector(b, res);
        lev.A->MultAdd(-1.0, x, res);

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
                    L1JacobiSmooth(level, b, x);
            }
            return;
        }

        // Pre-smooth: l1-Jacobi (fully parallel)
        for (int s = 0; s < num_smooth_; s++)
            L1JacobiSmooth(level, b, x);

        // Compute residual: r = b - A*x
        auto& res = *lev.residual;
        CopyVector(b, res);
        lev.A->MultAdd(-1.0, x, res);

        // Restrict to coarse: r_c = P^T * r
        auto& next = levels_[level + 1];
        auto& r_c = *next.tmp;
        lev.Pt->Mult(res, r_c);

        // Coarse solve: e_c = 0; VCycle(level+1, r_c, e_c)
        auto& e_c = *next.correction;
        e_c = 0;
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
