/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/// @file ams_wirebasket.hpp
/// @brief Compact AMS as the wirebasket (coarse) solver of NGSolve's BDDC.
///
/// For HCurl(order >= 1) the BDDC wirebasket holds every lowest-order edge
/// dof (dof e is the Whitney function of edge e) and, on badly shaped faces,
/// some face dofs as well. The wirebasket Schur complement restricted to the
/// edge block is a lowest-order Nedelec operator whose kernel still holds the
/// vertex gradients, so an AMS built on the mesh vertices can replace its
/// direct factorization. The extra face dofs enter the same AMS as rows with
/// an empty discrete gradient: its fine l1-Jacobi smoother covers them, the
/// auxiliary spaces do not. BDDC selects this through its own flag,
/// `coarsetype="sparsesolv_ams"`, once the module is imported.
///
/// A complex wirebasket matrix S (A = K + i w sigma M) is preconditioned by
/// the real AMS of Re S + Im S, the wirebasket counterpart of the usual
/// K + |w| sigma M surrogate: rotating A by exp(-i pi/4) keeps its Hermitian
/// part positive, Schur complements keep that property, so Re S + Im S is
/// positive semidefinite. `eps` adds eps * diag to the surrogate only.
/// `cycles` = k runs k stationary steps y += B (x - S y) on S itself.

#ifndef SPARSESOLV_AMS_WIREBASKET_HPP
#define SPARSESOLV_AMS_WIREBASKET_HPP

#include <comp.hpp>
#include <chrono>
#include <algorithm>
#include "sparsesolv/preconditioners/compact_ams.hpp"
#include "sparsesolv/preconditioners/complex_compact_ams.hpp"

namespace ngla {

/// Discrete gradient of the vertex space into the lowest-order edge dofs:
/// row e holds -1 / +1 at the lower / higher vertex number of edge e, the
/// global orientation NGSolve's HCurl spaces use. Rows past the edges (extra
/// dofs) stay empty.
inline shared_ptr<SparseMatrix<double>> BuildLowestOrderGradient(const ngcomp::MeshAccess& ma,
                                                                 size_t rows = 0) {
    const size_t nedge = ma.GetNEdges();
    rows = std::max(rows, nedge);
    Array<int> counts(rows);
    counts = 0;
    counts.Range(0, nedge) = 2;
    auto gradient = make_shared<SparseMatrix<double>>(counts, int(ma.GetNV()));
    ParallelFor(nedge, [&](size_t e) {
        auto pnums = ma.GetEdgePNums(e);
        const int first = pnums[0], second = pnums[1];
        auto columns = gradient->GetRowIndices(e);
        auto values = gradient->GetRowValues(e);
        if (first < second) {
            columns[0] = first;  values[0] = -1.0;
            columns[1] = second; values[1] = 1.0;
        } else {
            columns[0] = second; values[0] = 1.0;
            columns[1] = first;  values[1] = -1.0;
        }
    });
    return gradient;
}

/// Counters of the most recently built wirebasket AMS, read from Python.
struct AMSWirebasketStats {
    int builds = 0;
    int n_edges = 0;
    int n_extra = 0;
    int n_free = 0;
    int n_vertices = 0;
    int cycles = 1;
    bool complex = false;
    double extract_s = 0;
    double setup_s = 0;
    long applies = 0;
    double apply_s = 0;
    std::weak_ptr<ComplexHypreBasedAMS> complex_ams;  // stage timers of the complex cycle
};
inline AMSWirebasketStats& GetAMSWirebasketStats() {
    static AMSWirebasketStats stats;
    return stats;
}

/// Full-size operator: gathers the wirebasket block (edges, then the extra
/// dofs), applies the AMS there and zeroes every other dof, as the direct
/// wirebasket inverse does. With cycles = k > 1 it runs k steps of the
/// stationary iteration y += B (x - S y) on the wirebasket system S (B = one
/// AMS cycle): a fixed polynomial in B S times B, so the operator stays
/// (complex) symmetric and needs no flexible outer Krylov method.
template <typename SCAL>
class AMSWirebasketOperator : public BaseMatrix {
public:
    AMSWirebasketOperator(shared_ptr<BaseMatrix> ams, shared_ptr<SparseMatrix<SCAL>> system,
                          shared_ptr<BitArray> free, std::vector<int> extra,
                          size_t ndof, size_t nedge, int cycles)
        : ams_(ams), system_(system), free_(free), extra_(std::move(extra)),
          ndof_(ndof), nedge_(nedge), nloc_(nedge + extra_.size()),
          cycles_(std::max(cycles, 1)),
          xl_(nloc_), yl_(nloc_), res_(nloc_), cor_(nloc_) {}

    int VHeight() const override { return int(ndof_); }
    int VWidth() const override { return int(ndof_); }
    bool IsComplex() const override { return std::is_same_v<SCAL, Complex>; }
    AutoVector CreateRowVector() const override { return make_unique<VVector<SCAL>>(ndof_); }
    AutoVector CreateColVector() const override { return make_unique<VVector<SCAL>>(ndof_); }

    void Mult(const BaseVector& x, BaseVector& y) const override {
        Apply(x);
        auto fy = y.FV<SCAL>();
        auto yl = yl_.FV();
        ParallelFor(ndof_, [&](size_t i) { fy[i] = i < nedge_ ? yl[i] : SCAL(0); });
        for (size_t k = 0; k < extra_.size(); k++) fy[extra_[k]] = yl[nedge_ + k];
    }

    void MultAdd(double s, const BaseVector& x, BaseVector& y) const override {
        AddScaled(SCAL(s), x, y);
    }

    void MultAdd(Complex s, const BaseVector& x, BaseVector& y) const override {
        if constexpr (std::is_same_v<SCAL, Complex>) AddScaled(s, x, y);
        else BaseMatrix::MultAdd(s, x, y);
    }

    void MultTrans(const BaseVector& x, BaseVector& y) const override { Mult(x, y); }
    void MultTransAdd(double s, const BaseVector& x, BaseVector& y) const override { MultAdd(s, x, y); }

private:
    void AddScaled(SCAL s, const BaseVector& x, BaseVector& y) const {
        Apply(x);
        auto fy = y.FV<SCAL>();
        auto yl = yl_.FV();
        ParallelFor(nedge_, [&](size_t i) { fy[i] += s * yl[i]; });
        for (size_t k = 0; k < extra_.size(); k++) fy[extra_[k]] += s * yl[nedge_ + k];
    }

    /// yl = (k stationary AMS steps) applied to the gathered wirebasket part of x.
    void Apply(const BaseVector& x) const {
        auto t0 = std::chrono::steady_clock::now();
        auto fx = x.FV<SCAL>();
        auto xl = xl_.FV();
        ParallelFor(nedge_, [&](size_t i) { xl[i] = fx[i]; });
        for (size_t k = 0; k < extra_.size(); k++) xl[nedge_ + k] = fx[extra_[k]];
        ams_->Mult(xl_, yl_);
        if (cycles_ > 1) {
            auto y = yl_.FV();
            auto r = res_.FV();
            auto c = cor_.FV();
            const auto& S = *system_;
            for (int k = 1; k < cycles_; k++) {
                ParallelForRange(nloc_, [&](IntRange range) {
                    for (auto i : range) {
                        if (!free_->Test(i)) { r[i] = SCAL(0); continue; }
                        auto cols = S.GetRowIndices(i);
                        auto vals = S.GetRowValues(i);
                        SCAL sum = xl[i];
                        for (int j = 0; j < cols.Size(); j++) sum -= vals[j] * y[cols[j]];
                        r[i] = sum;
                    }
                });
                ams_->Mult(res_, cor_);
                ParallelFor(nloc_, [&](size_t i) { y[i] += c[i]; });
            }
        }
        auto& stats = GetAMSWirebasketStats();
        stats.applies++;
        stats.apply_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }

    shared_ptr<BaseMatrix> ams_;
    shared_ptr<SparseMatrix<SCAL>> system_;
    shared_ptr<BitArray> free_;
    std::vector<int> extra_;
    size_t ndof_, nedge_, nloc_;
    int cycles_;
    mutable VVector<SCAL> xl_, yl_, res_, cor_;
};

/// NGSolve preconditioner "sparsesolv_ams": BDDC hands it the wirebasket
/// free dofs (InitLevel) and the assembled wirebasket matrix (FinalizeLevel).
class AMSWirebasketPreconditioner : public ngcomp::Preconditioner {
public:
    AMSWirebasketPreconditioner(shared_ptr<ngcomp::BilinearForm> bfa,
                                const Flags& flags, const string name = "sparsesolv_ams")
        : ngcomp::Preconditioner(bfa, flags, name), bfa_(bfa)
    {
        cycle_type_ = int(flags.GetNumFlag("cycle_type", 1));
        cycles_ = int(flags.GetNumFlag("cycles", 1));
        lean_coarse_ = flags.GetNumFlag("lean_coarse", 1.0) != 0.0;
        num_smooth_ = int(flags.GetNumFlag("num_smooth", 1));
        print_level_ = int(flags.GetNumFlag("print_level", 0));
        eps_ = flags.GetNumFlag("eps", 0.0);
        beta_zero_ = flags.GetDefineFlag("beta_zero");
        mixed_precision_ = flags.GetDefineFlag("mixed_precision");
    }

    void InitLevel(shared_ptr<BitArray> freedofs) override { wb_free_ = freedofs; }

    void FinalizeLevel(const BaseMatrix* mat) override {
        if (!mat) throw Exception("sparsesolv_ams: no wirebasket matrix supplied");
        if (mat->IsComplex()) Build<Complex>(*mat);
        else Build<double>(*mat);
    }

    void Update() override {}

    const BaseMatrix& GetMatrix() const override {
        if (!op_) ThrowPreconditionerNotReady();
        return *op_;
    }
    shared_ptr<BaseMatrix> GetMatrixPtr() override {
        if (!op_) ThrowPreconditionerNotReady();
        return op_;
    }
    const BaseMatrix& GetAMatrix() const override { return bfa_->GetMatrix(); }
    const char* ClassName() const override { return "SparseSolv wirebasket AMS"; }

private:
    template <typename SCAL>
    void Build(const BaseMatrix& mat) {
        using clock = std::chrono::steady_clock;
        auto t0 = clock::now();
        auto fes = bfa_->GetFESpace();
        auto ma = fes->GetMeshAccess();
        if (ma->GetDimension() != 3)
            throw Exception("sparsesolv_ams: needs a three-dimensional mesh");
        const size_t nedge = ma->GetNEdges();
        const size_t nv = ma->GetNV();
        const size_t ndof = fes->GetNDof();

        // The edge block must be the lowest-order edge dofs, dof e on edge e.
        std::atomic<size_t> misnumbered{0};
        ParallelFor(nedge, [&](size_t e) {
            Array<ngcomp::DofId> dnums;
            fes->GetDofNrs(ngfem::NodeId(ngfem::NT_EDGE, e), dnums);
            if (dnums.Size() < 1 || size_t(dnums[0]) != e) misnumbered++;
        });
        if (misnumbered)
            throw Exception("sparsesolv_ams: the space does not number its lowest-order "
                            "edge dofs by edge (" + ToString(misnumbered.load()) + " edges differ); "
                            "it needs an HCurl space");

        // Local numbering: edges 0..nedge-1, then the free wirebasket dofs past them.
        auto free = wb_free_ ? wb_free_ : fes->GetFreeDofs();
        std::vector<int> extra;
        for (size_t i = nedge; i < free->Size(); i++)
            if (free->Test(i)) extra.push_back(int(i));
        const size_t nloc = nedge + extra.size();
        std::vector<int> local(ndof, -1);
        for (size_t e = 0; e < nedge; e++) local[e] = int(e);
        for (size_t k = 0; k < extra.size(); k++) local[extra[k]] = int(nedge + k);
        std::vector<int> global(nloc);
        for (size_t e = 0; e < nedge; e++) global[e] = int(e);
        for (size_t k = 0; k < extra.size(); k++) global[nedge + k] = extra[k];

        auto sym = dynamic_cast<const SparseMatrixSymmetric<SCAL>*>(&mat);
        auto gen = dynamic_cast<const SparseMatrix<SCAL>*>(&mat);
        if (!gen)
            throw Exception("sparsesolv_ams: the wirebasket matrix is not a SparseMatrix");
        if (size_t(gen->Height()) != ndof)
            throw Exception("sparsesolv_ams: the wirebasket matrix does not match the space");

        auto surrogate = [](SCAL v) -> double {
            if constexpr (std::is_same_v<SCAL, Complex>) return v.real() + v.imag();
            else return v;
        };

        // Full local pattern (symmetric storage keeps only the lower triangle).
        Array<int> counts(nloc);
        counts = 0;
        size_t coupled = 0;
        for (size_t li = 0; li < nloc; li++) {
            const int i = global[li];
            auto cols = gen->GetRowIndices(i);
            auto vals = gen->GetRowValues(i);
            for (int k = 0; k < cols.Size(); k++) {
                const int lc = local[cols[k]];
                if (lc < 0) {
                    if (free->Test(i) && free->Test(cols[k]) && vals[k] != SCAL(0)) coupled++;
                    continue;
                }
                counts[li]++;
                if (sym && size_t(lc) != li) counts[lc]++;
            }
        }
        if (coupled)
            throw Exception("sparsesolv_ams: the wirebasket matrix couples " + ToString(coupled)
                            + " free entries outside the wirebasket block");
        auto edge_mat = make_shared<SparseMatrix<double>>(counts, int(nloc));
        auto system = make_shared<SparseMatrix<SCAL>>(counts, int(nloc));  // S itself
        Array<int> fill(nloc);
        fill = 0;
        for (size_t li = 0; li < nloc; li++) {
            const int i = global[li];
            auto cols = gen->GetRowIndices(i);
            auto vals = gen->GetRowValues(i);
            for (int k = 0; k < cols.Size(); k++) {
                const int lc = local[cols[k]];
                if (lc < 0) continue;
                auto put = [&](size_t row, int col) {
                    edge_mat->GetRowIndices(row)[fill[row]] = col;
                    system->GetRowValues(row)[fill[row]] = vals[k];
                    edge_mat->GetRowValues(row)[fill[row]++] = surrogate(vals[k]);
                };
                put(li, lc);
                if (sym && size_t(lc) != li) put(lc, int(li));
            }
        }
        // Rows sorted by column (S shares the pattern); regularize the surrogate diagonal.
        ParallelFor(nloc, [&](size_t i) {
            auto cols = edge_mat->GetRowIndices(i);
            auto vals = edge_mat->GetRowValues(i);
            auto scols = system->GetRowIndices(i);
            auto svals = system->GetRowValues(i);
            const int n = cols.Size();
            std::vector<int> order(n);
            for (int k = 0; k < n; k++) order[k] = k;
            std::sort(order.begin(), order.end(), [&](int a, int b) { return cols[a] < cols[b]; });
            std::vector<int> c(n);
            std::vector<double> v(n);
            std::vector<SCAL> sv(n);
            for (int k = 0; k < n; k++) { c[k] = cols[order[k]]; v[k] = vals[order[k]]; sv[k] = svals[order[k]]; }
            for (int k = 0; k < n; k++) {
                cols[k] = scols[k] = c[k];
                vals[k] = v[k] * (c[k] == int(i) ? 1.0 + eps_ : 1.0);
                svals[k] = sv[k];
            }
        });

        auto local_free = make_shared<BitArray>(nloc);
        local_free->Clear();
        for (size_t li = 0; li < nloc; li++)
            if (free->Test(global[li])) local_free->SetBit(li);
        auto gradient = BuildLowestOrderGradient(*ma, nloc);
        std::vector<double> x(nv), y(nv), z(nv);
        for (size_t v = 0; v < nv; v++) {
            auto p = ma->GetPoint<3>(v);
            x[v] = p(0); y[v] = p(1); z[v] = p(2);
        }
        auto t1 = clock::now();

        const bool has_extra = !extra.empty();
        HypreBasedAMS::AllowActiveTaskManagerSetup allow;
        shared_ptr<BaseMatrix> ams;
        if constexpr (std::is_same_v<SCAL, Complex>) {
            if (beta_zero_ || mixed_precision_)
                throw Exception("sparsesolv_ams: beta_zero and mixed_precision are real-only");
            auto complex_ams = make_shared<ComplexHypreBasedAMS>(edge_mat, gradient, local_free, x, y, z,
                                                    int(nloc), cycle_type_, print_level_,
                                                    1.0, 0, num_smooth_, has_extra, lean_coarse_);
            GetAMSWirebasketStats().complex_ams = complex_ams;
            ams = complex_ams;
        } else {
            ams = make_shared<HypreBasedAMS>(edge_mat, gradient, local_free, x, y, z,
                                             cycle_type_, num_smooth_, 0.25, print_level_,
                                             1.0, 0, beta_zero_, false, mixed_precision_,
                                             has_extra, lean_coarse_);
        }
        auto& stats = GetAMSWirebasketStats();
        stats.builds++;
        stats.n_edges = int(nedge);
        stats.n_extra = int(extra.size());
        stats.n_free = int(local_free->NumSet());
        stats.n_vertices = int(nv);
        stats.cycles = cycles_;
        stats.complex = std::is_same_v<SCAL, Complex>;
        op_ = make_shared<AMSWirebasketOperator<SCAL>>(ams, system, local_free, std::move(extra),
                                                       ndof, nedge, cycles_);
        stats.extract_s = std::chrono::duration<double>(t1 - t0).count();
        stats.setup_s = std::chrono::duration<double>(clock::now() - t1).count();
        stats.applies = 0;
        stats.apply_s = 0;
    }

    shared_ptr<ngcomp::BilinearForm> bfa_;
    shared_ptr<BitArray> wb_free_;
    shared_ptr<BaseMatrix> op_;
    int cycle_type_ = 1, num_smooth_ = 1, print_level_ = 0, cycles_ = 1;
    double eps_ = 0.0;
    bool beta_zero_ = false, mixed_precision_ = false, lean_coarse_ = true;
};

/// Register "sparsesolv_ams" with NGSolve's preconditioner classes (once).
inline void RegisterAMSWirebasket() {
    static bool registered = false;
    if (registered) return;
    registered = true;
    ngcomp::GetPreconditionerClasses().AddPreconditioner(
        "sparsesolv_ams",
        [](shared_ptr<ngcomp::BilinearForm> bfa, const Flags& flags, const string name)
            -> shared_ptr<ngcomp::Preconditioner> {
            return make_shared<AMSWirebasketPreconditioner>(bfa, flags, name);
        });
}

}  // namespace ngla

#endif  // SPARSESOLV_AMS_WIREBASKET_HPP
