/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/// @file complex_compact_ams.hpp
/// @brief Complex AMS preconditioner with fused Re/Im SpMV operations
///
/// For complex eddy current problems (A = K + jw*sigma*M), applies CompactAMS
/// to real and imaginary parts with fused matrix operations.
///
/// Key optimization: fine-level SpMV is memory-bandwidth-bound. Loading matrix
/// row data once and processing both Re and Im halves the bandwidth cost.
/// AMG V-cycles (coarse levels) remain sequential per Re/Im.

#ifndef SPARSESOLV_COMPLEX_COMPACT_AMS_HPP
#define SPARSESOLV_COMPLEX_COMPACT_AMS_HPP

#include "compact_ams.hpp"
#include <comp.hpp>

namespace ngla {

/// Complex Compact AMS preconditioner with fused Re/Im SpMV.
///
/// Fine-level operations (l1-Jacobi, residual, restrict, prolongate) process
/// Re and Im simultaneously in a single pass over matrix rows. This halves
/// memory bandwidth for these operations. Coarse-level AMG V-cycles run
/// sequentially (Re then Im) using the shared CompactAMS hierarchy.
class ComplexCompactAMS : public BaseMatrix {
public:
    ComplexCompactAMS(
        shared_ptr<SparseMatrix<double>> a_real_mat,
        shared_ptr<SparseMatrix<double>> grad_mat,
        shared_ptr<BitArray> freedofs,
        const std::vector<double>& coord_x,
        const std::vector<double>& coord_y,
        const std::vector<double>& coord_z,
        int ndof_complex,
        int cycle_type = 1,
        int print_level = 0,
        double correction_weight = 1.0,
        int subspace_solver = 0,
        int num_smooth = 1)
        : ndof_complex_(ndof_complex > 0 ? ndof_complex : a_real_mat->Height())
    {
        // Build the real AMS hierarchy (shared for Re and Im)
        ams_ = make_shared<CompactAMS>(
            a_real_mat, grad_mat, freedofs,
            coord_x, coord_y, coord_z,
            cycle_type, num_smooth, 0.25, print_level, correction_weight,
            subspace_solver);

        // Cache accessors
        ndof_hc_ = ams_->GetNdofHC();
        ndof_h1_ = ams_->GetNdofH1();
        correction_weight_ = ams_->GetCorrectionWeight();
        num_smooth_ = ams_->GetNumSmooth();
        cycle_type_ = ams_->GetCycleType();

        // Allocate all work vectors (Re + Im pairs)
        b_re_ = std::make_unique<VVector<double>>(ndof_hc_);
        b_im_ = std::make_unique<VVector<double>>(ndof_hc_);
        x_re_ = std::make_unique<VVector<double>>(ndof_hc_);
        x_im_ = std::make_unique<VVector<double>>(ndof_hc_);
        res_re_ = std::make_unique<VVector<double>>(ndof_hc_);
        res_im_ = std::make_unique<VVector<double>>(ndof_hc_);

        rG_re_ = std::make_unique<VVector<double>>(ndof_h1_);
        rG_im_ = std::make_unique<VVector<double>>(ndof_h1_);
        gG_re_ = std::make_unique<VVector<double>>(ndof_h1_);
        gG_im_ = std::make_unique<VVector<double>>(ndof_h1_);

        rPx_re_ = std::make_unique<VVector<double>>(ndof_h1_);
        rPx_im_ = std::make_unique<VVector<double>>(ndof_h1_);
        gPx_re_ = std::make_unique<VVector<double>>(ndof_h1_);
        gPx_im_ = std::make_unique<VVector<double>>(ndof_h1_);
        rPy_re_ = std::make_unique<VVector<double>>(ndof_h1_);
        rPy_im_ = std::make_unique<VVector<double>>(ndof_h1_);
        gPy_re_ = std::make_unique<VVector<double>>(ndof_h1_);
        gPy_im_ = std::make_unique<VVector<double>>(ndof_h1_);
        rPz_re_ = std::make_unique<VVector<double>>(ndof_h1_);
        rPz_im_ = std::make_unique<VVector<double>>(ndof_h1_);
        gPz_re_ = std::make_unique<VVector<double>>(ndof_h1_);
        gPz_im_ = std::make_unique<VVector<double>>(ndof_h1_);

        if (print_level > 0)
            std::cout << "  ComplexCompactAMS: fused Re/Im SpMV "
                      << "(ndof=" << ndof_complex_ << ")" << std::endl;
    }

    /// Update preconditioner with current matrix values (re-reads stored matrix).
    void Update() {
        ams_->Update();
    }

    /// Update with a new real auxiliary matrix, then rebuild.
    void Update(shared_ptr<SparseMatrix<double>> new_a_real) {
        ams_->Update(new_a_real);
    }

    int VHeight() const override { return ndof_complex_; }
    int VWidth() const override { return ndof_complex_; }
    bool IsComplex() const override { return true; }

    AutoVector CreateRowVector() const override {
        return std::make_unique<VVector<Complex>>(ndof_complex_);
    }
    AutoVector CreateColVector() const override {
        return std::make_unique<VVector<Complex>>(ndof_complex_);
    }

    void Mult(const BaseVector& f, BaseVector& u) const override {
        auto fv_f = f.FVComplex();

        auto b_re = b_re_->FVDouble();
        auto b_im = b_im_->FVDouble();
        auto x_re = x_re_->FVDouble();
        auto x_im = x_im_->FVDouble();

        // Split complex -> Re/Im and initialize solution to zero
        auto fd = ams_->GetFreeDofs();
        ParallelFor(ndof_hc_, [&](size_t i) {
            bool free = !fd || fd->Test(i);
            b_re[i] = free ? fv_f[i].real() : 0.0;
            b_im[i] = free ? fv_f[i].imag() : 0.0;
            x_re[i] = 0.0;
            x_im[i] = 0.0;
        });

        // AMS V-cycle with fused Re/Im operations
        if (cycle_type_ == 7) {
            FusedFineSmooth(b_re, x_re, b_im, x_im);
            FusedGradientCorrect(b_re, x_re, b_im, x_im);
            FusedFineSmooth(b_re, x_re, b_im, x_im);
            FusedNodalCorrect(b_re, x_re, b_im, x_im);
            FusedFineSmooth(b_re, x_re, b_im, x_im);
            FusedGradientCorrect(b_re, x_re, b_im, x_im);
            FusedFineSmooth(b_re, x_re, b_im, x_im);
        } else {
            // cycle_type=1: "01210"
            FusedFineSmooth(b_re, x_re, b_im, x_im);
            FusedGradientCorrect(b_re, x_re, b_im, x_im);
            FusedNodalCorrect(b_re, x_re, b_im, x_im);
            FusedGradientCorrect(b_re, x_re, b_im, x_im);
            FusedFineSmooth(b_re, x_re, b_im, x_im);
        }

        // Combine Re/Im -> complex, zero constrained DOFs
        auto fv_u = u.FVComplex();
        ParallelFor(ndof_hc_, [&](size_t i) {
            bool free = !fd || fd->Test(i);
            fv_u[i] = free ? Complex(x_re[i], x_im[i]) : Complex(0.0, 0.0);
        });
    }

    void MultTrans(const BaseVector& f, BaseVector& u) const override {
        Mult(f, u);
    }

private:
    int ndof_complex_;
    int ndof_hc_, ndof_h1_;
    double correction_weight_;
    int num_smooth_;
    int cycle_type_;

    shared_ptr<CompactAMS> ams_;

    // Work vectors (Re + Im pairs)
    mutable std::unique_ptr<VVector<double>> b_re_, b_im_;
    mutable std::unique_ptr<VVector<double>> x_re_, x_im_;
    mutable std::unique_ptr<VVector<double>> res_re_, res_im_;
    mutable std::unique_ptr<VVector<double>> rG_re_, rG_im_, gG_re_, gG_im_;
    mutable std::unique_ptr<VVector<double>> rPx_re_, rPx_im_, gPx_re_, gPx_im_;
    mutable std::unique_ptr<VVector<double>> rPy_re_, rPy_im_, gPy_re_, gPy_im_;
    mutable std::unique_ptr<VVector<double>> rPz_re_, rPz_im_, gPz_re_, gPz_im_;

    // =====================================================================
    // Fused operations: process Re and Im in single pass over matrix rows
    // =====================================================================

    /// Fused l1-Jacobi smooth: x += (b - A*x) / l1_norm for Re and Im simultaneously.
    /// Two-phase: (1) fused residual computation, (2) Jacobi update.
    /// Phase 1 loads each matrix row ONCE for both Re and Im (halves bandwidth).
    /// Must be two phases to avoid reading partially-updated x (Jacobi, not GS).
    void FusedFineSmooth(FlatVector<double> b_re, FlatVector<double> x_re,
                         FlatVector<double> b_im, FlatVector<double> x_im) const {
        const auto& l1 = ams_->GetL1Norms();
        auto res_re = res_re_->FVDouble();
        auto res_im = res_im_->FVDouble();

        for (int s = 0; s < num_smooth_; s++) {
            // Phase 1: Fused residual (reads OLD x entirely, writes to res)
            FusedResidual(b_re, x_re, res_re, b_im, x_im, res_im);

            // Phase 2: Jacobi update from residual (no data dependency)
            ParallelFor(ndof_hc_, [&](size_t i) {
                double inv_l1 = 1.0 / l1[i];
                x_re[i] += res_re[i] * inv_l1;
                x_im[i] += res_im[i] * inv_l1;
            });
        }
    }

    /// Fused residual: res = b - A*x for Re and Im in single matrix pass.
    void FusedResidual(FlatVector<double> b_re, FlatVector<double> x_re,
                       FlatVector<double> res_re,
                       FlatVector<double> b_im, FlatVector<double> x_im,
                       FlatVector<double> res_im) const {
        const auto& A = ams_->GetAbc();
        ParallelFor(ndof_hc_, [&](size_t i) {
            auto cols = A.GetRowIndices(i);
            auto vals = A.GetRowValues(i);
            double dot_re = 0, dot_im = 0;
            for (int j = 0; j < cols.Size(); j++) {
                int c = cols[j];
                double v = vals[j];
                dot_re += v * x_re[c];
                dot_im += v * x_im[c];
            }
            res_re[i] = b_re[i] - dot_re;
            res_im[i] = b_im[i] - dot_im;
        });
    }

    /// Fused restrict: r_c = P^T * res for Re and Im in single matrix pass.
    void FusedRestrict(const SparseMatrix<double>& Pt,
                       FlatVector<double> res_re, FlatVector<double> r_c_re,
                       FlatVector<double> res_im, FlatVector<double> r_c_im) const {
        int nc = Pt.Height();
        ParallelFor(nc, [&](size_t i) {
            auto cols = Pt.GetRowIndices(i);
            auto vals = Pt.GetRowValues(i);
            double s_re = 0, s_im = 0;
            for (int j = 0; j < cols.Size(); j++) {
                int c = cols[j];
                double v = vals[j];
                s_re += v * res_re[c];
                s_im += v * res_im[c];
            }
            r_c_re[i] = s_re;
            r_c_im[i] = s_im;
        });
    }

    /// Fused prolongate: x += w * P * g for Re and Im in single matrix pass.
    void FusedProlongate(const SparseMatrix<double>& P, double w,
                         FlatVector<double> g_re, FlatVector<double> x_re,
                         FlatVector<double> g_im, FlatVector<double> x_im) const {
        int nrows = P.Height();
        ParallelFor(nrows, [&](size_t i) {
            auto cols = P.GetRowIndices(i);
            auto vals = P.GetRowValues(i);
            double s_re = 0, s_im = 0;
            for (int j = 0; j < cols.Size(); j++) {
                int c = cols[j];
                double v = vals[j];
                s_re += v * g_re[c];
                s_im += v * g_im[c];
            }
            x_re[i] += w * s_re;
            x_im[i] += w * s_im;
        });
    }

    // =====================================================================
    // Subspace corrections (fused fine-level, sequential AMG V-cycles)
    // =====================================================================

    void FusedGradientCorrect(FlatVector<double> b_re, FlatVector<double> x_re,
                              FlatVector<double> b_im, FlatVector<double> x_im) const {
        auto res_re = res_re_->FVDouble();
        auto res_im = res_im_->FVDouble();

        // Fused residual: 1 matrix pass for both Re and Im
        FusedResidual(b_re, x_re, res_re, b_im, x_im, res_im);

        // Fused restrict to gradient space
        auto rG_re = rG_re_->FVDouble();
        auto rG_im = rG_im_->FVDouble();
        FusedRestrict(ams_->GetGradT(), res_re, rG_re, res_im, rG_im);

        // AMG V-cycle: fused DualMult or sequential fallback
        auto* amg = ams_->GetBGAsAMG();
        if (amg) {
            amg->DualMult(*rG_re_, *gG_re_, *rG_im_, *gG_im_);
        } else {
            gG_re_->FVDouble() = 0;
            ams_->GetBG().Mult(*rG_re_, *gG_re_);
            gG_im_->FVDouble() = 0;
            ams_->GetBG().Mult(*rG_im_, *gG_im_);
        }

        // Fused prolongate
        FusedProlongate(ams_->GetGrad(), correction_weight_,
                        gG_re_->FVDouble(), x_re, gG_im_->FVDouble(), x_im);
    }

    void FusedNodalCorrect(FlatVector<double> b_re, FlatVector<double> x_re,
                           FlatVector<double> b_im, FlatVector<double> x_im) const {
        auto res_re = res_re_->FVDouble();
        auto res_im = res_im_->FVDouble();

        // Fused residual: 1 matrix pass for both Re and Im
        FusedResidual(b_re, x_re, res_re, b_im, x_im, res_im);

        // Fused restrict to all 3 Pi subspaces (from same residual)
        auto rPx_re = rPx_re_->FVDouble(); auto rPx_im = rPx_im_->FVDouble();
        auto rPy_re = rPy_re_->FVDouble(); auto rPy_im = rPy_im_->FVDouble();
        auto rPz_re = rPz_re_->FVDouble(); auto rPz_im = rPz_im_->FVDouble();

        FusedRestrict(ams_->GetPixT(), res_re, rPx_re, res_im, rPx_im);
        FusedRestrict(ams_->GetPiyT(), res_re, rPy_re, res_im, rPy_im);
        FusedRestrict(ams_->GetPizT(), res_re, rPz_re, res_im, rPz_im);

        // AMG V-cycles: fused DualMult (3 calls instead of 6)
        auto* amg_x = ams_->GetBPixAsAMG();
        auto* amg_y = ams_->GetBPiyAsAMG();
        auto* amg_z = ams_->GetBPizAsAMG();

        if (amg_x && amg_y && amg_z) {
            amg_x->DualMult(*rPx_re_, *gPx_re_, *rPx_im_, *gPx_im_);
            amg_y->DualMult(*rPy_re_, *gPy_re_, *rPy_im_, *gPy_im_);
            amg_z->DualMult(*rPz_re_, *gPz_re_, *rPz_im_, *gPz_im_);
        } else {
            gPx_re_->FVDouble() = 0; ams_->GetBPix().Mult(*rPx_re_, *gPx_re_);
            gPx_im_->FVDouble() = 0; ams_->GetBPix().Mult(*rPx_im_, *gPx_im_);
            gPy_re_->FVDouble() = 0; ams_->GetBPiy().Mult(*rPy_re_, *gPy_re_);
            gPy_im_->FVDouble() = 0; ams_->GetBPiy().Mult(*rPy_im_, *gPy_im_);
            gPz_re_->FVDouble() = 0; ams_->GetBPiz().Mult(*rPz_re_, *gPz_re_);
            gPz_im_->FVDouble() = 0; ams_->GetBPiz().Mult(*rPz_im_, *gPz_im_);
        }

        // Fused prolongate all 3 corrections
        FusedProlongate(ams_->GetPix(), correction_weight_,
                        gPx_re_->FVDouble(), x_re, gPx_im_->FVDouble(), x_im);
        FusedProlongate(ams_->GetPiy(), correction_weight_,
                        gPy_re_->FVDouble(), x_re, gPy_im_->FVDouble(), x_im);
        FusedProlongate(ams_->GetPiz(), correction_weight_,
                        gPz_re_->FVDouble(), x_re, gPz_im_->FVDouble(), x_im);
    }
};

}  // namespace ngla

#endif  // SPARSESOLV_COMPLEX_COMPACT_AMS_HPP
