/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/// @file compact_ams.hpp
/// @brief Compact HX (Compact Hiptmair-Xu) preconditioner using NGSolve TaskManager
///
/// Compact AMS implementation (~400 lines) built on CompactAMG.
/// Uses Hiptmair-Xu (2007) auxiliary space preconditioning:
///   - Gradient subspace: G^T * A_bc * G solved by CompactAMG
///   - Nodal subspace: Pi^T * A_bc * Pi solved by CompactAMG (component-wise)
///   - Fine-grid smoother: l1-Jacobi (fully TaskManager parallel)
///
/// All parallelism via NGSolve TaskManager (no OpenMP).
///
/// Reference: Hiptmair & Xu, SIAM J. Numer. Anal. 45(6), 2007.
///            Kolev & Vassilevski, J. Comput. Math. 27(5), 2009.

#ifndef SPARSESOLV_COMPACT_AMS_HPP
#define SPARSESOLV_COMPACT_AMS_HPP

#include "compact_amg.hpp"
#include <comp.hpp>
#include <vector>
#include <atomic>
#include <cmath>
#include <chrono>
#include <iostream>

namespace ngla {

/// Compact AMS preconditioner for real HCurl curl-curl + mass systems.
///
/// Implements the auxiliary space method: the HCurl system is preconditioned
/// by combining smoothing on the fine grid with corrections in the gradient
/// (H1) and Nedelec (H1^3) auxiliary spaces.
///
/// Usage:
///   auto ams = make_shared<CompactAMS>(mat, grad, freedofs,
///                                      coord_x, coord_y, coord_z);
///   // Use as preconditioner with COCR or CG
class CompactAMS : public BaseMatrix {
public:
    /// @param mat       HCurl system matrix (SparseMatrix<double>)
    /// @param grad      Discrete gradient G (H1 -> HCurl)
    /// @param freedofs  Free DOFs for HCurl space
    /// @param coord_x/y/z  Vertex coordinates (length = ndof_h1)
    /// @param cycle_type   1=01210 (default), 7=0201020
    /// @param num_smooth   Fine-grid smoother sweeps (default=1)
    /// @param amg_theta    AMG strength threshold (default=0.25)
    /// @param print_level  Verbosity (0=silent)
    /// @param subspace_solver  0=CompactAMG (default), 1=SparseCholesky (diagnostic)
    CompactAMS(shared_ptr<SparseMatrix<double>> mat,
               shared_ptr<SparseMatrix<double>> grad,
               shared_ptr<BitArray> freedofs,
               const std::vector<double>& coord_x,
               const std::vector<double>& coord_y,
               const std::vector<double>& coord_z,
               int cycle_type = 1,
               int num_smooth = 1,
               double amg_theta = 0.25,
               int print_level = 0,
               double correction_weight = 1.0,
               int subspace_solver = 0)
        : mat_(mat), grad_(grad), freedofs_(freedofs),
          ndof_hc_(mat->Height()), ndof_h1_(grad->Width()),
          cycle_type_(cycle_type), num_smooth_(num_smooth),
          print_level_(print_level), correction_weight_(correction_weight),
          subspace_solver_(subspace_solver), amg_theta_(amg_theta)
    {
        if ((int)coord_x.size() != ndof_h1_ ||
            (int)coord_y.size() != ndof_h1_ ||
            (int)coord_z.size() != ndof_h1_)
            throw std::runtime_error("CompactAMS: coordinate size mismatch with H1 DOFs");

        Setup(coord_x, coord_y, coord_z);
    }

    /// Update preconditioner with current matrix values.
    /// Geometry (G, Pi, transposes, work vectors) is preserved.
    /// Rebuilds: A_bc, Galerkin projections, AMG hierarchies, l1 norms.
    void Update() {
        mult_count_ = 0;
        t_smooth_ = t_grad_ = t_nodal_ = t_bc_ = 0;
        RebuildMatrix();
    }

    /// Update with a new system matrix, then rebuild.
    void Update(shared_ptr<SparseMatrix<double>> new_mat) {
        if (new_mat->Height() != mat_->Height() || new_mat->Width() != mat_->Width())
            throw std::invalid_argument("CompactAMS::Update: new matrix dimension ("
                + std::to_string(new_mat->Height()) + "x" + std::to_string(new_mat->Width())
                + ") does not match original ("
                + std::to_string(mat_->Height()) + "x" + std::to_string(mat_->Width()) + ")");
        mat_ = new_mat;
        Update();
    }

    // BaseMatrix interface
    int VHeight() const override { return ndof_hc_; }
    int VWidth() const override { return ndof_hc_; }
    bool IsComplex() const override { return false; }

    AutoVector CreateRowVector() const override { return mat_->CreateRowVector(); }
    AutoVector CreateColVector() const override { return mat_->CreateColVector(); }

    // Const accessors for fused Re/Im preconditioner (ComplexCompactAMS)
    const SparseMatrix<double>& GetAbc() const { return *A_bc_; }
    const SparseMatrix<double>& GetGrad() const { return *grad_; }
    const SparseMatrix<double>& GetGradT() const { return *grad_t_; }
    const SparseMatrix<double>& GetPix() const { return *Pix_; }
    const SparseMatrix<double>& GetPiy() const { return *Piy_; }
    const SparseMatrix<double>& GetPiz() const { return *Piz_; }
    const SparseMatrix<double>& GetPixT() const { return *Pix_t_; }
    const SparseMatrix<double>& GetPiyT() const { return *Piy_t_; }
    const SparseMatrix<double>& GetPizT() const { return *Piz_t_; }
    const BaseMatrix& GetBG() const { return *B_G_; }
    const BaseMatrix& GetBPix() const { return *B_Pix_; }
    const BaseMatrix& GetBPiy() const { return *B_Piy_; }
    const BaseMatrix& GetBPiz() const { return *B_Piz_; }

    // Typed accessors for fused DualMult (returns nullptr if not CompactAMG)
    CompactAMG* GetBGAsAMG() const { return dynamic_cast<CompactAMG*>(B_G_.get()); }
    CompactAMG* GetBPixAsAMG() const { return dynamic_cast<CompactAMG*>(B_Pix_.get()); }
    CompactAMG* GetBPiyAsAMG() const { return dynamic_cast<CompactAMG*>(B_Piy_.get()); }
    CompactAMG* GetBPizAsAMG() const { return dynamic_cast<CompactAMG*>(B_Piz_.get()); }
    const std::vector<double>& GetL1Norms() const { return fine_l1_norms_; }
    int GetNdofHC() const { return ndof_hc_; }
    int GetNdofH1() const { return ndof_h1_; }
    double GetCorrectionWeight() const { return correction_weight_; }
    int GetNumSmooth() const { return num_smooth_; }
    int GetCycleType() const { return cycle_type_; }
    shared_ptr<BitArray> GetFreeDofs() const { return freedofs_; }

    /// Apply one AMS V-cycle: cycle_type=1 -> "01210"
    void Mult(const BaseVector& b, BaseVector& x) const override {
        using clock = std::chrono::high_resolution_clock;
        auto tp = clock::now();
        auto elapsed = [&]() {
            auto now = clock::now();
            double dt = std::chrono::duration<double>(now - tp).count();
            tp = now;
            return dt;
        };

        x = 0;

        // Zero constrained DOFs in RHS
        auto& b0 = *b0_;
        CopyVector(b, b0);
        if (freedofs_) {
            auto fv = b0.FVDouble();
            ParallelFor(ndof_hc_, [&](size_t i) {
                if (!freedofs_->Test(i)) fv[i] = 0;
            });
        }
        t_bc_ += elapsed();

        if (cycle_type_ == 7) {
            FineSmooth(b0, x); t_smooth_ += elapsed();
            GradientCorrect(b0, x); t_grad_ += elapsed();
            FineSmooth(b0, x); t_smooth_ += elapsed();
            NodalCorrect(b0, x); t_nodal_ += elapsed();
            FineSmooth(b0, x); t_smooth_ += elapsed();
            GradientCorrect(b0, x); t_grad_ += elapsed();
            FineSmooth(b0, x); t_smooth_ += elapsed();
        } else {
            // 01210 (default, cycle_type=1)
            FineSmooth(b0, x); t_smooth_ += elapsed();
            GradientCorrect(b0, x); t_grad_ += elapsed();
            NodalCorrect(b0, x); t_nodal_ += elapsed();
            GradientCorrect(b0, x); t_grad_ += elapsed();
            FineSmooth(b0, x); t_smooth_ += elapsed();
        }

        // Zero constrained DOFs in output
        if (freedofs_) {
            auto fv = x.FVDouble();
            ParallelFor(ndof_hc_, [&](size_t i) {
                if (!freedofs_->Test(i)) fv[i] = 0;
            });
        }
        t_bc_ += elapsed();

        mult_count_++;
        if (print_level_ >= 1 && mult_count_ == 25) {
            std::cout << "\n  AMS Mult x" << mult_count_ << " breakdown:"
                      << " smooth=" << t_smooth_ << "s"
                      << " grad=" << t_grad_ << "s"
                      << " nodal=" << t_nodal_ << "s"
                      << " bc=" << t_bc_ << "s"
                      << " total=" << (t_smooth_+t_grad_+t_nodal_+t_bc_) << "s"
                      << std::endl;
        }
    }

    void MultTrans(const BaseVector& b, BaseVector& x) const override {
        Mult(b, x);
    }

private:
    shared_ptr<SparseMatrix<double>> mat_;
    shared_ptr<SparseMatrix<double>> A_bc_;  // BC-modified matrix (identity rows for constrained DOFs)
    shared_ptr<SparseMatrix<double>> grad_;
    shared_ptr<BitArray> freedofs_;
    int ndof_hc_, ndof_h1_;
    int cycle_type_;
    int num_smooth_;
    int print_level_;
    double correction_weight_;
    int subspace_solver_;  // 0=CompactAMG, 1=SparseCholesky
    double amg_theta_;     // AMG strength threshold (preserved for Update)

    // Gradient subspace
    shared_ptr<SparseMatrix<double>> grad_t_;  // G^T
    shared_ptr<SparseMatrix<double>> A_G_;     // G^T * A * G
    shared_ptr<BaseMatrix> B_G_;               // Solver for A_G (AMG or direct)

    // Nodal subspace (component-wise Pix, Piy, Piz)
    shared_ptr<SparseMatrix<double>> Pix_, Piy_, Piz_;
    shared_ptr<SparseMatrix<double>> Pix_t_, Piy_t_, Piz_t_;
    shared_ptr<SparseMatrix<double>> A_Pix_, A_Piy_, A_Piz_;
    shared_ptr<BaseMatrix> B_Pix_, B_Piy_, B_Piz_;

    // Fine-grid l1-Jacobi smoother norms
    std::vector<double> fine_l1_norms_;

    // Work vectors
    mutable std::unique_ptr<VVector<double>> b0_;  // Modified RHS (constrained DOFs zeroed)
    mutable std::unique_ptr<VVector<double>> r0_;  // Fine residual
    mutable std::unique_ptr<VVector<double>> r_G_, g_G_;      // Gradient space
    mutable std::unique_ptr<VVector<double>> r_Pix_, g_Pix_;  // Pix space
    mutable std::unique_ptr<VVector<double>> r_Piy_, g_Piy_;  // Piy space
    mutable std::unique_ptr<VVector<double>> r_Piz_, g_Piz_;  // Piz space

    // Accumulated timing (mutable for const Mult)
    mutable int mult_count_ = 0;
    mutable double t_smooth_ = 0, t_grad_ = 0, t_nodal_ = 0, t_bc_ = 0;

    // =====================================================================
    // Setup: split into geometry (one-time) + matrix (per-Update) + alloc
    // =====================================================================
    void Setup(const std::vector<double>& cx,
               const std::vector<double>& cy,
               const std::vector<double>& cz) {
        auto t_total = std::chrono::high_resolution_clock::now();
        if (print_level_ > 0)
            std::cout << "CompactAMS setup: HC=" << ndof_hc_
                      << " H1=" << ndof_h1_ << std::flush;

        // Phase 1: Geometry-dependent setup (one-time)
        SetupGeometry(cx, cy, cz);

        // Phase 2: Matrix-dependent setup (repeated on Update)
        RebuildMatrix();

        // Phase 3: Work vector allocation (one-time)
        AllocateWorkVectors();

        // Correction weight
        if (correction_weight_ <= 0.0)
            correction_weight_ = 1.0;
        if (print_level_ > 0)
            std::cout << "\n  correction_weight = " << correction_weight_;

        auto t_end = std::chrono::high_resolution_clock::now();
        if (print_level_ > 0)
            std::cout << "\n  CompactAMS setup complete: "
                      << std::chrono::duration<double>(t_end - t_total).count()
                      << "s total" << std::endl;
    }

    // =====================================================================
    // SetupGeometry: Pi matrices + transposes (geometry-only, one-time)
    // =====================================================================
    void SetupGeometry(const std::vector<double>& cx,
                       const std::vector<double>& cy,
                       const std::vector<double>& cz) {
        BuildPiComponents(cx, cy, cz);
        ParallelFor(4, [&](size_t d) {
            if (d == 0)
                grad_t_ = dynamic_pointer_cast<SparseMatrix<double>>(
                    grad_->CreateTranspose(true));
            else if (d == 1)
                Pix_t_ = dynamic_pointer_cast<SparseMatrix<double>>(
                    Pix_->CreateTranspose(true));
            else if (d == 2)
                Piy_t_ = dynamic_pointer_cast<SparseMatrix<double>>(
                    Piy_->CreateTranspose(true));
            else
                Piz_t_ = dynamic_pointer_cast<SparseMatrix<double>>(
                    Piz_->CreateTranspose(true));
        });
    }

    // =====================================================================
    // RebuildMatrix: A_bc, Galerkin projections, AMG, l1 norms
    // Called on every Update() — geometry is preserved.
    // =====================================================================
    void RebuildMatrix() {
        auto t_rebuild = std::chrono::high_resolution_clock::now();

        // 1. Create BC-modified matrix (identity rows for constrained DOFs).
        if (freedofs_) {
            A_bc_ = CreateBCModifiedMatrix(*mat_, *freedofs_);
        } else {
            A_bc_ = mat_;
        }

        // 2. Galerkin projection: A_G = G^T * A_bc * G, etc.
        auto t0 = std::chrono::high_resolution_clock::now();
        ParallelFor(4, [&](size_t d) {
            if (d == 0)
                A_G_ = dynamic_pointer_cast<SparseMatrix<double>>(A_bc_->Restrict(*grad_));
            else if (d == 1)
                A_Pix_ = dynamic_pointer_cast<SparseMatrix<double>>(A_bc_->Restrict(*Pix_));
            else if (d == 2)
                A_Piy_ = dynamic_pointer_cast<SparseMatrix<double>>(A_bc_->Restrict(*Piy_));
            else
                A_Piz_ = dynamic_pointer_cast<SparseMatrix<double>>(A_bc_->Restrict(*Piz_));
        });
        auto t1 = std::chrono::high_resolution_clock::now();
        double dt_restrict = std::chrono::duration<double>(t1 - t0).count();

        if (print_level_ > 0) {
            std::cout << "\n  Galerkin restrict: " << dt_restrict << "s"
                      << " A_G=" << A_G_->Height() << "x" << A_G_->Width()
                      << " nnz=" << A_G_->NZE() << std::flush;
        }

        // 3. Fix zero rows: set diag = 1.0 for truly zero rows.
        FixZeroRows(*A_G_);
        FixZeroRows(*A_Pix_);
        FixZeroRows(*A_Piy_);
        FixZeroRows(*A_Piz_);

        // 4. Build solvers for auxiliary spaces
        if (subspace_solver_ == 1) {
            if (print_level_ > 0)
                std::cout << "\n  Subspace solver: SparseCholesky (diagnostic)" << std::flush;

            t0 = std::chrono::high_resolution_clock::now();
            A_G_->SetInverseType(SPARSECHOLESKY);
            B_G_ = A_G_->InverseMatrix(shared_ptr<BitArray>(nullptr));
            t1 = std::chrono::high_resolution_clock::now();
            if (print_level_ > 0)
                std::cout << "\n  B_G setup: " << std::chrono::duration<double>(t1 - t0).count()
                          << "s (direct, n=" << A_G_->Height() << ")" << std::flush;

            t0 = std::chrono::high_resolution_clock::now();
            A_Pix_->SetInverseType(SPARSECHOLESKY);
            A_Piy_->SetInverseType(SPARSECHOLESKY);
            A_Piz_->SetInverseType(SPARSECHOLESKY);
            B_Pix_ = A_Pix_->InverseMatrix(shared_ptr<BitArray>(nullptr));
            B_Piy_ = A_Piy_->InverseMatrix(shared_ptr<BitArray>(nullptr));
            B_Piz_ = A_Piz_->InverseMatrix(shared_ptr<BitArray>(nullptr));
            t1 = std::chrono::high_resolution_clock::now();
            if (print_level_ > 0)
                std::cout << "\n  B_Pi setup: " << std::chrono::duration<double>(t1 - t0).count()
                          << "s (direct)" << std::flush;
        } else {
            int min_coarse_aux = 500;
            if (print_level_ > 0)
                std::cout << "\n  Subspace solver: CompactAMG (min_coarse="
                          << min_coarse_aux << ")" << std::flush;

            auto amg_G = make_shared<CompactAMG>(A_G_, nullptr, amg_theta_, 25, min_coarse_aux, 1, 0);
            auto amg_Pix = make_shared<CompactAMG>(A_Pix_, nullptr, amg_theta_, 25, min_coarse_aux, 1, 0);
            auto amg_Piy = make_shared<CompactAMG>(A_Piy_, nullptr, amg_theta_, 25, min_coarse_aux, 1, 0);
            auto amg_Piz = make_shared<CompactAMG>(A_Piz_, nullptr, amg_theta_, 25, min_coarse_aux, 1, 0);

            t0 = std::chrono::high_resolution_clock::now();
            ParallelFor(4, [&](size_t d) {
                if (d == 0) amg_G->Setup();
                else if (d == 1) amg_Pix->Setup();
                else if (d == 2) amg_Piy->Setup();
                else amg_Piz->Setup();
            });
            t1 = std::chrono::high_resolution_clock::now();
            if (print_level_ > 0)
                std::cout << "\n  AMG setup (4-way parallel): "
                          << std::chrono::duration<double>(t1 - t0).count()
                          << "s, levels: G=" << amg_G->NumLevels()
                          << " Px=" << amg_Pix->NumLevels()
                          << " Py=" << amg_Piy->NumLevels()
                          << " Pz=" << amg_Piz->NumLevels() << std::flush;

            B_G_ = amg_G;
            B_Pix_ = amg_Pix;
            B_Piy_ = amg_Piy;
            B_Piz_ = amg_Piz;
        }

        // 5. Compute truncated l1 norms for fine-grid smoother
        {
            fine_l1_norms_.resize(ndof_hc_);
            ParallelFor(ndof_hc_, [&](size_t i) {
                auto cols = A_bc_->GetRowIndices(i);
                auto vals = A_bc_->GetRowValues(i);
                double diag = 0, off_diag = 0;
                for (int j = 0; j < cols.Size(); j++) {
                    if (cols[j] == (int)i)
                        diag = std::abs(vals[j]);
                    else
                        off_diag += std::abs(vals[j]);
                }
                double temp = diag + 0.5 * off_diag;
                if (temp <= (4.0/3.0) * diag)
                    fine_l1_norms_[i] = (diag > 0) ? diag : 1.0;
                else
                    fine_l1_norms_[i] = (temp > 0) ? temp : 1.0;
            });
        }

        if (print_level_ > 0) {
            auto t_end = std::chrono::high_resolution_clock::now();
            std::cout << "\n  RebuildMatrix: "
                      << std::chrono::duration<double>(t_end - t_rebuild).count()
                      << "s" << std::flush;
        }
    }

    // =====================================================================
    // AllocateWorkVectors (one-time, sizes depend on DOF counts only)
    // =====================================================================
    void AllocateWorkVectors() {
        b0_ = std::make_unique<VVector<double>>(ndof_hc_);
        r0_ = std::make_unique<VVector<double>>(ndof_hc_);
        r_G_ = std::make_unique<VVector<double>>(ndof_h1_);
        g_G_ = std::make_unique<VVector<double>>(ndof_h1_);
        r_Pix_ = std::make_unique<VVector<double>>(ndof_h1_);
        g_Pix_ = std::make_unique<VVector<double>>(ndof_h1_);
        r_Piy_ = std::make_unique<VVector<double>>(ndof_h1_);
        g_Piy_ = std::make_unique<VVector<double>>(ndof_h1_);
        r_Piz_ = std::make_unique<VVector<double>>(ndof_h1_);
        g_Piz_ = std::make_unique<VVector<double>>(ndof_h1_);
    }

    // =====================================================================
    // Build Pi component matrices from G and coordinates
    // =====================================================================
    // Kolev-Vassilevski formula:
    //   Gd[e] = (G * coord_d)[e] = edge vector in d-direction
    //   Pix[e, v] = |G[e,v]| * 0.5 * Gx[e]
    //   Piy[e, v] = |G[e,v]| * 0.5 * Gy[e]
    //   Piz[e, v] = |G[e,v]| * 0.5 * Gz[e]
    //
    // For edge (v1,v2) with G[e,v1]=-1, G[e,v2]=+1:
    //   Gx[e] = x[v2] - x[v1] (edge vector x-component)
    //   Pix[e,v1] = Pix[e,v2] = 0.5 * (x[v2] - x[v1])
    //
    // Both vertices get the SAME value (half the edge vector).
    // This is the correct Nedelec interpolation operator.
    void BuildPiComponents(const std::vector<double>& cx,
                           const std::vector<double>& cy,
                           const std::vector<double>& cz) {
        int nedges = grad_->Height();

        // Step 1: Compute Gd[e] = (G * coord_d)[e] for d=x,y,z
        // This is the edge vector in each coordinate direction.
        std::vector<double> Gx(nedges, 0.0), Gy(nedges, 0.0), Gz(nedges, 0.0);
        ParallelFor(nedges, [&](size_t e) {
            auto cols = grad_->GetRowIndices(e);
            auto vals = grad_->GetRowValues(e);
            double gx = 0, gy = 0, gz = 0;
            for (int j = 0; j < cols.Size(); j++) {
                int v = cols[j];
                gx += vals[j] * cx[v];
                gy += vals[j] * cy[v];
                gz += vals[j] * cz[v];
            }
            Gx[e] = gx;
            Gy[e] = gy;
            Gz[e] = gz;
        });

        // Step 2: Build Pi matrices with same sparsity as G
        Array<int> cnt(nedges);
        for (int i = 0; i < nedges; i++)
            cnt[i] = grad_->GetRowIndices(i).Size();

        Pix_ = make_shared<SparseMatrix<double>>(cnt, ndof_h1_);
        Piy_ = make_shared<SparseMatrix<double>>(cnt, ndof_h1_);
        Piz_ = make_shared<SparseMatrix<double>>(cnt, ndof_h1_);

        double* Gd_ptrs[3] = {Gx.data(), Gy.data(), Gz.data()};
        SparseMatrix<double>* Pi_mats[3] = {Pix_.get(), Piy_.get(), Piz_.get()};

        ParallelFor(nedges, [&](size_t e) {
            auto g_cols = grad_->GetRowIndices(e);
            auto g_vals = grad_->GetRowValues(e);

            for (int d = 0; d < 3; d++) {
                auto pi_cols = Pi_mats[d]->GetRowIndices(e);
                auto pi_vals = Pi_mats[d]->GetRowValues(e);

                for (int j = 0; j < g_cols.Size(); j++) {
                    pi_cols[j] = g_cols[j];
                    // Pi formula: |G[e,v]| * 0.5 * Gd[e]
                    pi_vals[j] = std::abs(g_vals[j]) * 0.5 * Gd_ptrs[d][e];
                }
            }
        });
    }

    // =====================================================================
    // Fine-grid smoother: l1-Jacobi (fully TaskManager parallel)
    //
    // x += D_l1^{-1} * (b - A*x) per sweep, where D_l1 = truncated l1 norms.
    // Unlike GS, l1-Jacobi has NO data dependencies between rows, enabling
    // full parallelism via TaskManager. Each sweep requires an explicit
    // residual computation (A*x), which is also parallel (NGSolve SpMV).
    //
    // Reference: Baker et al., "Multigrid Smoothers for Ultraparallel Computing",
    //            SIAM J. Sci. Comput. 33(5), 2011, Remark 6.2.
    // =====================================================================

    /// Copy BaseVector data via FlatVector
    static void CopyVector(const BaseVector& src, BaseVector& dst) {
        auto fv_s = src.FVDouble();
        auto fv_d = dst.FVDouble();
        ParallelFor(fv_s.Size(), [&](size_t i) { fv_d[i] = fv_s[i]; });
    }

    void FineSmooth(const BaseVector& b, BaseVector& x) const {
        // l1-Jacobi: fully parallel, no data dependency between rows.
        // Each sweep: compute residual r = b - A*x, then x += r / l1_norm.
        auto& res = *r0_;

        for (int s = 0; s < num_smooth_; s++) {
            // Residual: r = b - A*x (NGSolve SpMV is TaskManager-parallel)
            CopyVector(b, res);
            A_bc_->MultAdd(-1.0, x, res);

            // Jacobi update: x[i] += r[i] / l1_norm[i] (fully parallel)
            auto fv_x = x.FVDouble();
            auto fv_r = res.FVDouble();
            ParallelFor(ndof_hc_, [&](size_t i) {
                fv_x[i] += fv_r[i] / fine_l1_norms_[i];
            });
        }
    }

    // =====================================================================
    // Subspace corrections
    // =====================================================================
    void SubspaceCorrect(const SparseMatrix<double>& P,
                         const SparseMatrix<double>& Pt,
                         const BaseMatrix& B,
                         const BaseVector& residual,
                         BaseVector& x,
                         BaseVector& r_c,
                         BaseVector& g_c,
                         const char* label = "") const {
        // Restrict: r_c = P^T * residual
        Pt.Mult(residual, r_c);

        // Coarse solve: g_c = B^{-1} * r_c (one AMG V-cycle)
        g_c.FVDouble() = 0;
        B.Mult(r_c, g_c);

        // Prolongate and add: x += omega * P * g_c
        P.MultAdd(correction_weight_, g_c, x);
    }

    void ComputeResidual(const BaseVector& b, const BaseVector& x, BaseVector& res) const {
        CopyVector(b, res);
        A_bc_->MultAdd(-1.0, x, res);
    }

    void GradientCorrect(const BaseVector& b, BaseVector& x) const {
        ComputeResidual(b, x, *r0_);
        SubspaceCorrect(*grad_, *grad_t_, *B_G_, *r0_, x, *r_G_, *g_G_, "G");
    }

    void NodalCorrect(const BaseVector& b, BaseVector& x) const {
        ComputeResidual(b, x, *r0_);

        // Additive Pi corrections: all 3 use the same residual (no recomputation).
        // Saves 2 fine-level SpMVs per AMS cycle vs multiplicative approach.

        // Restrict to all 3 subspaces from same residual
        Pix_t_->Mult(*r0_, *r_Pix_);
        Piy_t_->Mult(*r0_, *r_Piy_);
        Piz_t_->Mult(*r0_, *r_Piz_);

        // Solve all 3 (sequential: each AMG uses TaskManager internally)
        g_Pix_->FVDouble() = 0;
        B_Pix_->Mult(*r_Pix_, *g_Pix_);
        g_Piy_->FVDouble() = 0;
        B_Piy_->Mult(*r_Piy_, *g_Piy_);
        g_Piz_->FVDouble() = 0;
        B_Piz_->Mult(*r_Piz_, *g_Piz_);

        // Prolongate and add all 3 corrections
        Pix_->MultAdd(correction_weight_, *g_Pix_, x);
        Piy_->MultAdd(correction_weight_, *g_Piy_, x);
        Piz_->MultAdd(correction_weight_, *g_Piz_, x);
    }

    // =====================================================================
    // Utility
    // =====================================================================
    /// Create a copy of A with identity rows for constrained DOFs.
    shared_ptr<SparseMatrix<double>> CreateBCModifiedMatrix(
        const SparseMatrix<double>& A, const BitArray& freedofs) const
    {
        int n = A.Height();
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
                    dst_vals[j] = (src_cols[j] == (int)i) ? 1.0 : 0.0;
                } else if (!freedofs.Test(src_cols[j])) {
                    dst_vals[j] = 0.0;
                } else {
                    dst_vals[j] = src_vals[j];
                }
            }
        });
        return B;
    }

    /// Fix zero rows: set diagonal to 1.0 for truly zero rows.
    /// Only modifies rows where the l1 norm is exactly zero (truly zero rows).
    /// Sets diagonal to 1.0 for such rows. Does NOT add epsilon to all diagonals.
    void FixZeroRows(SparseMatrix<double>& A) const {
        int n = A.Height();
        std::atomic<int> n_fixed{0};
        ParallelFor(n, [&](size_t i) {
            auto cols = A.GetRowIndices(i);
            auto vals = A.GetRowValues(i);
            double l1_norm = 0;
            int diag_pos = -1;
            for (int j = 0; j < cols.Size(); j++) {
                l1_norm += std::abs(vals[j]);
                if (cols[j] == (int)i) diag_pos = j;
            }
            if (l1_norm == 0.0 && diag_pos >= 0) {
                vals[diag_pos] = 1.0;
                n_fixed.fetch_add(1, std::memory_order_relaxed);
            }
        });
        if (n_fixed.load() > 0 && print_level_ > 0) {
            std::cout << "\n  FixZeroRows: fixed " << n_fixed.load()
                      << " / " << n << " zero rows" << std::flush;
        }
    }
};

}  // namespace ngla

#endif  // SPARSESOLV_COMPACT_AMS_HPP
