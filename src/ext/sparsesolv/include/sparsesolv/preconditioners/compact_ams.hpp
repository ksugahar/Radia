/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/// @file compact_ams.hpp
/// @brief HypreBasedAMS -- a HYPRE-free reimplementation of the HYPRE AMS algorithm, on NGSolve TaskManager
///
/// This is NOT a new "compact" method: it is the SAME Auxiliary-space Maxwell Solver (AMS) algorithm that
/// HYPRE's AMS implements -- Hiptmair-Xu (2007) auxiliary-space preconditioning in the Kolev-Vassilevski
/// (2009) parallel form (Kolev & Vassilevski are the HYPRE AMS authors) -- reimplemented INDEPENDENTLY
/// (MPL-2.0, no HYPRE source, ~400 lines on CompactAMG) so it carries NO HYPRE dependency and runs natively
/// on the NGSolve TaskManager.  The name "HypreBasedAMS" states the lineage honestly: the algorithm is
/// HYPRE's AMS; only the implementation is ours/HYPRE-free.  (Renamed 2026-06-27 from the misleading
/// "CompactAMS", which read like a distinct reduced variant; CompactAMSPreconditioner remains a back-compat
/// alias.)  Requires a lowest-order HCurl space (order=1, nograds=True): the Pi interpolation below is
/// built from the Whitney edge-vertex incidence.  Components (identical to HYPRE AMS):
///   - Gradient subspace: G^T * A_bc * G solved by CompactAMG
///   - Nodal subspace: Pi^T * A_bc * Pi solved by CompactAMG (component-wise)
///   - Fine-grid smoother: l1-Jacobi (fully TaskManager parallel)
///
/// All parallelism via NGSolve TaskManager (no OpenMP).
///
/// Reference: Hiptmair & Xu, SIAM J. Numer. Anal. 45(6), 2007.
///            Kolev & Vassilevski, J. Comput. Math. 27(5), 2009 (the HYPRE AMS algorithm).

#ifndef SPARSESOLV_COMPACT_AMS_HPP
#define SPARSESOLV_COMPACT_AMS_HPP

#include "compact_amg.hpp"
#include <comp.hpp>
#include <core/taskmanager.hpp>
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
///   auto ams = make_shared<HypreBasedAMS>(mat, grad, freedofs,
///                                      coord_x, coord_y, coord_z);
///   // Use as preconditioner with COCR or CG
class HypreBasedAMS : public BaseMatrix {
public:
    // Reject caller-owned regions before touching matrix state. Setup owns
    // bounded internal regions, excluding coarse direct factorization.
    static void RequireSerialSetup() {
        if (ngcore::GetTaskManager() != nullptr)
            throw std::runtime_error(
                "AMS construction and Update must run outside ngsolve.TaskManager; "
                "leave the TaskManager context before setup, then re-enter for the solve.");
    }

    /// @param mat       HCurl system matrix (SparseMatrix<double>)
    /// @param grad      Discrete gradient G (H1 -> HCurl)
    /// @param freedofs  Free DOFs for HCurl space
    /// @param coord_x/y/z  Vertex coordinates (length = ndof_h1)
    /// @param cycle_type   1=01210 (default), 7=0201020
    /// @param num_smooth   Fine-grid smoother sweeps (default=1)
    /// @param amg_theta    AMG strength threshold (default=0.25)
    /// @param print_level  Verbosity (0=silent)
    /// @param subspace_solver  0=CompactAMG (default), 1=SparseCholesky (diagnostic)
    /// @param beta_zero Omit gradient correction for compatible pure curl-curl systems.
    HypreBasedAMS(shared_ptr<SparseMatrix<double>> mat,
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
               int subspace_solver = 0,
               bool beta_zero = false,
               bool reuse_hierarchy = false,
               bool mixed_precision = false)
        : mat_(mat), grad_(grad), freedofs_(freedofs),
          ndof_hc_(mat->Height()), ndof_h1_(grad->Width()),
          cycle_type_(cycle_type), num_smooth_(num_smooth),
          print_level_(print_level), correction_weight_(correction_weight),
          subspace_solver_(subspace_solver), amg_theta_(amg_theta), beta_zero_(beta_zero),
          reuse_hierarchy_(reuse_hierarchy), mixed_precision_(mixed_precision)
    {
        RequireSerialSetup();
        if (mixed_precision_ && !beta_zero_)
            throw std::invalid_argument("HypreBasedAMS: mixed_precision requires beta_zero (float32 matrix values "
                "lose the gauge-scale gradient components that the G correction amplifies)");
        if ((int)coord_x.size() != ndof_h1_ ||
            (int)coord_y.size() != ndof_h1_ ||
            (int)coord_z.size() != ndof_h1_)
            throw std::runtime_error("HypreBasedAMS: coordinate size mismatch with H1 DOFs");
        RequireLowestOrderGradient();

        Setup(coord_x, coord_y, coord_z);
    }

    /// The auxiliary spaces are built from the lowest-order (Whitney) edge-vertex
    /// incidence: BuildPiComponents treats every row of G as an edge with two
    /// vertices.  An order>=2 space with nograds=True still has one gradient
    /// column per vertex, so the coordinate-size check above cannot see it, but
    /// its non-edge rows are empty: those dofs would get smoothing only and the
    /// preconditioner would silently degrade instead of failing.  Reject it.
    void RequireLowestOrderGradient() const {
        if (grad_->Height() != ndof_hc_)
            throw std::runtime_error(
                "HypreBasedAMS: the discrete gradient has " + std::to_string(grad_->Height())
                + " rows but the HCurl matrix has " + std::to_string(ndof_hc_)
                + "; build both from the same HCurl space.");
        int bad = 0;
        for (int e = 0; e < grad_->Height(); e++) {
            auto vals = grad_->GetRowValues(e);
            int nonzeros = 0;
            for (int j = 0; j < vals.Size(); j++)
                if (vals[j] != 0.0) nonzeros++;
            if (nonzeros != 2) bad++;
        }
        if (bad > 0)
            throw std::runtime_error(
                "HypreBasedAMS requires a lowest-order HCurl space (order=1, nograds=True): "
                + std::to_string(bad) + " of " + std::to_string(grad_->Height())
                + " discrete-gradient rows are not edge-vertex pairs. "
                  "For order >= 2 use NGSolve's bddc preconditioner.");
    }

    /// Update preconditioner with current matrix values.
    /// Geometry (G, Pi, transposes, work vectors) is preserved.
    /// Rebuilds: A_bc, Galerkin projections, AMG hierarchies, l1 norms.
    void Update() {
        RequireSerialSetup();
        mult_count_ = 0;
        t_smooth_ = t_grad_ = t_nodal_ = t_bc_ = 0;
        t_residual_ = t_restrict_ = t_auxiliary_ = t_prolong_ = 0;
        RebuildMatrix();
    }

    /// Update with a new system matrix, then rebuild.
    void Update(shared_ptr<SparseMatrix<double>> new_mat) {
        RequireSerialSetup();
        if (new_mat->Height() != mat_->Height() || new_mat->Width() != mat_->Width())
            throw std::invalid_argument("HypreBasedAMS::Update: new matrix dimension ("
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

    // Const accessors for fused Re/Im preconditioner (ComplexHypreBasedAMS)
    const SparseMatrix<double>& GetAbc() const { return *A_bc_; }
    const SparseMatrix<double>& GetGrad() const { return *grad_; }
    const SparseMatrix<double>& GetGradT() const { return *grad_t_; }
    const SparseMatrix<double>& GetPix() const { return *Pix_; }
    const SparseMatrix<double>& GetPiy() const { return *Piy_; }
    const SparseMatrix<double>& GetPiz() const { return *Piz_; }
    const SparseMatrix<double>& GetPixT() const { return *Pix_t_; }
    const SparseMatrix<double>& GetPiyT() const { return *Piy_t_; }
    const SparseMatrix<double>& GetPizT() const { return *Piz_t_; }
    const BaseMatrix& GetBG() const {
        if (!B_G_) throw std::runtime_error("AMS beta_zero has no gradient solver");
        return *B_G_;
    }
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
    bool GetBetaZero() const { return beta_zero_; }
    bool GetReuseHierarchy() const { return reuse_hierarchy_; }
    bool GetMixedPrecision() const { return mixed_precision_; }
    int GetHierarchyRefreshes() const { return hierarchy_refreshes_; }
    int GetInPlaceUpdates() const { return in_place_updates_; }
    int GetSetupWorkers() const { return setup_workers_; }
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

        // The matrix gateway may supply a freshly allocated parallel vector.
        // Every entry of x is written before it is read: the first smoother
        // assigns x = b0 / l1 (see FineSmooth), and without smoothing x is zeroed
        // here, so the cycle does not depend on BaseVector::SetScalar behavior
        // across concrete parallel-vector implementations.
        if (num_smooth_ < 1) x.FVDouble() = 0.0;

        // Zero constrained DOFs in RHS (one fused copy-and-mask pass)
        auto& b0 = *b0_;
        {
            auto src = b.FVDouble();
            auto dst = b0.FVDouble();
            if (freedofs_)
                ParallelFor(ndof_hc_, [&](size_t i) { dst[i] = freedofs_->Test(i) ? src[i] : 0.0; });
            else
                ParallelFor(ndof_hc_, [&](size_t i) { dst[i] = src[i]; });
        }
        t_bc_ += elapsed();

        if (cycle_type_ == 7) {
            FineSmooth(b0, x, true); t_smooth_ += elapsed();
            GradientCorrect(b0, x); t_grad_ += elapsed();
            FineSmooth(b0, x); t_smooth_ += elapsed();
            NodalCorrect(b0, x); t_nodal_ += elapsed();
            FineSmooth(b0, x); t_smooth_ += elapsed();
            GradientCorrect(b0, x); t_grad_ += elapsed();
            FineSmooth(b0, x); t_smooth_ += elapsed();
        } else if (num_smooth_ == 1) {
            // 01210 with one sweep: run the cycle in the work vector and let the
            // last smoother write x = w + (b0 - A w) / l1, masked, in one pass.
            auto& w = *w_;
            FineSmooth(b0, w, true); t_smooth_ += elapsed();
            GradientCorrect(b0, w); t_grad_ += elapsed();
            NodalCorrect(b0, w); t_nodal_ += elapsed();
            GradientCorrect(b0, w); t_grad_ += elapsed();
            FinalSmoothInto(b0, w, x); t_smooth_ += elapsed();
            mult_count_++;
            PrintMultBreakdown();
            return;
        } else {
            // 01210 (default, cycle_type=1)
            FineSmooth(b0, x, true); t_smooth_ += elapsed();
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
        PrintMultBreakdown();
    }

    void PrintMultBreakdown() const {
        if (print_level_ >= 1 && mult_count_ == 25) {
            std::cout << "\n  AMS Mult x" << mult_count_ << " breakdown:"
                      << " smooth=" << t_smooth_ << "s"
                      << " grad=" << t_grad_ << "s"
                      << " nodal=" << t_nodal_ << "s"
                      << " bc=" << t_bc_ << "s"
                      << " total=" << (t_smooth_+t_grad_+t_nodal_+t_bc_) << "s"
                      << std::endl;
            std::cout << "  AMS correction subphases: residual=" << t_residual_
                      << "s restrict=" << t_restrict_ << "s auxiliary=" << t_auxiliary_
                      << "s prolong=" << t_prolong_ << "s" << std::endl;
        }
    }

    /// x = w + (b - A_bc w) / l1 on free dofs, 0 on constrained ones, in one
    /// pass (the last l1-Jacobi sweep of the 01210 cycle written to the output).
    void FinalSmoothInto(const BaseVector& b, const BaseVector& w, BaseVector& x) const {
        auto fb = b.FVDouble();
        auto fw = w.FVDouble();
        auto fx = x.FVDouble();
        const auto& A = *A_bc_;
        const float* vf = mixed_precision_ ? abc_values_f_.data() : nullptr;
        ParallelForRange(ndof_hc_, [&](IntRange range) {
            for (auto i : range) {
                if (freedofs_ && !freedofs_->Test(i)) { fx[i] = 0.0; continue; }
                auto cols = A.GetRowIndices(i);
                double s = fb[i];
                if (vf) {
                    const float* v = vf + A.First(i);
                    for (int j = 0; j < cols.Size(); j++) s -= double(v[j]) * fw[cols[j]];
                } else {
                    auto vals = A.GetRowValues(i);
                    for (int j = 0; j < cols.Size(); j++) s -= vals[j] * fw[cols[j]];
                }
                fx[i] = fw[i] + s / fine_l1_norms_[i];
            }
        });
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
    const bool beta_zero_; // Pure curl-curl: omit G correction and its hierarchy.
    const bool reuse_hierarchy_; // Update(): frozen AMG coarsening, refreshed Galerkin matrices.
    const bool mixed_precision_; // Residual SpMVs inside the cycle read float32 value mirrors.
    std::vector<float> abc_values_f_;  // float32 mirror of A_bc_ (mixed precision only)
    std::vector<float> pi_values_f_, pit_values_f_;  // interleaved Pi / Pi^T values (mixed, fused)
    int hierarchy_refreshes_ = 0;
    int in_place_updates_ = 0;   // Updates whose auxiliary products ran on fixed patterns.
    int pi_shared_pattern_ = -1; // Pi_x/y/z, transposes and products share patterns (-1 unknown).
    bool pi_fused_ = false;      // Pi_x/y/z and their transposes share patterns (set in SetupGeometry).
    shared_ptr<SparseMatrix<double>> pattern_source_;  // matrix A_bc_ was last built from
    int setup_workers_ = 0;

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
    mutable std::unique_ptr<VVector<double>> w_;   // Cycle iterate (01210, one sweep)
    mutable std::unique_ptr<VVector<double>> r_G_, g_G_;      // Gradient space
    mutable std::unique_ptr<VVector<double>> r_Pix_, g_Pix_;  // Pix space
    mutable std::unique_ptr<VVector<double>> r_Piy_, g_Piy_;  // Piy space
    mutable std::unique_ptr<VVector<double>> r_Piz_, g_Piz_;  // Piz space

    // Accumulated timing (mutable for const Mult)
    mutable int mult_count_ = 0;
    mutable double t_smooth_ = 0, t_grad_ = 0, t_nodal_ = 0, t_bc_ = 0;
    mutable double t_residual_ = 0, t_restrict_ = 0, t_auxiliary_ = 0, t_prolong_ = 0;

    template <class F> void ProfileCorrection(double& seconds, F&& action) const {
        if (print_level_ == 0) { action(); return; }
        auto start = std::chrono::steady_clock::now();
        action();
        seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }

    // =====================================================================
    // Setup: split into geometry (one-time) + matrix (per-Update) + alloc
    // =====================================================================
    void Setup(const std::vector<double>& cx,
               const std::vector<double>& cy,
               const std::vector<double>& cz) {
        auto t_total = std::chrono::high_resolution_clock::now();
        if (print_level_ > 0)
            std::cout << "HypreBasedAMS setup: HC=" << ndof_hc_
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
            std::cout << "\n  HypreBasedAMS setup complete: "
                      << std::chrono::duration<double>(t_end - t_total).count()
                      << "s total" << std::endl;
    }

    // =====================================================================
    // SetupGeometry: Pi matrices + transposes (geometry-only, one-time)
    // =====================================================================
    void SetupGeometry(const std::vector<double>& cx,
                       const std::vector<double>& cy,
                       const std::vector<double>& cz) {
        // Parallel region (setup runs outside any caller TaskManager).
        ngcore::RegionTaskManager tasks;
        auto t_geometry = std::chrono::steady_clock::now();
        BuildPiComponents(cx, cy, cz);
        grad_t_ = dynamic_pointer_cast<SparseMatrix<double>>(grad_->CreateTranspose(true));
        // Pi_d has G's pattern, so Pi_d^T has G^T's: copy it and read each
        // value from the edge row of Pi_d (an exact transpose, one transposition).
        Pix_t_ = make_shared<SparseMatrix<double>>(*grad_t_);
        Piy_t_ = make_shared<SparseMatrix<double>>(*grad_t_);
        Piz_t_ = make_shared<SparseMatrix<double>>(*grad_t_);
        ParallelFor(ndof_h1_, [&](size_t v) {
            auto edges = grad_t_->GetRowIndices(v);
            auto tx = Pix_t_->GetRowValues(v);
            auto ty = Piy_t_->GetRowValues(v);
            auto tz = Piz_t_->GetRowValues(v);
            for (int j = 0; j < edges.Size(); j++) {
                const int e = edges[j];
                auto ends = Pix_->GetRowIndices(e);
                int k = 0;
                while (k < ends.Size() && ends[k] != int(v)) k++;
                tx[j] = Pix_->GetRowValues(e)[k];
                ty[j] = Piy_->GetRowValues(e)[k];
                tz[j] = Piz_->GetRowValues(e)[k];
            }
        });
        // Pi_x/y/z are built on G's pattern; when their transposes share one
        // pattern too, restriction and prolongation run as single sweeps.
        pi_fused_ = SameSparsityPattern(*Pix_, *Piy_) && SameSparsityPattern(*Pix_, *Piz_)
            && SameSparsityPattern(*Pix_t_, *Piy_t_) && SameSparsityPattern(*Pix_t_, *Piz_t_);
        if (mixed_precision_ && pi_fused_) {
            // Geometry-only: interleaved (x, y, z) float32 values per entry.
            auto interleave = [](const SparseMatrix<double>& X, const SparseMatrix<double>& Y,
                                 const SparseMatrix<double>& Z, std::vector<float>& out) {
                auto vx = X.AsVector().FVDouble(), vy = Y.AsVector().FVDouble(), vz = Z.AsVector().FVDouble();
                out.resize(3 * vx.Size());
                ParallelFor(vx.Size(), [&](size_t k) {
                    out[3 * k] = float(vx[k]); out[3 * k + 1] = float(vy[k]); out[3 * k + 2] = float(vz[k]);
                });
            };
            interleave(*Pix_, *Piy_, *Piz_, pi_values_f_);
            interleave(*Pix_t_, *Piy_t_, *Piz_t_, pit_values_f_);
        }
        if (print_level_ > 0)
            std::cout << "\n  Geometry (Pi, transposes): "
                      << std::chrono::duration<double>(std::chrono::steady_clock::now() - t_geometry).count()
                      << "s" << std::flush;
    }

    // =====================================================================
    // RebuildMatrix: A_bc, Galerkin projections, AMG, l1 norms
    // Called on every Update() — geometry is preserved.
    // =====================================================================
    void RebuildMatrix() {
        setup_workers_ = 0;
        auto t_rebuild = std::chrono::high_resolution_clock::now();

        // reuse_hierarchy with an unchanged sparsity pattern: refresh A_bc and
        // the auxiliary Galerkin matrices in place (numeric products only).
        auto t0 = std::chrono::high_resolution_clock::now();
        bool in_place = false;
        {
            ngcore::RegionTaskManager tasks;
            // A SparseMatrix never changes its pattern, so the matrix object
            // A_bc was last built from needs no comparison.
            in_place = reuse_hierarchy_ && A_bc_ && A_Pix_ && A_Piy_ && A_Piz_
                && (beta_zero_ || A_G_)
                && (mat_ == pattern_source_ || SameSparsityPattern(*mat_, *A_bc_));
            if (in_place) {
                if (freedofs_) FillBCModifiedValues(*mat_, *freedofs_, *A_bc_);
                else A_bc_ = mat_;
                if (pi_shared_pattern_ < 0)
                    pi_shared_pattern_ = SameSparsityPattern(*Pix_, *Piy_) && SameSparsityPattern(*Pix_, *Piz_)
                        && SameSparsityPattern(*Pix_t_, *Piy_t_) && SameSparsityPattern(*Pix_t_, *Piz_t_)
                        && SameSparsityPattern(*A_Pix_, *A_Piy_) && SameSparsityPattern(*A_Pix_, *A_Piz_);
                const bool nodal = pi_shared_pattern_
                    ? FusedPiProducts()
                    : GalerkinProductOnPattern(*Pix_t_, *A_bc_, *Pix_, *A_Pix_)
                      && GalerkinProductOnPattern(*Piy_t_, *A_bc_, *Piy_, *A_Piy_)
                      && GalerkinProductOnPattern(*Piz_t_, *A_bc_, *Piz_, *A_Piz_);
                in_place = nodal
                    && (beta_zero_ || GalerkinProductOnPattern(*grad_t_, *A_bc_, *grad_, *A_G_));
                if (in_place) in_place_updates_++;
            }
        }
        if (!in_place) pi_shared_pattern_ = -1;  // the products below get new patterns
        pattern_source_ = mat_;
        if (!in_place) {
            ngcore::RegionTaskManager tasks;
            // 1. Create BC-modified matrix (identity rows for constrained DOFs).
            if (freedofs_) {
                A_bc_ = CreateBCModifiedMatrix(*mat_, *freedofs_);
            } else {
                A_bc_ = mat_;
            }

            // 2. Galerkin projection: A_G = G^T * A_bc * G, etc.
            if (reuse_hierarchy_ && pi_fused_) {
                // Native: one symbolic pass for the shared Pi pattern, then the
                // fused numeric pass for all three components.
                // For an element-graph A the pattern of Pi^T A Pi is the vertex
                // graph: every pair of tet vertices is an edge, so it is each
                // vertex and its edge neighbours (read off Pi^T). Any other A
                // falls back to the full symbolic product.
                A_Pix_ = VertexEdgePattern();
                A_Piy_ = make_shared<SparseMatrix<double>>(*A_Pix_);
                A_Piz_ = make_shared<SparseMatrix<double>>(*A_Pix_);
                pi_shared_pattern_ = 1;
                if (!FusedPiProducts()) {
                    A_Pix_ = GalerkinPattern(*Pix_t_, *A_bc_, *Pix_);
                    A_Piy_ = make_shared<SparseMatrix<double>>(*A_Pix_);
                    A_Piz_ = make_shared<SparseMatrix<double>>(*A_Pix_);
                    if (!FusedPiProducts())
                        throw std::runtime_error("HypreBasedAMS: Pi product outside its own pattern");
                }
                if (!beta_zero_) A_G_ = GalerkinProduct(*grad_t_, *A_bc_, *grad_);
            } else {
                if (!beta_zero_)
                    A_G_ = dynamic_pointer_cast<SparseMatrix<double>>(A_bc_->Restrict(*grad_));
                A_Pix_ = dynamic_pointer_cast<SparseMatrix<double>>(A_bc_->Restrict(*Pix_));
                A_Piy_ = dynamic_pointer_cast<SparseMatrix<double>>(A_bc_->Restrict(*Piy_));
                A_Piz_ = dynamic_pointer_cast<SparseMatrix<double>>(A_bc_->Restrict(*Piz_));
            }
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        double dt_restrict = std::chrono::duration<double>(t1 - t0).count();

        if (print_level_ > 0) {
            std::cout << (!in_place ? "\n  Galerkin restrict: "
                          : pi_shared_pattern_ == 1 ? "\n  Galerkin refresh (fixed pattern, fused Pi): "
                                                    : "\n  Galerkin refresh (fixed pattern): ")
                      << dt_restrict << "s";
            if (beta_zero_) std::cout << " beta_zero: gradient hierarchy omitted";
            else std::cout << " A_G=" << A_G_->Height() << "x" << A_G_->Width()
                           << " nnz=" << A_G_->NZE();
            std::cout << std::flush;
        }

        // 3. Fix zero rows: set diag = 1.0 for truly zero rows.
        if (!beta_zero_) FixZeroRows(*A_G_);
        FixZeroRows(*A_Pix_);
        FixZeroRows(*A_Piy_);
        FixZeroRows(*A_Piz_);

        // 4. Build solvers for auxiliary spaces
        if (subspace_solver_ == 1) {
            if (print_level_ > 0)
                std::cout << "\n  Subspace solver: SparseCholesky (diagnostic)" << std::flush;

            t0 = std::chrono::high_resolution_clock::now();
            if (!beta_zero_) {
                A_G_->SetInverseType("sparsecholesky");
                B_G_ = A_G_->InverseMatrix(shared_ptr<BitArray>(nullptr));
            }
            t1 = std::chrono::high_resolution_clock::now();
            if (print_level_ > 0 && !beta_zero_)
                std::cout << "\n  B_G setup: " << std::chrono::duration<double>(t1 - t0).count()
                          << "s (direct, n=" << A_G_->Height() << ")" << std::flush;

            t0 = std::chrono::high_resolution_clock::now();
            A_Pix_->SetInverseType("sparsecholesky");
            A_Piy_->SetInverseType("sparsecholesky");
            A_Piz_->SetInverseType("sparsecholesky");
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

            // reuse_hierarchy: keep the coarsening and interpolation of the
            // first build and refresh only the Galerkin coarse matrices.
            auto existing = [](const shared_ptr<BaseMatrix>& solver) {
                return dynamic_pointer_cast<CompactAMG>(solver);
            };
            const bool refresh = reuse_hierarchy_ && existing(B_Pix_) && existing(B_Piy_)
                && existing(B_Piz_) && (beta_zero_ || existing(B_G_));
            shared_ptr<CompactAMG> amg_G, amg_Pix, amg_Piy, amg_Piz;
            if (refresh) {
                if (!beta_zero_) amg_G = existing(B_G_);
                amg_Pix = existing(B_Pix_);
                amg_Piy = existing(B_Piy_);
                amg_Piz = existing(B_Piz_);
            } else {
                if (!beta_zero_)
                    amg_G = make_shared<CompactAMG>(A_G_, nullptr, amg_theta_, 25, min_coarse_aux, 1, 0);
                amg_Pix = make_shared<CompactAMG>(A_Pix_, nullptr, amg_theta_, 25, min_coarse_aux, 1, 0);
                amg_Piy = make_shared<CompactAMG>(A_Piy_, nullptr, amg_theta_, 25, min_coarse_aux, 1, 0);
                amg_Piz = make_shared<CompactAMG>(A_Piz_, nullptr, amg_theta_, 25, min_coarse_aux, 1, 0);
                // reuse_hierarchy: coarse matrices from the native Galerkin product
                // (the refresh path's numeric pass then runs on those patterns).
                for (auto& amg : {amg_G, amg_Pix, amg_Piy, amg_Piz})
                    if (amg) { amg->SetNativeGalerkin(reuse_hierarchy_); amg->SetMixedPrecision(mixed_precision_); }
            }

            t0 = std::chrono::high_resolution_clock::now();
            if (refresh) {
                if (!beta_zero_) amg_G->Refresh(A_G_);
                amg_Pix->Refresh(A_Pix_);
                amg_Piy->Refresh(A_Piy_);
                amg_Piz->Refresh(A_Piz_);
                hierarchy_refreshes_++;
            } else {
                if (!beta_zero_) amg_G->Setup();
                amg_Pix->Setup();
                amg_Piy->Setup();
                amg_Piz->Setup();
            }
            setup_workers_ = std::max({beta_zero_ ? 0 : amg_G->SetupWorkers(),
                amg_Pix->SetupWorkers(), amg_Piy->SetupWorkers(), amg_Piz->SetupWorkers()});
            t1 = std::chrono::high_resolution_clock::now();
            if (print_level_ > 0)
                std::cout << (refresh ? "\n  AMG refresh (frozen coarsening and interpolation): "
                                      : "\n  AMG setup (internally parallel hierarchy, serial coarse factorization): ")
                          << std::chrono::duration<double>(t1 - t0).count()
                          << "s, levels: G=" << (amg_G ? amg_G->NumLevels() : 0)
                          << " Px=" << amg_Pix->NumLevels()
                          << " Py=" << amg_Piy->NumLevels()
                          << " Pz=" << amg_Piz->NumLevels() << std::flush;

            if (print_level_ > 0 && !refresh) {
                std::array<double, 7> sum{};
                for (auto& amg : {amg_G, amg_Pix, amg_Piy, amg_Piz})
                    if (amg) for (int k = 0; k < 7; k++) sum[k] += amg->SetupPhases()[k];
                std::cout << "\n  AMG setup phases: strength=" << sum[0] << "s coarsen=" << sum[1]
                          << "s interp=" << sum[2] << "s transpose=" << sum[3] << "s galerkin=" << sum[4]
                          << "s l1+alloc=" << sum[5] << "s coarse_factor=" << sum[6] << "s" << std::flush;
            }
            B_G_ = amg_G;
            B_Pix_ = amg_Pix;
            B_Piy_ = amg_Piy;
            B_Piz_ = amg_Piz;
        }

        // 5. Float mirror and truncated l1 norms for the fine-grid smoother,
        // in a parallel region (the coarse factorizations above stay outside).
        {
            ngcore::RegionTaskManager tasks;
            if (mixed_precision_) MirrorValues(*A_bc_, abc_values_f_);
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
        w_ = std::make_unique<VVector<double>>(ndof_hc_);
        if (!beta_zero_) {
            r_G_ = std::make_unique<VVector<double>>(ndof_h1_);
            g_G_ = std::make_unique<VVector<double>>(ndof_h1_);
        }
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

    /// Zero matrix on the vertex graph read off Pi^T (= |G|^T): row v holds v
    /// and the other endpoint of each of its edges, sorted.
    shared_ptr<SparseMatrix<double>> VertexEdgePattern() const {
        Array<int> counts(ndof_h1_);
        ParallelFor(ndof_h1_, [&](size_t v) { counts[v] = int(Pix_t_->GetRowIndices(v).Size()) + 1; });
        auto C = make_shared<SparseMatrix<double>>(counts, ndof_h1_);
        ParallelFor(ndof_h1_, [&](size_t v) {
            auto cols = C->GetRowIndices(v);
            auto edges = Pix_t_->GetRowIndices(v);
            cols[0] = int(v);
            for (int a = 0; a < edges.Size(); a++) {
                auto ends = grad_->GetRowIndices(edges[a]);
                cols[a + 1] = ends[0] == int(v) ? ends[1] : ends[0];
            }
            std::sort(cols.Data(), cols.Data() + cols.Size());
            C->GetRowValues(v) = 0.0;
        });
        return C;
    }

    /// Pi_d^T A_bc Pi_d for d = x, y, z in one sweep of A_bc, on the existing
    /// patterns of A_Pix/y/z. The three Pi share G's pattern, so their
    /// transposes and products share theirs (checked by the caller); only the
    /// values differ. Returns false (values unspecified) when a nonzero
    /// contribution falls outside the pattern.
    bool FusedPiProducts() {
        auto& Cx = *A_Pix_;
        auto& Cy = *A_Piy_;
        auto& Cz = *A_Piz_;
        std::atomic<bool> inside{true};
        const int width = Cx.Width();
        ParallelForRange(Cx.Height(), [&](IntRange range) {
            std::vector<int> marker(width, -1);
            for (auto u : range) {
                auto ccols = Cx.GetRowIndices(u);
                auto cx = Cx.GetRowValues(u);
                auto cy = Cy.GetRowValues(u);
                auto cz = Cz.GetRowValues(u);
                for (int p = 0; p < ccols.Size(); p++) {
                    marker[ccols[p]] = p;
                    cx[p] = cy[p] = cz[p] = 0.0;
                }
                auto ecols = Pix_t_->GetRowIndices(u);
                auto tx = Pix_t_->GetRowValues(u);
                auto ty = Piy_t_->GetRowValues(u);
                auto tz = Piz_t_->GetRowValues(u);
                for (int a = 0; a < ecols.Size(); a++) {
                    const int e = ecols[a];
                    auto kcols = A_bc_->GetRowIndices(e);
                    auto kvals = A_bc_->GetRowValues(e);
                    for (int b = 0; b < kcols.Size(); b++) {
                        const double w = kvals[b];
                        if (w == 0.0) continue;
                        const double wx = tx[a] * w, wy = ty[a] * w, wz = tz[a] * w;
                        const int k = kcols[b];
                        auto jcols = Pix_->GetRowIndices(k);
                        auto px = Pix_->GetRowValues(k);
                        auto py = Piy_->GetRowValues(k);
                        auto pz = Piz_->GetRowValues(k);
                        for (int c = 0; c < jcols.Size(); c++) {
                            const int pos = marker[jcols[c]];
                            if (pos < 0) {
                                if (wx * px[c] != 0.0 || wy * py[c] != 0.0 || wz * pz[c] != 0.0)
                                    inside.store(false, std::memory_order_relaxed);
                                continue;
                            }
                            cx[pos] += wx * px[c];
                            cy[pos] += wy * py[c];
                            cz[pos] += wz * pz[c];
                        }
                    }
                }
                for (int p = 0; p < ccols.Size(); p++) marker[ccols[p]] = -1;
            }
        });
        return inside.load();
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

    /// res = b - A_bc x in one pass over the rows (no separate copy of b); with
    /// mixed precision A_bc's values come from its float32 mirror.
    void ResidualInto(const BaseVector& b, const BaseVector& x, BaseVector& res) const {
        ResidualWithValues(*A_bc_, mixed_precision_ ? abc_values_f_.data() : nullptr, b, x, res);
    }

    void FineSmooth(const BaseVector& b, BaseVector& x, bool initially_zero = false) const {
        // l1-Jacobi: fully parallel, no data dependency between rows.
        // Each sweep: compute residual r = b - A*x, then x += r / l1_norm.
        auto& res = *r0_;
        auto fv_x = x.FVDouble();

        for (int s = 0; s < num_smooth_; s++) {
            if (initially_zero && s == 0) {
                // x = 0 before the first sweep: the sweep is x = b / l1 (x is
                // assigned, so its previous content is never read).
                auto fv_b = b.FVDouble();
                ParallelFor(ndof_hc_, [&](size_t i) { fv_x[i] = fv_b[i] / fine_l1_norms_[i]; });
                continue;
            }
            ResidualInto(b, x, res);
            // Jacobi update: x[i] += r[i] / l1_norm[i] (fully parallel)
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
        ProfileCorrection(t_restrict_, [&]() { Pt.Mult(residual, r_c); });

        // Coarse solve: g_c = B^{-1} * r_c (one AMG V-cycle)
        g_c.FVDouble() = 0;
        ProfileCorrection(t_auxiliary_, [&]() { B.Mult(r_c, g_c); });

        // Prolongate and add: x += omega * P * g_c
        ProfileCorrection(t_prolong_, [&]() { P.MultAdd(correction_weight_, g_c, x); });
    }

    void ComputeResidual(const BaseVector& b, const BaseVector& x, BaseVector& res) const {
        ProfileCorrection(t_residual_, [&]() { ResidualInto(b, x, res); });
    }

    void GradientCorrect(const BaseVector& b, BaseVector& x) const {
        if (beta_zero_) return;
        ComputeResidual(b, x, *r0_);
        SubspaceCorrect(*grad_, *grad_t_, *B_G_, *r0_, x, *r_G_, *g_G_, "G");
    }

    void NodalCorrect(const BaseVector& b, BaseVector& x) const {
        ComputeResidual(b, x, *r0_);

        // Additive Pi corrections: all 3 use the same residual (no recomputation).
        // Saves 2 fine-level SpMVs per AMS cycle vs multiplicative approach.

        // Restrict to all 3 subspaces from same residual (one sweep when the
        // three Pi share a pattern, as built)
        ProfileCorrection(t_restrict_, [&]() {
        if (pi_fused_) {
            auto r = r0_->FVDouble();
            auto rx = r_Pix_->FVDouble();
            auto ry = r_Piy_->FVDouble();
            auto rz = r_Piz_->FVDouble();
            const float* tf = pit_values_f_.empty() ? nullptr : pit_values_f_.data();
            ParallelForRange(ndof_h1_, [&](IntRange range) {
                for (auto v : range) {
                    auto cols = Pix_t_->GetRowIndices(v);
                    double sx = 0.0, sy = 0.0, sz = 0.0;
                    if (tf) {
                        const float* t = tf + 3 * Pix_t_->First(v);
                        for (int j = 0; j < cols.Size(); j++) {
                            const double re = r[cols[j]];
                            sx += double(t[3 * j]) * re; sy += double(t[3 * j + 1]) * re;
                            sz += double(t[3 * j + 2]) * re;
                        }
                    } else {
                        auto tx = Pix_t_->GetRowValues(v);
                        auto ty = Piy_t_->GetRowValues(v);
                        auto tz = Piz_t_->GetRowValues(v);
                        for (int j = 0; j < cols.Size(); j++) {
                            const double re = r[cols[j]];
                            sx += tx[j] * re; sy += ty[j] * re; sz += tz[j] * re;
                        }
                    }
                    rx[v] = sx; ry[v] = sy; rz[v] = sz;
                }
            });
        } else {
            Pix_t_->Mult(*r0_, *r_Pix_);
            Piy_t_->Mult(*r0_, *r_Piy_);
            Piz_t_->Mult(*r0_, *r_Piz_);
        }
        });

        // Solve all 3 (sequential: each AMG uses TaskManager internally)
        ProfileCorrection(t_auxiliary_, [&]() {
        // CompactAMG::Mult initializes its output; other subspace solvers
        // (direct inverses) assign it as well.
        // Sequential on purpose: each AMG V-cycle is internally parallel, and
        // running the three inside one ParallelFor (nested jobs) livelocks.
        B_Pix_->Mult(*r_Pix_, *g_Pix_);
        B_Piy_->Mult(*r_Piy_, *g_Piy_);
        B_Piz_->Mult(*r_Piz_, *g_Piz_);
        });

        // Prolongate and add all 3 corrections (one sweep when fused)
        ProfileCorrection(t_prolong_, [&]() {
        if (pi_fused_) {
            auto fx = x.FVDouble();
            auto gx = g_Pix_->FVDouble();
            auto gy = g_Piy_->FVDouble();
            auto gz = g_Piz_->FVDouble();
            const double w = correction_weight_;
            const float* pf = pi_values_f_.empty() ? nullptr : pi_values_f_.data();
            ParallelForRange(ndof_hc_, [&](IntRange range) {
                for (auto e : range) {
                    auto cols = Pix_->GetRowIndices(e);
                    double s = 0.0;
                    if (pf) {
                        const float* p = pf + 3 * Pix_->First(e);
                        for (int j = 0; j < cols.Size(); j++) {
                            const int v = cols[j];
                            s += double(p[3 * j]) * gx[v] + double(p[3 * j + 1]) * gy[v]
                                 + double(p[3 * j + 2]) * gz[v];
                        }
                    } else {
                        auto px = Pix_->GetRowValues(e);
                        auto py = Piy_->GetRowValues(e);
                        auto pz = Piz_->GetRowValues(e);
                        for (int j = 0; j < cols.Size(); j++) {
                            const int v = cols[j];
                            s += px[j] * gx[v] + py[j] * gy[v] + pz[j] * gz[v];
                        }
                    }
                    fx[e] += w * s;
                }
            });
        } else {
            Pix_->MultAdd(correction_weight_, *g_Pix_, x);
            Piy_->MultAdd(correction_weight_, *g_Piy_, x);
            Piz_->MultAdd(correction_weight_, *g_Piz_, x);
        }
        });
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
