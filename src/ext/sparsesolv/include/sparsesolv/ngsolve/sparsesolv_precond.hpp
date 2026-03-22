/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/// @file sparsesolv_precond.hpp
/// @brief SparseSolv preconditioners (IC, SGS) as NGSolve BaseMatrix wrappers

#ifndef NGSOLVE_SPARSESOLV_PRECOND_HPP
#define NGSOLVE_SPARSESOLV_PRECOND_HPP

#include <basematrix.hpp>
#include <sparsematrix.hpp>

// SparseSolv headers
#include "../sparsesolv.hpp"

namespace ngla {

template<typename SCAL>
const SCAL* GetVectorData(const BaseVector& vec) {
    return vec.FV<SCAL>().Data();
}

template<typename SCAL>
SCAL* GetVectorData(BaseVector& vec) {
    return vec.FV<SCAL>().Data();
}

// ============================================================================
// Common: NGSolve SparseMatrix → SparseSolv view conversion
// ============================================================================

/// Convert NGSolve SparseMatrix to SparseSolv view (identity rows for constrained DOFs)
template<typename SCAL>
sparsesolv::SparseMatrixView<SCAL> BuildSparseMatrixView(
    const SparseMatrix<SCAL>& mat,
    const BitArray* freedofs,
    sparsesolv::index_t height, sparsesolv::index_t width,
    std::vector<sparsesolv::index_t>& row_ptr,
    std::vector<sparsesolv::index_t>& col_idx,
    std::vector<SCAL>& modified_values)
{
    const auto& firsti = mat.GetFirstArray();
    const auto& colnr = mat.GetColIndices();
    const auto& values = mat.GetValues();

    const size_t nrows = firsti.Size();
    const size_t nnz = colnr.Size();

    row_ptr.resize(nrows);
    sparsesolv::parallel_for(static_cast<sparsesolv::index_t>(nrows), [&](sparsesolv::index_t i) {
        row_ptr[i] = static_cast<sparsesolv::index_t>(firsti[i]);
    });

    col_idx.resize(nnz);
    sparsesolv::parallel_for(static_cast<sparsesolv::index_t>(nnz), [&](sparsesolv::index_t i) {
        col_idx[i] = static_cast<sparsesolv::index_t>(colnr[i]);
    });

    if (freedofs) {
        modified_values.resize(values.Size());
        sparsesolv::parallel_for(static_cast<sparsesolv::index_t>(values.Size()), [&](sparsesolv::index_t i) {
            modified_values[i] = values[i];
        });
        sparsesolv::parallel_for(height, [&](sparsesolv::index_t i) {
            for (sparsesolv::index_t k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                sparsesolv::index_t j = col_idx[k];
                if (!freedofs->Test(i)) {
                    modified_values[k] = (j == i) ? SCAL(1) : SCAL(0);
                } else if (!freedofs->Test(j)) {
                    modified_values[k] = SCAL(0);
                }
            }
        });
        return sparsesolv::SparseMatrixView<SCAL>(
            height, width, row_ptr.data(), col_idx.data(), modified_values.data());
    } else {
        return sparsesolv::SparseMatrixView<SCAL>(
            height, width, row_ptr.data(), col_idx.data(), values.Data());
    }
}

// ============================================================================
// Preconditioner Base Class
// ============================================================================

/// Base class: NGSolve SparseMatrix → SparseSolv view, FreeDofs handling
template<typename SCAL>
class SparseSolvPrecondBase : public BaseMatrix {
protected:
    shared_ptr<SparseMatrix<SCAL>> mat_;
    shared_ptr<BitArray> freedofs_;
    sparsesolv::index_t height_;
    sparsesolv::index_t width_;

    // CSR conversion buffers
    std::vector<sparsesolv::index_t> row_ptr_;
    std::vector<sparsesolv::index_t> col_idx_;
    std::vector<SCAL> modified_values_;

    // Cached temp buffer for MultAdd (avoids per-call allocation)
    mutable std::unique_ptr<SCAL[]> cached_temp_;
    mutable sparsesolv::index_t cached_temp_size_ = 0;

    SparseSolvPrecondBase(shared_ptr<SparseMatrix<SCAL>> mat,
                          shared_ptr<BitArray> freedofs)
        : mat_(mat)
        , freedofs_(freedofs ? make_shared<BitArray>(*freedofs) : nullptr)
        , height_(static_cast<sparsesolv::index_t>(mat->Height()))
        , width_(static_cast<sparsesolv::index_t>(mat->Width()))
    {}

    sparsesolv::SparseMatrixView<SCAL> prepare_matrix_view() {
        return BuildSparseMatrixView<SCAL>(
            *mat_, freedofs_.get(), height_, width_,
            row_ptr_, col_idx_, modified_values_);
    }

    /// Apply the underlying SparseSolv preconditioner (implemented by derived classes)
    virtual void apply_precond(const SCAL* x, SCAL* y) const = 0;

public:
    void Mult(const BaseVector& x, BaseVector& y) const override {
        const SCAL* x_data = GetVectorData<SCAL>(x);
        SCAL* y_data = GetVectorData<SCAL>(y);

        apply_precond(x_data, y_data);

        if (freedofs_) {
            sparsesolv::parallel_for(height_, [&](sparsesolv::index_t i) {
                if (!freedofs_->Test(i)) {
                    y_data[i] = SCAL(0);
                }
            });
        }
    }

    void MultAdd(double s, const BaseVector& x, BaseVector& y) const override {
        const SCAL* x_data = GetVectorData<SCAL>(x);
        SCAL* y_data = GetVectorData<SCAL>(y);

        if (!cached_temp_ || cached_temp_size_ != height_) {
            cached_temp_.reset(new SCAL[height_]);
            cached_temp_size_ = height_;
        }
        SCAL* temp = cached_temp_.get();
        apply_precond(x_data, temp);

        if (freedofs_) {
            sparsesolv::parallel_for(height_, [&](sparsesolv::index_t i) {
                if (freedofs_->Test(i)) {
                    y_data[i] += s * temp[i];
                }
            });
        } else {
            sparsesolv::parallel_for(height_, [&](sparsesolv::index_t i) {
                y_data[i] += s * temp[i];
            });
        }
    }

    int VHeight() const override { return static_cast<int>(height_); }
    int VWidth() const override { return static_cast<int>(width_); }

    AutoVector CreateRowVector() const override {
        return mat_->CreateRowVector();
    }

    AutoVector CreateColVector() const override {
        return mat_->CreateColVector();
    }
};

// ============================================================================
// IC Preconditioner
// ============================================================================
template<typename SCAL = double>
class SparseSolvICPreconditioner : public SparseSolvPrecondBase<SCAL> {
public:
    SparseSolvICPreconditioner(shared_ptr<SparseMatrix<SCAL>> mat,
                                shared_ptr<BitArray> freedofs = nullptr,
                                double shift = 1.05)
        : SparseSolvPrecondBase<SCAL>(mat, freedofs)
        , shift_(shift)
        , precond_(std::make_shared<sparsesolv::ICPreconditioner<SCAL>>(shift))
    {}

    void Update() {
        auto view = this->prepare_matrix_view();
        precond_->setup(view);
    }

    double GetShift() const { return shift_; }

    void SetShift(double shift) {
        shift_ = shift;
        precond_->set_shift_parameter(shift);
    }

    // ABMC accessors
    bool GetUseABMC() const { return config_.use_abmc; }
    void SetUseABMC(bool v) { config_.use_abmc = v; sync_config(); }
    int GetABMCBlockSize() const { return config_.abmc_block_size; }
    void SetABMCBlockSize(int v) { config_.abmc_block_size = v; sync_config(); }
    int GetABMCNumColors() const { return config_.abmc_num_colors; }
    void SetABMCNumColors(int v) { config_.abmc_num_colors = v; sync_config(); }
    bool GetDiagonalScaling() const { return config_.diagonal_scaling; }
    void SetDiagonalScaling(bool v) { config_.diagonal_scaling = v; sync_config(); }

protected:
    void apply_precond(const SCAL* x, SCAL* y) const override {
        precond_->apply(x, y, this->height_);
    }

private:
    double shift_;
    sparsesolv::SolverConfig config_;
    shared_ptr<sparsesolv::ICPreconditioner<SCAL>> precond_;

    void sync_config() {
        config_.shift_parameter = shift_;
        precond_->set_config(config_);
    }
};

// ============================================================================
// SGS Preconditioner
// ============================================================================
template<typename SCAL = double>
class SparseSolvSGSPreconditioner : public SparseSolvPrecondBase<SCAL> {
public:
    SparseSolvSGSPreconditioner(shared_ptr<SparseMatrix<SCAL>> mat,
                                 shared_ptr<BitArray> freedofs = nullptr)
        : SparseSolvPrecondBase<SCAL>(mat, freedofs)
        , precond_(std::make_shared<sparsesolv::SGSPreconditioner<SCAL>>())
    {}

    void Update() {
        auto view = this->prepare_matrix_view();
        precond_->setup(view);
    }

protected:
    void apply_precond(const SCAL* x, SCAL* y) const override {
        precond_->apply(x, y, this->height_);
    }

private:
    shared_ptr<sparsesolv::SGSPreconditioner<SCAL>> precond_;
};

// Type aliases for convenience
using ICPreconditioner = SparseSolvICPreconditioner<double>;
using SGSPreconditioner = SparseSolvSGSPreconditioner<double>;

} // namespace ngla

#endif // NGSOLVE_SPARSESOLV_PRECOND_HPP
