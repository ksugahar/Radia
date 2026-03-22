/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file sgs_preconditioner.hpp
 * @brief Symmetric Gauss-Seidel (SGS) preconditioner
 */

#ifndef SPARSESOLV_PRECONDITIONERS_SGS_PRECONDITIONER_HPP
#define SPARSESOLV_PRECONDITIONERS_SGS_PRECONDITIONER_HPP

#include "../core/preconditioner.hpp"
#include "../core/constants.hpp"
#include "../core/level_schedule.hpp"
#include <algorithm>
#include <cmath>

namespace sparsesolv {

/**
 * @brief Symmetric Gauss-Seidel (SGS) preconditioner
 *
 * Computes the SGS preconditioner: M = (D + L) * D^{-1} * (D + U)
 * where A = L + D + U (lower, diagonal, upper decomposition).
 *
 * The application M^{-1} * x is computed as:
 * 1. Forward sweep: solve (D + L) * y = x
 * 2. Diagonal scaling: z = D * y
 * 3. Backward sweep: solve (D + U) * result = z
 *
 * SGS is a good general-purpose preconditioner that:
 * - Requires no additional fill-in
 * - Works well with the SGS-MRTR solver
 * - Is symmetric, so can be used with CG
 *
 * @note This class stores non-owning pointers to the matrix CSR data.
 *       The caller MUST ensure the underlying data remains valid for the
 *       lifetime of this preconditioner (e.g., by holding a shared_ptr
 *       to the source matrix).
 * @note This class is NOT thread-safe for concurrent apply() calls.
 *       Use separate preconditioner instances for each thread.
 *
 * @tparam Scalar The scalar type (double or complex<double>)
 */
template<typename Scalar = double>
class SGSPreconditioner : public Preconditioner<Scalar> {
public:
    SGSPreconditioner() = default;

    // Disable copy (level schedules are non-trivial)
    SGSPreconditioner(const SGSPreconditioner&) = delete;
    SGSPreconditioner& operator=(const SGSPreconditioner&) = delete;

    // Enable move
    SGSPreconditioner(SGSPreconditioner&&) noexcept = default;
    SGSPreconditioner& operator=(SGSPreconditioner&&) noexcept = default;

    /**
     * @brief Setup the SGS preconditioner from matrix A
     *
     * Stores non-owning pointers to A's CSR data (no deep copy).
     * The caller must ensure the data outlives this preconditioner.
     *
     * @param A Sparse matrix view (CSR format, data must remain valid)
     */
    void setup(const SparseMatrixView<Scalar>& A) override {
        const index_t n = A.rows();
        size_ = n;

        // Store non-owning pointers (no deep copy — SGS does not modify values)
        row_ptr_ = A.row_ptr();
        col_idx_ = A.col_idx();
        values_ = A.values();

        // Extract and invert diagonal
        inv_diag_.resize(n);
        parallel_for(n, [&](index_t i) {
            Scalar d = A.diagonal(i);
            if (std::abs(d) > constants::MIN_DIAGONAL_TOLERANCE) {
                inv_diag_[i] = Scalar(1) / d;
            } else {
                inv_diag_[i] = Scalar(1);
            }
        });

        // Build level schedules for parallel triangular solves
        fwd_schedule_.build_from_lower(row_ptr_, col_idx_, n);
        bwd_schedule_.build_from_upper(row_ptr_, col_idx_, n);

        // Pre-allocate work vectors for apply()
        temp_.resize(size_);
        scaled_.resize(size_);

        this->is_setup_ = true;
    }

    /**
     * @brief Apply SGS preconditioner: y = M^{-1} * x
     *
     * M^{-1} = (D + U)^{-1} * D * (D + L)^{-1}
     *
     * @param x Input vector
     * @param y Output vector
     * @param size Vector size
     *
     * @note NOT thread-safe. Use separate instances for concurrent access.
     */
    void apply(const Scalar* x, Scalar* y, index_t size) const override {
        if (!this->is_setup_) {
            throw std::runtime_error("SGSPreconditioner::apply called before setup");
        }

        // Step 1: Forward Gauss-Seidel sweep (sequential - data dependency)
        // Solve (D + L) * temp = x  =>  temp = (D + L)^{-1} * x
        forward_sweep(x, temp_.data());

        // Step 2: Diagonal scaling (no data dependency)
        // scaled = D * temp
        parallel_for(size_, [&](index_t i) {
            // inv_diag_[i] = 1/D[i,i], so D[i,i] = 1/inv_diag_[i]
            scaled_[i] = temp_[i] / inv_diag_[i];
        });

        // Step 3: Backward Gauss-Seidel sweep (sequential - data dependency)
        // Solve (D + U) * y = scaled  =>  y = (D + U)^{-1} * D * (D + L)^{-1} * x
        backward_sweep(scaled_.data(), y);
    }

    std::string name() const override { return "SGS"; }

private:
    index_t size_ = 0;

    // Non-owning pointers to CSR data (caller ensures lifetime)
    const index_t* row_ptr_ = nullptr;
    const index_t* col_idx_ = nullptr;
    const Scalar* values_ = nullptr;

    std::vector<Scalar> inv_diag_;

    // Level schedules for parallel triangular solves
    LevelSchedule fwd_schedule_;     // For forward sweep
    LevelSchedule bwd_schedule_;     // For backward sweep

    // Mutable work vectors (allocated once in setup, reused in apply)
    mutable std::vector<Scalar> temp_;
    mutable std::vector<Scalar> scaled_;

    /**
     * @brief Forward Gauss-Seidel sweep: solve (D + L) * y = x
     *
     * Uses level scheduling for parallelism.
     */
    void forward_sweep(const Scalar* x, Scalar* y) const {
        const index_t* rp = row_ptr_;
        const index_t* ci = col_idx_;
        const Scalar* v = values_;

        for (const auto& level : fwd_schedule_.levels) {
            const index_t level_size = static_cast<index_t>(level.size());
            parallel_for(level_size, [&](index_t idx) {
                const index_t i = level[idx];
                Scalar sum = x[i];
                for (index_t k = rp[i]; k < rp[i + 1]; ++k) {
                    index_t j = ci[k];
                    if (j < i) {
                        sum -= v[k] * y[j];
                    }
                }
                y[i] = sum * inv_diag_[i];
            });
        }
    }

    /**
     * @brief Backward Gauss-Seidel sweep: solve (D + U) * y = x
     *
     * Uses level scheduling for parallelism.
     */
    void backward_sweep(const Scalar* x, Scalar* y) const {
        const index_t* rp = row_ptr_;
        const index_t* ci = col_idx_;
        const Scalar* v = values_;

        for (const auto& level : bwd_schedule_.levels) {
            const index_t level_size = static_cast<index_t>(level.size());
            parallel_for(level_size, [&](index_t idx) {
                const index_t i = level[idx];
                Scalar sum = x[i];
                for (index_t k = rp[i]; k < rp[i + 1]; ++k) {
                    index_t j = ci[k];
                    if (j > i) {
                        sum -= v[k] * y[j];
                    }
                }
                y[i] = sum * inv_diag_[i];
            });
        }
    }
};

} // namespace sparsesolv

#endif // SPARSESOLV_PRECONDITIONERS_SGS_PRECONDITIONER_HPP
