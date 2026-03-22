/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/// @file preconditioner.hpp
/// @brief Base class for preconditioners

#ifndef SPARSESOLV_CORE_PRECONDITIONER_HPP
#define SPARSESOLV_CORE_PRECONDITIONER_HPP

#include "types.hpp"
#include "constants.hpp"
#include "parallel.hpp"
#include "sparse_matrix_view.hpp"
#include <string>
#include <memory>

namespace sparsesolv {

/// Abstract base class for preconditioners: setup(A), then apply(x, y) = M^{-1}*x
template<typename Scalar = double>
class Preconditioner {
public:
    virtual ~Preconditioner() = default;

    virtual void setup(const SparseMatrixView<Scalar>& A) = 0;
    virtual void apply(const Scalar* x, Scalar* y, index_t size) const = 0;

    void apply(const std::vector<Scalar>& x, std::vector<Scalar>& y) const {
        assert(x.size() == y.size());
        apply(x.data(), y.data(), static_cast<index_t>(x.size()));
    }

    /**
     * @brief Apply preconditioner and compute dot product in one pass
     *
     * Computes y = M^{-1}*x and returns dot(r_for_dot, y).
     * Fusing avoids a separate pass over y[] for the dot product.
     *
     * Default implementation calls apply() then computes dot separately.
     * Derived classes (e.g., ICPreconditioner with ABMC) override this
     * to fuse the dot into the output phase.
     *
     * @param r_for_dot Vector to dot with the output (typically CG residual)
     * @param x Input vector
     * @param y Output vector (M^{-1}*x)
     * @param size Vector size
     * @param conjugate If true, compute conj(r_for_dot) . y (Hermitian)
     * @return dot(r_for_dot, y)
     */
    virtual Scalar apply_fused_dot(
        const Scalar* r_for_dot,
        const Scalar* x, Scalar* y,
        index_t size, bool conjugate) const
    {
        apply(x, y, size);
        return parallel_reduce_sum<Scalar>(size, [&](index_t i) -> Scalar {
            if constexpr (!std::is_same_v<Scalar, double>) {
                return conjugate ? std::conj(r_for_dot[i]) * y[i]
                                 : r_for_dot[i] * y[i];
            } else {
                return r_for_dot[i] * y[i];
            }
        });
    }

    virtual std::string name() const = 0;
    virtual bool is_ready() const { return is_setup_; }

protected:
    bool is_setup_ = false;
};

} // namespace sparsesolv

#endif // SPARSESOLV_CORE_PRECONDITIONER_HPP
