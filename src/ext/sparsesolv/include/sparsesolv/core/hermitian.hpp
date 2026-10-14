/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef SPARSESOLV_CORE_HERMITIAN_HPP
#define SPARSESOLV_CORE_HERMITIAN_HPP

#include "sparse_matrix_view.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

namespace sparsesolv {

// Structural validation is not a positive-definiteness certificate. The caller
// must supply an HPD matrix; CG also checks curvature along visited directions.
template<class Scalar>
void validate_hermitian(const SparseMatrixView<Scalar>& A) {
    if (A.rows() != A.cols() || A.rows() < 0 || !A.row_ptr()
        || A.row_ptr()[0] != 0)
        throw std::invalid_argument("Hermitian CG requires square CSR storage");
    const auto* rp = A.row_ptr();
    const auto* ci = A.col_idx();
    const auto* av = A.values();
    for (index_t i = 0; i < A.rows(); ++i) {
        if (rp[i] > rp[i + 1])
            throw std::invalid_argument("Hermitian CG requires ordered CSR rows");
        index_t previous = -1;
        for (index_t k = rp[i]; k < rp[i + 1]; ++k) {
            if (ci[k] <= previous || ci[k] >= A.cols()
                || !std::isfinite(std::real(av[k])) || !std::isfinite(std::imag(av[k])))
                throw std::invalid_argument("Hermitian CG requires finite sorted unique CSR entries");
            previous = ci[k];
        }
    }
    constexpr double tolerance = 64 * std::numeric_limits<double>::epsilon();
    auto component_scale = [](const Scalar& v) {
        return std::max(std::abs(std::real(v)), std::abs(std::imag(v)));
    };
    for (index_t i = 0; i < A.rows(); ++i) {
        const Scalar diagonal = A(i, i);
        if (!(std::real(diagonal) > 0)
            || std::abs(std::imag(diagonal)) > tolerance * std::real(diagonal))
            throw std::invalid_argument("Hermitian CG requires positive real diagonal entries");
        for (index_t k = rp[i]; k < rp[i + 1]; ++k) {
            const Scalar reverse = A(ci[k], i);
            const auto expected = std::conj(reverse);
            // Assembly can cancel almost all of an off-diagonal entry. Scale
            // symmetry error by the participating diagonals as well.
            const double scale = std::max({component_scale(av[k]), component_scale(reverse),
                                           component_scale(diagonal), component_scale(A(ci[k], ci[k]))});
            // Divide before subtracting to avoid overflow at extreme scales.
            if (scale > 0 && std::abs(av[k] / scale - expected / scale) > tolerance)
                throw std::invalid_argument("conjugate=True requires a Hermitian matrix");
        }
    }
}
} // namespace sparsesolv
#endif
