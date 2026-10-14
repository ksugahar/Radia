/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file norms.hpp
 * @brief Euclidean norms without avoidable intermediate underflow or overflow
 */

#ifndef SPARSESOLV_CORE_NORMS_HPP
#define SPARSESOLV_CORE_NORMS_HPP

#include "types.hpp"
#include "parallel.hpp"
#include <cmath>
#include <complex>
#include <limits>
#include <type_traits>

namespace sparsesolv {

/// v * 2^e computed on the exponent (exact unless the result leaves the range)
inline double scale_by_pow2(double v, int e) { return std::scalbn(v, e); }
inline std::complex<double> scale_by_pow2(const std::complex<double>& v, int e) {
    return {std::scalbn(v.real(), e), std::scalbn(v.imag(), e)};
}

/**
 * @brief Largest |Re| or |Im| component of v
 *
 * Returns 0 for a zero vector, NaN if any component is NaN and +inf if any
 * component is infinite (NaN takes precedence), so non-finite input is never
 * hidden by the maximum.
 */
template<typename Scalar>
double max_component(const Scalar* v, index_t n) {
    double m = 0.0;
    bool has_inf = false;
    auto visit = [&](double c) {
        if (std::isnan(c)) return false;
        if (std::isinf(c)) { has_inf = true; return true; }
        const double a = std::abs(c);
        if (a > m) m = a;
        return true;
    };
    for (index_t i = 0; i < n; ++i) {
        if constexpr (std::is_same_v<Scalar, double>) {
            if (!visit(v[i])) return std::numeric_limits<double>::quiet_NaN();
        } else {
            if (!visit(v[i].real()) || !visit(v[i].imag()))
                return std::numeric_limits<double>::quiet_NaN();
        }
    }
    return has_inf ? std::numeric_limits<double>::infinity() : m;
}

/// Exponent k with 2^k <= m < 2^(k+1) for finite m > 0 (subnormals included)
inline int scale_exponent(double m) { return std::ilogb(m); }

/**
 * @brief ||v||_2 from the components scaled by 2^-k on their exponents
 *
 * The scaling is exact, so for vectors of ordinary magnitude the result is
 * the plain sqrt(sum |v_i|^2); tiny (including subnormal) or huge vectors do
 * not underflow or overflow in the intermediate sum.  The returned norm can
 * still be +inf when the norm itself exceeds the double range.  Non-finite
 * components give NaN or +inf, as max_component does.
 */
template<typename Scalar>
double robust_norm(const Scalar* v, index_t n) {
    const double m = max_component(v, n);
    if (!(m > 0.0) || !std::isfinite(m)) return m;  // 0, inf or NaN
    const int k = scale_exponent(m);
    const double sum = parallel_reduce_sum<double>(n, [&](index_t i) {
        return std::norm(scale_by_pow2(v[i], -k));
    });
    return std::scalbn(std::sqrt(sum), k);
}

} // namespace sparsesolv

#endif // SPARSESOLV_CORE_NORMS_HPP
