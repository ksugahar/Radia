// precision_traits.hpp — type-templated precision support for radia_vim.
//
// The Newton-kernel volume Galerkin pipeline assembles a dense K matrix
// whose entries feed an eigenvalue + Cauer extraction pipeline. Empirical
// study (scripts/precision_sensitivity_study.py) shows:
//
//   precision         clean Cauer pairs (synthetic)
//   ----------------- -----------------------------
//   double (16 digits)        6-25 (mid-stage error grows beyond 30 modes)
//   __float128 (32 digits)    25-40 (mid-stage error ~6% at 50 modes)
//   mpfr  60+ digits          50+   (full reference accuracy)
//
// The K assembler is templated on the scalar type so the same code path
// can be compiled for double (fast), quad (medium), or mpfr (slow but
// extends Cauer extraction to many stages).
//
// Header-only Boost.Multiprecision provides cpp_bin_float_quad and
// mpfr_float; no link dependency for the headers themselves (mpfr_float
// requires linking against MPFR/GMP shared libraries).

#pragma once

#include <type_traits>
#include <string>

#ifdef RADIA_VIM_WITH_QUAD
  #include <boost/multiprecision/cpp_bin_float.hpp>
#endif

#ifdef RADIA_VIM_WITH_MPFR
  #include <boost/multiprecision/mpfr.hpp>
#endif

namespace radia_vim {

// Tag struct identifying which scalar types are compiled in.
struct PrecisionInfo {
    bool with_double = true;          // always available
    bool with_quad   = false;
    bool with_mpfr   = false;
    int  mpfr_digits = 0;             // 0 = not compiled in
    std::string boost_version_str;    // empty if neither QUAD nor MPFR
};

inline PrecisionInfo current_precision_info() {
    PrecisionInfo info;
#ifdef RADIA_VIM_WITH_QUAD
    info.with_quad = true;
#endif
#ifdef RADIA_VIM_WITH_MPFR
    info.with_mpfr = true;
    info.mpfr_digits = 60;            // default — overridden at runtime later
#endif
#if defined(RADIA_VIM_WITH_QUAD) || defined(RADIA_VIM_WITH_MPFR)
    info.boost_version_str =
        std::to_string(BOOST_VERSION / 100000) + "." +
        std::to_string((BOOST_VERSION / 100) % 1000) + "." +
        std::to_string(BOOST_VERSION % 100);
#endif
    return info;
}

// Type aliases. Always-present: Double. Conditionally-present: Quad, Mpfr.
using Double = double;

#ifdef RADIA_VIM_WITH_QUAD
  using Quad = boost::multiprecision::cpp_bin_float_quad;
#endif

#ifdef RADIA_VIM_WITH_MPFR
  using Mpfr = boost::multiprecision::mpfr_float;
#endif

// Compile-time predicate: does T have hardware-fast arithmetic?
template <typename T>
constexpr bool is_native_v = std::is_floating_point<T>::value;

}  // namespace radia_vim
