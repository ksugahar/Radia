// radia_vim.cpp — pybind11 module entry.
//
// Phase F-1 bootstrap: minimum viable build, exports
//   - version()          : version string
//   - hello()            : smoke test
//   - precision_info()   : reports which scalar types (double / quad / mpfr)
//                          are compiled in (controlled by RADIA_VIM_WITH_QUAD
//                          and RADIA_VIM_WITH_MPFR CMake options)
//   - newton_kernel(rA, rB) : direct evaluation of 1/|rA - rB|, demonstrating
//                          the precision plumbing
//
// The actual finite-element / quadrature / assembler implementations are
// gated behind RADIA_VIM_PHASE_F2..F4 CMake flags (off by default for the
// bootstrap build).

#include <comp.hpp>
#include <python_comp.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "precision_traits.hpp"

#include <array>
#include <cmath>

namespace py = pybind11;

#ifdef RADIA_VIM_PHASE_F2
namespace radia_vim {
    void ExportHDivDivFreeHexBasis(py::module& m);
}
#endif
#ifdef RADIA_VIM_PHASE_F3
namespace radia_vim {
    void ExportSauterSchwabHex(py::module& m);
}
#endif
#ifdef RADIA_VIM_PHASE_F4
namespace radia_vim {
    void ExportVolumeGalerkin(py::module& m);
}
#endif

namespace {

// Direct Newton-kernel evaluation, templated on scalar type.
// Used to verify that the precision plumbing compiles for each enabled type.
template <typename T>
T newton_kernel_impl(const std::array<T, 3>& rA,
                     const std::array<T, 3>& rB) {
    T dx = rA[0] - rB[0];
    T dy = rA[1] - rB[1];
    T dz = rA[2] - rB[2];
    T r2 = dx*dx + dy*dy + dz*dz;
    using std::sqrt;          // ADL picks up Boost.Multiprecision overloads
    return T(1) / sqrt(r2);
}

}  // anonymous namespace

PYBIND11_MODULE(radia_vim, m) {
    m.doc() = "radia_vim — Volume Integral Method (Newton-kernel volume "
              "Galerkin) for cuboid eddy currents, NGSolve extension.";
    m.attr("__version__") = "0.1.0";

    m.def("version", []() { return "0.1.0"; },
          "Return the radia_vim version string.");

    m.def("hello", []() {
        return std::string(
            "radia_vim: Newton-kernel volume Galerkin BEM for hex eddy currents");
    }, "Smoke test: verifies the C++ module loaded.");

    // --- Precision plumbing -------------------------------------------------
    py::class_<radia_vim::PrecisionInfo>(m, "PrecisionInfo")
        .def_readonly("with_double", &radia_vim::PrecisionInfo::with_double)
        .def_readonly("with_quad",   &radia_vim::PrecisionInfo::with_quad)
        .def_readonly("with_mpfr",   &radia_vim::PrecisionInfo::with_mpfr)
        .def_readonly("mpfr_digits", &radia_vim::PrecisionInfo::mpfr_digits)
        .def_readonly("boost_version", &radia_vim::PrecisionInfo::boost_version_str)
        .def("__repr__", [](const radia_vim::PrecisionInfo& p) {
            std::string s = "<PrecisionInfo ";
            s += std::string("double=") + (p.with_double ? "yes" : "no");
            s += std::string(" quad=")   + (p.with_quad   ? "yes" : "no");
            s += std::string(" mpfr=")   + (p.with_mpfr   ? "yes" : "no");
            if (p.with_mpfr) s += " (digits=" + std::to_string(p.mpfr_digits) + ")";
            if (!p.boost_version_str.empty())
                s += " boost=" + p.boost_version_str;
            s += ">";
            return s;
        });

    m.def("precision_info", &radia_vim::current_precision_info,
          "Report which scalar types (double / quad / mpfr) are compiled in.");

    // --- Newton kernel direct eval (precision smoke test) -------------------
    m.def("newton_kernel_double",
          [](std::array<double, 3> rA, std::array<double, 3> rB) {
              return newton_kernel_impl<double>(rA, rB);
          },
          py::arg("rA"), py::arg("rB"),
          "Evaluate 1/|rA - rB| at IEEE double precision. Verifies plumbing.");

#ifdef RADIA_VIM_WITH_QUAD
    m.def("newton_kernel_quad",
          [](std::array<double, 3> rA, std::array<double, 3> rB) {
              std::array<radia_vim::Quad, 3> rAq, rBq;
              for (int i = 0; i < 3; ++i) {
                  rAq[i] = radia_vim::Quad(rA[i]);
                  rBq[i] = radia_vim::Quad(rB[i]);
              }
              radia_vim::Quad result = newton_kernel_impl<radia_vim::Quad>(rAq, rBq);
              return result.str();   // return as decimal string (~32 digits)
          },
          py::arg("rA"), py::arg("rB"),
          "Evaluate 1/|rA - rB| at __float128 / cpp_bin_float_quad precision. "
          "Returns decimal string (~32 digits).");
#endif

#ifdef RADIA_VIM_WITH_MPFR
    m.def("newton_kernel_mpfr",
          [](std::array<double, 3> rA, std::array<double, 3> rB, unsigned digits) {
              radia_vim::Mpfr::default_precision(digits);
              std::array<radia_vim::Mpfr, 3> rAm, rBm;
              for (int i = 0; i < 3; ++i) {
                  rAm[i] = radia_vim::Mpfr(rA[i]);
                  rBm[i] = radia_vim::Mpfr(rB[i]);
              }
              radia_vim::Mpfr result = newton_kernel_impl<radia_vim::Mpfr>(rAm, rBm);
              return result.str(digits);
          },
          py::arg("rA"), py::arg("rB"), py::arg("digits") = 60u,
          "Evaluate 1/|rA - rB| at MPFR arbitrary precision. "
          "Returns decimal string with `digits` significant digits.");
#endif

#ifdef RADIA_VIM_PHASE_F2
    radia_vim::ExportHDivDivFreeHexBasis(m);
#endif
#ifdef RADIA_VIM_PHASE_F3
    radia_vim::ExportSauterSchwabHex(m);
#endif
#ifdef RADIA_VIM_PHASE_F4
    radia_vim::ExportVolumeGalerkin(m);
#endif
}
