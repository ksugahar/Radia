// sauter_schwab_hex.cpp — pybind11 bindings for spherical-Duffy quadrature.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "sauter_schwab_hex.hpp"

namespace py = pybind11;

namespace radia_vim {

void ExportSauterSchwabHex(py::module& m) {
    m.def("compute_K_entry_spherical_duffy",
          [](const HDivDivFreeHexBasis<double>& basis,
             int i, int j,
             double a, double b, double c,
             int n_omega_theta, int n_omega_phi,
             int n_rho, int n_c) {
              return compute_K_entry_spherical_duffy<double>(
                  basis, i, j, a, b, c,
                  n_omega_theta, n_omega_phi, n_rho, n_c);
          },
          py::arg("basis"), py::arg("i"), py::arg("j"),
          py::arg("a"), py::arg("b"), py::arg("c"),
          py::arg("n_omega_theta") = 16,
          py::arg("n_omega_phi") = 32,
          py::arg("n_rho") = 16,
          py::arg("n_c") = 8,
          "Compute K[i,j] entry via spherical-Duffy quadrature.\n"
          "Returns the bare integral; downstream code should multiply by\n"
          "mu0/(4*pi) for the physical K matrix.\n\n"
          "Args:\n"
          "  basis: HDivDivFreeHexBasis instance\n"
          "  i, j : DOF indices\n"
          "  a, b, c: cuboid dimensions [m]\n"
          "  n_omega_theta, n_omega_phi: GL order on (θ, φ) sphere\n"
          "  n_rho: GL order on radial ρ\n"
          "  n_c:   GL order on each c_x, c_y, c_z (3 dirs)");

    m.def("gauss_legendre_m1_1",
          [](int n) {
              std::vector<double> nodes, weights;
              gauss_legendre_m1_1<double>(n, nodes, weights);
              return std::make_pair(nodes, weights);
          },
          py::arg("n"),
          "Generate Gauss-Legendre nodes and weights on [-1, 1].\n"
          "Returns (nodes, weights) tuple.");
}

}  // namespace radia_vim
