// volume_galerkin.cpp — pybind11 bindings for the K matrix assembler.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "volume_galerkin.hpp"

namespace py = pybind11;

namespace radia_vim {

void ExportVolumeGalerkin(py::module& m) {
    m.def("assemble_K_bare",
          [](const HDivDivFreeHexBasis<double>& basis,
             double a, double b, double c,
             int n_omega_theta, int n_omega_phi,
             int n_rho, int n_c,
             int verbose) {
              const int n = basis.n_dofs();
              std::vector<double> Kflat;
              {
                  // Heavy lifting outside the GIL so OpenMP threads can
                  // parallel-compute without contention.
                  py::gil_scoped_release nogil;
                  Kflat = assemble_K_bare<double>(
                      basis, a, b, c,
                      n_omega_theta, n_omega_phi, n_rho, n_c, verbose);
              }
              // Now inside the GIL: create the NumPy array.
              py::array_t<double> arr({n, n});
              auto buf = arr.mutable_unchecked<2>();
              for (int i = 0; i < n; ++i)
                  for (int j = 0; j < n; ++j)
                      buf(i, j) = Kflat[static_cast<std::size_t>(i) * n + j];
              return arr;
          },
          py::arg("basis"),
          py::arg("a"), py::arg("b"), py::arg("c"),
          py::arg("n_omega_theta") = 16,
          py::arg("n_omega_phi")   = 32,
          py::arg("n_rho")         = 16,
          py::arg("n_c")           = 8,
          py::arg("verbose")       = 0,
          "Assemble the BARE K matrix (no mu0/(4 pi) prefactor) for the\n"
          "HDiv div-free hex BEM on cuboid (a x b x c).\n\n"
          "Returns a NumPy (n_dofs, n_dofs) double array. K is symmetric;\n"
          "we compute the upper triangle and mirror.\n\n"
          "Args:\n"
          "  basis            : HDivDivFreeHexBasis instance\n"
          "  a, b, c          : cuboid dimensions [m]\n"
          "  n_omega_theta etc: spherical-Duffy quadrature orders (see\n"
          "                     compute_K_entry_spherical_duffy)\n"
          "  verbose          : 0 = silent, k > 0 = stderr report every k pairs\n"
          "                     (use a large k like 100 for moderate n_dofs)\n\n"
          "OpenMP-parallel over (i, j) pairs if compiled with OpenMP.\n"
          "Multiply by mu0/(4 pi) downstream to recover the physical K.");
}

}  // namespace radia_vim
