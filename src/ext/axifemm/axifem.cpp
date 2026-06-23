// axifem.cpp — pybind11 module entry.
//
// The Radia build defines RADIA_AXIFEMM_PHASE_2B and exports the custom
// FiniteElement / FESpace / BFI implementation.  The version() and hello()
// helpers are kept as low-cost import smoke tests.

#include <comp.hpp>
#include <python_comp.hpp>

#ifdef RADIA_AXIFEMM_PHASE_2B
namespace axifem {
    void ExportAxiHenrotteFE(pybind11::module& m);
    void ExportAxiHenrotteFESpace(pybind11::module& m);
    void ExportAxiHenrotteIntegrators(pybind11::module& m);
}
#endif

PYBIND11_MODULE(axifem, m) {
    m.doc() = "axifem — FEMM/Henrotte axisymmetric finite elements for NGSolve.";
    m.attr("__version__") = "0.1.0";

    m.def("version", []() { return "0.1.0"; },
          "Return the axifem version string.");

    m.def("hello", []() {
        return std::string("axifem: Henrotte axisymmetric FE for NGSolve");
    }, "Smoke test: verifies the C++ module loaded.");

#ifdef RADIA_AXIFEMM_PHASE_2B
    axifem::ExportAxiHenrotteFE(m);
    axifem::ExportAxiHenrotteFESpace(m);
    axifem::ExportAxiHenrotteIntegrators(m);
#endif
}
