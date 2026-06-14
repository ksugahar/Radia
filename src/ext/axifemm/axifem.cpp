// radia_axifemm.cpp — pybind11 module entry.
//
// Phase 2-A bootstrap: minimum viable build, exports version() and hello()
// only. The custom FiniteElement / FESpace implementations are gated behind
// RADIA_AXIFEMM_PHASE_2B (off by default for the bootstrap build).

#include <comp.hpp>
#include <python_comp.hpp>

#ifdef RADIA_AXIFEMM_PHASE_2B
namespace radia_axifemm {
    void ExportAxiHenrotteFE(pybind11::module& m);
    void ExportAxiHenrotteFESpace(pybind11::module& m);
    void ExportAxiHenrotteIntegrators(pybind11::module& m);
}
#endif

PYBIND11_MODULE(radia_axifemm, m) {
    m.doc() = "radia_axifemm — FEMM/Henrotte axisymmetric finite elements for NGSolve.";
    m.attr("__version__") = "0.1.0";

    m.def("version", []() { return "0.1.0"; },
          "Return the radia_axifemm version string.");

    m.def("hello", []() {
        return std::string("radia_axifemm: Henrotte axisymmetric FE for NGSolve");
    }, "Smoke test: verifies the C++ module loaded.");

#ifdef RADIA_AXIFEMM_PHASE_2B
    radia_axifemm::ExportAxiHenrotteFE(m);
    radia_axifemm::ExportAxiHenrotteFESpace(m);
    radia_axifemm::ExportAxiHenrotteIntegrators(m);
#endif
}
