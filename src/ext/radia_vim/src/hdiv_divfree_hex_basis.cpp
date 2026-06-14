// hdiv_divfree_hex_basis.cpp — pybind11 bindings for HDiv div-free hex basis.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "hdiv_divfree_hex_basis.hpp"

namespace py = pybind11;

namespace radia_vim {

void ExportHDivDivFreeHexBasis(py::module& m) {
    py::enum_<HexBasisBlock>(m, "HexBasisBlock")
        .value("A", HexBasisBlock::A)
        .value("B", HexBasisBlock::B)
        .value("C", HexBasisBlock::C)
        .value("D", HexBasisBlock::D)
        .value("E", HexBasisBlock::E);

    py::class_<HexBasisDofKey>(m, "HexBasisDofKey")
        .def_readonly("block", &HexBasisDofKey::block)
        .def_readonly("i",     &HexBasisDofKey::i)
        .def_readonly("j",     &HexBasisDofKey::j)
        .def_readonly("k",     &HexBasisDofKey::k)
        .def("__repr__", [](const HexBasisDofKey& k) {
            std::string s = "<HexBasisDofKey block=";
            switch (k.block) {
              case HexBasisBlock::A: s += "A"; break;
              case HexBasisBlock::B: s += "B"; break;
              case HexBasisBlock::C: s += "C"; break;
              case HexBasisBlock::D: s += "D"; break;
              case HexBasisBlock::E: s += "E"; break;
            }
            s += " i=" + std::to_string(k.i);
            s += " j=" + std::to_string(k.j);
            s += " k=" + std::to_string(k.k);
            s += ">";
            return s;
        });

    py::class_<HDivDivFreeHexBasis<double>>(m, "HDivDivFreeHexBasis")
        .def(py::init<int>(), py::arg("order"),
             "HDiv div-free interior basis on the unit hex [0,1]^3. "
             "Total DOFs = 2*order^3 + 3*order^2.")
        .def_property_readonly("order",  &HDivDivFreeHexBasis<double>::order)
        .def_property_readonly("n_dofs", &HDivDivFreeHexBasis<double>::n_dofs)
        .def("dof_key", &HDivDivFreeHexBasis<double>::dof_key, py::arg("idx"),
             "Return the (block, i, j, k) index tuple for DOF idx.")
        .def("evaluate",
             [](const HDivDivFreeHexBasis<double>& self,
                int idx, double x, double y, double z) {
                 auto v = self.evaluate(idx, x, y, z);
                 return std::array<double, 3>{v.x, v.y, v.z};
             },
             py::arg("idx"), py::arg("x"), py::arg("y"), py::arg("z"),
             "Evaluate basis function `idx` at reference point (x, y, z) in [0,1]^3.\n"
             "Returns [vx, vy, vz] in IEEE double precision.");
}

}  // namespace radia_vim
