/* SPDX-License-Identifier: MPL-2.0 */
// Validation-only Python CPU route: normal NGSolve TaskManager dispatch,
// with no GPU snapshot/export work and no mexcuda/OpenMP dispatch macros.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "sparsesolv/preconditioners/compact_ams.hpp"
namespace py=pybind11;
PYBIND11_MODULE(_radia_ams_cpu_reference,m){
    py::module_::import("ngsolve.la");py::module_::import("ngsolve.comp");
    m.def("preconditioner",[](std::shared_ptr<ngla::BaseMatrix> matrix,
        std::shared_ptr<ngla::BaseMatrix> gradient,std::shared_ptr<ngcore::BitArray> free,
        std::vector<double> x,std::vector<double> y,std::vector<double> z,int cycle){
        if(ngcore::GetTaskManager()!=nullptr)throw std::invalid_argument("CPU AMS setup requires no active caller TaskManager");
        auto a=std::dynamic_pointer_cast<ngla::SparseMatrix<double>>(matrix);
        auto g=std::dynamic_pointer_cast<ngla::SparseMatrix<double>>(gradient);
        if(!a||!g)throw std::invalid_argument("real NGSolve sparse matrices required");
        return std::static_pointer_cast<ngla::BaseMatrix>(std::make_shared<ngla::HypreBasedAMS>(a,g,free,x,y,z,cycle));
    },py::arg("matrix"),py::arg("gradient"),py::arg("free"),py::arg("x"),py::arg("y"),py::arg("z"),py::arg("cycle")=1);
}
