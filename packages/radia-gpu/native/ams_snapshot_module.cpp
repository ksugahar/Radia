/* SPDX-License-Identifier: MPL-2.0 */
// Validation-only adapter; production AMS setup is compiled from its shared header.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "sparsesolv/preconditioners/compact_ams.hpp"
#include "sparsesolv/ngsolve/ams_accelerator_snapshot.hpp"
namespace py = pybind11;
PYBIND11_MODULE(_radia_ams_snapshot, m) {
    py::module_::import("ngsolve.la"); py::module_::import("ngsolve.comp");
    m.def("snapshot", [](std::shared_ptr<ngla::BaseMatrix> matrix,
        std::shared_ptr<ngla::BaseMatrix> gradient, std::shared_ptr<ngcore::BitArray> free,
        std::vector<double> x,std::vector<double> y,std::vector<double> z,
        py::array_t<double,py::array::c_style|py::array::forcecast> probes, int cycle) {
        auto a=std::dynamic_pointer_cast<ngla::SparseMatrix<double>>(matrix);
        auto g=std::dynamic_pointer_cast<ngla::SparseMatrix<double>>(gradient);
        if(!a || !g) throw std::invalid_argument("real nonsymmetric-storage NGSolve sparse matrices required");
        if(probes.ndim()!=2 || probes.shape(1)!=a->Height())
            throw std::invalid_argument("probe dimensions do not match matrix");
        auto ams=std::make_shared<ngla::HypreBasedAMS>(a,g,free,x,y,z,cycle);
        auto result=ngla::AcceleratorAMS(*ams);
        py::array_t<double> expected({probes.shape(0),probes.shape(1)});
        ngla::VVector<double> b(a->Height()), v(a->Height());
        auto p=probes.unchecked<2>(); auto q=expected.mutable_unchecked<2>();
        // CPU oracle uses the actual shared AMS Mult, not a rewritten algorithm.
        for(py::ssize_t k=0;k<probes.shape(0);++k) {
            for(int i=0;i<a->Height();++i) b.FVDouble()[i]=p(k,i);
            ams->Mult(b,v);
            for(int i=0;i<a->Height();++i) q(k,i)=v.FVDouble()[i];
        }
        result["cpu_probe_results"]=expected;
        result["cpu_preconditioner"]=std::static_pointer_cast<ngla::BaseMatrix>(ams);
        return result;
    },py::arg("matrix"),py::arg("gradient"),py::arg("free"),py::arg("x"),py::arg("y"),
      py::arg("z"),py::arg("probes"),py::arg("cycle")=1);
}
