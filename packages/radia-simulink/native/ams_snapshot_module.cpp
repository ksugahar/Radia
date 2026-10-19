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
    m.def("_check_openmp_contract", [](int threads) {
        const int before=ngla::CompactSetupThreads();
        {
            ngla::CompactSetupThreadScope scope(threads);
            std::vector<int> values(4096,-1);
            ngla::CompactParallelFor(values.size(),[&](size_t i) { values[i]=int(i); });
            for(size_t i=0;i<values.size();++i)
                if(values[i]!=int(i))throw std::runtime_error("OpenMP row coverage failed");
            bool propagated=false;
            try {
                ngla::CompactParallelFor(values.size(),[](size_t i) {
                    if(i==2048)throw std::runtime_error("worker exception probe");
                });
            } catch(const std::runtime_error& e) {
                propagated=std::string(e.what())=="worker exception probe";
            }
            if(!propagated)throw std::runtime_error("OpenMP exception propagation failed");
        }
        if(ngla::CompactSetupThreads()!=before)
            throw std::runtime_error("OpenMP setup thread scope leaked");
        return true;
    });
    m.def("snapshot", [](std::shared_ptr<ngla::BaseMatrix> matrix,
        std::shared_ptr<ngla::BaseMatrix> gradient, std::shared_ptr<ngcore::BitArray> free,
        std::vector<double> x,std::vector<double> y,std::vector<double> z,
        py::array_t<double,py::array::c_style|py::array::forcecast> probes, int cycle, int setup_threads) {
        if (ngcore::GetTaskManager() != nullptr)
            throw std::invalid_argument("mexcuda AMS preprocessing requires no active TaskManager");
        auto a=std::dynamic_pointer_cast<ngla::SparseMatrix<double>>(matrix);
        auto g=std::dynamic_pointer_cast<ngla::SparseMatrix<double>>(gradient);
        if(!a || !g) throw std::invalid_argument("real nonsymmetric-storage NGSolve sparse matrices required");
        if(probes.ndim()!=2 || probes.shape(1)!=a->Height())
            throw std::invalid_argument("probe dimensions do not match matrix");
        std::shared_ptr<ngla::HypreBasedAMS> ams;
        py::dict result;
        double hierarchy_seconds, export_seconds;
        {
            ngla::CompactSetupThreadScope scope(setup_threads);
            auto start=std::chrono::steady_clock::now();
            ams=std::make_shared<ngla::HypreBasedAMS>(a,g,free,x,y,z,cycle);
            hierarchy_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            start=std::chrono::steady_clock::now();
            result=ngla::AcceleratorAMS(*ams);
            export_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        }
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
        result["setup_diagnostics"]=py::dict(py::arg("hierarchy_seconds")=hierarchy_seconds,
            py::arg("export_seconds")=export_seconds,py::arg("requested_threads")=setup_threads,
            py::arg("strength_workers")=ams->GetSetupWorkers(),py::arg("taskmanager_used")=false);
        std::array<double,7> phases{};
        for(const auto* auxiliary : {ams->GetBGAsAMG(),ams->GetBPixAsAMG(),
                                    ams->GetBPiyAsAMG(),ams->GetBPizAsAMG()})
            if(auxiliary)for(size_t i=0;i<phases.size();++i)phases[i]+=auxiliary->SetupPhases()[i];
        result["setup_diagnostics"].cast<py::dict>()["amg_phase_seconds"]=py::cast(phases);
        return result;
    },py::arg("matrix"),py::arg("gradient"),py::arg("free"),py::arg("x"),py::arg("y"),
      py::arg("z"),py::arg("probes"),py::arg("cycle")=1,py::arg("setup_threads")=1);
}
