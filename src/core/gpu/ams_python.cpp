/* SPDX-License-Identifier: MPL-2.0 */
#include "snapshot.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <mutex>
#include <memory>
namespace py = pybind11;
namespace radia_gpu_ams_python {
using radia_sparse::Node;

Node snapshot(py::handle value, int depth=0) {
    if (depth > 64) throw py::value_error("snapshot nesting exceeds capacity");
    Node out;
    if (py::isinstance<py::dict>(value)) {
        out.kind=Node::record;
        for(auto item:py::reinterpret_borrow<py::dict>(value))
            out.fields[py::cast<std::string>(item.first)]=snapshot(item.second,depth+1);
    } else if (py::isinstance<py::list>(value) || py::isinstance<py::tuple>(value)) {
        out.kind=Node::list;
        for(auto item:py::reinterpret_borrow<py::sequence>(value))
            out.children.push_back(snapshot(item,depth+1));
    } else if (py::isinstance<py::bool_>(value)) {
        out.kind=Node::logical;out.boolean=py::cast<bool>(value);
    } else if (py::isinstance<py::str>(value)) {
        out.kind=Node::text;out.word=py::cast<std::string>(value);
    } else if (py::isinstance<py::int_>(value) || py::isinstance<py::float_>(value)) {
        out.kind=Node::number;out.scalar=py::cast<double>(value);
        out.real();
    } else if (py::isinstance<py::array>(value)) {
        auto a=py::reinterpret_borrow<py::array>(value);
        if(a.ndim()!=1 || !(a.flags() & py::array::c_style) || a.size()>INT_MAX-256)
            throw py::value_error("contiguous one-dimensional snapshot array required");
        out.count=size_t(a.size());
        if(a.dtype().is(py::dtype::of<int>())) {
            out.kind=Node::integers;auto p=static_cast<const int*>(a.data());
            if(out.count)out.int_values.assign(p,p+out.count);
        } else if(a.dtype().is(py::dtype::of<double>())) {
            out.kind=Node::reals;auto p=static_cast<const double*>(a.data());
            if(out.count)out.real_values.assign(p,p+out.count);
        } else throw py::value_error("snapshot arrays require exact int32 or float64 dtype");
    } else throw py::value_error("unsupported accelerator snapshot value");
    return out;
}

struct Resident {
    void* owner=nullptr;
    int n;
    std::string kind;
    std::mutex mutex;
    explicit Resident(py::dict input) {
        Node data=snapshot(input);
        n=data.at("A").at("rows").integer();
        kind=data.at("schema").string();
        { py::gil_scoped_release release;owner=radia_sparse::create(data); }
    }
    Resident(const Resident&)=delete;
    ~Resident(){radia_sparse::release(owner);}
    void close(){
        py::gil_scoped_release release;std::lock_guard<std::mutex> lock(mutex);
        radia_sparse::release(owner);owner=nullptr;
    }
    py::tuple run(py::array_t<double,py::array::c_style> rhs, bool solve,
                  double tolerance=1e-12,int iterations=10000){
        if(rhs.ndim()!=1||rhs.size()!=n||!std::isfinite(tolerance)||tolerance<=0
                ||tolerance>1e-6||iterations<1||iterations>100000)
            throw py::value_error("invalid sparse GPU operation inputs");
        std::vector<double> b(rhs.data(),rhs.data()+rhs.size());
        for(double value:b)if(!std::isfinite(value))throw py::value_error("nonfinite RHS");
        py::array_t<double> output(n);auto p=output.mutable_data();radia_sparse::Stats stats;
        {
            // Acquire/release the owner lock while the GIL is released so close
            // and concurrent calls cannot invalidate buffers or deadlock.
            py::gil_scoped_release release;std::lock_guard<std::mutex> lock(mutex);
            if(!owner)throw std::runtime_error("GPU hierarchy is closed");
            if(solve)stats=radia_sparse::solve(owner,b.data(),b.size(),p,tolerance,iterations);
            else radia_sparse::apply(owner,b.data(),b.size(),p);
        }
        py::dict info;info["backend"]="gpu";info["precision"]="float64";
        info["preconditioner"]=kind;
        info["operation"]=solve?"solve":"preconditioner_apply";
        info["residual_checked"]=solve;info["transfer_included"]=true;
        if(solve){info["iterations"]=stats.iterations;info["true_relative_residual"]=stats.residual;}
        return py::make_tuple(output,info);
    }
};
}  // namespace radia_gpu_ams_python
PYBIND11_MODULE(_gpu_ams,m){
    using Resident = radia_gpu_ams_python::Resident;
    m.doc()="Radia-owned FP64 AMG/AMS cycles; no MATLAB or NGSolve ABI";
    py::class_<Resident>(m,"Resident")
        .def(py::init<py::dict>())
        .def("close",&Resident::close)
        .def("solve",[](Resident& r,py::array_t<double,py::array::c_style> b,double t,int k){return r.run(b,true,t,k);},
             py::arg("rhs").noconvert(),py::arg("tolerance"),py::arg("max_iterations"))
        .def("apply",[](Resident& r,py::array_t<double,py::array::c_style> b){return r.run(b,false);},
             py::arg("rhs").noconvert());
}
