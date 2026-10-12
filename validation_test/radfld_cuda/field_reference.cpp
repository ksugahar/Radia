// Validation-only pybind11 adapter over the independently compiled pre-CUDA
// RadFld implementation. This module is never a GPU runtime dependency.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "rad_hdiv_field_evaluator.h"
namespace py=pybind11;
using Array=py::array_t<double,py::array::c_style|py::array::forcecast>;
using Field=rad_hdiv::HDivFieldEvaluator;
PYBIND11_MODULE(_radia_field_reference,m) {
    py::class_<Field,std::shared_ptr<Field>>(m,"Field")
        .def(py::init([](Array xyz,Array q,Array masks,Array signs,Array angles,int leaf,double theta) {
            if(xyz.ndim()!=2||xyz.shape(1)!=3||q.ndim()!=1||q.size()!=xyz.shape(0))
                throw std::invalid_argument("Invalid point cloud dimensions");
            std::vector<int> im;
            for(py::ssize_t i=0;i<masks.size();++i) im.push_back(int(masks.data()[i]));
            rad_hdiv::FieldEvaluatorOptions options;options.leaf_size=leaf;options.theta=theta;
            auto f=Field::FromCloud(std::vector<double>(xyz.data(),xyz.data()+xyz.size()),
                std::vector<double>(q.data(),q.data()+q.size()),im,
                std::vector<double>(signs.data(),signs.data()+signs.size()),options);
            f->SetImageRotations(std::vector<double>(angles.data(),angles.data()+angles.size()));
            return f;
        }))
        .def("field",[](Field& f,Array obs,const std::string& algorithm) {
            if(obs.ndim()!=2||obs.shape(1)!=3) throw std::invalid_argument("Expected (n,3) observations");
            py::array_t<double> out({obs.shape(0),py::ssize_t(3)});
            { py::gil_scoped_release release;
              f.Evaluate(obs.data(),std::size_t(obs.shape(0)),out.mutable_data(),Field::ParseAlgorithm(algorithm)); }
            return out;
        });
}
