// Validation-only pybind11 adapter over the independently compiled pre-CUDA
// RadFld implementation. This module is never a GPU runtime dependency.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <cmath>
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
        .def_static("from_elements",[](py::dict spec,Array masks,Array signs,Array angles,int leaf,double theta) {
            if(py::cast<double>(spec["version"])!=1.0) throw std::invalid_argument("Unsupported source version");
            auto values=[&](const char* key,int columns,bool optional=false) {
                if(optional&&!spec.contains(key)) return std::vector<double>{};
                Array a=py::cast<Array>(spec[key]);
                if(a.ndim()!=2||a.shape(1)!=columns) throw std::invalid_argument(std::string("Invalid shape: ")+key);
                std::vector<double> result(a.size());
                for(py::ssize_t i=0;i<a.size();++i) {
                    if(!std::isfinite(a.data()[i])) throw std::invalid_argument("Nonfinite source");
                    result[i]=a.data()[i];
                }
                return result;
            };
            std::vector<int> im;
            for(py::ssize_t i=0;i<masks.size();++i) {
                const double mask=masks.data()[i];
                if(!std::isfinite(mask)||mask!=std::floor(mask)||mask<0||mask>7) throw std::invalid_argument("Invalid image mask");
                im.push_back(int(mask));
            }
            std::vector<double> si(signs.size()),an(angles.size());
            for(py::ssize_t i=0;i<signs.size();++i) si[i]=signs.data()[i];
            for(py::ssize_t i=0;i<angles.size();++i) an[i]=angles.data()[i];
            rad_hdiv::FieldEvaluatorOptions options;options.leaf_size=leaf;options.theta=theta;
            const auto kind=py::cast<std::string>(spec["representation"]);
            int vc=0,sc=0;
            if(kind=="tet") {vc=16;sc=22;}
            else if(kind=="polynomial_tet") {vc=32;sc=22;}
            else if(kind=="cubic_tet") {vc=32;sc=29;}
            else if(kind=="sextic_quartic_tet") {vc=96;sc=44;}
            else if(kind=="curved_tet") {vc=34;sc=24;}
            else throw std::invalid_argument("Unsupported element representation");
            auto v=values("volume",vc),t=values("surface",sc);
            auto vo=values("volume_origins",3,true),so=values("surface_origins",3,true);
            std::shared_ptr<Field> f;
            if(kind=="curved_tet") {
                if(!vo.empty()||!so.empty()) throw std::invalid_argument("Curved origins unsupported");
                f=Field::FromCurvedTet(std::move(v),std::move(t),values("gauss_points",1),values("gauss_weights",1),im,si,options);
            } else if(kind=="tet") f=Field::FromTet(std::move(v),std::move(t),im,si,options,std::move(vo),std::move(so));
            else if(kind=="polynomial_tet") f=Field::FromPolynomialTet(std::move(v),std::move(t),im,si,options,std::move(vo),std::move(so));
            else if(kind=="cubic_tet") f=Field::FromCubicPolynomialTet(std::move(v),std::move(t),im,si,options,std::move(vo),std::move(so));
            else f=Field::FromSexticQuarticPolynomialTet(std::move(v),std::move(t),im,si,options,std::move(vo),std::move(so));
            f->SetImageRotations(an); return f;
        })
        .def("field",[](Field& f,Array obs,const std::string& algorithm) {
            if(obs.ndim()!=2||obs.shape(1)!=3) throw std::invalid_argument("Expected (n,3) observations");
            py::array_t<double> out({obs.shape(0),py::ssize_t(3)});
            { py::gil_scoped_release release;
              f.Evaluate(obs.data(),std::size_t(obs.shape(0)),out.mutable_data(),Field::ParseAlgorithm(algorithm)); }
            return out;
        });
}
