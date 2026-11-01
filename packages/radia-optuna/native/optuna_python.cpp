#include "optuna_cuda.hpp"
#include "optuna_runtime.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <cmath>
#include <stdexcept>
#include <limits>
namespace py = pybind11;
namespace {
struct Input {
    std::vector<double> weights, mu, sigma, probabilities;
    std::size_t choices;
    double low, high, step;
    bool categorical, log_scale;
    explicit Input(py::dict d) {
        weights=d["weights"].cast<std::vector<double>>();
        categorical=d["categorical"].cast<bool>();
        choices=categorical ? d["choices"].cast<std::size_t>() : 0;
        low=d["low"].cast<double>(); high=d["high"].cast<double>();
        step=d["step"].cast<double>(); log_scale=d["log_scale"].cast<bool>();
        if (weights.empty() || !std::isfinite(low) || !std::isfinite(high) || !(low < high) ||
            (!std::isnan(step) && (!std::isfinite(step) || step<=0)))
            throw std::invalid_argument("Estimator bounds/weights/step");
        if (categorical) {
            probabilities=d["probabilities"].cast<std::vector<double>>();
            if (!choices || weights.size()>std::numeric_limits<std::size_t>::max()/choices ||
                probabilities.size()!=weights.size()*choices)
                throw std::invalid_argument("Estimator categorical shape");
            for(double x:probabilities) if(!std::isfinite(x) || x<0)
                throw std::invalid_argument("Estimator probabilities");
        } else {
            mu=d["mu"].cast<std::vector<double>>();sigma=d["sigma"].cast<std::vector<double>>();
            if (mu.size()!=weights.size() || sigma.size()!=weights.size())
                throw std::invalid_argument("Estimator numerical shape");
            for (double x:sigma) if (!(x>0) || !std::isfinite(x))
                throw std::invalid_argument("Estimator sigma");
            for (double x:mu) if(!std::isfinite(x)) throw std::invalid_argument("Estimator mu");
        }
    }
    optuna_cuda::Estimator view() const {
        return {&weights,&mu,&sigma,&probabilities,choices,low,high,step,categorical,log_scale};
    }
};
struct Inputs {
    std::vector<Input> owners;
    std::vector<optuna_cuda::Estimator> views;
    explicit Inputs(py::list specs) {
        for (auto item:specs) owners.emplace_back(item.cast<py::dict>());
        if (owners.empty()) throw std::invalid_argument("Empty estimator dimensions");
        for (const auto& owner:owners) views.push_back(owner.view());
        for (const auto& owner:owners) {
            if (owner.weights!=owners.front().weights)
                throw std::invalid_argument("Shared mixture weights required");
            double sum=0;
            for(double x:owner.weights) {
                if(!std::isfinite(x) || x<0) throw std::invalid_argument("Mixture weights");
                sum+=x;
            }
            if(!(sum>0) || !std::isfinite(sum)) throw std::invalid_argument("Empty mixture weights");
        }
    }
    void validate(const std::vector<double>& samples) const {
        if(samples.size()%owners.size()) throw std::invalid_argument("Row-major samples shape");
        for(std::size_t i=0;i<samples.size();++i) {
            const auto& d=owners[i%owners.size()];const auto x=samples[i];
            if(!std::isfinite(x) || (d.log_scale && x<=0) ||
                (d.categorical && (x<1 || x>d.choices || x!=std::floor(x))))
                throw std::invalid_argument("Sample value / one-based categorical index");
        }
    }
};
struct Model {
    optuna_cuda::ModelConstruction result;
    std::vector<optuna_cuda::ModelDimension> dimensions;
    std::vector<double> weights, empty;
    std::vector<optuna_cuda::Estimator> views() const {
        std::vector<optuna_cuda::Estimator> out;
        for (const auto& d:dimensions)
            out.push_back({&weights,&empty,&empty,&empty,d.choices,d.low,d.high,d.step,
                           d.categorical,d.log_scale,result.device.get()});
        return out;
    }
};
void backend(const std::string& mode) {
    if(mode!="cpu" && mode!="gpu") throw std::invalid_argument("Explicit cpu/gpu backend required");
    optuna_cuda::install_runtime(optuna_cuda::cuda_runtime_hooks());
    optuna_cuda::set_gpu(mode=="gpu");
}
}
PYBIND11_MODULE(_optuna_kernels,m) {
    m.doc()="Private deterministic kernel boundary; host management/RNG are not implemented here.";
    m.add_object("_cache_lifetime",py::capsule(new int(0),[](void* token){
        optuna_cuda::release_cached_memory();delete static_cast<int*>(token);
    }));
    py::class_<optuna_cuda::ModelDimension>(m,"Dimension")
        .def(py::init<>()).def_readwrite("low",&optuna_cuda::ModelDimension::low)
        .def_readwrite("high",&optuna_cuda::ModelDimension::high)
        .def_readwrite("step",&optuna_cuda::ModelDimension::step)
        .def_readwrite("minimum",&optuna_cuda::ModelDimension::minimum)
        .def_readwrite("choices",&optuna_cuda::ModelDimension::choices)
        .def_readwrite("categorical",&optuna_cuda::ModelDimension::categorical)
        .def_readwrite("log_scale",&optuna_cuda::ModelDimension::log_scale);
    py::class_<Model,std::shared_ptr<Model>>(m,"Model")
        .def("close",[](Model& x){x.result.device.reset();})
        .def_property_readonly("mu",[](const Model& x){return x.result.mu;})
        .def_property_readonly("sigma",[](const Model& x){return x.result.sigma;})
        .def_property_readonly("probabilities",[](const Model& x){return x.result.probabilities;});
    m.def("log_pdf",[](std::string mode,std::vector<double> samples,py::list specs){
        backend(mode);Inputs inputs(specs);inputs.validate(samples);
        return optuna_cuda::log_pdf(samples,inputs.views);
    });
    m.def("best",[](std::string mode,std::vector<double> candidates,py::list below,py::list above){
        backend(mode);Inputs lower(below),upper(above);
        if(lower.views.size()!=upper.views.size()) throw std::invalid_argument("Candidate dimensions");
        if(candidates.empty()) throw std::invalid_argument("Empty candidates");
        lower.validate(candidates);upper.validate(candidates);
        return optuna_cuda::best(candidates,lower.views,upper.views);
    });
    m.def("build_model",[](std::string mode,std::vector<double> observations,
            std::vector<double> weights,std::vector<optuna_cuda::ModelDimension> dimensions,
            double prior,bool endpoints,bool host_sampling){
        backend(mode);auto model=std::make_shared<Model>();
        model->weights=weights;model->dimensions=dimensions;
        model->result=optuna_cuda::build_model(observations,weights,dimensions,prior,endpoints,host_sampling);
        return model;
    });
    m.def("best_prepared",[](std::string mode,std::vector<double> candidates,
            const Model& below,const Model& above){
        if(mode!="gpu") throw std::invalid_argument("Prepared device models require explicit GPU");
        if(!below.result.device || !above.result.device)
            throw std::invalid_argument("Prepared device model is closed");
        backend(mode);return optuna_cuda::best(candidates,below.views(),above.views());
    });
    m.def("release_cached_memory",&optuna_cuda::release_cached_memory);
    m.def("shutdown",&optuna_cuda::shutdown_runtime);
    m.def("statistics",[]{
        const auto s=optuna_cuda::statistics();py::dict d;
        d["build_route"]="nvcc-shared-library";d["cuda_device"]=s.device;
        d["density_launches"]=s.density_launches;d["selection_launches"]=s.selection_launches;
        d["model_build_launches"]=s.model_build_launches;d["model_reuses"]=s.model_reuses;
        d["host_to_device_bytes"]=s.host_to_device_bytes;d["device_to_host_bytes"]=s.device_to_host_bytes;
        d["cached_device_bytes"]=s.cached_device_bytes;d["taskmanager"]=false;
        return d;
    });
}
