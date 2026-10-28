#include "hdiv_resident_python_api.hpp"
#include "hdiv_inverse_bh.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <cstdint>
#include <algorithm>
namespace py=pybind11;
using namespace radia_hdiv_python;
template<class T>using Array=py::array_t<T,py::array::c_style>;
namespace {
#ifdef RADIA_HDIV_PYTHON_FAULT_INJECTION
thread_local std::size_t publicationFailAt=0,publicationAttempts=0;
#endif
void pythonPublicationProbe(){
#ifdef RADIA_HDIV_PYTHON_FAULT_INJECTION
    ++publicationAttempts;if(publicationFailAt && publicationAttempts==publicationFailAt)throw std::bad_alloc();
#endif
}
constexpr std::size_t hostSnapshotBudget=std::size_t(512)<<20;
template<class T,class Allocator=radia_hmatrix_cuda::HostFaultAllocator<T>>std::vector<T,Allocator> snapshot(Array<T> a,std::size_t& bytes) {
    if(a.ndim()!=1)throw py::value_error("Resident buffers must be one-dimensional and contiguous");
    const auto count=std::size_t(a.size());
    if(count>(hostSnapshotBudget-bytes)/sizeof(T))throw py::value_error("Resident host snapshot budget exceeded");
    bytes+=count*sizeof(T);if(!count)return {};
    return {a.data(),a.data()+count};
}
std::size_t integer(const py::dict& d,const char* key) {
    const auto value=d[py::str(key)];
    if(!py::isinstance<py::int_>(value)||py::isinstance<py::bool_>(value))throw py::value_error("Leaf dimensions require nonnegative integers");
    return py::cast<std::size_t>(value);
}
std::shared_ptr<DeviceOperator> construct(py::list input,Array<std::uint64_t> permutation,
    Array<double> scale,std::size_t dofs,Array<std::uint64_t> offsets,
    Array<std::uint64_t> indices,Array<double> values,const std::string& backend) {
    static_assert(sizeof(std::size_t)==sizeof(std::uint64_t));
    if(backend!="gpu")throw py::value_error("The resident CUDA operator requires explicit backend='gpu'");
    if(scale.ndim()!=1||scale.size()<1||scale.size()>4608||!dofs||dofs>16384||input.size()==0||input.size()>65536)
        throw py::value_error("Unsupported bounded resident dimensions");
    const auto n=std::size_t(scale.size());
    if(permutation.ndim()!=1||permutation.size()!=n||offsets.ndim()!=1||offsets.size()!=n+1||
       indices.ndim()!=1||values.ndim()!=1||indices.size()!=values.size())throw py::value_error("Resident buffer lengths disagree");
    std::size_t bytes=0;const auto perm64=snapshot(permutation,bytes);
    const auto physical=snapshot(scale,bytes),vv=snapshot(values,bytes);
    const auto rp64=snapshot(offsets,bytes),ci64=snapshot(indices,bytes);
    radia_hmatrix_cuda::HostVector<std::size_t> perm(perm64.begin(),perm64.end());
    MapInput map{n,dofs,{rp64.begin(),rp64.end()},{ci64.begin(),ci64.end()},vv};
    radia_hmatrix_cuda::HostVector<LeafInput> leaves;leaves.reserve(input.size());
    for(const auto item:input) {
        if(!py::isinstance<py::dict>(item))throw py::value_error("Each leaf must be a dictionary");
        auto d=py::reinterpret_borrow<py::dict>(item);
        if(d.size()!=7)throw py::value_error("Leaf schema requires row,column,rows,columns,rank,a1,a2");
        for(const char* key:{"row","column","rows","columns","rank","a1","a2"})
            if(!d.contains(key))throw py::value_error("Missing resident leaf field");
        if(!py::isinstance<Array<double>>(d["a1"])||!py::isinstance<Array<double>>(d["a2"]))
            throw py::value_error("Leaf factors require contiguous float64 arrays");
        LeafInput leaf{integer(d,"row"),integer(d,"column"),integer(d,"rows"),integer(d,"columns"),integer(d,"rank"),
            snapshot(py::cast<Array<double>>(d["a1"]),bytes),snapshot(py::cast<Array<double>>(d["a2"]),bytes)};
        leaves.push_back(std::move(leaf));
    }
    // Snapshots own all buffers while another Python thread can mutate inputs.
    py::gil_scoped_release release;auto owner=create(n,leaves,perm,physical,map);pythonPublicationProbe();return owner;
}
}
void registerMaterialStaging(py::module_&);
PYBIND11_MODULE(_hdiv_cuda,m) {
    registerMaterialStaging(m);
#ifdef RADIA_HDIV_PYTHON_FAULT_INJECTION
    m.def("_diagnostic_failure",[](const std::string& kind,std::size_t at){
        using namespace radia_hmatrix_cuda;
        if(kind=="host")configureHostFault(at);
        else if(kind=="device")configureDeviceFault(at);
        else if(kind=="resource")configureResourceFault(at);
        else if(kind=="publication"){publicationFailAt=at;publicationAttempts=0;}
        else throw py::value_error("Unknown diagnostic class");
    });
    m.def("_diagnostic_report",[](){
        using namespace radia_hmatrix_cuda;py::dict out;
        out["host_attempts"]=hostFaultState.attempts;out["host_live_bytes"]=hostFaultState.liveBytes;
        out["device_attempts"]=deviceFaultState.attempts;out["device_live_bytes"]=deviceFaultState.liveBytes;
        out["device_free_errors"]=deviceFaultState.freeErrors;out["resource_attempts"]=resourceFaultState.attempts;
        out["live_resources"]=resourceFaultState.live;out["resource_free_errors"]=resourceFaultState.freeErrors;
        out["publication_attempts"]=publicationAttempts;return out;
    });
#endif
    py::class_<radia_hdiv_bh::MaterialOperator,std::shared_ptr<radia_hdiv_bh::MaterialOperator>>(m,"_ResidentBHMaterial")
      .def(py::init([](Array<std::uint64_t> offsets,Array<std::uint64_t> indices,Array<double> values,Array<double> weights,
          std::size_t dofs,py::dict lawData,const std::string& backend){
          if(backend!="gpu")throw py::value_error("Material context requires explicit backend='gpu'");
          if(weights.ndim()!=1||weights.size()<1||weights.size()>65536||!dofs||dofs>16384||offsets.ndim()!=1||offsets.size()!=3*weights.size()+1||
             indices.ndim()!=1||values.ndim()!=1||indices.size()!=values.size()||values.size()>4000000)throw py::value_error("Mapped material payload shape or budget exceeded");
        for(const char* key:{"knots","coefficients","integral_base","cap","barrier_slope","last_field","initial_reluctivity"})if(!lawData.contains(key))throw py::value_error("Incomplete inverse BH law");
        if(!py::isinstance<Array<double>>(lawData["knots"])||!py::isinstance<Array<double>>(lawData["coefficients"])||!py::isinstance<Array<double>>(lawData["integral_base"]))throw py::value_error("Inverse BH arrays require contiguous float64");
        auto knots=py::cast<Array<double>>(lawData["knots"]),coefficients=py::cast<Array<double>>(lawData["coefficients"]),base=py::cast<Array<double>>(lawData["integral_base"]);
        if(knots.ndim()!=1||knots.size()<2||knots.size()>4096||coefficients.ndim()!=1||coefficients.size()!=4*(knots.size()-1)||base.ndim()!=1||base.size()!=knots.size())throw py::value_error("Inverse BH law shape exceeded");
        std::size_t bytes=0;const auto rp=snapshot(offsets,bytes),ci=snapshot(indices,bytes);
        radia_hdiv_bh::MaterialInput input{std::size_t(weights.size()),dofs,{rp.begin(),rp.end()},{ci.begin(),ci.end()},snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(values,bytes),snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(weights,bytes)};
        radia_hdiv_bh::Law law{snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(knots,bytes),snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(coefficients,bytes),snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(base,bytes),py::cast<double>(lawData["cap"]),py::cast<double>(lawData["barrier_slope"]),py::cast<double>(lawData["last_field"]),py::cast<double>(lawData["initial_reluctivity"])};

          py::gil_scoped_release release;auto owner=radia_hdiv_bh::createMaterial(input,law);pythonPublicationProbe();return owner;
      }),py::arg("offsets").noconvert(),py::arg("indices").noconvert(),py::arg("values").noconvert(),py::arg("weights").noconvert(),py::arg("dofs"),py::arg("law"),py::arg("backend"))
      .def("evaluate",[](radia_hdiv_bh::MaterialOperator& self,Array<double> state,Array<double> direction){
          const auto stats=self.stats();if(stats.closed)throw py::value_error("Material context is closed");
          if(state.ndim()!=1||std::size_t(state.size())!=stats.dofs||direction.ndim()!=1||(direction.size()!=0&&std::size_t(direction.size())!=stats.dofs))throw py::value_error("Material state or direction mismatch");
          std::size_t bytes=0;const auto x=snapshot(state,bytes),d=snapshot(direction,bytes);radia_hdiv_bh::MaterialEvaluation result;
          {py::gil_scoped_release release;result=self.evaluate(x,d);}
          auto array=[](const auto& v){pythonPublicationProbe();py::array_t<double> a(v.size());std::copy(v.begin(),v.end(),a.mutable_data());return a;};
          pythonPublicationProbe();py::dict out;out["gradient"]=array(result.gradient);out["tangent_action"]=array(result.tangentAction);out["energy"]=result.energy;out["maximum_magnetization"]=result.maximumMagnetization;out["device_payload_bound_bytes"]=result.ownedDeviceBytes;out["linearization_epoch"]=result.linearizationEpoch;return out;
      },py::arg("state").noconvert(),py::arg("direction").noconvert())
      .def("linearize",[](radia_hdiv_bh::MaterialOperator& self,Array<double> state){
          const auto stats=self.stats();if(stats.closed||state.ndim()!=1||std::size_t(state.size())!=stats.dofs)throw py::value_error("Material state dimension or liveness mismatch");
          std::size_t bytes=0;const auto x=snapshot(state,bytes);radia_hdiv_bh::MaterialEvaluation result;
          {py::gil_scoped_release release;result=self.linearize(x);}
          pythonPublicationProbe();py::array_t<double> gradient(result.gradient.size());std::copy(result.gradient.begin(),result.gradient.end(),gradient.mutable_data());pythonPublicationProbe();py::dict out;out["gradient"]=gradient;out["energy"]=result.energy;out["linearization_epoch"]=result.linearizationEpoch;return out;
      },py::arg("state").noconvert())
      .def("apply_tangent",[](radia_hdiv_bh::MaterialOperator& self,Array<double> direction,std::size_t epoch){
          const auto stats=self.stats();if(stats.closed||direction.ndim()!=1||std::size_t(direction.size())!=stats.dofs)throw py::value_error("Material direction dimension or liveness mismatch");
          std::size_t bytes=0;const auto d=snapshot(direction,bytes);radia_hmatrix_cuda::HostVector<double> result;{py::gil_scoped_release release;result=self.tangent(d,epoch);}
          pythonPublicationProbe();py::array_t<double> out(result.size());std::copy(result.begin(),result.end(),out.mutable_data());return out;
      },py::arg("direction").noconvert(),py::arg("linearization_epoch"))
      .def("stats",[](radia_hdiv_bh::MaterialOperator& self){const auto s=self.stats();pythonPublicationProbe();py::dict out;out["dofs"]=s.dofs;out["points"]=s.points;out["device_payload_bound_bytes"]=s.devicePayloadBound;out["setup_upload_bytes"]=s.setupUploadBytes;out["state_upload_bytes"]=s.stateUploadBytes;out["evaluations"]=s.evaluations;out["epoch"]=s.epoch;out["closed"]=s.closed;out["linearization_epoch"]=s.linearizationEpoch;out["linearization_valid"]=s.linearizationValid;out["tangent_applications"]=s.tangentApplications;return out;})
      .def("close",&radia_hdiv_bh::MaterialOperator::close,py::call_guard<py::gil_scoped_release>());

    m.def("_material_bh_evaluate",[](Array<std::uint64_t> offsets,Array<std::uint64_t> indices,Array<double> values,
        Array<double> weights,Array<double> state,Array<double> direction,py::dict lawData,const std::string& backend){
        if(backend!="gpu")throw py::value_error("Mapped material requires explicit backend='gpu'");
        if(weights.ndim()!=1||weights.size()<1||weights.size()>65536||state.ndim()!=1||state.size()<1||state.size()>16384||
           offsets.ndim()!=1||offsets.size()!=3*weights.size()+1||indices.ndim()!=1||values.ndim()!=1||indices.size()!=values.size()||values.size()>4000000||
           direction.ndim()!=1||(direction.size()!=0&&direction.size()!=state.size()))throw py::value_error("Mapped material payload shape or budget exceeded");
        for(const char* key:{"knots","coefficients","integral_base","cap","barrier_slope","last_field","initial_reluctivity"})if(!lawData.contains(key))throw py::value_error("Incomplete inverse BH law");
        if(!py::isinstance<Array<double>>(lawData["knots"])||!py::isinstance<Array<double>>(lawData["coefficients"])||!py::isinstance<Array<double>>(lawData["integral_base"]))throw py::value_error("Inverse BH arrays require contiguous float64");
        auto knots=py::cast<Array<double>>(lawData["knots"]),coefficients=py::cast<Array<double>>(lawData["coefficients"]),base=py::cast<Array<double>>(lawData["integral_base"]);
        if(knots.ndim()!=1||knots.size()<2||knots.size()>4096||coefficients.ndim()!=1||coefficients.size()!=4*(knots.size()-1)||base.ndim()!=1||base.size()!=knots.size())throw py::value_error("Inverse BH law shape exceeded");
        std::size_t bytes=0;const auto rp=snapshot(offsets,bytes),ci=snapshot(indices,bytes);
        radia_hdiv_bh::MaterialInput input{std::size_t(weights.size()),std::size_t(state.size()),{rp.begin(),rp.end()},{ci.begin(),ci.end()},snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(values,bytes),snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(weights,bytes)};
        radia_hdiv_bh::Law law{snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(knots,bytes),snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(coefficients,bytes),snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(base,bytes),py::cast<double>(lawData["cap"]),py::cast<double>(lawData["barrier_slope"]),py::cast<double>(lawData["last_field"]),py::cast<double>(lawData["initial_reluctivity"])};
        const auto x=snapshot(state,bytes),d=snapshot(direction,bytes);radia_hdiv_bh::MaterialEvaluation result;
        {py::gil_scoped_release release;result=radia_hdiv_bh::evaluateMaterial(input,law,x,d);}
        auto array=[](const auto& v){pythonPublicationProbe();py::array_t<double> a(v.size());std::copy(v.begin(),v.end(),a.mutable_data());return a;};
        pythonPublicationProbe();py::dict out;out["gradient"]=array(result.gradient);out["tangent_action"]=array(result.tangentAction);out["energy"]=result.energy;out["maximum_magnetization"]=result.maximumMagnetization;out["device_payload_bound_bytes"]=result.ownedDeviceBytes;out["linearization_epoch"]=result.linearizationEpoch;return out;
    },py::arg("offsets").noconvert(),py::arg("indices").noconvert(),py::arg("values").noconvert(),py::arg("weights").noconvert(),py::arg("state").noconvert(),py::arg("direction").noconvert(),py::arg("law"),py::arg("backend"));

    m.def("_inverse_bh_samples",[](Array<double> knots,Array<double> coefficients,Array<double> base,
        Array<double> amplitudes,double cap,double barrier,double last,double initial,const std::string& backend){
        if(backend!="gpu")throw py::value_error("Inverse BH requires explicit backend='gpu'");
        if(knots.ndim()!=1||knots.size()<2||knots.size()>4096||coefficients.ndim()!=1||coefficients.size()!=4*(knots.size()-1)||
           base.ndim()!=1||base.size()!=knots.size()||amplitudes.ndim()!=1||amplitudes.size()<1||amplitudes.size()>65536)
            throw py::value_error("Inverse BH payload shape or budget exceeded");
        std::size_t bytes=0;radia_hdiv_bh::Law law{snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(knots,bytes),snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(coefficients,bytes),snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(base,bytes),cap,barrier,last,initial};
        const auto a=snapshot(amplitudes,bytes);radia_hdiv_bh::Samples result;
        {py::gil_scoped_release release;result=radia_hdiv_bh::evaluate(law,a);}
        auto array=[](const auto& v){pythonPublicationProbe();py::array_t<double> a(v.size());std::copy(v.begin(),v.end(),a.mutable_data());return a;};
        pythonPublicationProbe();py::dict d;d["secant"]=array(result.secant);d["tangent"]=array(result.tangent);d["energy"]=array(result.energy);
        d["owned_device_payload_bytes"]=result.ownedDeviceBytes;d["backend"]="gpu";return d;
    },py::arg("knots").noconvert(),py::arg("coefficients").noconvert(),py::arg("integral_base").noconvert(),
      py::arg("amplitudes").noconvert(),py::arg("cap"),py::arg("barrier_slope"),py::arg("last_field"),py::arg("initial_reluctivity"),py::arg("backend"));
    m.doc()="Private shared B^T G B CUDA component; explicit shared material PCG; no FE assembly or nonlinear MMM solver";
    py::class_<DeviceOperator,std::shared_ptr<DeviceOperator>>(m,"_ResidentHDiv")
      .def(py::init(&construct),py::arg("leaves"),py::arg("permutation").noconvert(),py::arg("scale").noconvert(),
           py::arg("dofs"),py::arg("offsets").noconvert(),py::arg("indices").noconvert(),py::arg("values").noconvert(),py::arg("backend"))
      .def("apply",[](DeviceOperator& self,Array<double> x){
          std::size_t bytes=0;auto input=snapshot(x,bytes);radia_hmatrix_cuda::HostVector<double> output;
          {py::gil_scoped_release release;output=self.apply(input);}
          pythonPublicationProbe();py::array_t<double> result(output.size());std::copy(output.begin(),output.end(),result.mutable_data());return result;
      },py::arg("x").noconvert())
      .def("solve_linear_material",[](DeviceOperator& self,Array<std::uint64_t> offsets,
          Array<std::uint64_t> indices,Array<double> values,Array<double> rhs,double invChi,
          double tolerance,int maxIterations,const std::string& preconditioner,bool graph,int degree,double lower){
          if(preconditioner!="diagonal"&&preconditioner!="chebyshev")throw py::value_error("Unsupported explicit material preconditioner");
          if(graph&&preconditioner!="chebyshev")throw py::value_error("Material graph requires explicit Chebyshev preconditioner");
          const auto dimension=self.stats().dimension;
          if(offsets.ndim()!=1||offsets.size()!=dimension+1||indices.ndim()!=1||values.ndim()!=1||
              indices.size()!=values.size()||rhs.ndim()!=1||rhs.size()!=dimension)
              throw py::value_error("Material mass or RHS dimension mismatch");
          std::size_t bytes=0;const auto rp=snapshot(offsets,bytes),ci=snapshot(indices,bytes);
          MapInput mass{dimension,dimension,{rp.begin(),rp.end()},{ci.begin(),ci.end()},snapshot(values,bytes)};
          auto b=snapshot(rhs,bytes);MaterialResult result;
          {py::gil_scoped_release release;result=self.solve_linear(mass,b,invChi,tolerance,maxIterations,
              preconditioner=="chebyshev",graph,degree,lower);}
          pythonPublicationProbe();py::array_t<double> x(result.x.size());std::copy(result.x.begin(),result.x.end(),x.mutable_data());
          pythonPublicationProbe();py::dict d;d["x"]=x;d["iterations"]=result.iterations;d["relative_residual"]=result.relative_residual;
          d["owned_solver_payload_bytes"]=result.owned_bytes;d["preconditioner_graph_launches"]=result.graph_launches;
          d["backend"]="gpu";d["preconditioner"]=preconditioner;return d;
      },py::arg("mass_offsets").noconvert(),py::arg("mass_indices").noconvert(),py::arg("mass_values").noconvert(),
        py::arg("rhs").noconvert(),py::arg("inv_chi"),py::arg("tolerance"),py::arg("max_iterations"),
        py::arg("preconditioner"),py::arg("graph")=false,py::arg("degree")=32,py::arg("lower_fraction")=.01)
      .def("solve_nonlinear_material",[](DeviceOperator& self,radia_hdiv_bh::MaterialOperator& material,Array<std::uint64_t> offsets,Array<std::uint64_t> indices,Array<double> values,Array<double> rhs,Array<double> initial,double tolerance,double linearTolerance,int maximumIterations,int maximumLinearIterations,const std::string& preconditioner,int degree,double lower){
          const auto n=self.stats().dimension;
          if(offsets.ndim()!=1||std::size_t(offsets.size())!=n+1||indices.ndim()!=1||values.ndim()!=1||indices.size()!=values.size()||values.size()>4000000||rhs.ndim()!=1||std::size_t(rhs.size())!=n||initial.ndim()!=1||std::size_t(initial.size())!=n)throw py::value_error("Newton payload dimensions exceeded");
          int mode;if(preconditioner=="reference-diagonal")mode=0;else if(preconditioner=="reference-chebyshev")mode=1;else if(preconditioner=="reference-chebyshev-graph")mode=2;else throw py::value_error("Unsupported explicit Newton reference preconditioner");
          std::size_t bytes=0;const auto rp=snapshot(offsets,bytes),ci=snapshot(indices,bytes);MapInput reference{n,n,{rp.begin(),rp.end()},{ci.begin(),ci.end()},snapshot(values,bytes)};
          const auto b=snapshot(rhs,bytes),x=snapshot(initial,bytes);NewtonOptions options{tolerance,linearTolerance,maximumIterations,maximumLinearIterations,mode,degree,lower};NewtonResult result;
          {py::gil_scoped_release release;result=self.solve_nonlinear(material,reference,b,x,options);}
          pythonPublicationProbe();py::array_t<double> solution(result.x.size());std::copy(result.x.begin(),result.x.end(),solution.mutable_data());pythonPublicationProbe();py::dict out;out["x"]=solution;out["relative_residual"]=result.residual;out["energy"]=result.energy;out["newton_iterations"]=result.history.size();out["preconditioner_graph_launches"]=result.graphLaunches;
          pythonPublicationProbe();py::list history;for(const auto& s:result.history){pythonPublicationProbe();py::dict row;row["relative_residual"]=s.residual;row["energy"]=s.energy;row["alpha"]=s.alpha;row["linear_relative_residual"]=s.linearResidual;row["linear_iterations"]=s.linearIterations;row["backtracks"]=s.backtracks;row["original_residual_certified"]=s.residualCertified;history.append(row);}out["history"]=history;return out;
      },py::arg("material"),py::arg("reference_offsets").noconvert(),py::arg("reference_indices").noconvert(),py::arg("reference_values").noconvert(),py::arg("rhs").noconvert(),py::arg("initial").noconvert(),py::arg("nonlinear_tolerance"),py::arg("linear_tolerance"),py::arg("maximum_iterations"),py::arg("maximum_linear_iterations"),py::arg("preconditioner"),py::arg("degree")=32,py::arg("lower_fraction")=.01)
      .def("stats",[](DeviceOperator& self){Stats s;{py::gil_scoped_release release;s=self.stats();}
          pythonPublicationProbe();py::dict d;d["dimension"]=s.dimension;d["owned_device_bytes"]=s.owned_bytes;d["applications"]=s.applications;d["closed"]=s.closed;d["backend"]="gpu";return d;})
      .def("close",&DeviceOperator::close,py::call_guard<py::gil_scoped_release>());
}
