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
radia_hdiv_bh::Law lawFromDict(py::dict lawData,std::size_t& bytes){
    for(const char* key:{"knots","coefficients","integral_base","cap","barrier_slope","last_field","initial_reluctivity"})if(!lawData.contains(key))throw py::value_error("Incomplete inverse BH law");
    if(!py::isinstance<Array<double>>(lawData["knots"])||!py::isinstance<Array<double>>(lawData["coefficients"])||!py::isinstance<Array<double>>(lawData["integral_base"]))throw py::value_error("Inverse BH arrays require contiguous float64");
    auto knots=py::cast<Array<double>>(lawData["knots"]),coefficients=py::cast<Array<double>>(lawData["coefficients"]),base=py::cast<Array<double>>(lawData["integral_base"]);
    if(knots.ndim()!=1||knots.size()<2||knots.size()>4096||coefficients.ndim()!=1||coefficients.size()!=4*(knots.size()-1)||base.ndim()!=1||base.size()!=knots.size())throw py::value_error("Inverse BH law shape exceeded");
    return radia_hdiv_bh::Law{snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(knots,bytes),snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(coefficients,bytes),snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(base,bytes),py::cast<double>(lawData["cap"]),py::cast<double>(lawData["barrier_slope"]),py::cast<double>(lawData["last_field"]),py::cast<double>(lawData["initial_reluctivity"])};
}
template<class T>radia_hmatrix_cuda::HostVector<T> rows(Array<T> a,std::size_t width,std::size_t& bytes,const char* what){
    if(a.ndim()<1||(width&&(a.ndim()!=2||std::size_t(a.shape(1))!=width)))throw py::value_error(std::string("Flat evaluation table shape: ")+what);
    // C-contiguous tables of any rank are copied flat, row-major.
    const auto count=std::size_t(a.size());
    if(count>(hostSnapshotBudget-bytes)/sizeof(T))throw py::value_error("Resident host snapshot budget exceeded");
    bytes+=count*sizeof(T);return radia_hmatrix_cuda::HostVector<T>(a.data(),a.data()+count);
}
radia_hdiv_bh::FlatEvaluationInput flatEvaluationInput(Array<double> rule,Array<double> shapes,Array<double> vertices,
    Array<std::int32_t> elementClass,Array<std::int32_t> elementDofs,Array<std::int32_t> elementPoints,std::size_t dofs,std::size_t& bytes){
    if(shapes.ndim()!=4||shapes.shape(3)!=3||elementDofs.ndim()!=2||elementPoints.ndim()!=2||elementClass.ndim()!=1)
        throw py::value_error("Flat evaluation tables require (classes,points,local,3) shapes and (elements,n) topology");
    radia_hdiv_bh::FlatEvaluationInput input;
    input.elements=std::size_t(elementClass.size());input.local=std::size_t(elementDofs.shape(1));input.rulePoints=std::size_t(elementPoints.shape(1));
    input.classes=std::size_t(shapes.shape(0));input.points=input.elements*input.rulePoints;input.dofs=dofs;
    if(std::size_t(shapes.shape(1))!=input.rulePoints||std::size_t(shapes.shape(2))!=input.local||std::size_t(elementDofs.shape(0))!=input.elements||std::size_t(elementPoints.shape(0))!=input.elements)
        throw py::value_error("Flat evaluation table dimensions disagree");
    input.rule=rows(rule,4,bytes,"rule");input.shapes=rows(shapes,0,bytes,"shapes");input.vertices=rows(vertices,12,bytes,"vertices");
    input.elementClass=rows(elementClass,0,bytes,"class");input.elementDofs=rows(elementDofs,0,bytes,"dofs");input.elementPoints=rows(elementPoints,0,bytes,"points");
    return input;
}
template<bool Oracle>std::shared_ptr<radia_hdiv_bh::MaterialOperator> constructFlatMaterial(Array<double> rule,Array<double> shapes,Array<double> vertices,
    Array<std::int32_t> elementClass,Array<std::int32_t> elementDofs,Array<std::int32_t> elementPoints,std::size_t dofs,py::dict lawData,const std::string& backend){
    if(backend!="gpu")throw py::value_error("Material context requires explicit backend='gpu'");
    std::size_t bytes=0;auto law=lawFromDict(lawData,bytes);
    auto input=flatEvaluationInput(rule,shapes,vertices,elementClass,elementDofs,elementPoints,dofs,bytes);
    py::gil_scoped_release release;
    auto owner=Oracle?radia_hdiv_bh::createFlatMaterialOracle(input,law):radia_hdiv_bh::createFlatMaterial(input,law);pythonPublicationProbe();return owner;
}
py::dict validateFlatTables(Array<double> rule,Array<double> shapes,Array<double> vertices,
    Array<std::int32_t> elementClass,Array<std::int32_t> elementDofs,Array<std::int32_t> elementPoints,std::size_t dofs){
    std::size_t bytes=0;const auto budget=radia_hdiv_bh::validateFlatEvaluationScope(flatEvaluationInput(rule,shapes,vertices,elementClass,elementDofs,elementPoints,dofs,bytes));
    py::dict terms;for(std::size_t i=0;i<budget.count;++i)terms[radia_hdiv_bh::FlatEvaluationBudget::names[i]]=budget.bytes[i];
    py::dict out;out["schema"]="radia.hdiv-flat-material-budget.v1";out["terms"]=terms;out["total_bytes"]=budget.total();
    out["limiting_term"]=radia_hdiv_bh::FlatEvaluationBudget::names[budget.limiting()];
    out["excludes"]="CUB scratch and the inverse-BH law; both are added at construction";return out;
}
std::shared_ptr<DeviceOperator> constructFlat(Array<double> cells,Array<double> faces,Array<std::int32_t> host,
    Array<std::int32_t> kind,Array<std::int32_t> expo,Array<double> tet,Array<double> tri,Array<double> tetLow,
    Array<double> triLow,double farFactor,std::size_t dofs,Array<std::int32_t> cellClass,Array<std::int32_t> cellDofs,
    Array<double> cellTables,Array<std::int32_t> cellExponents,Array<std::int32_t> faceClass,Array<std::int32_t> faceDofs,
    Array<double> faceTables,Array<std::int32_t> faceExponents,const std::string& compression,py::object aca,const std::string& backend) {
    if(backend!="gpu")throw py::value_error("The device charge Gram requires explicit backend='gpu'");
    if(compression!="dense"&&compression!="aca")throw py::value_error("compression must be explicit 'dense' or 'aca'");
    if((compression=="dense")!=aca.is_none())throw py::value_error("ACA parameters are required for compression='aca' and rejected for 'dense'");
    const auto n=std::size_t(host.size());
    if(host.ndim()!=1||n<1||n>4608||!dofs||dofs>16384||kind.ndim()!=1||std::size_t(kind.size())!=n||
       expo.ndim()!=1||std::size_t(expo.size())!=3*n)
        throw py::value_error("Flat charge payload shape mismatch");
    for(auto* rule:{&tet,&tri,&tetLow,&triLow})
        if(rule->ndim()!=2||rule->shape(1)!=4)throw py::value_error("Reference rules require (points,4) float64 rows");
    if(cellClass.ndim()!=1||faceClass.ndim()!=1||cellDofs.ndim()!=2||faceDofs.ndim()!=2||cellTables.ndim()!=3||faceTables.ndim()!=3||
       cellExponents.ndim()!=2||cellExponents.shape(1)!=3||faceExponents.ndim()!=2||faceExponents.shape(1)!=2||
       cellDofs.shape(0)!=cellClass.size()||faceDofs.shape(0)!=faceClass.size()||cellTables.shape(1)!=cellExponents.shape(0)||
       faceTables.shape(1)!=faceExponents.shape(0)||cellTables.shape(2)!=cellDofs.shape(1)||faceTables.shape(2)!=faceDofs.shape(1))
        throw py::value_error("Flat charge tables require (classes,monomials,width) tables and (hosts,width) DOFs");
    std::size_t bytes=0;
    auto flat=[&](Array<double> a){
        const auto count=std::size_t(a.size());
        if(count>(hostSnapshotBudget-bytes)/sizeof(double))throw py::value_error("Resident host snapshot budget exceeded");
        bytes+=count*sizeof(double);return radia_hmatrix_cuda::HostVector<double>(a.data(),a.data()+count);};
    auto integers=[&](Array<std::int32_t> a){const auto v=snapshot(a,bytes);return radia_hmatrix_cuda::HostVector<int>(v.begin(),v.end());};
    auto table=[&](Array<std::int32_t> a,const char* what){const auto v=rows(a,0,bytes,what);return radia_hmatrix_cuda::HostVector<int>(v.begin(),v.end());};
    if(cells.ndim()!=1||faces.ndim()!=1)throw py::value_error("Geometry uses flat row-major coordinate vectors");
    FlatInput input{snapshot(cells,bytes),snapshot(faces,bytes),flat(tet),flat(tri),flat(tetLow),flat(triLow),
        integers(host),integers(kind),integers(expo),farFactor,compression};
    if(compression=="aca"){
        if(!py::isinstance<py::dict>(aca))throw py::value_error("aca must be a dictionary");
        auto d=py::reinterpret_borrow<py::dict>(aca);
        if(d.size()!=5)throw py::value_error("ACA schema requires eps,eta,leaf_size,max_rank,rank_cap");
        for(const char* key:{"eps","eta","leaf_size","max_rank","rank_cap"})if(!d.contains(key))throw py::value_error("Missing ACA parameter");
        if(!py::isinstance<py::float_>(d["eps"])||!py::isinstance<py::float_>(d["eta"]))throw py::value_error("ACA eps and eta require float values");
        if(!py::isinstance<py::str>(d["rank_cap"]))throw py::value_error("ACA rank_cap must be 'dense' or 'fail'");
        const auto cap=py::cast<std::string>(d["rank_cap"]);
        if(cap!="dense"&&cap!="fail")throw py::value_error("ACA rank_cap must be 'dense' or 'fail'");
        input.acaEps=py::cast<double>(d["eps"]);input.acaEta=py::cast<double>(d["eta"]);
        input.acaLeafSize=integer(d,"leaf_size");input.acaMaxRank=integer(d,"max_rank");input.acaRankCapDense=cap=="dense";
    }
    ChargeMapInput charges;charges.dofs=dofs;
    charges.cellMonomials=std::size_t(cellExponents.shape(0));charges.faceMonomials=std::size_t(faceExponents.shape(0));
    charges.cellWidth=std::size_t(cellDofs.shape(1));charges.faceWidth=std::size_t(faceDofs.shape(1));
    charges.cellTables=rows(cellTables,0,bytes,"cell_tables");charges.faceTables=rows(faceTables,0,bytes,"face_tables");
    charges.cellClass=integers(cellClass);charges.faceClass=integers(faceClass);
    charges.cellDofs=table(cellDofs,"cell_dofs");charges.faceDofs=table(faceDofs,"face_dofs");
    charges.cellExponents=table(cellExponents,"cell_exponents");charges.faceExponents=table(faceExponents,"face_exponents");
    py::gil_scoped_release release;auto owner=createFlat(input,charges);pythonPublicationProbe();return owner;
}
}
void registerMaterialStaging(py::module_&);
static int newtonMode(const std::string& preconditioner){
    if(preconditioner=="reference-diagonal")return 0;if(preconditioner=="reference-chebyshev")return 1;if(preconditioner=="reference-chebyshev-graph")return 2;
    throw py::value_error("Unsupported explicit Newton reference preconditioner");
}
static py::dict publishNewton(const NewtonResult& result){
    pythonPublicationProbe();py::array_t<double> solution(result.x.size());std::copy(result.x.begin(),result.x.end(),solution.mutable_data());pythonPublicationProbe();py::dict out;out["x"]=solution;out["relative_residual"]=result.residual;out["energy"]=result.energy;out["newton_iterations"]=result.history.size();out["preconditioner_graph_launches"]=result.graphLaunches;
    pythonPublicationProbe();py::list history;for(const auto& s:result.history){pythonPublicationProbe();py::dict row;row["relative_residual"]=s.residual;row["energy"]=s.energy;row["alpha"]=s.alpha;row["linear_relative_residual"]=s.linearResidual;row["linear_iterations"]=s.linearIterations;row["backtracks"]=s.backtracks;row["original_residual_certified"]=s.residualCertified;history.append(row);}out["history"]=history;return out;
}
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
        std::size_t bytes=0;auto law=lawFromDict(lawData,bytes);const auto rp=snapshot(offsets,bytes),ci=snapshot(indices,bytes);
        radia_hdiv_bh::MaterialInput input{std::size_t(weights.size()),dofs,{rp.begin(),rp.end()},{ci.begin(),ci.end()},snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(values,bytes),snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(weights,bytes)};

          py::gil_scoped_release release;auto owner=radia_hdiv_bh::createMaterial(input,law);pythonPublicationProbe();return owner;
      }),py::arg("offsets").noconvert(),py::arg("indices").noconvert(),py::arg("values").noconvert(),py::arg("weights").noconvert(),py::arg("dofs"),py::arg("law"),py::arg("backend"))
      .def_static("from_flat_tables",&constructFlatMaterial<false>,py::arg("rule").noconvert(),py::arg("reference_shapes").noconvert(),
           py::arg("vertices").noconvert(),py::arg("element_class").noconvert(),py::arg("element_dofs").noconvert(),
           py::arg("element_points").noconvert(),py::arg("dofs"),py::arg("law"),py::arg("backend"),
           "Matrix-free device E and E^T from resident NGSolve reference tables, with point measures and reference mass; no evaluation matrix exists or crosses.")
      .def_static("_explicit_evaluation_oracle",&constructFlatMaterial<true>,py::arg("rule").noconvert(),py::arg("reference_shapes").noconvert(),
           py::arg("vertices").noconvert(),py::arg("element_class").noconvert(),py::arg("element_dofs").noconvert(),
           py::arg("element_points").noconvert(),py::arg("dofs"),py::arg("law"),py::arg("backend"),
           "Diagnostic only, small cases: the same tables with an explicit device E/E^T, for matrix-free parity checks.")
      .def_static("validate_flat_tables",&validateFlatTables,py::arg("rule").noconvert(),py::arg("reference_shapes").noconvert(),
           py::arg("vertices").noconvert(),py::arg("element_class").noconvert(),py::arg("element_dofs").noconvert(),
           py::arg("element_points").noconvert(),py::arg("dofs"),
           "Host-only scope check of the flat evaluation tables (same checks as from_flat_tables) and their device byte budget; no device context.")
      .def("assembly",[](radia_hdiv_bh::MaterialOperator& self){std::string route;radia_hmatrix_cuda::HostVector<radia_hdiv_bh::LedgerRow> rows;radia_hdiv_bh::MaterialFootprint bytes;
          {py::gil_scoped_release release;route=self.route();rows=self.transfers();bytes=self.footprint();}
          pythonPublicationProbe();py::dict d;d["route"]=route;d["persistent_device_bytes"]=bytes.persistent;d["setup_transient_device_bytes"]=bytes.transient;py::list ledger;
          for(const auto& t:rows){py::dict row;row["payload"]=t.payload;row["direction"]=t.direction;row["bytes"]=t.bytes;ledger.append(row);}
          d["transfer_ledger"]=ledger;return d;})
      .def("reference_mass_apply",[](radia_hdiv_bh::MaterialOperator& self,Array<double> x){
          if(x.ndim()!=1)throw py::value_error("Reference mass input must be one-dimensional");
          std::size_t bytes=0;const auto input=snapshot(x,bytes);radia_hdiv_bh::MassAction result;
          {py::gil_scoped_release release;result=self.massAction(input);}
          pythonPublicationProbe();py::array_t<double> y(result.y.size());std::copy(result.y.begin(),result.y.end(),y.mutable_data());
          py::dict out;out["y"]=y;out["nonzeros"]=result.nonzeros;out["spectral_bound"]=result.bound;return out;
      },py::arg("x").noconvert(),"Diagnostic y=Mx of the device-built reference mass; vectors cross, the matrix never does.")
      .def("evaluation_apply",[](radia_hdiv_bh::MaterialOperator& self,Array<double> x,bool transpose){
          if(x.ndim()!=1)throw py::value_error("Evaluation action input must be one-dimensional");
          std::size_t bytes=0;const auto input=snapshot(x,bytes);radia_hmatrix_cuda::HostVector<double> y;
          {py::gil_scoped_release release;y=self.evaluationAction(input,transpose);}
          pythonPublicationProbe();py::array_t<double> out(y.size());std::copy(y.begin(),y.end(),out.mutable_data());return out;
      },py::arg("x").noconvert(),py::arg("transpose"),"Diagnostic E x (transpose=False) or E^T y (transpose=True); vectors cross, no evaluation matrix does.")
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
      .def_static("from_flat_charges",&constructFlat,py::arg("cell_vertices").noconvert(),py::arg("face_vertices").noconvert(),
           py::arg("charge_host").noconvert(),py::arg("charge_kind").noconvert(),py::arg("charge_exponents").noconvert(),
           py::arg("tet_rule").noconvert(),py::arg("tri_rule").noconvert(),py::arg("tet_low_rule").noconvert(),
           py::arg("tri_low_rule").noconvert(),py::arg("far_factor"),py::arg("dofs"),py::arg("cell_class").noconvert(),
           py::arg("cell_dofs").noconvert(),py::arg("cell_tables").noconvert(),py::arg("cell_exponents").noconvert(),
           py::arg("face_class").noconvert(),py::arg("face_dofs").noconvert(),py::arg("face_tables").noconvert(),
           py::arg("face_exponents").noconvert(),py::arg("compression"),py::arg("aca"),py::arg("backend"),
           "Generate the flat-TET BDM1/BDM2 charge Gram (explicit 'dense' block or device 'aca' leaves), B and B^T on the device; only geometry, monomials, rules, NGSolve tables and the geometric partition cross.")
      .def("assembly",[](DeviceOperator& self){Assembly a;{py::gil_scoped_release release;a=self.assembly();}
          pythonPublicationProbe();py::dict d;d["route"]=a.route;d["charges"]=a.charges;d["setup_peak_device_bytes"]=a.setupBytes;
          d["kernel_launches"]=a.launches;d["compression"]=a.compression;py::dict counters;
          for(const auto& c:a.counters)counters[py::str(c.first)]=c.second;d["hmatrix"]=counters;py::list ledger;
          for(const auto& t:a.transfers){py::dict row;row["payload"]=t.payload;row["direction"]=t.direction;row["bytes"]=t.bytes;ledger.append(row);}
          d["transfer_ledger"]=ledger;return d;})
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
      .def("solve_nonlinear_material",[](DeviceOperator& self,radia_hdiv_bh::MaterialOperator& material,Array<double> rhs,Array<double> initial,double tolerance,double linearTolerance,int maximumIterations,int maximumLinearIterations,const std::string& preconditioner,int degree,double lower){
          const auto n=self.stats().dimension;
          if(rhs.ndim()!=1||std::size_t(rhs.size())!=n||initial.ndim()!=1||std::size_t(initial.size())!=n)throw py::value_error("Newton payload dimensions exceeded");
          const int mode=newtonMode(preconditioner);
          std::size_t bytes=0;
          const auto b=snapshot(rhs,bytes),x=snapshot(initial,bytes);NewtonOptions options{tolerance,linearTolerance,maximumIterations,maximumLinearIterations,mode,degree,lower};NewtonResult result;
          {py::gil_scoped_release release;result=self.solve_nonlinear(material,b,x,options);}
          return publishNewton(result);
      },py::arg("material"),py::arg("rhs").noconvert(),py::arg("initial").noconvert(),py::arg("nonlinear_tolerance"),py::arg("linear_tolerance"),py::arg("maximum_iterations"),py::arg("maximum_linear_iterations"),py::arg("preconditioner"),py::arg("degree")=32,py::arg("lower_fraction")=.01)
      .def("solve_nonlinear_material_host_reference",[](DeviceOperator& self,radia_hdiv_bh::MaterialOperator& material,Array<std::uint64_t> offsets,Array<std::uint64_t> indices,Array<double> values,Array<double> rhs,Array<double> initial,double tolerance,double linearTolerance,int maximumIterations,int maximumLinearIterations,const std::string& preconditioner,int degree,double lower){
          const auto n=self.stats().dimension;
          if(offsets.ndim()!=1||std::size_t(offsets.size())!=n+1||indices.ndim()!=1||values.ndim()!=1||indices.size()!=values.size()||values.size()>4000000||rhs.ndim()!=1||std::size_t(rhs.size())!=n||initial.ndim()!=1||std::size_t(initial.size())!=n)throw py::value_error("Newton payload dimensions exceeded");
          const int mode=newtonMode(preconditioner);
          std::size_t bytes=0;const auto rp=snapshot(offsets,bytes),ci=snapshot(indices,bytes);MapInput reference{n,n,{rp.begin(),rp.end()},{ci.begin(),ci.end()},snapshot(values,bytes)};
          const auto b=snapshot(rhs,bytes),x=snapshot(initial,bytes);NewtonOptions options{tolerance,linearTolerance,maximumIterations,maximumLinearIterations,mode,degree,lower};NewtonResult result;
          {py::gil_scoped_release release;result=self.solve_nonlinear_host_reference(material,reference,b,x,options);}
          return publishNewton(result);
      },py::arg("material"),py::arg("reference_offsets").noconvert(),py::arg("reference_indices").noconvert(),py::arg("reference_values").noconvert(),py::arg("rhs").noconvert(),py::arg("initial").noconvert(),py::arg("nonlinear_tolerance"),py::arg("linear_tolerance"),py::arg("maximum_iterations"),py::arg("maximum_linear_iterations"),py::arg("preconditioner"),py::arg("degree")=32,py::arg("lower_fraction")=.01,
         "Host-evaluation material contexts only; a device-table context rejects a host reference mass.")
      .def("stats",[](DeviceOperator& self){Stats s;{py::gil_scoped_release release;s=self.stats();}
          pythonPublicationProbe();py::dict d;d["dimension"]=s.dimension;d["owned_device_bytes"]=s.owned_bytes;d["applications"]=s.applications;d["closed"]=s.closed;d["backend"]="gpu";return d;})
      .def("close",&DeviceOperator::close,py::call_guard<py::gil_scoped_release>());
}
