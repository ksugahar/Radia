#include "hdiv_resident_python_api.hpp"
#include "hdiv_inverse_bh.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <climits>
#include <cstdint>
#include <algorithm>
namespace py=pybind11;
using namespace radia_hdiv_python;
template<class T>using Array=py::array_t<T,py::array::c_style>;
namespace {
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
    // Geometry rows: 12 (affine vertices), 30 (curved P2 TET nodes) or 81 (HEX Q2 lattice) coordinates per element.
    if(vertices.ndim()!=2||(vertices.shape(1)!=12&&vertices.shape(1)!=30&&vertices.shape(1)!=81))
        throw py::value_error("Evaluation geometry requires (elements,12) affine vertices, (elements,30) P2 nodes or (elements,81) HEX Q2 nodes");
    input.nodes=std::size_t(vertices.shape(1))/3;
    input.rule=rows(rule,4,bytes,"rule");input.shapes=rows(shapes,0,bytes,"shapes");input.vertices=rows(vertices,0,bytes,"vertices");
    input.elementClass=rows(elementClass,0,bytes,"class");input.elementDofs=rows(elementDofs,0,bytes,"dofs");input.elementPoints=rows(elementPoints,0,bytes,"points");
    return input;
}
std::shared_ptr<radia_hdiv_bh::MaterialOperator> constructFlatMaterial(Array<double> rule,Array<double> shapes,Array<double> vertices,
    Array<std::int32_t> elementClass,Array<std::int32_t> elementDofs,Array<std::int32_t> elementPoints,std::size_t dofs,py::dict lawData,const std::string& backend){
    if(backend!="gpu")throw py::value_error("Material context requires explicit backend='gpu'");
    std::size_t bytes=0;auto law=lawFromDict(lawData,bytes);
    auto input=flatEvaluationInput(rule,shapes,vertices,elementClass,elementDofs,elementPoints,dofs,bytes);
    py::gil_scoped_release release;
    auto owner=radia_hdiv_bh::createFlatMaterial(input,law);return owner;
}
std::shared_ptr<radia_hdiv_bh::MassOperator> constructFlatMass(Array<double> rule,Array<double> shapes,Array<double> vertices,
    Array<std::int32_t> elementClass,Array<std::int32_t> elementDofs,Array<std::int32_t> elementPoints,std::size_t dofs,const std::string& backend){
    if(backend!="gpu")throw py::value_error("Mass context requires explicit backend='gpu'");
    std::size_t bytes=0;auto input=flatEvaluationInput(rule,shapes,vertices,elementClass,elementDofs,elementPoints,dofs,bytes);
    py::gil_scoped_release release;
    auto owner=radia_hdiv_bh::createFlatMass(input);return owner;
}
py::dict validateFlatTables(Array<double> rule,Array<double> shapes,Array<double> vertices,
    Array<std::int32_t> elementClass,Array<std::int32_t> elementDofs,Array<std::int32_t> elementPoints,std::size_t dofs){
    std::size_t bytes=0;const auto budget=radia_hdiv_bh::validateFlatEvaluationScope(flatEvaluationInput(rule,shapes,vertices,elementClass,elementDofs,elementPoints,dofs,bytes));
    py::dict terms;for(std::size_t i=0;i<budget.count;++i)terms[radia_hdiv_bh::FlatEvaluationBudget::names[i]]=budget.bytes[i];
    py::dict out;out["schema"]="radia.hdiv-flat-material-budget.v1";out["terms"]=terms;out["total_bytes"]=budget.total();
    out["limiting_term"]=radia_hdiv_bh::FlatEvaluationBudget::names[budget.limiting()];
    out["excludes"]="CUB scratch and the inverse-BH law; both are added at construction";return out;
}
std::shared_ptr<DeviceOperator> constructAnalyticGram(Array<double> cells,Array<double> faces,const std::string& backend) {
    if(backend!="gpu")throw py::value_error("The device analytic Gram requires explicit backend='gpu'");
    std::size_t bytes=0;const auto c=snapshot<double>(cells,bytes);const auto f=snapshot<double>(faces,bytes);
    py::gil_scoped_release release;return createAnalyticGram(c,f);
}
std::shared_ptr<DeviceOperator> constructFlat(Array<double> cells,Array<double> faces,Array<std::int32_t> host,
    Array<std::int32_t> kind,Array<std::int32_t> expo,Array<double> tet,Array<double> tri,Array<double> tetLow,
    Array<double> triLow,double farFactor,std::size_t dofs,Array<std::int32_t> cellClass,Array<std::int32_t> cellDofs,
    Array<double> cellTables,Array<std::int32_t> cellExponents,Array<std::int32_t> faceClass,Array<std::int32_t> faceDofs,
    Array<double> faceTables,Array<std::int32_t> faceExponents,const std::string& compression,py::object aca,const std::string& backend,
    const std::string& family,py::object cellNumbers,py::object faceNumbers,py::object curveGl,py::object curveGw,
    py::object cellMassRule,py::object faceMassRule,py::object hexRules,py::object hexBands) {
    if(backend!="gpu")throw py::value_error("The device charge Gram requires explicit backend='gpu'");
    if(compression!="dense"&&compression!="aca")throw py::value_error("compression must be explicit 'dense' or 'aca'");
    if((compression=="dense")!=aca.is_none())throw py::value_error("ACA parameters are required for compression='aca' and rejected for 'dense'");
    const auto n=std::size_t(host.size());
    if(host.ndim()!=1||n<1||n>std::size_t(INT_MAX)||!dofs||dofs>std::size_t(INT_MAX)||kind.ndim()!=1||std::size_t(kind.size())!=n||
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
    if(family!="affine-tet"&&family!="curved-tet"&&family!="hex")throw py::value_error("family must be explicit 'affine-tet', 'curved-tet' or 'hex'");
    const bool curved=family=="curved-tet",hex=family=="hex";
    for(const auto* value:{&cellNumbers,&faceNumbers,&curveGl,&curveGw,&cellMassRule,&faceMassRule})
        if(value->is_none()==curved)throw py::value_error("Curved geometry arrays are required for family='curved-tet' and rejected otherwise");
    for(const auto* value:{&hexRules,&hexBands})
        if(value->is_none()==hex)throw py::value_error("hex_rules and hex_bands are required for family='hex' and rejected otherwise");
    FlatInput input{snapshot(cells,bytes),snapshot(faces,bytes),flat(tet),flat(tri),flat(tetLow),flat(triLow),
        integers(host),integers(kind),integers(expo),farFactor,compression};
    input.family=family;
    auto typed=[](const py::object& value,const char* what,auto tag){using T=decltype(tag);
        if(!py::isinstance<Array<T>>(value))throw py::value_error(std::string(what)+" requires a contiguous array of the declared dtype");
        return py::cast<Array<T>>(value);};
    if(curved){
        input.cellNumbers=integers(typed(cellNumbers,"cell_numbers",std::int32_t{}));input.faceNumbers=integers(typed(faceNumbers,"face_numbers",std::int32_t{}));
        input.gl=flat(typed(curveGl,"curve_gl",double{}));input.gw=flat(typed(curveGw,"curve_gw",double{}));
    }
    if(hex){input.hexRules=flat(typed(hexRules,"hex_rules",double{}));input.hexBands=flat(typed(hexBands,"hex_bands",double{}));}
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
    charges.hex=hex;
    if(curved){charges.curved=true;
        auto rule=[&](const py::object& value,const char* what){auto a=typed(value,what,double{});
            if(a.ndim()!=2||a.shape(1)!=4)throw py::value_error(std::string(what)+" requires (points,4) float64 rows");return flat(a);};
        charges.cellRule=rule(cellMassRule,"cell_mass_rule");charges.faceRule=rule(faceMassRule,"face_mass_rule");}
    py::gil_scoped_release release;auto owner=createFlat(input,charges);return owner;
}
}
void registerMaterialStaging(py::module_&);
static int newtonMode(const std::string& preconditioner){
    if(preconditioner=="reference-diagonal")return 0;if(preconditioner=="reference-chebyshev")return 1;if(preconditioner=="reference-chebyshev-graph")return 2;
    throw py::value_error("Unsupported explicit Newton reference preconditioner");
}
static int massMode(const std::string& preconditioner){
    if(preconditioner=="mass-diagonal")return 0;if(preconditioner=="mass-cholesky")return 1;if(preconditioner=="mass-polynomial")return 2;
    if(preconditioner=="mass-polynomial-graph")return 3;if(preconditioner=="mass-chebyshev")return 4;if(preconditioner=="mass-chebyshev-graph")return 5;
    throw py::value_error("Unsupported material preconditioner");
}
static py::dict publishNewton(const NewtonResult& result){
    py::array_t<double> solution(result.x.size());std::copy(result.x.begin(),result.x.end(),solution.mutable_data());py::dict out;out["x"]=solution;out["relative_residual"]=result.residual;out["energy"]=result.energy;out["newton_iterations"]=result.history.size();out["preconditioner_graph_launches"]=result.graphLaunches;
    py::list history;for(const auto& s:result.history){py::dict row;row["relative_residual"]=s.residual;row["energy"]=s.energy;row["alpha"]=s.alpha;row["linear_relative_residual"]=s.linearResidual;row["linear_iterations"]=s.linearIterations;row["backtracks"]=s.backtracks;row["original_residual_certified"]=s.residualCertified;history.append(row);}out["history"]=history;return out;
}
PYBIND11_MODULE(_hdiv_cuda,m) {
    registerMaterialStaging(m);
    py::class_<radia_hdiv_bh::MaterialOperator,std::shared_ptr<radia_hdiv_bh::MaterialOperator>>(m,"_ResidentBHMaterial")
      .def_static("from_flat_tables",&constructFlatMaterial,py::arg("rule").noconvert(),py::arg("reference_shapes").noconvert(),
           py::arg("vertices").noconvert(),py::arg("element_class").noconvert(),py::arg("element_dofs").noconvert(),
           py::arg("element_points").noconvert(),py::arg("dofs"),py::arg("law"),py::arg("backend"),
           "Matrix-free device E and E^T from resident NGSolve reference tables, with point measures and reference mass; no evaluation matrix exists or crosses.")
      .def_static("validate_flat_tables",&validateFlatTables,py::arg("rule").noconvert(),py::arg("reference_shapes").noconvert(),
           py::arg("vertices").noconvert(),py::arg("element_class").noconvert(),py::arg("element_dofs").noconvert(),
           py::arg("element_points").noconvert(),py::arg("dofs"),
           "Host-only scope check of the flat evaluation tables (same checks as from_flat_tables) and their device byte budget; no device context.")
      .def("assembly",[](radia_hdiv_bh::MaterialOperator& self){std::string route;radia_hmatrix_cuda::HostVector<radia_hdiv_bh::LedgerRow> rows;radia_hdiv_bh::MaterialFootprint bytes;
          {py::gil_scoped_release release;route=self.route();rows=self.transfers();bytes=self.footprint();}
          py::dict d;d["route"]=route;d["persistent_device_bytes"]=bytes.persistent;d["setup_transient_device_bytes"]=bytes.transient;py::list ledger;
          for(const auto& t:rows){py::dict row;row["payload"]=t.payload;row["direction"]=t.direction;row["bytes"]=t.bytes;ledger.append(row);}
          d["transfer_ledger"]=ledger;return d;})
      .def("reference_mass_apply",[](radia_hdiv_bh::MaterialOperator& self,Array<double> x){
          if(x.ndim()!=1)throw py::value_error("Reference mass input must be one-dimensional");
          std::size_t bytes=0;const auto input=snapshot(x,bytes);radia_hdiv_bh::MassAction result;
          {py::gil_scoped_release release;result=self.massAction(input);}
          py::array_t<double> y(result.y.size());std::copy(result.y.begin(),result.y.end(),y.mutable_data());
          py::dict out;out["y"]=y;out["nonzeros"]=result.nonzeros;out["spectral_bound"]=result.bound;return out;
      },py::arg("x").noconvert(),"Diagnostic y=Mx of the device-built reference mass; vectors cross, the matrix never does.")
      .def("evaluation_apply",[](radia_hdiv_bh::MaterialOperator& self,Array<double> x,bool transpose){
          if(x.ndim()!=1)throw py::value_error("Evaluation action input must be one-dimensional");
          std::size_t bytes=0;const auto input=snapshot(x,bytes);radia_hmatrix_cuda::HostVector<double> y;
          {py::gil_scoped_release release;y=self.evaluationAction(input,transpose);}
          py::array_t<double> out(y.size());std::copy(y.begin(),y.end(),out.mutable_data());return out;
      },py::arg("x").noconvert(),py::arg("transpose"),"Diagnostic E x (transpose=False) or E^T y (transpose=True); vectors cross, no evaluation matrix does.")
      .def("evaluate",[](radia_hdiv_bh::MaterialOperator& self,Array<double> state,Array<double> direction){
          const auto stats=self.stats();if(stats.closed)throw py::value_error("Material context is closed");
          if(state.ndim()!=1||std::size_t(state.size())!=stats.dofs||direction.ndim()!=1||(direction.size()!=0&&std::size_t(direction.size())!=stats.dofs))throw py::value_error("Material state or direction mismatch");
          std::size_t bytes=0;const auto x=snapshot(state,bytes),d=snapshot(direction,bytes);radia_hdiv_bh::MaterialEvaluation result;
          {py::gil_scoped_release release;result=self.evaluate(x,d);}
          auto array=[](const auto& v){py::array_t<double> a(v.size());std::copy(v.begin(),v.end(),a.mutable_data());return a;};
          py::dict out;out["gradient"]=array(result.gradient);out["tangent_action"]=array(result.tangentAction);out["energy"]=result.energy;out["maximum_magnetization"]=result.maximumMagnetization;out["device_payload_bound_bytes"]=result.ownedDeviceBytes;out["linearization_epoch"]=result.linearizationEpoch;return out;
      },py::arg("state").noconvert(),py::arg("direction").noconvert())
      .def("linearize",[](radia_hdiv_bh::MaterialOperator& self,Array<double> state){
          const auto stats=self.stats();if(stats.closed||state.ndim()!=1||std::size_t(state.size())!=stats.dofs)throw py::value_error("Material state dimension or liveness mismatch");
          std::size_t bytes=0;const auto x=snapshot(state,bytes);radia_hdiv_bh::MaterialEvaluation result;
          {py::gil_scoped_release release;result=self.linearize(x);}
          py::array_t<double> gradient(result.gradient.size());std::copy(result.gradient.begin(),result.gradient.end(),gradient.mutable_data());py::dict out;out["gradient"]=gradient;out["energy"]=result.energy;out["linearization_epoch"]=result.linearizationEpoch;return out;
      },py::arg("state").noconvert())
      .def("apply_tangent",[](radia_hdiv_bh::MaterialOperator& self,Array<double> direction,std::size_t epoch){
          const auto stats=self.stats();if(stats.closed||direction.ndim()!=1||std::size_t(direction.size())!=stats.dofs)throw py::value_error("Material direction dimension or liveness mismatch");
          std::size_t bytes=0;const auto d=snapshot(direction,bytes);radia_hmatrix_cuda::HostVector<double> result;{py::gil_scoped_release release;result=self.tangent(d,epoch);}
          py::array_t<double> out(result.size());std::copy(result.begin(),result.end(),out.mutable_data());return out;
      },py::arg("direction").noconvert(),py::arg("linearization_epoch"))
      .def("stats",[](radia_hdiv_bh::MaterialOperator& self){const auto s=self.stats();py::dict out;out["dofs"]=s.dofs;out["points"]=s.points;out["device_payload_bound_bytes"]=s.devicePayloadBound;out["setup_upload_bytes"]=s.setupUploadBytes;out["state_upload_bytes"]=s.stateUploadBytes;out["evaluations"]=s.evaluations;out["epoch"]=s.epoch;out["closed"]=s.closed;out["linearization_epoch"]=s.linearizationEpoch;out["linearization_valid"]=s.linearizationValid;out["tangent_applications"]=s.tangentApplications;return out;})
      .def("close",&radia_hdiv_bh::MaterialOperator::close,py::call_guard<py::gil_scoped_release>());

    py::class_<radia_hdiv_bh::MassOperator,std::shared_ptr<radia_hdiv_bh::MassOperator>>(m,"_ResidentMass")
      .def_static("from_flat_tables",&constructFlatMass,py::arg("rule").noconvert(),py::arg("reference_shapes").noconvert(),
           py::arg("vertices").noconvert(),py::arg("element_class").noconvert(),py::arg("element_dofs").noconvert(),
           py::arg("element_points").noconvert(),py::arg("dofs"),py::arg("backend"),
           "Device-built reference mass M from resident NGSolve reference tables for the linear material solve; no mass entry crosses.")
      .def("assembly",[](radia_hdiv_bh::MassOperator& self){radia_hmatrix_cuda::HostVector<radia_hdiv_bh::LedgerRow> rows;radia_hdiv_bh::MaterialFootprint bytes;std::string route;
          {py::gil_scoped_release release;rows=self.transfers();bytes=self.footprint();route=self.route();}
          py::dict d;d["route"]=route;d["persistent_device_bytes"]=bytes.persistent;d["setup_transient_device_bytes"]=bytes.transient;py::list ledger;
          for(const auto& t:rows){py::dict row;row["payload"]=t.payload;row["direction"]=t.direction;row["bytes"]=t.bytes;ledger.append(row);}
          d["transfer_ledger"]=ledger;return d;})
      .def("reference_mass_apply",[](radia_hdiv_bh::MassOperator& self,Array<double> x){
          if(x.ndim()!=1)throw py::value_error("Reference mass input must be one-dimensional");
          std::size_t bytes=0;const auto input=snapshot(x,bytes);radia_hdiv_bh::MassAction result;
          {py::gil_scoped_release release;result=self.massAction(input);}
          py::array_t<double> y(result.y.size());std::copy(result.y.begin(),result.y.end(),y.mutable_data());
          py::dict out;out["y"]=y;out["nonzeros"]=result.nonzeros;out["spectral_bound"]=result.bound;return out;
      },py::arg("x").noconvert(),"Diagnostic y=Mx of the device-built reference mass; vectors cross, the matrix never does.")
      .def("stats",[](radia_hdiv_bh::MassOperator& self){std::size_t dofs;bool closed;{py::gil_scoped_release release;dofs=self.dofs();closed=self.closed();}
          py::dict d;d["dofs"]=dofs;d["closed"]=closed;d["backend"]="gpu";return d;})
      .def("close",&radia_hdiv_bh::MassOperator::close,py::call_guard<py::gil_scoped_release>());

    m.def("_inverse_bh_samples",[](Array<double> knots,Array<double> coefficients,Array<double> base,
        Array<double> amplitudes,double cap,double barrier,double last,double initial,const std::string& backend){
        if(backend!="gpu")throw py::value_error("Inverse BH requires explicit backend='gpu'");
        if(knots.ndim()!=1||knots.size()<2||knots.size()>4096||coefficients.ndim()!=1||coefficients.size()!=4*(knots.size()-1)||
           base.ndim()!=1||base.size()!=knots.size()||amplitudes.ndim()!=1||amplitudes.size()<1||amplitudes.size()>65536)
            throw py::value_error("Inverse BH payload shape or budget exceeded");
        std::size_t bytes=0;radia_hdiv_bh::Law law{snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(knots,bytes),snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(coefficients,bytes),snapshot<double,radia_hmatrix_cuda::HostFaultAllocator<double>>(base,bytes),cap,barrier,last,initial};
        const auto a=snapshot(amplitudes,bytes);radia_hdiv_bh::Samples result;
        {py::gil_scoped_release release;result=radia_hdiv_bh::evaluate(law,a);}
        auto array=[](const auto& v){py::array_t<double> a(v.size());std::copy(v.begin(),v.end(),a.mutable_data());return a;};
        py::dict d;d["secant"]=array(result.secant);d["tangent"]=array(result.tangent);d["energy"]=array(result.energy);
        d["owned_device_payload_bytes"]=result.ownedDeviceBytes;d["backend"]="gpu";return d;
    },py::arg("knots").noconvert(),py::arg("coefficients").noconvert(),py::arg("integral_base").noconvert(),
      py::arg("amplitudes").noconvert(),py::arg("cap"),py::arg("barrier_slope"),py::arg("last_field"),py::arg("initial_reluctivity"),py::arg("backend"));
    m.doc()="Private shared B^T G B CUDA component; explicit shared material PCG; no FE assembly or nonlinear MMM solver";
    m.def("_diagnostic_flat_charge_entries",[](Array<double> cells,Array<double> faces,Array<std::int32_t> host,Array<std::int32_t> kind,
        Array<std::int32_t> expo,Array<double> tet,Array<double> tri,Array<double> tetLow,Array<double> triLow,double farFactor,
        Array<std::uint64_t> rows,Array<std::uint64_t> cols,const std::string& backend,const std::string& family,
        py::object cellNumbers,py::object faceNumbers,py::object curveGl,py::object curveGw,py::object hexRules,py::object hexBands){
        if(backend!="gpu")throw py::value_error("Flat diagnostic entries require explicit backend='gpu'");
        if(family!="affine-tet"&&family!="curved-tet"&&family!="hex")throw py::value_error("family must be explicit 'affine-tet', 'curved-tet' or 'hex'");
        const bool curved=family=="curved-tet",hex=family=="hex";
        for(const auto* value:{&cellNumbers,&faceNumbers,&curveGl,&curveGw})
            if(value->is_none()==curved)throw py::value_error("Curved geometry arrays are required for family='curved-tet' and rejected otherwise");
        for(const auto* value:{&hexRules,&hexBands})
            if(value->is_none()==hex)throw py::value_error("hex_rules and hex_bands are required for family='hex' and rejected otherwise");
        const auto n=std::size_t(host.size());
        if(host.ndim()!=1||n<1||n>std::size_t(INT_MAX)||kind.ndim()!=1||std::size_t(kind.size())!=n||expo.ndim()!=1||std::size_t(expo.size())!=3*n||
           cells.ndim()!=1||faces.ndim()!=1||rows.ndim()!=1||cols.ndim()!=1||rows.size()!=cols.size()||rows.size()<1||rows.size()>65536)
            throw py::value_error("Flat diagnostic payload shape mismatch");
        for(auto* rule:{&tet,&tri,&tetLow,&triLow})
            if(rule->ndim()!=2||rule->shape(1)!=4)throw py::value_error("Reference rules require (points,4) float64 rows");
        std::size_t bytes=0;
        auto flat=[&](Array<double> a){
            const auto count=std::size_t(a.size());
            if(count>(hostSnapshotBudget-bytes)/sizeof(double))throw py::value_error("Resident host snapshot budget exceeded");
            bytes+=count*sizeof(double);return radia_hmatrix_cuda::HostVector<double>(a.data(),a.data()+count);};
        auto integers=[&](Array<std::int32_t> a){const auto v=snapshot(a,bytes);return radia_hmatrix_cuda::HostVector<int>(v.begin(),v.end());};
        FlatInput input{snapshot(cells,bytes),snapshot(faces,bytes),flat(tet),flat(tri),flat(tetLow),flat(triLow),
            integers(host),integers(kind),integers(expo),farFactor,"dense"};
        input.family=family;
        if(curved){input.cellNumbers=integers(py::cast<Array<std::int32_t>>(cellNumbers));input.faceNumbers=integers(py::cast<Array<std::int32_t>>(faceNumbers));
            input.gl=flat(py::cast<Array<double>>(curveGl));input.gw=flat(py::cast<Array<double>>(curveGw));}
        if(hex){input.hexRules=flat(py::cast<Array<double>>(hexRules));input.hexBands=flat(py::cast<Array<double>>(hexBands));}
        const auto r=snapshot(rows,bytes),c=snapshot(cols,bytes);
        radia_hmatrix_cuda::HostVector<std::size_t> rr(r.begin(),r.end()),cc(c.begin(),c.end());radia_hmatrix_cuda::HostVector<double> values;
        {py::gil_scoped_release release;values=flatDiagnosticEntries(input,rr,cc);}
        py::array_t<double> out(values.size());std::copy(values.begin(),values.end(),out.mutable_data());return out;
    },py::arg("cell_vertices").noconvert(),py::arg("face_vertices").noconvert(),py::arg("charge_host").noconvert(),
      py::arg("charge_kind").noconvert(),py::arg("charge_exponents").noconvert(),py::arg("tet_rule").noconvert(),
      py::arg("tri_rule").noconvert(),py::arg("tet_low_rule").noconvert(),py::arg("tri_low_rule").noconvert(),
      py::arg("far_factor"),py::arg("rows").noconvert(),py::arg("cols").noconvert(),py::arg("backend"),
      py::arg("family"),py::arg("cell_numbers")=py::none(),py::arg("face_numbers")=py::none(),
      py::arg("curve_gl")=py::none(),py::arg("curve_gw")=py::none(),py::arg("hex_rules")=py::none(),py::arg("hex_bands")=py::none(),
      "Diagnostic raw (unnormalized) device flat-Gram entries for explicit charge pairs; pairs and values cross, no matrix does.");
    py::class_<DeviceOperator,std::shared_ptr<DeviceOperator>>(m,"_ResidentHDiv")
      .def_static("from_analytic_constant_block",&constructAnalyticGram,py::arg("cell_vertices").noconvert(),py::arg("face_vertices").noconvert(),py::arg("backend"),
           "Fill a Gram-only constant-density block on the device from flat row-major linear tetrahedra and triangles; at most 128 charges, no images, no charge map, no material solve.")
      .def_static("from_flat_charges",&constructFlat,py::arg("cell_vertices").noconvert(),py::arg("face_vertices").noconvert(),
           py::arg("charge_host").noconvert(),py::arg("charge_kind").noconvert(),py::arg("charge_exponents").noconvert(),
           py::arg("tet_rule").noconvert(),py::arg("tri_rule").noconvert(),py::arg("tet_low_rule").noconvert(),
           py::arg("tri_low_rule").noconvert(),py::arg("far_factor"),py::arg("dofs"),py::arg("cell_class").noconvert(),
           py::arg("cell_dofs").noconvert(),py::arg("cell_tables").noconvert(),py::arg("cell_exponents").noconvert(),
           py::arg("face_class").noconvert(),py::arg("face_dofs").noconvert(),py::arg("face_tables").noconvert(),
           py::arg("face_exponents").noconvert(),py::arg("compression"),py::arg("aca"),py::arg("backend"),
           py::arg("family"),py::arg("cell_numbers")=py::none(),py::arg("face_numbers")=py::none(),
           py::arg("curve_gl")=py::none(),py::arg("curve_gw")=py::none(),py::arg("cell_mass_rule")=py::none(),
           py::arg("face_mass_rule")=py::none(),py::arg("hex_rules")=py::none(),py::arg("hex_bands")=py::none(),
           "Generate the flat-TET BDM1/BDM2 charge Gram (explicit 'dense' block or device 'aca' leaves), B and B^T on the device; only geometry, monomials, rules, NGSolve tables and the geometric partition cross.")
      .def("assembly",[](DeviceOperator& self){Assembly a;{py::gil_scoped_release release;a=self.assembly();}
          py::dict d;d["route"]=a.route;d["charges"]=a.charges;d["setup_peak_device_bytes"]=a.setupBytes;
          d["kernel_launches"]=a.launches;d["compression"]=a.compression;py::dict counters;
          for(const auto& c:a.counters)counters[py::str(c.first)]=c.second;d["hmatrix"]=counters;py::list ledger;
          if(a.budget.empty())d["charge_budget"]=py::none();
          else{py::dict terms;std::size_t total=0,top=0;std::string limiting;
              for(const auto& t:a.budget){terms[py::str(t.first)]=t.second;total+=t.second;if(limiting.empty()||t.second>top){top=t.second;limiting=t.first;}}
              py::dict b;b["schema"]="radia.hdiv-flat-charge-budget.v1";b["terms"]=terms;b["total_bytes"]=total;b["limiting_term"]=limiting;
              b["excludes"]="ACA factor storage and the rank-sized plan buffers; their allocation sites check the live reserve";d["charge_budget"]=b;}
          for(const auto& t:a.transfers){py::dict row;row["payload"]=t.payload;row["direction"]=t.direction;row["bytes"]=t.bytes;ledger.append(row);}
          d["transfer_ledger"]=ledger;return d;})
      .def("apply",[](DeviceOperator& self,Array<double> x){
          std::size_t bytes=0;auto input=snapshot(x,bytes);radia_hmatrix_cuda::HostVector<double> output;
          {py::gil_scoped_release release;output=self.apply(input);}
          py::array_t<double> result(output.size());std::copy(output.begin(),output.end(),result.mutable_data());return result;
      },py::arg("x").noconvert())
      .def("solve_nonlinear_material",[](DeviceOperator& self,radia_hdiv_bh::MaterialOperator& material,Array<double> rhs,Array<double> initial,double tolerance,double linearTolerance,int maximumIterations,int maximumLinearIterations,const std::string& preconditioner,int degree,double lower){
          const auto n=self.stats().dimension;
          if(rhs.ndim()!=1||std::size_t(rhs.size())!=n||initial.ndim()!=1||std::size_t(initial.size())!=n)throw py::value_error("Newton payload dimensions exceeded");
          const int mode=newtonMode(preconditioner);
          std::size_t bytes=0;
          const auto b=snapshot(rhs,bytes),x=snapshot(initial,bytes);NewtonOptions options{tolerance,linearTolerance,maximumIterations,maximumLinearIterations,mode,degree,lower};NewtonResult result;
          {py::gil_scoped_release release;result=self.solve_nonlinear(material,b,x,options);}
          return publishNewton(result);
      },py::arg("material"),py::arg("rhs").noconvert(),py::arg("initial").noconvert(),py::arg("nonlinear_tolerance"),py::arg("linear_tolerance"),py::arg("maximum_iterations"),py::arg("maximum_linear_iterations"),py::arg("preconditioner"),py::arg("degree")=32,py::arg("lower_fraction")=.01)
      .def("solve_linear",[](DeviceOperator& self,radia_hdiv_bh::MassOperator& mass,Array<double> rhs,double invChi,double tolerance,int maximumIterations,const std::string& preconditioner,py::object degree,py::object lower){
          const auto n=self.stats().dimension;
          if(rhs.ndim()!=1||std::size_t(rhs.size())!=n)throw py::value_error("Linear material RHS dimension mismatch");
          const int mode=massMode(preconditioner);LinearOptions options{invChi,tolerance,maximumIterations,mode};
          if(!degree.is_none()||!lower.is_none()){
              if(mode<4)throw py::value_error("Chebyshev settings require mass-chebyshev");
              if(!degree.is_none())options.degree=py::cast<int>(degree);if(!lower.is_none())options.lowerFraction=py::cast<double>(lower);}
          std::size_t bytes=0;const auto b=snapshot(rhs,bytes);LinearResult result;
          {py::gil_scoped_release release;result=self.solve_linear(mass,b,options);}
          py::array_t<double> x(result.x.size());std::copy(result.x.begin(),result.x.end(),x.mutable_data());
          py::dict out;out["x"]=x;out["iterations"]=result.iterations;out["relative_residual"]=result.relativeResidual;
          out["solver_owned_device_bytes"]=result.ownedBytes;out["preconditioner"]=preconditioner;out["preconditioner_graph_launches"]=result.graphLaunches;
          out["chebyshev_degree"]=mode>=4?options.degree:0;out["chebyshev_lower_fraction"]=mode>=4?options.lowerFraction:0.;return out;
      },py::arg("mass"),py::arg("rhs").noconvert(),py::arg("inv_chi"),py::arg("tolerance"),py::arg("maximum_iterations"),py::arg("preconditioner"),
        py::arg("degree")=py::none(),py::arg("lower_fraction")=py::none(),
        "Linear material PCG on (B^T G B + inv_chi M) x = rhs with the device-built mass context; only vectors cross.")
      .def("stats",[](DeviceOperator& self){Stats s;{py::gil_scoped_release release;s=self.stats();}
          py::dict d;d["dimension"]=s.dimension;d["owned_device_bytes"]=s.owned_bytes;d["applications"]=s.applications;d["closed"]=s.closed;d["backend"]="gpu";return d;})
      .def("close",&DeviceOperator::close,py::call_guard<py::gil_scoped_release>());
}
