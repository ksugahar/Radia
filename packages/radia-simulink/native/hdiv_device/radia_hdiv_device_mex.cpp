// Plain MEX adapter over the shared nvcc-built HDiv device kernels (radia_hdiv_kernels).
// The Gram, charge map, evaluation map and reference mass are generated on the
// device from primitive tables; only vectors and scalar results cross. Inputs are
// the flat C-order arrays of radia.hdiv-device-setup.v1; handles are checked uint64
// tokens, never pointers. No NGSolve, Python or mexcuda dependency.
#include "mex.h"
#include "mex_error_boundary.hpp"
#include "hdiv_resident_python_api.hpp"
#include <cmath>
#include <cstdint>
#include <cstring>
#include <map>
#include <memory>
#include <string>
using namespace radia_hdiv_python;
namespace {
using radia_hmatrix_cuda::HostVector;
struct Entry {
    std::shared_ptr<DeviceOperator> op;std::shared_ptr<radia_hdiv_bh::MassOperator> mass;
    std::shared_ptr<radia_hdiv_bh::MaterialOperator> material;
};
std::map<std::uint64_t,Entry>* registry=nullptr;
std::uint64_t counter=0;
constexpr std::uint64_t operatorTag=0x4844000000000000ULL,massTag=0x484d000000000000ULL,materialTag=0x4842000000000000ULL,tagMask=0xffff000000000000ULL;
void cleanup(){delete registry;registry=nullptr;}

std::string text(const mxArray* a,const char* what){
    if(!mxIsChar(a)||mxGetNumberOfElements(a)>256)throw std::invalid_argument(std::string(what)+" must be text");
    char* raw=mxArrayToString(a);if(!raw)throw std::bad_alloc();std::string s(raw);mxFree(raw);return s;
}
double scalar(const mxArray* a,const char* what){
    if(!mxIsDouble(a)||mxIsComplex(a)||mxIsSparse(a)||mxGetNumberOfElements(a)!=1||!std::isfinite(*mxGetPr(a)))
        throw std::invalid_argument(std::string(what)+" must be a finite real double scalar");
    return *mxGetPr(a);
}
std::size_t count(const mxArray* a,const char* what){
    const double v=scalar(a,what);if(v<0||v!=std::floor(v)||v>9007199254740992.)throw std::invalid_argument(std::string(what)+" must be a nonnegative integer");
    return std::size_t(v);
}
int integerScalar(const mxArray* a,const char* what){
    const double v=scalar(a,what);if(v!=std::floor(v)||std::abs(v)>2147483647.)throw std::invalid_argument(std::string(what)+" must be an int32 value");
    return int(v);
}
const mxArray* field(const mxArray* s,const char* name){
    if(!mxIsStruct(s)||mxGetNumberOfElements(s)!=1)throw std::invalid_argument("Device setup payload must be a scalar struct");
    const mxArray* f=mxGetField(s,0,name);if(!f)throw std::invalid_argument(std::string("Missing device setup field: ")+name);return f;
}
// Flat C-order payloads; MATLAB shape is irrelevant, the class and count are checked.
HostVector<double> reals(const mxArray* a,const char* what){
    if(!mxIsDouble(a)||mxIsComplex(a)||mxIsSparse(a))throw std::invalid_argument(std::string(what)+" must be a real dense double array");
    const auto n=mxGetNumberOfElements(a);const double* p=mxGetPr(a);return HostVector<double>(p,p+n);
}
HostVector<int> integers(const mxArray* a,const char* what){
    if(!mxIsInt32(a)||mxIsComplex(a))throw std::invalid_argument(std::string(what)+" must be an int32 array");
    const auto n=mxGetNumberOfElements(a);const auto* p=static_cast<const std::int32_t*>(mxGetData(a));return HostVector<int>(p,p+n);
}
std::size_t quotient(std::size_t total,std::size_t divisor,const char* what){
    if(!divisor||total%divisor)throw std::invalid_argument(std::string("Inconsistent flat table size: ")+what);return total/divisor;
}
template<class E>E& entry(const mxArray* a,std::uint64_t tag,std::shared_ptr<E> Entry::*member){
    if(!mxIsUint64(a)||mxIsComplex(a)||mxGetNumberOfElements(a)!=1)throw std::invalid_argument("Handle must be a scalar uint64");
    const auto id=*static_cast<const std::uint64_t*>(mxGetData(a));
    if((id&tagMask)!=tag||!registry)throw std::invalid_argument("Handle has the wrong kind or is not live");
    auto it=registry->find(id);if(it==registry->end()||!(it->second.*member))throw std::invalid_argument("Handle has the wrong kind or is not live");
    return *(it->second.*member);
}
mxArray* publish(Entry value,std::uint64_t tag){
    if(counter>=0x0000ffffffffffffULL)throw std::overflow_error("Handle counter exhausted");
    const auto id=tag|++counter;
    mxArray* out=mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);if(!out)throw std::bad_alloc();
    try{if(!registry)registry=new std::map<std::uint64_t,Entry>();if(!registry->emplace(id,std::move(value)).second)throw std::runtime_error("Handle collision");}
    catch(...){mxDestroyArray(out);if(registry&&registry->empty())cleanup();throw;}
    *static_cast<std::uint64_t*>(mxGetData(out))=id;mexLock();mexAtExit(cleanup);return out;
}
mxArray* column(const HostVector<double>& v){mxArray* a=mxCreateDoubleMatrix(v.size(),1,mxREAL);if(!a)throw std::bad_alloc();if(!v.empty())std::memcpy(mxGetPr(a),v.data(),v.size()*sizeof(double));return a;}
template<class V>mxArray* column(const V& v){HostVector<double> copy(v.begin(),v.end());return column(copy);}
mxArray* number(double v){mxArray* a=mxCreateDoubleScalar(v);if(!a)throw std::bad_alloc();return a;}
mxArray* string(const std::string& s){mxArray* a=mxCreateString(s.c_str());if(!a)throw std::bad_alloc();return a;}
// Struct assembled field by field; a failure destroys the partial root.
struct Builder {
    mxArray* root;
    Builder(){root=mxCreateStructMatrix(1,1,0,nullptr);if(!root)throw std::bad_alloc();}
    ~Builder(){if(root)mxDestroyArray(root);}
    Builder& set(const char* name,mxArray* value){if(mxAddField(root,name)<0){mxDestroyArray(value);throw std::bad_alloc();}mxSetField(root,0,name,value);return *this;}
    mxArray* release(){auto* r=root;root=nullptr;return r;}
};
mxArray* ledger(const HostVector<Transfer>& rows){
    const char* names[]={"payload","direction","bytes"};
    mxArray* a=mxCreateStructMatrix(rows.size(),1,3,names);if(!a)throw std::bad_alloc();
    for(std::size_t i=0;i<rows.size();++i){mxSetField(a,i,"payload",string(rows[i].payload));mxSetField(a,i,"direction",string(rows[i].direction));mxSetField(a,i,"bytes",number(double(rows[i].bytes)));}
    return a;
}

radia_hdiv_bh::FlatEvaluationInput evaluationInput(const mxArray* s,std::size_t dofs){
    radia_hdiv_bh::FlatEvaluationInput in;
    const auto rule=reals(field(s,"rule"),"rule"),shapes=reals(field(s,"reference_shapes"),"reference_shapes"),vertices=reals(field(s,"vertices"),"vertices");
    const auto cls=integers(field(s,"element_class"),"element_class"),dofTable=integers(field(s,"element_dofs"),"element_dofs"),points=integers(field(s,"element_points"),"element_points");
    in.elements=cls.size();in.rulePoints=quotient(rule.size(),4,"rule");in.local=quotient(dofTable.size(),in.elements,"element_dofs");
    in.classes=quotient(shapes.size(),in.rulePoints*in.local*3,"reference_shapes");in.points=in.elements*in.rulePoints;in.dofs=dofs;
    // Geometry: 4 affine vertices, 10 curved P2 nodes or a 27-node HEX Q2 lattice per element.
    in.nodes=quotient(vertices.size(),3*in.elements,"vertices");
    if(in.nodes!=4&&in.nodes!=10&&in.nodes!=27)throw std::invalid_argument("Evaluation geometry requires 4 affine vertices, 10 P2 nodes or 27 HEX lattice nodes per element");
    in.rule.assign(rule.begin(),rule.end());in.shapes.assign(shapes.begin(),shapes.end());in.vertices.assign(vertices.begin(),vertices.end());
    in.elementClass.assign(cls.begin(),cls.end());in.elementDofs.assign(dofTable.begin(),dofTable.end());in.elementPoints.assign(points.begin(),points.end());
    return in;
}
radia_hdiv_bh::Law law(const mxArray* s){
    auto vector=[&](const char* name){const auto v=reals(field(s,name),name);return radia_hdiv_bh::BoundaryVector<double>(v.begin(),v.end());};
    // makeHDivBHMaterial returns MATLAB's (knots-1)x4 pchip coefficients; the law is row-major.
    const auto knots=vector("knots");const mxArray* c=field(s,"coefficients");const auto values=reals(c,"coefficients");
    const std::size_t pieces=knots.empty()?0:knots.size()-1;
    if(mxGetNumberOfDimensions(c)!=2||mxGetM(c)!=pieces||mxGetN(c)!=4)throw std::invalid_argument("Inverse BH coefficients must be (knots-1)x4");
    radia_hdiv_bh::BoundaryVector<double> rows(4*pieces);
    for(std::size_t i=0;i<pieces;++i)for(std::size_t k=0;k<4;++k)rows[4*i+k]=values[i+pieces*k];
    return radia_hdiv_bh::Law{knots,rows,vector("integral_base"),scalar(field(s,"cap"),"cap"),
        scalar(field(s,"barrier_slope"),"barrier_slope"),scalar(field(s,"last_field"),"last_field"),scalar(field(s,"initial_reluctivity"),"initial_reluctivity")};
}
std::shared_ptr<DeviceOperator> createOperator(const mxArray* s){
    FlatInput input{reals(field(s,"cell_vertices"),"cell_vertices"),reals(field(s,"face_vertices"),"face_vertices"),
        reals(field(s,"tet_rule"),"tet_rule"),reals(field(s,"tri_rule"),"tri_rule"),reals(field(s,"tet_low_rule"),"tet_low_rule"),
        reals(field(s,"tri_low_rule"),"tri_low_rule"),integers(field(s,"charge_host"),"charge_host"),
        integers(field(s,"charge_kind"),"charge_kind"),integers(field(s,"charge_exponents"),"charge_exponents"),
        scalar(field(s,"far_factor"),"far_factor"),text(field(s,"compression"),"compression")};
    input.family=text(field(s,"family"),"family");
    if(input.family!="affine-tet"&&input.family!="curved-tet"&&input.family!="hex")
        throw std::invalid_argument("family must be explicit 'affine-tet', 'curved-tet' or 'hex'");
    const bool curved=input.family=="curved-tet",hex=input.family=="hex";
    for(const char* name:{"cell_numbers","face_numbers","curve_gl","curve_gw","cell_mass_rule","face_mass_rule"})
        if((mxGetField(s,0,name)!=nullptr)!=curved)throw std::invalid_argument(std::string(name)+" is required for family='curved-tet' and rejected otherwise");
    for(const char* name:{"hex_rules","hex_bands"})
        if((mxGetField(s,0,name)!=nullptr)!=hex)throw std::invalid_argument(std::string(name)+" is required for family='hex' and rejected otherwise");
    if(curved){
        input.cellNumbers=integers(field(s,"cell_numbers"),"cell_numbers");input.faceNumbers=integers(field(s,"face_numbers"),"face_numbers");
        input.gl=reals(field(s,"curve_gl"),"curve_gl");input.gw=reals(field(s,"curve_gw"),"curve_gw");
    }
    if(hex){input.hexRules=reals(field(s,"hex_rules"),"hex_rules");input.hexBands=reals(field(s,"hex_bands"),"hex_bands");}
    if(input.compression=="aca"){
        const mxArray* aca=field(s,"aca");
        input.acaEps=scalar(field(aca,"eps"),"eps");input.acaEta=scalar(field(aca,"eta"),"eta");
        input.acaLeafSize=count(field(aca,"leaf_size"),"leaf_size");input.acaMaxRank=count(field(aca,"max_rank"),"max_rank");
        const auto cap=text(field(aca,"rank_cap"),"rank_cap");
        if(cap!="dense"&&cap!="fail")throw std::invalid_argument("ACA rank_cap must be 'dense' or 'fail'");input.acaRankCapDense=cap=="dense";
    }else if(input.compression!="dense")throw std::invalid_argument("compression must be explicit 'dense' or 'aca'");
    else if(!mxIsEmpty(field(s,"aca")))throw std::invalid_argument("ACA parameters are rejected for compression='dense'");
    ChargeMapInput charges;charges.dofs=count(field(s,"dofs"),"dofs");
    charges.cellClass=integers(field(s,"cell_class"),"cell_class");charges.faceClass=integers(field(s,"face_class"),"face_class");
    charges.cellDofs=integers(field(s,"cell_dofs"),"cell_dofs");charges.faceDofs=integers(field(s,"face_dofs"),"face_dofs");
    charges.cellExponents=integers(field(s,"cell_exponents"),"cell_exponents");charges.faceExponents=integers(field(s,"face_exponents"),"face_exponents");
    charges.cellTables=reals(field(s,"cell_tables"),"cell_tables");charges.faceTables=reals(field(s,"face_tables"),"face_tables");
    charges.hex=hex;
    if(curved){charges.curved=true;
        charges.cellRule=reals(field(s,"cell_mass_rule"),"cell_mass_rule");charges.faceRule=reals(field(s,"face_mass_rule"),"face_mass_rule");
        quotient(charges.cellRule.size(),4,"cell_mass_rule");quotient(charges.faceRule.size(),4,"face_mass_rule");}
    charges.cellMonomials=quotient(charges.cellExponents.size(),3,"cell_exponents");charges.faceMonomials=quotient(charges.faceExponents.size(),2,"face_exponents");
    charges.cellWidth=quotient(charges.cellDofs.size(),charges.cellClass.size(),"cell_dofs");charges.faceWidth=quotient(charges.faceDofs.size(),charges.faceClass.size(),"face_dofs");
    return createFlat(input,charges);
}
int massMode(const std::string& p){
    const char* modes[]={"mass-diagonal","mass-cholesky","mass-polynomial","mass-polynomial-graph","mass-chebyshev","mass-chebyshev-graph"};
    for(int i=0;i<6;++i)if(p==modes[i])return i;throw std::invalid_argument("Unsupported material preconditioner");
}
int newtonMode(const std::string& p){
    if(p=="reference-diagonal")return 0;if(p=="reference-chebyshev")return 1;if(p=="reference-chebyshev-graph")return 2;
    throw std::invalid_argument("Unsupported explicit Newton reference preconditioner");
}
void arity(int nr,int expected,int nl,int outputs){if(nr!=expected||nl>outputs)throw std::invalid_argument("Invalid HDiv device command arity");}

mxArray* materialResult(const radia_hdiv_bh::MaterialEvaluation& r){
    Builder b;
    b.set("gradient",column(r.gradient)).set("tangent_action",column(r.tangentAction)).set("energy",number(r.energy))
     .set("maximum_magnetization",number(r.maximumMagnetization)).set("linearization_epoch",number(double(r.linearizationEpoch)));
    return b.release();
}

void dispatch(int nl,mxArray* out[],int nr,const mxArray* in[]){
    if(nr<1)throw std::invalid_argument("Missing HDiv device command");
    const auto command=text(in[0],"command");
    if(command=="create_operator"){arity(nr,2,nl,1);out[0]=publish(Entry{createOperator(in[1]),{},{}},operatorTag);return;}
    if(command=="create_analytic_gram"){arity(nr,3,nl,1);out[0]=publish(Entry{createAnalyticGram(reals(in[1],"cells"),reals(in[2],"faces")),{},{}},operatorTag);return;}
    if(command=="create_mass"){arity(nr,3,nl,1);out[0]=publish(Entry{{},radia_hdiv_bh::createFlatMass(evaluationInput(in[1],count(in[2],"dofs"))),{}},massTag);return;}
    if(command=="create_material"){arity(nr,4,nl,1);out[0]=publish(Entry{{},{},radia_hdiv_bh::createFlatMaterial(evaluationInput(in[1],count(in[2],"dofs")),law(in[3]))},materialTag);return;}
    if(command=="destroy"){
        arity(nr,2,nl,0);if(!mxIsUint64(in[1])||mxGetNumberOfElements(in[1])!=1||!registry)throw std::invalid_argument("Handle is not live");
        const auto id=*static_cast<const std::uint64_t*>(mxGetData(in[1]));auto it=registry->find(id);
        if(it==registry->end())throw std::invalid_argument("Handle is not live");
        registry->erase(it);if(registry->empty())cleanup();mexUnlock();return;
    }
    if(command=="live_count"){arity(nr,1,nl,1);out[0]=number(registry?double(registry->size()):0.);return;}
    if(command=="apply"){arity(nr,3,nl,1);out[0]=column(entry(in[1],operatorTag,&Entry::op).apply(reals(in[2],"x")));return;}
    if(command=="solve_linear"){
        if(nr!=8&&nr!=10)throw std::invalid_argument("Invalid HDiv device command arity");if(nl>1)throw std::invalid_argument("Invalid HDiv device command arity");
        auto& op=entry(in[1],operatorTag,&Entry::op);auto& mass=entry(in[2],massTag,&Entry::mass);
        const auto preconditioner=text(in[7],"preconditioner");LinearOptions options{scalar(in[4],"inv_chi"),scalar(in[5],"tolerance"),integerScalar(in[6],"maximum_iterations"),massMode(preconditioner)};
        if(nr==10){if(options.preconditioner<4)throw std::invalid_argument("Chebyshev settings require mass-chebyshev");
            options.degree=integerScalar(in[8],"degree");options.lowerFraction=scalar(in[9],"lower_fraction");}
        const auto r=op.solve_linear(mass,reals(in[3],"rhs"),options);Builder b;
        b.set("x",column(r.x)).set("iterations",number(r.iterations)).set("relative_residual",number(r.relativeResidual))
         .set("solver_owned_device_bytes",number(double(r.ownedBytes))).set("preconditioner",string(preconditioner))
         .set("preconditioner_graph_launches",number(double(r.graphLaunches)))
         .set("chebyshev_degree",number(options.preconditioner>=4?options.degree:0)).set("chebyshev_lower_fraction",number(options.preconditioner>=4?options.lowerFraction:0.));
        out[0]=b.release();return;
    }
    if(command=="solve_nonlinear"){
        arity(nr,6,nl,1);auto& op=entry(in[1],operatorTag,&Entry::op);auto& material=entry(in[2],materialTag,&Entry::material);const mxArray* o=in[5];
        const auto preconditioner=text(field(o,"preconditioner"),"preconditioner");
        NewtonOptions options{scalar(field(o,"nonlinear_tolerance"),"nonlinear_tolerance"),scalar(field(o,"linear_tolerance"),"linear_tolerance"),
            integerScalar(field(o,"maximum_iterations"),"maximum_iterations"),integerScalar(field(o,"maximum_linear_iterations"),"maximum_linear_iterations"),
            newtonMode(preconditioner),integerScalar(field(o,"degree"),"degree"),scalar(field(o,"lower_fraction"),"lower_fraction")};
        const auto r=op.solve_nonlinear(material,reals(in[3],"rhs"),reals(in[4],"initial"),options);
        const char* names[]={"relative_residual","energy","alpha","linear_relative_residual","linear_iterations","backtracks","original_residual_certified"};
        mxArray* history=mxCreateStructMatrix(r.history.size(),1,7,names);if(!history)throw std::bad_alloc();
        Builder b;b.set("history",history);
        for(std::size_t i=0;i<r.history.size();++i){const auto& h=r.history[i];
            mxSetField(history,i,names[0],number(h.residual));mxSetField(history,i,names[1],number(h.energy));mxSetField(history,i,names[2],number(h.alpha));
            mxSetField(history,i,names[3],number(h.linearResidual));mxSetField(history,i,names[4],number(h.linearIterations));
            mxSetField(history,i,names[5],number(h.backtracks));mxArray* c=mxCreateLogicalScalar(h.residualCertified);if(!c)throw std::bad_alloc();mxSetField(history,i,names[6],c);}
        b.set("x",column(r.x)).set("relative_residual",number(r.residual)).set("energy",number(r.energy))
         .set("newton_iterations",number(double(r.history.size()))).set("preconditioner_graph_launches",number(double(r.graphLaunches)));
        out[0]=b.release();return;
    }
    if(command=="assembly"){
        arity(nr,2,nl,1);const auto a=entry(in[1],operatorTag,&Entry::op).assembly();Builder b;
        b.set("route",string(a.route)).set("charges",number(double(a.charges))).set("setup_peak_device_bytes",number(double(a.setupBytes)))
         .set("kernel_launches",number(double(a.launches))).set("compression",string(a.compression)).set("transfer_ledger",ledger(a.transfers));
        Builder counters;for(const auto& c:a.counters)counters.set(c.first.c_str(),number(double(c.second)));b.set("hmatrix",counters.release());
        Builder budget;for(const auto& t:a.budget)budget.set(t.first.c_str(),number(double(t.second)));b.set("charge_budget",budget.release());
        out[0]=b.release();return;
    }
    if(command=="stats"){
        arity(nr,2,nl,1);const auto s=entry(in[1],operatorTag,&Entry::op).stats();Builder b;
        b.set("dimension",number(double(s.dimension))).set("owned_device_bytes",number(double(s.owned_bytes))).set("applications",number(double(s.applications)))
         .set("backend",string("gpu"));out[0]=b.release();return;
    }
    if(command=="mass_apply"){
        arity(nr,3,nl,1);const auto r=entry(in[1],massTag,&Entry::mass).massAction(reals(in[2],"x"));Builder b;
        b.set("y",column(r.y)).set("nonzeros",number(double(r.nonzeros))).set("spectral_bound",number(r.bound));out[0]=b.release();return;
    }
    if(command=="material_evaluate"){
        arity(nr,4,nl,1);out[0]=materialResult(entry(in[1],materialTag,&Entry::material).evaluate(reals(in[2],"state"),reals(in[3],"direction")));return;
    }
    if(command=="material_linearize"){arity(nr,3,nl,1);out[0]=materialResult(entry(in[1],materialTag,&Entry::material).linearize(reals(in[2],"state")));return;}
    if(command=="material_apply_tangent"){arity(nr,4,nl,1);out[0]=column(entry(in[1],materialTag,&Entry::material).tangent(reals(in[2],"direction"),count(in[3],"epoch")));return;}
    if(command=="inverse_bh_samples"){
        arity(nr,3,nl,1);const auto amplitudes=reals(in[2],"amplitudes");
        if(amplitudes.empty()||amplitudes.size()>65536)throw std::invalid_argument("Inverse BH amplitude budget exceeded");
        const auto r=radia_hdiv_bh::evaluate(law(in[1]),amplitudes);Builder b;
        b.set("secant",column(r.secant)).set("tangent",column(r.tangent)).set("energy",column(r.energy)).set("owned_device_payload_bytes",number(double(r.ownedDeviceBytes)));
        out[0]=b.release();return;
    }
    if(command=="material_stats"){
        arity(nr,2,nl,1);const auto s=entry(in[1],materialTag,&Entry::material).stats();Builder b;
        b.set("dofs",number(double(s.dofs))).set("points",number(double(s.points))).set("device_payload_bound_bytes",number(double(s.devicePayloadBound)))
         .set("setup_upload_bytes",number(double(s.setupUploadBytes))).set("state_upload_bytes",number(double(s.stateUploadBytes)))
         .set("evaluations",number(double(s.evaluations))).set("tangent_applications",number(double(s.tangentApplications)));
        out[0]=b.release();return;
    }
    throw std::invalid_argument("Unknown HDiv device command");
}
void withdrawn(int nl,mxArray** out,int nr,const mxArray** in){
    for(int i=0;i<nl;++i)out[i]=nullptr;
    try{dispatch(nl,out,nr,in);}catch(...){for(int i=0;i<nl;++i)if(out[i]){mxDestroyArray(out[i]);out[i]=nullptr;}throw;}
}
}
void mexFunction(int nl,mxArray* out[],int nr,const mxArray* in[]){
    char message[2048]={0};
    if(!radia_mex::checked(&withdrawn,nl,out,nr,in,message,sizeof(message)))mexErrMsgIdAndTxt("radia:hdiv:device","%s",message);
}
