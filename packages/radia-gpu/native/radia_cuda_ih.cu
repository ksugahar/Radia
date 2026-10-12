// Experimental persistent IH operator runtime. Physical assembly and dynamic
// ESIM are separate required frontends; unsupported modes fail explicitly.
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include "ih_coupled_state.cuh"
#include "ih_esim_coupling.cuh"
#include "rad_bem_galerkin.h"
#include "rad_biot_savart_filaments.h"
#include <map>
#include <string>
#include <cstdint>
#include <chrono>
using namespace radia_ih_cuda;
namespace {
const mxArray* field(const mxArray* a,const char* name){
 if(!a||!mxIsStruct(a)||mxGetNumberOfElements(a)!=1)throw std::invalid_argument("Expected scalar IH config");
 const auto* value=mxGetField(a,0,name);if(!value)throw std::invalid_argument(std::string("Missing IH field: ")+name);return value;
}
double number(const mxArray* a){if(!a||!mxIsDouble(a)||mxIsComplex(a)||mxIsSparse(a)||mxGetNumberOfElements(a)!=1||!std::isfinite(mxGetScalar(a)))throw std::invalid_argument("Expected finite scalar double");return mxGetScalar(a);}
std::string text(const mxArray* a){if(!a||!mxIsChar(a))throw std::invalid_argument("Expected character vector");char* c=mxArrayToString(a);if(!c)throw std::runtime_error("String conversion failed");std::string s(c);mxFree(c);return s;}
size_t dimension(size_t n){if(!n||n>INT_MAX)throw std::invalid_argument("IH dimension exceeds int32 indexing");return n;}
struct Buffer {
 void* p=nullptr;size_t bytes=0;
 ~Buffer(){if(p)cudaFree(p);}
 template<class T>void put(const std::vector<T>& values){bytes=values.size()*sizeof(T);if(bytes){cudaCheck(cudaMalloc(&p,bytes));cudaCheck(cudaMemcpy(p,values.data(),bytes,cudaMemcpyHostToDevice));}}
 template<class T>const T* get()const{return static_cast<const T*>(p);}
};
struct Matrix {
 size_t rows=0,columns=0;std::vector<int32_t> rowptr,col;
 std::vector<cuDoubleComplex> z;std::vector<double> real;
 Buffer dr,dc,dz,dreal,df;
 void read(const mxArray* a,bool complex){
  if(!mxIsDouble(a)||mxGetNumberOfDimensions(a)!=2||(!complex&&mxIsComplex(a)))throw std::invalid_argument("IH matrix must be double with declared complexity");
  rows=dimension(mxGetM(a));columns=dimension(mxGetN(a));
  rowptr.assign(rows+1,0);
  const bool sparse=mxIsSparse(a);const size_t count=sparse?mxGetJc(a)[columns]:mxGetNumberOfElements(a);
  if(count>INT_MAX)throw std::overflow_error("IH CSR nnz exceeds int32 indexing");
  if(sparse){for(size_t k=0;k<count;++k){if(mxGetIr(a)[k]>=rows)throw std::invalid_argument("Sparse row index out of range");++rowptr[mxGetIr(a)[k]+1];}}
  else for(size_t i=1;i<=rows;++i)rowptr[i]=int32_t(columns);
  for(size_t i=1;i<=rows;++i)rowptr[i]+=rowptr[i-1];
  col.resize(count);if(complex)z.resize(count);else real.resize(count);
  auto cursor=rowptr;
  const auto* re=mxIsComplex(a)?nullptr:mxGetDoubles(a);const auto* zi=mxIsComplex(a)?mxGetComplexDoubles(a):nullptr;
  for(size_t j=0;j<columns;++j){size_t begin=sparse?mxGetJc(a)[j]:j*rows,end=sparse?mxGetJc(a)[j+1]:(j+1)*rows;
   for(size_t k=begin;k<end;++k){const size_t i=sparse?mxGetIr(a)[k]:k-j*rows;const auto pos=cursor[i]++;col[pos]=int32_t(j);
    double r=zi?zi[k].real:re[k],im=zi?zi[k].imag:0.;if(!std::isfinite(r)||!std::isfinite(im))throw std::invalid_argument("Nonfinite IH matrix entry");
    if(complex)z[pos]=make_cuDoubleComplex(r,im);else real[pos]=r;
   }}
 }
 size_t bytes()const{return rowptr.size()*sizeof(int32_t)+col.size()*sizeof(int32_t)+z.size()*sizeof(cuDoubleComplex)+real.size()*sizeof(double);}
 void upload(){dr.put(rowptr);dc.put(col);dz.put(z);dreal.put(real);}
 void uploadFP32(){
  std::vector<cuFloatComplex> values;values.reserve(z.size());
  for(auto x:z){float r=float(cuCreal(x)),im=float(cuCimag(x));
   if(!std::isfinite(r)||!std::isfinite(im)||(cuCreal(x)!=0.&&r==0.f)||(cuCimag(x)!=0.&&im==0.f))
    throw std::invalid_argument("Preconditioner entry is outside FP32 representable range");
   values.push_back(make_cuFloatComplex(r,im));}
  dr.put(rowptr);dc.put(col);df.put(values);
 }
 ComplexCSRView complexView()const{return {rows,dr.get<int32_t>(),dc.get<int32_t>(),dz.get<cuDoubleComplex>()};}
 RealCSRView realView()const{return {rows,dr.get<int32_t>(),dc.get<int32_t>(),dreal.get<double>()};}
 void abandon(){dr.p=dc.p=dz.p=dreal.p=df.p=nullptr;}
};
std::vector<double> realVector(const mxArray* a,size_t count){
 if(!mxIsDouble(a)||mxIsComplex(a)||mxIsSparse(a)||mxGetNumberOfElements(a)!=count)throw std::invalid_argument("IH real vector shape mismatch");
 const auto* p=mxGetDoubles(a);std::vector<double> v(p,p+count);for(double x:v)if(!std::isfinite(x))throw std::invalid_argument("Nonfinite IH vector");return v;
}
#include "ih_peec_bem_assembly.cuh"
// Explicit bounded operator frontend. Geometry/reaction assembly is separate.
struct ElementMaterial {
 Matrix evaluation,entryMap;Buffer axis,rho,mu,saturation,geometry,fields,base,baseWeights;
 std::unique_ptr<ESIMCoupling> coupling;
 std::size_t owned=0;
 void abandon(){if(coupling)coupling->abandonAfterDeviceReset();evaluation.abandon();entryMap.abandon();for(auto* b:{&axis,&rho,&mu,&saturation,&geometry,&fields,&base,&baseWeights})b->p=nullptr;}
 explicit ElementMaterial(const mxArray* cfg,Matrix& a,Matrix& p,Buffer& weights,const std::vector<double>& w,std::size_t thermalDOF,unsigned components){
  const auto elements=p.rows/components;if(elements>4096||a.rows>4096)throw std::invalid_argument("Element ESIM research dimensions exceed 4096");
  evaluation.read(field(cfg,"temperature_evaluation"),false);entryMap.read(field(cfg,"operator_impedance_map"),true);
  if(evaluation.rows!=elements||evaluation.columns!=thermalDOF||entryMap.rows!=a.z.size()||entryMap.columns!=elements)
   throw std::invalid_argument("Element ESIM map dimensions disagree; impedance map rows follow eddy CSR entries");
  const auto* kelvin=field(cfg,"temperature_axis_K");const auto nt=mxGetNumberOfElements(kelvin);
  if(nt<2||nt>4096)throw std::invalid_argument("Temperature axis requires 2..4096 entries");
  auto temperatures=realVector(kelvin,nt);for(size_t j=0;j<nt;++j)if(temperatures[j]<=0.||(j&&temperatures[j]<=temperatures[j-1]))throw std::invalid_argument("Kelvin axis must be positive and strictly increasing");
  auto properties=[&](const char* name,bool positive){const auto* input=field(cfg,name);if(mxGetM(input)!=elements||mxGetN(input)!=nt)throw std::invalid_argument("Material table must be element by temperature");auto source=realVector(input,elements*nt);std::vector<double> out(elements*nt);for(size_t e=0;e<elements;++e)for(size_t t=0;t<nt;++t){auto v=source[e+elements*t];if(positive?v<=0.:v<0.)throw std::invalid_argument("Invalid material table value");out[e*nt+t]=v;}return out;};
  auto resistivity=properties("resistivity_ohm_m",true),sat=properties("saturation_magnetization_A_per_m",false);
  const auto* inputMu=field(cfg,"initial_permeability_H_per_m");
  if(!mxIsDouble(inputMu)||mxIsSparse(inputMu)||mxGetM(inputMu)!=elements||mxGetN(inputMu)!=nt)throw std::invalid_argument("Permeability table shape/type mismatch");
  std::vector<cuDoubleComplex> permeability(elements*nt);for(size_t e=0;e<elements;++e)for(size_t t=0;t<nt;++t){auto j=e+elements*t;auto z=mxIsComplex(inputMu)?make_cuDoubleComplex(mxGetComplexDoubles(inputMu)[j].real,mxGetComplexDoubles(inputMu)[j].imag):make_cuDoubleComplex(mxGetDoubles(inputMu)[j],0.);if(!cellFiniteHost(z)||cuCreal(z)<=0.||cuCimag(z)>0.)throw std::invalid_argument("Permeability must be finite and passive");permeability[e*nt+t]=z;}
  if(mxGetField(cfg,0,"bh_knots")||mxGetField(cfg,0,"bh_coefficients"))throw std::invalid_argument("Temperature-dependent nonlinear B-H gateway is not yet accepted");
  const auto* inputGeometry=field(cfg,"esim_cell_parameters");if(mxGetM(inputGeometry)!=elements||mxGetN(inputGeometry)!=7)throw std::invalid_argument("ESIM cells require element by 7 parameters");auto values=realVector(inputGeometry,elements*7);
  std::vector<ESIMCellParameters> cells(elements);int maxNodes=0;for(size_t e=0;e<elements;++e){auto at=[&](size_t k){return values[e+elements*k];};auto& cell=cells[e];
   if(at(0)<3||at(0)>4096||at(0)!=int(at(0))||at(1)<1||at(1)>10000||at(1)!=int(at(1))||(at(2)!=0.&&at(2)!=1.)||at(3)<=0.||at(4)<=0.||at(5)<=0.||at(5)>1e-6||at(6)<=0.||at(6)>1.)throw std::invalid_argument("Invalid ESIM cell geometry/iteration contract");
   cell.nodes=int(at(0));cell.maxIterations=int(at(1));cell.cylinder=at(2)!=0.;cell.halfThickness=at(3);cell.omega=at(4);cell.tolerance=at(5);cell.relaxation=at(6);maxNodes=std::max(maxNodes,cell.nodes);}
  auto guesses=realVector(field(cfg,"initial_field_A_per_m"),elements);for(auto x:guesses)if(x<=0.)throw std::invalid_argument("ESIM initial fields must be positive");
  double iterations=number(field(cfg,"max_constitutive_iterations"));if(iterations<2||iterations>10000||iterations!=int(iterations))throw std::invalid_argument("Invalid constitutive iteration limit");
  owned=evaluation.bytes()+entryMap.bytes()+temperatures.size()*sizeof(double)+resistivity.size()*sizeof(double)+permeability.size()*sizeof(cuDoubleComplex)+sat.size()*sizeof(double)+cells.size()*sizeof(ESIMCellParameters)+guesses.size()*sizeof(double)+a.z.size()*sizeof(cuDoubleComplex)+w.size()*sizeof(double);
  size_t available,total;cudaCheck(cudaMemGetInfo(&available,&total));auto reserve=std::max<size_t>(2ULL<<30,total/5);if(available<=reserve||owned>available-reserve)throw std::runtime_error("Element material upload exceeds GPU reserve");
  evaluation.upload();entryMap.upload();axis.put(temperatures);rho.put(resistivity);mu.put(permeability);saturation.put(sat);geometry.put(cells);fields.put(guesses);base.put(a.z);baseWeights.put(w);
  ESIMCouplingView v{};v.temperatureEvaluation=evaluation.realView();v.table={axis.get<double>(),rho.get<double>(),saturation.get<double>(),nullptr,nullptr,mu.get<cuDoubleComplex>(),int(nt),0};v.entryMap=entryMap.complexView();v.fieldProjection=p.complexView();v.matrix=a.complexView();v.baseValues=base.get<cuDoubleComplex>();v.mutableValues=static_cast<cuDoubleComplex*>(a.dz.p);v.geometry=geometry.get<ESIMCellParameters>();v.initialFields=fields.get<double>();v.baseLossWeights=baseWeights.get<double>();v.mutableLossWeights=static_cast<double*>(weights.p);v.elements=elements;v.components=components;v.workspaceNodes=maxNodes;v.maxIterations=int(iterations);v.tolerance=number(field(cfg,"constitutive_tolerance"));v.relaxation=number(field(cfg,"constitutive_relaxation"));v.fieldFloor=number(field(cfg,"field_floor_A_per_m"));coupling=std::make_unique<ESIMCoupling>(v);owned+=coupling->ownedBytes();
 }
 static bool cellFiniteHost(cuDoubleComplex z){return std::isfinite(cuCreal(z))&&std::isfinite(cuCimag(z));}
};
struct Context {
 int device=0;cudaStream_t guard=nullptr;
 Matrix a,p,system,k,c,projection,preconditioner;
 PreconditionerPrecision precision=PreconditionerPrecision::None;Buffer rhs,weights,ones,initial;
 std::unique_ptr<ElementMaterial> material;
 std::unique_ptr<CoupledState> state;size_t owned=0;
 ~Context(){
  int previous=device;cudaGetDevice(&previous);if(previous!=device)cudaSetDevice(device);
  const auto status=guard?cudaStreamQuery(guard):cudaSuccess;
  if(status!=cudaSuccess&&status!=cudaErrorNotReady){
   if(state)state->abandonAfterDeviceReset();if(material)material->abandon();
   for(auto* m:{&a,&p,&system,&k,&c,&projection,&preconditioner})m->abandon();
   rhs.p=weights.p=ones.p=initial.p=nullptr;guard=nullptr;
  }
  state.reset();material.reset();
  // Free operators while their owning device is selected.
  for(auto* m:{&a,&p,&system,&k,&c,&projection,&preconditioner})for(auto* b:{&m->dr,&m->dc,&m->dz,&m->dreal,&m->df}){if(b->p)cudaFree(b->p);b->p=nullptr;}
  for(auto* b:{&rhs,&weights,&ones,&initial}){if(b->p)cudaFree(b->p);b->p=nullptr;}
  if(guard)cudaStreamDestroy(guard);if(previous!=device)cudaSetDevice(previous);
 }
 void check(){int current;cudaCheck(cudaGetDevice(&current));if(current!=device)throw std::runtime_error("IH context belongs to another GPU");cudaCheck(cudaStreamQuery(guard));}
 explicit Context(const mxArray* cfg){
  if(number(field(cfg,"schema_version"))!=1.)throw std::invalid_argument("Unsupported IH operator schema");
  const auto materialMode=text(field(cfg,"material_mode"));if(materialMode!="fixed"&&materialMode!="element-esim-operators")throw std::invalid_argument("Unsupported IH material frontend");
  const auto* precisionField=mxGetField(cfg,0,"preconditioner_precision");
  if(precisionField){const auto mode=text(precisionField);
   if(mode=="fp32")precision=PreconditionerPrecision::FP32;
   else if(mode=="fp64")precision=PreconditionerPrecision::FP64;
   else if(mode!="none")throw std::invalid_argument("Unsupported IH preconditioner precision");
  }
  a.read(field(cfg,"eddy"),true);p.read(field(cfg,"loss_projection"),true);system.read(field(cfg,"thermal_system"),true);
  k.read(field(cfg,"stiffness"),false);c.read(field(cfg,"convection"),false);projection.read(field(cfg,"heat_projection"),false);
  if(precision!=PreconditionerPrecision::None){preconditioner.read(field(cfg,"preconditioner_inverse"),true);
   if(preconditioner.rows!=a.rows||preconditioner.columns!=a.columns)throw std::invalid_argument("Fixed preconditioner dimensions disagree");}
  const auto n=system.rows;double components=number(field(cfg,"loss_components"));
  if(components<1||components>INT_MAX||components!=int(components)||p.rows%int(components))throw std::invalid_argument("Invalid IH loss components");
  const size_t samples=p.rows/size_t(components);
  if(a.rows!=a.columns||p.columns!=a.rows||system.columns!=n||k.rows!=n||k.columns!=n||c.rows!=n||c.columns!=n||projection.rows!=n||projection.columns!=samples)throw std::invalid_argument("IH coupled matrix dimensions disagree");
  auto w=realVector(field(cfg,"loss_weights"),p.rows),o=realVector(field(cfg,"constant_coefficients"),n),t=realVector(field(cfg,"initial_temperature_K"),n);
  for(double x:w)if(x<0.)throw std::invalid_argument("Negative IH loss weight");
  const auto* r=field(cfg,"rhs_per_amp");if(mxIsSparse(r)||!mxIsDouble(r)||mxGetNumberOfElements(r)!=a.rows)throw std::invalid_argument("IH RHS shape mismatch");
  std::vector<cuDoubleComplex> b(a.rows);for(size_t i=0;i<a.rows;++i){double re=mxIsComplex(r)?mxGetComplexDoubles(r)[i].real:mxGetDoubles(r)[i],im=mxIsComplex(r)?mxGetComplexDoubles(r)[i].imag:0.;if(!std::isfinite(re)||!std::isfinite(im))throw std::invalid_argument("Nonfinite RHS");b[i]=make_cuDoubleComplex(re,im);}
  double dt=number(field(cfg,"sample_time_s")),conductivity=number(field(cfg,"conductivity_scale")),convection=number(field(cfg,"convection_scale"));
  if(dt<=0.||conductivity<0.||convection<0.)throw std::invalid_argument("Invalid IH thermal parameters");
  owned=a.bytes()+p.bytes()+system.bytes()+k.bytes()+c.bytes()+projection.bytes()+b.size()*sizeof(cuDoubleComplex)+(w.size()+o.size()+t.size())*sizeof(double);
  if(precision!=PreconditionerPrecision::None)owned+=preconditioner.rowptr.size()*sizeof(int32_t)+preconditioner.col.size()*sizeof(int32_t)+preconditioner.z.size()*(precision==PreconditionerPrecision::FP32?sizeof(cuFloatComplex):sizeof(cuDoubleComplex));
  const size_t workspace=(a.rows*(precision==PreconditionerPrecision::None?10:12)+n*9+2)*sizeof(cuDoubleComplex)+(2*samples+3*n+2)*sizeof(double)+2*sizeof(BiCGState)+sizeof(int);
  size_t available,total;cudaCheck(cudaMemGetInfo(&available,&total));const size_t reserve=std::max<size_t>(2ULL<<30,total/5);
  if(available<=reserve||owned>available-reserve||workspace>available-reserve-owned)throw std::runtime_error("IH context exceeds GPU memory budget");
  cudaCheck(cudaGetDevice(&device));cudaCheck(cudaStreamCreateWithFlags(&guard,cudaStreamNonBlocking));
  try {
  for(auto* m:{&a,&p,&system,&k,&c,&projection})m->upload();
  if(precision==PreconditionerPrecision::FP32)preconditioner.uploadFP32();
  else if(precision==PreconditionerPrecision::FP64)preconditioner.upload();
  rhs.put(b);weights.put(w);ones.put(o);initial.put(t);
  CoupledOperators op{};op.eddy=a.complexView();op.lossProjection=p.complexView();op.thermalSystem=system.complexView();op.stiffness=k.realView();op.convection=c.realView();op.heatProjection=projection.realView();
  if(precision!=PreconditionerPrecision::None)op.eddyPreconditioner={a.rows,preconditioner.dr.get<int32_t>(),preconditioner.dc.get<int32_t>(),preconditioner.dz.get<cuDoubleComplex>(),preconditioner.df.get<cuFloatComplex>(),precision};
  op.rhsPerAmp=rhs.get<cuDoubleComplex>();op.lossWeights=weights.get<double>();op.constantCoefficients=ones.get<double>();op.initialTemperature=initial.get<double>();
  op.heatSamples=samples;op.lossComponents=unsigned(components);op.dt=dt;op.conductivity=conductivity;op.convectionScale=convection;
  if(materialMode=="element-esim-operators"){
   material=std::make_unique<ElementMaterial>(cfg,a,p,weights,w,n,unsigned(components));owned+=material->owned;op.dynamicMaterial=true;op.materialOwner=material->coupling.get();
   op.updateMaterial=[](void* owner,cudaStream_t stream,const double* t,double current){static_cast<ESIMCoupling*>(owner)->prepare(stream,t,current);};
   op.solveMaterial=[](void* owner,ComplexSolver& solver,double tolerance,int iterations){return static_cast<ESIMCoupling*>(owner)->solve(solver,tolerance,iterations);};
   op.commitMaterial=[](void* owner){static_cast<ESIMCoupling*>(owner)->commit();};op.resetMaterial=[](void* owner){static_cast<ESIMCoupling*>(owner)->reset();};
  }
  state=std::make_unique<CoupledState>(op);owned+=state->ownedBytes();
  }catch(...){cudaStreamDestroy(guard);guard=nullptr;throw;}
 }
};
std::map<uint64_t,std::unique_ptr<Context>> contexts;uint64_t generation=uint64_t(std::chrono::steady_clock::now().time_since_epoch().count())&0x00007fffffffffffULL;
void cleanup(){contexts.clear();}
uint64_t handle(const mxArray* a){if(!mxIsUint64(a)||mxIsComplex(a)||mxGetNumberOfElements(a)!=1)throw std::invalid_argument("IH handle must be scalar uint64");return *static_cast<const uint64_t*>(mxGetData(a));}
Context& get(const mxArray* a){auto h=handle(a);auto it=contexts.find(h);if(it==contexts.end())throw std::invalid_argument("IH handle has wrong type, generation or liveness");it->second->check();return *it->second;}
mxArray* vectorOut(const std::vector<double>& v){auto* a=mxCreateDoubleMatrix(v.size(),1,mxREAL);std::copy(v.begin(),v.end(),mxGetDoubles(a));return a;}
}
void mexFunction(int nlhs,mxArray** out,int nrhs,const mxArray** in){
 try{
  mxInitGPU();if(!nrhs)throw std::invalid_argument("Missing IH command");const auto command=text(in[0]);
  if(command=="assemble_peec_bem"){if(nrhs!=2||nlhs!=1)throw std::invalid_argument("assemble_peec_bem requires geometry and one output");out[0]=assemblePEECBEM(in[1]);return;}
  if(command=="create"){
   if(nrhs!=2||nlhs!=1)throw std::invalid_argument("create requires config and one output");
   if(generation==0x0000ffffffffffffULL)throw std::overflow_error("IH handle generation exhausted");
   auto context=std::make_unique<Context>(in[1]);const uint64_t h=0x4948000000000000ULL|++generation;
   out[0]=mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);*static_cast<uint64_t*>(mxGetData(out[0]))=h;
   contexts.emplace(h,std::move(context));mexLock();mexAtExit(cleanup);return;
  }
  if(command=="live_count"){if(nrhs!=1||nlhs!=1)throw std::invalid_argument("live_count takes no arguments");out[0]=mxCreateDoubleScalar(double(contexts.size()));return;}
  if(command=="destroy"){
   if(nrhs!=2||nlhs)throw std::invalid_argument("destroy requires handle and no output");const auto h=handle(in[1]);
   if(!contexts.erase(h))throw std::invalid_argument("IH stale or wrong-type handle");mexUnlock();return;
  }
  if(nrhs<2)throw std::invalid_argument("Missing IH handle");auto& c=get(in[1]);
  if(command=="reset"){if(nrhs!=2||nlhs)throw std::invalid_argument("reset signature mismatch");c.state->reset();return;}
  if(command=="diagnostic_operator"){
   if(nrhs!=2||nlhs!=1)throw std::invalid_argument("diagnostic_operator returns candidate operator, not accepted state");
   std::vector<cuDoubleComplex> values(c.a.z.size());
   cudaCheck(cudaDeviceSynchronize());cudaCheck(cudaMemcpy(values.data(),c.a.dz.p,values.size()*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost));
   out[0]=mxCreateDoubleMatrix(c.a.rows,c.a.columns,mxCOMPLEX);
   for(size_t row=0;row<c.a.rows;++row)for(int32_t k=c.a.rowptr[row];k<c.a.rowptr[row+1];++k){auto& z=mxGetComplexDoubles(out[0])[row+c.a.rows*c.a.col[k]];z.real=cuCreal(values[k]);z.imag=cuCimag(values[k]);}return;
  }
  if(command=="update"){if(nrhs!=4||nlhs)throw std::invalid_argument("update requires current and ambient");c.state->update(number(in[2]),number(in[3]));return;}
  if(command=="electromagnetic_state"){
   if(nrhs!=2||nlhs!=1)throw std::invalid_argument("electromagnetic_state returns one complex vector");
   const auto v=c.state->electromagneticState();out[0]=mxCreateDoubleMatrix(v.size(),1,mxCOMPLEX);
   for(size_t i=0;i<v.size();++i){mxGetComplexDoubles(out[0])[i].real=cuCreal(v[i]);mxGetComplexDoubles(out[0])[i].imag=cuCimag(v[i]);}return;
  }
  if(command=="outputs"){
   if(nrhs!=2||nlhs!=3)throw std::invalid_argument("outputs returns temperature, heat and info");
   out[0]=vectorOut(c.state->temperature());out[1]=vectorOut(c.state->heat());const char* f[]={"steps","owned_device_bytes","eddy_iterations","eddy_true_residual","thermal_iterations","thermal_true_residual","preconditioner_precision","constitutive_iterations","constitutive_true_residual","constitutive_field_residual","material_mode"};out[2]=mxCreateStructMatrix(1,1,11,f);
   mxSetField(out[2],0,f[0],mxCreateDoubleScalar(double(c.state->steps())));mxSetField(out[2],0,f[1],mxCreateDoubleScalar(double(c.owned)));
   mxSetField(out[2],0,f[2],mxCreateDoubleScalar(c.state->eddyResult().iterations));
   mxSetField(out[2],0,f[3],mxCreateDoubleScalar(c.state->eddyResult().trueRelativeResidual));
   mxSetField(out[2],0,f[4],mxCreateDoubleScalar(c.state->thermalResult().iterations));
   mxSetField(out[2],0,f[5],mxCreateDoubleScalar(c.state->thermalResult().trueRelativeResidual));
   mxSetField(out[2],0,f[6],mxCreateString(c.precision==PreconditionerPrecision::FP32?"fp32":c.precision==PreconditionerPrecision::FP64?"fp64":"none"));mxSetField(out[2],0,f[7],mxCreateDoubleScalar(c.material?c.material->coupling->iterations():0));mxSetField(out[2],0,f[8],mxCreateDoubleScalar(c.material?c.material->coupling->residual():0));mxSetField(out[2],0,f[9],mxCreateDoubleScalar(c.material?c.material->coupling->fieldResidual():0));mxSetField(out[2],0,f[10],mxCreateString(c.material?"element-esim-operators":"fixed"));return;
  }
  throw std::invalid_argument("Unknown IH CUDA command");
 }catch(const std::exception& e){mexErrMsgIdAndTxt("radia:ih:cuda","%s",e.what());}
}
