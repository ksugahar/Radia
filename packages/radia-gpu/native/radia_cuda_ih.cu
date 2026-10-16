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
 Matrix evaluation,entryMap;Buffer axis,rho,mu,saturation,geometry,fields,base,baseWeights,bhKnots,bhCoefficients;
 std::unique_ptr<ESIMCoupling> coupling;
 std::size_t owned=0;
 void abandon(){if(coupling)coupling->abandonAfterDeviceReset();evaluation.abandon();entryMap.abandon();for(auto* b:{&axis,&rho,&mu,&saturation,&geometry,&fields,&base,&baseWeights,&bhKnots,&bhCoefficients})b->p=nullptr;}
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
  // Measured anhysteretic amplitude law, piecewise linear in H and T.
  // Build cubics locally (zero cubic/quadratic terms), never through Python.
  std::vector<double> knots,coefficients;
  const auto* inputKnots=mxGetField(cfg,0,"bh_knots");
  const auto* inputB=mxGetField(cfg,0,"bh_values_T");
  if(mxGetField(cfg,0,"bh_coefficients"))throw std::invalid_argument("Use checked bh_values_T, not unchecked polynomial coefficients");
  if(bool(inputKnots)!=bool(inputB))throw std::invalid_argument("B-H requires both bh_knots and bh_values_T");
  if(inputKnots){
   const auto nk=mxGetNumberOfElements(inputKnots);
   if(nk<2||nk>4096)throw std::invalid_argument("B-H requires 2..4096 common H knots");
   if(elements*nt*(nk-1)*4>(256ULL<<20)/sizeof(double))throw std::overflow_error("B-H coefficient table exceeds 256 MiB host budget");
   knots=realVector(inputKnots,nk);
   if(knots[0]!=0.)throw std::invalid_argument("B-H H axis must start at zero");
   for(size_t k=1;k<nk;++k)if(knots[k]<=knots[k-1])throw std::invalid_argument("B-H H axis must strictly increase");
   const auto* dims=mxGetDimensions(inputB);
   if(mxGetNumberOfDimensions(inputB)!=3||dims[0]!=elements||dims[1]!=nt||dims[2]!=nk)throw std::invalid_argument("bh_values_T must be element by temperature by H knot");
   auto values=realVector(inputB,elements*nt*nk);
   coefficients.resize(elements*nt*(nk-1)*4);
   constexpr double mu0=1.2566370614359172954e-6;
   for(size_t e=0;e<elements;++e)for(size_t t=0;t<nt;++t){
    auto at=[&](size_t k){return values[e+elements*(t+nt*k)];};
    if(at(0)!=0.||cuCimag(permeability[e*nt+t])!=0.)throw std::invalid_argument("Amplitude B-H starts at B=0 and requires real initial permeability; hysteresis is not represented");
    for(size_t k=0;k+1<nk;++k){
     const double slope=(at(k+1)-at(k))/(knots[k+1]-knots[k]);
     if(!std::isfinite(slope)||slope<mu0)throw std::invalid_argument("B-H differential permeability must be finite and at least mu0");
     const auto j=((e*nt+t)*(nk-1)+k)*4;
     coefficients[j]=coefficients[j+1]=0.;coefficients[j+2]=slope;coefficients[j+3]=at(k);
     if(k==0&&std::abs(slope-cuCreal(permeability[e*nt+t]))>1e-10*std::max(slope,cuCreal(permeability[e*nt+t])))throw std::invalid_argument("Initial permeability must match first B-H segment");
    }
    const double tail=mu0*(knots.back()+sat[e*nt+t]);
    if(std::abs(tail-at(nk-1))>1e-10*std::max(tail,at(nk-1)))throw std::invalid_argument("B-H saturation tail must be continuous");
   }
  }
  const auto* inputGeometry=field(cfg,"esim_cell_parameters");if(mxGetM(inputGeometry)!=elements||mxGetN(inputGeometry)!=7)throw std::invalid_argument("ESIM cells require element by 7 parameters");auto values=realVector(inputGeometry,elements*7);
  std::vector<ESIMCellParameters> cells(elements);int maxNodes=0;for(size_t e=0;e<elements;++e){auto at=[&](size_t k){return values[e+elements*k];};auto& cell=cells[e];
   if(at(0)<3||at(0)>4096||at(0)!=int(at(0))||at(1)<1||at(1)>10000||at(1)!=int(at(1))||(at(2)!=0.&&at(2)!=1.)||at(3)<=0.||at(4)<=0.||at(5)<=0.||at(5)>1e-6||at(6)<=0.||at(6)>1.)throw std::invalid_argument("Invalid ESIM cell geometry/iteration contract");
   cell.nodes=int(at(0));cell.maxIterations=int(at(1));cell.cylinder=at(2)!=0.;cell.halfThickness=at(3);cell.omega=at(4);cell.tolerance=at(5);cell.relaxation=at(6);maxNodes=std::max(maxNodes,cell.nodes);}
  auto guesses=realVector(field(cfg,"initial_field_A_per_m"),elements);for(auto x:guesses)if(x<=0.)throw std::invalid_argument("ESIM initial fields must be positive");
  double iterations=number(field(cfg,"max_constitutive_iterations"));if(iterations<2||iterations>10000||iterations!=int(iterations))throw std::invalid_argument("Invalid constitutive iteration limit");
  owned=evaluation.bytes()+entryMap.bytes()+temperatures.size()*sizeof(double)+resistivity.size()*sizeof(double)+permeability.size()*sizeof(cuDoubleComplex)+sat.size()*sizeof(double)+cells.size()*sizeof(ESIMCellParameters)+guesses.size()*sizeof(double)+a.z.size()*sizeof(cuDoubleComplex)+w.size()*sizeof(double);
  owned+=(knots.size()+coefficients.size())*sizeof(double);
  size_t available,total;cudaCheck(cudaMemGetInfo(&available,&total));auto reserve=std::max<size_t>(2ULL<<30,total/5);if(available<=reserve||owned>available-reserve)throw std::runtime_error("Element material upload exceeds GPU reserve");
  evaluation.upload();entryMap.upload();axis.put(temperatures);rho.put(resistivity);mu.put(permeability);saturation.put(sat);geometry.put(cells);fields.put(guesses);base.put(a.z);baseWeights.put(w);
  bhKnots.put(knots);bhCoefficients.put(coefficients);
  ESIMCouplingView v{};v.temperatureEvaluation=evaluation.realView();v.table={axis.get<double>(),rho.get<double>(),saturation.get<double>(),bhKnots.get<double>(),bhCoefficients.get<double>(),mu.get<cuDoubleComplex>(),int(nt),int(knots.size())};v.entryMap=entryMap.complexView();v.fieldProjection=p.complexView();v.matrix=a.complexView();v.baseValues=base.get<cuDoubleComplex>();v.mutableValues=static_cast<cuDoubleComplex*>(a.dz.p);v.geometry=geometry.get<ESIMCellParameters>();v.initialFields=fields.get<double>();v.baseLossWeights=baseWeights.get<double>();v.mutableLossWeights=static_cast<double*>(weights.p);v.elements=elements;v.components=components;v.workspaceNodes=maxNodes;v.maxIterations=int(iterations);v.tolerance=number(field(cfg,"constitutive_tolerance"));v.relaxation=number(field(cfg,"constitutive_relaxation"));v.fieldFloor=number(field(cfg,"field_floor_A_per_m"));coupling=std::make_unique<ESIMCoupling>(v);owned+=coupling->ownedBytes();
 }
 static bool cellFiniteHost(cuDoubleComplex z){return std::isfinite(cuCreal(z))&&std::isfinite(cuCimag(z));}
};
struct Context {
 int device=0;cudaStream_t guard=nullptr;bool poisoned=false;std::vector<int32_t> snapshotIndices,snapshotBlockDims;std::vector<std::string> snapshotUnits;cuDoubleComplex suppliedCoilImpedance{};
 double solveTolerance=1e-10;
 Matrix a,p,system,k,c,projection,preconditioner;
 PreconditionerPrecision precision=PreconditionerPrecision::None;Buffer rhs,weights,ones,initial;
 std::unique_ptr<ElementMaterial> material;
 std::unique_ptr<WorkSnapshot> snapshot;
 std::unique_ptr<CoupledState> state;size_t owned=0;
 ~Context(){
  int previous=device;cudaGetDevice(&previous);if(previous!=device)cudaSetDevice(device);
  const auto status=guard?cudaStreamQuery(guard):cudaSuccess;
  if(status!=cudaSuccess&&status!=cudaErrorNotReady){
   if(snapshot)snapshot->abandon();if(state)state->abandonAfterDeviceReset();if(material)material->abandon();
   for(auto* m:{&a,&p,&system,&k,&c,&projection,&preconditioner})m->abandon();
   rhs.p=weights.p=ones.p=initial.p=nullptr;guard=nullptr;
  }
  state.reset();snapshot.reset();material.reset();
  // Free operators while their owning device is selected.
  for(auto* m:{&a,&p,&system,&k,&c,&projection,&preconditioner})for(auto* b:{&m->dr,&m->dc,&m->dz,&m->dreal,&m->df}){if(b->p)cudaFree(b->p);b->p=nullptr;}
  for(auto* b:{&rhs,&weights,&ones,&initial}){if(b->p)cudaFree(b->p);b->p=nullptr;}
  if(guard)cudaStreamDestroy(guard);if(previous!=device)cudaSetDevice(previous);
 }
 static bool sticky(cudaError_t status)noexcept{return status==cudaErrorIllegalAddress||status==cudaErrorLaunchFailure||status==cudaErrorAssert||status==cudaErrorECCUncorrectable||status==cudaErrorLaunchTimeout||status==cudaErrorUnknown;}
 void noteFailure()noexcept{const auto last=cudaPeekAtLastError();const auto status=guard?cudaStreamQuery(guard):cudaSuccess;if(sticky(last)||sticky(status))poisoned=true;}
 void recoverFailure()noexcept{noteFailure();if(!poisoned&&state){const auto status=state->drainCandidateWork();if(sticky(status))poisoned=true;}}
 void check(){if(poisoned)throw std::runtime_error("IH context poisoned; destroy only");noteFailure();if(poisoned)throw std::runtime_error("IH context poisoned; destroy only");int current;cudaCheck(cudaGetDevice(&current));if(current!=device)throw std::runtime_error("IH context belongs to another GPU");cudaCheck(cudaStreamQuery(guard));}
 explicit Context(const mxArray* cfg){
  if(number(field(cfg,"schema_version"))!=1.)throw std::invalid_argument("Unsupported IH operator schema");
  const auto materialMode=text(field(cfg,"material_mode"));if(materialMode!="fixed"&&materialMode!="element-esim-operators")throw std::invalid_argument("Unsupported IH material frontend");
  const auto* precisionField=mxGetField(cfg,0,"preconditioner_precision");
  const auto* blocksField=mxGetField(cfg,0,"preconditioner_block_ends");std::vector<int> blockEnds;
  if(precisionField){const auto mode=text(precisionField);
   if(mode=="fp32")precision=PreconditionerPrecision::FP32;
   else if(mode=="fp64")precision=PreconditionerPrecision::FP64;
   else if(mode!="none")throw std::invalid_argument("Unsupported IH preconditioner precision");
  }
  a.read(field(cfg,"eddy"),true);p.read(field(cfg,"loss_projection"),true);system.read(field(cfg,"thermal_system"),true);
  k.read(field(cfg,"stiffness"),false);c.read(field(cfg,"convection"),false);projection.read(field(cfg,"heat_projection"),false);
  if(blocksField){
   if(precision!=PreconditionerPrecision::FP64||a.rows>512||mxGetField(cfg,0,"preconditioner_inverse"))throw std::invalid_argument("GPU block setup requires FP64, <=512 unknowns and no supplied inverse");
   auto ends=realVector(blocksField,mxGetNumberOfElements(blocksField));int previous=0;
   for(double x:ends){if(x<=previous||x>a.rows||x!=int(x))throw std::invalid_argument("Invalid block partition");blockEnds.push_back(int(x));previous=int(x);}if(previous!=a.rows)throw std::invalid_argument("Block partition must cover operator");
   preconditioner.rows=preconditioner.columns=a.rows;preconditioner.rowptr.resize(a.rows+1);preconditioner.col.resize(a.rows*a.rows);preconditioner.z.assign(a.rows*a.rows,make_cuDoubleComplex(0.,0.));
   for(size_t i=0;i<=a.rows;++i)preconditioner.rowptr[i]=int32_t(i*a.rows);for(size_t i=0;i<a.rows;++i)for(size_t j=0;j<a.rows;++j)preconditioner.col[i*a.rows+j]=int32_t(j);
  }else if(precision!=PreconditionerPrecision::None){preconditioner.read(field(cfg,"preconditioner_inverse"),true);
   if(preconditioner.rows!=a.rows||preconditioner.columns!=a.columns)throw std::invalid_argument("Fixed preconditioner dimensions disagree");}
  if(const auto* requested=mxGetField(cfg,0,"solver_relative_tolerance")){
   solveTolerance=number(requested);if(solveTolerance<=0.||solveTolerance>1e-6)throw std::invalid_argument("Invalid explicit IH solver tolerance");
  }
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
  const size_t workspace=(a.rows*(precision==PreconditionerPrecision::None?10:12)+n*9+2)*sizeof(cuDoubleComplex)+(2*samples+3*n+2)*sizeof(double)+2*sizeof(BiCGState)+sizeof(int)+(blocksField?a.rows*a.rows*sizeof(cuDoubleComplex)+blockEnds.size()*sizeof(int)+sizeof(int):0);
  size_t available,total;cudaCheck(cudaMemGetInfo(&available,&total));const size_t reserve=std::max<size_t>(2ULL<<30,total/5);
  if(available<=reserve||owned>available-reserve||workspace>available-reserve-owned)throw std::runtime_error("IH context exceeds GPU memory budget");
  cudaCheck(cudaGetDevice(&device));cudaCheck(cudaStreamCreateWithFlags(&guard,cudaStreamNonBlocking));
  try {
  for(auto* m:{&a,&p,&system,&k,&c,&projection})m->upload();
  if(precision==PreconditionerPrecision::FP32)preconditioner.uploadFP32();
  else if(precision==PreconditionerPrecision::FP64)preconditioner.upload();
  if(blocksField){Buffer work,ends,status;work.put(std::vector<cuDoubleComplex>(a.rows*a.rows));ends.put(blockEnds);status.put(std::vector<int>{0});
   buildBlockInverse<<<unsigned(blockEnds.size()),1>>>(a.complexView(),ends.get<int>(),int(blockEnds.size()),static_cast<cuDoubleComplex*>(work.p),static_cast<cuDoubleComplex*>(preconditioner.dz.p),static_cast<int*>(status.p));
   cudaCheck(cudaGetLastError());int failed=0;cudaCheck(cudaMemcpy(&failed,status.p,sizeof(int),cudaMemcpyDeviceToHost));if(failed)throw std::runtime_error("GPU block preconditioner is singular or nonfinite");
  }
  rhs.put(b);weights.put(w);ones.put(o);initial.put(t);
  CoupledOperators op{};op.eddy=a.complexView();op.lossProjection=p.complexView();op.thermalSystem=system.complexView();op.stiffness=k.realView();op.convection=c.realView();op.heatProjection=projection.realView();
  if(precision!=PreconditionerPrecision::None)op.eddyPreconditioner={a.rows,preconditioner.dr.get<int32_t>(),preconditioner.dc.get<int32_t>(),preconditioner.dz.get<cuDoubleComplex>(),preconditioner.df.get<cuFloatComplex>(),precision};
  op.rhsPerAmp=rhs.get<cuDoubleComplex>();op.lossWeights=weights.get<double>();op.constantCoefficients=ones.get<double>();op.initialTemperature=initial.get<double>();
  op.heatSamples=samples;op.lossComponents=unsigned(components);op.dt=dt;op.conductivity=conductivity;op.convectionScale=convection;
  if(materialMode=="element-esim-operators"){
   material=std::make_unique<ElementMaterial>(cfg,a,p,weights,w,n,unsigned(components));owned+=material->owned;op.dynamicMaterial=true;op.materialOwner=material->coupling.get();
   op.updateMaterial=[](void* owner,cudaStream_t stream,const double* t,double current){static_cast<ESIMCoupling*>(owner)->prepare(stream,t,current);};
   op.solveMaterial=[](void* owner,ComplexSolver& solver,double tolerance,int iterations){return static_cast<ESIMCoupling*>(owner)->solve(solver,tolerance,iterations);};
   static_assert(noexcept(std::declval<ESIMCoupling&>().commit()),"ESIM commit must be noexcept");
   op.stageMaterial=[](void* owner,cudaStream_t stream,bool reset){static_cast<ESIMCoupling*>(owner)->stageAccepted(stream,reset);};
   op.commitMaterial=[](void* owner)noexcept{static_cast<ESIMCoupling*>(owner)->commit();};
   op.validateCandidateTemperature=[](void* owner,cudaStream_t stream,const double* t){static_cast<ESIMCoupling*>(owner)->validateTemperature(stream,t);};
  }
  if(const auto* requested=mxGetField(cfg,0,"validation_work_snapshot")){
   if(!mxIsLogicalScalar(requested))throw std::invalid_argument("Snapshot option must be scalar logical");
   if(mxIsLogicalScalarTrue(requested)){
    if(!material||a.rows>4096||a.z.size()>=INT_MAX)throw std::invalid_argument("Snapshot requires bounded element ESIM and int32 CSR");
    if(a.rowptr.size()!=a.rows+1||a.rowptr.front()!=0||a.rowptr.back()!=int32_t(a.z.size()))throw std::invalid_argument("Snapshot CSR count mismatch");
    for(std::size_t row=0;row<a.rows;++row){if(a.rowptr[row]>a.rowptr[row+1])throw std::invalid_argument("Nonmonotone snapshot CSR");int32_t previous=-1;for(int32_t index=a.rowptr[row];index<a.rowptr[row+1];++index){if(a.col[index]<=previous||a.col[index]<0||std::size_t(a.col[index])>=a.columns)throw std::invalid_argument("Snapshot CSR must be sorted and unique");previous=a.col[index];}}
    auto indices=realVector(field(cfg,"validation_work_coil_port_indices"),2);
    for(double value:indices){if(value<1||value>double(a.rows)||value!=std::floor(value))throw std::invalid_argument("Checked one-based coil/port indices required");snapshotIndices.push_back(int32_t(value-1));}
    if(snapshotIndices[0]==snapshotIndices[1]||cuCreal(b[snapshotIndices[1]])!=1.||cuCimag(b[snapshotIndices[1]])!=0.)throw std::invalid_argument("Invalid port/current constraint mapping");
    for(int32_t index=a.rowptr[snapshotIndices[0]];index<a.rowptr[snapshotIndices[0]+1];++index)if(a.col[index]==snapshotIndices[0])suppliedCoilImpedance=a.z[index];
    const auto* dimsInput=field(cfg,"validation_work_block_dims");auto dims=realVector(dimsInput,mxGetNumberOfElements(dimsInput));std::size_t sum=0;for(double value:dims){if(value<1||value>double(a.rows)||value!=std::floor(value))throw std::invalid_argument("Positive block dimensions required");snapshotBlockDims.push_back(int32_t(value));sum+=std::size_t(value);}if(sum!=a.rows)throw std::invalid_argument("Snapshot block dimensions must partition unknowns");
    const auto* units=field(cfg,"validation_work_block_units");if(!mxIsCell(units)||mxGetNumberOfElements(units)!=snapshotBlockDims.size())throw std::invalid_argument("Block units must match checked partition");for(std::size_t index=0;index<snapshotBlockDims.size();++index){auto unit=text(mxGetCell(units,index));if(unit.empty()||unit.size()>128)throw std::invalid_argument("Explicit bounded physical block units required");snapshotUnits.push_back(unit);}
    for(std::size_t element=0;element<samples;++element){if(w[element*std::size_t(components)]<=0.)throw std::invalid_argument("Snapshot physical half-area must be positive");for(std::size_t component=1;component<std::size_t(components);++component)if(w[element*std::size_t(components)+component]!=w[element*std::size_t(components)])throw std::invalid_argument("Snapshot requires repeated physical half-area weights");}
    snapshot=std::make_unique<WorkSnapshot>();using K=WorkSnapshot::Kind;
    snapshot->add("certified_norms_b_r_x",3,K::Real);snapshot->add("x",a.rows,K::Complex);snapshot->add("rhs",a.rows,K::Complex);
    snapshot->add("T_material_eval",n,K::Real);snapshot->add("T_accepted_after",n,K::Real);snapshot->add("heat_Qe",samples,K::Real);
    material->coupling->configureSnapshot(*snapshot);snapshot->allocate(workspace+a.rows*sizeof(cuDoubleComplex),a.rows);owned+=snapshot->bytes();op.snapshot=snapshot.get();
    op.stageWorkMaterial=[](void* owner,WorkSnapshot& bank,cudaStream_t stream,WorkSnapshot::Meta& meta){static_cast<ESIMCoupling*>(owner)->stageSnapshot(bank,stream,meta);};
   }
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
  if(!nrhs)throw std::invalid_argument("Missing IH command");const auto command=text(in[0]);
  if(command=="destroy"){
   if(nrhs!=2||nlhs)throw std::invalid_argument("destroy requires handle and no output");const auto h=handle(in[1]);
   if(!contexts.erase(h))throw std::invalid_argument("IH stale or wrong-type handle");mexUnlock();return;
  }
  mxInitGPU();
  if(command=="assemble_peec_bem"){if(nrhs!=2||nlhs!=1)throw std::invalid_argument("assemble_peec_bem requires geometry and one output");out[0]=assemblePEECBEM(in[1]);return;}
  if(command=="create"){
   if(nrhs!=2||nlhs!=1)throw std::invalid_argument("create requires config and one output");
   if(generation==0x0000ffffffffffffULL)throw std::overflow_error("IH handle generation exhausted");
   auto context=std::make_unique<Context>(in[1]);const uint64_t h=0x4948000000000000ULL|++generation;
   out[0]=mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);*static_cast<uint64_t*>(mxGetData(out[0]))=h;
   contexts.emplace(h,std::move(context));mexLock();mexAtExit(cleanup);return;
  }
  if(command=="live_count"){if(nrhs!=1||nlhs!=1)throw std::invalid_argument("live_count takes no arguments");out[0]=mxCreateDoubleScalar(double(contexts.size()));return;}
  if(nrhs<2)throw std::invalid_argument("Missing IH handle");auto& c=get(in[1]);
  if(command=="reset"){if(nrhs!=2||nlhs)throw std::invalid_argument("reset signature mismatch");try{c.state->reset();}catch(...){c.recoverFailure();throw;}return;}
#ifdef RADIA_IH_TRANSACTION_TEST
  if(command=="transaction_test_failure"){
   if(nrhs!=3||nlhs)throw std::invalid_argument("Test failure stage signature");const auto stage=number(in[2]);if(stage<1||stage>7||stage!=std::floor(stage))throw std::invalid_argument("Test stage must be integer1..7");c.state->setTestFailure(unsigned(stage));return;
  }
  if(command=="transaction_test_poison"){
   if(nrhs!=2||nlhs)throw std::invalid_argument("Test poison signature");c.poisoned=true;return;
  }
  if(command=="transaction_test_capture"){
   if(nrhs!=2||nlhs!=1||!c.material)throw std::invalid_argument("Test capture needs element material");
   out[0]=mxCreateStructMatrix(1,1,0,nullptr);
   auto put=[&](const char* name,mxArray* value){mxAddField(out[0],name);mxSetField(out[0],0,name,value);};
   auto raw=[&](const char* name,const void* pointer,std::size_t bytes){auto* value=mxCreateNumericMatrix(bytes,1,mxUINT8_CLASS,mxREAL);if(bytes)cudaCheck(cudaMemcpy(mxGetData(value),pointer,bytes,cudaMemcpyDeviceToHost));put(name,value);};
   auto count=[&](const char* name,uint64_t number){auto* value=mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);*static_cast<uint64_t*>(mxGetData(value))=number;put(name,value);};
   raw("em",c.state->testAcceptedEM(),c.a.rows*sizeof(cuDoubleComplex));raw("temperature",c.state->testAcceptedTemperature(),c.system.rows*sizeof(double));raw("heat",c.state->testAcceptedHeat(),c.projection.columns*sizeof(double));
   raw("cells",c.material->coupling->acceptedCellsDevice(),c.material->evaluation.rows*sizeof(ESIMCellResult));raw("loss_weights",c.material->coupling->acceptedWeights(),c.p.rows*sizeof(double));raw("matrix_values",c.material->coupling->acceptedValues(),c.a.z.size()*sizeof(cuDoubleComplex));
   raw("snapshot_bank",c.snapshot?c.snapshot->testAcceptedBank():nullptr,c.snapshot?c.snapshot->testBankBytes():0);
   count("steps",c.state->steps());count("snapshot_step",c.snapshot&&c.snapshot->valid()?c.snapshot->meta().step:0);count("snapshot_generation",c.snapshot?c.snapshot->epoch():0);count("handle_generation",handle(in[1])&0x0000ffffffffffffULL);put("snapshot_valid",mxCreateLogicalScalar(c.snapshot&&c.snapshot->valid()));
   cudaDeviceProp properties{};cudaCheck(cudaGetDeviceProperties(&properties,c.device));std::string uuid;const char* hex="0123456789abcdef";for(unsigned char value:properties.uuid.bytes){uuid+=hex[value>>4];uuid+=hex[value&15];}put("device_uuid",mxCreateString(uuid.c_str()));count("handle_token",handle(in[1]));count("n_rows",c.a.rows);count("nnz",c.a.z.size());count("thermal_rows",c.system.rows);count("heat_samples",c.projection.columns);count("cell_count",c.material->evaluation.rows);count("cell_struct_size_bytes",sizeof(ESIMCellResult));count("complex_element_bytes",sizeof(cuDoubleComplex));
   put("schema_version",mxCreateString("radia.ih.transaction-test-capture.v1"));put("representation",mxCreateString("uint8 raw device bytes; little-endian Windows; complex interleaved FP64; cells include native padding; same MEX only; private capture; hash externally"));return;
  }
#endif
  if(command=="accepted_work_snapshot"){
   if(nrhs!=2||nlhs!=1||!c.snapshot||!c.snapshot->valid()||c.snapshot->meta().step!=c.state->steps())throw std::invalid_argument("Enabled valid same-step accepted snapshot required");
   const auto start=std::chrono::steady_clock::now();auto& bank=*c.snapshot;
   out[0]=mxCreateStructMatrix(1,1,0,nullptr);
   auto put=[&](const std::string& name,mxArray* value){mxAddField(out[0],name.c_str());mxSetField(out[0],0,name.c_str(),value);};
   auto real=[&](const std::string& name,const std::vector<double>& value){auto* a=mxCreateDoubleMatrix(value.size(),1,mxREAL);std::copy(value.begin(),value.end(),mxGetDoubles(a));put(name,a);};
   auto integer=[&](const std::string& name,const std::vector<int32_t>& value){auto* a=mxCreateNumericMatrix(value.size(),1,mxINT32_CLASS,mxREAL);std::copy(value.begin(),value.end(),static_cast<int32_t*>(mxGetData(a)));put(name,a);};
   auto count=[&](const char* name,uint64_t value){auto* a=mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);*static_cast<uint64_t*>(mxGetData(a))=value;put(name,a);};
   auto scalar=[&](const char* name,double value){put(name,mxCreateDoubleScalar(value));};
   auto complex=[&](const std::string& name,const std::vector<cuDoubleComplex>& value){std::vector<double> re,im;re.reserve(value.size());im.reserve(value.size());for(auto z:value){re.push_back(cuCreal(z));im.push_back(cuCimag(z));}real(name+"_re",re);real(name+"_im",im);};
   put("schema_version",mxCreateString("radia.ih.accepted-work.v2.1"));put("snapshot_valid",mxCreateLogicalScalar(true));count("snapshot_step",bank.meta().step);count("snapshot_generation",bank.epoch());count("handle_token",handle(in[1]));count("handle_generation",handle(in[1])&0x0000ffffffffffffULL);count("steps",c.state->steps());
   scalar("device_ordinal",c.device);cudaDeviceProp properties{};cudaCheck(cudaGetDeviceProperties(&properties,c.device));std::string uuid;const char* hex="0123456789abcdef";for(unsigned char value:properties.uuid.bytes){uuid+=hex[value>>4];uuid+=hex[value&15];}put("device_uuid",mxCreateString(uuid.c_str()));
   put("material_mode",mxCreateString("element-esim-operators"));put("phasor_convention",mxCreateString("peak exp(+i omega t)"));put("native_build_id",mxCreateString("external exact-source MEX receipt required"));
   count("n_rows",c.a.rows);count("n_cols",c.a.columns);count("nnz",c.a.z.size());integer("csr_rowptr",c.a.rowptr);integer("csr_col",c.a.col);
   put("csr_layout",mxCreateString("zero-based sorted unique row-major; split real/imag FP64"));put("last_linear_operator_saved",mxCreateLogicalScalar(false));put("condition_estimate",mxCreateString("none; independent CPU validation only"));
   scalar("linear_tolerance",bank.meta().linearTolerance);scalar("constitutive_tolerance",bank.meta().constitutiveTolerance);scalar("field_floor_A_per_m",bank.meta().fieldFloor);scalar("constitutive_relaxation",bank.meta().relaxation);count("linear_max_iterations",bank.meta().linearMaxIterations);count("constitutive_max_iterations",bank.meta().constitutiveMaxIterations);count("linear_iterations",bank.meta().linearIterations);count("constitutive_iterations",bank.meta().constitutiveIterations);scalar("controls_current",bank.meta().current);scalar("controls_ambient",bank.meta().ambient);scalar("linear_true_residual",bank.meta().linear);scalar("certified_true_residual_native",bank.meta().certified);scalar("field_fixed_point_residual",bank.meta().field);
   put("residual_definition",mxCreateString("certified norm(Ax-b)/norm(b), absolute norm(Ax) if b=0; linear residual belongs to unsaved lagged operator"));
   for(const auto& slot:bank.slots()){
    if(slot.kind==WorkSnapshot::Kind::Real)real(slot.name,bank.read<double>(slot));
    else if(slot.kind==WorkSnapshot::Kind::Complex)complex(slot.name,bank.read<cuDoubleComplex>(slot));
    else if(slot.kind==WorkSnapshot::Kind::Cells){auto values=bank.read<ESIMCellResult>(slot);std::vector<cuDoubleComplex> z;std::vector<double> r;std::vector<int32_t> status;for(const auto& cell:values){z.push_back(cell.impedance);r.push_back(cell.trueResidual);status.push_back(cell.status);}complex("Zs",z);real("cell_true_residual",r);integer("cell_status",status);}
    else{auto values=bank.read<ESIMCellParameters>(slot);std::vector<double> rho,omega,thickness,tolerance,relaxation,saturation;std::vector<int32_t> nodes,cylinder;std::vector<cuDoubleComplex> mu;for(const auto& cell:values){rho.push_back(cell.resistivity);omega.push_back(cell.omega);thickness.push_back(cell.halfThickness);tolerance.push_back(cell.tolerance);relaxation.push_back(cell.relaxation);saturation.push_back(cell.law.saturationMagnetization);nodes.push_back(cell.nodes);cylinder.push_back(cell.cylinder?1:0);mu.push_back(cell.law.initialMu);}real("evaluated_resistivity_ohm_m",rho);real("omega",omega);real("cell_half_thickness_m",thickness);real("cell_tolerance",tolerance);real("cell_relaxation",relaxation);real("saturation_magnetization_A_per_m",saturation);integer("cell_nodes",nodes);integer("cell_cylinder",cylinder);complex("initial_mu",mu);}
   }
   integer("loss_projection_rowptr",c.p.rowptr);integer("loss_projection_col",c.p.col);complex("loss_projection",c.p.z);count("loss_projection_rows",c.p.rows);count("loss_projection_cols",c.p.columns);count("loss_components",c.p.rows/c.material->evaluation.rows);integer("port_indices",c.snapshotIndices);integer("block_dims",c.snapshotBlockDims);auto* units=mxCreateCellMatrix(c.snapshotUnits.size(),1);for(std::size_t index=0;index<c.snapshotUnits.size();++index)mxSetCell(units,index,mxCreateString(c.snapshotUnits[index].c_str()));put("block_units",units);complex("supplied_coil_impedance",std::vector<cuDoubleComplex>{c.suppliedCoilImpedance});
   integer("temperature_evaluation_rowptr",c.material->evaluation.rowptr);integer("temperature_evaluation_col",c.material->evaluation.col);real("temperature_evaluation",c.material->evaluation.real);
   std::vector<double> base(c.p.rows);cudaCheck(cudaMemcpy(base.data(),c.material->baseWeights.p,base.size()*sizeof(double),cudaMemcpyDeviceToHost));real("base_loss_weights",base);std::vector<double> areas;const auto components=c.p.rows/c.material->evaluation.rows;for(std::size_t element=0;element<c.material->evaluation.rows;++element)areas.push_back(2.*base[element*components]);real("face_areas_m2",areas);count("temperature_evaluation_rows",c.material->evaluation.rows);count("temperature_evaluation_cols",c.material->evaluation.columns);
   count("snapshot_bytes",bank.bytes());count("atomicity_fix_bytes",c.state->atomicityBytes()+c.material->coupling->atomicityBytes());count("owned_device_bytes",c.owned);put("owned_setup_peak_bytes",mxCreateDoubleMatrix(0,0,mxREAL));put("setup_peak_applicability",mxCreateString("not yet instrumented; owned_device_bytes excludes freed setup temporaries"));scalar("t_d2d_sync_s",bank.meta().stageSeconds);scalar("t_d2d_gpu_s",bank.meta().gpuCopySeconds);std::size_t transferCount=1,transferBytes=base.size()*sizeof(double);for(const auto& slot:bank.slots())if(slot.bytes){++transferCount;transferBytes+=slot.bytes;}count("d2h_transfer_count",transferCount);count("d2h_transfer_bytes",transferBytes);scalar("t_d2h_s",std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());return;
  }
  if(command=="diagnostic_operator"){
   if(nrhs!=2||nlhs!=1)throw std::invalid_argument("diagnostic_operator returns candidate operator, not accepted state");
   std::vector<cuDoubleComplex> values(c.a.z.size());
   cudaCheck(cudaDeviceSynchronize());cudaCheck(cudaMemcpy(values.data(),c.a.dz.p,values.size()*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost));
   out[0]=mxCreateDoubleMatrix(c.a.rows,c.a.columns,mxCOMPLEX);
   for(size_t row=0;row<c.a.rows;++row)for(int32_t k=c.a.rowptr[row];k<c.a.rowptr[row+1];++k){auto& z=mxGetComplexDoubles(out[0])[row+c.a.rows*c.a.col[k]];z.real=cuCreal(values[k]);z.imag=cuCimag(values[k]);}return;
  }
  if(command=="update"){if(nrhs!=4||nlhs)throw std::invalid_argument("update requires current and ambient");try{c.state->update(number(in[2]),number(in[3]),c.solveTolerance);}catch(...){c.recoverFailure();throw;}return;}
  if(command=="element_impedance"){
   if(nrhs!=2||nlhs!=1||!c.material)throw std::invalid_argument("element_impedance requires element material context and one output");
   const auto values=c.material->coupling->acceptedCells();out[0]=mxCreateDoubleMatrix(values.size(),1,mxCOMPLEX);
   for(size_t i=0;i<values.size();++i){mxGetComplexDoubles(out[0])[i].real=cuCreal(values[i].impedance);mxGetComplexDoubles(out[0])[i].imag=cuCimag(values[i].impedance);}return;
  }
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
