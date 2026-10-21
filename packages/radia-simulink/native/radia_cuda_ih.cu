#include "mex_error_boundary.hpp"
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
 Matrix evaluation,entryMap;Buffer axis,rho,mu,saturation,knots,coefficients,geometry,fields,base,baseWeights;
 std::unique_ptr<ESIMCoupling> coupling;
 std::size_t owned=0;
 void abandon(){if(coupling)coupling->abandonAfterDeviceReset();evaluation.abandon();entryMap.abandon();for(auto* b:{&axis,&rho,&mu,&saturation,&knots,&coefficients,&geometry,&fields,&base,&baseWeights})b->p=nullptr;}
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
  std::vector<double> bhKnots,bhCoefficients;int knotCount=0;
  if(mxGetField(cfg,0,"bh_knots")||mxGetField(cfg,0,"bh_coefficients")){
   const auto* inputKnots=field(cfg,"bh_knots");const auto count=mxGetNumberOfElements(inputKnots);
   if(count<2||count>256)throw std::invalid_argument("BH requires 2..256 H knots");
   bhKnots=realVector(inputKnots,count);if(bhKnots[0]!=0.)throw std::invalid_argument("BH H knots must start at zero");
   for(size_t i=1;i<count;++i)if(bhKnots[i]<=bhKnots[i-1])throw std::invalid_argument("BH H knots must increase");
   const size_t intervals=count-1,stride=intervals*4;
   if(elements*nt>(2ULL<<30)/sizeof(double)/stride)throw std::overflow_error("BH coefficient payload exceeds 2 GiB");
   const auto* input=field(cfg,"bh_coefficients");const auto* shape=mxGetDimensions(input);
   if(mxGetNumberOfDimensions(input)!=4||shape[0]!=elements||shape[1]!=nt||shape[2]!=intervals||shape[3]!=4)
    throw std::invalid_argument("BH coefficients require [element,temperature,interval,4] descending powers");
   auto raw=realVector(input,elements*nt*stride);bhCoefficients.resize(raw.size());
   for(size_t e=0;e<elements;++e)for(size_t t=0;t<nt;++t)for(size_t i=0;i<intervals;++i)for(size_t c=0;c<4;++c)
    bhCoefficients[(e*nt+t)*stride+4*i+c]=raw[e+elements*(t+nt*(i+intervals*c))];
   for(size_t row=0;row<elements*nt;++row){
    const auto* coeff=bhCoefficients.data()+row*stride;
    if(fabs(coeff[3])>1e-12||fabs(coeff[2]-cuCreal(permeability[row]))>1e-8*fmax(1e-12,fabs(coeff[2]))||cuCimag(permeability[row])!=0.)
     throw std::invalid_argument("Nonlinear BH must start at B=0 and match real initial permeability");
    for(size_t i=0;i<intervals;++i){const auto* c=coeff+4*i;double h=bhKnots[i+1]-bhKnots[i];
     auto derivative=[&](double x){return (3*c[0]*x+2*c[1])*x+c[2];};
     double smallest=std::min(derivative(0.),derivative(h));
     if(c[0]!=0.){double x=-c[1]/(3*c[0]);if(x>0.&&x<h)smallest=std::min(smallest,derivative(x));}
     if(!std::isfinite(smallest)||smallest<=0.)throw std::invalid_argument("BH differential permeability must remain finite and positive");
     double end=((c[0]*h+c[1])*h+c[2])*h+c[3];
     double next=i+1<intervals?coeff[4*(i+1)+3]:1.2566370614359172954e-6*(bhKnots.back()+sat[row]);
     if(!std::isfinite(end)||!std::isfinite(next)||fabs(end-next)>1e-8*fmax(1e-12,std::max(fabs(end),fabs(next))))throw std::invalid_argument("BH interval/saturation continuity mismatch");
    }
   }
   knotCount=int(count);
  }
  const auto* inputGeometry=field(cfg,"esim_cell_parameters");if(mxGetM(inputGeometry)!=elements||mxGetN(inputGeometry)!=7)throw std::invalid_argument("ESIM cells require element by 7 parameters");auto values=realVector(inputGeometry,elements*7);
  std::vector<ESIMCellParameters> cells(elements);int maxNodes=0;for(size_t e=0;e<elements;++e){auto at=[&](size_t k){return values[e+elements*k];};auto& cell=cells[e];
   if(at(0)<3||at(0)>4096||at(0)!=int(at(0))||at(1)<1||at(1)>10000||at(1)!=int(at(1))||(at(2)!=0.&&at(2)!=1.)||at(3)<=0.||at(4)<=0.||at(5)<=0.||at(5)>1e-6||at(6)<=0.||at(6)>1.)throw std::invalid_argument("Invalid ESIM cell geometry/iteration contract");
   cell.nodes=int(at(0));cell.maxIterations=int(at(1));cell.cylinder=at(2)!=0.;cell.halfThickness=at(3);cell.omega=at(4);cell.tolerance=at(5);cell.relaxation=at(6);maxNodes=std::max(maxNodes,cell.nodes);}
  auto guesses=realVector(field(cfg,"initial_field_A_per_m"),elements);for(auto x:guesses)if(x<=0.)throw std::invalid_argument("ESIM initial fields must be positive");
  double iterations=number(field(cfg,"max_constitutive_iterations"));if(iterations<2||iterations>10000||iterations!=int(iterations))throw std::invalid_argument("Invalid constitutive iteration limit");
  owned=evaluation.bytes()+entryMap.bytes()+temperatures.size()*sizeof(double)+resistivity.size()*sizeof(double)+permeability.size()*sizeof(cuDoubleComplex)+sat.size()*sizeof(double)+cells.size()*sizeof(ESIMCellParameters)+guesses.size()*sizeof(double)+a.z.size()*sizeof(cuDoubleComplex)+w.size()*sizeof(double)+(bhKnots.size()+bhCoefficients.size())*sizeof(double);
  size_t available,total;cudaCheck(cudaMemGetInfo(&available,&total));auto reserve=std::max<size_t>(2ULL<<30,total/5);if(available<=reserve||owned>available-reserve)throw std::runtime_error("Element material upload exceeds GPU reserve");
  evaluation.upload();entryMap.upload();axis.put(temperatures);rho.put(resistivity);mu.put(permeability);saturation.put(sat);knots.put(bhKnots);coefficients.put(bhCoefficients);geometry.put(cells);fields.put(guesses);base.put(a.z);baseWeights.put(w);
  ESIMCouplingView v{};v.temperatureEvaluation=evaluation.realView();v.table={axis.get<double>(),rho.get<double>(),saturation.get<double>(),knots.get<double>(),coefficients.get<double>(),mu.get<cuDoubleComplex>(),int(nt),knotCount};v.entryMap=entryMap.complexView();v.fieldProjection=p.complexView();v.matrix=a.complexView();v.baseValues=base.get<cuDoubleComplex>();v.mutableValues=static_cast<cuDoubleComplex*>(a.dz.p);v.geometry=geometry.get<ESIMCellParameters>();v.initialFields=fields.get<double>();v.baseLossWeights=baseWeights.get<double>();v.mutableLossWeights=static_cast<double*>(weights.p);v.elements=elements;v.components=components;v.workspaceNodes=maxNodes;v.maxIterations=int(iterations);v.tolerance=number(field(cfg,"constitutive_tolerance"));v.relaxation=number(field(cfg,"constitutive_relaxation"));v.fieldFloor=number(field(cfg,"field_floor_A_per_m"));coupling=std::make_unique<ESIMCoupling>(v);owned+=coupling->ownedBytes();
 }
 static bool cellFiniteHost(cuDoubleComplex z){return std::isfinite(cuCreal(z))&&std::isfinite(cuCimag(z));}
};
// IH-only nonlinear heat owner. Caller/NGSolve supplies all FE evaluation
// maps. Bounded dense Newton updates do not use shared AMG/IC kernels.
struct QuadratureThermal {
 ThermalQuadratureView v{};cuDoubleComplex* matrix=nullptr;
 Buffer axis,basis,gradient,volume,density,cp,kappa,boundaryBasis,area,h,emissivity;
 Buffer current,volumeWork,boundaryWork,stats,invalid,constant;
 double tolerance=1e-10,energyLimit=1e-6,floor=1e-12;int maxIterations=30;
 double candidateStats[6]{},acceptedStats[6]{};unsigned candidateIterations=0,acceptedIterations=0;
 size_t owned=0;
 std::vector<double> rows(const mxArray* a,size_t nr,size_t nc){
  if(mxGetNumberOfDimensions(a)!=2||mxGetM(a)!=nr||mxGetN(a)!=nc)throw std::invalid_argument("IH quadrature table dimensions disagree");
  auto raw=realVector(a,nr*nc);std::vector<double> out(raw.size());
  for(size_t i=0;i<nr;++i)for(size_t j=0;j<nc;++j)out[i*nc+j]=raw[i+nr*j];return out;
 }
 explicit QuadratureThermal(const mxArray* cfg,Matrix& system,const std::vector<double>& one){
  v.n=system.rows;if(v.n>256||system.z.size()!=v.n*v.n)throw std::invalid_argument("Dynamic thermal reference requires complete dense pattern, <=256 coefficients");
  for(size_t i=0;i<v.n;++i)for(size_t j=0;j<v.n;++j)if(system.col[i*v.n+j]!=int(j))throw std::invalid_argument("Dynamic thermal dense CSR pattern required");
  const auto* e=field(cfg,"thermal_basis");v.q=mxGetM(e);
  const auto* b=field(cfg,"thermal_boundary_basis");v.b=mxGetM(b);
  double dim=number(field(cfg,"thermal_dimension"));if(dim<1||dim>3||dim!=int(dim)||!v.q||v.q>4096||!v.b||v.b>4096)throw std::invalid_argument("IH quadrature dimensions out of range");v.dimension=int(dim);
  const auto* ts=field(cfg,"thermal_temperature_axis_K");const auto nt=mxGetNumberOfElements(ts);if(nt<2||nt>1024)throw std::invalid_argument("Thermal temperature axis requires2..1024 entries");v.nt=int(nt);
  auto t=realVector(ts,v.nt);for(int j=0;j<v.nt;++j)if(t[j]<=0.||(j&&t[j]<=t[j-1]))throw std::invalid_argument("Thermal Kelvin axis must increase");
  auto ev=rows(e,v.q,v.n),grad=rows(field(cfg,"thermal_gradient"),v.q*v.dimension,v.n),bv=rows(b,v.b,v.n);
  auto vw=realVector(field(cfg,"thermal_volume_weights_m3"),v.q),rho=realVector(field(cfg,"thermal_density_kg_per_m3"),v.q);
  auto cv=rows(field(cfg,"thermal_cp_J_per_kgK"),v.q,v.nt),kv=rows(field(cfg,"thermal_kappa_W_per_mK"),v.q,v.nt);
  auto aw=realVector(field(cfg,"thermal_boundary_weights_m2"),v.b),hv=realVector(field(cfg,"thermal_convection_W_per_m2K"),v.b),eps=rows(field(cfg,"thermal_emissivity"),v.b,v.nt);
  for(auto x:vw)if(x<=0.)throw std::invalid_argument("Positive volume weights required");for(auto x:rho)if(x<=0.)throw std::invalid_argument("Positive mass density required");
  for(auto x:cv)if(x<=0.)throw std::invalid_argument("Positive cp required");for(auto x:kv)if(x<=0.)throw std::invalid_argument("Positive kappa required");
  for(auto x:aw)if(x<=0.)throw std::invalid_argument("Positive surface weights required");for(auto x:hv)if(x<0.)throw std::invalid_argument("Nonnegative convection required");for(auto x:eps)if(x<0.||x>1.)throw std::invalid_argument("Emissivity outside[0,1]");
  auto preserves=[&](const std::vector<double>& map,size_t nr,double expected){for(size_t i=0;i<nr;++i){double s=0.;for(size_t j=0;j<v.n;++j)s+=map[i*v.n+j]*one[j];if(fabs(s-expected)>1e-12)throw std::invalid_argument("Quadrature maps do not preserve the constant field");}};
  preserves(ev,v.q,1.);preserves(bv,v.b,1.);preserves(grad,v.q*v.dimension,0.);
  tolerance=number(field(cfg,"thermal_nonlinear_tolerance"));energyLimit=number(field(cfg,"thermal_energy_relative_limit"));floor=number(field(cfg,"thermal_energy_power_floor_W"));
  double count=number(field(cfg,"thermal_max_iterations"));if(tolerance<=0.||tolerance>1e-6||energyLimit<=0.||energyLimit>1e-6||floor<=0.||floor>1e-12||count<1||count>100||count!=int(count))throw std::invalid_argument("Invalid thermal acceptance limits");maxIterations=int(count);
  const size_t need=(t.size()+ev.size()+grad.size()+vw.size()+rho.size()+cv.size()+kv.size()+bv.size()+aw.size()+hv.size()+eps.size()+2*v.n+v.q*(4+v.dimension)+2*v.b+6)*sizeof(double)+sizeof(int);
  size_t free,total;cudaCheck(cudaMemGetInfo(&free,&total));const auto reserve=std::max<size_t>(2ULL<<30,total/5);if(free<=reserve||need>free-reserve)throw std::runtime_error("IH thermal quadrature exceeds device reserve");
  axis.put(t);basis.put(ev);gradient.put(grad);volume.put(vw);density.put(rho);cp.put(cv);kappa.put(kv);boundaryBasis.put(bv);area.put(aw);h.put(hv);emissivity.put(eps);constant.put(one);
  current.put(std::vector<double>(v.n));volumeWork.put(std::vector<double>(v.q*(4+v.dimension)));boundaryWork.put(std::vector<double>(2*v.b));stats.put(std::vector<double>(6));invalid.put(std::vector<int>{0});owned=need;
  v.axis=axis.get<double>();v.basis=basis.get<double>();v.gradient=gradient.get<double>();v.volume=volume.get<double>();v.density=density.get<double>();v.cp=cp.get<double>();v.kappa=kappa.get<double>();v.boundaryBasis=boundaryBasis.get<double>();v.area=area.get<double>();v.convection=h.get<double>();v.emissivity=emissivity.get<double>();matrix=static_cast<cuDoubleComplex*>(system.dz.p);
 }
 void abandon(){for(auto* p:{&axis,&basis,&gradient,&volume,&density,&cp,&kappa,&boundaryBasis,&area,&h,&emissivity,&current,&volumeWork,&boundaryWork,&stats,&invalid,&constant})p->p=nullptr;matrix=nullptr;}
 ComplexSolver::Result solve(ComplexSolver& solver,const double* old,const double* load,double ambient,double dt,double* candidate){
  if(!std::isfinite(ambient)||ambient<=0.)throw std::invalid_argument("Ambient must be positive Kelvin");
  auto stream=solver.stream();cudaCheck(cudaMemcpyAsync(current.p,old,v.n*sizeof(double),cudaMemcpyDeviceToDevice,stream));
  ComplexSolver::Result last{};unsigned total=0;
  for(int iteration=0;iteration<maxIterations;++iteration){
   cudaCheck(cudaMemsetAsync(invalid.p,0,sizeof(int),stream));
   thermalEvaluateQuadrature<<<unsigned((std::max(v.q,v.b)+127)/128),128,0,stream>>>(v,current.get<double>(),old,ambient,dt,static_cast<double*>(volumeWork.p),static_cast<double*>(boundaryWork.p),static_cast<int*>(invalid.p));
   int bad=0;cudaCheck(cudaMemcpyAsync(&bad,invalid.p,sizeof(int),cudaMemcpyDeviceToHost,stream));cudaCheck(cudaStreamSynchronize(stream));if(bad)throw std::runtime_error("Thermal table range/finite candidate rejected; old state preserved");
   thermalAssembleQuadrature<<<unsigned((v.n+127)/128),128,0,stream>>>(v,volumeWork.get<double>(),boundaryWork.get<double>(),load,matrix,solver.rhs(),static_cast<int*>(invalid.p));
   thermalQuadratureBalance<<<1,1,0,stream>>>(v,volumeWork.get<double>(),boundaryWork.get<double>(),load,constant.get<double>(),solver.rhs(),floor,static_cast<double*>(stats.p),static_cast<int*>(invalid.p));
   cudaCheck(cudaMemcpyAsync(candidateStats,stats.p,sizeof(candidateStats),cudaMemcpyDeviceToHost,stream));cudaCheck(cudaMemcpyAsync(&bad,invalid.p,sizeof(int),cudaMemcpyDeviceToHost,stream));cudaCheck(cudaStreamSynchronize(stream));
   if(bad)throw std::runtime_error("Nonfinite nonlinear thermal residual/Jacobian");
   if(candidateStats[0]<=tolerance&&candidateStats[1]<=energyLimit){cudaCheck(cudaMemcpyAsync(candidate,current.p,v.n*sizeof(double),cudaMemcpyDeviceToDevice,stream));candidateIterations=unsigned(iteration);return {total,last.trueRelativeResidual};}
   last=solver.solve(tolerance,2000);total+=last.iterations;
   thermalCandidate<<<unsigned((v.n+127)/128),128,0,stream>>>(current.get<double>(),solver.solution(),candidate,static_cast<int*>(invalid.p),v.n);
   cudaCheck(cudaMemcpyAsync(&bad,invalid.p,sizeof(int),cudaMemcpyDeviceToHost,stream));cudaCheck(cudaStreamSynchronize(stream));if(bad)throw std::runtime_error("Invalid thermal Newton increment");
   cudaCheck(cudaMemcpyAsync(current.p,candidate,v.n*sizeof(double),cudaMemcpyDeviceToDevice,stream));
  }
  throw std::runtime_error("Nonlinear thermal/energy acceptance failed; old state preserved");
 }
 void commit(){std::copy(candidateStats,candidateStats+6,acceptedStats);acceptedIterations=candidateIterations;}
};
struct Context {
 int device=0;cudaStream_t guard=nullptr;
 Matrix a,p,system,k,c,projection,preconditioner;
 PreconditionerPrecision precision=PreconditionerPrecision::None;Buffer rhs,weights,ones,initial,voltage;
 std::unique_ptr<ElementMaterial> material;
 std::unique_ptr<QuadratureThermal> thermalMaterial;
 std::unique_ptr<CoupledState> state;size_t owned=0;
 ~Context(){
  int previous=device;cudaGetDevice(&previous);if(previous!=device)cudaSetDevice(device);
  const auto status=guard?cudaStreamQuery(guard):cudaSuccess;
  if(status!=cudaSuccess&&status!=cudaErrorNotReady){
   if(state)state->abandonAfterDeviceReset();if(material)material->abandon();if(thermalMaterial)thermalMaterial->abandon();
   for(auto* m:{&a,&p,&system,&k,&c,&projection,&preconditioner})m->abandon();
   rhs.p=weights.p=ones.p=initial.p=voltage.p=nullptr;guard=nullptr;
  }
  state.reset();material.reset();thermalMaterial.reset();
  // Free operators while their owning device is selected.
  for(auto* m:{&a,&p,&system,&k,&c,&projection,&preconditioner})for(auto* b:{&m->dr,&m->dc,&m->dz,&m->dreal,&m->df}){if(b->p)cudaFree(b->p);b->p=nullptr;}
  for(auto* b:{&rhs,&weights,&ones,&initial,&voltage}){if(b->p)cudaFree(b->p);b->p=nullptr;}
  if(guard)cudaStreamDestroy(guard);if(previous!=device)cudaSetDevice(previous);
 }
 void check(){int current;cudaCheck(cudaGetDevice(&current));if(current!=device)throw std::runtime_error("IH context belongs to another GPU");cudaCheck(cudaStreamQuery(guard));}
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
  for(const auto value:system.z)if(cuCimag(value)!=0.)
   throw std::invalid_argument("Time-domain IH thermal system must be real");
  k.read(field(cfg,"stiffness"),false);c.read(field(cfg,"convection"),false);projection.read(field(cfg,"heat_projection"),false);
  if(blocksField){
   if(precision!=PreconditionerPrecision::FP64||a.rows>512||mxGetField(cfg,0,"preconditioner_inverse"))throw std::invalid_argument("GPU block setup requires FP64, <=512 unknowns and no supplied inverse");
   auto ends=realVector(blocksField,mxGetNumberOfElements(blocksField));int previous=0;
   for(double x:ends){if(x<=previous||x>a.rows||x!=int(x))throw std::invalid_argument("Invalid block partition");blockEnds.push_back(int(x));previous=int(x);}if(previous!=a.rows)throw std::invalid_argument("Block partition must cover operator");
   preconditioner.rows=preconditioner.columns=a.rows;preconditioner.rowptr.resize(a.rows+1);preconditioner.col.resize(a.rows*a.rows);preconditioner.z.assign(a.rows*a.rows,make_cuDoubleComplex(0.,0.));
   for(size_t i=0;i<=a.rows;++i)preconditioner.rowptr[i]=int32_t(i*a.rows);for(size_t i=0;i<a.rows;++i)for(size_t j=0;j<a.rows;++j)preconditioner.col[i*a.rows+j]=int32_t(j);
  }else if(precision!=PreconditionerPrecision::None){preconditioner.read(field(cfg,"preconditioner_inverse"),true);
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
  const size_t workspace=(a.rows*(precision==PreconditionerPrecision::None?10:12)+n*9+2)*sizeof(cuDoubleComplex)+(2*samples+3*n+7)*sizeof(double)+2*sizeof(BiCGState)+sizeof(int)+(blocksField?a.rows*a.rows*sizeof(cuDoubleComplex)+blockEnds.size()*sizeof(int)+sizeof(int):0);
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
   op.commitMaterial=[](void* owner){static_cast<ESIMCoupling*>(owner)->commit();};op.resetMaterial=[](void* owner){static_cast<ESIMCoupling*>(owner)->reset();};
  }
  if(const auto* thermalMode=mxGetField(cfg,0,"thermal_material_mode")){
   const auto mode=text(thermalMode);
   if(mode=="quadrature-enthalpy-v1"){
    thermalMaterial=std::make_unique<QuadratureThermal>(cfg,system,o);owned+=thermalMaterial->owned;
    op.thermalOwner=thermalMaterial.get();
    op.solveThermal=[](void* owner,ComplexSolver& solver,const double* old,const double* load,double ambient,double dt,double* candidate){return static_cast<QuadratureThermal*>(owner)->solve(solver,old,load,ambient,dt,candidate);};
    op.commitThermal=[](void* owner){static_cast<QuadratureThermal*>(owner)->commit();};
    op.resetThermal=[](void* owner){auto* t=static_cast<QuadratureThermal*>(owner);std::fill(t->acceptedStats,t->acceptedStats+6,0.);t->acceptedIterations=0;};
   }else if(mode!="fixed")throw std::invalid_argument("Unsupported thermal material mode; no fixed substitution");
  }
  if(const auto* pv=mxGetField(cfg,0,"port_voltage_projection")){
   if(!mxIsDouble(pv)||mxIsSparse(pv)||mxGetNumberOfElements(pv)!=a.rows)throw std::invalid_argument("Independent voltage projection shape/type mismatch");
   std::vector<cuDoubleComplex> vp(a.rows);
   for(size_t i=0;i<a.rows;++i){const auto z=mxIsComplex(pv)?make_cuDoubleComplex(mxGetComplexDoubles(pv)[i].real,mxGetComplexDoubles(pv)[i].imag):make_cuDoubleComplex(mxGetDoubles(pv)[i],0.);if(!(std::isfinite(cuCreal(z))&&std::isfinite(cuCimag(z))))throw std::invalid_argument("Nonfinite voltage projection");vp[i]=z;}
   op.coilResistance=number(field(cfg,"coil_resistance_ohm"));op.energyRelativeLimit=number(field(cfg,"energy_input_relative_limit"));op.energyPowerFloor=number(field(cfg,"energy_power_floor_W"));
   if(op.coilResistance<0.||op.energyRelativeLimit<=0.||op.energyRelativeLimit>1e-6||op.energyPowerFloor<=0.||op.energyPowerFloor>1e-12)throw std::invalid_argument("Invalid independent IH power limits");
   voltage.put(vp);owned+=vp.size()*sizeof(cuDoubleComplex);op.voltageProjection=voltage.get<cuDoubleComplex>();
  }
  if(const auto* required=mxGetField(cfg,0,"require_four_property_coupling")){
   if(number(required)!=1.||!material||!thermalMaterial||!mxGetField(cfg,0,"bh_knots")||!op.voltageProjection)
    throw std::invalid_argument("Four-property IH requires nonlinear BH/rho tables, dynamic enthalpy/kappa and independent port power; no fixed fallback");
  }
  cudaCheck(cudaMemGetInfo(&available,&total));if(available<=reserve||workspace>available-reserve)throw std::runtime_error("IH remaining solver workspace exceeds GPU reserve");
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
void dispatchIh(int nlhs,mxArray** out,int nrhs,const mxArray** in){
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
   out[0]=vectorOut(c.state->temperature());out[1]=vectorOut(c.state->heat());const char* f[]={"steps","owned_device_bytes","eddy_iterations","eddy_true_residual","thermal_iterations","thermal_true_residual","preconditioner_precision","constitutive_iterations","constitutive_true_residual","constitutive_field_residual","material_mode","thermal_material_mode","thermal_nonlinear_residual","thermal_energy_relative_defect","thermal_input_power_W","thermal_storage_power_W","thermal_boundary_loss_W","thermal_energy_defect_W","port_input_power_W","coil_loss_power_W","integrated_heat_power_W","port_heat_defect_W","port_heat_relative_defect"};out[2]=mxCreateStructMatrix(1,1,23,f);
   mxSetField(out[2],0,f[0],mxCreateDoubleScalar(double(c.state->steps())));mxSetField(out[2],0,f[1],mxCreateDoubleScalar(double(c.owned)));
   mxSetField(out[2],0,f[2],mxCreateDoubleScalar(c.state->eddyResult().iterations));
   mxSetField(out[2],0,f[3],mxCreateDoubleScalar(c.state->eddyResult().trueRelativeResidual));
   mxSetField(out[2],0,f[4],mxCreateDoubleScalar(c.state->thermalResult().iterations));
   mxSetField(out[2],0,f[5],mxCreateDoubleScalar(c.state->thermalResult().trueRelativeResidual));
   mxSetField(out[2],0,f[6],mxCreateString(c.precision==PreconditionerPrecision::FP32?"fp32":c.precision==PreconditionerPrecision::FP64?"fp64":"none"));mxSetField(out[2],0,f[7],mxCreateDoubleScalar(c.material?c.material->coupling->iterations():0));mxSetField(out[2],0,f[8],mxCreateDoubleScalar(c.material?c.material->coupling->residual():0));mxSetField(out[2],0,f[9],mxCreateDoubleScalar(c.material?c.material->coupling->fieldResidual():0));mxSetField(out[2],0,f[10],mxCreateString(c.material?"element-esim-operators":"fixed"));
   mxSetField(out[2],0,f[11],mxCreateString(c.thermalMaterial?"quadrature-enthalpy-v1":"fixed"));
   for(int i=0;i<6;++i)mxSetField(out[2],0,f[12+i],mxCreateDoubleScalar(c.thermalMaterial?c.thermalMaterial->acceptedStats[i]:mxGetNaN()));
   for(int i=0;i<5;++i)mxSetField(out[2],0,f[18+i],mxCreateDoubleScalar(c.voltage.p?c.state->powerBalance()[i]:mxGetNaN()));
   return;
  }
  throw std::invalid_argument("Unknown IH CUDA command");

}
void mexFunction(int nl,mxArray** out,int nr,const mxArray** in){

    char errorMessage[2048]={0};
    if(!radia_mex::checked(dispatchIh,nl,out,nr,in,errorMessage,sizeof(errorMessage))){
        mexErrMsgIdAndTxt("radia:ih:cuda","%s",errorMessage);
    }
}
