// Validation-only native ESIM cell entry point; not the physical IH frontend.
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include "ih_esim_cell.cuh"
#include <vector>
#include <stdexcept>
#include <limits>
using namespace radia_ih_cuda;
namespace {
void check(cudaError_t e){if(e!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(e));}
struct Buffer {void* p=nullptr;~Buffer(){if(p)cudaFree(p);}void allocate(size_t n){check(cudaMalloc(&p,n));}
 template<class T> void put(const std::vector<T>& a){allocate(a.size()*sizeof(T));check(cudaMemcpy(p,a.data(),a.size()*sizeof(T),cudaMemcpyHostToDevice));}
 template<class T> T* get(){return static_cast<T*>(p);}};
std::vector<double> real(const mxArray* a){
 if(!mxIsDouble(a)||mxIsComplex(a)||mxIsSparse(a))throw std::invalid_argument("Expected full real double input");
 const auto* p=mxGetDoubles(a);std::vector<double> v(p,p+mxGetNumberOfElements(a));
 for(double x:v)if(!std::isfinite(x))throw std::invalid_argument("Nonfinite cell input");return v;
}
}
void mexFunction(int nlhs,mxArray** out,int nrhs,const mxArray** in) {
 try{
  if(nrhs!=4||nlhs!=3)throw std::invalid_argument("Expected fields, parameters, knots, interval-major coefficient rows; three outputs");
  mxInitGPU();auto fields=real(in[0]),v=real(in[1]),knots=real(in[2]),coeff=real(in[3]);
  if(fields.empty()||v.size()!=11||v[0]<3||v[0]>4096||v[0]!=int(v[0])||v[1]<1||v[1]>10000||v[1]!=int(v[1])||
    (v[2]!=0.&&v[2]!=1.)||v[3]<=0.||v[4]<=0.||v[5]<0.||v[6]<=0.||v[6]>1e-6||v[7]<=0.||v[7]>1.||v[8]<=0.||v[9]>0.)
      throw std::invalid_argument("Invalid finite ESIM cell parameters");
  for(double f:fields)if(f<=0.)throw std::invalid_argument("Surface fields must be positive; apply declared outer field floor explicitly");
  if(!knots.empty()){
   if(knots.size()<2||knots.size()>INT_MAX||knots[0]!=0.||v[9]!=0.||mxGetM(in[3])!=knots.size()-1||mxGetN(in[3])!=4)
    throw std::invalid_argument("Invalid real BH piecewise polynomial");
   for(size_t i=1;i<knots.size();++i)if(knots[i]<=knots[i-1])throw std::invalid_argument("BH knots must increase");
  }else if(!coeff.empty())throw std::invalid_argument("Unexpected BH coefficients");
  const size_t n=size_t(v[0]),panels=fields.size();
  if(panels>std::numeric_limits<size_t>::max()/(4*n*sizeof(cuDoubleComplex)))throw std::overflow_error("ESIM workspace size overflow");
  size_t available,total;check(cudaMemGetInfo(&available,&total));size_t reserve=total/5;if(reserve<(2ULL<<30))reserve=2ULL<<30;
  size_t need=4*n*panels*sizeof(cuDoubleComplex)+panels*(sizeof(double)+sizeof(ESIMCellResult))+knots.size()*sizeof(double)+coeff.size()*sizeof(double);
  if(available<=reserve||need>available-reserve)throw std::runtime_error("ESIM cell workspace exceeds GPU budget");
  Buffer df,dk,dc,work,dr;df.put(fields);if(!knots.empty()){dk.put(knots);std::vector<double> row(coeff.size());
    for(size_t i=0;i+1<knots.size();++i)for(int j=0;j<4;++j)row[4*i+j]=coeff[i+(knots.size()-1)*j];dc.put(row);}
  work.allocate(4*n*panels*sizeof(cuDoubleComplex));dr.allocate(panels*sizeof(ESIMCellResult));
  ESIMCellParameters params{};params.nodes=int(n);params.maxIterations=int(v[1]);params.cylinder=v[2]!=0.;
  params.halfThickness=v[3];params.resistivity=v[4];params.omega=v[5];params.tolerance=v[6];params.relaxation=v[7];
  params.law={dk.get<double>(),dc.get<double>(),int(knots.size()),make_cuDoubleComplex(v[8],v[9]),v[10]};
  finiteESIMCells<<<unsigned((panels+31)/32),32>>>(params,df.get<double>(),panels,work.get<cuDoubleComplex>(),dr.get<ESIMCellResult>());
  check(cudaGetLastError());std::vector<ESIMCellResult> result(panels);check(cudaMemcpy(result.data(),dr.p,panels*sizeof(ESIMCellResult),cudaMemcpyDeviceToHost));
  for(const auto& r:result)if(r.status!=1)throw std::runtime_error("ESIM cell failed convergence, passivity or true residual gate; status="+std::to_string(r.status));
  std::vector<cuDoubleComplex> solution(4*n*panels);check(cudaMemcpy(solution.data(),work.p,solution.size()*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost));
  out[0]=mxCreateDoubleMatrix(panels,1,mxCOMPLEX);out[1]=mxCreateDoubleMatrix(panels,6,mxREAL);out[2]=mxCreateDoubleMatrix(n,panels,mxCOMPLEX);
  for(size_t p=0;p<panels;++p){const auto& r=result[p];auto& z=mxGetComplexDoubles(out[0])[p];z.real=cuCreal(r.impedance);z.imag=cuCimag(r.impedance);
   double values[]={r.activePower,r.reactivePower,r.magneticPower,r.relativeChange,r.trueResidual,double(r.iterations)};
   for(int j=0;j<6;++j)mxGetDoubles(out[1])[p+panels*j]=values[j];
   for(size_t i=0;i<n;++i){auto h=solution[p*4*n+i];auto& target=mxGetComplexDoubles(out[2])[i+n*p];target.real=cuCreal(h);target.imag=cuCimag(h);}}
 }catch(const std::exception& e){mexErrMsgIdAndTxt("radia:ih:esimCell","%s",e.what());}
}
