// Generic three-cell material plumbing diagnostic, not a PEEC-BEM case.
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include "ih_temperature_material.cuh"
#include <vector>
#include <stdexcept>
using namespace radia_ih_cuda;
namespace {
void check(cudaError_t e){if(e!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(e));}
struct Buffer {void* p=nullptr;~Buffer(){if(p)cudaFree(p);}
 void alloc(size_t n){check(cudaMalloc(&p,n));}
 template<class T>void put(const std::vector<T>& v){alloc(v.size()*sizeof(T));check(cudaMemcpy(p,v.data(),v.size()*sizeof(T),cudaMemcpyHostToDevice));}
 template<class T>T* get(){return static_cast<T*>(p);}};
}
void mexFunction(int nlhs,mxArray** out,int nrhs,const mxArray** in){
 try{
  if(nrhs!=1||nlhs!=1||!mxIsDouble(in[0])||mxIsComplex(in[0])||mxIsSparse(in[0])||mxGetNumberOfElements(in[0])!=3)
      throw std::invalid_argument("Expected three real temperature coefficients and one output");
  mxInitGPU();constexpr int nodes=32;constexpr size_t cells=3;
  std::vector<double> t(mxGetDoubles(in[0]),mxGetDoubles(in[0])+3);
  for(double v:t)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite temperature");
  Buffer temperature,row,column,value,axis,rho,mu,saturation,geometry,candidate,local,invalid,fields,work,result;
  temperature.put(t);row.put(std::vector<int32_t>{0,2,3,4});column.put(std::vector<int32_t>{0,1,1,2});value.put(std::vector<double>{.5,.5,1.,1.});
  axis.put(std::vector<double>{300.,600.});rho.put(std::vector<double>{1e-7,4e-7,2e-7,8e-7,3e-7,12e-7});
  std::vector<cuDoubleComplex> mus;for(int e=0;e<3;++e){mus.push_back(make_cuDoubleComplex(1.2566370614359173e-4,0));mus.push_back(make_cuDoubleComplex(6.283185307179586e-5,0));}
  mu.put(mus);saturation.put(std::vector<double>(6,0.));
  std::vector<ESIMCellParameters> g(cells);for(size_t e=0;e<cells;++e){g[e].nodes=nodes;g[e].maxIterations=50;g[e].cylinder=e==1;g[e].halfThickness=.001;g[e].omega=6283.185307179586;g[e].tolerance=1e-6;g[e].relaxation=.5;}
  geometry.put(g);candidate.alloc(cells*sizeof(ESIMCellParameters));local.alloc(cells*sizeof(double));invalid.put(std::vector<int>{0});
  TemperatureMaterialTable table{axis.get<double>(),rho.get<double>(),saturation.get<double>(),nullptr,nullptr,mu.get<cuDoubleComplex>(),2,0};
  updateElementESIMMaterials<<<1,32>>>({cells,row.get<int32_t>(),column.get<int32_t>(),value.get<double>()},temperature.get<double>(),table,geometry.get<ESIMCellParameters>(),candidate.get<ESIMCellParameters>(),local.get<double>(),nullptr,invalid.get<int>());
  check(cudaGetLastError());int bad;check(cudaMemcpy(&bad,invalid.p,sizeof(int),cudaMemcpyDeviceToHost));
  if(bad)throw std::runtime_error("Material temperature outside declared table or invalid candidate; accepted state untouched");
  fields.put(std::vector<double>{100.,200.,300.});work.alloc(cells*4*nodes*sizeof(cuDoubleComplex));result.alloc(cells*sizeof(ESIMCellResult));
  finiteElementESIMCells<<<1,32>>>(candidate.get<ESIMCellParameters>(),fields.get<double>(),cells,nodes,work.get<cuDoubleComplex>(),result.get<ESIMCellResult>());
  check(cudaGetLastError());std::vector<ESIMCellResult> r(cells);check(cudaMemcpy(r.data(),result.p,cells*sizeof(ESIMCellResult),cudaMemcpyDeviceToHost));
  std::vector<ESIMCellParameters> parameters(cells);check(cudaMemcpy(parameters.data(),candidate.p,cells*sizeof(ESIMCellParameters),cudaMemcpyDeviceToHost));
  std::vector<double> mapped(cells);check(cudaMemcpy(mapped.data(),local.p,cells*sizeof(double),cudaMemcpyDeviceToHost));
  out[0]=mxCreateDoubleMatrix(cells,7,mxREAL);
  for(size_t e=0;e<cells;++e){if(r[e].status!=1)throw std::runtime_error("Element ESIM cell rejected");
   const double values[]={mapped[e],parameters[e].resistivity,cuCreal(parameters[e].law.initialMu),cuCreal(r[e].impedance),cuCimag(r[e].impedance),r[e].trueResidual,double(r[e].iterations)};
   for(size_t k=0;k<7;++k)mxGetDoubles(out[0])[e+cells*k]=values[k];}
 }catch(const std::exception& e){mexErrMsgIdAndTxt("radia:ih:temperatureSmoke","%s",e.what());}
}
