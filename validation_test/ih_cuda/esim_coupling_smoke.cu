// Generic three-cell constitutive block diagnostic, not a physical PEEC-BEM case.
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include "ih_esim_coupling.cuh"
#include "ih_coupled_state.cuh"
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
  Buffer diagonalRows,diagonalColumns,base,entries,entryValues,fieldValues,baseWeights,weights,rhs;
  diagonalRows.put(std::vector<int32_t>{0,1,2,3});diagonalColumns.put(std::vector<int32_t>{0,1,2});
  base.put(std::vector<cuDoubleComplex>(cells,make_cuDoubleComplex(2.,.1)));
  entries.put(std::vector<cuDoubleComplex>(cells,make_cuDoubleComplex(2.,.1)));
  entryValues.put(std::vector<cuDoubleComplex>(cells,make_cuDoubleComplex(10.,0.)));
  fieldValues.put(std::vector<cuDoubleComplex>(cells,make_cuDoubleComplex(1.,0.)));
  baseWeights.put(std::vector<double>(cells,.5));weights.put(std::vector<double>(cells,0.));
  ComplexCSRView matrix{cells,diagonalRows.get<int32_t>(),diagonalColumns.get<int32_t>(),entries.get<cuDoubleComplex>()};
  ESIMCouplingView v{};v.temperatureEvaluation={cells,row.get<int32_t>(),column.get<int32_t>(),value.get<double>()};v.table=table;
  v.entryMap={cells,diagonalRows.get<int32_t>(),diagonalColumns.get<int32_t>(),entryValues.get<cuDoubleComplex>()};
  v.fieldProjection={cells,diagonalRows.get<int32_t>(),diagonalColumns.get<int32_t>(),fieldValues.get<cuDoubleComplex>()};
  v.matrix=matrix;v.baseValues=base.get<cuDoubleComplex>();v.mutableValues=entries.get<cuDoubleComplex>();v.geometry=geometry.get<ESIMCellParameters>();
  v.initialFields=fields.get<double>();v.baseLossWeights=baseWeights.get<double>();v.mutableLossWeights=weights.get<double>();
  v.elements=cells;v.components=1;v.workspaceNodes=nodes;v.maxIterations=100;v.tolerance=1e-8;v.relaxation=.8;v.fieldFloor=1e-6;
  ESIMCoupling coupling(v);ComplexSolver solver(matrix);
  std::vector<cuDoubleComplex> drive(cells,make_cuDoubleComplex(100.,0.));
  for(int repetition=0;repetition<5;++repetition){
   coupling.prepare(solver.stream(),temperature.get<double>(),1.);
   check(cudaMemcpyAsync(solver.rhs(),drive.data(),cells*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice,solver.stream()));
   const auto solved=coupling.solve(solver,1e-10,200);coupling.commit();
   if(solved.trueRelativeResidual>1e-10||coupling.residual()>1e-8||coupling.fieldResidual()>1e-8)throw std::runtime_error("Constitutive block residual gate failed");
  }
  Buffer thermalValues,zeroValues,heatValues,constantValues,driveValues;
  thermalValues.put(std::vector<cuDoubleComplex>(cells,make_cuDoubleComplex(1.,0.)));
  zeroValues.put(std::vector<double>(cells,0.));heatValues.put(std::vector<double>(cells,.001));
  constantValues.put(std::vector<double>(cells,1.));driveValues.put(drive);
  CoupledOperators op{};op.eddy=matrix;op.lossProjection=v.fieldProjection;
  op.thermalSystem={cells,diagonalRows.get<int32_t>(),diagonalColumns.get<int32_t>(),thermalValues.get<cuDoubleComplex>()};
  op.stiffness=op.convection={cells,diagonalRows.get<int32_t>(),diagonalColumns.get<int32_t>(),zeroValues.get<double>()};
  op.heatProjection={cells,diagonalRows.get<int32_t>(),diagonalColumns.get<int32_t>(),heatValues.get<double>()};
  op.rhsPerAmp=driveValues.get<cuDoubleComplex>();op.lossWeights=weights.get<double>();op.constantCoefficients=constantValues.get<double>();op.initialTemperature=temperature.get<double>();
  op.heatSamples=cells;op.lossComponents=1;op.dt=.001;op.conductivity=1.;op.convectionScale=0.;op.dynamicMaterial=true;op.materialOwner=&coupling;
  op.updateMaterial=[](void* owner,cudaStream_t stream,const double* t,double current){static_cast<ESIMCoupling*>(owner)->prepare(stream,t,current);};
  op.solveMaterial=[](void* owner,ComplexSolver& linear,double tolerance,int iterations){return static_cast<ESIMCoupling*>(owner)->solve(linear,tolerance,iterations);};
  op.commitMaterial=[](void* owner){static_cast<ESIMCoupling*>(owner)->commit();};op.resetMaterial=[](void* owner){static_cast<ESIMCoupling*>(owner)->reset();};
  CoupledState state(op);for(int step=0;step<5;++step)state.update(1.,300.);
  const auto accepted=state.electromagneticState();const auto acceptedTemperature=state.temperature();const auto acceptedHeat=state.heat();
  bool rejected=false;try{state.update(1.7976931348623157e308,300.);}catch(const std::exception&){rejected=true;}
  const auto after=state.electromagneticState();
  if(!rejected||state.steps()!=5||acceptedTemperature!=state.temperature()||acceptedHeat!=state.heat())throw std::runtime_error("Candidate failure changed accepted coupled state");
  for(size_t i=0;i<cells;++i)if(cuCreal(after[i])!=cuCreal(accepted[i])||cuCimag(after[i])!=cuCimag(accepted[i]))throw std::runtime_error("Candidate failure changed accepted electromagnetic state");
  out[0]=mxCreateDoubleMatrix(1,7,mxREAL);auto* values=mxGetDoubles(out[0]);
  values[0]=coupling.residual();values[1]=coupling.fieldResidual();values[2]=coupling.iterations();values[3]=coupling.ownedBytes()+state.ownedBytes();
  values[4]=state.eddyResult().trueRelativeResidual;values[5]=state.thermalResult().trueRelativeResidual;values[6]=state.steps();
  state.reset();for(auto z:state.electromagneticState())if(cuCreal(z)!=0.||cuCimag(z)!=0.)throw std::runtime_error("Reset retained electromagnetic state");
 }catch(const std::exception& e){mexErrMsgIdAndTxt("radia:ih:couplingSmoke","%s",e.what());}
}
