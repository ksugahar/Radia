// Validation-only gateway for resident algebra; not an IH solver.
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include "ih_resident_kernels.cuh"
#include "ih_bicgstab_kernels.cuh"
#include "ih_complex_solver.cuh"
#include "ih_thermal_kernels.cuh"
#include "ih_coupled_state.cuh"
#include <vector>
#include <stdexcept>
#include <cmath>
using namespace radia_ih_cuda;
void check(cudaError_t e){if(e!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(e));}
struct Buffer {
    void* p=nullptr;
    ~Buffer(){if(p)cudaFree(p);}
    template<class T> void put(const std::vector<T>& a){check(cudaMalloc(&p,a.size()*sizeof(T)));check(cudaMemcpy(p,a.data(),a.size()*sizeof(T),cudaMemcpyHostToDevice));}
    template<class T> T* get(){return static_cast<T*>(p);}
};
std::vector<cuDoubleComplex> complexValues(const mxArray* a,size_t count){
    if(!mxIsDouble(a)||mxIsSparse(a)||mxGetNumberOfElements(a)!=count)throw std::runtime_error("Invalid complex array");
    std::vector<cuDoubleComplex> result(count);
    const auto* z=mxIsComplex(a)?mxGetComplexDoubles(a):nullptr;
    const auto* r=mxIsComplex(a)?nullptr:mxGetDoubles(a);
    for(size_t i=0;i<count;++i){double re=z?z[i].real:r[i],im=z?z[i].imag:0.;if(!std::isfinite(re)||!std::isfinite(im))throw std::runtime_error("Nonfinite input");result[i]=make_cuDoubleComplex(re,im);}
    return result;
}
void mexFunction(int nlhs,mxArray** out,int nrhs,const mxArray** in){
 try{
    mxInitGPU();
    if(nrhs!=4||nlhs!=4)throw std::runtime_error("Expected A(3x3),P(6x3),x(3),current scalar and four outputs");
    if(mxGetM(in[0])!=3||mxGetN(in[0])!=3||mxGetM(in[1])!=6||mxGetN(in[1])!=3)throw std::runtime_error("Invalid matrix shape");
    if(!mxIsDouble(in[3])||mxIsComplex(in[3])||mxGetNumberOfElements(in[3])!=1||!std::isfinite(mxGetScalar(in[3])))throw std::runtime_error("Invalid drive");
    auto ac=complexValues(in[0],9),pc=complexValues(in[1],18),x=complexValues(in[2],3);
    std::vector<cuDoubleComplex> a(9),p(18);
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)a[3*i+j]=ac[i+3*j];
    for(int i=0;i<6;++i)for(int j=0;j<3;++j)p[3*i+j]=pc[i+6*j];
    std::vector<int32_t> ar={0,3,6,9},pr={0,3,6,9,12,15,18},cols={0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2};
    Buffer da,dp,dx,drhs,dy,dq,dar,dpr,dc,di,dw;
    da.put(a);dp.put(p);dx.put(x);dar.put(ar);dpr.put(pr);dc.put(cols);
    di.put(std::vector<double>{mxGetScalar(in[3])});dw.put(std::vector<double>(6,.5));
    drhs.put(std::vector<cuDoubleComplex>(3));dy.put(std::vector<cuDoubleComplex>(3));dq.put(std::vector<double>(2));
    ComplexCSRView av{3,dar.get<int32_t>(),dc.get<int32_t>(),da.get<cuDoubleComplex>()};
    ComplexCSRView pv{6,dpr.get<int32_t>(),dc.get<int32_t>(),dp.get<cuDoubleComplex>()};
    // No host vector exchange or explicit synchronization between stages.
    for(int repeat=0;repeat<10;++repeat){
        scaleDrive<<<1,32>>>(dx.get<cuDoubleComplex>(),di.get<double>(),3,drhs.get<cuDoubleComplex>());
        complexMatvec<<<1,32>>>(av,drhs.get<cuDoubleComplex>(),dy.get<cuDoubleComplex>());
        quadraticLoss<<<1,32>>>(pv,dy.get<cuDoubleComplex>(),dw.get<double>(),2,3,dq.get<double>());
    }
    check(cudaGetLastError());check(cudaDeviceSynchronize());
    ComplexSolver eddy(av);
    check(cudaMemcpyAsync(eddy.rhs(),dy.p,3*sizeof(cuDoubleComplex),cudaMemcpyDeviceToDevice,eddy.stream()));
    auto eddyResult=eddy.solve(1e-12,200);
    quadraticLoss<<<1,32,0,eddy.stream()>>>(pv,eddy.solution(),dw.get<double>(),2,3,dq.get<double>());
    check(cudaStreamSynchronize(eddy.stream()));
    Buffer dk,dcv,dmap,dmr,dmc,dload,dtemp,done,dambient,dsystem,dcandidate,dbad;
    std::vector<double> k={2.,-1.,0.,-1.,2.,-1.,0.,-1.,2.},cv(9,0.),sys(9);
    std::vector<cuDoubleComplex> system(9);
    for(int i=0;i<3;++i)cv[3*i+i]=.1;
    for(int i=0;i<9;++i)system[i]=make_cuDoubleComplex((i%4==0?1.:0.)+.01*k[i]+.02*cv[i],0.);
    dk.put(k);dcv.put(cv);dsystem.put(system);dmap.put(std::vector<double>{1.,0.,.5,.5,0.,1.});
    dmr.put(std::vector<int32_t>{0,2,4,6});dmc.put(std::vector<int32_t>{0,1,0,1,0,1});
    dload.put(std::vector<double>(3));dtemp.put(std::vector<double>{300.,301.,302.});
    done.put(std::vector<double>(3,1.));dambient.put(std::vector<double>{293.15});
    dcandidate.put(std::vector<double>(3));dbad.put(std::vector<int>{0});
    ComplexCSRView thermalMatrix{3,dar.get<int32_t>(),dc.get<int32_t>(),dsystem.get<cuDoubleComplex>()};
    ComplexSolver thermal(thermalMatrix);
    RealCSRView kv{3,dar.get<int32_t>(),dc.get<int32_t>(),dk.get<double>()};
    RealCSRView cvv{3,dar.get<int32_t>(),dc.get<int32_t>(),dcv.get<double>()};
    RealCSRView loadMap{3,dmr.get<int32_t>(),dmc.get<int32_t>(),dmap.get<double>()};
    thermalLoad<<<1,32,0,thermal.stream()>>>(loadMap,dq.get<double>(),dload.get<double>());
    thermalIncrementRHS<<<1,32,0,thermal.stream()>>>(kv,cvv,dload.get<double>(),dtemp.get<double>(),done.get<double>(),dambient.get<double>(),.01,1.,2.,thermal.rhs());
    auto thermalResult=thermal.solve(1e-12,200);
    thermalCandidate<<<1,32,0,thermal.stream()>>>(dtemp.get<double>(),thermal.solution(),dcandidate.get<double>(),dbad.get<int>(),3);
    int invalid=0;check(cudaMemcpyAsync(&invalid,dbad.p,sizeof(int),cudaMemcpyDeviceToHost,thermal.stream()));
    check(cudaStreamSynchronize(thermal.stream()));if(invalid)throw std::runtime_error("Invalid thermal candidate");
    out[0]=mxCreateDoubleMatrix(3,1,mxCOMPLEX);out[1]=mxCreateDoubleMatrix(2,1,mxREAL);out[3]=mxCreateDoubleMatrix(3,1,mxREAL);
    std::vector<cuDoubleComplex> y(3);check(cudaMemcpy(y.data(),eddy.solution(),3*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost));
    check(cudaMemcpy(mxGetDoubles(out[1]),dq.p,2*sizeof(double),cudaMemcpyDeviceToHost));
    check(cudaMemcpy(mxGetDoubles(out[3]),dcandidate.p,3*sizeof(double),cudaMemcpyDeviceToHost));
    auto* result=mxGetComplexDoubles(out[0]);for(int i=0;i<3;++i){result[i].real=cuCreal(y[i]);result[i].imag=cuCimag(y[i]);}
    // Identity-system first-correction convergence: omega is deliberately
    // not formed when s=0. This checks the resident BiCGStab stage contract.
    double norm=0.;for(const auto z:x)norm+=cuCreal(z)*cuCreal(z)+cuCimag(z)*cuCimag(z);
    BiCGState state{};state.rho=state.denominator=make_cuDoubleComplex(norm,0.);
    state.rhsNormSquared=norm;state.relativeToleranceSquared=1e-24;
    Buffer ds,dsolution,dscratch,dr,dt,drr;
    ds.put(std::vector<BiCGState>{state});dsolution.put(std::vector<cuDoubleComplex>(3));
    dscratch.put(std::vector<cuDoubleComplex>(3));dr.put(x);dt.put(std::vector<cuDoubleComplex>(3));
    drr.put(std::vector<cuDoubleComplex>{make_cuDoubleComplex(0.,0.)});
    bicgAlpha<<<1,1>>>(ds.get<BiCGState>());
    bicgIntermediate<<<1,32>>>(ds.get<BiCGState>(),dx.get<cuDoubleComplex>(),dx.get<cuDoubleComplex>(),dscratch.get<cuDoubleComplex>(),3);
    bicgCheckIntermediate<<<1,1>>>(ds.get<BiCGState>(),drr.get<cuDoubleComplex>());
    bicgOmega<<<1,1>>>(ds.get<BiCGState>());
    bicgUpdate<<<1,32>>>(ds.get<BiCGState>(),dx.get<cuDoubleComplex>(),dscratch.get<cuDoubleComplex>(),dt.get<cuDoubleComplex>(),dsolution.get<cuDoubleComplex>(),dr.get<cuDoubleComplex>(),3);
    bicgFinish<<<1,1>>>(ds.get<BiCGState>(),drr.get<cuDoubleComplex>());
    check(cudaGetLastError());check(cudaDeviceSynchronize());
    std::vector<cuDoubleComplex> solution(3);
    check(cudaMemcpy(solution.data(),dsolution.p,3*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost));
    check(cudaMemcpy(&state,ds.p,sizeof(state),cudaMemcpyDeviceToHost));
    if(state.status!=1||state.iteration!=1)throw std::runtime_error("BiCGStab first-correction transition failed");
    for(int i=0;i<3;++i)if(cuCabs(cuCsub(solution[i],x[i]))>1e-12)throw std::runtime_error("BiCGStab identity correction failed");
    // Persistent coupled state: operators remain uploaded across all updates.
    // Use A*x as the per-amp RHS, matching the independent CPU algebra oracle.
    complexMatvec<<<1,32>>>(av,dx.get<cuDoubleComplex>(),drhs.get<cuDoubleComplex>());
    check(cudaDeviceSynchronize());
    CoupledOperators coupled{};
    coupled.eddy=av;coupled.lossProjection=pv;coupled.thermalSystem=thermalMatrix;
    coupled.stiffness=kv;coupled.convection=cvv;coupled.heatProjection=loadMap;
    coupled.rhsPerAmp=drhs.get<cuDoubleComplex>();coupled.lossWeights=dw.get<double>();
    coupled.constantCoefficients=done.get<double>();coupled.initialTemperature=dtemp.get<double>();
    coupled.heatSamples=2;coupled.lossComponents=3;coupled.dt=.01;coupled.conductivity=1.;coupled.convectionScale=2.;
    // A missing dynamic material implementation must never freeze ESIM silently.
    coupled.dynamicMaterial=true;
    bool rejected=false;try{CoupledState unsupported(coupled);}catch(const std::invalid_argument&){rejected=true;}
    if(!rejected)throw std::runtime_error("Missing dynamic material update accepted");
    coupled.dynamicMaterial=false;
    CoupledState persistent(coupled);
    for(int step=0;step<10;++step)persistent.update(mxGetScalar(in[3]),293.15,1e-12,200);
    if(persistent.steps()!=10)throw std::runtime_error("Resident update count mismatch");
    auto trajectory=persistent.temperature();auto acceptedHeat=persistent.heat();
    bool invalidControlRejected=false;
    try{persistent.update(std::numeric_limits<double>::quiet_NaN(),293.15);}catch(const std::invalid_argument&){invalidControlRejected=true;}
    if(!invalidControlRejected||persistent.steps()!=10||persistent.temperature()!=trajectory||persistent.heat()!=acceptedHeat)
        throw std::runtime_error("Invalid update changed accepted resident state");
    persistent.reset();
    const auto resetTemperature=persistent.temperature();const auto resetHeat=persistent.heat();
    if(persistent.steps()!=0||resetTemperature!=std::vector<double>({300.,301.,302.})||resetHeat!=std::vector<double>({0.,0.}))
        throw std::runtime_error("Resident reset mismatch");
    const char* fields[]={"resident_repeats","eddy_true_residual","thermal_true_residual","resident_temperature","resident_steps","reset_checked","invalid_update_checked","missing_material_rejected"};
    out[2]=mxCreateStructMatrix(1,1,8,fields);
    mxSetField(out[2],0,fields[0],mxCreateDoubleScalar(10));
    mxSetField(out[2],0,fields[1],mxCreateDoubleScalar(eddyResult.trueRelativeResidual));
    mxSetField(out[2],0,fields[2],mxCreateDoubleScalar(thermalResult.trueRelativeResidual));
    auto* residentTemperature=mxCreateDoubleMatrix(3,1,mxREAL);
    std::copy(trajectory.begin(),trajectory.end(),mxGetDoubles(residentTemperature));
    mxSetField(out[2],0,fields[3],residentTemperature);
    mxSetField(out[2],0,fields[4],mxCreateDoubleScalar(10));
    for(int f=5;f<8;++f)mxSetField(out[2],0,fields[f],mxCreateLogicalScalar(true));
 }catch(const std::exception& e){mexErrMsgIdAndTxt("radia:ih:cudaSmoke","%s",e.what());}
}
