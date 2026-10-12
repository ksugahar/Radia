// Validation-only gateway for resident algebra; not an IH solver.
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include "ih_resident_kernels.cuh"
#include "ih_bicgstab_kernels.cuh"
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
    if(nrhs!=4||nlhs!=3)throw std::runtime_error("Expected A(3x3),P(6x3),x(3),current scalar and three outputs");
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
    out[0]=mxCreateDoubleMatrix(3,1,mxCOMPLEX);out[1]=mxCreateDoubleMatrix(2,1,mxREAL);
    std::vector<cuDoubleComplex> y(3);check(cudaMemcpy(y.data(),dy.p,3*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost));
    check(cudaMemcpy(mxGetDoubles(out[1]),dq.p,2*sizeof(double),cudaMemcpyDeviceToHost));
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
    out[2]=mxCreateDoubleScalar(10);
 }catch(const std::exception& e){mexErrMsgIdAndTxt("radia:ih:cudaSmoke","%s",e.what());}
}
