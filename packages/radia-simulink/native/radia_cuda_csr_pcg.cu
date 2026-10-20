#include "mex_error_boundary.hpp"
// CUDA-only SPD CSR solve. NGSolve owns finite-element assembly and numbering.
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <cuda_runtime.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr int threads = 256;
void check(cudaError_t code) {
    if (code != cudaSuccess) throw std::runtime_error(cudaGetErrorString(code));
}
struct Inputs {
    std::vector<const mxGPUArray*> arrays;
    ~Inputs() { for (auto a : arrays) mxGPUDestroyGPUArray(a); }
    const mxGPUArray* get(const mxArray* arg, mxClassID type) {
        if (!mxIsGPUArray(arg)) throw std::runtime_error("All system arrays must be gpuArray; no CPU fallback");
        auto a = mxGPUCreateFromMxArray(arg); arrays.push_back(a);
        if (mxGPUGetClassID(a) != type || mxGPUGetComplexity(a) != mxREAL)
            throw std::runtime_error("Expected real int32 CSR indices and real double numerical arrays");
        return a;
    }
};
struct Buffer {
    double* data = nullptr;
    explicit Buffer(size_t n) { check(cudaMalloc(&data, n*sizeof(double))); }
    ~Buffer() { if (data) cudaFree(data); }
};
struct GPUOutput {
    mxGPUArray* array;
    explicit GPUOutput(mwSize n) : array(mxGPUCreateGPUArray(1, &n, mxDOUBLE_CLASS, mxREAL, MX_GPU_INITIALIZE_VALUES)) {}
    ~GPUOutput() { mxGPUDestroyGPUArray(array); }
};
__device__ double reduce(double v) {
    __shared__ double s[threads]; s[threadIdx.x] = v; __syncthreads();
    for (int offset=threads/2; offset; offset/=2) {
        if (threadIdx.x < offset) s[threadIdx.x] += s[threadIdx.x+offset];
        __syncthreads();
    }
    return s[0];
}
// state: rho, pAp, rho_new, rr, bnorm2, done, breakdown, iterations, nonzero_rhs.
__global__ void validate(int n, int nnz, const int* rp, const int* ci,
                         const double* av, const double* b, const double* di, double* state) {
    size_t i = size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if (i<n && (rp[i]<0 || rp[i+1]<rp[i] || rp[i+1]>nnz ||
                !isfinite(b[i]) || !isfinite(di[i]) || di[i]<=0)) atomicAdd(state+6,1.0);
    if (i<n && b[i]!=0) atomicAdd(state+8,1.0);
    if (i<nnz && (ci[i]<0 || ci[i]>=n || !isfinite(av[i]))) atomicAdd(state+6,1.0);
    if (i==0 && (rp[0]!=0 || rp[n]!=nnz)) atomicAdd(state+6,1.0);
}
__global__ void initialize(int n, const double* b, const double* di,
                           double* r, double* p, double* state) {
    size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    double rho=0, rr=0;
    if (i<n) { r[i]=b[i]; p[i]=di[i]*b[i]; rho=b[i]*p[i]; rr=b[i]*b[i]; }
    double sum=reduce(rho); if(threadIdx.x==0) atomicAdd(state,sum);
    sum=reduce(rr); if(threadIdx.x==0) atomicAdd(state+4,sum);
}
__global__ void spmv(int n, const int* rp, const int* ci, const double* av,
                     const double* x, double* y, const double* state, bool guarded) {
    if (guarded && state[5]) return;
    size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if (i<n) { double sum=0; for(int j=rp[i];j<rp[i+1];++j) sum+=av[j]*x[ci[j]]; y[i]=sum; }
}
__global__ void dot_pap(int n, const double* p, const double* q, double* state) {
    if(state[5]) return;
    size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    double sum=reduce(i<n ? p[i]*q[i] : 0); if(threadIdx.x==0) atomicAdd(state+1,sum);
}
__global__ void update_xr(int n, const double* di, const double* p, const double* q,
                          double* x, double* r, double* state) {
    if(state[5]) return;
    if(!(state[1]>0) || !isfinite(state[1]) || !(state[0]>0)) {
        if(blockIdx.x==0 && threadIdx.x==0) { state[6]=1; state[5]=1; } return;
    }
    size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    double rn=0, rr=0, alpha=state[0]/state[1];
    if(i<n) { x[i]+=alpha*p[i]; r[i]-=alpha*q[i]; rr=r[i]*r[i]; rn=rr*di[i]; }
    double sum=reduce(rn); if(threadIdx.x==0) atomicAdd(state+2,sum);
    sum=reduce(rr); if(threadIdx.x==0) atomicAdd(state+3,sum);
}
__global__ void update_p(int n, const double* di, const double* r, double* p, const double* state) {
    if(state[5]) return;
    size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(i<n) p[i]=di[i]*r[i]+(state[2]/state[0])*p[i];
}
__global__ void advance(double* state, double tol2) {
    if(state[5]) return;
    state[7]+=1; state[0]=state[2];
    if(!isfinite(state[3]) || !isfinite(state[0])) { state[6]=1; state[5]=1; }
    else if(state[3]<=tol2*state[4]) state[5]=1;
    state[1]=state[2]=state[3]=0;
}
__global__ void residual(int n, const double* b, const double* ax, double* state) {
    size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    double r=i<n ? b[i]-ax[i] : 0;
    double scaled=state[4]>0 ? r/sqrt(state[4]) : r;
    double sum=reduce(scaled*scaled); if(threadIdx.x==0) atomicAdd(state+3,sum);
}
double scalar(const mxArray* a) {
    if(!mxIsDouble(a) || mxIsComplex(a) || mxGetNumberOfElements(a)!=1)
        throw std::runtime_error("Tolerance and iteration limit must be real double scalars");
    return mxGetScalar(a);
}
}

void dispatchCsr(int nlhs,mxArray** plhs,int nrhs,const mxArray** prhs){
        mxInitGPU();
        Inputs inputs;
        auto rp=inputs.get(prhs[0],mxINT32_CLASS), ci=inputs.get(prhs[1],mxINT32_CLASS);
        auto av=inputs.get(prhs[2],mxDOUBLE_CLASS), b=inputs.get(prhs[3],mxDOUBLE_CLASS);
        auto di=inputs.get(prhs[4],mxDOUBLE_CLASS);
        size_t n=mxGPUGetNumberOfElements(b), nnz=mxGPUGetNumberOfElements(av);
        double tol=scalar(prhs[5]), maxiter=scalar(prhs[6]);
        if(n==0 || n>INT32_MAX || nnz>INT32_MAX || mxGPUGetNumberOfElements(rp)!=n+1 ||
           mxGPUGetNumberOfElements(ci)!=nnz || mxGPUGetNumberOfElements(di)!=n ||
           !(tol>0 && tol<=1e-6) || !(maxiter>=1 && maxiter<=100000 && maxiter==floor(maxiter)))
            throw std::runtime_error("Invalid dimensions, tolerance (0,1e-6], or integer iteration limit [1,100000]");
        const int *row=static_cast<const int*>(mxGPUGetDataReadOnly(rp)), *col=static_cast<const int*>(mxGPUGetDataReadOnly(ci));
        const double *val=static_cast<const double*>(mxGPUGetDataReadOnly(av)), *rhs=static_cast<const double*>(mxGPUGetDataReadOnly(b)), *dinv=static_cast<const double*>(mxGPUGetDataReadOnly(di));
        GPUOutput output(n); auto x=static_cast<double*>(mxGPUGetData(output.array));
        Buffer r(n), p(n), q(n), state(9);
        check(cudaMemset(state.data,0,9*sizeof(double)));
        validate<<<(std::max(n,nnz)+threads-1)/threads,threads>>>((int)n,(int)nnz,row,col,val,rhs,dinv,state.data);
        double host[9]; check(cudaMemcpy(host,state.data,9*sizeof(double),cudaMemcpyDeviceToHost));
        if(host[6]) throw std::runtime_error("Invalid CSR indices or nonfinite values / nonpositive diagonal inverse");
        int blocks=(int)((n+threads-1)/threads);
        initialize<<<blocks,threads>>>((int)n,rhs,dinv,r.data,p.data,state.data);
        check(cudaMemcpy(host,state.data,9*sizeof(double),cudaMemcpyDeviceToHost));
        if(!std::isfinite(host[4]) || (host[4]==0 && host[8]>0))
            throw std::runtime_error("RHS norm is outside supported FP64 reduction range; rescale the system explicitly");
        if(host[4]>0) {
            for(int k=0;k<(int)maxiter;++k) {
                spmv<<<blocks,threads>>>((int)n,row,col,val,p.data,q.data,state.data,true);
                dot_pap<<<blocks,threads>>>((int)n,p.data,q.data,state.data);
                update_xr<<<blocks,threads>>>((int)n,dinv,p.data,q.data,x,r.data,state.data);
                update_p<<<blocks,threads>>>((int)n,dinv,r.data,p.data,state.data);
                advance<<<1,1>>>(state.data,tol*tol);
                if(k%10==9 || k+1==(int)maxiter) {
                    check(cudaMemcpy(host,state.data,9*sizeof(double),cudaMemcpyDeviceToHost));
                    if(host[5]) break;
                }
            }
        }
        spmv<<<blocks,threads>>>((int)n,row,col,val,x,q.data,state.data,false);
        // A breakdown may leave an accumulated recursive residual; always reset.
        check(cudaMemset(state.data+3,0,sizeof(double)));
        residual<<<blocks,threads>>>((int)n,rhs,q.data,state.data);
        check(cudaGetLastError());
        check(cudaMemcpy(host,state.data,9*sizeof(double),cudaMemcpyDeviceToHost));
        double relative=host[4]>0 ? std::sqrt(host[3]) : 0;
        if(host[6] || !std::isfinite(relative) || relative>1e-6 || relative>std::max(10*tol,1e-12))
            throw std::runtime_error("CG breakdown or true residual exceeds acceptance: "+std::to_string(relative));
        const char* names[]={"iterations","true_relative_residual","precision","backend","preconditioner"};
        plhs[1]=mxCreateStructMatrix(1,1,5,names);
        mxSetField(plhs[1],0,names[0],mxCreateDoubleScalar(host[7]));
        mxSetField(plhs[1],0,names[1],mxCreateDoubleScalar(relative));
        mxSetField(plhs[1],0,names[2],mxCreateString("float64"));
        mxSetField(plhs[1],0,names[3],mxCreateString("CUDA custom CSR kernels; NGSolve CPU assembly"));
        mxSetField(plhs[1],0,names[4],mxCreateString("Jacobi; not IC or AMS"));
        plhs[0]=mxGPUCreateMxArrayOnGPU(output.array);

}
void mexFunction(int nl,mxArray** out,int nr,const mxArray** in){
    if(nr!=7 || nl!=2) mexErrMsgIdAndTxt("radia:cuda:arity", "[x,info]=radia_cuda_csr_pcg(int32_rowptr,int32_columns,double_values,b,diagonal_inverse,tolerance,maxiter)");
    char errorMessage[2048]={0};
    if(!radia_mex::checked(dispatchCsr,nl,out,nr,in,errorMessage,sizeof(errorMessage))){
        mexErrMsgIdAndTxt("radia:cuda:solve","%s",errorMessage);
    }
}
