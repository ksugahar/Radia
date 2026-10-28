#pragma once
// Shared domain implementation: no MATLAB or Python runtime dependencies.
#include "csr_pcg_host.hpp"
#include <cuda_runtime.h>
namespace radia_csr {
constexpr int threads=256;
inline void check(cudaError_t code) {
    if(code!=cudaSuccess) throw std::runtime_error(cudaGetErrorString(code));
}
template<class T> struct DeviceBuffer {
    T* data=nullptr;
    explicit DeviceBuffer(size_t count) { check(cudaMalloc(reinterpret_cast<void**>(&data),count*sizeof(T))); }
    DeviceBuffer(const DeviceBuffer&)=delete;
    DeviceBuffer& operator=(const DeviceBuffer&)=delete;
    ~DeviceBuffer(){if(data) cudaFree(data);}
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
inline Result solve_device(int n,int nnz,const int* row,const int* col,
    const double* val,const double* rhs,const double* dinv,double tol,int maxiter,double* x) {
    validate_dimensions(n,nnz,tol,maxiter);
        DeviceBuffer<double> r(n), p(n), q(n), state(9);
        check(cudaMemset(state.data,0,9*sizeof(double)));
        check(cudaMemset(x,0,n*sizeof(double)));
        validate<<<(size_t(std::max(n,nnz))+threads-1)/threads,threads>>>((int)n,(int)nnz,row,col,val,rhs,dinv,state.data);
        double host[9]; check(cudaMemcpy(host,state.data,9*sizeof(double),cudaMemcpyDeviceToHost));
        if(host[6]) throw std::runtime_error("Invalid CSR indices or nonfinite values / nonpositive diagonal inverse");
        int blocks=(int)((size_t(n)+threads-1)/threads);
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

    return {static_cast<int>(host[7]),relative};
}
} // namespace radia_csr
