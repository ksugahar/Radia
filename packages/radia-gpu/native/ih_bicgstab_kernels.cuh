// Device-scalar BiCGStab stages for the IH complex operator.
// cuBLAS reductions write directly to these resident scalars. Host orchestration
// and true-residual certification are mandatory; these kernels alone are not a solver.
#pragma once
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cmath>
#include <cstddef>
namespace radia_ih_cuda {
struct BiCGState {
    cuDoubleComplex rho, previousRho, alpha, omega, denominator, ts, tt;
    double rhsNormSquared, relativeToleranceSquared;
    unsigned iteration;
    int status; // 0=active, 1=recursive convergence, -1=breakdown/nonfinite
};
__device__ inline bool usable(cuDoubleComplex z) {
    return isfinite(cuCreal(z))&&isfinite(cuCimag(z)) &&
           (cuCreal(z)!=0. || cuCimag(z)!=0.);
}
__global__ void bicgDirection(BiCGState* scalar,const cuDoubleComplex* r,
        const cuDoubleComplex* v,cuDoubleComplex* p,std::size_t n) {
    const std::size_t i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(i>=n||scalar->status!=0) return;
    if(scalar->iteration==0) {p[i]=r[i];return;}
    const cuDoubleComplex beta=cuCmul(cuCdiv(scalar->rho,scalar->previousRho),
                                     cuCdiv(scalar->alpha,scalar->omega));
    p[i]=cuCadd(r[i],cuCmul(beta,cuCsub(p[i],cuCmul(scalar->omega,v[i]))));
}
__global__ void bicgAlpha(BiCGState* scalar) {
    if(scalar->status!=0)return;
    if(!usable(scalar->rho)||!usable(scalar->denominator)) {scalar->status=-1;return;}
    scalar->alpha=cuCdiv(scalar->rho,scalar->denominator);
    if(!usable(scalar->alpha))scalar->status=-1;
}
__global__ void bicgIntermediate(const BiCGState* scalar,
        const cuDoubleComplex* r,const cuDoubleComplex* v,cuDoubleComplex* s,std::size_t n) {
    const std::size_t i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(i<n&&scalar->status==0)s[i]=cuCsub(r[i],cuCmul(scalar->alpha,v[i]));
}
// rr is a device reduction of s^H*s. An exact/near-exact first correction
// must be accepted before computing omega, whose denominator would be zero.
__global__ void bicgCheckIntermediate(BiCGState* scalar,const cuDoubleComplex* rr) {
    if(scalar->status!=0)return;
    const double norm=cuCreal(*rr);
    if(!isfinite(norm)||norm<0.)scalar->status=-1;
    else if(norm<=scalar->relativeToleranceSquared*scalar->rhsNormSquared)scalar->status=2;
}
__global__ void bicgOmega(BiCGState* scalar) {
    if(scalar->status!=0)return;
    if(!usable(scalar->tt)||!usable(scalar->ts)) {scalar->status=-1;return;}
    scalar->omega=cuCdiv(scalar->ts,scalar->tt);
    if(!usable(scalar->omega))scalar->status=-1;
}
__global__ void bicgUpdate(const BiCGState* scalar,const cuDoubleComplex* p,
        const cuDoubleComplex* s,const cuDoubleComplex* t,
        cuDoubleComplex* x,cuDoubleComplex* r,std::size_t n) {
    const std::size_t i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(i>=n||scalar->status<0||scalar->status==1)return;
    if(scalar->status==2) {
        x[i]=cuCadd(x[i],cuCmul(scalar->alpha,p[i]));r[i]=s[i];return;
    }
    x[i]=cuCadd(x[i],cuCadd(cuCmul(scalar->alpha,p[i]),cuCmul(scalar->omega,s[i])));
    r[i]=cuCsub(s[i],cuCmul(scalar->omega,t[i]));
}
__global__ void bicgFinish(BiCGState* scalar,const cuDoubleComplex* rr) {
    if(scalar->status<0||scalar->status==1)return;
    const double norm=cuCreal(*rr);
    if(!isfinite(norm)||norm<0.) {scalar->status=-1;return;}
    scalar->previousRho=scalar->rho;++scalar->iteration;
    if(scalar->status==2||norm<=scalar->relativeToleranceSquared*scalar->rhsNormSquared)scalar->status=1;
}
__global__ void trueResidual(const cuDoubleComplex* rhs,const cuDoubleComplex* ax,
        cuDoubleComplex* residual,std::size_t n) {
    const std::size_t i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(i<n)residual[i]=cuCsub(rhs[i],ax[i]);
}
} // namespace radia_ih_cuda
