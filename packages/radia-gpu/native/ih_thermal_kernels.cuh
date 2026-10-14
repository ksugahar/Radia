// Backward-Euler increment form shared with the CPU IH physical contract.
#pragma once
#include "ih_resident_kernels.cuh"
namespace radia_ih_cuda {
struct RealCSRView {std::size_t rows;const std::int32_t* rowptr;const std::int32_t* columns;const double* values;};
__global__ void thermalLoad(RealCSRView projection,const double* heat,double* load) {
    const std::size_t i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(i>=projection.rows)return;
    double value=0.;for(std::int32_t p=projection.rowptr[i];p<projection.rowptr[i+1];++p)value+=projection.values[p]*heat[projection.columns[p]];
    load[i]=value;
}
// (M+dt*k*K+dt*h*C)*dT = dt*(Q-k*K*T+h*C*(Tamb*one-T)).
// Q is the assembled load, not heat density multiplied by weights again.
__global__ void thermalIncrementRHS(RealCSRView stiffness,RealCSRView convection,
        const double* load,const double* temperature,const double* one,
        const double* ambient,double dt,double conductivity,double h,
        cuDoubleComplex* rhs) {
    const std::size_t i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(i>=stiffness.rows)return;
    double value=load[i];
    for(std::int32_t p=stiffness.rowptr[i];p<stiffness.rowptr[i+1];++p)
        value-=conductivity*stiffness.values[p]*temperature[stiffness.columns[p]];
    if(convection.values)for(std::int32_t p=convection.rowptr[i];p<convection.rowptr[i+1];++p){
        const auto j=convection.columns[p];value+=h*convection.values[p]*(ambient[0]*one[j]-temperature[j]);
    }
    rhs[i]=make_cuDoubleComplex(dt*value,0.);
}
// Write a candidate state; the owner commits it only after solve/residual and
// finite/imaginary-state checks succeed. Failed updates preserve old state.
__global__ void thermalCandidate(const double* accepted,const cuDoubleComplex* increment,
        double* candidate,int* invalid,std::size_t n) {
    const std::size_t i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(i>=n)return;
    const double value=accepted[i]+cuCreal(increment[i]);
    if(!isfinite(value)||fabs(cuCimag(increment[i]))>1e-12*fmax(1.,fabs(cuCreal(increment[i]))))atomicExch(invalid,1);
    candidate[i]=value;
}
} // namespace radia_ih_cuda
