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
// NGSolve/caller-owned basis and gradient evaluation, row-major on device.
// No FE geometry, quadrature generation or basis convention is reimplemented.
struct ThermalQuadratureView {
    std::size_t n,q,b;int dimension,nt;
    const double *axis,*basis,*gradient,*volume,*density,*cp,*kappa;
    const double *boundaryBasis,*area,*convection,*emissivity;
};
__device__ inline bool thermalTable(const double* axis,const double* table,int nt,
        std::size_t row,double t,double& value,double& slope,double* integral=nullptr) {
    if(!isfinite(t)||t<axis[0]||t>axis[nt-1])return false;
    int lo=0,hi=nt-1;while(hi-lo>1){const int m=(lo+hi)/2;if(t<axis[m])hi=m;else lo=m;}
    const auto* v=table+row*nt;const double h=axis[hi]-axis[lo];
    slope=(v[hi]-v[lo])/h;value=v[lo]+slope*(t-axis[lo]);
    if(integral){double sum=0.;for(int j=0;j<lo;++j)sum+=.5*(v[j]+v[j+1])*(axis[j+1]-axis[j]);
        *integral=sum+.5*(v[lo]+value)*(t-axis[lo]);}
    return isfinite(value)&&isfinite(slope)&&(!integral||isfinite(*integral));
}
// Volume samples: [storage W, volumetric capacity/dt, k, dk/dT, grad(T)...].
__global__ void thermalEvaluateQuadrature(ThermalQuadratureView v,const double* x,
        const double* old,double ambient,double dt,double* volume,double* boundary,int* invalid) {
    const auto q=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(q<v.q){double t=0.,before=0.;for(std::size_t j=0;j<v.n;++j){t+=v.basis[q*v.n+j]*x[j];before+=v.basis[q*v.n+j]*old[j];}
        double cp,dcp,H,oldcp,oldslope,oldH,k,dk;
        if(!thermalTable(v.axis,v.cp,v.nt,q,t,cp,dcp,&H)||!thermalTable(v.axis,v.cp,v.nt,q,before,oldcp,oldslope,&oldH)||
           !thermalTable(v.axis,v.kappa,v.nt,q,t,k,dk)||cp<=0.||k<=0.){atomicExch(invalid,1);return;}
        auto* a=volume+q*(4+v.dimension);a[0]=v.volume[q]*v.density[q]*(H-oldH)/dt;a[1]=v.density[q]*cp/dt;a[2]=k;a[3]=dk;
        for(int d=0;d<v.dimension;++d){double g=0.;for(std::size_t j=0;j<v.n;++j)g+=v.gradient[(q*v.dimension+d)*v.n+j]*x[j];a[4+d]=g;}
    }
    if(q<v.b){double t=0.;for(std::size_t j=0;j<v.n;++j)t+=v.boundaryBasis[q*v.n+j]*x[j];double eps,deps;
        if(!thermalTable(v.axis,v.emissivity,v.nt,q,t,eps,deps)||eps<0.||eps>1.){atomicExch(invalid,1);return;}
        constexpr double sigma=5.670374419e-8;const double fourth=t*t*t*t-ambient*ambient*ambient*ambient;
        boundary[2*q]=v.convection[q]*(t-ambient)+sigma*eps*fourth;
        boundary[2*q+1]=v.convection[q]+sigma*(deps*fourth+4.*eps*t*t*t);
        if(!isfinite(boundary[2*q])||!isfinite(boundary[2*q+1]))atomicExch(invalid,1);
    }
}
// Bounded dense Newton matrix. The complete supplied FE maps are contracted;
// diagonal mass or nodal-temperature substitutions are not made.
__global__ void thermalAssembleQuadrature(ThermalQuadratureView v,const double* volume,
        const double* boundary,const double* load,cuDoubleComplex* matrix,cuDoubleComplex* rhs,int* invalid) {
    const auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=v.n)return;
    double residual=-load[i];
    for(std::size_t j=0;j<v.n;++j)matrix[i*v.n+j]=make_cuDoubleComplex(0.,0.);
    for(std::size_t q=0;q<v.q;++q){const auto* a=volume+q*(4+v.dimension);const double phi=v.basis[q*v.n+i];double grad=0.;
        for(int d=0;d<v.dimension;++d)grad+=v.gradient[(q*v.dimension+d)*v.n+i]*a[4+d];
        residual+=phi*a[0]+v.volume[q]*a[2]*grad;
        for(std::size_t j=0;j<v.n;++j){double gg=0.;for(int d=0;d<v.dimension;++d)gg+=v.gradient[(q*v.dimension+d)*v.n+i]*v.gradient[(q*v.dimension+d)*v.n+j];
            matrix[i*v.n+j].x+=v.volume[q]*(a[1]*phi*v.basis[q*v.n+j]+a[2]*gg+a[3]*grad*v.basis[q*v.n+j]);}
    }
    for(std::size_t q=0;q<v.b;++q){const double phi=v.boundaryBasis[q*v.n+i];residual+=v.area[q]*phi*boundary[2*q];
        for(std::size_t j=0;j<v.n;++j)matrix[i*v.n+j].x+=v.area[q]*phi*boundary[2*q+1]*v.boundaryBasis[q*v.n+j];}
    rhs[i]=make_cuDoubleComplex(-residual,0.);if(!isfinite(residual))atomicExch(invalid,1);
    for(std::size_t j=0;j<v.n;++j)if(!isfinite(matrix[i*v.n+j].x))atomicExch(invalid,1);
}
__global__ void thermalQuadratureBalance(ThermalQuadratureView v,const double* volume,
        const double* boundary,const double* load,const double* constant,const cuDoubleComplex* rhs,
        double floor,double* result,int* invalid) {
    if(blockIdx.x||threadIdx.x)return;double input=0.,storage=0.,loss=0.,rr=0.,ff=0.;
    for(std::size_t i=0;i<v.n;++i){input+=constant[i]*load[i];rr+=cuCreal(rhs[i])*cuCreal(rhs[i]);ff+=load[i]*load[i];}
    for(std::size_t q=0;q<v.q;++q)storage+=volume[q*(4+v.dimension)];
    for(std::size_t q=0;q<v.b;++q)loss+=v.area[q]*boundary[2*q];
    const double scale=fmax(fmax(fabs(input),fabs(storage)+fabs(loss)),floor);
    result[0]=sqrt(rr)/fmax(sqrt(ff),floor);result[1]=fabs(input-storage-loss)/scale;
    result[2]=input;result[3]=storage;result[4]=loss;result[5]=input-storage-loss;
    for(int i=0;i<6;++i)if(!isfinite(result[i]))atomicExch(invalid,1);
}
} // namespace radia_ih_cuda
