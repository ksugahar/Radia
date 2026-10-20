// Finite slab/cylinder ESIM cells: same uniform-grid finite differences,
// Picard relaxation and midpoint power integration as ESIMFiniteSlabSolver.
// This is a finite-difference cell, not an NGSolve FE replacement.
#pragma once
#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cmath>
namespace radia_ih_cuda {
struct ESIMLaw {
    const double* knots; // H in A/m, starts at zero
    const double* coefficients; // interval-major cubic B(H), descending powers
    int count; // zero: constant permeability, otherwise at least two
    cuDoubleComplex initialMu;
    double saturationMagnetization; // constant M beyond final knot
};
struct ESIMCellParameters {
    int nodes,maxIterations;
    bool cylinder;
    double halfThickness,resistivity,omega,tolerance,relaxation;
    ESIMLaw law;
};
struct ESIMCellResult {
    cuDoubleComplex impedance;
    double activePower,reactivePower,magneticPower,relativeChange,trueResidual;
    int iterations,status; // 1 converged; -1 invalid; -2 breakdown; -3 no convergence
};
__device__ inline bool cellFinite(cuDoubleComplex z){return isfinite(cuCreal(z))&&isfinite(cuCimag(z));}
__device__ inline cuDoubleComplex cellMu(const ESIMLaw& law,double h) {
    if(!law.count||h<1e-10)return law.initialMu;
    constexpr double mu0=1.2566370614359172954e-6;
    if(h>law.knots[law.count-1])return make_cuDoubleComplex(mu0*(h+law.saturationMagnetization)/h,0.);
    int lo=0,hi=law.count-1;
    while(hi-lo>1){int m=(lo+hi)/2;if(h<law.knots[m])hi=m;else lo=m;}
    const double x=h-law.knots[lo];const double* c=law.coefficients+4*lo;
    return make_cuDoubleComplex(((c[0]*x+c[1])*x+c[2])*x/h+c[3]/h,0.);
}
__device__ inline void cellRow(const ESIMCellParameters& p,int i,cuDoubleComplex mu,double surface,
        cuDoubleComplex& lower,cuDoubleComplex& diagonal,cuDoubleComplex& upper,cuDoubleComplex& rhs) {
    const double h=p.halfThickness/(p.nodes-1),coef=p.resistivity/(h*h);
    lower=upper=rhs=make_cuDoubleComplex(0.,0.);
    if((!p.cylinder&&i==0)||(p.cylinder&&i==p.nodes-1)) {
        diagonal=make_cuDoubleComplex(1.,0.);rhs=make_cuDoubleComplex(surface,0.);return;
    }
    double l=coef,u=coef;
    if(p.cylinder&&i==0){l=0.;u=4.*coef;}
    else if(p.cylinder){l=coef*(1.-.5/i);u=coef*(1.+.5/i);}
    else if(i==p.nodes-1){l=2.*coef;u=0.;}
    lower=make_cuDoubleComplex(l,0.);upper=make_cuDoubleComplex(u,0.);
    diagonal=cuCadd(make_cuDoubleComplex(-l-u,0.),cuCmul(make_cuDoubleComplex(0.,p.omega),mu));
}
// One cell arithmetic is shared by uniform and element-specific material paths.
__device__ inline ESIMCellResult solveFiniteESIMCell(ESIMCellParameters p,double surface,cuDoubleComplex* H) {
    ESIMCellResult result{};result.status=-1;
    if(p.nodes<3||p.maxIterations<1||!isfinite(surface)||surface<=0.||
       !isfinite(p.halfThickness)||p.halfThickness<=0.||!isfinite(p.resistivity)||p.resistivity<=0.||
       !isfinite(p.omega)||p.omega<0.||!isfinite(p.tolerance)||p.tolerance<=0.||
       !isfinite(p.relaxation)||p.relaxation<=0.||p.relaxation>1.||!cellFinite(p.law.initialMu)||
       cuCreal(p.law.initialMu)<=0.||cuCimag(p.law.initialMu)>0.) {return result;}
    auto* mu=H+p.nodes;auto* cp=mu+p.nodes;auto* dp=cp+p.nodes;
    for(int i=0;i<p.nodes;++i)mu[i]=p.law.initialMu;
    for(int iteration=0;iteration<p.maxIterations;++iteration) {
        for(int i=0;i<p.nodes;++i) {
            cuDoubleComplex l,d,u,b;cellRow(p,i,mu[i],surface,l,d,u,b);
            if(i){d=cuCsub(d,cuCmul(l,cp[i-1]));b=cuCsub(b,cuCmul(l,dp[i-1]));}
            if(!cellFinite(d)||cuCabs(d)==0.){result.status=-2;return result;}
            cp[i]=cuCdiv(u,d);dp[i]=cuCdiv(b,d);
        }
        H[p.nodes-1]=dp[p.nodes-1];
        for(int i=p.nodes-2;i>=0;--i)H[i]=cuCsub(dp[i],cuCmul(cp[i],H[i+1]));
        double residualSquared=0.,change=0.,muScale=0.;
        for(int i=0;i<p.nodes;++i) {
            cuDoubleComplex l,d,u,b;cellRow(p,i,mu[i],surface,l,d,u,b);
            auto residual=cuCsub(cuCmul(d,H[i]),b);
            if(i)residual=cuCadd(residual,cuCmul(l,H[i-1]));
            if(i+1<p.nodes)residual=cuCadd(residual,cuCmul(u,H[i+1]));
            const double r=cuCabs(residual);residualSquared+=r*r;
            auto next=cellMu(p.law,cuCabs(H[i]));
            if(!cellFinite(H[i])||!cellFinite(next)||cuCreal(next)<=0.||cuCimag(next)>0.){result.status=-2;return result;}
            change=fmax(change,cuCabs(cuCsub(next,mu[i])));muScale=fmax(muScale,cuCabs(mu[i]));
            // cp is no longer needed after back substitution.
            cp[i]=next;
        }
        result.trueResidual=sqrt(residualSquared)/surface;
        if(!isfinite(result.trueResidual)||result.trueResidual>1e-6){result.status=-2;return result;}
        result.relativeChange=change/muScale;result.iterations=iteration+1;
        for(int i=0;i<p.nodes;++i)mu[i]=cuCadd(cuCmul(make_cuDoubleComplex(1.-p.relaxation,0.),mu[i]),cuCmul(make_cuDoubleComplex(p.relaxation,0.),cp[i]));
        if(result.relativeChange<p.tolerance){result.status=1;break;}
    }
    if(result.status!=1){result.status=-3;return result;}
    const double h=p.halfThickness/(p.nodes-1);
    for(int i=0;i+1<p.nodes;++i) {
        const auto derivative=cuCdiv(cuCsub(H[i+1],H[i]),make_cuDoubleComplex(h,0.));
        const auto midpoint=cuCmul(make_cuDoubleComplex(.5,0.),cuCadd(H[i],H[i+1]));
        const auto muMid=cuCmul(make_cuDoubleComplex(.5,0.),cuCadd(mu[i],mu[i+1]));
        const double weight=p.cylinder?(i+.5)*h/p.halfThickness:1.;
        const double h2=cuCreal(midpoint)*cuCreal(midpoint)+cuCimag(midpoint)*cuCimag(midpoint);
        result.activePower+=.5*p.resistivity*pow(cuCabs(derivative),2)*weight*h;
        result.magneticPower+=-.5*p.omega*cuCimag(muMid)*h2*weight*h;
        result.reactivePower+=.5*p.omega*cuCreal(muMid)*h2*weight*h;
    }
    result.activePower+=result.magneticPower;
    result.impedance=make_cuDoubleComplex(2.*result.activePower/(surface*surface),2.*result.reactivePower/(surface*surface));
    if(!cellFinite(result.impedance)||result.activePower<0.)result.status=-2;
    return result;
}
// Uniform-law entry remains available for existing differential tests.
__global__ void finiteESIMCells(ESIMCellParameters p,const double* surfaceFields,std::size_t panels,
        cuDoubleComplex* workspace,ESIMCellResult* results) {
    const auto panel=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(panel>=panels)return;
    results[panel]=solveFiniteESIMCell(p,surfaceFields[panel],workspace+panel*std::size_t(p.nodes)*4);
}
// Element laws are staged on device from the accepted local temperature.
// A fixed stride is checked by the owner; no shared temperature or element-0 law.
__global__ void finiteElementESIMCells(const ESIMCellParameters* parameters,const double* surfaceFields,
        std::size_t panels,int workspaceNodes,cuDoubleComplex* workspace,ESIMCellResult* results) {
    const auto panel=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(panel>=panels)return;
    const auto p=parameters[panel];
    if(p.nodes<3||p.nodes>workspaceNodes){ESIMCellResult r{};r.status=-1;results[panel]=r;return;}
    results[panel]=solveFiniteESIMCell(p,surfaceFields[panel],workspace+panel*std::size_t(workspaceNodes)*4);
}
} // namespace radia_ih_cuda
