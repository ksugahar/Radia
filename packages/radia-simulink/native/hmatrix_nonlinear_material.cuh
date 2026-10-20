// Shared host Newton/globalization control; all domain operator/tangent work is CUDA.
#pragma once
#include "hmatrix_bh_material.cuh"
#include "hmatrix_material_solver.cuh"
#include "hdiv_resident_python_api.hpp"
namespace radia_hmatrix_cuda {
template<class A,class B>inline double newtonDot(const std::vector<double,A>& a,const std::vector<double,B>& b){
    double sum=0,correction=0;for(std::size_t i=0;i<a.size();++i){const double product=a[i]*b[i],value=sum+product;
        correction+=(std::abs(sum)>=std::abs(product)?(sum-value)+product:(product-value)+sum)+std::fma(a[i],b[i],-product);sum=value;}return sum+correction;
}
template<class A>inline double newtonNorm(const std::vector<double,A>& a){double norm=0;for(double x:a)norm=std::hypot(norm,x);return norm;}
template<class RhsAllocator,class StateAllocator>inline radia_hdiv_python::NewtonResult solveNonlinearMaterial(Resident& op,radia_hdiv_bh::MappedMaterial& material,
    const SparseMap& referenceMass,const std::vector<double,RhsAllocator>& rhs,const std::vector<double,StateAllocator>& initial,const radia_hdiv_python::NewtonOptions& options){
    const auto n=op.dimension();
    if(!op.isDemag()||material.dofs()!=n||rhs.size()!=n||initial.size()!=n||referenceMass.rows!=n||referenceMass.columns!=n||
       !std::isfinite(options.tolerance)||options.tolerance<=0||options.tolerance>1e-6||!std::isfinite(options.linearTolerance)||options.linearTolerance<=0||options.linearTolerance>1e-6||
       options.maximumIterations<1||options.maximumIterations>100||options.maximumLinearIterations<1||options.maximumLinearIterations>100000||options.preconditioner<0||options.preconditioner>2||options.degree<1||options.degree>64||!std::isfinite(options.lowerFraction)||options.lowerFraction<=0||options.lowerFraction>=1)
        throw std::invalid_argument("Invalid shared Newton options or dimensions");
    referenceMass.validate();
    for(double v:rhs)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite Newton RHS");
    for(double v:initial)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite Newton initial state");
    radia_hdiv_python::NewtonResult result;result.x.assign(initial.begin(),initial.end());
    const double sourceNorm=std::max(newtonNorm(rhs),std::numeric_limits<double>::min());
    auto demag=[&](const auto& x){op.upload(x);op.apply(1);return op.download();};
    auto residual=[&](const auto& mat,const auto& nm){HostVector<double> r(n,HostFaultAllocator<double>("newton.residual"));for(std::size_t i=0;i<n;++i)r[i]=mat.gradient[i]+nm[i]-rhs[i];return r;};
    auto energy=[&](const auto& mat,const auto& x,const auto& nm){return mat.energy+.5*newtonDot(x,nm)-newtonDot(rhs,x);};
    for(int iteration=0;iteration<=options.maximumIterations;++iteration){
        const auto mat=material.linearize(result.x);const auto nm=demag(result.x);auto r=residual(mat,nm);
        result.residual=newtonNorm(r)/sourceNorm;result.energy=energy(mat,result.x,nm);
        if(!std::isfinite(result.residual)||!std::isfinite(result.energy))throw std::runtime_error("Nonfinite Newton original residual or energy");
        // Re-evaluate the original equation at every accepted state; no trial certificate substitutes for this check.
        if(result.residual<=options.tolerance)return result;
        if(iteration==options.maximumIterations)throw std::runtime_error("Shared Newton failed original residual convergence");
        for(double& v:r)v=-v;
        const auto epoch=mat.linearizationEpoch;
        auto action=[&](const double* x,double* y,cudaStream_t stream){material.applyFrozenDevice(x,y,stream,epoch);};
        const bool chebyshev=options.preconditioner!=0,graph=options.preconditioner==2;
        // referenceMass is an explicitly chosen SPD preconditioner, never the nonlinear Jacobian.
        const auto linear=solveLinearMaterial(op,referenceMass,r,1.,options.linearTolerance,options.maximumLinearIterations,false,false,graph,chebyshev,options.degree,options.lowerFraction,action);
        result.graphLaunches+=linear.preconditionerGraphLaunches;HostVector<double> delta(linear.x.begin(),linear.x.end(),HostFaultAllocator<double>("newton.direction"));
        const double decrement=newtonDot(r,delta);if(!std::isfinite(decrement)||decrement<=0)throw std::runtime_error("Newton direction is not descending");
        double alpha=1;int backtracks=0;bool accepted=false,certified=false;HostVector<double> candidate(n,HostFaultAllocator<double>("newton.trial"));
        while(alpha>1e-10){
            for(std::size_t i=0;i<n;++i)candidate[i]=result.x[i]+alpha*delta[i];
            const auto trial=material.evaluate(candidate,{});const auto trialN=demag(candidate);const double trialEnergy=energy(trial,candidate,trialN);
            const double trialResidual=newtonNorm(residual(trial,trialN))/sourceNorm;
            if(std::isfinite(trialEnergy)&&std::isfinite(trialResidual)&&trialResidual<=options.tolerance){accepted=true;certified=true;break;}
            if(std::isfinite(trialEnergy)&&trialEnergy<=result.energy-1e-4*alpha*decrement){accepted=true;break;}
            alpha*=.5;++backtracks;
        }
        if(!accepted)throw std::runtime_error("Shared Newton Armijo line search failed; last accepted state unchanged");
        result.history.push_back({result.residual,result.energy,alpha,linear.relativeResidual,linear.iterations,backtracks,certified});result.x=std::move(candidate);
    }
    throw std::runtime_error("Shared Newton unreachable convergence state");
}
}
