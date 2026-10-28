// FP64 material solve A = B^T G B + inv_chi M, explicit diagonal PCG.
// This is not the default CPU mass-Riesz preconditioner, nor an FE assembler.
#pragma once
#include "hmatrix_resident.cuh"
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <sstream>
#include <functional>
namespace radia_hmatrix_cuda {
inline void blasCheck(cublasStatus_t s){if(s!=CUBLAS_STATUS_SUCCESS)throw std::runtime_error("HDiv PCG cuBLAS failure");}
inline void solverCheck(cusolverStatus_t status){if(status!=CUSOLVER_STATUS_SUCCESS)throw std::runtime_error("HDiv mass Cholesky cuSOLVER failure");}
struct SolverHandle {
    cusolverDnHandle_t value=nullptr;
    explicit SolverHandle(cudaStream_t stream){resourceFaultPoint("cusolver.create");solverCheck(cusolverDnCreate(&value));resourceAcquired();try{resourceFaultPoint("cusolver.after-create");solverCheck(cusolverDnSetStream(value,stream));}catch(...){resourceReleased(cusolverDnDestroy(value)==CUSOLVER_STATUS_SUCCESS);value=nullptr;throw;}}
    ~SolverHandle(){if(value)resourceReleased(cusolverDnDestroy(value)==CUSOLVER_STATUS_SUCCESS);}
};
struct BlasHandle {
    cublasHandle_t value=nullptr;
    explicit BlasHandle(cudaStream_t stream){resourceFaultPoint("cublas.create");blasCheck(cublasCreate(&value));resourceAcquired();try{resourceFaultPoint("cublas.after-create");blasCheck(cublasSetStream(value,stream));}catch(...){resourceReleased(cublasDestroy(value)==CUBLAS_STATUS_SUCCESS);value=nullptr;throw;}}
    ~BlasHandle(){if(value)resourceReleased(cublasDestroy(value)==CUBLAS_STATUS_SUCCESS);}
};
__global__ void materialSum(std::size_t n,double invChi,const double* mass,double* y){auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n)y[i]+=invChi*mass[i];}
__global__ void pcgPrecondition(std::size_t n,const double* diagonal,const double* r,double* z){auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n)z[i]=r[i]/diagonal[i];}
__global__ void pcgUpdate(std::size_t n,double alpha,const double* p,const double* ap,double* x,double* r){auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n){x[i]+=alpha*p[i];r[i]-=alpha*ap[i];}}
__global__ void pcgDirection(std::size_t n,double beta,const double* z,double* p){auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n)p[i]=z[i]+beta*p[i];}
__global__ void trueResidual(std::size_t n,const double* b,const double* ax,double* r){auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n)r[i]=b[i]-ax[i];}
__global__ void pcgPolynomialStep(std::size_t n,double omega,double invChi,const double* diagonal,const double* r,const double* mz,double* z){auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n)z[i]+=omega*(r[i]-invChi*mz[i])/diagonal[i];}
__global__ void pcgChebyshevStep(std::size_t n,double alpha,double beta,double invChi,const double* diagonal,const double* r,const double* mz,double* direction,double* z){auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n){double v=alpha*(r[i]-invChi*mz[i])/diagonal[i];if(beta!=0.)v+=beta*direction[i];direction[i]=v;z[i]+=v;}}
struct MaterialPreconditionGraph {
    cudaGraph_t graph=nullptr;cudaGraphExec_t executable=nullptr;
    MaterialPreconditionGraph()=default;
    MaterialPreconditionGraph(const MaterialPreconditionGraph&)=delete;
    MaterialPreconditionGraph& operator=(const MaterialPreconditionGraph&)=delete;
    ~MaterialPreconditionGraph(){if(executable)resourceReleased(cudaGraphExecDestroy(executable)==cudaSuccess);if(graph)resourceReleased(cudaGraphDestroy(graph)==cudaSuccess);}
    template<class Submit> void capture(cudaStream_t stream,Submit submit){
        cudaStreamCaptureStatus status;check(cudaStreamIsCapturing(stream,&status));
        if(status!=cudaStreamCaptureStatusNone)throw std::runtime_error("Material graph requires an uncaptured owned stream");
        resourceFaultPoint("graph.begin-capture");check(cudaStreamBeginCapture(stream,cudaStreamCaptureModeThreadLocal));
        try{submit();resourceFaultPoint("graph.after-submit");}catch(...){cudaGraph_t failed=nullptr;cudaStreamEndCapture(stream,&failed);if(failed){resourceAcquired();resourceReleased(cudaGraphDestroy(failed)==cudaSuccess);}throw;}
        check(cudaStreamEndCapture(stream,&graph));if(graph)resourceAcquired();
        std::size_t nodes=0;check(cudaGraphGetNodes(graph,nullptr,&nodes));
        if(!nodes||nodes>1024)throw std::runtime_error("Material graph node budget exceeded");
        resourceFaultPoint("graph.instantiate");check(cudaGraphInstantiate(&executable,graph,nullptr,nullptr,0));resourceAcquired();resourceFaultPoint("graph.after-instantiate");
        std::size_t freeBytes=0,totalBytes=0;check(cudaMemGetInfo(&freeBytes,&totalBytes));
        if(freeBytes<std::max(std::size_t(2)<<30,totalBytes/5))throw std::runtime_error("Material graph violates device reserve");
    }
    void launch(cudaStream_t stream){check(cudaGraphLaunch(executable,stream));}
};
struct MaterialSolution {HostVector<double> x;int iterations=0;double relativeResidual=0;std::size_t ownedBytes=0,preconditionerGraphLaunches=0;};
template<class Allocator>
inline MaterialSolution solveLinearMaterial(Resident& op,const SparseMap& mass,const std::vector<double,Allocator>& rhs,double invChi,double tolerance,int maxIterations,bool massCholesky=false,bool massPolynomial=false,bool polynomialGraph=false,bool massChebyshev=false,int chebyshevDegree=8,double chebyshevLowerFraction=.1,const std::function<void(const double*,double*,cudaStream_t)>& materialAction={}){
    const auto n=op.dimension();
    if(!op.isDemag()||mass.rows!=n||mass.columns!=n||rhs.size()!=n||n>2147483647u||!std::isfinite(invChi)||invChi<=0||!std::isfinite(tolerance)||tolerance<=0||tolerance>1.e-6||maxIterations<1||maxIterations>100000)
        throw std::invalid_argument("Invalid linear material PCG configuration");
    if(polynomialGraph&&!massPolynomial&&!massChebyshev)throw std::invalid_argument("Material graph requires fixed polynomial preconditioner");
    if(chebyshevDegree<1||chebyshevDegree>64||!std::isfinite(chebyshevLowerFraction)||chebyshevLowerFraction<=0||chebyshevLowerFraction>=1)throw std::invalid_argument("Invalid fixed Chebyshev configuration");
    mass.validate();const auto transpose=mass.transpose();HostVector<double> diag(n,0.);
    // PCG requires a symmetric mass. No silent symmetrization is permitted.
    for(std::size_t row=0;row<n;++row){
        if(mass.offsets[row+1]-mass.offsets[row]!=transpose.offsets[row+1]-transpose.offsets[row])throw std::invalid_argument("Material mass is not symmetric");
        for(auto k=mass.offsets[row],t=transpose.offsets[row];k<mass.offsets[row+1];++k,++t){
            if(mass.indices[k]!=transpose.indices[t]||std::abs(mass.values[k]-transpose.values[t])>1.e-13*std::max(1.,std::max(std::abs(mass.values[k]),std::abs(transpose.values[t]))))throw std::invalid_argument("Material mass is not symmetric");
            if(mass.indices[k]==row)diag[row]=invChi*mass.values[k];
        }
        if(!std::isfinite(rhs[row])||!std::isfinite(diag[row])||diag[row]<=0)throw std::invalid_argument("Invalid RHS or material diagonal");
    }
    // Eight fixed damped-Jacobi steps from zero are a symmetric polynomial
    // inverse action: D^-1/2 omega sum_{j=0}^7(I-omega Ahat)^j D^-1/2.
    // Ahat=invChi D^-1/2 M D^-1/2. Its spectral maximum is bounded by
    // the infinity norm of the similar D^-1 invChi M. omega=1/bound
    // makes I-omega Ahat positive semidefinite, hence this inverse action
    // is strictly SPD even without claiming an exact mass inverse.
    double polynomialOmega=1.,spectralBound=1.;
    if(massPolynomial||massChebyshev){
        double bound=0.;
        for(std::size_t row=0;row<n;++row){double sum=0.;for(auto k=mass.offsets[row];k<mass.offsets[row+1];++k)sum+=std::abs(invChi*mass.values[k]/diag[row]);bound=std::max(bound,sum);}
        if(!std::isfinite(bound)||bound<1.)throw std::invalid_argument("Invalid sparse polynomial spectral bound");
        polynomialOmega=1./bound;spectralBound=bound;
    }
    MaterialSolution result;result.ownedBytes=checkedProduct(n,8*sizeof(double))+(mass.offsets.size()+mass.indices.size())*sizeof(std::size_t)+mass.values.size()*sizeof(double);
    if(massChebyshev)result.ownedBytes+=checkedProduct(n,sizeof(double));
    std::size_t free,total;check(cudaMemGetInfo(&free,&total));const auto reserve=std::max(std::size_t(2)<<30,total/5);
    if(free<=reserve||result.ownedBytes>free-reserve)throw std::runtime_error("Material PCG GPU memory reserve exceeded");
    auto stream=op.executionStream();BlasHandle blas(stream);DeviceMap m;m.upload(mass);
    Buffer<double> b,x,r,z,p,ap,mx,d;b.upload(rhs);d.upload(diag);
    for(auto* v:{&x,&r,&z,&p,&ap,&mx})v->allocate(n);
    Buffer<double> chebyshevDirection;if(massChebyshev)chebyshevDirection.allocate(n);
    Buffer<double> massFactor,factorWorkspace;Buffer<int> factorInfo;
    if(massCholesky){
        const auto factorBytes=checkedProduct(checkedProduct(n,n),sizeof(double));
        if(factorBytes>(std::size_t(2)<<30))throw std::invalid_argument("Mass Cholesky factor exceeds bounded 2 GiB");
        check(cudaMemGetInfo(&free,&total));
        if(free<=reserve||factorBytes>free-reserve)throw std::runtime_error("Mass Cholesky factor reserve exceeded");
        HostVector<double> dense(checkedProduct(n,n),0.);
        for(std::size_t row=0;row<n;++row)for(auto k=mass.offsets[row];k<mass.offsets[row+1];++k)dense[row*n+mass.indices[k]]=mass.values[k];
        massFactor.upload(dense);SolverHandle solver(stream);int workspace=0;
        solverCheck(cusolverDnDpotrf_bufferSize(solver.value,CUBLAS_FILL_MODE_LOWER,int(n),massFactor.data,int(n),&workspace));
        if(workspace<0)throw std::runtime_error("Invalid mass factor workspace");
        const auto workspaceBytes=checkedProduct(std::size_t(workspace),sizeof(double))+sizeof(int);
        check(cudaMemGetInfo(&free,&total));
        if(free<=reserve||workspaceBytes>free-reserve)throw std::runtime_error("Mass Cholesky workspace reserve exceeded");
        factorWorkspace.allocate(std::size_t(workspace));factorInfo.allocate(1);
        solverCheck(cusolverDnDpotrf(solver.value,CUBLAS_FILL_MODE_LOWER,int(n),massFactor.data,int(n),factorWorkspace.data,workspace,factorInfo.data));
        int info=0;check(cudaMemcpyAsync(&info,factorInfo.data,sizeof(int),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));
        if(info!=0)throw std::runtime_error("Mass Cholesky rejected non-SPD mass; no fallback");
        result.ownedBytes+=factorBytes+workspaceBytes;
    }
    check(cudaMemsetAsync(x.data,0,n*sizeof(double),stream));check(cudaMemcpyAsync(r.data,b.data,n*sizeof(double),cudaMemcpyDeviceToDevice,stream));
    auto dot=[&](const double* a,const double* b){double v;blasCheck(cublasDdot(blas.value,int(n),a,1,b,1,&v));if(!std::isfinite(v))throw std::runtime_error("Nonfinite PCG reduction");return v;};
    auto apply=[&](const double* v){op.applyDevice(v,ap.data);if(materialAction)materialAction(v,mx.data,stream);else m.apply(v,mx.data,stream);materialSum<<<unsigned((n+255)/256),256,0,stream>>>(n,invChi,mx.data,ap.data);check(cudaGetLastError());};
    auto polynomialSubmit=[&](){
        check(cudaMemsetAsync(z.data,0,n*sizeof(double),stream));
        for(int step=0;step<8;++step){m.apply(z.data,mx.data,stream);pcgPolynomialStep<<<unsigned((n+255)/256),256,0,stream>>>(n,polynomialOmega,invChi,d.data,r.data,mx.data,z.data);check(cudaGetLastError());}
    };
    auto chebyshevSubmit=[&](){
            // Explicit fixed degree; no RHS-dependent stopping. Upper bound is
            // the infinity norm of D^-1 invChi M with a rounding guard.
            // q(lambda)=(1-Tm((theta-lambda)/delta)/Tm(theta/delta))/lambda
            // is positive for every 0<lambda<=upper, including below lower.
            // Thus D^-1/2 q(Ahat) D^-1/2 is SPD for SPD material mass.
            // The declared lower fraction is a smoothing target, not lambda_min.
            const double upper=spectralBound*(1.+1.e-12),lower=upper*chebyshevLowerFraction;
            const double theta=(upper+lower)*.5,delta=(upper-lower)*.5,sigma=theta/delta;
            double rho=1./sigma;check(cudaMemsetAsync(z.data,0,n*sizeof(double),stream));
            for(int step=0;step<chebyshevDegree;++step){
                double alpha=1./theta,beta=0.;
                if(step){const double next=1./(2.*sigma-rho);beta=next*rho;alpha=2.*next/delta;rho=next;}
                m.apply(z.data,mx.data,stream);
                pcgChebyshevStep<<<unsigned((n+255)/256),256,0,stream>>>(n,alpha,beta,invChi,d.data,r.data,mx.data,chebyshevDirection.data,z.data);check(cudaGetLastError());
            }
    };
    // Graph captures the same fixed arithmetic only. Buffers outlive graph,
    // and every launch uses this solve's immutable mass/diagonal/fixed polynomial parameters.
    // Driver-owned graph allocations are not included in ownedBytes.
    MaterialPreconditionGraph preconditionGraph;
    if(polynomialGraph){
        if(massChebyshev)preconditionGraph.capture(stream,chebyshevSubmit);
        else preconditionGraph.capture(stream,polynomialSubmit);
    }
    auto precondition=[&](){
        if(polynomialGraph){preconditionGraph.launch(stream);++result.preconditionerGraphLaunches;
        }else if(massChebyshev){
            chebyshevSubmit();
        }else if(massPolynomial){
            polynomialSubmit();
        }else if(massCholesky){
            check(cudaMemcpyAsync(z.data,r.data,n*sizeof(double),cudaMemcpyDeviceToDevice,stream));
            blasCheck(cublasDtrsv(blas.value,CUBLAS_FILL_MODE_LOWER,CUBLAS_OP_N,CUBLAS_DIAG_NON_UNIT,int(n),massFactor.data,int(n),z.data,1));
            blasCheck(cublasDtrsv(blas.value,CUBLAS_FILL_MODE_LOWER,CUBLAS_OP_T,CUBLAS_DIAG_NON_UNIT,int(n),massFactor.data,int(n),z.data,1));
        }else{pcgPrecondition<<<unsigned((n+255)/256),256,0,stream>>>(n,d.data,r.data,z.data);check(cudaGetLastError());}
    };
    const double b2=dot(b.data,b.data);if(b2<0)throw std::runtime_error("Invalid RHS norm");
    if(b2>0){
        precondition();check(cudaMemcpyAsync(p.data,z.data,n*sizeof(double),cudaMemcpyDeviceToDevice,stream));double rz=dot(r.data,z.data);
        bool converged=false;
        for(int k=0;k<maxIterations;++k){
            apply(p.data);double pap=dot(p.data,ap.data);if(pap<=0||rz<=0){std::ostringstream message;message.precision(17);message<<"Material PCG lost positive definiteness: iteration="<<k<<" pAp="<<pap<<" rP r="<<rz<<" b2="<<b2<<" r2="<<dot(r.data,r.data);throw std::runtime_error(message.str());}
            double alpha=rz/pap;if(!std::isfinite(alpha))throw std::runtime_error("Invalid PCG step");
            pcgUpdate<<<unsigned((n+255)/256),256,0,stream>>>(n,alpha,p.data,ap.data,x.data,r.data);check(cudaGetLastError());result.iterations=k+1;
            if(dot(r.data,r.data)<=tolerance*tolerance*b2){
                apply(x.data);trueResidual<<<unsigned((n+255)/256),256,0,stream>>>(n,b.data,ap.data,r.data);check(cudaGetLastError());
                if(dot(r.data,r.data)<=tolerance*tolerance*b2){converged=true;break;}
                // Residual replacement restarts CG; never report recursive
                // convergence when the original equation fails its gate.
                precondition();check(cudaMemcpyAsync(p.data,z.data,n*sizeof(double),cudaMemcpyDeviceToDevice,stream));rz=dot(r.data,z.data);continue;
            }
            precondition();double next=dot(r.data,z.data);double beta=next/rz;if(!std::isfinite(beta)||next<=0)throw std::runtime_error("Invalid PCG recurrence");
            pcgDirection<<<unsigned((n+255)/256),256,0,stream>>>(n,beta,z.data,p.data);check(cudaGetLastError());rz=next;
        }
        if(!converged){
            apply(x.data);trueResidual<<<unsigned((n+255)/256),256,0,stream>>>(n,b.data,ap.data,r.data);check(cudaGetLastError());
            const double actual=std::sqrt(dot(r.data,r.data)/b2);
            std::ostringstream message;message.precision(17);
            message<<"Material PCG failed true residual convergence: iterations="<<result.iterations
                   <<" requested="<<tolerance<<" actual="<<actual;
            throw std::runtime_error(message.str());
        }
    }
    apply(x.data);trueResidual<<<unsigned((n+255)/256),256,0,stream>>>(n,b.data,ap.data,r.data);check(cudaGetLastError());
    result.relativeResidual=b2>0?std::sqrt(dot(r.data,r.data)/b2):std::sqrt(dot(r.data,r.data));
    if(result.relativeResidual>tolerance||!std::isfinite(result.relativeResidual))throw std::runtime_error("Material PCG true residual rejected");
    result.x.resize(n);check(cudaMemcpyAsync(result.x.data(),x.data,n*sizeof(double),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));return result;
}
}
