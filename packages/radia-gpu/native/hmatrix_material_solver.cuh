// FP64 material solve A = B^T G B + inv_chi M, explicit diagonal PCG.
// This is not the default CPU mass-Riesz preconditioner, nor an FE assembler.
#pragma once
#include "hmatrix_resident.cuh"
#include <cublas_v2.h>
namespace radia_hmatrix_cuda {
inline void blasCheck(cublasStatus_t s){if(s!=CUBLAS_STATUS_SUCCESS)throw std::runtime_error("HDiv PCG cuBLAS failure");}
struct BlasHandle {
    cublasHandle_t value=nullptr;
    explicit BlasHandle(cudaStream_t stream){blasCheck(cublasCreate(&value));try{blasCheck(cublasSetStream(value,stream));}catch(...){cublasDestroy(value);throw;}}
    ~BlasHandle(){if(value)cublasDestroy(value);}
};
__global__ void materialSum(std::size_t n,double invChi,const double* mass,double* y){auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n)y[i]+=invChi*mass[i];}
__global__ void pcgPrecondition(std::size_t n,const double* diagonal,const double* r,double* z){auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n)z[i]=r[i]/diagonal[i];}
__global__ void pcgUpdate(std::size_t n,double alpha,const double* p,const double* ap,double* x,double* r){auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n){x[i]+=alpha*p[i];r[i]-=alpha*ap[i];}}
__global__ void pcgDirection(std::size_t n,double beta,const double* z,double* p){auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n)p[i]=z[i]+beta*p[i];}
__global__ void trueResidual(std::size_t n,const double* b,const double* ax,double* r){auto i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n)r[i]=b[i]-ax[i];}
struct MaterialSolution {std::vector<double> x;int iterations=0;double relativeResidual=0;std::size_t ownedBytes=0;};
inline MaterialSolution solveLinearMaterial(Resident& op,const SparseMap& mass,const std::vector<double>& rhs,double invChi,double tolerance,int maxIterations){
    const auto n=op.dimension();
    if(!op.isDemag()||mass.rows!=n||mass.columns!=n||rhs.size()!=n||n>2147483647u||!std::isfinite(invChi)||invChi<=0||!std::isfinite(tolerance)||tolerance<=0||tolerance>1.e-6||maxIterations<1||maxIterations>100000)
        throw std::invalid_argument("Invalid linear material PCG configuration");
    mass.validate();const auto transpose=mass.transpose();std::vector<double> diag(n,0.);
    // PCG requires a symmetric mass. No silent symmetrization is permitted.
    for(std::size_t row=0;row<n;++row){
        if(mass.offsets[row+1]-mass.offsets[row]!=transpose.offsets[row+1]-transpose.offsets[row])throw std::invalid_argument("Material mass is not symmetric");
        for(auto k=mass.offsets[row],t=transpose.offsets[row];k<mass.offsets[row+1];++k,++t){
            if(mass.indices[k]!=transpose.indices[t]||std::abs(mass.values[k]-transpose.values[t])>1.e-13*std::max(1.,std::max(std::abs(mass.values[k]),std::abs(transpose.values[t]))))throw std::invalid_argument("Material mass is not symmetric");
            if(mass.indices[k]==row)diag[row]=invChi*mass.values[k];
        }
        if(!std::isfinite(rhs[row])||!std::isfinite(diag[row])||diag[row]<=0)throw std::invalid_argument("Invalid RHS or material diagonal");
    }
    MaterialSolution result;result.ownedBytes=checkedProduct(n,8*sizeof(double))+(mass.offsets.size()+mass.indices.size())*sizeof(std::size_t)+mass.values.size()*sizeof(double);
    std::size_t free,total;check(cudaMemGetInfo(&free,&total));const auto reserve=std::max(std::size_t(2)<<30,total/5);
    if(free<=reserve||result.ownedBytes>free-reserve)throw std::runtime_error("Material PCG GPU memory reserve exceeded");
    auto stream=op.executionStream();BlasHandle blas(stream);DeviceMap m;m.upload(mass);
    Buffer<double> b,x,r,z,p,ap,mx,d;b.upload(rhs);d.upload(diag);
    for(auto* v:{&x,&r,&z,&p,&ap,&mx})v->allocate(n);
    check(cudaMemsetAsync(x.data,0,n*sizeof(double),stream));check(cudaMemcpyAsync(r.data,b.data,n*sizeof(double),cudaMemcpyDeviceToDevice,stream));
    auto dot=[&](const double* a,const double* b){double v;blasCheck(cublasDdot(blas.value,int(n),a,1,b,1,&v));if(!std::isfinite(v))throw std::runtime_error("Nonfinite PCG reduction");return v;};
    auto apply=[&](const double* v){op.applyDevice(v,ap.data);m.apply(v,mx.data,stream);materialSum<<<unsigned((n+255)/256),256,0,stream>>>(n,invChi,mx.data,ap.data);check(cudaGetLastError());};
    auto precondition=[&](){pcgPrecondition<<<unsigned((n+255)/256),256,0,stream>>>(n,d.data,r.data,z.data);check(cudaGetLastError());};
    const double b2=dot(b.data,b.data);if(b2<0)throw std::runtime_error("Invalid RHS norm");
    if(b2>0){
        precondition();check(cudaMemcpyAsync(p.data,z.data,n*sizeof(double),cudaMemcpyDeviceToDevice,stream));double rz=dot(r.data,z.data);
        bool converged=false;
        for(int k=0;k<maxIterations;++k){
            apply(p.data);double pap=dot(p.data,ap.data);if(pap<=0||rz<=0)throw std::runtime_error("Material PCG lost positive definiteness");
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
        if(!converged)throw std::runtime_error("Material PCG failed true residual convergence");
    }
    apply(x.data);trueResidual<<<unsigned((n+255)/256),256,0,stream>>>(n,b.data,ap.data,r.data);check(cudaGetLastError());
    result.relativeResidual=b2>0?std::sqrt(dot(r.data,r.data)/b2):std::sqrt(dot(r.data,r.data));
    if(result.relativeResidual>tolerance||!std::isfinite(result.relativeResidual))throw std::runtime_error("Material PCG true residual rejected");
    result.x.resize(n);check(cudaMemcpyAsync(result.x.data(),x.data,n*sizeof(double),cudaMemcpyDeviceToHost,stream));check(cudaStreamSynchronize(stream));return result;
}
}
