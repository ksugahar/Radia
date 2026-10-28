#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>
namespace radia_csr {
struct Result { int iterations; double true_relative_residual; };
inline void validate_dimensions(size_t n,size_t nnz,double tol,int maxiter) {
    if(n==0 || n>INT32_MAX || nnz>INT32_MAX || !(tol>0 && tol<=1e-6)
       || maxiter<1 || maxiter>100000)
        throw std::runtime_error("Invalid CSR dimensions, tolerance or iteration limit");
    // Bounded candidate: conservative per-phase matrix plus vector accounting.
    // int32 dimensions and a 64-bit size_t make these products non-overflowing.
    static_assert(sizeof(size_t)>=8,"CSR requires a 64-bit process");
    const size_t bound=(n+1)*sizeof(int)+nnz*(sizeof(int)+sizeof(double))+n*8*sizeof(double)+9*sizeof(double);
    if(bound>size_t(8)*1024*1024*1024)
        throw std::runtime_error("CSR phase exceeds the 8 GiB candidate memory bound");
}
inline void validate_host(int n,int nnz,const int* row,const int* col,
                          const double* val,const double* rhs,const double* dinv) {
    if(row[0]!=0 || row[n]!=nnz) throw std::runtime_error("Invalid CSR endpoints");
    for(int i=0;i<n;++i) if(row[i]<0 || row[i+1]<row[i] || row[i+1]>nnz
        || !std::isfinite(rhs[i]) || !std::isfinite(dinv[i]) || !(dinv[i]>0))
        throw std::runtime_error("Invalid CSR row, RHS or inverse diagonal");
    for(int j=0;j<nnz;++j) if(col[j]<0 || col[j]>=n || !std::isfinite(val[j]))
        throw std::runtime_error("Invalid CSR column or value");
}
inline void multiply(int n,const int* row,const int* col,const double* val,
                     const double* x,double* y) {
    for(int i=0;i<n;++i){double sum=0;for(int j=row[i];j<row[i+1];++j)sum+=val[j]*x[col[j]];y[i]=sum;}
}
inline double dot(int n,const double* a,const double* b) {
    double sum=0;for(int i=0;i<n;++i)sum+=a[i]*b[i];return sum;
}
// Serial host execution deliberately has no CUDA or TaskManager dependency.
// The caller supplies an SPD system and a positive Jacobi inverse diagonal.
inline Result solve_cpu(int n,int nnz,const int* row,const int* col,const double* val,
    const double* rhs,const double* dinv,double tol,int maxiter,double* x) {
    validate_dimensions(n,nnz,tol,maxiter);validate_host(n,nnz,row,col,val,rhs,dinv);
    std::fill(x,x+n,0.0);
    std::vector<double> r(rhs,rhs+n),p(n),q(n),z(n);
    const double bnorm2=dot(n,rhs,rhs);
    bool nonzero=false;for(int i=0;i<n;++i)nonzero=nonzero || rhs[i]!=0;
    if(!std::isfinite(bnorm2) || (bnorm2==0 && nonzero))
        throw std::runtime_error("RHS norm outside supported FP64 range; rescale explicitly");
    int iterations=0;
    if(bnorm2>0){
        for(int i=0;i<n;++i)p[i]=z[i]=dinv[i]*r[i];
        double rho=dot(n,r.data(),z.data());
        for(int k=0;k<maxiter;++k){
            multiply(n,row,col,val,p.data(),q.data());
            const double pap=dot(n,p.data(),q.data());
            if(!(pap>0) || !std::isfinite(pap) || !(rho>0) || !std::isfinite(rho))
                throw std::runtime_error("PCG requires a finite SPD system");
            const double alpha=rho/pap;
            for(int i=0;i<n;++i){x[i]+=alpha*p[i];r[i]-=alpha*q[i];}
            iterations=k+1;
            const double rr=dot(n,r.data(),r.data());
            if(!std::isfinite(rr))throw std::runtime_error("PCG residual overflow");
            if(rr<=tol*tol*bnorm2)break;
            for(int i=0;i<n;++i)z[i]=dinv[i]*r[i];
            const double next=dot(n,r.data(),z.data()),beta=next/rho;
            for(int i=0;i<n;++i)p[i]=z[i]+beta*p[i];
            rho=next;
        }
    }
    multiply(n,row,col,val,x,q.data());
    double relative2=0;
    for(int i=0;i<n;++i){double e=bnorm2>0 ? (rhs[i]-q[i])/std::sqrt(bnorm2) : rhs[i]-q[i];relative2+=e*e;}
    const double relative=std::sqrt(relative2);
    if(!std::isfinite(relative) || relative>1e-6 || relative>std::max(10*tol,1e-12))
        throw std::runtime_error("PCG true original-system residual exceeds acceptance");
    return {iterations,relative};
}
// Defined by the CUDA translation unit; uses host arrays with explicit transfers.
Result solve_gpu_host(int n,int nnz,const int*,const int*,const double*,const double*,
                      const double*,double,int,double*);
// Internal AMG smoothing phase, not a complete preconditioner or solver.
void amg_smooth_gpu_host(int n,int nnz,const int*,const int*,const double*,
                         const double*,const double*,int,double*);
} // namespace radia_csr
