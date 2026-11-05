// Laplace kernel expansion of the far blocks of the mixed charge Gram: hybrid cross approximation. The charge kernel
// is 1/(4 pi r), so a block between two separated clusters factors through tensor Chebyshev interpolation of r^-1 on
// their tight boxes (Borm and Grasedyck, Numer. Math. 101, 2005; the interpolation of Fong and Darve, JCP 228, 2009):
//   G(i,j) = sum_kl M_t(i,k) K(k,l) M_s(j,l), M(i,k) = int rho_i(x) L_k(x) dx / sigma_i, K(k,l) = 1/(4 pi |x_k - y_l|).
// K costs one kernel per entry where a Gram entry costs a double sum, so setup needs no per-entry quadrature. The
// moments, K, its factors and the block factors are built and kept on the device. Only cluster boxes and work lists
// cross, as host metadata. The factors go into the resident H-matrix: this is a setup route, not a matrix-free FMM.
#pragma once
#include "hdiv_flat_gram.cuh"
namespace radia_hmatrix_cuda {
constexpr int hcaOrderCap=8;
// K is cut below the block tolerance, because the moments and the recompression still follow.
constexpr double hcaKernelEps=.25;
struct HcaNodes {double x[hcaOrderCap],scale[hcaOrderCap];int n;};
inline HcaNodes hcaNodes(int n){HcaNodes e{};e.n=n;
    for(int k=0;k<n;++k)e.x[k]=std::cos(3.14159265358979323846*(2*k+1)/(2*n));
    for(int k=0;k<n;++k){double v=1.0;for(int j=0;j<n;++j)if(j!=k)v*=e.x[k]-e.x[j];e.scale[k]=1.0/v;}
    return e;}
// Tight box of one cluster, the first of its moments and its charges in permuted order.
struct HcaBox {double lo[3],hi[3];std::size_t table,charge,charges;};
struct HcaMomentItem {std::size_t box,start;int host;};
// One factor side: rows x rank from the moments of a cluster and q x rank of K.
struct HcaProductSpan {std::size_t rows,rank,table,factor,out;};
__device__ inline double hcaLagrange(const HcaNodes& e,int k,double x){double v=e.scale[k];for(int j=0;j<e.n;++j)if(j!=k)v*=x-e.x[j];return v;}
__device__ inline void hcaNode(const HcaNodes& e,const HcaBox& b,int k,double* xi){
    const int n=e.n,axis[3]={k%n,(k/n)%n,k/(n*n)};
    for(int a=0;a<3;++a)xi[a]=b.lo[a]+.5*(b.hi[a]-b.lo[a])*(e.x[axis[a]]+1.0);
}
// A host's moments on the Chebyshev nodes of one cluster holding it, on the far rule of the host.
__global__ void __launch_bounds__(hexThreads) hcaMoments(HexDevice d,HcaNodes e,const HcaBox* boxes,const HcaMomentItem* items,const double* sigma,
        std::size_t first,std::size_t count,double* scratch,std::size_t stride,double* moments,int* invalid){
    const int n1=e.n,q=n1*n1*n1;
    for(std::size_t item=first+blockIdx.x;item<first+count;item+=gridDim.x){
        const HcaMomentItem it=items[item];const HcaBox b=boxes[it.box];const int g=it.host,nC=d.count[g];
        const int N=d.out.n,np=hexCell(d,g)?N*N*N:N*N;
        double* D=scratch+std::size_t(blockIdx.x)*stride;double* L=D+std::size_t(np)*nC;
        hexFarDensities(d,g,nC,D);
        const double* X=d.farX+3*std::size_t(g)*d.farStride;const double* W=d.farW+std::size_t(g)*d.farStride;
        for(int p=threadIdx.x;p<np;p+=blockDim.x)for(int a=0;a<3;++a){const double w=b.hi[a]-b.lo[a],x=w>0.0?2.0*(X[3*p+a]-b.lo[a])/w-1.0:0.0;
            for(int k=0;k<n1;++k)L[(std::size_t(p)*3+a)*n1+k]=hcaLagrange(e,k,x);}
        __syncthreads();
        const std::size_t row=it.start-b.charge;
        for(int m=threadIdx.x;m<nC*q;m+=blockDim.x){const int c=m/q,k=m%q,k0=k%n1,k1=(k/n1)%n1,k2=k/(n1*n1);double sum=0.0;
            for(int p=0;p<np;++p){const double* l=L+std::size_t(p)*3*n1;sum+=W[p]*D[std::size_t(p)*nC+c]*l[k0]*l[n1+k1]*l[2*n1+k2];}
            const double value=sum/sigma[d.first[g]+c];if(!isfinite(value))atomicExch(invalid,1);
            moments[b.table+(row+c)*q+k]=value;}
        __syncthreads();
    }
}
__global__ void __launch_bounds__(hexThreads) hcaProduct(const HcaProductSpan* spans,int q,const double* moments,const double* factors,double* out,std::size_t first,std::size_t count){
    for(std::size_t item=first+blockIdx.x;item<first+count;item+=gridDim.x){const HcaProductSpan s=spans[item];
        for(std::size_t m=threadIdx.x;m<s.rows*s.rank;m+=blockDim.x){const std::size_t i=m%s.rows,l=m/s.rows;double sum=0.0;
            const double* row=moments+s.table+i*q;const double* column=factors+s.factor+std::size_t(q)*l;
            for(int k=0;k<q;++k)sum+=row[k]*column[k];
            out[s.out+m]=sum;}
        __syncthreads();
    }
}
} // namespace radia_hmatrix_cuda
