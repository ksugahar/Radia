/* SPDX-License-Identifier: MPL-2.0 */
#pragma once
// Shared AMS/CompactAMG device primitives. No MATLAB or Python ABI.
// Include inside the adapter's kernel namespace after CUDA headers.
// Keep row reduction, smoothing, coarse factor and warp dependency order intact.
__global__ void verify_csr(int n,int m,int nz,const int* rp,const int* ci,const double* a,int* err) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i<n && (rp[i]<0||rp[i+1]<rp[i]||rp[i+1]>nz)) atomicExch(err,1);
    if(i<nz && (ci[i]<0||ci[i]>=m||!isfinite(a[i]))) atomicExch(err,1);
    if(i==0 && (rp[0]!=0||rp[n]!=nz)) atomicExch(err,1);
}
__global__ void verify_vec(int n,const double* a,bool positive,int* err) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i<n&&(!isfinite(a[i])||(positive&&a[i]<=0))) atomicExch(err,1);
}
__global__ void verify_mask(int n,const int* free,const double* b,int* err) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i<n&&((free[i]!=0&&free[i]!=1)||(!free[i]&&b[i]!=0))) atomicExch(err,1);
}
__global__ void multiply(int n,const int* rp,const int* ci,const double* av,const double* x,double* y) {
    int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=n)return;
    double v=0;for(int j=rp[i];j<rp[i+1];++j)v+=av[j]*x[ci[j]];y[i]=v;
}
// Keep the CSR reduction order; only combine independent row-local work.
__global__ void multiply_residual(int n,const int* rp,const int* ci,const double* av,
    const double* x,const double* b,double* r) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
    double value=0;for(int j=rp[i];j<rp[i+1];++j)value+=av[j]*x[ci[j]];
    r[i]=b[i]-value;
}
__global__ void initial_relax(int n,const double* b,const double* l1,double* x) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)x[i]=b[i]/l1[i];
}
// Read the old iterate and write a separate buffer: in-place CSR Jacobi races.
__global__ void multiply_relax(int n,const int* rp,const int* ci,const double* av,
    const double* x,const double* b,const double* l1,double* next) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
    double value=0;for(int j=rp[i];j<rp[i+1];++j)value+=av[j]*x[ci[j]];
    next[i]=x[i]+(b[i]-value)/l1[i];
}
__global__ void multiply_add(int n,const int* rp,const int* ci,const double* av,
    const double* x,double alpha,double* y) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
    double value=0;for(int j=rp[i];j<rp[i+1];++j)value+=av[j]*x[ci[j]];
    y[i]+=alpha*value;
}
__global__ void pcg_update(int n,double alpha,const double* p,const double* ap,double* x,double* r) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n){x[i]+=alpha*p[i];r[i]-=alpha*ap[i];}
}
__global__ void pcg_update_device(int n,const double* scalars,const double* p,const double* ap,double* x,double* r) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
    double rho=scalars[0],pap=scalars[1];
    // The immediately following host residual norm reports breakdown loudly.
    if(!(rho>0)||!(pap>0)||!isfinite(rho)||!isfinite(pap)){r[i]=nan("");return;}
    double alpha=rho/pap;x[i]+=alpha*p[i];r[i]-=alpha*ap[i];
}
__global__ void combine(int n,const double* z,const double* scalars,double* p) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)p[i]=z[i]+(scalars[2]/scalars[0])*p[i];
}
__global__ void mask(int n,const int* free,double* x) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n&&!free[i])x[i]=0;
}
__global__ void dense_coarse(int n,const double* a,const double* b,double* x) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
    double s=0;for(int j=0;j<n;++j)s+=a[size_t(i)*n+j]*b[j];x[i]=s;
}
// Preserve NGSolve's sparse LDL factor and its dependency order. A warp
// shares each row's products; row dependencies remain ordered on device.
__global__ void sparse_coarse(int n,int used,const int* order,const int* rp,const int* ci,
    const double* a,const double* dinv,const double* b,double* x,double* w) {
    int lane=threadIdx.x;
    for(int i=lane;i<used;i+=32)w[i]=0;
    __syncwarp();
    for(int i=lane;i<n;i+=32)if(order[i]>=0)w[order[i]]=b[i];
    __syncwarp();
    for(int i=0;i<used;++i){double value=w[i];
        for(int j=rp[i]+lane;j<rp[i+1];j+=32)w[ci[j]]-=a[j]*value;
        __syncwarp();
    }
    for(int i=lane;i<used;i+=32)w[i]*=dinv[i];
    __syncwarp();
    for(int i=used-1;i>=0;--i){double sum=0;
        for(int j=rp[i]+lane;j<rp[i+1];j+=32)sum+=a[j]*w[ci[j]];
        for(int offset=16;offset;offset/=2)sum+=__shfl_down_sync(0xffffffff,sum,offset);
        if(lane==0)w[i]-=sum;
        __syncwarp();
    }
    for(int i=lane;i<n;i+=32)x[i]=order[i]<0?0:w[order[i]];
}
__global__ void verify_factor(int n,int used,int nz,const int* order,const int* rp,const int* ci,int* seen,int* err) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i<n&&(order[i]<-1||order[i]>=used))atomicExch(err,1);
    if(i<n&&order[i]>=0&&order[i]<used)atomicAdd(seen+order[i],1);
    if(i<used){int first=rp[i],last=rp[i+1];
        // This kernel must be safe even while the independent CSR check fails.
        if(first<0||last<first||last>nz){atomicExch(err,1);return;}
        for(int j=first;j<last;++j)
            if(ci[j]<=i||(j>first&&ci[j]<=ci[j-1]))atomicExch(err,1);
    }
}
__global__ void verify_permutation(int used,const int* seen,int* err){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<used&&seen[i]!=1)atomicExch(err,1);}
// Cache the original sparse LDLt action in Radia-owned row-major storage.
// Each warp computes one column using the original factor dependency order.
__global__ void inverse_sparse_coarse(int used,const int* rp,const int* ci,const double* a,const double* di,double* inverse){
    int column=blockIdx.x,lane=threadIdx.x;double* w=inverse+column;
    for(int i=lane;i<used;i+=32)w[size_t(i)*used]=i==column?1.:0.;__syncwarp();
    for(int i=0;i<used;++i){double value=w[size_t(i)*used];for(int j=rp[i]+lane;j<rp[i+1];j+=32)w[size_t(ci[j])*used]-=a[j]*value;__syncwarp();}
    for(int i=lane;i<used;i+=32)w[size_t(i)*used]*=di[i];__syncwarp();
    for(int i=used-1;i>=0;--i){double sum=0;for(int j=rp[i]+lane;j<rp[i+1];j+=32)sum+=a[j]*w[size_t(ci[j])*used];
        for(int offset=16;offset;offset/=2)sum+=__shfl_down_sync(0xffffffff,sum,offset);
        if(lane==0)w[size_t(i)*used]-=sum;__syncwarp();}
}
__global__ void cache_inverse_order(int n,const int* order,int* reverse){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n&&order[i]>=0)reverse[order[i]]=i;}
__global__ void cached_sparse_coarse(int n,int used,const int* order,const int* reverse,const double* inverse,const double* b,double* x){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;int row=order[i];double value=0;
    if(row>=0)for(int j=0;j<used;++j)value+=inverse[size_t(row)*used+j]*b[reverse[j]];x[i]=value;
}
