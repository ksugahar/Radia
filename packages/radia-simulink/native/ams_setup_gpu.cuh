/* SPDX-License-Identifier: MPL-2.0 */
// Included by radia_cuda_ams.cu. Deterministic row-owned P^T A P generation.
#include "ams_setup_gpu.hpp"
#include <chrono>
#include <limits>
namespace setup_gpu {
struct Context { AMSSetupGPUStats stats; size_t live=0; const void* identity=nullptr; };
struct Allocation {
    Context& c; void* p=nullptr; size_t bytes;
    Allocation(Context& ctx,size_t n):c(ctx),bytes(n) {
        size_t available,total;ck(cudaMemGetInfo(&available,&total));
        const size_t reserve=std::max<size_t>(size_t(2)<<30,total/5);
        if(n>available||available-n<reserve)throw std::runtime_error("AMS CUDA setup exceeds GPU memory budget");
        ck(cudaMalloc(&p,std::max<size_t>(n,1)));c.live+=n;c.stats.peak_owned_bytes=std::max(c.stats.peak_owned_bytes,c.live);
    }
    ~Allocation(){if(p)cudaFree(p);c.live-=bytes;}
    Allocation(const Allocation&)=delete;
    template<class V> V* as(){return static_cast<V*>(p);}
};
struct Matrix {
    int rows,cols,nnz;Allocation rp,ci,av;
    Matrix(Context& c,int n,int m,int nz):rows(n),cols(m),nnz(nz),
      rp(c,(size_t(rows)+1)*sizeof(int)),ci(c,size_t(nnz)*sizeof(int)),av(c,size_t(nnz)*sizeof(double)) {}
    Matrix(Context& c,const AMSSetupCSR& h):rows(h.rows),cols(h.cols),nnz(h.nnz),
      rp(c,(size_t(rows)+1)*sizeof(int)),ci(c,size_t(nnz)*sizeof(int)),av(c,size_t(nnz)*sizeof(double)) {
        ck(cudaMemcpy(rp.p,h.rowptr,rp.bytes,cudaMemcpyHostToDevice));
        ck(cudaMemcpy(ci.p,h.columns,ci.bytes,cudaMemcpyHostToDevice));
        ck(cudaMemcpy(av.p,h.values,av.bytes,cudaMemcpyHostToDevice));
        c.stats.upload_bytes+=rp.bytes+ci.bytes+av.bytes;
    }
};
struct Session:Context {std::unique_ptr<Matrix> a;};
__global__ void multiply_rows(int n,int cap,const int* ar,const int* ac,const double* av,
    const int* br,const int* bc,const double* bv,int* keys,double* vals,int* counts,int* failed) {
    int row=blockIdx.x*blockDim.x+threadIdx.x;if(row>=n)return;
    size_t offset=size_t(row)*cap;int* key=keys+offset;double* val=vals+offset;
    int count=0;
    // Match NGSolve's two sorted CSR products: AP followed by Pt(AP).
    // Each output row has one owner, with no unordered floating-point atomics.
    for(int e=ar[row];e<ar[row+1];++e) {
        int i=ac[e];double a=av[e];
        for(int q=br[i];q<br[i+1];++q) {
            int col=bc[q];unsigned slot=(unsigned(col)*2654435761u)&unsigned(cap-1);
            int step=0;
            while(key[slot]!=-1&&key[slot]!=col&&step<cap){slot=(slot+1)&(cap-1);++step;}
            if(step==cap){atomicExch(failed,1);return;}
            if(key[slot]==-1){key[slot]=col;val[slot]=0;count++;}
            val[slot]=__dadd_rn(val[slot],__dmul_rn(a,bv[q]));
        }
    }
    // Pack then sort columns without changing any accumulated value.
    int used=0;
    for(int j=0;j<cap;++j)if(key[j]>=0){key[used]=key[j];val[used]=val[j];++used;}
    for(int j=1;j<used;++j){int col=key[j];double v=val[j];int k=j-1;
        while(k>=0&&key[k]>col){key[k+1]=key[k];val[k+1]=val[k];--k;}key[k+1]=col;val[k+1]=v;}
    for(int j=0;j<used;++j)if(!isfinite(val[j]))atomicExch(failed,2);
    counts[row]=count;
}
__global__ void compact_rows(int n,int cap,const int* rp,const int* keys,const double* vals,int* ci,double* av) {
    int row=blockIdx.x;if(row>=n)return;
    for(int k=threadIdx.x;k<rp[row+1]-rp[row];k+=blockDim.x){ci[rp[row]+k]=keys[size_t(row)*cap+k];av[rp[row]+k]=vals[size_t(row)*cap+k];}
}
std::unique_ptr<Matrix> multiply(Context& c,Matrix& left,Matrix& right) {
    if(left.cols!=right.rows)throw std::runtime_error("CUDA CSR multiply dimensions mismatch");
    Allocation counts(c,size_t(left.rows)*sizeof(int)),failed(c,sizeof(int));
    int cap=128;
    for(;;) {
        size_t slots=size_t(left.rows)*size_t(cap);
        if(slots>std::numeric_limits<size_t>::max()/sizeof(double))throw std::overflow_error("AMS CUDA hash capacity overflow");
        Allocation keys(c,slots*sizeof(int)),vals(c,slots*sizeof(double));
        ck(cudaMemset(keys.p,0xff,keys.bytes));ck(cudaMemset(failed.p,0,sizeof(int)));
        multiply_rows<<<(left.rows+127)/128,128>>>(left.rows,cap,left.rp.as<int>(),left.ci.as<int>(),left.av.as<double>(),
            right.rp.as<int>(),right.ci.as<int>(),right.av.as<double>(),keys.as<int>(),vals.as<double>(),counts.as<int>(),failed.as<int>());
        ck(cudaGetLastError());int status;ck(cudaMemcpy(&status,failed.p,sizeof(int),cudaMemcpyDeviceToHost));
        if(status==2)throw std::runtime_error("Nonfinite CUDA Galerkin product");
        if(status){if(cap>INT_MAX/2)throw std::overflow_error("AMS CUDA row hash capacity overflow");cap*=2;continue;}
        std::vector<int> sizes(left.rows),rp(size_t(left.rows)+1);
        ck(cudaMemcpy(sizes.data(),counts.p,counts.bytes,cudaMemcpyDeviceToHost));c.stats.download_bytes+=counts.bytes;
        size_t nz=0;
        for(int i=0;i<left.rows;++i){rp[i]=int(nz);nz+=sizes[i];if(nz>size_t(INT_MAX-256))throw std::overflow_error("AMS CUDA product exceeds CSR capacity");}
        rp[left.rows]=int(nz);
        auto out=std::make_unique<Matrix>(c,left.rows,right.cols,int(nz));
        ck(cudaMemcpy(out->rp.p,rp.data(),out->rp.bytes,cudaMemcpyHostToDevice));c.stats.upload_bytes+=out->rp.bytes;
        compact_rows<<<left.rows,128>>>(left.rows,cap,out->rp.as<int>(),keys.as<int>(),vals.as<double>(),out->ci.as<int>(),out->av.as<double>());
        ck(cudaGetLastError());
        return out;
    }
}
}
void* AMSSetupGPUCreate(){mxInitGPU();return new setup_gpu::Session;}
void AMSSetupGPUDestroy(void* p) noexcept {delete static_cast<setup_gpu::Session*>(p);}
AMSSetupGPUStats AMSSetupGPUStatistics(void* p){return static_cast<setup_gpu::Session*>(p)->stats;}
AMSSetupProduct AMSSetupGPURestrict(void* context,const void* identity,const AMSSetupCSR& a,
    const AMSSetupCSR& p,const AMSSetupCSR& pt) {
    using namespace setup_gpu;
    auto& c=*static_cast<Session*>(context);auto start=std::chrono::steady_clock::now();
    if(a.rows!=a.cols||p.rows!=a.rows||pt.rows!=p.cols||pt.cols!=p.rows)
        throw std::runtime_error("AMS CUDA Galerkin dimensions mismatch");
    if(c.identity!=identity){c.a.reset();c.a=std::make_unique<Matrix>(c,a);c.identity=identity;}
    Matrix prolong(c,p),transpose(c,pt);
    auto ap=multiply(c,*c.a,prolong);
    auto product=multiply(c,transpose,*ap);
    AMSSetupProduct out;out.rowptr.resize(size_t(product->rows)+1);out.columns.resize(product->nnz);out.values.resize(product->nnz);
    ck(cudaMemcpy(out.rowptr.data(),product->rp.p,product->rp.bytes,cudaMemcpyDeviceToHost));
    ck(cudaMemcpy(out.columns.data(),product->ci.p,product->ci.bytes,cudaMemcpyDeviceToHost));
    ck(cudaMemcpy(out.values.data(),product->av.p,product->av.bytes,cudaMemcpyDeviceToHost));
    c.stats.download_bytes+=product->rp.bytes+product->ci.bytes+product->av.bytes;
    c.stats.products++;c.stats.seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return out;
}
