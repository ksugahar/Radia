/* SPDX-License-Identifier: MPL-2.0 */
// Experimental FP64 AMS execution of an existing NGSolve/Radia hierarchy.
// Host: validation/control only. Device: all vector, matrix and coarse solve work.
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cmath>
#include <algorithm>
#include <climits>
#include <memory>
#include <string>
#include <stdexcept>
#include <vector>
namespace {
constexpr int T=256;
void ck(cudaError_t s) { if(s!=cudaSuccess) throw std::runtime_error(cudaGetErrorString(s)); }
void cb(cublasStatus_t s) { if(s!=CUBLAS_STATUS_SUCCESS) throw std::runtime_error("cuBLAS failed: "+std::to_string(s)); }
const mxArray* field(const mxArray* a,const char* name) {
    if(!mxIsStruct(a)||mxGetNumberOfElements(a)!=1) throw std::runtime_error("scalar struct required");
    auto v=mxGetField(a,0,name); if(!v) throw std::runtime_error(std::string("missing field: ")+name); return v;
}
double scalar(const mxArray* a) {
    if(!mxIsDouble(a)||mxIsComplex(a)||mxGetNumberOfElements(a)!=1||!std::isfinite(mxGetScalar(a)))
        throw std::runtime_error("finite real double metadata scalar required");
    return mxGetScalar(a);
}
int integer(const mxArray* a) { double x=scalar(a); if(x<0||x>INT_MAX||x!=floor(x)) throw std::runtime_error("invalid int32 dimension"); return int(x); }
std::string word(const mxArray* a) {
    if(!mxIsChar(a)) throw std::runtime_error("character metadata required");
    char* p=mxArrayToString(a); if(!p) throw std::runtime_error("invalid text");
    std::string s(p); mxFree(p); return s;
}
struct Input {
    std::vector<const mxGPUArray*> refs;
    ~Input(){ for(auto p:refs) mxGPUDestroyGPUArray(p); }
    template<class V> const V* get(const mxArray* a,mxClassID type,size_t count) {
        if(!mxIsGPUArray(a)) throw std::runtime_error("numerical arrays must be gpuArray; CPU fallback is disabled");
        auto p=mxGPUCreateFromMxArray(a); refs.push_back(p);
        if(mxGPUGetClassID(p)!=type||mxGPUGetComplexity(p)!=mxREAL||mxGPUGetNumberOfElements(p)!=count)
            throw std::runtime_error("GPU array type/length mismatch");
        return static_cast<const V*>(mxGPUGetDataReadOnly(p));
    }
};
struct Buffer {
    double* p=nullptr; size_t count;
    explicit Buffer(size_t n):count(n){ ck(cudaMalloc(&p,std::max<size_t>(1,n)*sizeof(double))); }
    ~Buffer(){if(p)cudaFree(p);}
    Buffer(const Buffer&)=delete; Buffer& operator=(const Buffer&)=delete;
    void zero(){ck(cudaMemset(p,0,count*sizeof(double)));}
};
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
struct Flag {
    int* p=nullptr; explicit Flag(size_t count=1){size_t bytes=std::max<size_t>(1,count)*sizeof(int);ck(cudaMalloc(&p,bytes));ck(cudaMemset(p,0,bytes));}
    ~Flag(){if(p)cudaFree(p);}
    void check(){int value;ck(cudaMemcpy(&value,p,sizeof(int),cudaMemcpyDeviceToHost));if(value)throw std::runtime_error("invalid GPU numeric data or CSR indices");}
};
__global__ void multiply(int n,const int* rp,const int* ci,const double* av,const double* x,double* y) {
    int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=n)return;
    double v=0;for(int j=rp[i];j<rp[i+1];++j)v+=av[j]*x[ci[j]];y[i]=v;
}
__global__ void residual(int n,const double* b,const double* ax,double* r) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)r[i]=b[i]-ax[i];
}
__global__ void smooth_update(int n,const double* r,const double* l1,double* x) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)x[i]+=r[i]/l1[i];
}
__global__ void scaled_add(int n,double alpha,const double* a,double* b) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)b[i]+=alpha*a[i];
}
__global__ void combine(int n,const double* z,double beta,double* p) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)p[i]=z[i]+beta*p[i];
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
__global__ void verify_factor(int n,int used,const int* order,const int* rp,const int* ci,int* seen,int* err) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i<n&&(order[i]<-1||order[i]>=used))atomicExch(err,1);
    if(i<n&&order[i]>=0&&order[i]<used)atomicAdd(seen+order[i],1);
    if(i<used)for(int j=rp[i];j<rp[i+1];++j)
        if(ci[j]<=i||(j>rp[i]&&ci[j]<=ci[j-1]))atomicExch(err,1);
}
__global__ void verify_permutation(int used,const int* seen,int* err){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<used&&seen[i]!=1)atomicExch(err,1);}
struct CSR {
    int n,m,nz;const int *rp,*ci;const double* a;
    CSR(const mxArray* v,Input& in):n(integer(field(v,"rows"))),m(integer(field(v,"cols"))) {
        if(n<1||m<1||n>INT_MAX-T||m>INT_MAX-T)throw std::runtime_error("operator dimensions exceed supported capacity");
        auto vals=field(v,"values");if(!mxIsGPUArray(vals))throw std::runtime_error("CSR values require GPU");
        auto temp=mxGPUCreateFromMxArray(vals);size_t count=mxGPUGetNumberOfElements(temp);mxGPUDestroyGPUArray(temp);
        if(count>INT_MAX-T)throw std::runtime_error("CSR index capacity exceeded");nz=int(count);
        rp=in.get<int>(field(v,"rowptr"),mxINT32_CLASS,size_t(n)+1);
        ci=in.get<int>(field(v,"columns"),mxINT32_CLASS,nz);a=in.get<double>(vals,mxDOUBLE_CLASS,nz);
        Flag flag;verify_csr<<<(std::max(n,nz)+T-1)/T,T>>>(n,m,nz,rp,ci,a,flag.p);flag.check();
    }
    void mult(const double* x,double* y)const{multiply<<<(n+T-1)/T,T>>>(n,rp,ci,a,x,y);}
};
void finite_vec(Input& in,const mxArray* a,int n,const double*& result,bool positive=false){
    result=in.get<double>(a,mxDOUBLE_CLASS,n);Flag f;verify_vec<<<(n+T-1)/T,T>>>(n,result,positive,f.p);f.check();
}
struct Coarse {
    std::string kind;int n,used=0;const double *a=nullptr,*di=nullptr;const int *order=nullptr,*rp=nullptr,*ci=nullptr;
    std::unique_ptr<Buffer> work;
    Coarse(const mxArray* v,Input& in,int rows):kind(word(field(v,"kind"))),n(rows){
        if(integer(field(v,"rows"))!=n)throw std::runtime_error("coarse size mismatch");
        if(kind=="dense_inverse") {
            a=in.get<double>(field(v,"values"),mxDOUBLE_CLASS,size_t(n)*n);
            if(size_t(n)*n>INT_MAX-T)throw std::runtime_error("dense coarse capacity exceeded");
            Flag f;verify_vec<<<(n*n+T-1)/T,T>>>(n*n,a,false,f.p);f.check();
        } else if(kind=="ngsolve_sparse_ldlt") {
            auto d=field(v,"diagonal_inverse");if(!mxIsGPUArray(d))throw std::runtime_error("coarse diagonal requires GPU");
            auto tmp=mxGPUCreateFromMxArray(d);size_t count=mxGPUGetNumberOfElements(tmp);mxGPUDestroyGPUArray(tmp);
            if(count>size_t(n))throw std::runtime_error("invalid used coarse rows");used=int(count);
            di=in.get<double>(d,mxDOUBLE_CLASS,used);order=in.get<int>(field(v,"order"),mxINT32_CLASS,n);
            auto val=field(v,"values");if(!mxIsGPUArray(val))throw std::runtime_error("coarse factors require GPU arrays");
            tmp=mxGPUCreateFromMxArray(val);count=mxGPUGetNumberOfElements(tmp);mxGPUDestroyGPUArray(tmp);
            if(count>INT_MAX-T)throw std::runtime_error("coarse capacity exceeded");
            rp=in.get<int>(field(v,"rowptr"),mxINT32_CLASS,size_t(used)+1);
            ci=in.get<int>(field(v,"columns"),mxINT32_CLASS,count);a=in.get<double>(val,mxDOUBLE_CLASS,count);
            Flag f,seen(used);verify_csr<<<std::max(1,(std::max(used,int(count))+T-1)/T),T>>>(used,used,int(count),rp,ci,a,f.p);f.check();
            verify_factor<<<(n+T-1)/T,T>>>(n,used,order,rp,ci,seen.p,f.p);
            verify_permutation<<<std::max(1,(used+T-1)/T),T>>>(used,seen.p,f.p);
            verify_vec<<<std::max(1,(used+T-1)/T),T>>>(used,di,false,f.p);f.check();work=std::make_unique<Buffer>(used);
        } else if(kind!="ten_l1_jacobi_sweeps")throw std::runtime_error("unknown coarse solver");
    }
    void apply(const double* b,double* x)const{
        if(kind=="dense_inverse")dense_coarse<<<(n+T-1)/T,T>>>(n,a,b,x);
        else sparse_coarse<<<1,32>>>(n,used,order,rp,ci,a,di,b,x,work->p);
    }
};
struct Level {
    CSR a;std::unique_ptr<CSR> p,pt;const double* l1;int smooth;
    Buffer x,b,r,tmp;std::unique_ptr<Coarse> coarse;
    Level(const mxArray* v,Input& in,bool last):a(field(v,"A"),in),smooth(integer(field(v,"num_smooth"))),x(a.n),b(a.n),r(a.n),tmp(a.n){
        if(a.n!=a.m||smooth<1)throw std::runtime_error("invalid AMG level");finite_vec(in,field(v,"l1"),a.n,l1,true);
        if(last)coarse=std::make_unique<Coarse>(field(v,"coarse"),in,a.n);
        else {p=std::make_unique<CSR>(field(v,"P"),in);pt=std::make_unique<CSR>(field(v,"Pt"),in);}
    }
    void relax(){a.mult(x.p,tmp.p);residual<<<(a.n+T-1)/T,T>>>(a.n,b.p,tmp.p,r.p);smooth_update<<<(a.n+T-1)/T,T>>>(a.n,r.p,l1,x.p);}
};
struct Space {
    CSR p,pt;std::vector<std::unique_ptr<Level>> levels;Buffer rhs,out;
    Space(const mxArray* v,Input& in):p(field(v,"P"),in),pt(field(v,"Pt"),in),rhs(p.m),out(p.n){
        auto list=field(v,"levels");if(!mxIsCell(list)||mxGetNumberOfElements(list)<1)throw std::runtime_error("AMG levels require nonempty cell array");
        size_t count=mxGetNumberOfElements(list);
        for(size_t k=0;k<count;++k)levels.push_back(std::make_unique<Level>(mxGetCell(list,k),in,k+1==count));
        if(pt.n!=p.m||pt.m!=p.n||levels[0]->a.n!=p.m)throw std::runtime_error("auxiliary dimensions mismatch");
        for(size_t k=0;k+1<count;++k){auto& l=*levels[k];int nc=levels[k+1]->a.n;
            if(l.p->n!=l.a.n||l.p->m!=nc||l.pt->n!=nc||l.pt->m!=l.a.n)throw std::runtime_error("AMG transfer dimensions mismatch");}
    }
    void cycle(size_t k){
        auto& l=*levels[k];l.x.zero();
        if(l.coarse){if(l.coarse->kind=="ten_l1_jacobi_sweeps")for(int s=0;s<10;++s)l.relax();else l.coarse->apply(l.b.p,l.x.p);return;}
        for(int s=0;s<l.smooth;++s)l.relax();
        l.a.mult(l.x.p,l.tmp.p);residual<<<(l.a.n+T-1)/T,T>>>(l.a.n,l.b.p,l.tmp.p,l.r.p);
        l.pt->mult(l.r.p,levels[k+1]->b.p);cycle(k+1);l.p->mult(levels[k+1]->x.p,l.tmp.p);
        scaled_add<<<(l.a.n+T-1)/T,T>>>(l.a.n,1,l.tmp.p,l.x.p);
        for(int s=0;s<l.smooth;++s)l.relax();
    }
    void apply(const double* residual){pt.mult(residual,levels[0]->b.p);cycle(0);p.mult(levels[0]->x.p,out.p);}
};
struct AMS {
    CSR a;const double* l1;const int* free;int smooth,cycle,gradient;double omega;
    std::vector<std::unique_ptr<Space>> spaces;Buffer b,x,r,tmp;
    AMS(const mxArray* v,Input& in):a(field(v,"A"),in),smooth(integer(field(v,"num_smooth"))),cycle(integer(field(v,"cycle_type"))),
        gradient(0),omega(scalar(field(v,"correction_weight"))),b(a.n),x(a.n),r(a.n),tmp(a.n){
        if(word(field(v,"schema"))!="radia.ams-accelerator.v1"||a.n!=a.m||smooth<1||(cycle!=1&&cycle!=7))throw std::runtime_error("unsupported AMS snapshot");
        auto beta=field(v,"beta_zero");if(!mxIsLogicalScalar(beta))throw std::runtime_error("logical beta_zero required");gradient=!mxIsLogicalScalarTrue(beta);
        finite_vec(in,field(v,"l1"),a.n,l1,true);free=in.get<int>(field(v,"free"),mxINT32_CLASS,a.n);
        auto list=field(v,"spaces");if(!mxIsCell(list)||mxGetNumberOfElements(list)!=size_t(3+gradient))throw std::runtime_error("invalid AMS spaces");
        for(size_t k=0;k<mxGetNumberOfElements(list);++k){spaces.push_back(std::make_unique<Space>(mxGetCell(list,k),in));
            if(spaces.back()->p.n!=a.n)throw std::runtime_error("AMS interpolation dimension mismatch");}
    }
    void res(){a.mult(x.p,tmp.p);residual<<<(a.n+T-1)/T,T>>>(a.n,b.p,tmp.p,r.p);}
    void relax(){for(int s=0;s<smooth;++s){res();smooth_update<<<(a.n+T-1)/T,T>>>(a.n,r.p,l1,x.p);}}
    void grad(){if(!gradient)return;res();spaces[0]->apply(r.p);scaled_add<<<(a.n+T-1)/T,T>>>(a.n,omega,spaces[0]->out.p,x.p);}
    void nodal(){res();for(int k=gradient;k<gradient+3;++k){spaces[k]->apply(r.p);scaled_add<<<(a.n+T-1)/T,T>>>(a.n,omega,spaces[k]->out.p,x.p);}}
    void apply(const double* rhs,double* output){ck(cudaMemcpy(b.p,rhs,a.n*sizeof(double),cudaMemcpyDeviceToDevice));mask<<<(a.n+T-1)/T,T>>>(a.n,free,b.p);x.zero();
        if(cycle==1){relax();grad();nodal();grad();relax();}
        else {relax();grad();relax();nodal();relax();grad();relax();}
        mask<<<(a.n+T-1)/T,T>>>(a.n,free,x.p);ck(cudaMemcpy(output,x.p,a.n*sizeof(double),cudaMemcpyDeviceToDevice));}
};
struct Blas {cublasHandle_t h;Blas(){cb(cublasCreate(&h));}~Blas(){cublasDestroy(h);}double dot(int n,const double* a,const double* b){double d;cb(cublasDdot(h,n,a,1,b,1,&d));return d;}double norm(int n,const double* a){double d;cb(cublasDnrm2(h,n,a,1,&d));return d;}};
struct Output {mxGPUArray* a;explicit Output(int n){mwSize s=n;a=mxGPUCreateGPUArray(1,&s,mxDOUBLE_CLASS,mxREAL,MX_GPU_INITIALIZE_VALUES);}~Output(){mxGPUDestroyGPUArray(a);}double* data(){return static_cast<double*>(mxGPUGetData(a));}};
void run(int nlhs,mxArray** out,int nrhs,const mxArray** args){
    if(nrhs!=5||nlhs!=2)throw std::runtime_error("[x,info]=radia_cuda_ams(snapshot,b,tol,maxiter,apply_only)");
    double tol=scalar(args[2]);int maxiter=integer(args[3]);
    if(!(tol>0&&tol<=1e-6)||maxiter<1||maxiter==INT_MAX||!mxIsLogicalScalar(args[4]))throw std::runtime_error("invalid tolerance, iteration limit or apply flag");
    mxInitGPU();size_t memory_before,total,available;ck(cudaMemGetInfo(&memory_before,&total));
    Input in;AMS ams(args[0],in);int n=ams.a.n;const double* b;
    ck(cudaMemGetInfo(&available,&total));size_t workspace_sample=memory_before>available?memory_before-available:0;
    finite_vec(in,args[1],n,b);bool apply_only=mxIsLogicalScalarTrue(args[4]);
    Buffer zero(n);zero.zero();Flag f;verify_mask<<<(n+T-1)/T,T>>>(n,ams.free,apply_only?zero.p:b,f.p);f.check();
    Output solution(n);double rel=0;int it=0;
    if(apply_only)ams.apply(b,solution.data());
    else {
        Buffer r(n),p(n),z(n),ap(n);Blas blas;
        ck(cudaMemGetInfo(&available,&total));workspace_sample=std::max(workspace_sample,memory_before>available?memory_before-available:0);
        ck(cudaMemcpy(r.p,b,n*sizeof(double),cudaMemcpyDeviceToDevice));
        double bnorm=blas.norm(n,b);if(!std::isfinite(bnorm)||bnorm<0)throw std::runtime_error("invalid RHS norm");
        if(bnorm>0){ams.apply(r.p,z.p);ck(cudaMemcpy(p.p,z.p,n*sizeof(double),cudaMemcpyDeviceToDevice));double rho=blas.dot(n,r.p,z.p);rel=1;
            for(it=1;it<=maxiter;++it){ams.a.mult(p.p,ap.p);double pap=blas.dot(n,p.p,ap.p);
                if(!(pap>0)||!(rho>0)||!std::isfinite(pap)||!std::isfinite(rho))throw std::runtime_error("AMS PCG breakdown: nonpositive curvature");
                double alpha=rho/pap;scaled_add<<<(n+T-1)/T,T>>>(n,alpha,p.p,solution.data());scaled_add<<<(n+T-1)/T,T>>>(n,-alpha,ap.p,r.p);
                rel=blas.norm(n,r.p)/bnorm;if(!std::isfinite(rel))throw std::runtime_error("nonfinite residual");
                if(rel<=tol){ams.a.mult(solution.data(),ap.p);residual<<<(n+T-1)/T,T>>>(n,b,ap.p,r.p);rel=blas.norm(n,r.p)/bnorm;if(rel<=tol)break;
                    ams.apply(r.p,z.p);ck(cudaMemcpy(p.p,z.p,n*sizeof(double),cudaMemcpyDeviceToDevice));rho=blas.dot(n,r.p,z.p);continue;}
                ams.apply(r.p,z.p);double next=blas.dot(n,r.p,z.p);combine<<<(n+T-1)/T,T>>>(n,z.p,next/rho,p.p);rho=next;
            }
            ams.a.mult(solution.data(),ap.p);residual<<<(n+T-1)/T,T>>>(n,b,ap.p,r.p);rel=blas.norm(n,r.p)/bnorm;
            if(!std::isfinite(rel)||rel>tol||it>maxiter)throw std::runtime_error("AMS PCG failed its FP64 true residual gate");
        }
    }
    ck(cudaGetLastError());ck(cudaDeviceSynchronize());out[0]=mxGPUCreateMxArrayOnGPU(solution.a);
    const char* names[]={"iterations","true_relative_residual","coarse_execution","setup_execution","workspace_sample_bytes"};
    out[1]=mxCreateStructMatrix(1,1,5,names);mxSetField(out[1],0,names[0],mxCreateDoubleScalar(it));mxSetField(out[1],0,names[1],mxCreateDoubleScalar(rel));
    mxSetField(out[1],0,names[2],mxCreateString("GPU: existing sparse LDL factors with warp-parallel row products, or existing dense inverse"));
    mxSetField(out[1],0,names[3],mxCreateString("CPU Python/NGSolve: FE assembly, AMS coarsening, interpolation, Galerkin products and factorization"));
    mxSetField(out[1],0,names[4],mxCreateDoubleScalar(double(workspace_sample)));
}
}
void mexFunction(int nlhs,mxArray** out,int nrhs,const mxArray** args){try{run(nlhs,out,nrhs,args);}catch(const std::exception& e){mexErrMsgIdAndTxt("radia:cuda:ams","%s",e.what());}}
