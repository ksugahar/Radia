#include "mex_error_boundary.hpp"
/* SPDX-License-Identifier: MPL-2.0 */
// FP64 natural-order IC(0). Host performs symbolic indexing only; CUDA owns
// numerical factorization, triangular solves, PCG and projected-operator residuals.
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
constexpr int T=256;
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point start){return std::chrono::duration<double>(Clock::now()-start).count();}
void ck(cudaError_t s){if(s!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(s));}
void cs(cusparseStatus_t s){if(s!=CUSPARSE_STATUS_SUCCESS)throw std::runtime_error("cuSPARSE error "+std::to_string(s));}
void cb(cublasStatus_t s){if(s!=CUBLAS_STATUS_SUCCESS)throw std::runtime_error("cuBLAS error "+std::to_string(s));}
template<class V> struct Buffer {
    V* p=nullptr;size_t n=0;
    ~Buffer(){if(p)cudaFree(p);}
    void alloc(size_t count){if(p)throw std::runtime_error("buffer already allocated");if(count>SIZE_MAX/sizeof(V))throw std::overflow_error("allocation overflow");n=count;ck(cudaMalloc(&p,std::max<size_t>(1,n)*sizeof(V)));}
    void upload(const std::vector<V>& v){alloc(v.size());ck(cudaMemcpy(p,v.data(),n*sizeof(V),cudaMemcpyHostToDevice));}
    size_t bytes()const{return n*sizeof(V);}
    Buffer()=default;Buffer(const Buffer&)=delete;
};
struct Inputs {
    std::vector<const mxGPUArray*> a;
    ~Inputs(){for(auto p:a)mxGPUDestroyGPUArray(p);}
    template<class V> const V* get(const mxArray* value,mxClassID type,size_t count){
        if(!mxIsGPUArray(value))throw std::runtime_error("gpuArray required; no CPU solver fallback");
        auto p=mxGPUCreateFromMxArray(value);a.push_back(p);
        if(mxGPUGetClassID(p)!=type||mxGPUGetComplexity(p)!=mxREAL||mxGPUGetNumberOfElements(p)!=count)throw std::runtime_error("GPU type or dimension mismatch");
        return static_cast<const V*>(mxGPUGetDataReadOnly(p));
    }
};
const mxArray* field(const mxArray* a,const char* name){if(!mxIsStruct(a)||mxGetNumberOfElements(a)!=1)throw std::runtime_error("scalar CSR struct required");auto p=mxGetField(a,0,name);if(!p)throw std::runtime_error(std::string("missing ")+name);return p;}
double scalar(const mxArray* a){if(!mxIsDouble(a)||mxIsComplex(a)||mxGetNumberOfElements(a)!=1||!std::isfinite(mxGetScalar(a)))throw std::runtime_error("finite double scalar required");return mxGetScalar(a);}
int dimension(const mxArray* a){double v=scalar(a);if(v<1||v>INT_MAX||v!=floor(v))throw std::runtime_error("positive int32 dimension required");return int(v);}
int blocks(size_t n){return int((n+T-1)/T);}
__global__ void validate(int n,const int* rp,const int* ci,const int* free,const double* av,int* failure){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
    if(free[i]!=0&&free[i]!=1){atomicExch(failure,1);return;}
    for(int p=rp[i];p<rp[i+1];++p){int j=ci[p];double value=av[p];if(!isfinite(value)){atomicExch(failure,1);continue;}
        int lo=rp[j],hi=rp[j+1];while(lo<hi){int m=lo+(hi-lo)/2;if(ci[m]<i)lo=m+1;else hi=m;}
        if(lo>=rp[j+1]||ci[lo]!=i||!isfinite(av[lo])||fabs(value-av[lo])>64*2.220446049250313e-16*fmax(1.,fmax(fabs(value),fabs(av[lo]))))atomicExch(failure,1);
        if(i==j&&free[i]&&!(value>0))atomicExch(failure,1);
        if((!free[i]||!free[j])&&((i!=j&&value!=0)||(i==j&&value!=0&&value!=1)))atomicExch(failure,1);
    }
}
__global__ void identityRows(int n,const int* rp,const int* ci,const int* free,double* av){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n&&!free[i])for(int p=rp[i];p<rp[i+1];++p)if(ci[p]==i)av[p]=1.;}
__global__ void lowerValues(int nz,const int* source,const double* original,double* lower){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<nz)lower[i]=original[source[i]];}
__global__ void factorLevel(int count,const int* rows,const int* rp,const int* ci,double* ld,double shift,int* failure){
    int t=blockIdx.x*blockDim.x+threadIdx.x;if(t>=count)return;int i=rows[t],begin=rp[i],diag=rp[i+1]-1;
    for(int p=begin;p<diag;++p){int j=ci[p],jdiag=rp[j+1]-1;double sum=ld[p];int a=begin,b=rp[j];
        while(a<p&&b<jdiag){int ka=ci[a],kb=ci[b];if(ka<kb)++a;else if(kb<ka)++b;else{sum-=ld[a]*ld[b]/ld[rp[ka+1]-1];++a;++b;}}
        ld[p]=sum;
    }
    double original=ld[diag],value=shift*original;
    for(int p=begin;p<diag;++p)value-=ld[p]*ld[p]/ld[rp[ci[p]+1]-1];
    if(!isfinite(value)||!(value>0)||value<1e-14*original){atomicExch(failure,1);ld[diag]=1.;return;}
    ld[diag]=value;
}
__global__ void toCholesky(int n,const int* rp,const int* ci,const double* ld,double* l){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;for(int p=rp[i];p<rp[i+1];++p)l[p]=(ci[p]==i)?sqrt(ld[p]):ld[p]/sqrt(ld[rp[ci[p]+1]-1]);}
__global__ void checkVector(int n,const double* x,const int* free,int* failure){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n){if(!isfinite(x[i])||(!free[i]&&x[i]!=0))atomicExch(failure,1);if(x[i]!=0)atomicExch(failure+1,1);}}
__global__ void spmv(int n,const int* rp,const int* ci,const double* av,const double* x,double* y){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n){double s=0;for(int p=rp[i];p<rp[i+1];++p)s+=av[p]*x[ci[p]];y[i]=s;}}
__global__ void update(int n,double alpha,const double* p,const double* ap,double* x,double* r){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n){x[i]+=alpha*p[i];r[i]-=alpha*ap[i];}}
__global__ void direction(int n,double beta,const double* z,double* p){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)p[i]=z[i]+beta*p[i];}
__global__ void subtract(int n,const double* b,double* ax,double bnorm){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n){double residual=b[i]-ax[i];ax[i]=bnorm>0?residual/bnorm:residual;}}

struct IterationState {double rho,alpha,beta,bnorm,tolerance;int iteration,limit,active,error;};
constexpr int GRAPH_BATCH=8;
__global__ void graphAlpha(IterationState* s,const double* scalars){
    if(!s->active)return;
    double pap=scalars[0];
    if(!isfinite(pap)||!(pap>0)||!isfinite(s->rho)||!(s->rho>0)){s->error=1;s->active=0;return;}
    s->alpha=s->rho/pap;if(!isfinite(s->alpha)){s->error=2;s->active=0;}
}
__global__ void graphUpdate(int n,const IterationState* s,const double* p,const double* ap,double* x,double* r){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n&&s->active){x[i]+=s->alpha*p[i];r[i]-=s->alpha*ap[i];}
}
__global__ void graphConvergence(IterationState* s,const double* scalars){
    if(!s->active)return;++s->iteration;double norm=scalars[1];
    if(!isfinite(norm)){s->error=3;s->active=0;return;}
    if(norm/s->bnorm<=s->tolerance){s->active=0;return;}
    if(s->iteration>=s->limit){s->error=4;s->active=0;}
}
__global__ void graphBeta(IterationState* s,const double* scalars){
    if(!s->active)return;double next=scalars[2];
    if(!isfinite(next)||!(next>0)){s->error=1;s->active=0;return;}
    s->beta=next/s->rho;s->rho=next;if(!isfinite(s->beta)){s->error=2;s->active=0;}
}
__global__ void graphDirection(int n,const IterationState* s,const double* z,double* p){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n&&s->active)p[i]=z[i]+s->beta*p[i];}

struct Context {
    int n=0,nz=0,device=-1,levels=0;unsigned solves=0;double shift=1.05;
    size_t symbolicDownload=0,symbolicUpload=0,setupPeakOwned=0;
    double symbolicSeconds=0,factorSeconds=0,analysisSeconds=0,graphSeconds=0;int iterationHostSyncs=0,pruned=0;
    Buffer<int> rp,ci,lrp,lci,urp,uci,free,failure;bool explicitFree=false;
    Buffer<double> av,l,u,x,r,p,ap,tmp,z,scalars;
    Buffer<IterationState> state;
    cudaStream_t stream=nullptr;cudaGraph_t graph=nullptr;cudaGraphExec_t graphExec=nullptr;
    // cudaMalloc alignment meets cuBLAS workspace alignment. Retain the buffer
    // until graph and cuBLAS handle destruction (before member destructors).
    Buffer<char> forwardBuffer,backwardBuffer,blasWorkspace;
    cusparseHandle_t sparse=nullptr;cublasHandle_t blas=nullptr;
    cusparseSpMatDescr_t lower=nullptr,upper=nullptr;cusparseDnVecDescr_t vr=nullptr,vt=nullptr,vz=nullptr;
    cusparseSpSVDescr_t forward=nullptr,backward=nullptr;
    ~Context(){if(device>=0)cudaSetDevice(device);if(graphExec)cudaGraphExecDestroy(graphExec);if(graph)cudaGraphDestroy(graph);if(forward)cusparseSpSV_destroyDescr(forward);if(backward)cusparseSpSV_destroyDescr(backward);if(vr)cusparseDestroyDnVec(vr);if(vt)cusparseDestroyDnVec(vt);if(vz)cusparseDestroyDnVec(vz);if(lower)cusparseDestroySpMat(lower);if(upper)cusparseDestroySpMat(upper);if(sparse)cusparseDestroy(sparse);if(blas)cublasDestroy(blas);if(stream)cudaStreamDestroy(stream);}
#ifdef RADIA_IC_LIFECYCLE_TEST
    bool simulateDeviceMismatch=false;
#endif
    void selected(){
#ifdef RADIA_IC_LIFECYCLE_TEST
        if(simulateDeviceMismatch)throw std::runtime_error("IC context belongs to another CUDA device");
#endif
        int d;ck(cudaGetDevice(&d));if(d!=device)throw std::runtime_error("IC context belongs to another CUDA device");}
    void flag(){int v;ck(cudaMemcpy(&v,failure.p,sizeof(v),cudaMemcpyDeviceToHost));if(v)throw std::runtime_error("Invalid/nonfinite/asymmetric matrix, RHS, or nonpositive IC pivot");}
    size_t bytes()const{return rp.bytes()+ci.bytes()+lrp.bytes()+lci.bytes()+urp.bytes()+uci.bytes()+free.bytes()+failure.bytes()+av.bytes()+l.bytes()+u.bytes()+x.bytes()+r.bytes()+p.bytes()+ap.bytes()+tmp.bytes()+z.bytes()+scalars.bytes()+state.bytes()+forwardBuffer.bytes()+backwardBuffer.bytes()+blasWorkspace.bytes();}
    void setup(const mxArray* op,double s){
        auto phase=Clock::now();
        shift=s;if(!(shift>=1&&shift<=100))throw std::runtime_error("IC fixed shift must be in [1,100]");ck(cudaGetDevice(&device));ck(cudaDeviceSynchronize());
        n=dimension(field(op,"rows"));if(dimension(field(op,"cols"))!=n)throw std::runtime_error("square SPD operator required");
        auto values=field(op,"values");if(!mxIsGPUArray(values))throw std::runtime_error("GPU CSR required");auto check=mxGPUCreateFromMxArray(values);size_t count=mxGPUGetNumberOfElements(check);mxGPUDestroyGPUArray(check);if(count>INT_MAX||count<size_t(n))throw std::runtime_error("int32 CSR storage limit");nz=int(count);
        Inputs input;auto rows=input.get<int>(field(op,"rowptr"),mxINT32_CLASS,size_t(n)+1);auto cols=input.get<int>(field(op,"columns"),mxINT32_CLASS,nz);auto vals=input.get<double>(values,mxDOUBLE_CLASS,nz);
        std::vector<int> hrp(size_t(n)+1),hci(nz);ck(cudaMemcpy(hrp.data(),rows,hrp.size()*sizeof(int),cudaMemcpyDeviceToHost));ck(cudaMemcpy(hci.data(),cols,hci.size()*sizeof(int),cudaMemcpyDeviceToHost));
        if(hrp[0]!=0||hrp[n]!=nz)throw std::runtime_error("CSR endpoints invalid");
        auto mask=mxGetField(op,0,"free");explicitFree=mask!=nullptr;std::vector<int> hostFree(n,1);
        if(mask){auto f=input.get<int>(mask,mxINT32_CLASS,n);ck(cudaMemcpy(hostFree.data(),f,size_t(n)*sizeof(int),cudaMemcpyDeviceToHost));for(int value:hostFree)if(value!=0&&value!=1)throw std::runtime_error("Invalid free mask");}
        std::vector<int> lp(size_t(n)+1),lc,source,level(n);lc.reserve(size_t(nz)/2+n);source.reserve(size_t(nz)/2+n);int maximum=0;
        for(int i=0;i<n;++i){if(hrp[i]<0||hrp[i+1]<hrp[i]||hrp[i+1]>nz)throw std::runtime_error("CSR row range invalid");int previous=-1;bool diagonal=false;
            for(int k=hrp[i];k<hrp[i+1];++k){int j=hci[k];if(j<0||j>=n||j<=previous)throw std::runtime_error("CSR columns must be sorted, unique and in range");previous=j;if(j<=i){if(j<i&&(!hostFree[i]||!hostFree[j])){++pruned;continue;}lc.push_back(j);source.push_back(k);if(j<i)level[i]=std::max(level[i],level[j]+1);else diagonal=true;}}
            if(!diagonal)throw std::runtime_error("Missing CSR diagonal");lp[i+1]=int(lc.size());maximum=std::max(maximum,level[i]);
        }
        symbolicDownload=(size_t(n)+1+size_t(nz)+(explicitFree?size_t(n):0))*sizeof(int);free.upload(hostFree);
        symbolicUpload=(3*(size_t(n)+1)+size_t(nz)+2*size_t(n)+4*lc.size()+2)*sizeof(int);
        levels=maximum+1;rp.upload(hrp);ci.upload(hci);lrp.upload(lp);lci.upload(lc);av.alloc(nz);ck(cudaMemcpy(av.p,vals,av.bytes(),cudaMemcpyDeviceToDevice));failure.upload(std::vector<int>{0,0});
        validate<<<blocks(n),T>>>(n,rp.p,ci.p,free.p,av.p,failure.p);flag();identityRows<<<blocks(n),T>>>(n,rp.p,ci.p,free.p,av.p);
        Buffer<int> mapping,ordered;mapping.upload(source);std::vector<int> offsets(size_t(levels)+1),order(n);
        for(int v:level)++offsets[v+1];for(int k=0;k<levels;++k)offsets[k+1]+=offsets[k];auto cursor=offsets;
        for(int i=0;i<n;++i)order[cursor[level[i]]++]=i;ordered.upload(order);
        std::vector<int> up(size_t(n)+1),uc(lc.size()),transposeMap(lc.size());
        for(int j:lc)++up[j+1];for(int i=0;i<n;++i)up[i+1]+=up[i];auto next=up;
        for(int i=0;i<n;++i)for(int k=lp[i];k<lp[i+1];++k){int dest=next[lc[k]]++;uc[dest]=i;transposeMap[dest]=k;}
        Buffer<int> transpose;transpose.upload(transposeMap);urp.upload(up);uci.upload(uc);u.alloc(lc.size());
        ck(cudaDeviceSynchronize());symbolicSeconds=elapsed(phase);phase=Clock::now();
        Buffer<double> ld;ld.alloc(lc.size());l.alloc(lc.size());lowerValues<<<blocks(lc.size()),T>>>(int(lc.size()),mapping.p,av.p,ld.p);
        for(int k=0;k<levels;++k){int size=offsets[k+1]-offsets[k];factorLevel<<<blocks(size),T>>>(size,ordered.p+offsets[k],lrp.p,lci.p,ld.p,shift,failure.p);}
        flag();toCholesky<<<blocks(n),T>>>(n,lrp.p,lci.p,ld.p,l.p);ck(cudaGetLastError());
        // Only already-validated constrained zero couplings were omitted. Free
        // IC(0) entries/order and the original CSR operator remain unchanged.
        lowerValues<<<blocks(lc.size()),T>>>(int(lc.size()),transpose.p,l.p,u.p);
        ck(cudaDeviceSynchronize());factorSeconds=elapsed(phase);phase=Clock::now();
        x.alloc(n);r.alloc(n);p.alloc(n);ap.alloc(n);tmp.alloc(n);z.alloc(n);
        cs(cusparseCreate(&sparse));cb(cublasCreate(&blas));ck(cudaStreamCreate(&stream));cs(cusparseSetStream(sparse,stream));cb(cublasSetStream(blas,stream));cs(cusparseCreateCsr(&lower,n,n,lc.size(),lrp.p,lci.p,l.p,CUSPARSE_INDEX_32I,CUSPARSE_INDEX_32I,CUSPARSE_INDEX_BASE_ZERO,CUDA_R_64F));
        cs(cusparseCreateCsr(&upper,n,n,lc.size(),urp.p,uci.p,u.p,CUSPARSE_INDEX_32I,CUSPARSE_INDEX_32I,CUSPARSE_INDEX_BASE_ZERO,CUDA_R_64F));
        cusparseFillMode_t fill=CUSPARSE_FILL_MODE_LOWER;cusparseDiagType_t diag=CUSPARSE_DIAG_TYPE_NON_UNIT;
        cs(cusparseSpMatSetAttribute(lower,CUSPARSE_SPMAT_FILL_MODE,&fill,sizeof(fill)));cs(cusparseSpMatSetAttribute(lower,CUSPARSE_SPMAT_DIAG_TYPE,&diag,sizeof(diag)));
        fill=CUSPARSE_FILL_MODE_UPPER;cs(cusparseSpMatSetAttribute(upper,CUSPARSE_SPMAT_FILL_MODE,&fill,sizeof(fill)));cs(cusparseSpMatSetAttribute(upper,CUSPARSE_SPMAT_DIAG_TYPE,&diag,sizeof(diag)));
        cs(cusparseCreateDnVec(&vr,n,r.p,CUDA_R_64F));cs(cusparseCreateDnVec(&vt,n,tmp.p,CUDA_R_64F));cs(cusparseCreateDnVec(&vz,n,z.p,CUDA_R_64F));cs(cusparseSpSV_createDescr(&forward));cs(cusparseSpSV_createDescr(&backward));
        double one=1;size_t a=0,b=0;
        cs(cusparseSpSV_bufferSize(sparse,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,lower,vr,vt,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,forward,&a));
        cs(cusparseSpSV_bufferSize(sparse,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,upper,vt,vz,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,backward,&b));forwardBuffer.alloc(a);backwardBuffer.alloc(b);
        cs(cusparseSpSV_analysis(sparse,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,lower,vr,vt,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,forward,forwardBuffer.p));
        cs(cusparseSpSV_analysis(sparse,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,upper,vt,vz,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,backward,backwardBuffer.p));ck(cudaDeviceSynchronize());
        analysisSeconds=elapsed(phase);phase=Clock::now();scalars.alloc(3);state.alloc(1);
        // SetStream resets a handle's workspace, so assign this only after the
        // stream is fixed and before warmup/capture. No later SetStream calls.
        blasWorkspace.alloc(4ULL<<20);cb(cublasSetWorkspace(blas,blasWorkspace.p,blasWorkspace.bytes()));
        captureGraph();graphSeconds=elapsed(phase);setupPeakOwned=bytes()+mapping.bytes()+ordered.bytes()+ld.bytes()+transpose.bytes();
    }
    void apply(){double one=1;cs(cusparseSpSV_solve(sparse,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,lower,vr,vt,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,forward));cs(cusparseSpSV_solve(sparse,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,upper,vt,vz,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,backward));}
    void captureGraph(){
        cb(cublasSetPointerMode(blas,CUBLAS_POINTER_MODE_DEVICE));ck(cudaMemsetAsync(r.p,0,r.bytes(),stream));ck(cudaMemsetAsync(p.p,0,p.bytes(),stream));
        cb(cublasDdot(blas,n,r.p,1,r.p,1,scalars.p));cb(cublasDnrm2(blas,n,r.p,1,scalars.p+1));ck(cudaStreamSynchronize(stream));
        ck(cudaStreamBeginCapture(stream,cudaStreamCaptureModeThreadLocal));
        try{
            for(int k=0;k<GRAPH_BATCH;++k){
                spmv<<<blocks(n),T,0,stream>>>(n,rp.p,ci.p,av.p,p.p,ap.p);cb(cublasDdot(blas,n,p.p,1,ap.p,1,scalars.p));
                graphAlpha<<<1,1,0,stream>>>(state.p,scalars.p);graphUpdate<<<blocks(n),T,0,stream>>>(n,state.p,p.p,ap.p,x.p,r.p);
                cb(cublasDnrm2(blas,n,r.p,1,scalars.p+1));graphConvergence<<<1,1,0,stream>>>(state.p,scalars.p);
                apply();cb(cublasDdot(blas,n,r.p,1,z.p,1,scalars.p+2));graphBeta<<<1,1,0,stream>>>(state.p,scalars.p);graphDirection<<<blocks(n),T,0,stream>>>(n,state.p,z.p,p.p);
            }
            ck(cudaStreamEndCapture(stream,&graph));ck(cudaGraphInstantiate(&graphExec,graph,nullptr,nullptr,0));
        }catch(...){cudaGraph_t aborted=nullptr;cudaStreamEndCapture(stream,&aborted);if(aborted)cudaGraphDestroy(aborted);throw;}
        cb(cublasSetPointerMode(blas,CUBLAS_POINTER_MODE_HOST));
    }
    double norm(const double* v){double value;cb(cublasDnrm2(blas,n,v,1,&value));if(!std::isfinite(value))throw std::runtime_error("nonfinite FP64 norm");return value;}
    double dot(const double* a,const double* b){double value;cb(cublasDdot(blas,n,a,1,b,1,&value));if(!std::isfinite(value))throw std::runtime_error("nonfinite PCG inner product");return value;}
    mxArray* profile(const double* rhs,int count){
        ck(cudaMemset(failure.p,0,2*sizeof(int)));checkVector<<<blocks(n),T>>>(n,rhs,free.p,failure.p);flag();
        ck(cudaMemcpy(r.p,rhs,r.bytes(),cudaMemcpyDeviceToDevice));apply();
        auto measure=[&](auto operation){ck(cudaDeviceSynchronize());auto start=Clock::now();for(int i=0;i<count;++i)operation();ck(cudaDeviceSynchronize());return elapsed(start)/count;};
        double one=1;
        const char* names[]={"forward_seconds","backward_seconds","spmv_seconds","dot_host_seconds","norm_host_seconds","repetitions"};
        double values[]={measure([&]{cs(cusparseSpSV_solve(sparse,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,lower,vr,vt,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,forward));}),
            measure([&]{cs(cusparseSpSV_solve(sparse,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,upper,vt,vz,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,backward));}),
            measure([&]{spmv<<<blocks(n),T>>>(n,rp.p,ci.p,av.p,r.p,ap.p);}),measure([&]{dot(r.p,r.p);}),measure([&]{norm(r.p);}),double(count)};
        auto result=mxCreateStructMatrix(1,1,6,names);for(int i=0;i<6;++i)mxSetField(result,0,names[i],mxCreateDoubleScalar(values[i]));return result;
    }
    std::pair<int,double> solve(const double* rhs,double tol,int limit){
        selected();iterationHostSyncs=0;ck(cudaMemset(failure.p,0,2*sizeof(int)));checkVector<<<blocks(n),T>>>(n,rhs,free.p,failure.p);flag();ck(cudaMemset(x.p,0,x.bytes()));ck(cudaMemcpy(r.p,rhs,r.bytes(),cudaMemcpyDeviceToDevice));double bnorm=norm(rhs);int iterations=0;
        if(bnorm==0){int nonzero;ck(cudaMemcpy(&nonzero,failure.p+1,sizeof(int),cudaMemcpyDeviceToHost));if(nonzero)throw std::runtime_error("RHS norm outside supported FP64 range; explicitly rescale the system");}
        bool converged=bnorm==0;
        if(bnorm>0){apply();ck(cudaMemcpy(p.p,z.p,p.bytes(),cudaMemcpyDeviceToDevice));double rho=dot(r.p,z.p);
            if(!(rho>0))throw std::runtime_error("IC PCG nonpositive preconditioned curvature");
            IterationState status{rho,0,0,bnorm,tol,0,limit,1,0};ck(cudaMemcpyAsync(state.p,&status,sizeof(status),cudaMemcpyHostToDevice,stream));
            for(int k=0;k<limit;k+=GRAPH_BATCH){ck(cudaGraphLaunch(graphExec,stream));ck(cudaMemcpyAsync(&status,state.p,sizeof(status),cudaMemcpyDeviceToHost,stream));ck(cudaStreamSynchronize(stream));++iterationHostSyncs;if(!status.active)break;}
            if(status.error==1)throw std::runtime_error("IC PCG nonpositive/nonfinite curvature");
            if(status.error==2)throw std::runtime_error("IC PCG invalid step or direction");
            if(status.error==3)throw std::runtime_error("IC PCG nonfinite residual norm");
            iterations=status.iteration;converged=!status.active&&status.error==0;
        }
        if(!converged)throw std::runtime_error("IC PCG iteration limit exhausted");
        spmv<<<blocks(n),T>>>(n,rp.p,ci.p,av.p,x.p,ap.p);subtract<<<blocks(n),T>>>(n,rhs,ap.p,bnorm);double relative=norm(ap.p);
        if(!std::isfinite(relative)||relative>std::min(1e-6,std::max(10*tol,1e-12)))throw std::runtime_error("IC PCG projected original operator true residual failed");++solves;return {iterations,relative};
    }
};
using Handle=uint64_t;constexpr Handle TAG=0x1CC0;
std::unordered_map<Handle,std::unique_ptr<Context>> contexts;
Handle id(const mxArray* a){if(mxGetClassID(a)!=mxUINT64_CLASS||mxIsComplex(a)||mxGetNumberOfElements(a)!=1)throw std::runtime_error("typed uint64 IC handle required");Handle h=*static_cast<const Handle*>(mxGetData(a));if((h>>48)!=TAG)throw std::runtime_error("IC handle type mismatch");return h;}
auto lookup(const mxArray* a){auto it=contexts.find(id(a));if(it==contexts.end())throw std::runtime_error("stale/released IC handle");return it;}
Context& context(const mxArray* a){auto it=lookup(a);it->second->selected();return *it->second;}
void release(const mxArray* a){
    auto it=lookup(a); // Validate type/liveness without requiring the owning device.
    int caller;ck(cudaGetDevice(&caller));
    contexts.erase(it); // Context and buffer destructors free on their owning device.
    mexUnlock(); // Balance the erased context even if restoring the device fails.
    ck(cudaSetDevice(caller));
}
void atExit(){int d=-1;cudaGetDevice(&d);contexts.clear();if(d>=0)cudaSetDevice(d);}
mxArray* copyGPU(const double* p,int n){mwSize size=n;auto a=mxGPUCreateGPUArray(1,&size,mxDOUBLE_CLASS,mxREAL,MX_GPU_DO_NOT_INITIALIZE);try{ck(cudaMemcpy(mxGPUGetData(a),p,size_t(n)*sizeof(double),cudaMemcpyDeviceToDevice));auto out=mxGPUCreateMxArrayOnGPU(a);mxGPUDestroyGPUArray(a);return out;}catch(...){mxGPUDestroyGPUArray(a);throw;}}
mxArray* info(Context& c,int iterations,double residual,bool checked=true){
    const char* names[]={"iterations","true_relative_residual","shift","factor_levels","gpu_owned_bytes","solve_count","backend","preconditioner","factorization_backend","triangular_backend","residual_checked","symbolic_download_bytes","symbolic_upload_bytes","gpu_setup_peak_owned_bytes","memory_definition","explicit_free_mask","symbolic_and_input_seconds","numeric_factor_seconds","triangular_analysis_seconds","graph_capture_seconds","iteration_host_syncs","graph_batch_iterations","graph_cached","factor_nonzeros","projected_zero_entries_pruned"};
    auto a=mxCreateStructMatrix(1,1,25,names);double v[]={double(iterations),residual,c.shift,double(c.levels),double(c.bytes()),double(c.solves)};
    for(int i=0;i<6;++i)mxSetField(a,0,names[i],i==1&&!checked?mxCreateDoubleMatrix(0,0,mxREAL):mxCreateDoubleScalar(v[i]));
    const char* text[]={"MATLAB/mexcuda CUDA FP64","natural-order fixed-shift IC(0)","CUDA level-scheduled incomplete LDLT converted to Cholesky","cuSPARSE SpSV explicit lower/upper non-transpose CSR"};for(int i=0;i<4;++i)mxSetField(a,0,names[i+6],mxCreateString(text[i]));
    mxSetField(a,0,names[10],mxCreateLogicalScalar(checked));mxSetField(a,0,names[11],mxCreateDoubleScalar(double(c.symbolicDownload)));mxSetField(a,0,names[12],mxCreateDoubleScalar(double(c.symbolicUpload)));mxSetField(a,0,names[13],mxCreateDoubleScalar(double(c.setupPeakOwned)));mxSetField(a,0,names[14],mxCreateString("Explicit owned CUDA buffers and queried SpSV workspaces; excludes opaque library/runtime allocations and caller inputs/outputs"));mxSetField(a,0,names[15],mxCreateLogicalScalar(c.explicitFree));
    mxSetField(a,0,names[16],mxCreateDoubleScalar(c.symbolicSeconds));mxSetField(a,0,names[17],mxCreateDoubleScalar(c.factorSeconds));mxSetField(a,0,names[18],mxCreateDoubleScalar(c.analysisSeconds));
    mxSetField(a,0,names[19],mxCreateDoubleScalar(c.graphSeconds));mxSetField(a,0,names[20],mxCreateDoubleScalar(c.iterationHostSyncs));mxSetField(a,0,names[21],mxCreateDoubleScalar(GRAPH_BATCH));mxSetField(a,0,names[22],mxCreateLogicalScalar(c.graphExec!=nullptr));mxSetField(a,0,names[23],mxCreateDoubleScalar(c.l.n));mxSetField(a,0,names[24],mxCreateDoubleScalar(c.pruned));return a;
}
}
void dispatchIc(int nlhs,mxArray** out,int nrhs,const mxArray** args){mxInitGPU();if(nrhs<1||!mxIsChar(args[0]))throw std::runtime_error("IC command required");char* chars=mxArrayToString(args[0]);if(!chars)throw std::runtime_error("invalid command");std::string command(chars);mxFree(chars);
        if(command=="create"){
            if(nrhs!=3||nlhs!=2)throw std::runtime_error("[handle,info]=radia_cuda_ic('create',GPU_CSR,shift)");auto c=std::make_unique<Context>();c->setup(args[1],scalar(args[2]));
            static Handle generation=Handle(std::chrono::high_resolution_clock::now().time_since_epoch().count())&0x0000ffffffffffffULL;if(generation==0x0000ffffffffffffULL)throw std::runtime_error("IC generation exhausted");Handle h=(TAG<<48)|++generation;
            out[1]=info(*c,0,0,false);out[0]=mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);*static_cast<Handle*>(mxGetData(out[0]))=h;contexts.emplace(h,std::move(c));static bool registered=false;if(!registered){mexAtExit(atExit);registered=true;}mexLock();return;
        }
        if(command=="release"){if(nrhs!=2||nlhs!=0)throw std::runtime_error("release requires IC handle");release(args[1]);return;}
#ifdef RADIA_IC_LIFECYCLE_TEST
        if(command=="test_device_mismatch"){
            if(nrhs!=2||nlhs!=0)throw std::runtime_error("test_device_mismatch arity");
            lookup(args[1])->second->simulateDeviceMismatch=true;return;
        }
#endif
        if(command=="live_count"){if(nrhs!=1||nlhs!=1)throw std::runtime_error("live_count arity");out[0]=mxCreateDoubleScalar(double(contexts.size()));return;}
        if(command=="profile"){if(nrhs!=4||nlhs!=1)throw std::runtime_error("profile(handle,gpuRhs,repetitions)");auto& c=context(args[1]);Inputs input;auto b=input.get<double>(args[2],mxDOUBLE_CLASS,c.n);int count=dimension(args[3]);if(count>500)throw std::runtime_error("Profile repetitions exceed bound");out[0]=c.profile(b,count);return;}
        if(command=="solve"||command=="apply"){
            bool solving=command=="solve";if(nrhs!=(solving?5:3)||nlhs!=(solving?2:1))throw std::runtime_error("IC solve/apply arity");auto& c=context(args[1]);Inputs input;auto b=input.get<double>(args[2],mxDOUBLE_CLASS,c.n);ck(cudaDeviceSynchronize());
            if(solving){double tol=scalar(args[3]),lim=scalar(args[4]);if(!(tol>0&&tol<=1e-6)||lim<1||lim>100000||lim!=floor(lim))throw std::runtime_error("IC tolerance or iteration limit invalid");auto result=c.solve(b,tol,int(lim));out[0]=copyGPU(c.x.p,c.n);out[1]=info(c,result.first,result.second);}
            else{ck(cudaMemset(c.failure.p,0,2*sizeof(int)));checkVector<<<blocks(c.n),T>>>(c.n,b,c.free.p,c.failure.p);c.flag();ck(cudaMemcpy(c.r.p,b,c.r.bytes(),cudaMemcpyDeviceToDevice));c.apply();checkVector<<<blocks(c.n),T>>>(c.n,c.z.p,c.free.p,c.failure.p);c.flag();out[0]=copyGPU(c.z.p,c.n);}return;
        }
        throw std::runtime_error("unknown IC command");

}
void mexFunction(int nl,mxArray** out,int nr,const mxArray** in){

    char errorMessage[2048]={0};
    if(!radia_mex::checked(dispatchIc,nl,out,nr,in,errorMessage,sizeof(errorMessage))){
        mexErrMsgIdAndTxt("radia:cuda:ic","%s",errorMessage);
    }
}
