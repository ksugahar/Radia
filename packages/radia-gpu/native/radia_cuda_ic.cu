/* SPDX-License-Identifier: MPL-2.0 */
// FP64 natural-order IC(0). Host performs symbolic indexing only; CUDA owns
// numerical factorization, triangular solves, PCG and original-A residuals.
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
__global__ void checkVector(int n,const double* x,const int* free,int* failure){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n&&(!isfinite(x[i])||(!free[i]&&x[i]!=0)))atomicExch(failure,1);}
__global__ void spmv(int n,const int* rp,const int* ci,const double* av,const double* x,double* y){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n){double s=0;for(int p=rp[i];p<rp[i+1];++p)s+=av[p]*x[ci[p]];y[i]=s;}}
__global__ void update(int n,double alpha,const double* p,const double* ap,double* x,double* r){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n){x[i]+=alpha*p[i];r[i]-=alpha*ap[i];}}
__global__ void direction(int n,double beta,const double* z,double* p){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)p[i]=z[i]+beta*p[i];}
__global__ void subtract(int n,const double* b,double* ax){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)ax[i]=b[i]-ax[i];}

struct Context {
    int n=0,nz=0,device=-1,levels=0;unsigned solves=0;double shift=1.05;
    size_t symbolicDownload=0,symbolicUpload=0,setupPeakOwned=0;
    Buffer<int> rp,ci,lrp,lci,free,failure;bool explicitFree=false;
    Buffer<double> av,l,x,r,p,ap,tmp,z;
    Buffer<char> forwardBuffer,backwardBuffer;
    cusparseHandle_t sparse=nullptr;cublasHandle_t blas=nullptr;
    cusparseSpMatDescr_t lower=nullptr;cusparseDnVecDescr_t vr=nullptr,vt=nullptr,vz=nullptr;
    cusparseSpSVDescr_t forward=nullptr,backward=nullptr;
    ~Context(){if(device>=0)cudaSetDevice(device);if(forward)cusparseSpSV_destroyDescr(forward);if(backward)cusparseSpSV_destroyDescr(backward);if(vr)cusparseDestroyDnVec(vr);if(vt)cusparseDestroyDnVec(vt);if(vz)cusparseDestroyDnVec(vz);if(lower)cusparseDestroySpMat(lower);if(sparse)cusparseDestroy(sparse);if(blas)cublasDestroy(blas);}
    void selected(){int d;ck(cudaGetDevice(&d));if(d!=device)throw std::runtime_error("IC context belongs to another CUDA device");}
    void flag(){int v;ck(cudaMemcpy(&v,failure.p,sizeof(v),cudaMemcpyDeviceToHost));if(v)throw std::runtime_error("Invalid/nonfinite/asymmetric matrix, RHS, or nonpositive IC pivot");}
    size_t bytes()const{return rp.bytes()+ci.bytes()+lrp.bytes()+lci.bytes()+free.bytes()+failure.bytes()+av.bytes()+l.bytes()+x.bytes()+r.bytes()+p.bytes()+ap.bytes()+tmp.bytes()+z.bytes()+forwardBuffer.bytes()+backwardBuffer.bytes();}
    void setup(const mxArray* op,double s){
        shift=s;if(!(shift>=1&&shift<=100))throw std::runtime_error("IC fixed shift must be in [1,100]");ck(cudaGetDevice(&device));ck(cudaDeviceSynchronize());
        n=dimension(field(op,"rows"));if(dimension(field(op,"cols"))!=n)throw std::runtime_error("square SPD operator required");
        auto values=field(op,"values");if(!mxIsGPUArray(values))throw std::runtime_error("GPU CSR required");auto check=mxGPUCreateFromMxArray(values);size_t count=mxGPUGetNumberOfElements(check);mxGPUDestroyGPUArray(check);if(count>INT_MAX||count<size_t(n))throw std::runtime_error("int32 CSR storage limit");nz=int(count);
        Inputs input;auto rows=input.get<int>(field(op,"rowptr"),mxINT32_CLASS,size_t(n)+1);auto cols=input.get<int>(field(op,"columns"),mxINT32_CLASS,nz);auto vals=input.get<double>(values,mxDOUBLE_CLASS,nz);
        std::vector<int> hrp(size_t(n)+1),hci(nz);ck(cudaMemcpy(hrp.data(),rows,hrp.size()*sizeof(int),cudaMemcpyDeviceToHost));ck(cudaMemcpy(hci.data(),cols,hci.size()*sizeof(int),cudaMemcpyDeviceToHost));
        if(hrp[0]!=0||hrp[n]!=nz)throw std::runtime_error("CSR endpoints invalid");
        std::vector<int> lp(size_t(n)+1),lc,source,level(n);int maximum=0;
        for(int i=0;i<n;++i){if(hrp[i]<0||hrp[i+1]<hrp[i]||hrp[i+1]>nz)throw std::runtime_error("CSR row range invalid");int previous=-1;bool diagonal=false;
            for(int k=hrp[i];k<hrp[i+1];++k){int j=hci[k];if(j<0||j>=n||j<=previous)throw std::runtime_error("CSR columns must be sorted, unique and in range");previous=j;if(j<=i){lc.push_back(j);source.push_back(k);if(j<i)level[i]=std::max(level[i],level[j]+1);else diagonal=true;}}
            if(!diagonal)throw std::runtime_error("Missing CSR diagonal");lp[i+1]=int(lc.size());maximum=std::max(maximum,level[i]);
        }
        symbolicDownload=(size_t(n)+1+size_t(nz))*sizeof(int);
        auto mask=mxGetField(op,0,"free");explicitFree=mask!=nullptr;
        if(mask){auto f=input.get<int>(mask,mxINT32_CLASS,n);free.alloc(n);ck(cudaMemcpy(free.p,f,free.bytes(),cudaMemcpyDeviceToDevice));}
        else free.upload(std::vector<int>(n,1));
        symbolicUpload=(size_t(n)+1+size_t(nz)+size_t(n)+1+2*lc.size()+size_t(n)+1+(explicitFree?0:size_t(n)))*sizeof(int);
        levels=maximum+1;rp.upload(hrp);ci.upload(hci);lrp.upload(lp);lci.upload(lc);av.alloc(nz);ck(cudaMemcpy(av.p,vals,av.bytes(),cudaMemcpyDeviceToDevice));failure.upload(std::vector<int>{0});
        validate<<<blocks(n),T>>>(n,rp.p,ci.p,free.p,av.p,failure.p);flag();identityRows<<<blocks(n),T>>>(n,rp.p,ci.p,free.p,av.p);
        Buffer<int> mapping,ordered;mapping.upload(source);std::vector<int> offsets(size_t(levels)+1),order(n);
        for(int v:level)++offsets[v+1];for(int k=0;k<levels;++k)offsets[k+1]+=offsets[k];auto cursor=offsets;
        for(int i=0;i<n;++i)order[cursor[level[i]]++]=i;ordered.upload(order);
        Buffer<double> ld;ld.alloc(lc.size());l.alloc(lc.size());lowerValues<<<blocks(lc.size()),T>>>(int(lc.size()),mapping.p,av.p,ld.p);
        for(int k=0;k<levels;++k){int size=offsets[k+1]-offsets[k];factorLevel<<<blocks(size),T>>>(size,ordered.p+offsets[k],lrp.p,lci.p,ld.p,shift,failure.p);}
        flag();toCholesky<<<blocks(n),T>>>(n,lrp.p,lci.p,ld.p,l.p);ck(cudaGetLastError());
        x.alloc(n);r.alloc(n);p.alloc(n);ap.alloc(n);tmp.alloc(n);z.alloc(n);
        cs(cusparseCreate(&sparse));cb(cublasCreate(&blas));cs(cusparseCreateCsr(&lower,n,n,lc.size(),lrp.p,lci.p,l.p,CUSPARSE_INDEX_32I,CUSPARSE_INDEX_32I,CUSPARSE_INDEX_BASE_ZERO,CUDA_R_64F));
        cusparseFillMode_t fill=CUSPARSE_FILL_MODE_LOWER;cusparseDiagType_t diag=CUSPARSE_DIAG_TYPE_NON_UNIT;
        cs(cusparseSpMatSetAttribute(lower,CUSPARSE_SPMAT_FILL_MODE,&fill,sizeof(fill)));cs(cusparseSpMatSetAttribute(lower,CUSPARSE_SPMAT_DIAG_TYPE,&diag,sizeof(diag)));
        cs(cusparseCreateDnVec(&vr,n,r.p,CUDA_R_64F));cs(cusparseCreateDnVec(&vt,n,tmp.p,CUDA_R_64F));cs(cusparseCreateDnVec(&vz,n,z.p,CUDA_R_64F));cs(cusparseSpSV_createDescr(&forward));cs(cusparseSpSV_createDescr(&backward));
        double one=1;size_t a=0,b=0;
        cs(cusparseSpSV_bufferSize(sparse,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,lower,vr,vt,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,forward,&a));
        cs(cusparseSpSV_bufferSize(sparse,CUSPARSE_OPERATION_TRANSPOSE,&one,lower,vt,vz,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,backward,&b));forwardBuffer.alloc(a);backwardBuffer.alloc(b);
        cs(cusparseSpSV_analysis(sparse,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,lower,vr,vt,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,forward,forwardBuffer.p));
        cs(cusparseSpSV_analysis(sparse,CUSPARSE_OPERATION_TRANSPOSE,&one,lower,vt,vz,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,backward,backwardBuffer.p));ck(cudaDeviceSynchronize());
        setupPeakOwned=bytes()+mapping.bytes()+ordered.bytes()+ld.bytes();
    }
    void apply(){double one=1;cs(cusparseSpSV_solve(sparse,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,lower,vr,vt,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,forward));cs(cusparseSpSV_solve(sparse,CUSPARSE_OPERATION_TRANSPOSE,&one,lower,vt,vz,CUDA_R_64F,CUSPARSE_SPSV_ALG_DEFAULT,backward));}
    double norm(const double* v){double value;cb(cublasDnrm2(blas,n,v,1,&value));if(!std::isfinite(value))throw std::runtime_error("nonfinite FP64 norm");return value;}
    double dot(const double* a,const double* b){double value;cb(cublasDdot(blas,n,a,1,b,1,&value));if(!std::isfinite(value))throw std::runtime_error("nonfinite PCG inner product");return value;}
    std::pair<int,double> solve(const double* rhs,double tol,int limit){
        selected();ck(cudaMemset(failure.p,0,sizeof(int)));checkVector<<<blocks(n),T>>>(n,rhs,free.p,failure.p);flag();ck(cudaMemset(x.p,0,x.bytes()));ck(cudaMemcpy(r.p,rhs,r.bytes(),cudaMemcpyDeviceToDevice));double bnorm=norm(rhs);int iterations=0;
        bool converged=bnorm==0;
        if(bnorm>0){apply();ck(cudaMemcpy(p.p,z.p,p.bytes(),cudaMemcpyDeviceToDevice));double rho=dot(r.p,z.p);
            for(int k=0;k<limit;++k){if(!(rho>0))throw std::runtime_error("IC PCG nonpositive preconditioned curvature");spmv<<<blocks(n),T>>>(n,rp.p,ci.p,av.p,p.p,ap.p);double pap=dot(p.p,ap.p);if(!(pap>0))throw std::runtime_error("IC PCG nonpositive operator curvature");double alpha=rho/pap;if(!std::isfinite(alpha))throw std::runtime_error("IC PCG invalid step");update<<<blocks(n),T>>>(n,alpha,p.p,ap.p,x.p,r.p);iterations=k+1;
                if(norm(r.p)/bnorm<=tol){converged=true;break;}
                apply();double next=dot(r.p,z.p);direction<<<blocks(n),T>>>(n,next/rho,z.p,p.p);rho=next;
            }
        }
        if(!converged)throw std::runtime_error("IC PCG iteration limit exhausted");
        spmv<<<blocks(n),T>>>(n,rp.p,ci.p,av.p,x.p,ap.p);subtract<<<blocks(n),T>>>(n,rhs,ap.p);double error=norm(ap.p);double relative=bnorm>0?error/bnorm:error;
        if(!std::isfinite(relative)||relative>std::min(1e-6,std::max(10*tol,1e-12)))throw std::runtime_error("IC PCG original-A true residual failed");++solves;return {iterations,relative};
    }
};
using Handle=uint64_t;constexpr Handle TAG=0x1CC0;
std::unordered_map<Handle,std::unique_ptr<Context>> contexts;
Handle id(const mxArray* a){if(mxGetClassID(a)!=mxUINT64_CLASS||mxIsComplex(a)||mxGetNumberOfElements(a)!=1)throw std::runtime_error("typed uint64 IC handle required");Handle h=*static_cast<const Handle*>(mxGetData(a));if((h>>48)!=TAG)throw std::runtime_error("IC handle type mismatch");return h;}
Context& context(const mxArray* a){auto it=contexts.find(id(a));if(it==contexts.end())throw std::runtime_error("stale/released IC handle");it->second->selected();return *it->second;}
void atExit(){int d=-1;cudaGetDevice(&d);contexts.clear();if(d>=0)cudaSetDevice(d);}
mxArray* copyGPU(const double* p,int n){mwSize size=n;auto a=mxGPUCreateGPUArray(1,&size,mxDOUBLE_CLASS,mxREAL,MX_GPU_DO_NOT_INITIALIZE);try{ck(cudaMemcpy(mxGPUGetData(a),p,size_t(n)*sizeof(double),cudaMemcpyDeviceToDevice));auto out=mxGPUCreateMxArrayOnGPU(a);mxGPUDestroyGPUArray(a);return out;}catch(...){mxGPUDestroyGPUArray(a);throw;}}
mxArray* info(Context& c,int iterations,double residual,bool checked=true){
    const char* names[]={"iterations","true_relative_residual","shift","factor_levels","gpu_owned_bytes","solve_count","backend","preconditioner","factorization_backend","triangular_backend","residual_checked","symbolic_download_bytes","symbolic_upload_bytes","gpu_setup_peak_owned_bytes","memory_definition","explicit_free_mask"};
    auto a=mxCreateStructMatrix(1,1,16,names);double v[]={double(iterations),residual,c.shift,double(c.levels),double(c.bytes()),double(c.solves)};
    for(int i=0;i<6;++i)mxSetField(a,0,names[i],i==1&&!checked?mxCreateDoubleMatrix(0,0,mxREAL):mxCreateDoubleScalar(v[i]));
    const char* text[]={"MATLAB/mexcuda CUDA FP64","natural-order fixed-shift IC(0)","CUDA level-scheduled incomplete LDLT converted to Cholesky","cuSPARSE SpSV forward and transpose"};for(int i=0;i<4;++i)mxSetField(a,0,names[i+6],mxCreateString(text[i]));
    mxSetField(a,0,names[10],mxCreateLogicalScalar(checked));mxSetField(a,0,names[11],mxCreateDoubleScalar(double(c.symbolicDownload)));mxSetField(a,0,names[12],mxCreateDoubleScalar(double(c.symbolicUpload)));mxSetField(a,0,names[13],mxCreateDoubleScalar(double(c.setupPeakOwned)));mxSetField(a,0,names[14],mxCreateString("Explicit owned CUDA buffers and queried SpSV workspaces; excludes opaque library/runtime allocations and caller inputs/outputs"));mxSetField(a,0,names[15],mxCreateLogicalScalar(c.explicitFree));return a;
}
}
void mexFunction(int nlhs,mxArray** out,int nrhs,const mxArray** args){
    try{mxInitGPU();if(nrhs<1||!mxIsChar(args[0]))throw std::runtime_error("IC command required");char* chars=mxArrayToString(args[0]);if(!chars)throw std::runtime_error("invalid command");std::string command(chars);mxFree(chars);
        if(command=="create"){
            if(nrhs!=3||nlhs!=2)throw std::runtime_error("[handle,info]=radia_cuda_ic('create',GPU_CSR,shift)");auto c=std::make_unique<Context>();c->setup(args[1],scalar(args[2]));
            static Handle generation=Handle(std::chrono::high_resolution_clock::now().time_since_epoch().count())&0x0000ffffffffffffULL;if(generation==0x0000ffffffffffffULL)throw std::runtime_error("IC generation exhausted");Handle h=(TAG<<48)|++generation;
            out[1]=info(*c,0,0,false);out[0]=mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);*static_cast<Handle*>(mxGetData(out[0]))=h;contexts.emplace(h,std::move(c));static bool registered=false;if(!registered){mexAtExit(atExit);registered=true;}mexLock();return;
        }
        if(command=="release"){if(nrhs!=2||nlhs!=0)throw std::runtime_error("release requires IC handle");context(args[1]);contexts.erase(id(args[1]));mexUnlock();return;}
        if(command=="live_count"){if(nrhs!=1||nlhs!=1)throw std::runtime_error("live_count arity");out[0]=mxCreateDoubleScalar(double(contexts.size()));return;}
        if(command=="solve"||command=="apply"){
            bool solving=command=="solve";if(nrhs!=(solving?5:3)||nlhs!=(solving?2:1))throw std::runtime_error("IC solve/apply arity");auto& c=context(args[1]);Inputs input;auto b=input.get<double>(args[2],mxDOUBLE_CLASS,c.n);ck(cudaDeviceSynchronize());
            if(solving){double tol=scalar(args[3]),lim=scalar(args[4]);if(!(tol>0&&tol<=1e-6)||lim<1||lim>100000||lim!=floor(lim))throw std::runtime_error("IC tolerance or iteration limit invalid");auto result=c.solve(b,tol,int(lim));out[0]=copyGPU(c.x.p,c.n);out[1]=info(c,result.first,result.second);}
            else{ck(cudaMemset(c.failure.p,0,sizeof(int)));checkVector<<<blocks(c.n),T>>>(c.n,b,c.free.p,c.failure.p);c.flag();ck(cudaMemcpy(c.r.p,b,c.r.bytes(),cudaMemcpyDeviceToDevice));c.apply();checkVector<<<blocks(c.n),T>>>(c.n,c.z.p,c.free.p,c.failure.p);c.flag();out[0]=copyGPU(c.z.p,c.n);}return;
        }
        throw std::runtime_error("unknown IC command");
    }catch(const std::exception& ex){mexErrMsgIdAndTxt("radia:cuda:ic", "%s",ex.what());}
}
