#include "mex_error_boundary.hpp"
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
#include <unordered_map>
#include <chrono>
#include <cstdint>
#include "ams_fem_host.hpp"
void RadiaAMSSetup(int,mxArray**,int,const mxArray**);
namespace {
constexpr int T=256;
// Nonzero only while recording a preconditioner graph; normal work uses MATLAB's legacy stream.
thread_local cudaStream_t launch_stream=nullptr;
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
struct Flag {
    int* p=nullptr; bool owned=true;
    explicit Flag(int* borrowed):p(borrowed),owned(false){}
    explicit Flag(size_t count=1){size_t bytes=std::max<size_t>(1,count)*sizeof(int);ck(cudaMalloc(&p,bytes));ck(cudaMemset(p,0,bytes));}
    ~Flag(){if(p&&owned)cudaFree(p);}
    void check(){int value;ck(cudaMemcpy(&value,p,sizeof(int),cudaMemcpyDeviceToHost));if(value)throw std::runtime_error("invalid GPU numeric data or CSR indices");}
};
struct Input {
    Flag validation;
    Input()=default;explicit Input(int* shared):validation(shared){}
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
uint64_t nativeBufferAllocations=0;
struct Buffer {
    double* p=nullptr; size_t count;
    explicit Buffer(size_t n):count(n){ ck(cudaMalloc(&p,std::max<size_t>(1,n)*sizeof(double)));++nativeBufferAllocations; }
    ~Buffer(){if(p)cudaFree(p);}
    Buffer(const Buffer&)=delete; Buffer& operator=(const Buffer&)=delete;
    void zero(){ck(cudaMemset(p,0,count*sizeof(double)));}
};
#include "ams_amg_kernels.cuh"
struct CSR {
    int n,m,nz;const int *rp,*ci;const double* a;
    CSR(const mxArray* v,Input& in):n(integer(field(v,"rows"))),m(integer(field(v,"cols"))) {
        if(n<1||m<1||n>INT_MAX-T||m>INT_MAX-T)throw std::runtime_error("operator dimensions exceed supported capacity");
        auto vals=field(v,"values");if(!mxIsGPUArray(vals))throw std::runtime_error("CSR values require GPU");
        auto temp=mxGPUCreateFromMxArray(vals);size_t count=mxGPUGetNumberOfElements(temp);mxGPUDestroyGPUArray(temp);
        if(count>INT_MAX-T)throw std::runtime_error("CSR index capacity exceeded");nz=int(count);
        rp=in.get<int>(field(v,"rowptr"),mxINT32_CLASS,size_t(n)+1);
        ci=in.get<int>(field(v,"columns"),mxINT32_CLASS,nz);a=in.get<double>(vals,mxDOUBLE_CLASS,nz);
        verify_csr<<<(std::max(n,nz)+T-1)/T,T,0,launch_stream>>>(n,m,nz,rp,ci,a,in.validation.p);
    }
    void mult(const double* x,double* y)const{multiply<<<(n+T-1)/T,T,0,launch_stream>>>(n,rp,ci,a,x,y);}
    void add(const double* x,double alpha,double* y)const{multiply_add<<<(n+T-1)/T,T,0,launch_stream>>>(n,rp,ci,a,x,alpha,y);}
    void relax(Buffer& x,const double* b,const double* l1,Buffer& tmp)const{
        multiply_relax<<<(n+T-1)/T,T,0,launch_stream>>>(n,rp,ci,a,x.p,b,l1,tmp.p);
        std::swap(x.p,tmp.p);
    }
    void res(const double* x,const double* b,double* r)const{
        multiply_residual<<<(n+T-1)/T,T,0,launch_stream>>>(n,rp,ci,a,x,b,r);
    }
};
void finite_vec(Input& in,const mxArray* a,int n,const double*& result,bool positive=false){
    result=in.get<double>(a,mxDOUBLE_CLASS,n);verify_vec<<<(n+T-1)/T,T,0,launch_stream>>>(n,result,positive,in.validation.p);
}
struct Coarse {
    std::string kind;int n,used=0;const double *a=nullptr,*di=nullptr;const int *order=nullptr,*rp=nullptr,*ci=nullptr;
    std::unique_ptr<Buffer> work;
    std::unique_ptr<Flag> seen;
    std::unique_ptr<Buffer> inverse;
    std::unique_ptr<Flag> reverse;
    Coarse(const mxArray* v,Input& in,int rows):kind(word(field(v,"kind"))),n(rows){
        if(integer(field(v,"rows"))!=n)throw std::runtime_error("coarse size mismatch");
        if(kind=="dense_inverse") {
            a=in.get<double>(field(v,"values"),mxDOUBLE_CLASS,size_t(n)*n);
            if(size_t(n)*n>INT_MAX-T)throw std::runtime_error("dense coarse capacity exceeded");
            verify_vec<<<(n*n+T-1)/T,T,0,launch_stream>>>(n*n,a,false,in.validation.p);
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
            seen=std::make_unique<Flag>(used);
            verify_csr<<<std::max(1,(std::max(used,int(count))+T-1)/T),T,0,launch_stream>>>(used,used,int(count),rp,ci,a,in.validation.p);
            verify_factor<<<(n+T-1)/T,T,0,launch_stream>>>(n,used,int(count),order,rp,ci,seen->p,in.validation.p);
            verify_permutation<<<std::max(1,(used+T-1)/T),T,0,launch_stream>>>(used,seen->p,in.validation.p);
            verify_vec<<<std::max(1,(used+T-1)/T),T,0,launch_stream>>>(used,di,false,in.validation.p);work=std::make_unique<Buffer>(used);
            // Validate every factor/index before construction dereferences it.
            in.validation.check();
            if(used>0&&used<=512){size_t available,total;ck(cudaMemGetInfo(&available,&total));size_t bytes=size_t(used)*used*sizeof(double)+size_t(used)*sizeof(int);
                size_t reserve=std::max<size_t>(2ULL<<30,total/5);if(available<=reserve||bytes>available-reserve)throw std::runtime_error("coarse action cache exceeds CUDA memory reserve");
                inverse=std::make_unique<Buffer>(size_t(used)*used);reverse=std::make_unique<Flag>(used);
                cache_inverse_order<<<(n+T-1)/T,T>>>(n,order,reverse->p);
                inverse_sparse_coarse<<<used,32>>>(used,rp,ci,a,di,inverse->p);
                verify_vec<<<(used*used+T-1)/T,T>>>(used*used,inverse->p,false,in.validation.p);in.validation.check();
            }
        } else if(kind!="ten_l1_jacobi_sweeps")throw std::runtime_error("unknown coarse solver");
    }
    void apply(const double* b,double* x)const{
        if(kind=="dense_inverse")dense_coarse<<<(n+T-1)/T,T,0,launch_stream>>>(n,a,b,x);
        else if(inverse)cached_sparse_coarse<<<(n+T-1)/T,T,0,launch_stream>>>(n,used,order,reverse->p,inverse->p,b,x);
        else sparse_coarse<<<1,32,0,launch_stream>>>(n,used,order,rp,ci,a,di,b,x,work->p);
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
    void relax(){a.relax(x,b.p,l1,tmp);}
    void initial(){initial_relax<<<(a.n+T-1)/T,T,0,launch_stream>>>(a.n,b.p,l1,x.p);}
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
        auto& l=*levels[k];
        if(l.coarse){if(l.coarse->kind=="ten_l1_jacobi_sweeps"){l.initial();for(int s=1;s<10;++s)l.relax();}else l.coarse->apply(l.b.p,l.x.p);return;}
        l.initial();for(int s=1;s<l.smooth;++s)l.relax();
        l.a.res(l.x.p,l.b.p,l.r.p);
        l.pt->mult(l.r.p,levels[k+1]->b.p);cycle(k+1);l.p->add(levels[k+1]->x.p,1,l.x.p);
        for(int s=0;s<l.smooth;++s)l.relax();
    }
    void apply(const double* residual){pt.mult(residual,levels[0]->b.p);cycle(0);p.mult(levels[0]->x.p,out.p);}
    void add(const double* residual,double alpha,double* output){pt.mult(residual,levels[0]->b.p);cycle(0);p.add(levels[0]->x.p,alpha,output);}
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
    void res(){a.res(x.p,b.p,r.p);}
    size_t coarseCacheBytes()const{size_t bytes=0;for(const auto& s:spaces)for(const auto& l:s->levels)if(l->coarse&&l->coarse->inverse)bytes+=l->coarse->inverse->count*sizeof(double)+size_t(l->coarse->used)*sizeof(int);return bytes;}
    void relax(){for(int s=0;s<smooth;++s)a.relax(x,b.p,l1,tmp);}
    void initial(){initial_relax<<<(a.n+T-1)/T,T,0,launch_stream>>>(a.n,b.p,l1,x.p);for(int s=1;s<smooth;++s)a.relax(x,b.p,l1,tmp);}
    void grad(){if(!gradient)return;res();spaces[0]->add(r.p,omega,x.p);}
    void nodal(){res();for(int k=gradient;k<gradient+3;++k)spaces[k]->add(r.p,omega,x.p);}
    void apply(const double* rhs,double* output){ck(cudaMemcpyAsync(b.p,rhs,a.n*sizeof(double),cudaMemcpyDeviceToDevice,launch_stream));mask<<<(a.n+T-1)/T,T,0,launch_stream>>>(a.n,free,b.p);
        if(cycle==1){initial();grad();nodal();grad();relax();}
        else {initial();grad();relax();nodal();relax();grad();relax();}
        mask<<<(a.n+T-1)/T,T,0,launch_stream>>>(a.n,free,x.p);ck(cudaMemcpyAsync(output,x.p,a.n*sizeof(double),cudaMemcpyDeviceToDevice,launch_stream));}
};
// Record the fixed AMS cycle once per solve. All graph buffers remain owned by
// AMS and the PCG workspace until the graph is destroyed. Capture does not run
// kernels; each replay initializes its iterates before consuming them.
struct AMSGraph {
    cudaStream_t stream=nullptr;
    cudaGraph_t graph=nullptr;
    cudaGraphExec_t executable=nullptr;
    AMSGraph(AMS& ams,const double* rhs,double* output) {
        ck(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
        try {
            ck(cudaStreamBeginCapture(stream,cudaStreamCaptureModeThreadLocal));
            launch_stream=stream;
            try { ams.apply(rhs,output); }
            catch(...) { launch_stream=nullptr; cudaStreamEndCapture(stream,&graph); throw; }
            launch_stream=nullptr;
            ck(cudaStreamEndCapture(stream,&graph));
            ck(cudaGraphInstantiate(&executable,graph,nullptr,nullptr,0));
        } catch(...) { release(); throw; }
    }
    void release() noexcept {
        if(executable)cudaGraphExecDestroy(executable);
        if(graph)cudaGraphDestroy(graph);
        if(stream)cudaStreamDestroy(stream);
    }
    ~AMSGraph(){release();}
    AMSGraph(const AMSGraph&)=delete;
    AMSGraph& operator=(const AMSGraph&)=delete;
    void apply(){ck(cudaGraphLaunch(executable,nullptr));}
};
struct Blas {
    cublasHandle_t h;
    Blas(){cb(cublasCreate(&h));}
    ~Blas(){cublasDestroy(h);}
    void dot(int n,const double* a,const double* b,double* device_result){
        cb(cublasSetPointerMode(h,CUBLAS_POINTER_MODE_DEVICE));
        cb(cublasDdot(h,n,a,1,b,1,device_result));
    }
    double norm(int n,const double* a){
        cb(cublasSetPointerMode(h,CUBLAS_POINTER_MODE_HOST));
        double d;cb(cublasDnrm2(h,n,a,1,&d));return d;
    }
};
struct Output {mxGPUArray* a;explicit Output(int n){mwSize s=n;a=mxGPUCreateGPUArray(1,&s,mxDOUBLE_CLASS,mxREAL,MX_GPU_INITIALIZE_VALUES);}~Output(){mxGPUDestroyGPUArray(a);}double* data(){return static_cast<double*>(mxGPUGetData(a));}};
#include "ams_resident.cuh"
#include "ams_nonlinear.cuh"
#include "ams_fem.cuh"
void run(int nlhs,mxArray** out,int nrhs,const mxArray** args){
    if(nrhs>0 && mxIsChar(args[0])) {
        auto command=word(args[0]);
        if(command=="setup"){RadiaAMSSetup(nlhs,out,nrhs,args);return;}
        mxInitGPU();
        if(femCommand(command,nlhs,out,nrhs,args)||amsResidentCommand(command,nlhs,out,nrhs,args)||nonlinearCommand(command,nlhs,out,nrhs,args))return;
        throw std::runtime_error("Unknown AMS command");
    }
    if((nrhs<5||nrhs>7)||nlhs!=2)throw std::runtime_error("[x,info]=radia_cuda_ams(snapshot,b,tol,maxiter,apply_only[,use_graph])");
    if(nrhs>=6&&!mxIsLogicalScalar(args[5]))throw std::runtime_error("logical use_graph required");
    const bool use_graph=nrhs>=6&&mxIsLogicalScalarTrue(args[5]);
    double tol=scalar(args[2]);int maxiter=integer(args[3]);
    if(!(tol>0&&tol<=1e-6)||maxiter<1||maxiter==INT_MAX||!mxIsLogicalScalar(args[4]))throw std::runtime_error("invalid tolerance, iteration limit or apply flag");
    mxInitGPU();size_t memory_before,total,available;ck(cudaMemGetInfo(&memory_before,&total));
    Input in;AMS ams(args[0],in);int n=ams.a.n;const double* b;
    // Frozen AMS preconditions the current tangent; residuals use that tangent.
    std::unique_ptr<CSR> current;
    if(nrhs==7){current=std::make_unique<CSR>(args[6],in);
        if(current->n!=n||current->m!=n)throw std::runtime_error("current operator dimensions mismatch");}
    const CSR& op=current?*current:ams.a;
    ck(cudaMemGetInfo(&available,&total));size_t workspace_sample=memory_before>available?memory_before-available:0;
    finite_vec(in,args[1],n,b);bool apply_only=mxIsLogicalScalarTrue(args[4]);
    Buffer zero(n);zero.zero();verify_mask<<<(n+T-1)/T,T,0,launch_stream>>>(n,ams.free,apply_only?zero.p:b,in.validation.p);
    in.validation.check(); // No solve may consume unvalidated numerical data.
    Output solution(n);double rel=apply_only?mxGetNaN():0;int it=0;
    if(apply_only){
        ams.apply(b,solution.data());
        verify_vec<<<(n+T-1)/T,T,0,launch_stream>>>(n,solution.data(),false,in.validation.p);
        in.validation.check(); // Finite inputs can still overflow during an AMS cycle.
    }
    else {
        ResidentWork work(n);std::unique_ptr<AMSGraph> graph;
        if(use_graph)graph=std::make_unique<AMSGraph>(ams,work.r.p,work.z.p);
        auto stats=residentPCG(ams,op,b,work,tol,maxiter,graph.get());
        it=stats.iterations;rel=stats.residual;
        workspace_sample=std::max(workspace_sample,work.bytes());
        ck(cudaMemcpy(solution.data(),work.x.p,size_t(n)*sizeof(double),cudaMemcpyDeviceToDevice));
    }
    ck(cudaGetLastError());ck(cudaDeviceSynchronize());out[0]=mxGPUCreateMxArrayOnGPU(solution.a);
    const char* names[]={"iterations","true_relative_residual","coarse_execution","setup_execution","workspace_sample_bytes","residual_checked","preconditioner_execution"};
    out[1]=mxCreateStructMatrix(1,1,7,names);mxSetField(out[1],0,names[0],mxCreateDoubleScalar(it));mxSetField(out[1],0,names[1],mxCreateDoubleScalar(rel));
    mxSetField(out[1],0,names[2],mxCreateString("GPU: existing sparse LDL factors with warp-parallel row products, or existing dense inverse"));
    mxSetField(out[1],0,names[3],mxCreateString("Input hierarchy: MATLAB/mexcuda host setup or explicit validation fixture; no runtime Python fallback"));
    mxSetField(out[1],0,names[4],mxCreateDoubleScalar(double(workspace_sample)));
    mxSetField(out[1],0,names[5],mxCreateLogicalScalar(!apply_only));
    mxSetField(out[1],0,names[6],mxCreateString((!apply_only&&use_graph)?"CUDA Graph replay of the unchanged AMS cycle":"direct CUDA kernels"));
}
}
// These definitions are linked by the separately compiled MATLAB setup
// adapter. Keep external linkage, outside the private solver namespace.
#include "ams_setup_gpu.cuh"
void mexFunction(int nl,mxArray** out,int nr,const mxArray** in){

    char errorMessage[2048]={0};
    if(!radia_mex::checked(run,nl,out,nr,in,errorMessage,sizeof(errorMessage))){
        mexErrMsgIdAndTxt("radia:cuda:ams","%s",errorMessage);
    }
}
