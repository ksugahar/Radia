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
#include <atomic>
#include <climits>
#include <memory>
#include <string>
#include <stdexcept>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <cstdint>
#include "ams_fem_host.hpp"
#include "../../../src/core/gpu/snapshot.hpp"
void RadiaAMSSetup(int,mxArray**,int,const mxArray**);
namespace {
#include "../../../src/core/gpu/ams_core.cuh"
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
Node device_node(const mxArray* a,Input& in){
    Node v;
    if(mxIsGPUArray(a)){
        auto p=mxGPUCreateFromMxArray(a);in.refs.push_back(p);
        if(mxGPUGetComplexity(p)!=mxREAL)throw std::runtime_error("real GPU snapshot array required");
        auto type=mxGPUGetClassID(p);
        if(type==mxINT32_CLASS)v.kind=Node::integers;
        else if(type==mxDOUBLE_CLASS)v.kind=Node::reals;
        else throw std::runtime_error("typed int32/FP64 GPU snapshot required");
        v.device=mxGPUGetDataReadOnly(p);v.count=mxGPUGetNumberOfElements(p);
    }else if(mxIsStruct(a)){
        if(mxGetNumberOfElements(a)!=1)throw std::runtime_error("scalar snapshot struct required");
        v.kind=Node::record;
        for(int i=0;i<mxGetNumberOfFields(a);++i){auto name=mxGetFieldNameByNumber(a,i);v.fields[name]=device_node(field(a,name),in);}
    }else if(mxIsCell(a)){
        v.kind=Node::list;for(size_t i=0;i<mxGetNumberOfElements(a);++i)v.children.push_back(device_node(mxGetCell(a,i),in));
    }else if(mxIsLogicalScalar(a)){v.kind=Node::logical;v.boolean=mxIsLogicalScalarTrue(a);}
    else if(mxIsChar(a)){v.kind=Node::text;v.word=word(a);}
    else{v.kind=Node::number;v.scalar=scalar(a);}
    return v;
}
void finite_vec(Input& in,const mxArray* a,int n,const double*& result,bool positive=false){
    result=in.get<double>(a,mxDOUBLE_CLASS,n);
    verify_vec<<<(n+T-1)/T,T,0,launch_stream>>>(n,result,positive,in.validation.p);
}
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
    Input in;AMS ams(device_node(args[0],in),in.validation);int n=ams.a.n;const double* b;
    // Frozen AMS preconditions the current tangent; residuals use that tangent.
    std::unique_ptr<CSR> current;
    if(nrhs==7){current=std::make_unique<CSR>(device_node(args[6],in),in.validation);
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
