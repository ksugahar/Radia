#include "mex_error_boundary.hpp"
// CUDA-only SPD CSR solve. NGSolve owns finite-element assembly and numbering.
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <cuda_runtime.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../../src/core/gpu/csr_pcg_cuda.cuh"
namespace {
struct Inputs {
    std::vector<const mxGPUArray*> arrays;
    ~Inputs() { for (auto a : arrays) mxGPUDestroyGPUArray(a); }
    const mxGPUArray* get(const mxArray* arg, mxClassID type) {
        if (!mxIsGPUArray(arg)) throw std::runtime_error("All system arrays must be gpuArray; no CPU fallback");
        auto a = mxGPUCreateFromMxArray(arg);
        try { arrays.push_back(a); } catch(...) { mxGPUDestroyGPUArray(a); throw; }
        if (mxGPUGetClassID(a) != type || mxGPUGetComplexity(a) != mxREAL)
            throw std::runtime_error("Expected real int32 CSR indices and real double numerical arrays");
        return a;
    }
};
struct GPUOutput {
    mxGPUArray* array;
    explicit GPUOutput(mwSize n) : array(mxGPUCreateGPUArray(1, &n, mxDOUBLE_CLASS, mxREAL, MX_GPU_INITIALIZE_VALUES)) {}
    ~GPUOutput() { mxGPUDestroyGPUArray(array); }
};
double scalar(const mxArray* a) {
    if(!mxIsDouble(a) || mxIsComplex(a) || mxGetNumberOfElements(a)!=1)
        throw std::runtime_error("Tolerance and iteration limit must be real double scalars");
    return mxGetScalar(a);
}
}

void dispatchCsr(int nl,mxArray** out,int nr,const mxArray** in) {
    if(nl!=2 || (nr!=7 && nr!=8))throw std::runtime_error("CSR solve requires seven arrays/options and optional backend");
    bool cpu=false;
    if(nr==8){
        if(!mxIsChar(in[7]) || mxGetM(in[7])!=1)throw std::runtime_error("Backend must be cpu or gpu character row");
        char backend[8]={0};
        if(mxGetString(in[7],backend,sizeof(backend))!=0)throw std::runtime_error("Invalid Backend");
        if(std::string(backend)=="cpu")cpu=true;
        else if(std::string(backend)!="gpu")throw std::runtime_error("Backend must be cpu or gpu");
    }
    const double tol=scalar(in[5]),limit=scalar(in[6]);
    if(!std::isfinite(limit) || limit<1 || limit>100000 || limit!=std::floor(limit))
        throw std::runtime_error("Invalid iteration limit");
    radia_csr::Result result{};
    if(cpu){
        for(int i=0;i<5;++i){
            const mxClassID type=i<2?mxINT32_CLASS:mxDOUBLE_CLASS;
            if(mxGetClassID(in[i])!=type || mxIsSparse(in[i]) || mxIsComplex(in[i]))
                throw std::runtime_error("CPU arrays must be full real int32 indices or double values");
        }
        const size_t n=mxGetNumberOfElements(in[3]),nnz=mxGetNumberOfElements(in[2]);
        radia_csr::validate_dimensions(n,nnz,tol,static_cast<int>(limit));
        if(mxGetNumberOfElements(in[0])!=n+1 || mxGetNumberOfElements(in[1])!=nnz || mxGetNumberOfElements(in[4])!=n)
            throw std::runtime_error("CSR dimensions disagree");
        // All C++ allocations unwind before the POD-only mexFunction reports errors.
        out[0]=mxCreateDoubleMatrix(n,1,mxREAL);
        result=radia_csr::solve_cpu(static_cast<int>(n),static_cast<int>(nnz),
            static_cast<const int*>(mxGetData(in[0])),static_cast<const int*>(mxGetData(in[1])),
            mxGetPr(in[2]),mxGetPr(in[3]),mxGetPr(in[4]),tol,static_cast<int>(limit),mxGetPr(out[0]));
    }else{
        mxInitGPU();
        Inputs inputs;
        auto rp=inputs.get(in[0],mxINT32_CLASS),ci=inputs.get(in[1],mxINT32_CLASS);
        auto av=inputs.get(in[2],mxDOUBLE_CLASS),b=inputs.get(in[3],mxDOUBLE_CLASS),di=inputs.get(in[4],mxDOUBLE_CLASS);
        const size_t n=mxGPUGetNumberOfElements(b),nnz=mxGPUGetNumberOfElements(av);
        radia_csr::validate_dimensions(n,nnz,tol,static_cast<int>(limit));
        if(mxGPUGetNumberOfElements(rp)!=n+1 || mxGPUGetNumberOfElements(ci)!=nnz || mxGPUGetNumberOfElements(di)!=n)
            throw std::runtime_error("CSR dimensions disagree");
        GPUOutput output(n);
        result=radia_csr::solve_device(static_cast<int>(n),static_cast<int>(nnz),
            static_cast<const int*>(mxGPUGetDataReadOnly(rp)),static_cast<const int*>(mxGPUGetDataReadOnly(ci)),
            static_cast<const double*>(mxGPUGetDataReadOnly(av)),static_cast<const double*>(mxGPUGetDataReadOnly(b)),
            static_cast<const double*>(mxGPUGetDataReadOnly(di)),tol,static_cast<int>(limit),static_cast<double*>(mxGPUGetData(output.array)));
        out[0]=mxGPUCreateMxArrayOnGPU(output.array);
    }
    const char* names[]={"iterations","true_relative_residual","precision","backend","preconditioner","host_threads"};
    out[1]=mxCreateStructMatrix(1,1,6,names);
    mxSetField(out[1],0,names[0],mxCreateDoubleScalar(result.iterations));
    mxSetField(out[1],0,names[1],mxCreateDoubleScalar(result.true_relative_residual));
    mxSetField(out[1],0,names[2],mxCreateString("float64"));
    mxSetField(out[1],0,names[3],mxCreateString(cpu?"cpu":"gpu"));
    mxSetField(out[1],0,names[4],mxCreateString("Jacobi; not IC or AMS"));
    mxSetField(out[1],0,names[5],mxCreateDoubleScalar(1));
}
void mexFunction(int nl,mxArray** out,int nr,const mxArray** in){
    char errorMessage[2048]={0};
    if(!radia_mex::checked(dispatchCsr,nl,out,nr,in,errorMessage,sizeof(errorMessage))){
        mexErrMsgIdAndTxt("radia:cuda:solve","%s",errorMessage);
    }
}
