// Resident complex BiCGStab with device reductions and batched graph replay.
#pragma once
#include "ih_resident_kernels.cuh"
#include "ih_bicgstab_kernels.cuh"
#include <cublas_v2.h>
#include <stdexcept>
#include <limits>
#include <algorithm>
namespace radia_ih_cuda {
inline void cudaCheck(cudaError_t e){if(e!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(e));}
inline void blasCheck(cublasStatus_t e){if(e!=CUBLAS_STATUS_SUCCESS)throw std::runtime_error("IH cuBLAS operation failed");}
__global__ void initializeBiCG(BiCGState* scalar,const cuDoubleComplex* norm,double tolerance) {
    *scalar=BiCGState{};
    scalar->rho=scalar->previousRho=scalar->alpha=scalar->omega=make_cuDoubleComplex(1.,0.);
    scalar->rhsNormSquared=cuCreal(*norm);scalar->relativeToleranceSquared=tolerance*tolerance;
    if(!isfinite(scalar->rhsNormSquared)||scalar->rhsNormSquared<0.)scalar->status=-1;
    else if(scalar->rhsNormSquared==0.)scalar->status=1;
}
class ComplexSolver {
    cudaStream_t stream_=nullptr;
    cublasHandle_t blas_=nullptr;
    cudaGraph_t graph_=nullptr;
    cudaGraphExec_t executable_=nullptr;
    cuDoubleComplex *workspace_=nullptr,*rr_=nullptr;
    BiCGState* scalar_=nullptr;
    ComplexCSRView matrix_{};
    cuDoubleComplex *rhs_,*solution_,*r_,*r0_,*p_,*v_,*s_,*t_,*ax_;
    int n_=0;
    void dot(const cuDoubleComplex* a,const cuDoubleComplex* b,cuDoubleComplex* out) {
        blasCheck(cublasZdotc(blas_,n_,a,1,b,1,out));
    }
    void iteration() {
        const unsigned blocks=(unsigned(n_)+127)/128;
        dot(r0_,r_,&scalar_->rho);
        bicgDirection<<<blocks,128,0,stream_>>>(scalar_,r_,v_,p_,n_);
        complexMatvec<<<blocks,128,0,stream_>>>(matrix_,p_,v_);
        dot(r0_,v_,&scalar_->denominator);
        bicgAlpha<<<1,1,0,stream_>>>(scalar_);
        bicgIntermediate<<<blocks,128,0,stream_>>>(scalar_,r_,v_,s_,n_);
        dot(s_,s_,rr_);bicgCheckIntermediate<<<1,1,0,stream_>>>(scalar_,rr_);
        complexMatvec<<<blocks,128,0,stream_>>>(matrix_,s_,t_);
        dot(t_,s_,&scalar_->ts);dot(t_,t_,&scalar_->tt);
        bicgOmega<<<1,1,0,stream_>>>(scalar_);
        bicgUpdate<<<blocks,128,0,stream_>>>(scalar_,p_,s_,t_,solution_,r_,n_);
        dot(r_,r_,rr_);bicgFinish<<<1,1,0,stream_>>>(scalar_,rr_);
    }
    void release() noexcept {
        if(stream_)cudaStreamSynchronize(stream_);
        if(executable_)cudaGraphExecDestroy(executable_);
        if(graph_)cudaGraphDestroy(graph_);
        if(blas_)cublasDestroy(blas_);
        if(scalar_)cudaFree(scalar_);if(workspace_)cudaFree(workspace_);
        if(stream_)cudaStreamDestroy(stream_);
    }
public:
    ComplexSolver(const ComplexSolver&)=delete;
    ComplexSolver& operator=(const ComplexSolver&)=delete;
    explicit ComplexSolver(ComplexCSRView matrix):matrix_(matrix) {
        if(!matrix.rows||matrix.rows>std::size_t(std::numeric_limits<int>::max()))throw std::invalid_argument("IH solver dimension out of range");
        n_=int(matrix.rows);
        try {
            cudaCheck(cudaStreamCreateWithFlags(&stream_,cudaStreamNonBlocking));
            blasCheck(cublasCreate(&blas_));blasCheck(cublasSetStream(blas_,stream_));
            blasCheck(cublasSetPointerMode(blas_,CUBLAS_POINTER_MODE_DEVICE));
            cudaCheck(cudaMalloc(&workspace_,(std::size_t(n_)*9+1)*sizeof(cuDoubleComplex)));
            cudaCheck(cudaMalloc(&scalar_,sizeof(BiCGState)));
            rhs_=workspace_;solution_=rhs_+n_;r_=solution_+n_;r0_=r_+n_;
            p_=r0_+n_;v_=p_+n_;s_=v_+n_;t_=s_+n_;ax_=t_+n_;rr_=ax_+n_;
            cudaCheck(cudaMemsetAsync(workspace_,0,(std::size_t(n_)*9+1)*sizeof(cuDoubleComplex),stream_));
            dot(rhs_,rhs_,rr_);cudaCheck(cudaStreamSynchronize(stream_));
            cudaCheck(cudaStreamBeginCapture(stream_,cudaStreamCaptureModeThreadLocal));
            for(int i=0;i<8;++i)iteration();
            cudaCheck(cudaStreamEndCapture(stream_,&graph_));
            cudaCheck(cudaGraphInstantiate(&executable_,graph_,nullptr,nullptr,0));
        } catch(...) {release();throw;}
    }
    ~ComplexSolver(){release();}
    cudaStream_t stream()const{return stream_;}
    cuDoubleComplex* rhs(){return rhs_;}
    const cuDoubleComplex* solution()const{return solution_;}
    std::size_t ownedBytes()const{return (std::size_t(n_)*9+1)*sizeof(cuDoubleComplex)+sizeof(BiCGState);}
    struct Result {unsigned iterations;double trueRelativeResidual;};
    Result solve(double tolerance,int maxIterations) {
        if(!std::isfinite(tolerance)||tolerance<=0.||tolerance>1e-6||maxIterations<1)throw std::invalid_argument("Invalid IH solve limits");
        cudaCheck(cudaMemsetAsync(solution_,0,std::size_t(n_)*sizeof(cuDoubleComplex),stream_));
        cudaCheck(cudaMemcpyAsync(r_,rhs_,std::size_t(n_)*sizeof(cuDoubleComplex),cudaMemcpyDeviceToDevice,stream_));
        cudaCheck(cudaMemcpyAsync(r0_,rhs_,std::size_t(n_)*sizeof(cuDoubleComplex),cudaMemcpyDeviceToDevice,stream_));
        dot(rhs_,rhs_,rr_);initializeBiCG<<<1,1,0,stream_>>>(scalar_,rr_,tolerance);
        BiCGState state{};
        for(int i=0;i<maxIterations;i+=8) {
            if(maxIterations-i>=8)cudaCheck(cudaGraphLaunch(executable_,stream_));
            else for(int j=i;j<maxIterations;++j)iteration();
            cudaCheck(cudaMemcpyAsync(&state,scalar_,sizeof(state),cudaMemcpyDeviceToHost,stream_));
            cudaCheck(cudaStreamSynchronize(stream_));
            if(state.status!=0)break;
        }
        if(state.status!=1)throw std::runtime_error("IH complex solve did not converge or broke down");
        const unsigned blocks=(unsigned(n_)+127)/128;
        complexMatvec<<<blocks,128,0,stream_>>>(matrix_,solution_,ax_);
        trueResidual<<<blocks,128,0,stream_>>>(rhs_,ax_,r_,n_);dot(r_,r_,rr_);
        cuDoubleComplex norm{};cudaCheck(cudaMemcpyAsync(&norm,rr_,sizeof(norm),cudaMemcpyDeviceToHost,stream_));
        cudaCheck(cudaStreamSynchronize(stream_));
        const double residual=state.rhsNormSquared==0.?std::sqrt(cuCreal(norm)):std::sqrt(cuCreal(norm)/state.rhsNormSquared);
        if(!std::isfinite(residual)||residual>tolerance)throw std::runtime_error("IH complex true residual exceeds tolerance");
        return {state.iteration,residual};
    }
};
} // namespace radia_ih_cuda
