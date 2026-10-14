// Resident complex BiCGStab with device reductions and batched graph replay.
#pragma once
#include "ih_resident_kernels.cuh"
#include "ih_bicgstab_kernels.cuh"
#include "ih_preconditioner.cuh"
#include <cublas_v2.h>
#include <stdexcept>
#include <limits>
#include <algorithm>
namespace radia_ih_cuda {
// A saddle-point RHS can give r^H A r == 0 at the first BiCG step.
// BiCGStab permits any fixed shadow with nonzero shadow^H r. Use a
// deterministic dense shadow perturbation; certify the original equation.
__global__ void initializeShadow(const cuDoubleComplex* rhs,cuDoubleComplex* shadow,
        const cuDoubleComplex* squaredNorm,int n) {
    const int i=int(blockIdx.x*blockDim.x+threadIdx.x);if(i>=n)return;
    const double scale=sqrt(cuCreal(*squaredNorm)/double(n));
    const double fraction=double(i+1)/double(n+1);
    shadow[i]=cuCadd(rhs[i],make_cuDoubleComplex(scale*(1.+fraction),scale*.25*fraction));
}
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
    ComplexPreconditioner preconditioner_{};
    cuDoubleComplex *phat_=nullptr,*shat_=nullptr;
    cuDoubleComplex *rhs_,*solution_,*r_,*r0_,*p_,*v_,*s_,*t_,*ax_;
    int n_=0;
    void dot(const cuDoubleComplex* a,const cuDoubleComplex* b,cuDoubleComplex* out) {
        blasCheck(cublasZdotc(blas_,n_,a,1,b,1,out));
    }
    void iteration() {
        const unsigned blocks=(unsigned(n_)+127)/128;
        dot(r0_,r_,&scalar_->rho);
        bicgDirection<<<blocks,128,0,stream_>>>(scalar_,r_,v_,p_,n_);
        const cuDoubleComplex* actionP=p_;
        if(preconditioner_.precision!=PreconditionerPrecision::None){
            applyFixedPreconditioner<<<blocks,128,0,stream_>>>(preconditioner_,p_,phat_,scalar_);actionP=phat_;
        }
        complexMatvec<<<blocks,128,0,stream_>>>(matrix_,actionP,v_);
        dot(r0_,v_,&scalar_->denominator);
        bicgAlpha<<<1,1,0,stream_>>>(scalar_);
        bicgIntermediate<<<blocks,128,0,stream_>>>(scalar_,r_,v_,s_,n_);
        dot(s_,s_,rr_);bicgCheckIntermediate<<<1,1,0,stream_>>>(scalar_,rr_);
        const cuDoubleComplex* actionS=s_;
        if(preconditioner_.precision!=PreconditionerPrecision::None){
            applyFixedPreconditioner<<<blocks,128,0,stream_>>>(preconditioner_,s_,shat_,scalar_);actionS=shat_;
        }
        complexMatvec<<<blocks,128,0,stream_>>>(matrix_,actionS,t_);
        dot(t_,s_,&scalar_->ts);dot(t_,t_,&scalar_->tt);
        bicgOmega<<<1,1,0,stream_>>>(scalar_);
        bicgUpdateRight<<<blocks,128,0,stream_>>>(scalar_,actionP,actionS,s_,t_,solution_,r_,n_);
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
    explicit ComplexSolver(ComplexCSRView matrix,ComplexPreconditioner preconditioner={}):matrix_(matrix),preconditioner_(preconditioner) {
        if(!matrix.rows||matrix.rows>std::size_t(std::numeric_limits<int>::max()))throw std::invalid_argument("IH solver dimension out of range");
        n_=int(matrix.rows);
        if(preconditioner.precision!=PreconditionerPrecision::None&&
            (preconditioner.rows!=matrix.rows||!preconditioner.rowptr||!preconditioner.columns||
             (preconditioner.precision==PreconditionerPrecision::FP32?!preconditioner.values32:!preconditioner.values64)))
            throw std::invalid_argument("Invalid fixed IH preconditioner");
        try {
            cudaCheck(cudaStreamCreateWithFlags(&stream_,cudaStreamNonBlocking));
            blasCheck(cublasCreate(&blas_));blasCheck(cublasSetStream(blas_,stream_));
            blasCheck(cublasSetPointerMode(blas_,CUBLAS_POINTER_MODE_DEVICE));
            cudaCheck(cudaMalloc(&workspace_,workspaceBytes()));
            cudaCheck(cudaMalloc(&scalar_,sizeof(BiCGState)));
            rhs_=workspace_;solution_=rhs_+n_;r_=solution_+n_;r0_=r_+n_;
            p_=r0_+n_;v_=p_+n_;s_=v_+n_;t_=s_+n_;ax_=t_+n_;rr_=ax_+n_;
            cudaCheck(cudaMemsetAsync(workspace_,0,workspaceBytes(),stream_));
            if(preconditioner_.precision!=PreconditionerPrecision::None){phat_=rr_+1;shat_=phat_+n_;}
            dot(rhs_,rhs_,rr_);cudaCheck(cudaStreamSynchronize(stream_));
            cudaCheck(cudaStreamBeginCapture(stream_,cudaStreamCaptureModeThreadLocal));
            for(int i=0;i<8;++i)iteration();
            cudaCheck(cudaStreamEndCapture(stream_,&graph_));
            cudaCheck(cudaGraphInstantiate(&executable_,graph_,nullptr,nullptr,0));
        } catch(...) {release();throw;}
    }
    // Only the checked owner may call this after CUDA confirms a device reset.
    // Stale device addresses must never be freed into a new CUDA context.
    void abandonAfterDeviceReset() noexcept {
        stream_=nullptr;blas_=nullptr;graph_=nullptr;executable_=nullptr;
        scalar_=nullptr;workspace_=nullptr;
    }
    ~ComplexSolver(){release();}
    cudaStream_t stream()const{return stream_;}
    cuDoubleComplex* rhs(){return rhs_;}
    const cuDoubleComplex* solution()const{return solution_;}
    std::size_t workspaceBytes()const{return (std::size_t(n_)*(preconditioner_.precision==PreconditionerPrecision::None?9:11)+1)*sizeof(cuDoubleComplex);}
    std::size_t ownedBytes()const{return workspaceBytes()+sizeof(BiCGState);}
    struct Result {unsigned iterations;double trueRelativeResidual;};
    Result solve(double tolerance,int maxIterations) {
        if(!std::isfinite(tolerance)||tolerance<=0.||tolerance>1e-6||maxIterations<1)throw std::invalid_argument("Invalid IH solve limits");
        cudaCheck(cudaMemsetAsync(solution_,0,std::size_t(n_)*sizeof(cuDoubleComplex),stream_));
        cudaCheck(cudaMemcpyAsync(r_,rhs_,std::size_t(n_)*sizeof(cuDoubleComplex),cudaMemcpyDeviceToDevice,stream_));
        cudaCheck(cudaMemcpyAsync(r0_,rhs_,std::size_t(n_)*sizeof(cuDoubleComplex),cudaMemcpyDeviceToDevice,stream_));
        dot(rhs_,rhs_,rr_);
        initializeShadow<<<(unsigned(n_)+127)/128,128,0,stream_>>>(rhs_,r0_,rr_,n_);
        initializeBiCG<<<1,1,0,stream_>>>(scalar_,rr_,tolerance);
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
