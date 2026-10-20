// Fixed right preconditioners only. Outer vectors/operator/residual stay FP64.
#pragma once
#include "ih_resident_kernels.cuh"
#include "ih_bicgstab_kernels.cuh"
namespace radia_ih_cuda {
// Bounded frozen block-Jacobi setup on device. Pivoting changes only the
// preconditioner; the outer operator and true residual remain original FP64.
__global__ void buildBlockInverse(ComplexCSRView a,const int* ends,int blocks,
        cuDoubleComplex* work,cuDoubleComplex* inverse,int* failure) {
    const int block=int(blockIdx.x);if(block>=blocks||threadIdx.x)return;
    const int begin=block?ends[block-1]:0,end=ends[block],n=int(a.rows),m=end-begin;
    for(int i=begin;i<end;++i){
        for(int j=begin;j<end;++j){work[i*n+j]=make_cuDoubleComplex(0.,0.);inverse[i*n+j]=make_cuDoubleComplex(i==j?1.:0.,0.);}
        for(int p=a.rowptr[i];p<a.rowptr[i+1];++p){int j=a.columns[p];if(j>=begin&&j<end)work[i*n+j]=a.values[p];}
    }
    for(int k=begin;k<end;++k){
        int pivot=k;double best=0.;for(int i=k;i<end;++i){double v=cuCabs(work[i*n+k]);if(v>best){best=v;pivot=i;}}
        if(!isfinite(best)||best==0.){atomicExch(failure,1);return;}
        if(pivot!=k)for(int j=begin;j<end;++j){auto v=work[k*n+j];work[k*n+j]=work[pivot*n+j];work[pivot*n+j]=v;v=inverse[k*n+j];inverse[k*n+j]=inverse[pivot*n+j];inverse[pivot*n+j]=v;}
        auto diagonal=work[k*n+k];for(int j=begin;j<end;++j){work[k*n+j]=cuCdiv(work[k*n+j],diagonal);inverse[k*n+j]=cuCdiv(inverse[k*n+j],diagonal);}
        for(int i=begin;i<end;++i)if(i!=k){auto factor=work[i*n+k];for(int j=begin;j<end;++j){work[i*n+j]=cuCsub(work[i*n+j],cuCmul(factor,work[k*n+j]));inverse[i*n+j]=cuCsub(inverse[i*n+j],cuCmul(factor,inverse[k*n+j]));}}
    }
    for(int i=begin;i<end;++i)for(int j=begin;j<end;++j)if(!isfinite(cuCreal(inverse[i*n+j]))||!isfinite(cuCimag(inverse[i*n+j])))atomicExch(failure,1);
}
enum class PreconditionerPrecision { None, FP64, FP32 };
struct ComplexPreconditioner {
    std::size_t rows=0;
    const std::int32_t *rowptr=nullptr,*columns=nullptr;
    const cuDoubleComplex* values64=nullptr;
    const cuFloatComplex* values32=nullptr;
    PreconditionerPrecision precision=PreconditionerPrecision::None;
};
__global__ void applyFixedPreconditioner(ComplexPreconditioner preconditioner,
        const cuDoubleComplex* input,cuDoubleComplex* output,BiCGState* scalar) {
    const std::size_t row=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(row>=preconditioner.rows)return;
    // No early return based on status: other threads can report failure here.
    cuDoubleComplex value=make_cuDoubleComplex(0.,0.);
    if(preconditioner.precision==PreconditionerPrecision::FP32) {
        cuFloatComplex sum=make_cuFloatComplex(0.f,0.f);
        for(std::int32_t p=preconditioner.rowptr[row];p<preconditioner.rowptr[row+1];++p) {
            auto x=input[preconditioner.columns[p]];
            const auto xf=make_cuFloatComplex(float(cuCreal(x)),float(cuCimag(x)));
            sum=cuCaddf(sum,cuCmulf(preconditioner.values32[p],xf));
        }
        value=make_cuDoubleComplex(double(cuCrealf(sum)),double(cuCimagf(sum)));
    }else{
        for(std::int32_t p=preconditioner.rowptr[row];p<preconditioner.rowptr[row+1];++p)
            value=cuCadd(value,cuCmul(preconditioner.values64[p],input[preconditioner.columns[p]]));
    }
    if(!isfinite(cuCreal(value))||!isfinite(cuCimag(value)))atomicExch(&scalar->status,-1);
    output[row]=value;
}
} // namespace radia_ih_cuda
