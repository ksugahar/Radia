// Fixed right preconditioners only. Outer vectors/operator/residual stay FP64.
#pragma once
#include "ih_resident_kernels.cuh"
#include "ih_bicgstab_kernels.cuh"
namespace radia_ih_cuda {
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
