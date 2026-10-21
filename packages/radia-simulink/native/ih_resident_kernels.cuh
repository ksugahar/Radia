// Resident IH algebra primitives. Not a complete solver or release artifact.
// NGSolve/Radia retain ownership of operator assembly and physical contracts.
#pragma once
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cstddef>
#include <cstdint>

namespace radia_ih_cuda {
// Checked CSR buffers and state are owned by the future mexcuda context.
// Complex values contain the complete PEEC/BEM/ESIM coupled operator; no
// diagonal, sampled, real-only or frozen-mode substitution is made here.
struct ComplexCSRView {
    std::size_t rows;
    const std::int32_t* rowptr;
    const std::int32_t* columns;
    const cuDoubleComplex* values;
};

__global__ void complexMatvec(ComplexCSRView matrix,
        const cuDoubleComplex* state, cuDoubleComplex* result) {
    const std::size_t row=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(row>=matrix.rows) return;
    cuDoubleComplex value=make_cuDoubleComplex(0.,0.);
    for(std::int32_t p=matrix.rowptr[row];p<matrix.rowptr[row+1];++p)
        value=cuCadd(value,cuCmul(matrix.values[p],state[matrix.columns[p]]));
    result[row]=value;
}

// Q = sum_c weight_qc * |sum_j P_qcj*x_j|^2.
// P contains the complete physical loss projection. Squaring follows the
// coherent field/current sum, preserving every basis cross term. Weights
// include the peak-phasor factor 1/2, constitutive coefficient and declared
// integration/projection measure; units are owned by the checked assembler.
// Volume and surface loss must use distinct, explicitly declared measures.
__global__ void quadraticLoss(ComplexCSRView projection,
        const cuDoubleComplex* state, const double* weights,
        std::size_t samples, unsigned components, double* loss) {
    const std::size_t q=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(q>=samples) return;
    double total=0.;
    for(unsigned component=0;component<components;++component) {
        const std::size_t row=q*components+component;
        cuDoubleComplex field=make_cuDoubleComplex(0.,0.);
        for(std::int32_t p=projection.rowptr[row];p<projection.rowptr[row+1];++p)
            field=cuCadd(field,cuCmul(projection.values[p],state[projection.columns[p]]));
        total+=weights[row]*(cuCreal(field)*cuCreal(field)+cuCimag(field)*cuCimag(field));
    }
    loss[q]=total;
}

// Simulink uploads drive scalars; vector RHS generation stays on device.
__global__ void scaleDrive(const cuDoubleComplex* rhsPerAmp,
        const double* currentAmp, std::size_t n, cuDoubleComplex* rhs) {
    const std::size_t i=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(i<n) rhs[i]=make_cuDoubleComplex(cuCreal(rhsPerAmp[i])*currentAmp[0],
                                       cuCimag(rhsPerAmp[i])*currentAmp[0]);
}
// Voltage observation is supplied independently by the physical PEEC/BEM
// assembler, never reconstructed from the heat weights being checked.
__global__ void certifyIHPower(const cuDoubleComplex* voltageProjection,
        const cuDoubleComplex* state,const double* load,const double* constant,
        std::size_t emSize,std::size_t thermalSize,double current,double coilResistance,
        double relativeLimit,double floor,double* report,int* invalid) {
    if(blockIdx.x||threadIdx.x)return;
    double voltage=coilResistance*current,heat=0.;
    for(std::size_t j=0;j<emSize;++j)voltage+=cuCreal(cuCmul(voltageProjection[j],state[j]));
    for(std::size_t i=0;i<thermalSize;++i)heat+=constant[i]*load[i];
    const double input=.5*current*voltage,coil=.5*coilResistance*current*current;
    const double defect=input-coil-heat,scale=fmax(fmax(fabs(input),fabs(coil)+fabs(heat)),floor);
    report[0]=input;report[1]=coil;report[2]=heat;report[3]=defect;report[4]=fabs(defect)/scale;
    if(!isfinite(report[4])||report[4]>relativeLimit)atomicExch(invalid,1);
}
} // namespace radia_ih_cuda
