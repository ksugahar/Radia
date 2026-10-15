// Research material law: explicit piecewise-linear interpolation in kelvin.
// NGSolve temperature coefficients are evaluated with an owner-supplied map.
#pragma once
#include "ih_esim_cell.cuh"
#include "ih_thermal_kernels.cuh"
namespace radia_ih_cuda {
struct TemperatureMaterialTable {
    const double *kelvin,*resistivity,*saturation,*knots,*coefficients;
    const cuDoubleComplex* initialMu;
    int temperatures,knotsCount;
    // Properties are [element][temperature]; cubics add [interval][4].
};
__global__ void validateElementTemperatureRange(RealCSRView map,const double* temperature,
        const double* axis,int count,int* invalid){
    const auto e=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(e>=map.rows)return;
    double t=0.;for(int k=map.rowptr[e];k<map.rowptr[e+1];++k)t+=map.values[k]*temperature[map.columns[k]];
    if(!isfinite(t)||t<axis[0]||t>axis[count-1])atomicExch(invalid,1);
}
__global__ void updateElementESIMMaterials(RealCSRView temperatureEvaluation,
        const double* acceptedTemperature,TemperatureMaterialTable table,
        const ESIMCellParameters* geometry,ESIMCellParameters* candidate,
        double* localTemperature,double* candidateCoefficients,int* invalid) {
    const auto e=std::size_t(blockIdx.x)*blockDim.x+threadIdx.x;
    if(e>=temperatureEvaluation.rows)return;
    double t=0.;for(int k=temperatureEvaluation.rowptr[e];k<temperatureEvaluation.rowptr[e+1];++k)
        t+=temperatureEvaluation.values[k]*acceptedTemperature[temperatureEvaluation.columns[k]];
    localTemperature[e]=t;
    if(!isfinite(t)||t<=0.||table.temperatures<2||
       t<table.kelvin[0]||t>table.kelvin[table.temperatures-1]){atomicExch(invalid,1);return;}
    int lo=0,hi=table.temperatures-1;
    while(hi-lo>1){const int m=(lo+hi)/2;if(t<table.kelvin[m])hi=m;else lo=m;}
    const double delta=table.kelvin[hi]-table.kelvin[lo];
    if(!isfinite(delta)||delta<=0.){atomicExch(invalid,1);return;}
    const double f=(t-table.kelvin[lo])/delta;
    const auto a=e*std::size_t(table.temperatures)+lo,b=a+(hi-lo);
    auto p=geometry[e];
    p.resistivity=(1.-f)*table.resistivity[a]+f*table.resistivity[b];
    p.law.initialMu=cuCadd(cuCmul(make_cuDoubleComplex(1.-f,0.),table.initialMu[a]),cuCmul(make_cuDoubleComplex(f,0.),table.initialMu[b]));
    p.law.saturationMagnetization=(1.-f)*table.saturation[a]+f*table.saturation[b];
    p.law.count=table.knotsCount;p.law.knots=table.knots;
    if(table.knotsCount){
        const auto n=std::size_t(table.knotsCount-1)*4;
        auto* coefficients=candidateCoefficients+e*n;
        for(std::size_t k=0;k<n;++k){const double c=(1.-f)*table.coefficients[a*n+k]+f*table.coefficients[b*n+k];
            coefficients[k]=c;if(!isfinite(c))atomicExch(invalid,1);}
        p.law.coefficients=coefficients;
    }else p.law.coefficients=nullptr;
    if(!isfinite(p.resistivity)||p.resistivity<=0.||!cellFinite(p.law.initialMu)||
       cuCreal(p.law.initialMu)<=0.||cuCimag(p.law.initialMu)>0.||
       !isfinite(p.law.saturationMagnetization)||p.law.saturationMagnetization<0.)atomicExch(invalid,1);
    candidate[e]=p;
}
} // namespace radia_ih_cuda
