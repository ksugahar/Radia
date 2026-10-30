// Device instantiation of the production scalar RadFld formulas.
#pragma once
#include <array>
#include <cmath>
#include <limits>
#include "rad_field_cloud_kernel.h"
#include "rad_hdiv_vim.h"
#define RAD_FIELD_NAMESPACE rad_hdiv_cuda
#define RAD_FIELD_KERNEL __device__
#define RAD_FIELD_DEVICE_KERNELS
#define RAD_FIELD_POLY_CAP 6
#include "rad_hdiv_field_math.inc"
#include "rad_field_source_kernel.inc"
#undef RAD_FIELD_POLY_CAP
#undef RAD_FIELD_DEVICE_KERNELS
#undef RAD_FIELD_KERNEL
#undef RAD_FIELD_NAMESPACE

inline cudaError_t UploadFieldFarGauss() {
    rad_hdiv_cuda::GaussLegendreTable table{};
    rad_hdiv::CopyFieldFarGauss(&table.x[0][0],&table.w[0][0]);
    return cudaMemcpyToSymbol(rad_hdiv_cuda::DeviceFarGauss,&table,sizeof(table));
}
