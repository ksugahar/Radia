#pragma once
#if defined(RADIA_BLAS_MKL) && !defined(RADIA_BLAS_NVPL)
#include "fftw_mkl.h"
#elif defined(RADIA_BLAS_NVPL) && !defined(RADIA_BLAS_MKL)
#include "fftw_nvpl.h"
#else
#error "Select exactly one FFT provider with the BLAS provider"
#endif
