#include "mex.h"
#include "gpu/mxGPUArray.h"
#include "optuna_runtime.hpp"
#include <stdexcept>

namespace optuna_cuda {
namespace {
void initialize() { mxInitGPU(); }
void* allocate(std::size_t bytes, void** owner) {
    const mwSize size = bytes;
    auto* array = mxGPUCreateGPUArray(1, &size, mxUINT8_CLASS, mxREAL,
                                     MX_GPU_DO_NOT_INITIALIZE);
    void* data = mxGPUGetData(array);
    *owner = mxGPUCreateMxArrayOnGPU(array);
    mexMakeArrayPersistent(static_cast<mxArray*>(*owner));
    mxGPUDestroyGPUArray(array);
    return data;
}
bool valid(void* owner) noexcept {
    return owner && mxGPUIsValidGPUData(static_cast<mxArray*>(owner));
}
void release(void* owner) noexcept { mxDestroyArray(static_cast<mxArray*>(owner)); }
}
void install_matlab_runtime() {
    install_runtime({initialize, allocate, valid, release});
}
}
