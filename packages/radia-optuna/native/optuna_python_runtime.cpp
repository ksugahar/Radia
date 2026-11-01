#include "optuna_runtime.hpp"
#include <cuda.h>
#include <cuda_runtime.h>
#include <windows.h>
#include <memory>
#include <stdexcept>

namespace optuna_cuda {
namespace {
struct Allocation { void* data = nullptr; unsigned long long context_id = 0; };
using GetContext = CUresult (CUDAAPI*)(CUcontext*);
using GetId = CUresult (CUDAAPI*)(CUcontext,unsigned long long*);
HMODULE driver = nullptr; // retained for process lifetime, only explicit GPU loads it
GetContext get_context = nullptr;
GetId get_id = nullptr;
unsigned long long context_id() noexcept {
    CUcontext context = nullptr;
    unsigned long long id = 0;
    if (!get_context || !get_id || get_context(&context) != CUDA_SUCCESS || !context ||
        get_id(context, &id) != CUDA_SUCCESS) return 0;
    return id;
}
void initialize() {
    if (!driver) {
        driver=LoadLibraryExW(L"nvcuda.dll",nullptr,LOAD_LIBRARY_SEARCH_SYSTEM32);
        if (!driver) throw std::runtime_error("Optuna GPU requires the NVIDIA system driver");
        get_context=reinterpret_cast<GetContext>(GetProcAddress(driver,"cuCtxGetCurrent"));
        get_id=reinterpret_cast<GetId>(GetProcAddress(driver,"cuCtxGetId"));
    }
    if (!get_context || !get_id)
        throw std::runtime_error("Optuna GPU driver lacks context lifetime identity API");
    if (cudaFree(nullptr) != cudaSuccess || !context_id())
        throw std::runtime_error("Optuna explicit GPU context initialization failed");
}
bool valid(void* owner) noexcept {
    const auto* allocation = static_cast<Allocation*>(owner);
    return allocation && allocation->context_id && allocation->context_id == context_id();
}
void* allocate(std::size_t bytes, void** owner) {
    auto allocation = std::make_unique<Allocation>();
    allocation->context_id = context_id();
    if (!allocation->context_id || cudaMalloc(&allocation->data, bytes) != cudaSuccess)
        throw std::runtime_error("Optuna CUDA workspace allocation failed");
    void* result = allocation->data;
    *owner = allocation.release();
    return result;
}
void release(void* owner) noexcept {
    auto* allocation = static_cast<Allocation*>(owner);
    // Device reset already freed the allocation; never free a reused address.
    if (valid(owner)) cudaFree(allocation->data);
    delete allocation;
}
}
RuntimeHooks cuda_runtime_hooks() { return {initialize, allocate, valid, release}; }
}
