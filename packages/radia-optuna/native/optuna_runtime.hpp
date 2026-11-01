#pragma once
#include <cstddef>
#include <cstdint>

namespace optuna_cuda {
// Module-local adapter functions. No MATLAB, Python or solver types cross here.
// valid must identify allocation/context lifetime, never just pointer attributes.
struct RuntimeHooks {
    void (*initialize)();
    void* (*allocate)(std::size_t bytes, void** owner);
    bool (*valid)(void* owner) noexcept;
    void (*release)(void* owner) noexcept;
};
void install_runtime(RuntimeHooks hooks);
void shutdown_runtime(); // rejects live buffers before discarding callbacks
RuntimeHooks cuda_runtime_hooks();
void install_matlab_runtime();
}
