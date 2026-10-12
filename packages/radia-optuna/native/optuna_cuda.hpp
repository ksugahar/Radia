#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

// Host-owned estimator views. CUDA owns only per-call temporary buffers.
namespace optuna_cuda {
struct Estimator {
    const std::vector<double>* weights;
    const std::vector<double>* mu;
    const std::vector<double>* sigma;
    const std::vector<double>* probabilities;
    std::size_t choices;
    double low, high, step;
    bool categorical, log_scale;
};
struct Statistics {
    std::uint64_t density_launches = 0, selection_launches = 0;
    std::uint64_t host_to_device_bytes = 0, device_to_host_bytes = 0;
    std::size_t peak_temporary_bytes = 0;
    std::uint64_t cpu_density_calls = 0, cpu_selection_calls = 0;
    int runtime_version = 0, driver_version = 0, device = -1;
    bool stream_ordered_allocator = false;
};
std::vector<double> log_pdf(const std::vector<double>& samples,
                            const std::vector<Estimator>& estimators);
std::size_t best(const std::vector<double>& candidates,
                 const std::vector<Estimator>& below,
                 const std::vector<Estimator>& above);
bool use_gpu();
void set_gpu(bool enabled);
Statistics statistics();
}
