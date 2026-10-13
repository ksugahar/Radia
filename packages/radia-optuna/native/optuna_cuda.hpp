#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>
#include <memory>

// Host-owned RNG inputs and optional device-resident estimator construction.
namespace optuna_cuda {
class PreparedModel;
struct ModelDimension {
    double low, high, step, minimum;
    std::size_t choices;
    bool categorical, log_scale;
};
struct ModelConstruction {
    std::shared_ptr<PreparedModel> device;
    std::vector<std::vector<double>> mu, sigma, probabilities;
};
ModelConstruction build_model(const std::vector<double>& observations,
    const std::vector<double>& weights, const std::vector<ModelDimension>& dimensions,
    double prior_weight, bool endpoints, bool host_sampling);
struct Estimator {
    const std::vector<double>* weights;
    const std::vector<double>* mu;
    const std::vector<double>* sigma;
    const std::vector<double>* probabilities;
    std::size_t choices;
    double low, high, step;
    bool categorical, log_scale;
    const PreparedModel* prepared = nullptr;
    const std::vector<double>* denominator = nullptr;
};
struct Statistics {
    std::uint64_t density_launches = 0, selection_launches = 0;
    std::uint64_t host_to_device_bytes = 0, device_to_host_bytes = 0;
    std::size_t peak_temporary_bytes = 0;
    std::uint64_t cpu_density_calls = 0, cpu_selection_calls = 0;
    std::uint64_t cpu_parallel_density_calls = 0;
    int cpu_max_density_threads = 0;
    std::uint64_t device_allocation_calls = 0, workspace_reuses = 0;
    std::uint64_t avoided_upload_bytes = 0;
    std::size_t cached_device_bytes = 0;
    std::uint64_t model_build_launches = 0, model_reuses = 0;
    std::uint64_t preparation_launches = 0, model_download_bytes = 0;
    std::uint64_t cpu_preparation_calls = 0, cpu_cached_normalizer_values = 0;
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
void release_cached_memory() noexcept;
}
