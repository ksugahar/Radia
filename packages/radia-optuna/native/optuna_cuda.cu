#include "optuna_cuda.hpp"
#include "optuna_runtime.hpp"
#include <cuda_runtime.h>
#include <cub/device/device_segmented_radix_sort.cuh>
#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <map>
#include <omp.h>

namespace optuna_cuda {
namespace {
Statistics stats;
bool gpu_enabled = true;
RuntimeHooks runtime{};
std::size_t active_buffers = 0;
void check(cudaError_t code, const char* operation) {
    if (code != cudaSuccess)
        throw std::runtime_error(std::string("Optuna CUDA ") + operation +
                                 ": " + cudaGetErrorString(code));
}
std::size_t product(std::size_t a, std::size_t b) {
    if (b && a > std::numeric_limits<std::size_t>::max() / b)
        throw std::overflow_error("Optuna CUDA allocation size overflow");
    return a * b;
}
std::size_t sum(std::size_t a, std::size_t b) {
    if (a > std::numeric_limits<std::size_t>::max() - b)
        throw std::overflow_error("Optuna CUDA arena size overflow");
    return a + b;
}
// Pack input into one range. Device workspaces retain an exact host shadow
// so unchanged bytes can remain resident without relying on object identity.
struct Arena {
    std::vector<unsigned char> input;
    std::size_t bytes = 0;
    template<class T> std::size_t reserve(std::size_t count) {
        const auto padding = (alignof(T) - bytes % alignof(T)) % alignof(T);
        const auto offset = sum(bytes, padding);
        bytes = sum(offset, product(count, sizeof(T)));
        return offset;
    }
    template<class T> std::size_t append(const std::vector<T>& values) {
        const auto offset = reserve<T>(values.size());
        input.resize(bytes);
        if (!values.empty()) std::memcpy(input.data() + offset, values.data(), product(values.size(), sizeof(T)));
        return offset;
    }
};
struct Workspace {
    void* data = nullptr;
    void* owner = nullptr;
    std::size_t capacity = 0;
    std::vector<unsigned char> input;
    bool model_ready = false;
    ~Workspace() { if (owner) runtime.release(owner); }
};
std::map<std::pair<int, int>, std::shared_ptr<Workspace>> workspaces;
constexpr std::size_t max_cached_bytes = 64 * 1024 * 1024;
void update_cached_bytes() {
    stats.cached_device_bytes = 0;
    for (const auto& item : workspaces) stats.cached_device_bytes += item.second->capacity;
}
struct Buffer {
    void* data = nullptr;
    void* owner = nullptr;
    std::shared_ptr<Workspace> cached;
    explicit Buffer(std::size_t bytes, int slot = 0) {
        std::size_t capacity = 4096;
        while (capacity < bytes && capacity < max_cached_bytes) capacity *= 2;
        std::size_t other_bytes = 0;
        for (const auto& item : workspaces)
            if (item.first.first == stats.device && item.first.second != slot)
                other_bytes += item.second->capacity;
        if (bytes && bytes <= max_cached_bytes && capacity <= max_cached_bytes - other_bytes) {
            auto& workspace = workspaces[{stats.device, slot}];
            if (workspace && workspace->owner && !runtime.valid(workspace->owner)) {
                // An invalid context epoch is not repaired by pointer reuse.
                workspace.reset();
                throw std::runtime_error("Optuna device allocation epoch expired");
            }
            // A live PreparedModel pins its outputs. A later build in the same
            // slot must not overwrite it, even if the address/capacity matches.
            if (workspace && workspace.use_count() > 1) workspace.reset();
            if (!workspace || workspace->capacity < bytes) {
                auto next = std::make_shared<Workspace>();
                next->data = runtime.allocate(capacity, &next->owner);
                if (!next->data || !next->owner)
                    throw std::runtime_error("Optuna runtime returned an invalid allocation");
                next->capacity = capacity;
                workspace = std::move(next);
                ++stats.device_allocation_calls;
            } else {
                ++stats.workspace_reuses;
            }
            cached = workspace;
            data = workspace->data;
            update_cached_bytes();
            ++active_buffers;
            return;
        }
        if (bytes) {
            data = runtime.allocate(bytes, &owner);
            if (!data || !owner) {
                if (owner) runtime.release(owner);
                throw std::runtime_error("Optuna runtime returned an invalid temporary allocation");
            }
        }
        if (bytes) ++stats.device_allocation_calls;
        ++active_buffers;
    }
    ~Buffer() {
        if (data && !cached) {
            runtime.release(owner);
        }
        --active_buffers;
    }
    Buffer(const Buffer&) = delete;
    Buffer& operator=(const Buffer&) = delete;
    template<class T> T* as(std::size_t offset = 0) {
        if (cached && !runtime.valid(cached->owner))
            throw std::runtime_error("Optuna prepared allocation epoch expired");
        if (owner && !runtime.valid(owner))
            throw std::runtime_error("Optuna temporary allocation epoch expired");
        return reinterpret_cast<T*>(static_cast<unsigned char*>(data) + offset);
    }
};
bool upload(Buffer& target, const void* source, std::size_t bytes) {
    const auto* input = static_cast<const unsigned char*>(source);
    std::size_t first = 0, last = bytes;
    std::vector<unsigned char> next;
    if (target.cached) {
        const auto& previous = target.cached->input;
        if (previous.size() == bytes) {
            while (first < bytes && previous[first] == input[first]) ++first;
            while (last > first && previous[last - 1] == input[last - 1]) --last;
        }
        next.assign(input, input + bytes);
    }
    try {
        if (last > first) check(cudaMemcpy(static_cast<unsigned char*>(target.data) + first,
            input + first, last - first, cudaMemcpyHostToDevice), "upload");
    } catch (...) {
        if (target.cached) target.cached->input.clear();
        throw;
    }
    stats.host_to_device_bytes += last - first;
    stats.avoided_upload_bytes += bytes - (last - first);
    if (target.cached) target.cached->input.swap(next);
    return last > first;
}
void download(void* target, const void* source, std::size_t bytes) {
    if (bytes) check(cudaMemcpy(target, source, bytes, cudaMemcpyDeviceToHost), "download");
    stats.device_to_host_bytes += bytes;
}
struct DeviceEstimator {
    std::size_t weights, mu, sigma, probabilities, kernels, choices;
    double low, high, step;
    bool categorical, log_scale;
    std::size_t denominator = static_cast<std::size_t>(-1);
    std::size_t log_sigma = static_cast<std::size_t>(-1);
};
struct Packed {
    std::vector<DeviceEstimator> estimators;
    std::vector<double> values;
    std::size_t normalization_count = 0;
    std::size_t append(const std::vector<double>* value) {
        const auto offset = values.size();
        if (value) values.insert(values.end(), value->begin(), value->end());
        return offset;
    }
    explicit Packed(const std::vector<Estimator>& source, bool prepare = false) {
        if (source.empty()) throw std::invalid_argument("Optuna CUDA requires estimators");
        const auto kernels = source.front().weights->size();
        estimators.resize(source.size());
        std::memset(estimators.data(), 0, product(estimators.size(), sizeof(DeviceEstimator)));
        for (std::size_t d = 0; d < source.size(); ++d) {
            const auto& e = source[d]; auto& descriptor = estimators[d];
            if (!kernels || e.weights->size() != kernels)
                throw std::invalid_argument("Optuna CUDA mixture sizes disagree");
            descriptor.weights = append(e.weights); descriptor.mu = append(e.mu);
            descriptor.sigma = append(e.sigma); descriptor.probabilities = append(e.probabilities);
            descriptor.kernels = kernels; descriptor.choices = e.choices;
            descriptor.low = e.low; descriptor.high = e.high; descriptor.step = e.step;
            descriptor.categorical = e.categorical; descriptor.log_scale = e.log_scale;
            descriptor.denominator = descriptor.log_sigma = static_cast<std::size_t>(-1);
            if (prepare && !e.categorical) {
                descriptor.denominator = normalization_count;
                normalization_count = sum(normalization_count, kernels);
                descriptor.log_sigma = normalization_count;
                normalization_count = sum(normalization_count, kernels);
            }
        }
    }
};
__host__ __device__ double log_mass(double lower, double upper) {
    constexpr double inv = 0.70710678118654752440;
    double mass;
    if (upper <= 0.0) mass = 0.5 * (erfc(-upper * inv) - erfc(-lower * inv));
    else if (lower >= 0.0) mass = 0.5 * (erfc(lower * inv) - erfc(upper * inv));
    else mass = 0.5 * (erf(upper * inv) - erf(lower * inv));
    return log(fmax(mass, DBL_MIN));
}
__host__ __device__ double component(const DeviceEstimator& e, const double* data,
                            double sample, std::size_t kernel, const double* normalization = nullptr) {
    if (e.categorical) {
        const auto choice = static_cast<std::size_t>(sample) - 1;
        return log(fmax(data[e.probabilities + kernel * e.choices + choice], DBL_MIN));
    }
    const double mu = data[e.mu + kernel], sigma = data[e.sigma + kernel];
    const bool prepared = e.denominator != static_cast<std::size_t>(-1);
    const double denominator = prepared ? normalization[e.denominator + kernel] :
        log_mass((e.low - mu) / sigma, (e.high - mu) / sigma);
    if (isfinite(e.step)) {
        double lo = sample - 0.5 * e.step, hi = sample + 0.5 * e.step;
        if (e.log_scale) { lo = log(lo); hi = log(hi); }
        return log_mass((lo - mu) / sigma, (hi - mu) / sigma) - denominator;
    }
    const double x = e.log_scale ? log(sample) : sample;
    const double z = (x - mu) / sigma;
    return -0.5 * z * z - 0.91893853320467274178 -
        (prepared ? normalization[e.log_sigma + kernel] : log(sigma)) - denominator;
}
__global__ void prepare_density(const DeviceEstimator* estimators, const double* data, double* normalization) {
    const auto& e = estimators[blockIdx.x];
    if (e.categorical) return;
    for (std::size_t k = threadIdx.x; k < e.kernels; k += blockDim.x) {
        const double mu = data[e.mu + k], sigma = data[e.sigma + k];
        normalization[e.denominator + k] = log_mass((e.low - mu) / sigma, (e.high - mu) / sigma);
        normalization[e.log_sigma + k] = log(sigma);
    }
}
struct SortValue { double key; std::size_t index; };
// One block owns a dimension. Read both lanes before writing so duplicate-key
// choices exactly follow the existing NumPy four-lane network (<=256 rows).
__device__ void network_stage(SortValue* values, std::size_t padded,
    int permutation, unsigned maximum_mask) {
    const auto i = static_cast<std::size_t>(threadIdx.x);
    SortValue selected{};
    if (i < padded) {
        const auto lane = i % 4;
        const auto other = i - lane + ((permutation >> (lane * 2)) & 3);
        const auto a = values[i], b = values[other];
        const bool maximum = (maximum_mask & (1U << lane)) != 0;
        selected = (maximum ? a.key >= b.key : a.key <= b.key) ? a : b;
    }
    __syncthreads();
    if (i < padded) values[i] = selected;
    __syncthreads();
}
__device__ void network_blocks(SortValue* values, std::size_t padded,
    std::size_t width, bool reverse) {
    const auto i = static_cast<std::size_t>(threadIdx.x);
    SortValue selected{};
    if (i < padded) {
        const auto block = i / 4, lane = i % 4;
        const auto region = block / width * width, offset = block % width;
        const auto other_block = reverse ? region + width - offset - 1 :
            region + (offset < width / 2 ? offset + width / 2 : offset - width / 2);
        const auto a = values[i], b = values[other_block * 4 + (reverse ? 3 - lane : lane)];
        const bool maximum = offset >= width / 2;
        selected = (maximum ? a.key >= b.key : a.key <= b.key) ? a : b;
    }
    __syncthreads();
    if (i < padded) values[i] = selected;
    __syncthreads();
}
__global__ void initialize_large_model(const double* observations, std::size_t n,
    std::size_t dimensions, std::size_t padded, const ModelDimension* settings,
    const DeviceEstimator* estimators, double* data, double* keys, std::size_t* indices) {
    const auto i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= dimensions * padded) return;
    const auto dim = i / padded, row = i % padded;
    const auto& e = estimators[dim];
    const double raw = row < n && !e.categorical ? observations[row * dimensions + dim] : INFINITY;
    const double mu = row < n && !e.categorical && settings[dim].log_scale ? log(raw) : raw;
    keys[i] = mu; indices[i] = row;
    if (row < n && !e.categorical) data[e.mu + row] = mu;
}
__global__ void construct_model(const double* observations, std::size_t n,
    std::size_t dimensions, const ModelDimension* settings,
    const DeviceEstimator* estimators, double* data, SortValue* scratch,
    std::size_t padded, const double* sorted_keys, const std::size_t* sorted_indices, double prior_weight, bool endpoints) {
    const auto dim = static_cast<std::size_t>(blockIdx.x);
    const auto& e = estimators[dim]; const auto& setting = settings[dim];
    if (e.categorical) {
        for (std::size_t row = threadIdx.x; row <= n; row += blockDim.x) {
            data[e.mu + row] = 0.0; data[e.sigma + row] = 0.0;
            const double initial = n == 0 ? 1.0 / static_cast<double>(e.choices) :
                prior_weight / static_cast<double>(n + 1);
            double total = 0.0;
            const auto observed = row < n ? static_cast<std::size_t>(observations[row * dimensions + dim] - 1) : e.choices;
            for (std::size_t choice = 0; choice < e.choices; ++choice) {
                const double value = initial + (n && choice == observed ? 1.0 : 0.0);
                data[e.probabilities + row * e.choices + choice] = value;
                total += value;
            }
            if (n) for (std::size_t choice = 0; choice < e.choices; ++choice)
                data[e.probabilities + row * e.choices + choice] = total > 0.0 ?
                    data[e.probabilities + row * e.choices + choice] / total :
                    1.0 / static_cast<double>(e.choices);
        }
        return;
    }
    auto* sorted = scratch + dim * padded;
    if (n <= 256) for (std::size_t i = threadIdx.x; i < padded; i += blockDim.x) {
        const double raw = i < n ? observations[i * dimensions + dim] : INFINITY;
        const double mu = i < n && setting.log_scale ? log(raw) : raw;
        sorted[i] = {mu, i < n ? i : static_cast<std::size_t>(-1)};
        if (i < n) data[e.mu + i] = mu;
    }
    __syncthreads();
    if (n > 1 && n <= 256) {
        network_stage(sorted, padded, 0xB1, 0xA); // 1,0,3,2
        network_stage(sorted, padded, 0x1B, 0xC); // 3,2,1,0
        network_stage(sorted, padded, 0xB1, 0xA);
        for (std::size_t width = 2; width <= padded / 4; width *= 2) {
            network_blocks(sorted, padded, width, true);
            for (std::size_t inner = width / 2; inner >= 2; inner /= 2)
                network_blocks(sorted, padded, inner, false);
            network_stage(sorted, padded, 0x4E, 0xC); // 2,3,0,1
            network_stage(sorted, padded, 0xB1, 0xA);
        }
    }
    for (std::size_t i = threadIdx.x; i < n; i += blockDim.x) {
        const auto position = dim * padded + i;
        const double value = n > 256 ? sorted_keys[position] : sorted[i].key;
        const double left = i == 0 ? e.low : (n > 256 ? sorted_keys[position - 1] : sorted[i - 1].key);
        const double right = i + 1 == n ? e.high : (n > 256 ? sorted_keys[position + 1] : sorted[i + 1].key);
        double sigma = fmax(value - left, right - value);
        if (!endpoints && n >= 2) {
            if (i == 0) sigma = (n > 256 ? sorted_keys[position + 1] : sorted[1].key) - value;
            if (i + 1 == n) sigma = value - (n > 256 ? sorted_keys[position - 1] : sorted[i - 1].key);
        }
        data[e.sigma + (n > 256 ? sorted_indices[position] : sorted[i].index)] = fmin(fmax(sigma, setting.minimum), e.high - e.low);
    }
    if (threadIdx.x == 0) {
        data[e.mu + n] = 0.5 * (e.low + e.high);
        data[e.sigma + n] = e.high - e.low;
    }
    __syncthreads();
    for (std::size_t k = threadIdx.x; k <= n; k += blockDim.x) {
        const double mu = data[e.mu + k], sigma = data[e.sigma + k];
        data[e.denominator + k] = log_mass((e.low - mu) / sigma, (e.high - mu) / sigma);
        data[e.log_sigma + k] = log(sigma);
    }
}
// A block per candidate evaluates mixture components in parallel. The final
// ordered reduction retains the oracle's summation and tie-breaking contract.
__global__ void density(const double* samples, std::size_t count,
                         const DeviceEstimator* estimators, std::size_t dimensions,
                         const double* data, const double* normalization, double* terms, double* output) {
    const std::size_t row = blockIdx.x;
    if (row >= count) return;
    const auto kernels = estimators[0].kernels;
    double* components = terms + row * kernels;
    for (std::size_t k = threadIdx.x; k < kernels; k += blockDim.x) {
        double value = log(data[estimators[0].weights + k]);
        for (std::size_t d = 0; d < dimensions; ++d)
            value += component(estimators[d], data, samples[row * dimensions + d], k, normalization);
        components[k] = value;
    }
    __syncthreads();
    __shared__ double maximum;
    if (threadIdx.x == 0) {
        maximum = -INFINITY;
        for (std::size_t k = 0; k < kernels; ++k) maximum = fmax(maximum, components[k]);
    }
    if (kernels < 128) {
        if (threadIdx.x != 0) return;
        if (maximum == -INFINITY) { output[row] = maximum; return; }
        double total = 0.0;
        for (std::size_t k = 0; k < kernels; ++k) total += exp(components[k] - maximum);
        output[row] = maximum + log(total);
        return;
    }
    __syncthreads();
    if (maximum == -INFINITY) {
        if (threadIdx.x == 0) output[row] = maximum;
        return;
    }
    // Independent exponentials can run in parallel. Keep both maximum and
    // final addition ordered, including the first-maximum candidate contract.
    for (std::size_t k = threadIdx.x; k < kernels; k += blockDim.x)
        components[k] = exp(components[k] - maximum);
    __syncthreads();
    if (threadIdx.x != 0) return;
    double total = 0.0;
    for (std::size_t k = 0; k < kernels; ++k) total += components[k];
    output[row] = maximum + log(total);
}
__global__ void select_best(const double* below, const double* above,
                            std::size_t count, std::size_t* result) {
    std::size_t best = 0;
    double maximum = -INFINITY;
    for (std::size_t row = 0; row < count; ++row) {
        const double acquisition = below[row] - above[row];
        if (!isfinite(acquisition)) { *result = count; return; }
        if (acquisition > maximum) { maximum = acquisition; best = row; }
    }
    *result = best;
}
// Share component equations with CUDA; parallelize candidates only. Ordered
// per-candidate summation and first-maximum selection preserve oracle behavior.
std::vector<double> host_pdf(const std::vector<double>& samples,
                             const std::vector<Estimator>& estimators) {
    Packed packed(estimators, true);
    std::vector<double> normalization(packed.normalization_count);
    for (std::size_t d = 0; d < estimators.size(); ++d) {
        const auto& e = packed.estimators[d];
        if (e.categorical) continue;
        const auto* cached = estimators[d].denominator;
        const bool reuse = cached && cached->size() == e.kernels;
        for (std::size_t k = 0; k < e.kernels; ++k) {
            const double mu = packed.values[e.mu + k], sigma = packed.values[e.sigma + k];
            normalization[e.denominator + k] = reuse ? (*cached)[k] :
                log_mass((e.low - mu) / sigma, (e.high - mu) / sigma);
            normalization[e.log_sigma + k] = std::log(sigma);
        }
        if (reuse) stats.cpu_cached_normalizer_values += e.kernels;
    }
    ++stats.cpu_preparation_calls;
    if (samples.size() % estimators.size())
        throw std::invalid_argument("Optuna CPU sample shape");
    const auto count = samples.size() / estimators.size();
    const auto kernels = packed.estimators.front().kernels;
    if (count > static_cast<std::size_t>(LLONG_MAX))
        throw std::overflow_error("Optuna CPU loop size overflow");
    const auto component_count = product(count, kernels);
    const auto work = product(component_count, estimators.size());
    std::vector<double> terms(component_count), result(count);
    // Default TPE has only 24 candidates but can have thousands of components.
    // Keep short densities serial and bound the team by independent rows.
    const int threads = static_cast<int>(std::min<std::size_t>(count, omp_get_max_threads()));
    const bool parallel = threads > 1 && (count >= 256 || work >= 16384);
    int used_threads = 1;
    #pragma omp parallel if(parallel) num_threads(threads > 0 ? threads : 1)
    {
    #pragma omp single
    used_threads = omp_get_num_threads();
    #pragma omp for schedule(static)
    for (long long row = 0; row < static_cast<long long>(count); ++row) {
        double* values = terms.data() + static_cast<std::size_t>(row) * kernels;
        double maximum = -INFINITY;
        for (std::size_t k = 0; k < kernels; ++k) {
            double value = std::log(packed.values[packed.estimators[0].weights + k]);
            for (std::size_t d = 0; d < estimators.size(); ++d)
                value += component(packed.estimators[d], packed.values.data(),
                                   samples[row * estimators.size() + d], k, normalization.data());
            values[k] = value;
            maximum = std::fmax(maximum, value);
        }
        double total = 0.0;
        for (std::size_t k = 0; k < kernels; ++k)
            total += std::exp(values[k] - maximum);
        result[row] = maximum + std::log(total);
    }
    }
    for (double value : result)
        if (!std::isfinite(value)) throw std::runtime_error("Optuna CPU density is nonfinite");
    ++stats.cpu_density_calls;
    if (used_threads > 1) ++stats.cpu_parallel_density_calls;
    stats.cpu_max_density_threads = std::max(stats.cpu_max_density_threads, used_threads);
    return result;
}
void initialize() {
    if (!runtime.initialize) throw std::logic_error("Optuna runtime adapter not installed");
    runtime.initialize();
    check(cudaGetDevice(&stats.device), "get device");
    check(cudaRuntimeGetVersion(&stats.runtime_version), "runtime version");
    check(cudaDriverGetVersion(&stats.driver_version), "driver version");
    int pool_supported = 0;
    check(cudaDeviceGetAttribute(&pool_supported, cudaDevAttrMemoryPoolsSupported,
                                 stats.device), "memory pool support");
    stats.stream_ordered_allocator = false; // allocator ownership belongs to adapter
}
struct DeviceMixture {
    std::size_t descriptors, data, normalization = 0;
    DeviceMixture(Arena& arena, const Packed& p)
        : descriptors(arena.append(p.estimators)), data(arena.append(p.values)) {}
};
void launch(Buffer& arena, std::size_t samples, std::size_t count, std::size_t dimensions,
            const DeviceMixture& mixture, std::size_t terms, std::size_t output) {
    if (count > static_cast<std::size_t>(INT_MAX))
        throw std::overflow_error("Optuna CUDA grid size overflow");
    density<<<static_cast<unsigned>(count), 128>>>(
        arena.as<double>(samples), count, arena.as<DeviceEstimator>(mixture.descriptors),
        dimensions, arena.as<double>(mixture.data), arena.as<double>(mixture.normalization), arena.as<double>(terms), arena.as<double>(output));
    check(cudaGetLastError(), "launch density");
    ++stats.density_launches;
}
}
class PreparedModel {
public:
    std::unique_ptr<Buffer> storage;
    std::size_t descriptors = 0, dimensions = 0, kernels = 0;
    mutable bool completed = false;
    ~PreparedModel() {
        // A launch, allocation, download or scoring failure must never leave
        // a matching input shadow pointing at incomplete model outputs.
        if (!completed && storage && storage->cached) {
            storage->cached->input.clear(); storage->cached->model_ready = false;
        }
    }
    void accept() const noexcept {
        completed = true;
        if (storage->cached) storage->cached->model_ready = true;
    }
};
ModelConstruction build_model(const std::vector<double>& observations,
    const std::vector<double>& weights, const std::vector<ModelDimension>& dimensions,
    double prior_weight, bool endpoints, bool host_sampling) {
    if (!gpu_enabled) throw std::logic_error("GPU model construction requires explicit GPU mode");
    initialize();
    const auto dims = dimensions.size();
    if (!dims || dims > static_cast<std::size_t>(INT_MAX) || observations.size() % dims ||
        weights.size() != observations.size() / dims + 1)
        throw std::invalid_argument("Optuna CUDA model shape");
    const auto n = observations.size() / dims, kernels = n + 1;
    std::size_t padded = 4;
    while (padded < n) padded = product(padded, 2);
    Arena layout;
    const auto input = layout.append(observations);
    const auto weight_offset = layout.append(weights);
    const auto settings_offset = layout.append(dimensions);
    layout.append(std::vector<double>{prior_weight, endpoints ? 1.0 : 0.0});
    const auto items = product(dims, padded);
    if (n > 256 && items > static_cast<std::size_t>(INT_MAX))
        throw std::overflow_error("Optuna CUDA segmented sort size overflow");
    std::vector<int> offsets;
    if (n > 256) for (std::size_t d = 0; d <= dims; ++d) offsets.push_back(static_cast<int>(d * padded));
    const auto offsets_offset = layout.append(offsets);
    const auto descriptors_offset = layout.reserve<DeviceEstimator>(dims);
    const auto input_bytes = layout.bytes;
    const auto model_count = product(dims, kernels);
    const auto mu_offset = layout.reserve<double>(model_count);
    const auto sigma_offset = layout.reserve<double>(model_count);
    std::size_t probability_count = 0;
    for (const auto& d : dimensions) if (d.categorical)
        probability_count = sum(probability_count, product(kernels, d.choices));
    const auto probabilities_offset = layout.reserve<double>(probability_count);
    const auto denominator_offset = layout.reserve<double>(model_count);
    const auto log_sigma_offset = layout.reserve<double>(model_count);
    const auto scratch_offset = layout.reserve<SortValue>(n <= 256 ? items : 0);
    const auto keys_in = layout.reserve<double>(n > 256 ? items : 0);
    const auto keys_out = layout.reserve<double>(n > 256 ? items : 0);
    const auto indices_in = layout.reserve<std::size_t>(n > 256 ? items : 0);
    const auto indices_out = layout.reserve<std::size_t>(n > 256 ? items : 0);
    std::vector<DeviceEstimator> descriptors(dims);
    std::memset(descriptors.data(), 0, product(dims, sizeof(DeviceEstimator)));
    std::size_t probability_position = 0;
    for (std::size_t d = 0; d < dims; ++d) {
        auto& e = descriptors[d]; const auto& setting = dimensions[d];
        e.weights = weight_offset / sizeof(double);
        e.mu = mu_offset / sizeof(double) + d * kernels;
        e.sigma = sigma_offset / sizeof(double) + d * kernels;
        e.probabilities = probabilities_offset / sizeof(double) + probability_position;
        e.denominator = denominator_offset / sizeof(double) + d * kernels;
        e.log_sigma = log_sigma_offset / sizeof(double) + d * kernels;
        e.kernels = kernels; e.choices = setting.choices;
        e.low = setting.low; e.high = setting.high;
        e.step = setting.step;
        e.categorical = setting.categorical; e.log_scale = setting.log_scale;
        if (setting.categorical) probability_position += kernels * setting.choices;
    }
    layout.input.resize(input_bytes);
    std::memcpy(layout.input.data() + descriptors_offset, descriptors.data(), product(dims, sizeof(DeviceEstimator)));
    ModelConstruction result;
    result.device = std::make_shared<PreparedModel>();
    auto& model = *result.device;
    model.storage = std::make_unique<Buffer>(layout.bytes, host_sampling ? 1 : 2);
    model.descriptors = descriptors_offset; model.dimensions = dims; model.kernels = kernels;
    auto& arena = *model.storage;
    stats.peak_temporary_bytes = std::max(stats.peak_temporary_bytes, layout.bytes);
    const bool changed = upload(arena, layout.input.data(), layout.input.size());
    if (changed || !arena.cached || !arena.cached->model_ready) {
        // Keep readiness in the workspace too: MATLAB allocation errors can
        // leave the C++ stack without running ordinary exception destructors.
        if (arena.cached) arena.cached->model_ready = false;
        if (n > 256) {
            initialize_large_model<<<static_cast<unsigned>((items + 255) / 256), 256>>>(
                arena.as<double>(input), n, dims, padded, arena.as<ModelDimension>(settings_offset),
                arena.as<DeviceEstimator>(descriptors_offset), arena.as<double>(),
                arena.as<double>(keys_in), arena.as<std::size_t>(indices_in));
            check(cudaGetLastError(), "initialize model sort"); ++stats.model_build_launches;
            std::size_t bytes = 0;
            check(cub::DeviceSegmentedRadixSort::SortPairs(nullptr, bytes,
                arena.as<double>(keys_in), arena.as<double>(keys_out),
                arena.as<std::size_t>(indices_in), arena.as<std::size_t>(indices_out),
                static_cast<int>(items), static_cast<int>(dims), arena.as<int>(offsets_offset),
                arena.as<int>(offsets_offset) + 1), "query segmented model sort");
            stats.peak_sort_temporary_bytes = std::max(stats.peak_sort_temporary_bytes, bytes);
            Buffer temporary(std::max<std::size_t>(bytes, 1), 3);
            check(cub::DeviceSegmentedRadixSort::SortPairs(temporary.data, bytes,
                arena.as<double>(keys_in), arena.as<double>(keys_out),
                arena.as<std::size_t>(indices_in), arena.as<std::size_t>(indices_out),
                static_cast<int>(items), static_cast<int>(dims), arena.as<int>(offsets_offset),
                arena.as<int>(offsets_offset) + 1), "sort model observations");
            ++stats.model_radix_sorts;
        }
        construct_model<<<static_cast<unsigned>(dims), 256>>>(arena.as<double>(input), n,
            dims, arena.as<ModelDimension>(settings_offset), arena.as<DeviceEstimator>(descriptors_offset),
            arena.as<double>(), arena.as<SortValue>(scratch_offset), padded, arena.as<double>(keys_out),
            arena.as<std::size_t>(indices_out), prior_weight, endpoints);
        check(cudaGetLastError(), "construct TPE model");
        ++stats.model_build_launches;
    } else ++stats.model_reuses;
    result.mu.resize(dims); result.sigma.resize(dims); result.probabilities.resize(dims);
    if (host_sampling) {
        // Only the good model returns to the ordered host RNG. The bad model
        // remains entirely on device through candidate scoring.
        const auto host_count = sum(product(model_count, 2), probability_count);
        std::vector<double> host(host_count);
        const auto bytes = product(host_count, sizeof(double));
        download(host.data(), arena.as<double>(mu_offset), bytes);
        stats.model_download_bytes += bytes;
        probability_position = 2 * model_count;
        for (std::size_t d = 0; d < dims; ++d) {
            if (dimensions[d].categorical) {
                const auto count = kernels * dimensions[d].choices;
                result.probabilities[d].assign(host.begin() + probability_position,
                    host.begin() + probability_position + count);
                probability_position += count;
            } else {
                result.mu[d].assign(host.begin() + d * kernels, host.begin() + (d + 1) * kernels);
                result.sigma[d].assign(host.begin() + model_count + d * kernels,
                    host.begin() + model_count + (d + 1) * kernels);
            }
        }
    }
    return result;
}
std::vector<double> log_pdf(const std::vector<double>& samples,
                            const std::vector<Estimator>& estimators) {
    if (!gpu_enabled) return host_pdf(samples, estimators);
    initialize();
    Packed packed(estimators, true);
    if (samples.size() % estimators.size()) throw std::invalid_argument("Optuna CUDA sample shape");
    const auto count = samples.size() / estimators.size();
    std::vector<double> result(count);
    if (!count) return result;
    const auto output_bytes = product(count, sizeof(double));
    Arena layout;
    const auto input = layout.append(samples);
    DeviceMixture mixture(layout, packed);
    mixture.normalization = layout.reserve<double>(packed.normalization_count);
    const auto output = layout.reserve<double>(count);
    const auto terms = layout.reserve<double>(product(count, packed.estimators.front().kernels));
    Buffer arena(layout.bytes);
    stats.peak_temporary_bytes = std::max(stats.peak_temporary_bytes, layout.bytes);
    upload(arena, layout.input.data(), layout.input.size());
    prepare_density<<<static_cast<unsigned>(estimators.size()), 128>>>(
        arena.as<DeviceEstimator>(mixture.descriptors), arena.as<double>(mixture.data), arena.as<double>(mixture.normalization));
    check(cudaGetLastError(), "prepare density"); ++stats.preparation_launches;
    launch(arena, input, count, estimators.size(), mixture, terms, output);
    download(result.data(), arena.as<double>(output), output_bytes);
    for (double value : result)
        if (!std::isfinite(value)) throw std::runtime_error("Optuna CUDA density is nonfinite");
    return result;
}
std::size_t best(const std::vector<double>& candidates,
                 const std::vector<Estimator>& below, const std::vector<Estimator>& above) {
    if (!gpu_enabled) {
        if (below.empty() || below.size() != above.size() || candidates.empty() || candidates.size() % below.size())
            throw std::invalid_argument("Optuna CPU candidate shape");
        const auto lower = host_pdf(candidates, below), upper = host_pdf(candidates, above);
        std::size_t result = 0;
        double maximum = -INFINITY;
        for (std::size_t row = 0; row < lower.size(); ++row) {
            const double acquisition = lower[row] - upper[row];
            if (!std::isfinite(acquisition))
                throw std::runtime_error("Optuna CPU candidate acquisition is nonfinite");
            if (acquisition > maximum) { maximum = acquisition; result = row; }
        }
        ++stats.cpu_selection_calls;
        return result;
    }
    initialize();
    if (below.empty() || below.size() != above.size() || candidates.empty() || candidates.size() % below.size())
        throw std::invalid_argument("Optuna CUDA candidate shape");
    const auto* lower_model = below.front().prepared;
    const auto* upper_model = above.front().prepared;
    if ((lower_model == nullptr) != (upper_model == nullptr))
        throw std::invalid_argument("Optuna CUDA mixed model ownership");
    std::unique_ptr<Packed> lower, upper;
    if (!lower_model) {
        lower = std::make_unique<Packed>(below, true);
        upper = std::make_unique<Packed>(above, true);
    } else {
        if (lower_model->dimensions != below.size() || upper_model->dimensions != above.size())
            throw std::invalid_argument("Optuna CUDA prepared model dimensions");
        for (std::size_t d = 0; d < below.size(); ++d)
            if (below[d].prepared != lower_model || above[d].prepared != upper_model)
                throw std::invalid_argument("Optuna CUDA prepared model ownership");
    }
    const auto count = candidates.size() / below.size();
    Arena layout;
    const auto input = layout.append(candidates);
    std::unique_ptr<DeviceMixture> lower_mix, upper_mix;
    if (!lower_model) {
        lower_mix = std::make_unique<DeviceMixture>(layout, *lower);
        upper_mix = std::make_unique<DeviceMixture>(layout, *upper);
        lower_mix->normalization = layout.reserve<double>(lower->normalization_count);
        upper_mix->normalization = layout.reserve<double>(upper->normalization_count);
    }
    const auto lower_pdf = layout.reserve<double>(count), upper_pdf = layout.reserve<double>(count);
    const auto index = layout.reserve<std::size_t>(1);
    const auto terms = layout.reserve<double>(product(count,
        std::max(below.front().weights->size(), above.front().weights->size())));
    Buffer arena(layout.bytes);
    stats.peak_temporary_bytes = std::max(stats.peak_temporary_bytes, layout.bytes);
    upload(arena, layout.input.data(), layout.input.size());
    if (lower_model) {
        if (count > static_cast<std::size_t>(INT_MAX))
            throw std::overflow_error("Optuna CUDA grid size overflow");
        const PreparedModel* models[] = {lower_model, upper_model};
        const std::size_t outputs[] = {lower_pdf, upper_pdf};
        for (int i = 0; i < 2; ++i) {
            density<<<static_cast<unsigned>(count), 128>>>(arena.as<double>(input), count,
                models[i]->storage->as<DeviceEstimator>(models[i]->descriptors), below.size(),
                models[i]->storage->as<double>(), models[i]->storage->as<double>(), arena.as<double>(terms), arena.as<double>(outputs[i]));
            check(cudaGetLastError(), "launch prepared density"); ++stats.density_launches;
        }
    } else {
        const DeviceMixture* mixtures[] = {lower_mix.get(), upper_mix.get()};
        const std::size_t outputs[] = {lower_pdf, upper_pdf};
        for (int i = 0; i < 2; ++i) {
            prepare_density<<<static_cast<unsigned>(below.size()), 128>>>(
                arena.as<DeviceEstimator>(mixtures[i]->descriptors), arena.as<double>(mixtures[i]->data), arena.as<double>(mixtures[i]->normalization));
            check(cudaGetLastError(), "prepare candidate density"); ++stats.preparation_launches;
            launch(arena, input, count, below.size(), *mixtures[i], terms, outputs[i]);
        }
    }
    select_best<<<1, 1>>>(arena.as<double>(lower_pdf), arena.as<double>(upper_pdf), count, arena.as<std::size_t>(index));
    check(cudaGetLastError(), "launch selection");
    ++stats.selection_launches;
    std::size_t result;
    download(&result, arena.as<std::size_t>(index), sizeof(result));

    if (result >= count)
        throw std::runtime_error("Optuna CUDA candidate acquisition is nonfinite");
    if (lower_model) {
        lower_model->accept();
        upper_model->accept();
    }
    return result;
}
bool use_gpu() { return gpu_enabled; }
void set_gpu(bool enabled) {
    if (!enabled) release_cached_memory();
    gpu_enabled = enabled;
}
Statistics statistics() { return stats; }
void release_cached_memory() noexcept {
    if (workspaces.empty()) return; // CPU-only use never enters CUDA.
    workspaces.clear();
    stats.cached_device_bytes = 0;
}
void install_runtime(RuntimeHooks hooks) {
    if (!hooks.initialize || !hooks.allocate || !hooks.valid || !hooks.release)
        throw std::invalid_argument("Optuna incomplete runtime hooks");
    if (runtime.initialize) {
        if (runtime.initialize != hooks.initialize || runtime.allocate != hooks.allocate ||
            runtime.valid != hooks.valid || runtime.release != hooks.release)
            throw std::logic_error("Optuna runtime adapter cannot be replaced");
        return;
    }
    runtime = hooks; // registration alone makes no CUDA API call
}
void shutdown_runtime() {
    if (active_buffers) throw std::logic_error("Optuna live buffers prevent runtime unload");
    release_cached_memory();
    runtime = {};
}
}
