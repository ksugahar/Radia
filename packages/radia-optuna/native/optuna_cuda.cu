#include "optuna_cuda.hpp"
#include <cuda_runtime.h>
#include "gpu/mxGPUArray.h"
#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace optuna_cuda {
namespace {
Statistics stats;
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
struct Buffer {
    void* data = nullptr;
    bool stream_ordered = stats.stream_ordered_allocator;
    explicit Buffer(std::size_t bytes) {
        if (bytes) check(stream_ordered ? cudaMallocAsync(&data, bytes, 0) :
                                         cudaMalloc(&data, bytes), "allocate");
    }
    ~Buffer() {
        if (data) {
            if (stream_ordered) cudaFreeAsync(data, 0);
            else cudaFree(data);
        }
    }
    Buffer(const Buffer&) = delete;
    Buffer& operator=(const Buffer&) = delete;
    template<class T> T* as() { return static_cast<T*>(data); }
};
void upload(void* target, const void* source, std::size_t bytes) {
    if (bytes) check(cudaMemcpy(target, source, bytes, cudaMemcpyHostToDevice), "upload");
    stats.host_to_device_bytes += bytes;
}
void download(void* target, const void* source, std::size_t bytes) {
    if (bytes) check(cudaMemcpy(target, source, bytes, cudaMemcpyDeviceToHost), "download");
    stats.device_to_host_bytes += bytes;
}
struct DeviceEstimator {
    std::size_t weights, mu, sigma, probabilities, kernels, choices;
    double low, high, step;
    bool categorical, log_scale;
};
struct Packed {
    std::vector<DeviceEstimator> estimators;
    std::vector<double> values;
    std::size_t append(const std::vector<double>* value) {
        const auto offset = values.size();
        if (value) values.insert(values.end(), value->begin(), value->end());
        return offset;
    }
    explicit Packed(const std::vector<Estimator>& source) {
        if (source.empty()) throw std::invalid_argument("Optuna CUDA requires estimators");
        const auto kernels = source.front().weights->size();
        for (const auto& e : source) {
            if (!kernels || e.weights->size() != kernels)
                throw std::invalid_argument("Optuna CUDA mixture sizes disagree");
            estimators.push_back({append(e.weights), append(e.mu), append(e.sigma),
                append(e.probabilities), kernels, e.choices, e.low, e.high,
                e.step, e.categorical, e.log_scale});
        }
    }
};
__device__ double log_mass(double lower, double upper) {
    constexpr double inv = 0.70710678118654752440;
    double mass;
    if (upper <= 0.0) mass = 0.5 * (erfc(-upper * inv) - erfc(-lower * inv));
    else if (lower >= 0.0) mass = 0.5 * (erfc(lower * inv) - erfc(upper * inv));
    else mass = 0.5 * (erf(upper * inv) - erf(lower * inv));
    return log(fmax(mass, DBL_MIN));
}
__device__ double component(const DeviceEstimator& e, const double* data,
                            double sample, std::size_t kernel) {
    if (e.categorical) {
        const auto choice = static_cast<std::size_t>(sample) - 1;
        return log(fmax(data[e.probabilities + kernel * e.choices + choice], DBL_MIN));
    }
    const double mu = data[e.mu + kernel], sigma = data[e.sigma + kernel];
    const double denominator = log_mass((e.low - mu) / sigma, (e.high - mu) / sigma);
    if (isfinite(e.step)) {
        double lo = sample - 0.5 * e.step, hi = sample + 0.5 * e.step;
        if (e.log_scale) { lo = log(lo); hi = log(hi); }
        return log_mass((lo - mu) / sigma, (hi - mu) / sigma) - denominator;
    }
    const double x = e.log_scale ? log(sample) : sample;
    const double z = (x - mu) / sigma;
    return -0.5 * z * z - 0.91893853320467274178 - log(sigma) - denominator;
}
// A block per candidate evaluates mixture components in parallel. The final
// ordered reduction retains the oracle's summation and tie-breaking contract.
__global__ void density(const double* samples, std::size_t count,
                         const DeviceEstimator* estimators, std::size_t dimensions,
                         const double* data, double* terms, double* output) {
    const std::size_t row = blockIdx.x;
    if (row >= count) return;
    const auto kernels = estimators[0].kernels;
    double* components = terms + row * kernels;
    for (std::size_t k = threadIdx.x; k < kernels; k += blockDim.x) {
        double value = log(data[estimators[0].weights + k]);
        for (std::size_t d = 0; d < dimensions; ++d)
            value += component(estimators[d], data, samples[row * dimensions + d], k);
        components[k] = value;
    }
    __syncthreads();
    if (threadIdx.x != 0) return;
    double maximum = -INFINITY;
    for (std::size_t k = 0; k < kernels; ++k) maximum = fmax(maximum, components[k]);
    if (maximum == -INFINITY) { output[row] = maximum; return; }
    double total = 0.0;
    for (std::size_t k = 0; k < kernels; ++k) total += exp(components[k] - maximum);
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
void initialize() {
    mxInitGPU();
    check(cudaGetDevice(&stats.device), "get device");
    check(cudaRuntimeGetVersion(&stats.runtime_version), "runtime version");
    check(cudaDriverGetVersion(&stats.driver_version), "driver version");
    int pool_supported = 0;
    check(cudaDeviceGetAttribute(&pool_supported, cudaDevAttrMemoryPoolsSupported,
                                 stats.device), "memory pool support");
    stats.stream_ordered_allocator = pool_supported != 0;
}
struct DeviceMixture {
    Buffer descriptors, data;
    explicit DeviceMixture(const Packed& p)
        : descriptors(product(p.estimators.size(), sizeof(DeviceEstimator))),
          data(product(p.values.size(), sizeof(double))) {
        upload(descriptors.data, p.estimators.data(), product(p.estimators.size(), sizeof(DeviceEstimator)));
        upload(data.data, p.values.data(), product(p.values.size(), sizeof(double)));
    }
};
void launch(Buffer& samples, std::size_t count, std::size_t dimensions,
            DeviceMixture& mixture, Buffer& terms, Buffer& output) {
    if (count > static_cast<std::size_t>(INT_MAX))
        throw std::overflow_error("Optuna CUDA grid size overflow");
    density<<<static_cast<unsigned>(count), 128>>>(
        samples.as<double>(), count, mixture.descriptors.as<DeviceEstimator>(),
        dimensions, mixture.data.as<double>(), terms.as<double>(), output.as<double>());
    check(cudaGetLastError(), "launch density");
    ++stats.density_launches;
}
}
std::vector<double> log_pdf(const std::vector<double>& samples,
                            const std::vector<Estimator>& estimators) {
    initialize();
    Packed packed(estimators);
    if (samples.size() % estimators.size()) throw std::invalid_argument("Optuna CUDA sample shape");
    const auto count = samples.size() / estimators.size();
    std::vector<double> result(count);
    if (!count) return result;
    const auto input_bytes = product(samples.size(), sizeof(double));
    const auto output_bytes = product(count, sizeof(double));
    const auto term_bytes = product(product(count, packed.estimators.front().kernels), sizeof(double));
    Buffer input(input_bytes), output(output_bytes), terms(term_bytes);
    DeviceMixture mixture(packed);
    stats.peak_temporary_bytes = std::max(stats.peak_temporary_bytes, input_bytes + output_bytes + term_bytes +
        product(packed.values.size(), sizeof(double)) + product(packed.estimators.size(), sizeof(DeviceEstimator)));
    upload(input.data, samples.data(), input_bytes);
    launch(input, count, estimators.size(), mixture, terms, output);
    download(result.data(), output.data, output_bytes);
    for (double value : result)
        if (!std::isfinite(value)) throw std::runtime_error("Optuna CUDA density is nonfinite");
    return result;
}
std::size_t best(const std::vector<double>& candidates,
                 const std::vector<Estimator>& below, const std::vector<Estimator>& above) {
    initialize();
    if (below.empty() || below.size() != above.size() || candidates.empty() || candidates.size() % below.size())
        throw std::invalid_argument("Optuna CUDA candidate shape");
    Packed lower(below), upper(above);
    const auto count = candidates.size() / below.size();
    const auto input_bytes = product(candidates.size(), sizeof(double));
    const auto output_bytes = product(count, sizeof(double));
    const auto term_bytes = product(product(count,
        std::max(lower.estimators.front().kernels, upper.estimators.front().kernels)), sizeof(double));
    Buffer input(input_bytes), lower_pdf(output_bytes), upper_pdf(output_bytes), index(sizeof(std::size_t));
    Buffer terms(term_bytes);
    DeviceMixture lower_mix(lower), upper_mix(upper);
    stats.peak_temporary_bytes = std::max(stats.peak_temporary_bytes, input_bytes + 2 * output_bytes + term_bytes + sizeof(std::size_t) +
        product(lower.values.size() + upper.values.size(), sizeof(double)) +
        product(lower.estimators.size() + upper.estimators.size(), sizeof(DeviceEstimator)));
    upload(input.data, candidates.data(), input_bytes);
    launch(input, count, below.size(), lower_mix, terms, lower_pdf);
    launch(input, count, above.size(), upper_mix, terms, upper_pdf);
    select_best<<<1, 1>>>(lower_pdf.as<double>(), upper_pdf.as<double>(), count, index.as<std::size_t>());
    check(cudaGetLastError(), "launch selection");
    ++stats.selection_launches;
    std::size_t result;
    download(&result, index.data, sizeof(result));
    if (result >= count)
        throw std::runtime_error("Optuna CUDA candidate acquisition is nonfinite");
    return result;
}
Statistics statistics() { return stats; }
}
