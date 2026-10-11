#pragma once
#include <vector>

namespace optuna_cuda {
struct Dimension {
    std::size_t mu, sigma, denominator, probabilities, choices;
    bool categorical, log_scale;
    double step;
};
std::vector<double> joint(const std::vector<double>& samples, std::size_t count,
    const std::vector<double>& weights, const std::vector<double>& parameters,
    const std::vector<Dimension>& dimensions);
// Checked host inputs, row-major categorical probabilities. No solver runtime.
std::vector<double> numerical(const std::vector<double>& samples,
    const std::vector<double>& weights, const std::vector<double>& mu,
    const std::vector<double>& sigma, const std::vector<double>& denominator,
    bool log_scale, double step);
std::vector<double> categorical(const std::vector<double>& samples,
    const std::vector<double>& weights, const std::vector<double>& probabilities,
    std::size_t choices);
}
