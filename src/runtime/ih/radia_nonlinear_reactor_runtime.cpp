#include "radia_nonlinear_reactor_runtime.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace radia { namespace reactor {
namespace {

constexpr double kMu0 = 4.0e-7 * 3.141592653589793238462643383279502884;

bool finite_values(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(),
                       [](double value) { return std::isfinite(value); });
}

double dot(const std::vector<double>& left,
           const std::vector<double>& right) {
    double result = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index)
        result += left[index] * right[index];
    return result;
}

double norm(const std::vector<double>& values) {
    return std::sqrt(std::max(dot(values, values), 0.0));
}

std::vector<double> matvec(const std::vector<double>& matrix,
                           const std::vector<double>& vector, int n) {
    std::vector<double> result(static_cast<std::size_t>(n), 0.0);
    for (int row = 0; row < n; ++row)
        for (int column = 0; column < n; ++column)
            result[static_cast<std::size_t>(row)] +=
                matrix[static_cast<std::size_t>(row * n + column)] *
                vector[static_cast<std::size_t>(column)];
    return result;
}

std::vector<double> solve_linear(std::vector<double> matrix,
                                 std::vector<double> rhs, int n) {
    double matrix_scale = 0.0;
    for (double value : matrix) matrix_scale = std::max(matrix_scale, std::abs(value));
    const double pivot_floor = std::max(
        matrix_scale * 100.0 * std::numeric_limits<double>::epsilon(),
        std::numeric_limits<double>::min());
    for (int column = 0; column < n; ++column) {
        int pivot = column;
        double best = std::abs(
            matrix[static_cast<std::size_t>(column * n + column)]);
        for (int row = column + 1; row < n; ++row) {
            const double candidate = std::abs(
                matrix[static_cast<std::size_t>(row * n + column)]);
            if (candidate > best) {
                best = candidate;
                pivot = row;
            }
        }
        if (!(best > pivot_floor) || !std::isfinite(best))
            throw std::runtime_error(
                "nonlinear reactor tangent is singular or ill-conditioned");
        if (pivot != column) {
            for (int entry = column; entry < n; ++entry)
                std::swap(
                    matrix[static_cast<std::size_t>(column * n + entry)],
                    matrix[static_cast<std::size_t>(pivot * n + entry)]);
            std::swap(rhs[static_cast<std::size_t>(column)],
                      rhs[static_cast<std::size_t>(pivot)]);
        }
        for (int row = column + 1; row < n; ++row) {
            const double factor =
                matrix[static_cast<std::size_t>(row * n + column)] /
                matrix[static_cast<std::size_t>(column * n + column)];
            matrix[static_cast<std::size_t>(row * n + column)] = 0.0;
            for (int entry = column + 1; entry < n; ++entry)
                matrix[static_cast<std::size_t>(row * n + entry)] -=
                    factor * matrix[static_cast<std::size_t>(column * n + entry)];
            rhs[static_cast<std::size_t>(row)] -=
                factor * rhs[static_cast<std::size_t>(column)];
        }
    }
    std::vector<double> solution(static_cast<std::size_t>(n), 0.0);
    for (int row = n - 1; row >= 0; --row) {
        double value = rhs[static_cast<std::size_t>(row)];
        for (int column = row + 1; column < n; ++column)
            value -= matrix[static_cast<std::size_t>(row * n + column)] *
                     solution[static_cast<std::size_t>(column)];
        solution[static_cast<std::size_t>(row)] =
            value / matrix[static_cast<std::size_t>(row * n + row)];
    }
    return solution;
}

struct ConstitutiveValue {
    double field = 0.0;
    double differential_reluctivity = 0.0;
    double coenergy = 0.0;
};

ConstitutiveValue constitutive(const Config& config, double magnetization) {
    const auto& m = config.magnetization_table_A_per_m;
    const auto& h = config.field_table_A_per_m;
    const double magnitude = std::max(magnetization, 0.0);
    std::size_t upper = 1;
    if (magnitude >= m.back()) {
        upper = m.size() - 1;
    } else {
        upper = static_cast<std::size_t>(
            std::upper_bound(m.begin(), m.end(), magnitude) - m.begin());
        upper = std::max<std::size_t>(upper, 1);
    }
    const std::size_t lower = upper - 1;
    const double slope = (h[upper] - h[lower]) / (m[upper] - m[lower]);
    double field = h[lower] + slope * (magnitude - m[lower]);
    if (magnitude >= m.back())
        field = h.back() + slope * (magnitude - m.back());

    double integral = 0.0;
    for (std::size_t index = 1; index <= lower; ++index)
        integral += 0.5 * (h[index] + h[index - 1]) *
                    (m[index] - m[index - 1]);
    const double delta = magnitude - m[lower];
    integral += h[lower] * delta + 0.5 * slope * delta * delta;
    return {field, slope, integral};
}

struct Evaluation {
    std::vector<double> residual;
    std::vector<double> tangent;
    std::vector<double> sample_magnitude;
    std::vector<double> sample_field;
    double objective = 0.0;
    double constitutive_coenergy = 0.0;
};

Evaluation evaluate(const Config& config, const std::vector<double>& coefficients,
                    double current_A) {
    const int n_modes = config.n_modes;
    const int n_samples = config.n_samples;
    Evaluation result;
    result.residual = matvec(config.demag, coefficients, n_modes);
    result.tangent = config.demag;
    result.sample_magnitude.assign(static_cast<std::size_t>(n_samples), 0.0);
    result.sample_field.assign(static_cast<std::size_t>(n_samples), 0.0);

    for (int mode = 0; mode < n_modes; ++mode)
        result.residual[static_cast<std::size_t>(mode)] -=
            current_A * config.excitation_per_amp[static_cast<std::size_t>(mode)];

    for (int sample = 0; sample < n_samples; ++sample) {
        double vector[3] = {0.0, 0.0, 0.0};
        for (int mode = 0; mode < n_modes; ++mode)
            for (int component = 0; component < 3; ++component)
                vector[component] +=
                    coefficients[static_cast<std::size_t>(mode)] *
                    config.magnetization_modes[static_cast<std::size_t>(
                        (mode * n_samples + sample) * 3 + component)];
        const double magnitude = std::sqrt(
            vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
        const ConstitutiveValue material = constitutive(config, magnitude);
        const double secant = magnitude > 1.0e-30
            ? material.field / magnitude
            : material.differential_reluctivity;
        const double direction[3] = {
            magnitude > 1.0e-30 ? vector[0] / magnitude : 0.0,
            magnitude > 1.0e-30 ? vector[1] / magnitude : 0.0,
            magnitude > 1.0e-30 ? vector[2] / magnitude : 0.0};
        const double weight = config.sample_weights[static_cast<std::size_t>(sample)];
        result.sample_magnitude[static_cast<std::size_t>(sample)] = magnitude;
        result.sample_field[static_cast<std::size_t>(sample)] = material.field;
        result.constitutive_coenergy += weight * material.coenergy;

        std::vector<double> parallel(static_cast<std::size_t>(n_modes), 0.0);
        for (int mode = 0; mode < n_modes; ++mode) {
            double projection = 0.0;
            double field_projection = 0.0;
            for (int component = 0; component < 3; ++component) {
                const double basis = config.magnetization_modes[
                    static_cast<std::size_t>((mode * n_samples + sample) * 3 + component)];
                projection += basis * direction[component];
                field_projection += basis * vector[component];
            }
            parallel[static_cast<std::size_t>(mode)] = projection;
            result.residual[static_cast<std::size_t>(mode)] +=
                weight * secant * field_projection;
        }
        for (int row = 0; row < n_modes; ++row) {
            for (int column = 0; column < n_modes; ++column) {
                double basis_dot = 0.0;
                for (int component = 0; component < 3; ++component)
                    basis_dot += config.magnetization_modes[static_cast<std::size_t>(
                                     (row * n_samples + sample) * 3 + component)] *
                                 config.magnetization_modes[static_cast<std::size_t>(
                                     (column * n_samples + sample) * 3 + component)];
                result.tangent[static_cast<std::size_t>(row * n_modes + column)] +=
                    weight * (secant * basis_dot +
                              (material.differential_reluctivity - secant) *
                                  parallel[static_cast<std::size_t>(row)] *
                                  parallel[static_cast<std::size_t>(column)]);
            }
        }
    }

    const std::vector<double> demag_m =
        matvec(config.demag, coefficients, n_modes);
    result.objective = result.constitutive_coenergy +
        0.5 * dot(coefficients, demag_m) -
        current_A * dot(config.excitation_per_amp, coefficients);
    return result;
}

}  // namespace

Runtime::Runtime(Config config) : config_(std::move(config)) {
    validate();
    previous_magnetization_.assign(static_cast<std::size_t>(config_.n_modes), 0.0);
    reset();
}

void Runtime::validate() const {
    const std::size_t n_modes = static_cast<std::size_t>(config_.n_modes);
    const std::size_t n_samples = static_cast<std::size_t>(config_.n_samples);
    if (config_.n_modes <= 0 || config_.n_samples <= 0 ||
        config_.demag.size() != n_modes * n_modes ||
        config_.magnetization_modes.size() != n_modes * n_samples * 3 ||
        config_.sample_weights.size() != n_samples ||
        config_.excitation_per_amp.size() != n_modes ||
        config_.magnetization_table_A_per_m.size() < 2 ||
        config_.magnetization_table_A_per_m.size() !=
            config_.field_table_A_per_m.size())
        throw std::invalid_argument("invalid nonlinear reactor dimensions");
    if (!finite_values(config_.demag) ||
        !finite_values(config_.magnetization_modes) ||
        !finite_values(config_.sample_weights) ||
        !finite_values(config_.excitation_per_amp) ||
        !finite_values(config_.magnetization_table_A_per_m) ||
        !finite_values(config_.field_table_A_per_m))
        throw std::invalid_argument(
            "nonlinear reactor configuration must contain finite values");
    if (!std::all_of(config_.sample_weights.begin(), config_.sample_weights.end(),
                     [](double value) { return value > 0.0; }))
        throw std::invalid_argument("sample weights must be positive");
    if (norm(config_.excitation_per_amp) == 0.0)
        throw std::invalid_argument("reactor excitation must be nonzero");
    for (int row = 0; row < config_.n_modes; ++row)
        for (int column = 0; column < config_.n_modes; ++column) {
            const double left = config_.demag[static_cast<std::size_t>(
                row * config_.n_modes + column)];
            const double right = config_.demag[static_cast<std::size_t>(
                column * config_.n_modes + row)];
            if (std::abs(left - right) >
                1.0e-10 * std::max({1.0, std::abs(left), std::abs(right)}))
                throw std::invalid_argument("demag operator must be symmetric");
        }
    const auto& m = config_.magnetization_table_A_per_m;
    const auto& h = config_.field_table_A_per_m;
    if (m.front() != 0.0 || h.front() != 0.0)
        throw std::invalid_argument("inverse BH table must start at zero");
    for (std::size_t index = 1; index < m.size(); ++index)
        if (!(m[index] > m[index - 1]) || !(h[index] > h[index - 1]))
            throw std::invalid_argument(
                "inverse BH magnetization and field values must be strictly increasing");
    if (!(config_.sample_time_s > 0.0) ||
        !(config_.residual_tolerance > 0.0) || config_.max_iterations < 1 ||
        !(config_.line_search_minimum > 0.0 && config_.line_search_minimum <= 1.0) ||
        config_.air_inductance_H < 0.0 || config_.winding_resistance_Ohm < 0.0 ||
        !std::isfinite(config_.sample_time_s) ||
        !std::isfinite(config_.initial_current_A) ||
        !std::isfinite(config_.air_inductance_H) ||
        !std::isfinite(config_.winding_resistance_Ohm))
        throw std::invalid_argument("invalid nonlinear reactor scalar configuration");
}

Output Runtime::solve(double current_A,
                      const std::vector<double>& initial) const {
    if (!std::isfinite(current_A))
        throw std::invalid_argument("reactor current must be finite");
    std::vector<double> coefficients = initial;
    Evaluation state = evaluate(config_, coefficients, current_A);
    const double rhs_scale = std::max(
        std::abs(current_A) * norm(config_.excitation_per_amp), 1.0);
    double relative = norm(state.residual) / rhs_scale;
    int iterations = 0;
    for (; relative > config_.residual_tolerance &&
           iterations < config_.max_iterations; ++iterations) {
        std::vector<double> rhs = state.residual;
        for (double& value : rhs) value = -value;
        const std::vector<double> step =
            solve_linear(state.tangent, std::move(rhs), config_.n_modes);
        const double directional = dot(state.residual, step);
        double factor = 1.0;
        bool accepted = false;
        while (factor >= config_.line_search_minimum) {
            std::vector<double> candidate = coefficients;
            for (int mode = 0; mode < config_.n_modes; ++mode)
                candidate[static_cast<std::size_t>(mode)] +=
                    factor * step[static_cast<std::size_t>(mode)];
            Evaluation trial = evaluate(config_, candidate, current_A);
            if (trial.objective <=
                    state.objective + 1.0e-4 * factor * directional ||
                norm(trial.residual) < norm(state.residual)) {
                coefficients = std::move(candidate);
                state = std::move(trial);
                accepted = true;
                break;
            }
            factor *= 0.5;
        }
        if (!accepted)
            throw std::runtime_error(
                "nonlinear reactor Newton line search failed");
        relative = norm(state.residual) / rhs_scale;
    }
    if (relative > config_.residual_tolerance)
        throw std::runtime_error(
            "nonlinear reactor solve did not converge within max_iterations");

    Output output;
    output.magnetization_coefficients = coefficients;
    output.residual_relative_norm = relative;
    output.nonlinear_iterations = iterations;
    output.flux_linkage_Wb_turn =
        config_.air_inductance_H * current_A +
        kMu0 * dot(config_.excitation_per_amp, coefficients);
    output.voltage_V = config_.winding_resistance_Ohm * current_A +
        (output.flux_linkage_Wb_turn - previous_flux_linkage_) /
            config_.sample_time_s;

    const std::vector<double> sensitivity = solve_linear(
        state.tangent, config_.excitation_per_amp, config_.n_modes);
    output.differential_inductance_H = config_.air_inductance_H +
        kMu0 * dot(config_.excitation_per_amp, sensitivity);
    const std::vector<double> demag_m =
        matvec(config_.demag, coefficients, config_.n_modes);
    output.magnetic_energy_J = 0.5 * config_.air_inductance_H *
        current_A * current_A + kMu0 *
        (state.constitutive_coenergy + 0.5 * dot(coefficients, demag_m));
    output.flux_density_T.resize(static_cast<std::size_t>(config_.n_samples));
    for (int sample = 0; sample < config_.n_samples; ++sample) {
        const double value = kMu0 *
            (state.sample_magnitude[static_cast<std::size_t>(sample)] +
             state.sample_field[static_cast<std::size_t>(sample)]);
        output.flux_density_T[static_cast<std::size_t>(sample)] = value;
        output.peak_flux_density_T = std::max(output.peak_flux_density_T, value);
    }
    return output;
}

const Output& Runtime::output(double current_A) {
    if (!cache_valid_ || current_A != cached_current_A_) {
        cached_output_ = solve(current_A, previous_magnetization_);
        cached_current_A_ = current_A;
        cache_valid_ = true;
    }
    return cached_output_;
}

void Runtime::update(double current_A) {
    const Output& trial = output(current_A);
    previous_current_A_ = current_A;
    previous_flux_linkage_ = trial.flux_linkage_Wb_turn;
    previous_magnetization_ = trial.magnetization_coefficients;
    ++accepted_steps_;
    cache_valid_ = false;
}

void Runtime::reset() {
    previous_magnetization_.assign(static_cast<std::size_t>(config_.n_modes), 0.0);
    previous_current_A_ = config_.initial_current_A;
    previous_flux_linkage_ = 0.0;
    accepted_steps_ = 0;
    cache_valid_ = false;
    Output initial = solve(config_.initial_current_A, previous_magnetization_);
    previous_magnetization_ = initial.magnetization_coefficients;
    previous_flux_linkage_ = initial.flux_linkage_Wb_turn;
    cache_valid_ = false;
}

Snapshot Runtime::snapshot() const {
    return {previous_current_A_, previous_flux_linkage_, accepted_steps_,
            previous_magnetization_};
}

void Runtime::restore(const Snapshot& state) {
    if (!std::isfinite(state.previous_current_A) ||
        !std::isfinite(state.previous_flux_linkage_Wb_turn) ||
        state.magnetization_coefficients.size() !=
            static_cast<std::size_t>(config_.n_modes) ||
        !finite_values(state.magnetization_coefficients))
        throw std::invalid_argument("invalid nonlinear reactor snapshot");
    previous_current_A_ = state.previous_current_A;
    previous_flux_linkage_ = state.previous_flux_linkage_Wb_turn;
    accepted_steps_ = state.accepted_steps;
    previous_magnetization_ = state.magnetization_coefficients;
    cache_valid_ = false;
}

}}  // namespace radia::reactor
