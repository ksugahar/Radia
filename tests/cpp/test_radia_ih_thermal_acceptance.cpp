#include "radia_ih_thermal.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using radia::ih::CSRMatrix;
using radia::ih::ThermalState;
using radia::ih::ThermalStepOptions;

CSRMatrix Diagonal(int n, double value) {
    CSRMatrix result;
    result.n = n;
    result.row_ptr.resize(static_cast<std::size_t>(n + 1));
    result.col.resize(static_cast<std::size_t>(n));
    result.value.assign(static_cast<std::size_t>(n), value);
    for (int index = 0; index < n; ++index) {
        result.row_ptr[static_cast<std::size_t>(index)] = index;
        result.col[static_cast<std::size_t>(index)] = index;
    }
    result.row_ptr[static_cast<std::size_t>(n)] = n;
    return result;
}

CSRMatrix DenseTwoByTwo(double a00, double a01, double a10, double a11) {
    return CSRMatrix{2, {0, 2, 4}, {0, 1, 0, 1}, {a00, a01, a10, a11}};
}

void RequireEqual(double actual, double expected, const std::string& name) {
    if (actual != expected)
        throw std::runtime_error(name + " changed after rejected solve");
}

void TestTinyLoadCannotUseAnAbsoluteResidualFloor() {
    ThermalState state{{300.0}, 7.0, 0.25};
    ThermalStepOptions options;
    options.dt_s = 1.0;
    options.tolerance = 1.0e-6;
    options.max_iterations = 0;
    bool rejected = false;
    try {
        radia::ih::advance_thermal(
            Diagonal(1, 1.0), Diagonal(1, 0.0), nullptr,
            {1.0e-10}, {1.0}, 0.75, options, state);
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    if (!rejected)
        throw std::runtime_error("tiny unresolved thermal load was accepted");
    RequireEqual(state.temperature_K[0], 300.0, "temperature");
    RequireEqual(state.time_s, 7.0, "time");
    RequireEqual(state.previous_angle_rad, 0.25, "angle");
}

void TestFailedIterationDoesNotPublishPartialState() {
    ThermalState state{{300.0, 310.0}, 2.0, 0.1};
    ThermalStepOptions options;
    options.dt_s = 1.0;
    options.tolerance = 1.0e-12;
    options.max_iterations = 1;
    bool rejected = false;
    try {
        radia::ih::advance_thermal(
            DenseTwoByTwo(1.0, 0.0, 0.0, 1.0),
            DenseTwoByTwo(2.0, -1.0, -1.0, 2.0), nullptr,
            {4.0, -3.0}, {1.0, 1.0}, 0.9, options, state);
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    if (!rejected)
        throw std::runtime_error("under-iterated thermal solve was accepted");
    RequireEqual(state.temperature_K[0], 300.0, "first temperature");
    RequireEqual(state.temperature_K[1], 310.0, "second temperature");
    RequireEqual(state.time_s, 2.0, "time");
    RequireEqual(state.previous_angle_rad, 0.1, "angle");
}

void TestJacobiPcgAcceptsAResolvedSystem() {
    ThermalState state{{300.0, 310.0}, 0.0, 0.0};
    ThermalStepOptions options;
    options.dt_s = 1.0;
    options.tolerance = 1.0e-12;
    options.max_iterations = 10;
    radia::ih::advance_thermal(
        DenseTwoByTwo(1.0, 0.0, 0.0, 1.0),
        DenseTwoByTwo(2.0, -1.0, -1.0, 2.0), nullptr,
        {4.0, -3.0}, {1.0, 1.0}, 0.9, options, state);
    const std::vector<double> rhs = {304.0, 307.0};
    const std::vector<double> residual = {
        rhs[0] - (3.0 * state.temperature_K[0] - state.temperature_K[1]),
        rhs[1] - (-state.temperature_K[0] + 3.0 * state.temperature_K[1])};
    const double relative = std::hypot(residual[0], residual[1]) /
                            std::hypot(rhs[0], rhs[1]);
    if (!(relative <= 1.0e-12))
        throw std::runtime_error("accepted solve failed the independent true residual");
}
}  // namespace

int main() {
    try {
        TestTinyLoadCannotUseAnAbsoluteResidualFloor();
        TestFailedIterationDoesNotPublishPartialState();
        TestJacobiPcgAcceptsAResolvedSystem();
        std::cout << "radia_ih_thermal: all tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "radia_ih_thermal: FAILED: " << error.what() << '\n';
        return 1;
    }
}
