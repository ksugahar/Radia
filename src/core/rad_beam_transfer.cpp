#include "rad_beam_transfer.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

namespace radia::beam {
namespace {

template <typename Value>
void AddScaled(Value& destination, const Value& source, double scale) {
    for (std::size_t index = 0; index < destination.values.size(); ++index)
        destination.values[index] += scale * source.values[index];
}

template <typename Value>
Value Difference(const Value& left, const Value& right) {
    Value result;
    for (std::size_t index = 0; index < result.values.size(); ++index)
        result.values[index] = left.values[index] - right.values[index];
    return result;
}

template <typename Value>
void RequireFinite(const Value& value, const char* name) {
    for (double item : value.values) {
        if (!std::isfinite(item))
            throw std::invalid_argument(std::string(name) +
                                        " must contain finite values");
    }
}

Tensor3Map6 LeftMultiply(const Matrix6& left, const Tensor3Map6& right) {
    Tensor3Map6 result;
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t a = 0; a < 6; ++a) {
            const double coefficient = left(i, a);
            for (std::size_t j = 0; j < 6; ++j)
                for (std::size_t k = 0; k < 6; ++k)
                    result(i, j, k) += coefficient * right(a, j, k);
        }
    return result;
}

Tensor4Map6 LeftMultiply(const Matrix6& left, const Tensor4Map6& right) {
    Tensor4Map6 result;
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t a = 0; a < 6; ++a) {
            const double coefficient = left(i, a);
            for (std::size_t j = 0; j < 6; ++j)
                for (std::size_t k = 0; k < 6; ++k)
                    for (std::size_t l = 0; l < 6; ++l)
                        result(i, j, k, l) +=
                            coefficient * right(a, j, k, l);
        }
    return result;
}

Tensor3Map6 TransformInputs(const Tensor3Map6& tensor,
                            const Matrix6& transform) {
    Tensor3Map6 first;
    Tensor3Map6 result;
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t a = 0; a < 6; ++a)
            for (std::size_t b = 0; b < 6; ++b)
                for (std::size_t j = 0; j < 6; ++j)
                    first(i, j, b) += tensor(i, a, b) * transform(a, j);
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            for (std::size_t b = 0; b < 6; ++b)
                for (std::size_t k = 0; k < 6; ++k)
                    result(i, j, k) += first(i, j, b) * transform(b, k);
    return result;
}

Tensor4Map6 TransformInputs(const Tensor4Map6& tensor,
                            const Matrix6& transform) {
    Tensor4Map6 first;
    Tensor4Map6 second;
    Tensor4Map6 result;
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t a = 0; a < 6; ++a)
            for (std::size_t b = 0; b < 6; ++b)
                for (std::size_t c = 0; c < 6; ++c)
                    for (std::size_t j = 0; j < 6; ++j)
                        first(i, j, b, c) +=
                            tensor(i, a, b, c) * transform(a, j);
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            for (std::size_t b = 0; b < 6; ++b)
                for (std::size_t c = 0; c < 6; ++c)
                    for (std::size_t k = 0; k < 6; ++k)
                        second(i, j, k, c) +=
                            first(i, j, b, c) * transform(b, k);
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            for (std::size_t k = 0; k < 6; ++k)
                for (std::size_t c = 0; c < 6; ++c)
                    for (std::size_t l = 0; l < 6; ++l)
                        result(i, j, k, l) +=
                            second(i, j, k, c) * transform(c, l);
    return result;
}

// Returns the complete 3 * outer[R,T] term under the factorial convention.
Tensor4Map6 CrossSecondOrder(const Tensor3Map6& outer,
                             const Matrix6& inner_r,
                             const Tensor3Map6& inner_t) {
    Tensor3Map6 first;
    Tensor4Map6 result;
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t a = 0; a < 6; ++a)
            for (std::size_t b = 0; b < 6; ++b)
                for (std::size_t j = 0; j < 6; ++j)
                    first(i, j, b) += outer(i, a, b) * inner_r(a, j);

    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            for (std::size_t k = 0; k < 6; ++k)
                for (std::size_t l = 0; l < 6; ++l)
                    for (std::size_t b = 0; b < 6; ++b)
                        result(i, j, k, l) +=
                            first(i, j, b) * inner_t(b, k, l) +
                            first(i, k, b) * inner_t(b, j, l) +
                            first(i, l, b) * inner_t(b, j, k);
    return result;
}

TaylorMap6 MapDerivative(const DynamicsJet6& jet, const TaylorMap6& map,
                         unsigned maximum_order) {
    TaylorMap6 result;
    result.r = Multiply(jet.a_per_m, map.r);
    if (maximum_order >= 2) {
        result.t = LeftMultiply(jet.a_per_m, map.t);
        AddScaled(result.t, TransformInputs(jet.f2_per_m, map.r), 1.0);
    }
    if (maximum_order >= 3) {
        result.u = LeftMultiply(jet.a_per_m, map.u);
        AddScaled(result.u,
                  CrossSecondOrder(jet.f2_per_m, map.r, map.t), 1.0);
        AddScaled(result.u, TransformInputs(jet.f3_per_m, map.r), 1.0);
    }
    return result;
}

TaylorMap6 StatePlus(const TaylorMap6& state, const TaylorMap6& derivative,
                     double scale, unsigned maximum_order) {
    TaylorMap6 result = state;
    AddScaled(result.r, derivative.r, scale);
    if (maximum_order >= 2) AddScaled(result.t, derivative.t, scale);
    if (maximum_order >= 3) AddScaled(result.u, derivative.u, scale);
    return result;
}

TaylorMap6 IntegrateConstantJetUnchecked(const DynamicsJet6& jet,
                                         double length_m,
                                         unsigned maximum_order) {
    TaylorMap6 state = IdentityTaylorMap6();
    const TaylorMap6 k1 = MapDerivative(jet, state, maximum_order);
    const TaylorMap6 k2 = MapDerivative(
        jet, StatePlus(state, k1, 0.5 * length_m, maximum_order),
        maximum_order);
    const TaylorMap6 k3 = MapDerivative(
        jet, StatePlus(state, k2, 0.5 * length_m, maximum_order),
        maximum_order);
    const TaylorMap6 k4 = MapDerivative(
        jet, StatePlus(state, k3, length_m, maximum_order), maximum_order);
    AddScaled(state.r, k1.r, length_m / 6.0);
    AddScaled(state.r, k2.r, length_m / 3.0);
    AddScaled(state.r, k3.r, length_m / 3.0);
    AddScaled(state.r, k4.r, length_m / 6.0);
    if (maximum_order >= 2) {
        AddScaled(state.t, k1.t, length_m / 6.0);
        AddScaled(state.t, k2.t, length_m / 3.0);
        AddScaled(state.t, k3.t, length_m / 3.0);
        AddScaled(state.t, k4.t, length_m / 6.0);
    }
    if (maximum_order >= 3) {
        AddScaled(state.u, k1.u, length_m / 6.0);
        AddScaled(state.u, k2.u, length_m / 3.0);
        AddScaled(state.u, k3.u, length_m / 3.0);
        AddScaled(state.u, k4.u, length_m / 6.0);
    }
    return state;
}

struct StepDescriptor {
    std::size_t region = 0;
    double path_begin_m = 0.0;
    double length_m = 0.0;
};

double SmallFactorial(std::size_t value) {
    if (value < 2) return 1.0;
    if (value == 2) return 2.0;
    if (value == 3) return 6.0;
    if (value == 4) return 24.0;
    if (value == 5) return 120.0;
    throw std::invalid_argument("beam Taylor monomial degree exceeds five");
}

void AddHamiltonianMonomial(HamiltonianJet6& jet,
                            const std::array<std::size_t, 6>& powers,
                            double coefficient) {
    std::size_t degree = 0;
    double derivative = coefficient;
    for (std::size_t power : powers) {
        degree += power;
        derivative *= SmallFactorial(power);
    }
    if (coefficient == 0.0) return;
    if (degree < 2 || degree > 5)
        throw std::invalid_argument(
            "Hamiltonian monomial degree must be between two and five");
    const auto matches = [&](const std::array<std::size_t, 5>& inputs,
                             std::size_t count) {
        std::array<std::size_t, 6> actual{};
        for (std::size_t index = 0; index < count; ++index)
            ++actual[inputs[index]];
        return actual == powers;
    };
    if (degree == 2) {
        for (std::size_t i = 0; i < 6; ++i)
            for (std::size_t j = 0; j < 6; ++j)
                if (matches({i, j, 0, 0, 0}, 2))
                    jet.h2_per_m(i, j) += derivative;
        return;
    }
    if (degree == 3) {
        for (std::size_t i = 0; i < 6; ++i)
            for (std::size_t j = 0; j < 6; ++j)
                for (std::size_t k = 0; k < 6; ++k)
                    if (matches({i, j, k, 0, 0}, 3))
                        jet.h3_per_m(i, j, k) += derivative;
        return;
    }
    if (degree == 4) {
        for (std::size_t i = 0; i < 6; ++i)
            for (std::size_t j = 0; j < 6; ++j)
                for (std::size_t k = 0; k < 6; ++k)
                    for (std::size_t l = 0; l < 6; ++l)
                        if (matches({i, j, k, l, 0}, 4))
                            jet.h4_per_m(i, j, k, l) += derivative;
        return;
    }
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            for (std::size_t k = 0; k < 6; ++k)
                for (std::size_t l = 0; l < 6; ++l)
                    for (std::size_t m = 0; m < 6; ++m)
                        if (matches({i, j, k, l, m}, 5))
                            jet.h5_per_m(i, j, k, l, m) += derivative;
}

double CanonicalPoissonEntry(std::size_t row, std::size_t column) {
    if (row == 0 && column == 1) return 1.0;
    if (row == 1 && column == 0) return -1.0;
    if (row == 2 && column == 3) return 1.0;
    if (row == 3 && column == 2) return -1.0;
    if (row == 4 && column == 5) return -1.0;
    if (row == 5 && column == 4) return 1.0;
    return 0.0;
}

void AddMonomialDerivative(DynamicsJet6& jet, std::size_t output,
                           const std::array<std::size_t, 6>& powers,
                           double coefficient, unsigned maximum_order) {
    std::size_t degree = 0;
    double derivative = coefficient;
    for (std::size_t power : powers) {
        degree += power;
        derivative *= SmallFactorial(power);
    }
    if (degree == 0 || degree > maximum_order || coefficient == 0.0)
        return;
    const auto matches = [&](const std::array<std::size_t, 3>& inputs,
                             std::size_t count) {
        std::array<std::size_t, 6> actual{};
        for (std::size_t index = 0; index < count; ++index)
            ++actual[inputs[index]];
        return actual == powers;
    };
    if (degree == 1) {
        for (std::size_t i = 0; i < 6; ++i)
            if (powers[i] == 1) jet.a_per_m(output, i) += derivative;
        return;
    }
    if (degree == 2) {
        for (std::size_t i = 0; i < 6; ++i)
            for (std::size_t j = 0; j < 6; ++j)
                if (matches({i, j, 0}, 2))
                    jet.f2_per_m(output, i, j) += derivative;
        return;
    }
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            for (std::size_t k = 0; k < 6; ++k)
                if (matches({i, j, k}, 3))
                    jet.f3_per_m(output, i, j, k) += derivative;
}

double Binomial(std::size_t order, std::size_t selected) {
    static constexpr double values[6][6] = {
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 1.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 2.0, 1.0, 0.0, 0.0, 0.0},
        {1.0, 3.0, 3.0, 1.0, 0.0, 0.0},
        {1.0, 4.0, 6.0, 4.0, 1.0, 0.0},
        {1.0, 5.0, 10.0, 10.0, 5.0, 1.0},
    };
    return values[order][selected];
}

void ValidateJet(const DynamicsJet6& jet, unsigned maximum_order,
                 double symmetry_tolerance) {
    RequireFinite(jet.a_per_m, "A");
    if (maximum_order >= 2) {
        RequireFinite(jet.f2_per_m, "F2");
        if (InputSymmetryDefect(jet.f2_per_m) > symmetry_tolerance)
            throw std::invalid_argument(
                "F2 input indices must be symmetric within tolerance");
    }
    if (maximum_order >= 3) {
        RequireFinite(jet.f3_per_m, "F3");
        if (InputSymmetryDefect(jet.f3_per_m) > symmetry_tolerance)
            throw std::invalid_argument(
                "F3 input indices must be symmetric within tolerance");
    }
}

std::vector<StepDescriptor> BuildSteps(
        const std::vector<DynamicsSegment6>& segments,
        const VariationalOptions& options,
        std::vector<std::size_t>& boundary_steps,
        std::vector<double>& boundary_paths) {
    std::vector<StepDescriptor> steps;
    boundary_steps.clear();
    boundary_paths.clear();
    boundary_steps.push_back(0);
    boundary_paths.push_back(0.0);
    double path = 0.0;
    for (std::size_t region = 0; region < segments.size(); ++region) {
        const auto& segment = segments[region];
        if (!std::isfinite(segment.length_m) || segment.length_m <= 0.0)
            throw std::invalid_argument(
                "segment lengths must be finite and positive");
        ValidateJet(segment.jet, options.maximum_order,
                    options.input_symmetry_tolerance);
        const double raw_steps = std::ceil(
            segment.length_m / options.maximum_step_m);
        if (!std::isfinite(raw_steps) || raw_steps < 1.0 ||
            raw_steps > static_cast<double>(options.maximum_steps))
            throw std::invalid_argument("invalid integration step count");
        const std::size_t count = static_cast<std::size_t>(raw_steps);
        if (steps.size() > options.maximum_steps - count)
            throw std::invalid_argument(
                "variational integration exceeds maximum_steps");
        const double step_length = segment.length_m / count;
        for (std::size_t index = 0; index < count; ++index) {
            steps.push_back({region, path + index * step_length, step_length});
        }
        path += segment.length_m;
        boundary_steps.push_back(steps.size());
        boundary_paths.push_back(path);
    }
    return steps;
}

}  // namespace

Matrix6 IdentityMatrix6() {
    Matrix6 result;
    for (std::size_t index = 0; index < 6; ++index)
        result(index, index) = 1.0;
    return result;
}

TaylorMap6 IdentityTaylorMap6() {
    TaylorMap6 result;
    result.r = IdentityMatrix6();
    return result;
}

Matrix6 Multiply(const Matrix6& left, const Matrix6& right) {
    Matrix6 result;
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t k = 0; k < 6; ++k) {
            const double coefficient = left(i, k);
            for (std::size_t j = 0; j < 6; ++j)
                result(i, j) += coefficient * right(k, j);
        }
    return result;
}

TaylorMap6 ComposeTaylorMaps(const TaylorMap6& outer,
                             const TaylorMap6& inner,
                             unsigned maximum_order) {
    if (maximum_order < 1 || maximum_order > 3)
        throw std::invalid_argument("maximum_order must be 1, 2, or 3");
    TaylorMap6 result;
    result.r = Multiply(outer.r, inner.r);
    if (maximum_order >= 2) {
        result.t = LeftMultiply(outer.r, inner.t);
        AddScaled(result.t, TransformInputs(outer.t, inner.r), 1.0);
    }
    if (maximum_order >= 3) {
        result.u = LeftMultiply(outer.r, inner.u);
        AddScaled(result.u,
                  CrossSecondOrder(outer.t, inner.r, inner.t), 1.0);
        AddScaled(result.u, TransformInputs(outer.u, inner.r), 1.0);
    }
    return result;
}

TaylorMap6 IntegrateConstantJet(const DynamicsJet6& jet, double length_m,
                                unsigned maximum_order) {
    if (!std::isfinite(length_m) || length_m <= 0.0)
        throw std::invalid_argument("length_m must be finite and positive");
    if (maximum_order < 1 || maximum_order > 3)
        throw std::invalid_argument("maximum_order must be 1, 2, or 3");
    ValidateJet(jet, maximum_order, 1.0e-12);
    return IntegrateConstantJetUnchecked(jet, length_m, maximum_order);
}

DynamicsJet6 BuildParaxialMagneticDynamicsJet(
        const TransverseMagneticMultipoleExpansion& expansion,
        double magnetic_rigidity_t_m, double curvature_sign,
        double gradient_sign, unsigned maximum_order) {
    if (expansion.order > 4)
        throw std::invalid_argument("multipole order must not exceed four");
    if (maximum_order < 1 || maximum_order > 3)
        throw std::invalid_argument("maximum_order must be 1, 2, or 3");
    if (!std::isfinite(magnetic_rigidity_t_m) ||
        magnetic_rigidity_t_m == 0.0)
        throw std::invalid_argument(
            "magnetic_rigidity_t_m must be finite and nonzero");
    if (!std::isfinite(curvature_sign) || !std::isfinite(gradient_sign))
        throw std::invalid_argument(
            "curvature_sign and gradient_sign must be finite");
    for (double coefficient : expansion.normal_t_per_m_power)
        if (!std::isfinite(coefficient))
            throw std::invalid_argument(
                "normal multipole coefficients must be finite");
    for (double coefficient : expansion.skew_t_per_m_power)
        if (!std::isfinite(coefficient))
            throw std::invalid_argument(
                "skew multipole coefficients must be finite");

    DynamicsJet6 jet;
    const double curvature = curvature_sign *
        expansion.normal_t_per_m_power[0] / magnetic_rigidity_t_m;
    jet.a_per_m(0, 1) = 1.0;
    jet.a_per_m(1, 0) = -curvature * curvature;
    jet.a_per_m(1, 5) = curvature;
    jet.a_per_m(2, 3) = 1.0;
    jet.a_per_m(4, 0) = curvature;

    // x' = px/(1+delta), y' = py/(1+delta), through cubic order.
    AddMonomialDerivative(jet, 0, {0, 1, 0, 0, 0, 1}, -1.0,
                          maximum_order);
    AddMonomialDerivative(jet, 0, {0, 1, 0, 0, 0, 2}, 1.0,
                          maximum_order);
    AddMonomialDerivative(jet, 2, {0, 0, 0, 1, 0, 1}, -1.0,
                          maximum_order);
    AddMonomialDerivative(jet, 2, {0, 0, 0, 1, 0, 2}, 1.0,
                          maximum_order);

    const std::complex<double> imaginary(0.0, 1.0);
    const unsigned field_order = std::min(expansion.order, maximum_order);
    for (unsigned order = 1; order <= field_order; ++order) {
        const std::complex<double> multipole(
            expansion.normal_t_per_m_power[order],
            expansion.skew_t_per_m_power[order]);
        for (unsigned y_power = 0; y_power <= order; ++y_power) {
            const std::complex<double> polynomial =
                multipole * Binomial(order, y_power) *
                std::pow(imaginary, static_cast<int>(y_power));
            for (unsigned delta_power = 0;
                 delta_power + order <= maximum_order; ++delta_power) {
                std::array<std::size_t, 6> powers{};
                powers[0] = order - y_power;
                powers[2] = y_power;
                powers[5] = delta_power;
                const double chromatic = delta_power % 2 == 0 ? 1.0 : -1.0;
                AddMonomialDerivative(
                    jet, 1, powers,
                    -gradient_sign * chromatic * polynomial.real() /
                        magnetic_rigidity_t_m,
                    maximum_order);
                AddMonomialDerivative(
                    jet, 3, powers,
                    gradient_sign * chromatic * polynomial.imag() /
                        magnetic_rigidity_t_m,
                    maximum_order);
            }
        }
    }
    ValidateJet(jet, maximum_order, 1.0e-12);
    return jet;
}

HamiltonianJet6 BuildCanonicalBodyHamiltonianJet(
        const TransverseMagneticMultipoleExpansion& expansion,
        double magnetic_rigidity_t_m, double curvature_sign,
        double gradient_sign, double reference_beta,
        std::optional<double> reference_curvature_per_m) {
    if (expansion.order > 4)
        throw std::invalid_argument("multipole order must not exceed four");
    if (!std::isfinite(magnetic_rigidity_t_m) ||
        magnetic_rigidity_t_m == 0.0)
        throw std::invalid_argument(
            "magnetic_rigidity_t_m must be finite and nonzero");
    if (!std::isfinite(curvature_sign) || !std::isfinite(gradient_sign))
        throw std::invalid_argument(
            "curvature_sign and gradient_sign must be finite");
    if (!std::isfinite(reference_beta) || reference_beta <= 0.0 ||
        reference_beta > 1.0)
        throw std::invalid_argument(
            "reference_beta must be finite and in (0, 1]");
    if (reference_curvature_per_m.has_value() &&
        !std::isfinite(*reference_curvature_per_m))
        throw std::invalid_argument(
            "reference_curvature_per_m must be finite when supplied");
    for (double coefficient : expansion.normal_t_per_m_power)
        if (!std::isfinite(coefficient))
            throw std::invalid_argument(
                "normal multipole coefficients must be finite");
    for (double coefficient : expansion.skew_t_per_m_power)
        if (!std::isfinite(coefficient))
            throw std::invalid_argument(
                "skew multipole coefficients must be finite");

    HamiltonianJet6 result;
    result.reference_beta = reference_beta;
    const double field_curvature = curvature_sign *
        expansion.normal_t_per_m_power[0] / magnetic_rigidity_t_m;
    const double curvature = reference_curvature_per_m.value_or(
        field_curvature);
    result.reference_curvature_per_m = curvature;
    result.field_curvature_per_m = field_curvature;

    // Exact parent:
    // -(1+h*x)*sqrt((1+delta)^2-px^2-py^2) - a_s + H_ref.
    AddHamiltonianMonomial(result, {0, 2, 0, 0, 0, 0}, 0.5);
    AddHamiltonianMonomial(result, {0, 0, 0, 2, 0, 0}, 0.5);
    AddHamiltonianMonomial(result, {2, 0, 0, 0, 0, 0},
                           0.5 * curvature * field_curvature);
    AddHamiltonianMonomial(result, {1, 0, 0, 0, 0, 1}, -curvature);
    AddHamiltonianMonomial(
        result, {0, 0, 0, 0, 0, 2},
        0.5 * (1.0 - reference_beta * reference_beta));

    for (std::size_t momentum : {std::size_t{1}, std::size_t{3}}) {
        std::array<std::size_t, 6> powers{};
        powers[momentum] = 2;
        powers[5] = 1;
        AddHamiltonianMonomial(result, powers, -0.5);
        powers[5] = 0;
        powers[0] = 1;
        AddHamiltonianMonomial(result, powers, 0.5 * curvature);
    }
    AddHamiltonianMonomial(
        result, {0, 0, 0, 0, 0, 3},
        -0.5 * reference_beta * reference_beta *
            (1.0 - reference_beta * reference_beta));

    for (std::size_t momentum : {std::size_t{1}, std::size_t{3}}) {
        std::array<std::size_t, 6> powers{};
        powers[momentum] = 2;
        powers[5] = 2;
        AddHamiltonianMonomial(result, powers, 0.5);
        powers[5] = 1;
        powers[0] = 1;
        AddHamiltonianMonomial(result, powers, -0.5 * curvature);
    }
    AddHamiltonianMonomial(result, {0, 4, 0, 0, 0, 0}, 0.125);
    AddHamiltonianMonomial(result, {0, 0, 0, 4, 0, 0}, 0.125);
    AddHamiltonianMonomial(result, {0, 2, 0, 2, 0, 0}, 0.25);
    const double beta2 = reference_beta * reference_beta;
    AddHamiltonianMonomial(
        result, {0, 0, 0, 0, 0, 4},
        beta2 * (1.0 - beta2) * (5.0 * beta2 - 1.0) / 8.0);

    // Fifth-degree kinematic terms from the same exact square roots.
    for (std::size_t momentum : {std::size_t{1}, std::size_t{3}}) {
        std::array<std::size_t, 6> powers{};
        powers[momentum] = 2;
        powers[5] = 3;
        AddHamiltonianMonomial(result, powers, -0.5);
        powers[5] = 2;
        powers[0] = 1;
        AddHamiltonianMonomial(result, powers, 0.5 * curvature);

        powers = {};
        powers[momentum] = 4;
        powers[5] = 1;
        AddHamiltonianMonomial(result, powers, -3.0 / 8.0);
        powers[5] = 0;
        powers[0] = 1;
        AddHamiltonianMonomial(result, powers, curvature / 8.0);
    }
    AddHamiltonianMonomial(result, {0, 2, 0, 2, 0, 1}, -3.0 / 4.0);
    AddHamiltonianMonomial(
        result, {1, 2, 0, 2, 0, 0}, curvature / 4.0);
    AddHamiltonianMonomial(
        result, {0, 0, 0, 0, 0, 5},
        beta2 * beta2 * (1.0 - beta2) * (3.0 - 7.0 * beta2) / 8.0);

    const std::complex<double> imaginary(0.0, 1.0);
    for (unsigned order = 1; order <= expansion.order; ++order) {
        const unsigned degree = order + 1;
        const std::complex<double> normalized = gradient_sign *
            std::complex<double>(
                expansion.normal_t_per_m_power[order],
                expansion.skew_t_per_m_power[order]) /
            magnetic_rigidity_t_m;
        for (unsigned y_power = 0; y_power <= degree; ++y_power) {
            std::array<std::size_t, 6> powers{};
            powers[0] = degree - y_power;
            powers[2] = y_power;
            const std::complex<double> polynomial = normalized *
                Binomial(degree, y_power) *
                std::pow(imaginary, static_cast<int>(y_power)) /
                static_cast<double>(degree);
            AddHamiltonianMonomial(result, powers, polynomial.real());
        }
    }

    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t a = 0; a < 6; ++a) {
            const double poisson = CanonicalPoissonEntry(i, a);
            for (std::size_t j = 0; j < 6; ++j) {
                result.dynamics.a_per_m(i, j) +=
                    poisson * result.h2_per_m(a, j);
                for (std::size_t k = 0; k < 6; ++k) {
                    result.dynamics.f2_per_m(i, j, k) +=
                        poisson * result.h3_per_m(a, j, k);
                    for (std::size_t l = 0; l < 6; ++l) {
                        result.dynamics.f3_per_m(i, j, k, l) +=
                            poisson * result.h4_per_m(a, j, k, l);
                        for (std::size_t m = 0; m < 6; ++m)
                            result.dynamics.f4_per_m(i, j, k, l, m) +=
                                poisson * result.h5_per_m(a, j, k, l, m);
                    }
                }
            }
        }
    ValidateJet(result.dynamics, 3, 1.0e-12);
    return result;
}

VariationalReport6 PropagateVariationalMap(
        const std::vector<DynamicsSegment6>& segments,
        const VariationalOptions& options) {
    if (segments.empty())
        throw std::invalid_argument("at least one dynamics segment is required");
    if (options.maximum_order < 1 || options.maximum_order > 3)
        throw std::invalid_argument("maximum_order must be 1, 2, or 3");
    if (!std::isfinite(options.maximum_step_m) ||
        options.maximum_step_m <= 0.0)
        throw std::invalid_argument(
            "maximum_step_m must be finite and positive");
    if (options.maximum_steps == 0 || options.maximum_region_pairs == 0)
        throw std::invalid_argument(
            "maximum_steps and maximum_region_pairs must be positive");
    if (!std::isfinite(options.input_symmetry_tolerance) ||
        options.input_symmetry_tolerance < 0.0)
        throw std::invalid_argument(
            "input_symmetry_tolerance must be finite and nonnegative");

    std::vector<std::size_t> boundary_steps;
    std::vector<double> boundary_paths;
    const std::vector<StepDescriptor> steps = BuildSteps(
        segments, options, boundary_steps, boundary_paths);
    const std::size_t step_count = steps.size();

    std::vector<Matrix6> prefix_r(step_count + 1);
    std::vector<Matrix6> local_r(step_count);
    prefix_r[0] = IdentityMatrix6();
    TaylorMap6 endpoint = IdentityTaylorMap6();
    std::vector<TaylorMap6> boundary_maps;
    boundary_maps.reserve(segments.size() + 1);
    boundary_maps.push_back(endpoint);
    std::size_t next_boundary = 1;
    for (std::size_t step = 0; step < step_count; ++step) {
        const auto& descriptor = steps[step];
        const TaylorMap6 local = IntegrateConstantJetUnchecked(
            segments[descriptor.region].jet, descriptor.length_m,
            options.maximum_order);
        local_r[step] = local.r;
        endpoint = ComposeTaylorMaps(local, endpoint, options.maximum_order);
        prefix_r[step + 1] = endpoint.r;
        if (next_boundary < boundary_steps.size() &&
            step + 1 == boundary_steps[next_boundary]) {
            boundary_maps.push_back(endpoint);
            ++next_boundary;
        }
    }

    std::vector<Matrix6> suffix_r(step_count + 1);
    suffix_r[step_count] = IdentityMatrix6();
    for (std::size_t reverse = step_count; reverse > 0; --reverse) {
        const std::size_t step = reverse - 1;
        suffix_r[step] = Multiply(suffix_r[step + 1], local_r[step]);
    }

    VariationalReport6 report;
    report.maximum_order = options.maximum_order;
    report.endpoint_map = endpoint;
    report.diagnostics.integration_steps = step_count;
    report.diagnostics.t_input_symmetry_defect =
        InputSymmetryDefect(endpoint.t);
    report.diagnostics.u_input_symmetry_defect =
        InputSymmetryDefect(endpoint.u);

    report.stations.reserve(segments.size() + 1);
    for (std::size_t boundary = 0; boundary < boundary_steps.size(); ++boundary) {
        TransferStation6 station;
        station.path_length_m = boundary_paths[boundary];
        station.boundary_index = boundary;
        station.map_from_start = boundary_maps[boundary];
        station.r_to_end = suffix_r[boundary_steps[boundary]];
        report.stations.push_back(std::move(station));
        const Matrix6 recomposed = Multiply(
            report.stations.back().r_to_end,
            report.stations.back().map_from_start.r);
        report.diagnostics.r_composition_error = std::max(
            report.diagnostics.r_composition_error,
            MaximumAbsoluteDifference(recomposed, endpoint.r));
    }

    report.regions.resize(segments.size());
    double path = 0.0;
    for (std::size_t region = 0; region < segments.size(); ++region) {
        auto& output = report.regions[region];
        output.region_index = region;
        output.name = segments[region].name.empty()
            ? "segment_" + std::to_string(region)
            : segments[region].name;
        output.s_begin_m = path;
        path += segments[region].length_m;
        output.s_end_m = path;
    }

    if (options.maximum_order >= 2) {
        std::vector<Tensor3Map6> source_t_at_station(segments.size());
        std::map<std::pair<std::size_t, std::size_t>, Tensor4Map6>
            pair_contributions;

        for (std::size_t step = 0; step < step_count; ++step) {
            const auto& descriptor = steps[step];
            const std::size_t downstream = descriptor.region;
            const DynamicsJet6& jet = segments[downstream].jet;
            const TaylorMap6 local = IntegrateConstantJetUnchecked(
                jet, descriptor.length_m, options.maximum_order);
            const Matrix6& before = prefix_r[step];
            const Matrix6& after_to_end = suffix_r[step + 1];

            const Tensor3Map6 local_t_at_station =
                TransformInputs(local.t, before);
            AddScaled(report.regions[downstream].t_at_end,
                      LeftMultiply(after_to_end, local_t_at_station), 1.0);

            if (options.maximum_order >= 3) {
                DynamicsJet6 direct_jet = jet;
                direct_jet.f2_per_m = Tensor3Map6{};
                const TaylorMap6 direct = IntegrateConstantJetUnchecked(
                    direct_jet, descriptor.length_m, 3);
                const Tensor4Map6 local_cascade =
                    Difference(local.u, direct.u);
                AddScaled(report.regions[downstream].u_direct_at_end,
                          LeftMultiply(after_to_end,
                              TransformInputs(direct.u, before)), 1.0);
                AddScaled(report.regions[downstream].u_local_cascade_at_end,
                          LeftMultiply(after_to_end,
                              TransformInputs(local_cascade, before)), 1.0);

                if (MaximumAbsoluteEntry(local.t) > 0.0) {
                    for (std::size_t upstream = 0;
                         upstream < source_t_at_station.size(); ++upstream) {
                        if (MaximumAbsoluteEntry(
                                source_t_at_station[upstream]) == 0.0)
                            continue;
                        const Tensor4Map6 pair_at_output = CrossSecondOrder(
                            local.t, before, source_t_at_station[upstream]);
                        const Tensor4Map6 pair_at_end = LeftMultiply(
                            after_to_end, pair_at_output);
                        if (MaximumAbsoluteEntry(pair_at_end) == 0.0)
                            continue;
                        if (upstream == downstream) {
                            AddScaled(report.regions[downstream]
                                          .u_local_cascade_at_end,
                                      pair_at_end, 1.0);
                            continue;
                        }
                        const auto key = std::make_pair(upstream, downstream);
                        auto found = pair_contributions.find(key);
                        if (found == pair_contributions.end()) {
                            if (pair_contributions.size() >=
                                options.maximum_region_pairs)
                                throw std::runtime_error(
                                    "nonlinear attribution exceeds "
                                    "maximum_region_pairs");
                            found = pair_contributions.emplace(
                                key, Tensor4Map6{}).first;
                        }
                        AddScaled(found->second, pair_at_end, 1.0);
                    }
                }
            }

            for (auto& source : source_t_at_station)
                source = LeftMultiply(local.r, source);
            AddScaled(source_t_at_station[downstream],
                      local_t_at_station, 1.0);
        }

        report.region_pairs.reserve(pair_contributions.size());
        for (const auto& item : pair_contributions) {
            RegionPairNonlinearContribution6 output;
            output.upstream_region = item.first.first;
            output.downstream_region = item.first.second;
            output.u_cascade_at_end = item.second;
            output.maximum_absolute_entry =
                MaximumAbsoluteEntry(item.second);
            report.region_pairs.push_back(std::move(output));
        }

        Tensor3Map6 reconstructed_t;
        Tensor4Map6 reconstructed_u;
        for (const auto& region : report.regions) {
            AddScaled(reconstructed_t, region.t_at_end, 1.0);
            AddScaled(reconstructed_u, region.u_direct_at_end, 1.0);
            AddScaled(reconstructed_u,
                      region.u_local_cascade_at_end, 1.0);
        }
        for (const auto& pair : report.region_pairs)
            AddScaled(reconstructed_u, pair.u_cascade_at_end, 1.0);
        report.diagnostics.t_reconstruction_error =
            MaximumAbsoluteDifference(reconstructed_t, endpoint.t);
        if (options.maximum_order >= 3)
            report.diagnostics.u_reconstruction_error =
                MaximumAbsoluteDifference(reconstructed_u, endpoint.u);
    }

    return report;
}

double MaximumAbsoluteEntry(const Matrix6& value) {
    double result = 0.0;
    for (double item : value.values) result = std::max(result, std::abs(item));
    return result;
}

double MaximumAbsoluteEntry(const Tensor3Map6& value) {
    double result = 0.0;
    for (double item : value.values) result = std::max(result, std::abs(item));
    return result;
}

double MaximumAbsoluteEntry(const Tensor4Map6& value) {
    double result = 0.0;
    for (double item : value.values) result = std::max(result, std::abs(item));
    return result;
}

double MaximumAbsoluteDifference(const Matrix6& left, const Matrix6& right) {
    return MaximumAbsoluteEntry(Difference(left, right));
}

double MaximumAbsoluteDifference(const Tensor3Map6& left,
                                 const Tensor3Map6& right) {
    return MaximumAbsoluteEntry(Difference(left, right));
}

double MaximumAbsoluteDifference(const Tensor4Map6& left,
                                 const Tensor4Map6& right) {
    return MaximumAbsoluteEntry(Difference(left, right));
}

double InputSymmetryDefect(const Tensor3Map6& value) {
    double result = 0.0;
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            for (std::size_t k = 0; k < 6; ++k)
                result = std::max(
                    result, std::abs(value(i, j, k) - value(i, k, j)));
    return result;
}

double InputSymmetryDefect(const Tensor4Map6& value) {
    double result = 0.0;
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            for (std::size_t k = 0; k < 6; ++k)
                for (std::size_t l = 0; l < 6; ++l) {
                    result = std::max(result,
                        std::abs(value(i, j, k, l) -
                                 value(i, k, j, l)));
                    result = std::max(result,
                        std::abs(value(i, j, k, l) -
                                 value(i, j, l, k)));
                }
    return result;
}

}  // namespace radia::beam
