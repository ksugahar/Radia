#include "rad_beam_transfer.h"

#include <algorithm>
#include <cmath>
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
