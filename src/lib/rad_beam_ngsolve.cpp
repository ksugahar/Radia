#include "rad_beam_ngsolve.h"

#include <core/taskmanager.hpp>
#include <elementtransformation.hpp>
#include <meshaccess.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace radia::beam {
namespace {

constexpr std::size_t kSampleCount = 9;

bool Finite(const Vector3& value) {
    return std::all_of(value.begin(), value.end(), [](double item) {
        return std::isfinite(item);
    });
}

double Dot(const Vector3& left, const Vector3& right) {
    return left[0] * right[0] + left[1] * right[1] +
           left[2] * right[2];
}

Vector3 Scale(const Vector3& value, double factor) {
    return {factor * value[0], factor * value[1], factor * value[2]};
}

Vector3 Add(const Vector3& left, const Vector3& right) {
    return {left[0] + right[0], left[1] + right[1],
            left[2] + right[2]};
}

Vector3 Subtract(const Vector3& left, const Vector3& right) {
    return {left[0] - right[0], left[1] - right[1],
            left[2] - right[2]};
}

Vector3 Cross(const Vector3& left, const Vector3& right) {
    return {
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    };
}

double Norm(const Vector3& value) {
    return std::sqrt(Dot(value, value));
}

Vector3 Normalize(const Vector3& value, const char* name) {
    const double norm = Norm(value);
    if (!std::isfinite(norm) || norm <= 64.0 *
            std::numeric_limits<double>::epsilon())
        throw std::invalid_argument(std::string(name) +
                                    " must have nonzero finite norm");
    return Scale(value, 1.0 / norm);
}

Vector3 ProjectNormal(const Vector3& value, const Vector3& tangent,
                      const char* name) {
    return Normalize(Subtract(value, Scale(tangent, Dot(value, tangent))),
                     name);
}

Vector3 TransportHorizontal(const Vector3& preferred,
                            const Vector3& initial,
                            const Vector3& tangent) {
    const auto projected_norm = [&](const Vector3& candidate) {
        return Norm(Subtract(
            candidate, Scale(tangent, Dot(candidate, tangent))));
    };
    if (projected_norm(preferred) > 64.0 *
            std::numeric_limits<double>::epsilon())
        return ProjectNormal(preferred, tangent,
                             "transported horizontal axis");
    if (projected_norm(initial) > 64.0 *
            std::numeric_limits<double>::epsilon())
        return ProjectNormal(initial, tangent,
                             "initial horizontal axis");

    const std::array<Vector3, 3> coordinate_axes{{
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
    }};
    const auto least_aligned = std::min_element(
        coordinate_axes.begin(), coordinate_axes.end(),
        [&](const Vector3& left, const Vector3& right) {
            return std::abs(Dot(left, tangent)) <
                   std::abs(Dot(right, tangent));
        });
    return ProjectNormal(*least_aligned, tangent,
                         "fallback horizontal axis");
}

Vector3 EvaluateField(
        const ngcomp::GridFunctionCoefficientFunction& coefficient,
        const std::shared_ptr<ngcomp::MeshAccess>& mesh,
        const Vector3& point, ngstd::LocalHeap& local_heap,
        std::size_t segment, std::size_t sample) {
    local_heap.CleanUp();
    std::array<double, 3> physical_storage = point;
    ngbla::FlatVector<double> physical_point(
        3, physical_storage.data());
    ngfem::IntegrationPoint integration_point;
    const ngfem::ElementId element = mesh->FindElementOfPoint(
        physical_point, integration_point, true);
    if (element.IsInvalid() || !element.IsVolume())
        throw std::runtime_error(
            "beam GridFunction sample is outside the volume mesh at segment " +
            std::to_string(segment) + ", sample " +
            std::to_string(sample));
    auto& transformation = mesh->GetTrafo(element, local_heap);
    ngfem::MappedIntegrationPoint<3, 3> mapped_point(
        integration_point, transformation);
    ngbla::FlatVector<double> result(3, local_heap);
    coefficient.Evaluate(mapped_point, result);
    Vector3 output{result[0], result[1], result[2]};
    if (!Finite(output))
        throw std::runtime_error(
            "beam GridFunction evaluation returned a non-finite field at "
            "segment " + std::to_string(segment));
    return output;
}

Vector3 LocalComponents(const Vector3& field, const Vector3& horizontal,
                        const Vector3& vertical, const Vector3& tangent) {
    return {Dot(field, horizontal), Dot(field, vertical),
            Dot(field, tangent)};
}

void BuildGenerator(GridFunctionSegmentLinearization& output,
                    double magnetic_rigidity_t_m,
                    double curvature_sign, double gradient_sign) {
    const double d_bx_dx = output.field_gradient_local_t_per_m[0];
    const double d_bx_dy = output.field_gradient_local_t_per_m[1];
    const double d_by_dx = output.field_gradient_local_t_per_m[2];
    const double d_by_dy = output.field_gradient_local_t_per_m[3];
    const double normal_gradient_t_per_m =
        0.5 * (d_by_dx + d_bx_dy);
    const double skew_gradient_t_per_m =
        0.5 * (d_bx_dx - d_by_dy);

    output.curvature_per_m = curvature_sign *
        output.center_field_local_t[1] / magnetic_rigidity_t_m;
    output.normal_gradient_per_m2 = gradient_sign *
        normal_gradient_t_per_m / magnetic_rigidity_t_m;
    output.skew_gradient_per_m2 = gradient_sign *
        skew_gradient_t_per_m / magnetic_rigidity_t_m;
    output.transverse_divergence_t_per_m = d_bx_dx + d_by_dy;
    output.transverse_curl_mismatch_t_per_m =
        0.5 * (d_by_dx - d_bx_dy);

    const double h = output.curvature_per_m;
    const double k1 = output.normal_gradient_per_m2;
    const double k1s = output.skew_gradient_per_m2;
    Matrix6& a = output.a_per_m;
    a(0, 1) = 1.0;
    a(1, 0) = -(h * h + k1);
    a(1, 2) = k1s;
    a(1, 5) = h;
    a(2, 3) = 1.0;
    a(3, 0) = k1s;
    a(3, 2) = k1;
    a(4, 0) = h;
}

GridFunctionSegmentLinearization LinearizeSegment(
        const ngcomp::GridFunctionCoefficientFunction& coefficient,
        const std::shared_ptr<ngcomp::MeshAccess>& mesh,
        const Vector3& position, const Vector3& tangent,
        const Vector3& horizontal, const Vector3& vertical,
        const GridFunctionLinearizationOptions& options,
        ngstd::LocalHeap& local_heap, std::size_t segment) {
    GridFunctionSegmentLinearization output;
    output.reference_position_m = position;
    output.horizontal_axis = horizontal;
    output.vertical_axis = vertical;
    output.tangent_axis = tangent;

    const double radius = options.sample_radius_m;
    const double diagonal = radius / std::sqrt(2.0);
    const std::array<std::array<double, 2>, kSampleCount> offsets{{
        {0.0, 0.0},
        {radius, 0.0}, {-radius, 0.0},
        {0.0, radius}, {0.0, -radius},
        {diagonal, diagonal}, {diagonal, -diagonal},
        {-diagonal, diagonal}, {-diagonal, -diagonal},
    }};

    std::array<Vector3, kSampleCount> local_fields{};
    for (std::size_t sample = 0; sample < kSampleCount; ++sample) {
        const Vector3 sample_position = Add(
            position, Add(Scale(horizontal, offsets[sample][0]),
                          Scale(vertical, offsets[sample][1])));
        const Vector3 global_field = EvaluateField(
            coefficient, mesh, sample_position, local_heap, segment, sample);
        local_fields[sample] = LocalComponents(
            global_field, horizontal, vertical, tangent);
        for (std::size_t component = 0; component < 3; ++component) {
            output.sample_positions_m[3 * sample + component] =
                sample_position[component];
            output.sample_fields_local_t[3 * sample + component] =
                local_fields[sample][component];
        }
    }
    output.center_field_local_t = local_fields[0];

    constexpr double sample_count = static_cast<double>(kSampleCount);
    const double slope_denominator = 4.0 * radius * radius;
    for (std::size_t component = 0; component < 3; ++component) {
        double mean = 0.0;
        double x_moment = 0.0;
        double y_moment = 0.0;
        for (std::size_t sample = 0; sample < kSampleCount; ++sample) {
            mean += local_fields[sample][component];
            x_moment += offsets[sample][0] * local_fields[sample][component];
            y_moment += offsets[sample][1] * local_fields[sample][component];
        }
        output.fitted_center_field_local_t[component] = mean / sample_count;
        output.field_gradient_local_t_per_m[2 * component] =
            x_moment / slope_denominator;
        output.field_gradient_local_t_per_m[2 * component + 1] =
            y_moment / slope_denominator;
    }

    double residual_squared = 0.0;
    double maximum_residual = 0.0;
    for (std::size_t sample = 0; sample < kSampleCount; ++sample) {
        double sample_residual_squared = 0.0;
        for (std::size_t component = 0; component < 3; ++component) {
            const double predicted =
                output.fitted_center_field_local_t[component] +
                output.field_gradient_local_t_per_m[2 * component] *
                    offsets[sample][0] +
                output.field_gradient_local_t_per_m[2 * component + 1] *
                    offsets[sample][1];
            const double residual = local_fields[sample][component] - predicted;
            sample_residual_squared += residual * residual;
        }
        residual_squared += sample_residual_squared;
        maximum_residual = std::max(
            maximum_residual, std::sqrt(sample_residual_squared));
    }
    output.rms_fit_residual_t =
        std::sqrt(residual_squared / sample_count);
    output.maximum_fit_residual_t = maximum_residual;
    output.center_fit_bias_t = Norm(Subtract(
        output.fitted_center_field_local_t,
        output.center_field_local_t));
    BuildGenerator(output, options.magnetic_rigidity_t_m,
                   options.curvature_sign, options.gradient_sign);
    return output;
}

}  // namespace

GridFunctionTransferReport6 PropagateGridFunctionLinearMap(
        const std::shared_ptr<ngcomp::GridFunction>& field,
        const std::vector<double>& segment_lengths_m,
        const std::vector<Vector3>& reference_positions_m,
        const std::vector<Vector3>& reference_tangents,
        const std::vector<std::string>& names,
        const GridFunctionLinearizationOptions& options) {
    if (!field)
        throw std::invalid_argument("field GridFunction must not be null");
    const std::size_t count = segment_lengths_m.size();
    if (count == 0 || reference_positions_m.size() != count ||
        reference_tangents.size() != count)
        throw std::invalid_argument(
            "segment lengths, reference positions, and tangents must have "
            "the same nonzero length");
    if (!names.empty() && names.size() != count)
        throw std::invalid_argument(
            "names must contain one entry per reference segment");
    if (!std::isfinite(options.magnetic_rigidity_t_m) ||
        options.magnetic_rigidity_t_m == 0.0)
        throw std::invalid_argument(
            "magnetic_rigidity_t_m must be finite and nonzero");
    if (!std::isfinite(options.sample_radius_m) ||
        options.sample_radius_m <= 0.0)
        throw std::invalid_argument(
            "sample_radius_m must be finite and positive");
    if (!Finite(options.initial_horizontal) ||
        !std::isfinite(options.curvature_sign) ||
        !std::isfinite(options.gradient_sign) ||
        !std::isfinite(options.maximum_step_m) ||
        options.maximum_step_m <= 0.0 || options.maximum_steps == 0)
        throw std::invalid_argument(
            "GridFunction linearization options must be finite and valid");
    for (std::size_t index = 0; index < count; ++index) {
        if (!std::isfinite(segment_lengths_m[index]) ||
            segment_lengths_m[index] <= 0.0 ||
            !Finite(reference_positions_m[index]) ||
            !Finite(reference_tangents[index]))
            throw std::invalid_argument(
                "reference segment data must be finite and lengths positive");
    }

    const auto mesh = field->GetMeshAccess();
    if (!mesh || mesh->GetDimension() != 3)
        throw std::invalid_argument(
            "beam GridFunction must belong to a three-dimensional mesh");
    ngcomp::GridFunctionCoefficientFunction coefficient(field);
    if (coefficient.IsComplex())
        throw std::invalid_argument(
            "beam GridFunction must be real-valued");
    if (coefficient.Dimension() != 3)
        throw std::invalid_argument(
            "beam GridFunction must have exactly three field components");

    GridFunctionTransferReport6 output;
    output.magnetic_rigidity_t_m = options.magnetic_rigidity_t_m;
    output.sample_radius_m = options.sample_radius_m;
    output.linearizations.reserve(count);
    std::vector<DynamicsSegment6> dynamics(count);

    ngcore::RegionTaskManager task_manager;
    ngstd::LocalHeap local_heap(1 << 20, "radia_beam_gridfunction");
    Vector3 previous_horizontal{};
    for (std::size_t index = 0; index < count; ++index) {
        const Vector3 tangent = Normalize(
            reference_tangents[index], "reference tangent");
        const Vector3 horizontal = TransportHorizontal(
            index == 0 ? options.initial_horizontal : previous_horizontal,
            options.initial_horizontal, tangent);
        const Vector3 vertical = Normalize(
            Cross(tangent, horizontal), "vertical axis");
        previous_horizontal = horizontal;
        output.linearizations.push_back(LinearizeSegment(
            coefficient, mesh, reference_positions_m[index], tangent,
            horizontal, vertical, options, local_heap, index));

        dynamics[index].length_m = segment_lengths_m[index];
        dynamics[index].jet.a_per_m = output.linearizations.back().a_per_m;
        dynamics[index].name = names.empty()
            ? "segment_" + std::to_string(index + 1)
            : names[index];
    }

    VariationalOptions variational;
    variational.maximum_order = 1;
    variational.maximum_step_m = options.maximum_step_m;
    variational.maximum_steps = options.maximum_steps;
    output.transfer = PropagateVariationalMap(dynamics, variational);
    return output;
}

}  // namespace radia::beam
