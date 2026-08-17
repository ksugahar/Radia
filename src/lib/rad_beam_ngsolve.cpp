#include "rad_beam_ngsolve.h"

#include <core/taskmanager.hpp>
#include <elementtransformation.hpp>
#include <hcurlhofespace.hpp>
#include <meshaccess.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
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

Vector3 SeedHorizontal(const Vector3& initial, const Vector3& tangent) {
    const auto projected_norm = [&](const Vector3& candidate) {
        return Norm(Subtract(
            candidate, Scale(tangent, Dot(candidate, tangent))));
    };
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

Vector3 ReflectInPlane(const Vector3& value, const Vector3& plane_normal,
                       double normal_squared) {
    return Subtract(value, Scale(
        plane_normal, 2.0 * Dot(plane_normal, value) / normal_squared));
}

Vector3 TransportHorizontalDoubleReflection(
        const Vector3& previous_position, const Vector3& position,
        const Vector3& previous_tangent, const Vector3& tangent,
        const Vector3& previous_horizontal) {
    // Wang, Juttler, Zheng, and Liu, ACM TOG 27(1), 2008,
    // DOI 10.1145/1330511.1330513.  The two Householder reflections give a
    // fourth-order discrete approximation of the Bishop/RMF initial-value
    // problem without differentiating the tracked orbit.  In particular it
    // remains defined through zero-curvature stations where Frenet framing
    // has no normal.
    const Vector3 chord = Subtract(position, previous_position);
    const double chord_squared = Dot(chord, chord);
    const double position_scale = std::max({
        1.0, Norm(previous_position), Norm(position)});
    const double chord_tolerance = 64.0 *
        std::numeric_limits<double>::epsilon() * position_scale;
    if (chord_squared <= chord_tolerance * chord_tolerance)
        throw std::invalid_argument(
            "consecutive reference positions must be distinct for the "
            "Bishop double-reflection frame");

    const Vector3 reflected_horizontal = ReflectInPlane(
        previous_horizontal, chord, chord_squared);
    const Vector3 reflected_tangent = ReflectInPlane(
        previous_tangent, chord, chord_squared);
    const Vector3 tangent_difference = Subtract(tangent, reflected_tangent);
    const double difference_squared = Dot(
        tangent_difference, tangent_difference);
    Vector3 horizontal = reflected_horizontal;
    const double tangent_tolerance = 64.0 *
        std::numeric_limits<double>::epsilon();
    if (difference_squared > tangent_tolerance * tangent_tolerance)
        horizontal = ReflectInPlane(
            reflected_horizontal, tangent_difference, difference_squared);

    // Reflections preserve the ideal constraints.  Reproject only to remove
    // floating-point drift accumulated along long RK trajectories.
    return ProjectNormal(horizontal, tangent,
                         "Bishop transported horizontal axis");
}

Vector3 RotateNormalAboutTangent(const Vector3& normal,
                                 const Vector3& tangent, double angle) {
    return ProjectNormal(Add(
        Scale(normal, std::cos(angle)),
        Scale(Cross(tangent, normal), std::sin(angle))),
        tangent, "minimal-twist horizontal axis");
}

std::shared_ptr<ngcomp::GridFunctionCoefficientFunction>
MagneticCoefficient(
        const std::shared_ptr<ngcomp::GridFunction>& field,
        GridFunctionMagneticInput magnetic_input) {
    if (!field)
        throw std::invalid_argument("field GridFunction must not be null");
    const auto space = field->GetFESpace();
    if (!space)
        throw std::invalid_argument(
            "beam GridFunction must own an NGSolve FESpace");

    std::shared_ptr<ngcomp::GridFunctionCoefficientFunction> coefficient;
    if (magnetic_input ==
            GridFunctionMagneticInput::HCurlVectorPotential) {
        if (!std::dynamic_pointer_cast<ngcomp::HCurlHighOrderFESpace>(space))
            throw std::invalid_argument(
                "HCurl vector-potential input requires a GridFunction on "
                "an NGSolve HCurl space");
        coefficient = field->Deriv();
    } else {
        coefficient = std::make_shared<
            ngcomp::GridFunctionCoefficientFunction>(field);
    }
    if (!coefficient || coefficient->IsComplex())
        throw std::invalid_argument(
            "beam GridFunction magnetic input must be real-valued");
    if (coefficient->Dimension() != 3)
        throw std::invalid_argument(
            "beam GridFunction magnetic input must produce exactly three "
            "components");
    return coefficient;
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

void BuildDynamicsJet(GridFunctionSegmentLinearization& output,
                      const GridFunctionLinearizationOptions& options) {
    const double d_bx_dx = output.field_gradient_local_t_per_m[0];
    const double d_bx_dy = output.field_gradient_local_t_per_m[1];
    const double d_by_dx = output.field_gradient_local_t_per_m[2];
    const double d_by_dy = output.field_gradient_local_t_per_m[3];
    const double normal_gradient_t_per_m =
        0.5 * (d_by_dx + d_bx_dy);
    const double skew_gradient_t_per_m =
        0.5 * (d_bx_dx - d_by_dy);

    output.curvature_per_m = options.curvature_sign *
        output.center_field_local_t[1] /
        options.magnetic_rigidity_t_m;
    output.normal_gradient_per_m2 = options.gradient_sign *
        output.multipoles.normal_t_per_m_power[1] /
        options.magnetic_rigidity_t_m;
    output.skew_gradient_per_m2 = options.gradient_sign *
        output.multipoles.skew_t_per_m_power[1] /
        options.magnetic_rigidity_t_m;
    output.transverse_divergence_t_per_m = d_bx_dx + d_by_dy;
    output.transverse_curl_mismatch_t_per_m =
        0.5 * (d_by_dx - d_bx_dy);

    // Keep the affine Maxwell diagnostics independently inspectable. The map
    // itself is generated from the harmonic multipole fit below.
    (void)normal_gradient_t_per_m;
    (void)skew_gradient_t_per_m;
    output.dynamics_jet = BuildParaxialMagneticDynamicsJet(
        output.multipoles, options.magnetic_rigidity_t_m,
        options.curvature_sign, options.gradient_sign,
        options.maximum_map_order);
    output.a_per_m = output.dynamics_jet.a_per_m;
}

void FitTransverseMultipoles(
        GridFunctionSegmentLinearization& output,
        const std::array<std::array<double, 2>, kSampleCount>& offsets,
        const std::array<Vector3, kSampleCount>& local_fields,
        unsigned order, double radius) {
    output.multipoles.order = order;
    output.multipoles.normal_t_per_m_power[0] = local_fields[0][1];
    output.multipoles.skew_t_per_m_power[0] = local_fields[0][0];
    for (unsigned degree = 1; degree <= order; ++degree) {
        std::complex<double> moment{};
        for (std::size_t sample = 1; sample < kSampleCount; ++sample) {
            const std::complex<double> normalized(
                offsets[sample][0] / radius,
                offsets[sample][1] / radius);
            const std::complex<double> field(
                local_fields[sample][1], local_fields[sample][0]);
            moment += field * std::conj(std::pow(
                normalized, static_cast<int>(degree)));
        }
        const std::complex<double> coefficient =
            moment / (8.0 * std::pow(radius, static_cast<int>(degree)));
        output.multipoles.normal_t_per_m_power[degree] =
            coefficient.real();
        output.multipoles.skew_t_per_m_power[degree] =
            coefficient.imag();
    }

    double residual_squared = 0.0;
    double maximum_residual = 0.0;
    for (std::size_t sample = 0; sample < kSampleCount; ++sample) {
        const std::complex<double> coordinate(
            offsets[sample][0], offsets[sample][1]);
        std::complex<double> predicted{};
        for (unsigned degree = 0; degree <= order; ++degree) {
            const std::complex<double> coefficient(
                output.multipoles.normal_t_per_m_power[degree],
                output.multipoles.skew_t_per_m_power[degree]);
            predicted += coefficient * std::pow(
                coordinate, static_cast<int>(degree));
        }
        const std::complex<double> actual(
            local_fields[sample][1], local_fields[sample][0]);
        const double residual = std::abs(actual - predicted);
        residual_squared += residual * residual;
        maximum_residual = std::max(maximum_residual, residual);
    }
    output.multipole_rms_fit_residual_t = std::sqrt(
        residual_squared / static_cast<double>(kSampleCount));
    output.multipole_maximum_fit_residual_t = maximum_residual;
    output.multipole_fit_rank = order + 1;
    output.multipole_scaled_design_condition = 1.0;
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
    FitTransverseMultipoles(output, offsets, local_fields,
                            options.multipole_order, radius);
    BuildDynamicsJet(output, options);
    return output;
}

}  // namespace

GridFunctionTransferReport6 PropagateGridFunctionMultipoleMap(
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
    if (options.multipole_order < 1 || options.multipole_order > 4 ||
        options.maximum_map_order < 1 || options.maximum_map_order > 3 ||
        !Finite(options.initial_horizontal) ||
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
    auto coefficient = MagneticCoefficient(field, options.magnetic_input);

    GridFunctionTransferReport6 output;
    output.magnetic_input = options.magnetic_input;
    output.grid_function_space_class = field->GetFESpace()->GetClassName();
    output.grid_function_space_order = field->GetFESpace()->GetOrder();
    output.magnetic_rigidity_t_m = options.magnetic_rigidity_t_m;
    output.sample_radius_m = options.sample_radius_m;
    output.multipole_order = options.multipole_order;
    output.periodic_frame = options.periodic_frame;
    output.linearizations.reserve(count);
    std::vector<DynamicsSegment6> dynamics(count);

    std::vector<Vector3> tangents(count);
    std::vector<Vector3> horizontal(count);
    std::vector<Vector3> vertical(count);
    for (std::size_t index = 0; index < count; ++index) {
        tangents[index] = Normalize(
            reference_tangents[index], "reference tangent");
        horizontal[index] = index == 0
            ? SeedHorizontal(options.initial_horizontal, tangents[index])
            : TransportHorizontalDoubleReflection(
                reference_positions_m[index - 1],
                reference_positions_m[index], tangents[index - 1],
                tangents[index], horizontal[index - 1]);
    }

    if (options.periodic_frame) {
        if (count < 3)
            throw std::invalid_argument(
                "periodic Bishop frame needs at least three reference "
                "stations");
        const Vector3 closure_horizontal =
            TransportHorizontalDoubleReflection(
                reference_positions_m.back(), reference_positions_m.front(),
                tangents.back(), tangents.front(), horizontal.back());
        const double correction = std::atan2(
            Dot(tangents.front(), Cross(
                closure_horizontal, horizontal.front())),
            Dot(closure_horizontal, horizontal.front()));
        std::vector<double> cumulative(count, 0.0);
        for (std::size_t index = 1; index < count; ++index)
            cumulative[index] = cumulative[index - 1] + Norm(Subtract(
                reference_positions_m[index],
                reference_positions_m[index - 1]));
        const double total_length = cumulative.back() + Norm(Subtract(
            reference_positions_m.front(), reference_positions_m.back()));
        if (!std::isfinite(total_length) || total_length <= 0.0)
            throw std::invalid_argument(
                "periodic Bishop frame needs a nondegenerate closed path");
        for (std::size_t index = 1; index < count; ++index)
            horizontal[index] = RotateNormalAboutTangent(
                horizontal[index], tangents[index],
                correction * cumulative[index] / total_length);
        output.frame_holonomy_correction_rad = correction;
    }
    for (std::size_t index = 0; index < count; ++index)
        vertical[index] = Normalize(
            Cross(tangents[index], horizontal[index]), "vertical axis");

    ngcore::RegionTaskManager task_manager;
    ngstd::LocalHeap local_heap(1 << 20, "radia_beam_gridfunction");
    for (std::size_t index = 0; index < count; ++index) {
        output.linearizations.push_back(LinearizeSegment(
            *coefficient, mesh, reference_positions_m[index], tangents[index],
            horizontal[index], vertical[index], options, local_heap, index));

        dynamics[index].length_m = segment_lengths_m[index];
        dynamics[index].jet = output.linearizations.back().dynamics_jet;
        dynamics[index].name = names.empty()
            ? "segment_" + std::to_string(index + 1)
            : names[index];
    }

    VariationalOptions variational;
    variational.maximum_order = options.maximum_map_order;
    variational.maximum_step_m = options.maximum_step_m;
    variational.maximum_steps = options.maximum_steps;
    output.transfer = PropagateVariationalMap(dynamics, variational);
    return output;
}

GridFunctionTransferReport6 PropagateGridFunctionLinearMap(
        const std::shared_ptr<ngcomp::GridFunction>& field,
        const std::vector<double>& segment_lengths_m,
        const std::vector<Vector3>& reference_positions_m,
        const std::vector<Vector3>& reference_tangents,
        const std::vector<std::string>& names,
        const GridFunctionLinearizationOptions& options) {
    GridFunctionLinearizationOptions linear_options = options;
    linear_options.multipole_order = 1;
    linear_options.maximum_map_order = 1;
    return PropagateGridFunctionMultipoleMap(
        field, segment_lengths_m, reference_positions_m,
        reference_tangents, names, linear_options);
}

NGSolveGridFunctionField::NGSolveGridFunctionField(
        std::shared_ptr<ngcomp::GridFunction> field,
        GridFunctionMagneticInput magnetic_input)
    : field_(std::move(field)), magnetic_input_(magnetic_input) {
    if (!field_)
        throw std::invalid_argument("field GridFunction must not be null");
    const auto mesh = field_->GetMeshAccess();
    if (!mesh || mesh->GetDimension() != 3)
        throw std::invalid_argument(
            "beam GridFunction must belong to a three-dimensional mesh");
    magnetic_coefficient_ = MagneticCoefficient(field_, magnetic_input_);
}

FieldSample NGSolveGridFunctionField::Evaluate(
        const Vec3& position_m, double time_s,
        const FieldRequest& request) const {
    if (!std::isfinite(position_m.x) || !std::isfinite(position_m.y) ||
        !std::isfinite(position_m.z))
        throw std::invalid_argument(
            "position_m must contain finite values");
    if (!std::isfinite(time_s))
        throw std::invalid_argument("time_s must be finite");
    FieldSample output;
    if (!request.magnetic) return output;
    ngcore::RegionTaskManager task_manager;
    ngstd::LocalHeap local_heap(1 << 16, "radia_beam_direct_gridfunction");
    const Vector3 point{position_m.x, position_m.y, position_m.z};
    const Vector3 field = EvaluateField(
        *magnetic_coefficient_, field_->GetMeshAccess(), point,
        local_heap, 0, 0);
    output.magnetic_t = {field[0], field[1], field[2]};
    return output;
}

std::string NGSolveGridFunctionField::TypeName() const {
    return "ngsolve-grid-function-magnetic";
}

const std::shared_ptr<ngcomp::GridFunction>&
NGSolveGridFunctionField::GridFunction() const {
    return field_;
}

GridFunctionMagneticInput
NGSolveGridFunctionField::MagneticInput() const {
    return magnetic_input_;
}

}  // namespace radia::beam
