#include "rad_beam_transfer_pybind.h"

#include "rad_beam_ngsolve.h"
#include "rad_beam_transfer.h"

#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;

namespace {

using F64Array = py::array_t<double, py::array::c_style |
                                     py::array::forcecast>;

void RequireShape(const py::buffer_info& buffer,
                  const std::vector<py::ssize_t>& expected,
                  const char* name) {
    if (buffer.ndim != static_cast<py::ssize_t>(expected.size()))
        throw std::invalid_argument(std::string(name) +
                                    " has the wrong number of dimensions");
    for (py::ssize_t axis = 0; axis < buffer.ndim; ++axis) {
        if (buffer.shape[axis] != expected[static_cast<std::size_t>(axis)])
            throw std::invalid_argument(std::string(name) +
                                        " has the wrong shape");
    }
}

template <typename Value>
void CopyValue(Value& destination, const double* source) {
    std::copy(source, source + destination.values.size(),
              destination.values.begin());
}

std::vector<radia::beam::DynamicsSegment6> ReadSegments(
        const F64Array& lengths, const F64Array& a,
        const py::object& f2_object, const py::object& f3_object,
        const py::object& names_object) {
    const auto lengths_buffer = lengths.request();
    if (lengths_buffer.ndim != 1 || lengths_buffer.shape[0] <= 0)
        throw std::invalid_argument("lengths must have shape (n_segment,)");
    const py::ssize_t count = lengths_buffer.shape[0];
    const auto a_buffer = a.request();
    RequireShape(a_buffer, {count, 6, 6}, "A");

    F64Array f2;
    F64Array f3;
    const double* f2_data = nullptr;
    const double* f3_data = nullptr;
    if (!f2_object.is_none()) {
        f2 = F64Array::ensure(f2_object);
        if (!f2)
            throw std::invalid_argument("F2 must be a real C-contiguous array");
        const auto buffer = f2.request();
        RequireShape(buffer, {count, 6, 6, 6}, "F2");
        f2_data = static_cast<const double*>(buffer.ptr);
    }
    if (!f3_object.is_none()) {
        f3 = F64Array::ensure(f3_object);
        if (!f3)
            throw std::invalid_argument("F3 must be a real C-contiguous array");
        const auto buffer = f3.request();
        RequireShape(buffer, {count, 6, 6, 6, 6}, "F3");
        f3_data = static_cast<const double*>(buffer.ptr);
    }

    std::vector<std::string> names(static_cast<std::size_t>(count));
    if (!names_object.is_none()) {
        names = names_object.cast<std::vector<std::string>>();
        if (names.size() != static_cast<std::size_t>(count))
            throw std::invalid_argument(
                "names must contain one string per segment");
    }

    const double* length_data =
        static_cast<const double*>(lengths_buffer.ptr);
    const double* a_data = static_cast<const double*>(a_buffer.ptr);
    std::vector<radia::beam::DynamicsSegment6> segments(
        static_cast<std::size_t>(count));
    for (py::ssize_t index = 0; index < count; ++index) {
        auto& segment = segments[static_cast<std::size_t>(index)];
        segment.length_m = length_data[index];
        segment.name = names[static_cast<std::size_t>(index)];
        CopyValue(segment.jet.a_per_m, a_data + index * 36);
        if (f2_data)
            CopyValue(segment.jet.f2_per_m, f2_data + index * 216);
        if (f3_data)
            CopyValue(segment.jet.f3_per_m, f3_data + index * 1296);
    }
    return segments;
}

template <typename Value>
py::array_t<double> ValueArray(const Value& value,
                               const std::vector<py::ssize_t>& shape) {
    py::array_t<double> output(shape);
    std::memcpy(output.mutable_data(), value.values.data(),
                value.values.size() * sizeof(double));
    return output;
}

template <typename Value, typename Selector>
py::array_t<double> StackValues(std::size_t count,
                                const std::vector<py::ssize_t>& item_shape,
                                Selector selector) {
    std::vector<py::ssize_t> shape;
    shape.reserve(item_shape.size() + 1);
    shape.push_back(static_cast<py::ssize_t>(count));
    shape.insert(shape.end(), item_shape.begin(), item_shape.end());
    py::array_t<double> output(shape);
    double* destination = output.mutable_data();
    for (std::size_t index = 0; index < count; ++index) {
        const Value& value = selector(index);
        std::memcpy(destination + index * value.values.size(),
                    value.values.data(), value.values.size() * sizeof(double));
    }
    return output;
}

py::dict ReportDictionary(const radia::beam::VariationalReport6& report) {
    py::dict output;
    output["schema"] = "radia.beam.variational-map.result.v1";
    output["backend"] = "native-cpp";
    output["factorial_convention"] =
        "u_out = R*u + 1/2*T[u,u] + 1/6*U[u,u,u]";
    output["coordinate_order"] = py::make_tuple(
        "x", "px_over_p0", "y", "py_over_p0", "sigma", "delta");
    output["maximum_order"] = report.maximum_order;
    output["R"] = ValueArray(report.endpoint_map.r, {6, 6});
    output["T"] = ValueArray(report.endpoint_map.t, {6, 6, 6});
    output["U"] = ValueArray(report.endpoint_map.u, {6, 6, 6, 6});

    py::array_t<double> station_s(report.stations.size());
    for (std::size_t index = 0; index < report.stations.size(); ++index)
        station_s.mutable_data()[index] = report.stations[index].path_length_m;
    output["station_s_m"] = std::move(station_s);
    output["station_R"] = StackValues<radia::beam::Matrix6>(
        report.stations.size(), {6, 6},
        [&](std::size_t index) -> const radia::beam::Matrix6& {
            return report.stations[index].map_from_start.r;
        });
    output["station_T"] = StackValues<radia::beam::Tensor3Map6>(
        report.stations.size(), {6, 6, 6},
        [&](std::size_t index) -> const radia::beam::Tensor3Map6& {
            return report.stations[index].map_from_start.t;
        });
    output["station_U"] = StackValues<radia::beam::Tensor4Map6>(
        report.stations.size(), {6, 6, 6, 6},
        [&](std::size_t index) -> const radia::beam::Tensor4Map6& {
            return report.stations[index].map_from_start.u;
        });
    output["station_R_to_end"] = StackValues<radia::beam::Matrix6>(
        report.stations.size(), {6, 6},
        [&](std::size_t index) -> const radia::beam::Matrix6& {
            return report.stations[index].r_to_end;
        });

    py::list region_names;
    py::array_t<double> region_bounds(
        {static_cast<py::ssize_t>(report.regions.size()),
         static_cast<py::ssize_t>(2)});
    auto bounds = region_bounds.mutable_unchecked<2>();
    for (std::size_t index = 0; index < report.regions.size(); ++index) {
        region_names.append(report.regions[index].name);
        bounds(static_cast<py::ssize_t>(index), 0) =
            report.regions[index].s_begin_m;
        bounds(static_cast<py::ssize_t>(index), 1) =
            report.regions[index].s_end_m;
    }
    output["region_names"] = std::move(region_names);
    output["region_bounds_m"] = std::move(region_bounds);
    output["region_T"] = StackValues<radia::beam::Tensor3Map6>(
        report.regions.size(), {6, 6, 6},
        [&](std::size_t index) -> const radia::beam::Tensor3Map6& {
            return report.regions[index].t_at_end;
        });
    output["region_U_direct"] = StackValues<radia::beam::Tensor4Map6>(
        report.regions.size(), {6, 6, 6, 6},
        [&](std::size_t index) -> const radia::beam::Tensor4Map6& {
            return report.regions[index].u_direct_at_end;
        });
    output["region_U_local_cascade"] =
        StackValues<radia::beam::Tensor4Map6>(
            report.regions.size(), {6, 6, 6, 6},
            [&](std::size_t index) -> const radia::beam::Tensor4Map6& {
                return report.regions[index].u_local_cascade_at_end;
            });

    py::array_t<std::int64_t> pair_indices(
        {static_cast<py::ssize_t>(report.region_pairs.size()),
         static_cast<py::ssize_t>(2)});
    auto pairs = pair_indices.mutable_unchecked<2>();
    for (std::size_t index = 0; index < report.region_pairs.size(); ++index) {
        pairs(static_cast<py::ssize_t>(index), 0) =
            static_cast<std::int64_t>(
                report.region_pairs[index].upstream_region);
        pairs(static_cast<py::ssize_t>(index), 1) =
            static_cast<std::int64_t>(
                report.region_pairs[index].downstream_region);
    }
    output["pair_regions"] = std::move(pair_indices);
    output["pair_U_cascade"] = StackValues<radia::beam::Tensor4Map6>(
        report.region_pairs.size(), {6, 6, 6, 6},
        [&](std::size_t index) -> const radia::beam::Tensor4Map6& {
            return report.region_pairs[index].u_cascade_at_end;
        });

    py::dict diagnostics;
    diagnostics["integration_steps"] =
        report.diagnostics.integration_steps;
    diagnostics["R_composition_error"] =
        report.diagnostics.r_composition_error;
    diagnostics["T_reconstruction_error"] =
        report.diagnostics.t_reconstruction_error;
    diagnostics["U_reconstruction_error"] =
        report.diagnostics.u_reconstruction_error;
    diagnostics["T_symmetry_defect"] =
        report.diagnostics.t_input_symmetry_defect;
    diagnostics["U_symmetry_defect"] =
        report.diagnostics.u_input_symmetry_defect;
    output["diagnostics"] = std::move(diagnostics);
    return output;
}

std::vector<radia::beam::Vector3> Vector3Rows(
        const F64Array& values, py::ssize_t count, const char* name) {
    const auto buffer = values.request();
    RequireShape(buffer, {count, 3}, name);
    const double* data = static_cast<const double*>(buffer.ptr);
    std::vector<radia::beam::Vector3> output(
        static_cast<std::size_t>(count));
    for (py::ssize_t row = 0; row < count; ++row)
        for (py::ssize_t component = 0; component < 3; ++component)
            output[static_cast<std::size_t>(row)][
                static_cast<std::size_t>(component)] =
                data[row * 3 + component];
    return output;
}

py::dict GridFunctionReportDictionary(
        const radia::beam::GridFunctionTransferReport6& report) {
    py::dict output = ReportDictionary(report.transfer);
    output["schema"] =
        "radia.beam.grid-function-linear-map.result.v1";
    output["backend"] = "native-cpp-ngsolve-gridfunction";
    output["field_source"] = "ngsolve.GridFunction";
    output["linearization_order"] = 1;
    output["magnetic_rigidity_t_m"] = report.magnetic_rigidity_t_m;
    output["sample_radius_m"] = report.sample_radius_m;
    output["fit_model"] =
        "nine-point transverse least-squares affine field jet";
    output["frame_convention"] =
        "right-handed parallel transport seeded by initial_horizontal";

    const py::ssize_t count = static_cast<py::ssize_t>(
        report.linearizations.size());
    py::array_t<double> positions(
        std::vector<py::ssize_t>{count, 3});
    py::array_t<double> horizontal(
        std::vector<py::ssize_t>{count, 3});
    py::array_t<double> vertical(
        std::vector<py::ssize_t>{count, 3});
    py::array_t<double> tangent(
        std::vector<py::ssize_t>{count, 3});
    py::array_t<double> center_field(
        std::vector<py::ssize_t>{count, 3});
    py::array_t<double> fitted_center_field(
        std::vector<py::ssize_t>{count, 3});
    py::array_t<double> gradients(
        std::vector<py::ssize_t>{count, 3, 2});
    py::array_t<double> raw_positions(
        std::vector<py::ssize_t>{count, 9, 3});
    py::array_t<double> raw_fields(
        std::vector<py::ssize_t>{count, 9, 3});
    py::array_t<double> curvature(count);
    py::array_t<double> normal_gradient(count);
    py::array_t<double> skew_gradient(count);
    py::array_t<double> divergence(count);
    py::array_t<double> curl_mismatch(count);
    py::array_t<double> center_bias(count);
    py::array_t<double> rms_residual(count);
    py::array_t<double> maximum_residual(count);
    py::array_t<std::int64_t> fit_rank(count);
    py::array_t<double> fit_condition(count);

    for (py::ssize_t index = 0; index < count; ++index) {
        const auto& value = report.linearizations[
            static_cast<std::size_t>(index)];
        for (py::ssize_t component = 0; component < 3; ++component) {
            positions.mutable_data()[index * 3 + component] =
                value.reference_position_m[component];
            horizontal.mutable_data()[index * 3 + component] =
                value.horizontal_axis[component];
            vertical.mutable_data()[index * 3 + component] =
                value.vertical_axis[component];
            tangent.mutable_data()[index * 3 + component] =
                value.tangent_axis[component];
            center_field.mutable_data()[index * 3 + component] =
                value.center_field_local_t[component];
            fitted_center_field.mutable_data()[index * 3 + component] =
                value.fitted_center_field_local_t[component];
            for (py::ssize_t derivative = 0; derivative < 2; ++derivative)
                gradients.mutable_data()[
                    (index * 3 + component) * 2 + derivative] =
                    value.field_gradient_local_t_per_m[
                        component * 2 + derivative];
        }
        for (py::ssize_t sample = 0; sample < 9; ++sample)
            for (py::ssize_t component = 0; component < 3; ++component) {
                raw_positions.mutable_data()[
                    (index * 9 + sample) * 3 + component] =
                    value.sample_positions_m[sample * 3 + component];
                raw_fields.mutable_data()[
                    (index * 9 + sample) * 3 + component] =
                    value.sample_fields_local_t[sample * 3 + component];
            }
        curvature.mutable_data()[index] = value.curvature_per_m;
        normal_gradient.mutable_data()[index] =
            value.normal_gradient_per_m2;
        skew_gradient.mutable_data()[index] = value.skew_gradient_per_m2;
        divergence.mutable_data()[index] =
            value.transverse_divergence_t_per_m;
        curl_mismatch.mutable_data()[index] =
            value.transverse_curl_mismatch_t_per_m;
        center_bias.mutable_data()[index] = value.center_fit_bias_t;
        rms_residual.mutable_data()[index] = value.rms_fit_residual_t;
        maximum_residual.mutable_data()[index] =
            value.maximum_fit_residual_t;
        fit_rank.mutable_data()[index] = static_cast<std::int64_t>(
            value.fit_rank);
        fit_condition.mutable_data()[index] =
            value.scaled_design_condition;
    }
    output["reference_positions_m"] = std::move(positions);
    output["frame_horizontal"] = std::move(horizontal);
    output["frame_vertical"] = std::move(vertical);
    output["frame_tangent"] = std::move(tangent);
    output["center_field_local_t"] = std::move(center_field);
    output["fitted_center_field_local_t"] =
        std::move(fitted_center_field);
    output["field_gradient_local_t_per_m"] = std::move(gradients);
    output["field_sample_positions_m"] = std::move(raw_positions);
    output["field_samples_local_t"] = std::move(raw_fields);
    output["curvature_per_m"] = std::move(curvature);
    output["normal_gradient_per_m2"] = std::move(normal_gradient);
    output["skew_gradient_per_m2"] = std::move(skew_gradient);
    output["transverse_divergence_t_per_m"] = std::move(divergence);
    output["transverse_curl_mismatch_t_per_m"] =
        std::move(curl_mismatch);
    output["center_fit_bias_t"] = std::move(center_bias);
    output["rms_fit_residual_t"] = std::move(rms_residual);
    output["maximum_fit_residual_t"] = std::move(maximum_residual);
    output["fit_rank"] = std::move(fit_rank);
    output["scaled_design_condition"] = std::move(fit_condition);
    output["local_A_per_m"] = StackValues<radia::beam::Matrix6>(
        report.linearizations.size(), {6, 6},
        [&](std::size_t index) -> const radia::beam::Matrix6& {
            return report.linearizations[index].a_per_m;
        });
    return output;
}

}  // namespace

void ExportBeamTransfer(py::module_& module) {
    module.def(
        "_beam_variational_map",
        [](F64Array lengths, F64Array a, py::object f2, py::object f3,
           py::object names, unsigned maximum_order, double maximum_step_m,
           std::size_t maximum_steps, std::size_t maximum_region_pairs,
           double input_symmetry_tolerance) {
            auto segments = ReadSegments(lengths, a, f2, f3, names);
            radia::beam::VariationalOptions options;
            options.maximum_order = maximum_order;
            options.maximum_step_m = maximum_step_m;
            options.maximum_steps = maximum_steps;
            options.maximum_region_pairs = maximum_region_pairs;
            options.input_symmetry_tolerance = input_symmetry_tolerance;
            radia::beam::VariationalReport6 report;
            {
                py::gil_scoped_release release;
                report = radia::beam::PropagateVariationalMap(
                    segments, options);
            }
            return ReportDictionary(report);
        },
        py::arg("lengths_m"), py::arg("A_per_m"),
        py::arg("F2_per_m") = py::none(),
        py::arg("F3_per_m") = py::none(),
        py::arg("names") = py::none(), py::arg("maximum_order") = 3,
        py::arg("maximum_step_m") = 1.0e-3,
        py::arg("maximum_steps") = 1000000,
        py::arg("maximum_region_pairs") = 100000,
        py::arg("input_symmetry_tolerance") = 1.0e-12,
        R"pbdoc(
Propagate a canonical six-dimensional Taylor map through constant jet regions.

The native C++ kernel integrates R/T/U through third order and returns
region-resolved T, direct/local U, and ordered upstream-to-downstream U
cascade attribution. Arrays are C-order with shapes (n,6,6), (n,6,6,6),
and (n,6,6,6,6). The map convention is
u_out = R*u + 1/2*T[u,u] + 1/6*U[u,u,u].
)pbdoc");

    module.def(
        "_beam_grid_function_linear_map",
        [](py::object field_object, F64Array lengths,
           F64Array reference_positions, F64Array reference_tangents,
           double magnetic_rigidity_t_m, F64Array initial_horizontal,
           double sample_radius_m, py::object names_object,
           double curvature_sign, double gradient_sign,
           double maximum_step_m, std::size_t maximum_steps) {
            const auto lengths_buffer = lengths.request();
            if (lengths_buffer.ndim != 1 || lengths_buffer.shape[0] <= 0)
                throw std::invalid_argument(
                    "lengths_m must have shape (n_segment,)");
            const py::ssize_t count = lengths_buffer.shape[0];
            const double* length_data = static_cast<const double*>(
                lengths_buffer.ptr);
            std::vector<double> segment_lengths(
                length_data, length_data + count);
            auto positions = Vector3Rows(
                reference_positions, count, "reference_positions_m");
            auto tangents = Vector3Rows(
                reference_tangents, count, "reference_tangents");
            const auto horizontal_buffer = initial_horizontal.request();
            RequireShape(horizontal_buffer, {3}, "initial_horizontal");
            const double* horizontal_data = static_cast<const double*>(
                horizontal_buffer.ptr);

            std::vector<std::string> names;
            if (!names_object.is_none()) {
                names = names_object.cast<std::vector<std::string>>();
                if (names.size() != static_cast<std::size_t>(count))
                    throw std::invalid_argument(
                        "names must contain one string per segment");
            }
            auto field = field_object.cast<
                std::shared_ptr<ngcomp::GridFunction>>();
            radia::beam::GridFunctionLinearizationOptions options;
            options.magnetic_rigidity_t_m = magnetic_rigidity_t_m;
            options.sample_radius_m = sample_radius_m;
            std::copy(horizontal_data, horizontal_data + 3,
                      options.initial_horizontal.begin());
            options.curvature_sign = curvature_sign;
            options.gradient_sign = gradient_sign;
            options.maximum_step_m = maximum_step_m;
            options.maximum_steps = maximum_steps;

            radia::beam::GridFunctionTransferReport6 report;
            {
                py::gil_scoped_release release;
                report = radia::beam::PropagateGridFunctionLinearMap(
                    field, segment_lengths, positions, tangents, names,
                    options);
            }
            return GridFunctionReportDictionary(report);
        },
        py::arg("field"), py::arg("lengths_m"),
        py::arg("reference_positions_m"),
        py::arg("reference_tangents"),
        py::arg("magnetic_rigidity_t_m"),
        py::arg("initial_horizontal"),
        py::arg("sample_radius_m") = 1.0e-3,
        py::arg("names") = py::none(),
        py::arg("curvature_sign") = 1.0,
        py::arg("gradient_sign") = 1.0,
        py::arg("maximum_step_m") = 1.0e-3,
        py::arg("maximum_steps") = 1000000,
        R"pbdoc(
Build a first-order transfer map by evaluating an NGSolve GridFunction directly.

The C++ adapter samples nine transverse points at each reference station via
NGSolve's mapped GridFunction evaluation. It reports the local field fit,
normal/skew quadrupole profiles, residual diagnostics, and the accumulated
six-dimensional R map without constructing a regular-grid field map.
)pbdoc");
}
