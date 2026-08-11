#include "radia_beam_mex_commands.h"

#include "rad_beam_ngsolve.h"
#include "rad_beam_transfer.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

const mxArray* Field(const mxArray* value, const char* name) {
    return value && mxIsStruct(value) ? mxGetField(value, 0, name) : nullptr;
}

void RequireScalarStruct(const mxArray* value) {
    if (!value || !mxIsStruct(value) || mxGetNumberOfElements(value) != 1)
        throw std::invalid_argument(
            "beam variational map requires a scalar configuration struct");
}

std::string Text(const mxArray* value, const char* name) {
    if (!value || !mxIsChar(value))
        throw std::invalid_argument(std::string(name) + " must be text");
    char* buffer = mxArrayToUTF8String(value);
    if (!buffer)
        throw std::invalid_argument(std::string(name) + " is invalid");
    std::string result(buffer);
    mxFree(buffer);
    return result;
}

double Scalar(const mxArray* value, const char* name, double fallback,
              bool optional) {
    if (!value && optional) return fallback;
    if (!value || !mxIsDouble(value) || mxIsComplex(value) ||
        mxGetNumberOfElements(value) != 1)
        throw std::invalid_argument(std::string(name) +
                                    " must be a real double scalar");
    const double result = mxGetScalar(value);
    if (!std::isfinite(result))
        throw std::invalid_argument(std::string(name) + " must be finite");
    return result;
}

std::size_t PositiveInteger(const mxArray* value, const char* name,
                            std::size_t fallback, bool optional) {
    const double number = Scalar(
        value, name, static_cast<double>(fallback), optional);
    if (number < 1.0 || number != std::floor(number) ||
        number > static_cast<double>(
            std::numeric_limits<std::size_t>::max()))
        throw std::invalid_argument(std::string(name) +
                                    " must be a positive integer");
    return static_cast<std::size_t>(number);
}

const double* RealData(const mxArray* value, const char* name) {
    if (!value || !mxIsDouble(value) || mxIsComplex(value))
        throw std::invalid_argument(std::string(name) +
                                    " must be a real double array");
    return mxGetPr(value);
}

void RequireLeadingShape(const mxArray* value, const char* name,
                         std::size_t leading_rank, std::size_t count) {
    (void)RealData(value, name);
    const mwSize rank = mxGetNumberOfDimensions(value);
    const mwSize* dimensions = mxGetDimensions(value);
    if (rank < leading_rank)
        throw std::invalid_argument(std::string(name) +
                                    " has too few dimensions");
    std::size_t expected = count;
    for (std::size_t axis = 0; axis < leading_rank; ++axis) {
        if (dimensions[axis] != 6)
            throw std::invalid_argument(std::string(name) +
                                        " must have leading dimensions of 6");
        expected *= 6;
    }
    if (mxGetNumberOfElements(value) != expected)
        throw std::invalid_argument(std::string(name) +
                                    " has the wrong segment count");
}

std::vector<std::string> Names(const mxArray* value, std::size_t count) {
    std::vector<std::string> result(count);
    if (!value) return result;
    if (!mxIsCell(value) || mxGetNumberOfElements(value) != count)
        throw std::invalid_argument(
            "names must be a cell array with one entry per segment");
    for (std::size_t index = 0; index < count; ++index)
        result[index] = Text(mxGetCell(value, index), "names entry");
    return result;
}

std::vector<radia::beam::DynamicsSegment6> ParseSegments(
        const mxArray* config) {
    const mxArray* lengths_value = Field(config, "lengths_m");
    const double* lengths = RealData(lengths_value, "lengths_m");
    const std::size_t count = mxGetNumberOfElements(lengths_value);
    if (count == 0)
        throw std::invalid_argument("lengths_m must not be empty");
    const mxArray* a_value = Field(config, "A_per_m");
    RequireLeadingShape(a_value, "A_per_m", 2, count);
    const double* a = mxGetPr(a_value);
    const mxArray* f2_value = Field(config, "F2_per_m");
    const mxArray* f3_value = Field(config, "F3_per_m");
    const double* f2 = nullptr;
    const double* f3 = nullptr;
    if (f2_value) {
        RequireLeadingShape(f2_value, "F2_per_m", 3, count);
        f2 = mxGetPr(f2_value);
    }
    if (f3_value) {
        RequireLeadingShape(f3_value, "F3_per_m", 4, count);
        f3 = mxGetPr(f3_value);
    }
    const std::vector<std::string> names = Names(Field(config, "names"), count);

    std::vector<radia::beam::DynamicsSegment6> segments(count);
    for (std::size_t segment = 0; segment < count; ++segment) {
        auto& output = segments[segment];
        output.length_m = lengths[segment];
        output.name = names[segment];
        for (std::size_t i = 0; i < 6; ++i)
            for (std::size_t j = 0; j < 6; ++j)
                output.jet.a_per_m(i, j) =
                    a[i + 6 * (j + 6 * segment)];
        if (f2) {
            for (std::size_t i = 0; i < 6; ++i)
                for (std::size_t j = 0; j < 6; ++j)
                    for (std::size_t k = 0; k < 6; ++k)
                        output.jet.f2_per_m(i, j, k) =
                            f2[i + 6 * (j + 6 * (k + 6 * segment))];
        }
        if (f3) {
            for (std::size_t i = 0; i < 6; ++i)
                for (std::size_t j = 0; j < 6; ++j)
                    for (std::size_t k = 0; k < 6; ++k)
                        for (std::size_t l = 0; l < 6; ++l)
                            output.jet.f3_per_m(i, j, k, l) =
                                f3[i + 6 * (j + 6 * (k + 6 *
                                    (l + 6 * segment)))];
        }
    }
    return segments;
}

std::vector<radia::beam::Vector3> Rows3(
        const mxArray* value, std::size_t count, const char* name) {
    const double* data = RealData(value, name);
    if (mxGetNumberOfDimensions(value) != 2 ||
        mxGetM(value) != count || mxGetN(value) != 3)
        throw std::invalid_argument(
            std::string(name) + " must have shape (n_segment,3)");
    std::vector<radia::beam::Vector3> output(count);
    for (std::size_t row = 0; row < count; ++row)
        for (std::size_t component = 0; component < 3; ++component)
            output[row][component] = data[row + count * component];
    return output;
}

radia::beam::Vector3 OptionalVector3(
        const mxArray* value, const char* name,
        radia::beam::Vector3 fallback) {
    if (!value) return fallback;
    const double* data = RealData(value, name);
    if (mxGetNumberOfElements(value) != 3)
        throw std::invalid_argument(std::string(name) +
                                    " must contain three entries");
    return {data[0], data[1], data[2]};
}

struct GridFunctionInput {
    std::vector<double> lengths_m;
    std::vector<radia::beam::Vector3> positions_m;
    std::vector<radia::beam::Vector3> tangents;
    std::vector<std::string> names;
    radia::beam::GridFunctionLinearizationOptions options;
};

GridFunctionInput ParseGridFunctionInput(const mxArray* config) {
    RequireScalarStruct(config);
    const mxArray* schema = Field(config, "schema");
    if (schema && Text(schema, "schema") !=
            "radia.beam.grid-function-linear-map.v1")
        throw std::invalid_argument(
            "unsupported beam GridFunction linear-map schema");
    const mxArray* lengths_value = Field(config, "lengths_m");
    const double* lengths = RealData(lengths_value, "lengths_m");
    const std::size_t count = mxGetNumberOfElements(lengths_value);
    if (count == 0)
        throw std::invalid_argument("lengths_m must not be empty");

    GridFunctionInput output;
    output.lengths_m.assign(lengths, lengths + count);
    output.positions_m = Rows3(
        Field(config, "reference_positions_m"), count,
        "reference_positions_m");
    output.tangents = Rows3(
        Field(config, "reference_tangents"), count,
        "reference_tangents");
    output.names = Names(Field(config, "names"), count);
    if (!Field(config, "names")) output.names.clear();
    output.options.magnetic_rigidity_t_m = Scalar(
        Field(config, "magnetic_rigidity_t_m"),
        "magnetic_rigidity_t_m", 0.0, false);
    output.options.sample_radius_m = Scalar(
        Field(config, "sample_radius_m"), "sample_radius_m", 1.0e-3, true);
    output.options.initial_horizontal = OptionalVector3(
        Field(config, "initial_horizontal"), "initial_horizontal",
        {1.0, 0.0, 0.0});
    output.options.curvature_sign = Scalar(
        Field(config, "curvature_sign"), "curvature_sign", 1.0, true);
    output.options.gradient_sign = Scalar(
        Field(config, "gradient_sign"), "gradient_sign", 1.0, true);
    output.options.maximum_step_m = Scalar(
        Field(config, "maximum_step_m"), "maximum_step_m", 1.0e-3, true);
    output.options.maximum_steps = PositiveInteger(
        Field(config, "maximum_steps"), "maximum_steps", 1000000, true);
    return output;
}

mxArray* MatrixArray(const radia::beam::Matrix6& value) {
    mxArray* output = mxCreateDoubleMatrix(6, 6, mxREAL);
    double* data = mxGetPr(output);
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            data[i + 6 * j] = value(i, j);
    return output;
}

mxArray* Tensor3Array(const radia::beam::Tensor3Map6& value) {
    const mwSize dimensions[] = {6, 6, 6};
    mxArray* output = mxCreateNumericArray(3, dimensions, mxDOUBLE_CLASS,
                                           mxREAL);
    double* data = mxGetPr(output);
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            for (std::size_t k = 0; k < 6; ++k)
                data[i + 6 * (j + 6 * k)] = value(i, j, k);
    return output;
}

mxArray* Tensor4Array(const radia::beam::Tensor4Map6& value) {
    const mwSize dimensions[] = {6, 6, 6, 6};
    mxArray* output = mxCreateNumericArray(4, dimensions, mxDOUBLE_CLASS,
                                           mxREAL);
    double* data = mxGetPr(output);
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            for (std::size_t k = 0; k < 6; ++k)
                for (std::size_t l = 0; l < 6; ++l)
                    data[i + 6 * (j + 6 * (k + 6 * l))] = value(i, j, k, l);
    return output;
}

template <typename Selector>
mxArray* MatrixStack(std::size_t count, Selector selector) {
    const mwSize dimensions[] = {6, 6, static_cast<mwSize>(count)};
    mxArray* output = mxCreateNumericArray(3, dimensions, mxDOUBLE_CLASS,
                                           mxREAL);
    double* data = mxGetPr(output);
    for (std::size_t item = 0; item < count; ++item) {
        const auto& value = selector(item);
        for (std::size_t i = 0; i < 6; ++i)
            for (std::size_t j = 0; j < 6; ++j)
                data[i + 6 * (j + 6 * item)] = value(i, j);
    }
    return output;
}

template <typename Selector>
mxArray* Tensor3Stack(std::size_t count, Selector selector) {
    const mwSize dimensions[] = {6, 6, 6, static_cast<mwSize>(count)};
    mxArray* output = mxCreateNumericArray(4, dimensions, mxDOUBLE_CLASS,
                                           mxREAL);
    double* data = mxGetPr(output);
    for (std::size_t item = 0; item < count; ++item) {
        const auto& value = selector(item);
        for (std::size_t i = 0; i < 6; ++i)
            for (std::size_t j = 0; j < 6; ++j)
                for (std::size_t k = 0; k < 6; ++k)
                    data[i + 6 * (j + 6 * (k + 6 * item))] = value(i, j, k);
    }
    return output;
}

template <typename Selector>
mxArray* Tensor4Stack(std::size_t count, Selector selector) {
    const mwSize dimensions[] = {6, 6, 6, 6, static_cast<mwSize>(count)};
    mxArray* output = mxCreateNumericArray(5, dimensions, mxDOUBLE_CLASS,
                                           mxREAL);
    double* data = mxGetPr(output);
    for (std::size_t item = 0; item < count; ++item) {
        const auto& value = selector(item);
        for (std::size_t i = 0; i < 6; ++i)
            for (std::size_t j = 0; j < 6; ++j)
                for (std::size_t k = 0; k < 6; ++k)
                    for (std::size_t l = 0; l < 6; ++l)
                        data[i + 6 * (j + 6 * (k + 6 *
                            (l + 6 * item)))] = value(i, j, k, l);
    }
    return output;
}

mxArray* StringCell(const std::vector<std::string>& values) {
    mxArray* output = mxCreateCellMatrix(values.size(), 1);
    for (std::size_t index = 0; index < values.size(); ++index)
        mxSetCell(output, index, mxCreateString(values[index].c_str()));
    return output;
}

mxArray* Diagnostics(const radia::beam::VariationalDiagnostics& value) {
    const char* fields[] = {
        "integration_steps", "R_composition_error",
        "T_reconstruction_error", "U_reconstruction_error",
        "T_symmetry_defect", "U_symmetry_defect"};
    mxArray* output = mxCreateStructMatrix(1, 1, 6, fields);
    mxSetField(output, 0, fields[0], mxCreateDoubleScalar(
        static_cast<double>(value.integration_steps)));
    mxSetField(output, 0, fields[1],
               mxCreateDoubleScalar(value.r_composition_error));
    mxSetField(output, 0, fields[2],
               mxCreateDoubleScalar(value.t_reconstruction_error));
    mxSetField(output, 0, fields[3],
               mxCreateDoubleScalar(value.u_reconstruction_error));
    mxSetField(output, 0, fields[4],
               mxCreateDoubleScalar(value.t_input_symmetry_defect));
    mxSetField(output, 0, fields[5],
               mxCreateDoubleScalar(value.u_input_symmetry_defect));
    return output;
}

mxArray* Report(const radia::beam::VariationalReport6& value) {
    const char* fields[] = {
        "schema", "backend", "factorial_convention", "maximum_order",
        "R", "T", "U", "station_s_m", "station_R",
        "station_R_to_end", "region_names", "region_bounds_m", "region_T",
        "region_U_direct", "region_U_local_cascade", "pair_regions",
        "pair_U_cascade", "diagnostics", "coordinate_order"};
    mxArray* output = mxCreateStructMatrix(1, 1, 19, fields);
    mxSetField(output, 0, fields[0],
               mxCreateString("radia.beam.variational-map.result.v1"));
    mxSetField(output, 0, fields[1], mxCreateString("native-cpp-mex"));
    mxSetField(output, 0, fields[2], mxCreateString(
        "u_out = R*u + 1/2*T[u,u] + 1/6*U[u,u,u]"));
    mxSetField(output, 0, fields[3],
               mxCreateDoubleScalar(value.maximum_order));
    mxSetField(output, 0, fields[4], MatrixArray(value.endpoint_map.r));
    mxSetField(output, 0, fields[5], Tensor3Array(value.endpoint_map.t));
    mxSetField(output, 0, fields[6], Tensor4Array(value.endpoint_map.u));

    mxArray* station_s = mxCreateDoubleMatrix(value.stations.size(), 1, mxREAL);
    for (std::size_t index = 0; index < value.stations.size(); ++index)
        mxGetPr(station_s)[index] = value.stations[index].path_length_m;
    mxSetField(output, 0, fields[7], station_s);
    mxSetField(output, 0, fields[8], MatrixStack(
        value.stations.size(), [&](std::size_t index) -> const auto& {
            return value.stations[index].map_from_start.r;
        }));
    mxSetField(output, 0, fields[9], MatrixStack(
        value.stations.size(), [&](std::size_t index) -> const auto& {
            return value.stations[index].r_to_end;
        }));

    std::vector<std::string> names;
    names.reserve(value.regions.size());
    mxArray* bounds = mxCreateDoubleMatrix(value.regions.size(), 2, mxREAL);
    for (std::size_t index = 0; index < value.regions.size(); ++index) {
        names.push_back(value.regions[index].name);
        mxGetPr(bounds)[index] = value.regions[index].s_begin_m;
        mxGetPr(bounds)[index + value.regions.size()] =
            value.regions[index].s_end_m;
    }
    mxSetField(output, 0, fields[10], StringCell(names));
    mxSetField(output, 0, fields[11], bounds);
    mxSetField(output, 0, fields[12], Tensor3Stack(
        value.regions.size(), [&](std::size_t index) -> const auto& {
            return value.regions[index].t_at_end;
        }));
    mxSetField(output, 0, fields[13], Tensor4Stack(
        value.regions.size(), [&](std::size_t index) -> const auto& {
            return value.regions[index].u_direct_at_end;
        }));
    mxSetField(output, 0, fields[14], Tensor4Stack(
        value.regions.size(), [&](std::size_t index) -> const auto& {
            return value.regions[index].u_local_cascade_at_end;
        }));

    mxArray* pair_regions = mxCreateDoubleMatrix(
        value.region_pairs.size(), 2, mxREAL);
    for (std::size_t index = 0; index < value.region_pairs.size(); ++index) {
        mxGetPr(pair_regions)[index] =
            static_cast<double>(value.region_pairs[index].upstream_region + 1);
        mxGetPr(pair_regions)[index + value.region_pairs.size()] =
            static_cast<double>(value.region_pairs[index].downstream_region + 1);
    }
    mxSetField(output, 0, fields[15], pair_regions);
    mxSetField(output, 0, fields[16], Tensor4Stack(
        value.region_pairs.size(), [&](std::size_t index) -> const auto& {
            return value.region_pairs[index].u_cascade_at_end;
        }));
    mxSetField(output, 0, fields[17], Diagnostics(value.diagnostics));
    mxSetField(output, 0, fields[18], StringCell(
        {"x", "px_over_p0", "y", "py_over_p0", "sigma", "delta"}));
    return output;
}

void ReplaceStructField(mxArray* output, const char* name, mxArray* value) {
    mxArray* previous = mxGetField(output, 0, name);
    mxSetField(output, 0, name, value);
    if (previous) mxDestroyArray(previous);
}

void AddStructField(mxArray* output, const char* name, mxArray* value) {
    if (mxAddField(output, name) < 0) {
        mxDestroyArray(value);
        throw std::runtime_error(std::string("failed to add result field ") +
                                 name);
    }
    mxSetField(output, 0, name, value);
}

template <typename Selector>
mxArray* Vector3Rows(std::size_t count, Selector selector) {
    mxArray* output = mxCreateDoubleMatrix(count, 3, mxREAL);
    double* data = mxGetPr(output);
    for (std::size_t row = 0; row < count; ++row) {
        const auto& value = selector(row);
        for (std::size_t component = 0; component < 3; ++component)
            data[row + count * component] = value[component];
    }
    return output;
}

template <typename Selector>
mxArray* ScalarColumn(std::size_t count, Selector selector) {
    mxArray* output = mxCreateDoubleMatrix(count, 1, mxREAL);
    double* data = mxGetPr(output);
    for (std::size_t row = 0; row < count; ++row)
        data[row] = selector(row);
    return output;
}

mxArray* FieldGradientStack(
        const std::vector<radia::beam::GridFunctionSegmentLinearization>&
            values) {
    const mwSize dimensions[] = {
        static_cast<mwSize>(values.size()), 3, 2};
    mxArray* output = mxCreateNumericArray(3, dimensions, mxDOUBLE_CLASS,
                                           mxREAL);
    double* data = mxGetPr(output);
    const std::size_t count = values.size();
    for (std::size_t row = 0; row < count; ++row)
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t derivative = 0; derivative < 2; ++derivative)
                data[row + count * (component + 3 * derivative)] =
                    values[row].field_gradient_local_t_per_m[
                        2 * component + derivative];
    return output;
}

template <typename Selector>
mxArray* SampleStack(std::size_t count, Selector selector) {
    const mwSize dimensions[] = {9, 3, static_cast<mwSize>(count)};
    mxArray* output = mxCreateNumericArray(3, dimensions, mxDOUBLE_CLASS,
                                           mxREAL);
    double* data = mxGetPr(output);
    for (std::size_t item = 0; item < count; ++item) {
        const auto& value = selector(item);
        for (std::size_t sample = 0; sample < 9; ++sample)
            for (std::size_t component = 0; component < 3; ++component)
                data[sample + 9 * (component + 3 * item)] =
                    value[3 * sample + component];
    }
    return output;
}

mxArray* GridFunctionReport(
        const radia::beam::GridFunctionTransferReport6& value) {
    mxArray* output = Report(value.transfer);
    ReplaceStructField(output, "schema", mxCreateString(
        "radia.beam.grid-function-linear-map.result.v1"));
    ReplaceStructField(output, "backend", mxCreateString(
        "native-cpp-ngsolve-gridfunction-mex"));
    AddStructField(output, "field_source",
                   mxCreateString("ngsolve.GridFunction"));
    AddStructField(output, "linearization_order", mxCreateDoubleScalar(1.0));
    AddStructField(output, "magnetic_rigidity_t_m",
                   mxCreateDoubleScalar(value.magnetic_rigidity_t_m));
    AddStructField(output, "sample_radius_m",
                   mxCreateDoubleScalar(value.sample_radius_m));
    AddStructField(output, "fit_model", mxCreateString(
        "nine-point transverse least-squares affine field jet"));
    AddStructField(output, "frame_convention", mxCreateString(
        "right-handed parallel transport seeded by initial_horizontal"));

    const auto& samples = value.linearizations;
    const std::size_t count = samples.size();
    AddStructField(output, "reference_positions_m", Vector3Rows(
        count, [&](std::size_t index) -> const auto& {
            return samples[index].reference_position_m;
        }));
    AddStructField(output, "frame_horizontal", Vector3Rows(
        count, [&](std::size_t index) -> const auto& {
            return samples[index].horizontal_axis;
        }));
    AddStructField(output, "frame_vertical", Vector3Rows(
        count, [&](std::size_t index) -> const auto& {
            return samples[index].vertical_axis;
        }));
    AddStructField(output, "frame_tangent", Vector3Rows(
        count, [&](std::size_t index) -> const auto& {
            return samples[index].tangent_axis;
        }));
    AddStructField(output, "center_field_local_t", Vector3Rows(
        count, [&](std::size_t index) -> const auto& {
            return samples[index].center_field_local_t;
        }));
    AddStructField(output, "fitted_center_field_local_t", Vector3Rows(
        count, [&](std::size_t index) -> const auto& {
            return samples[index].fitted_center_field_local_t;
        }));
    AddStructField(output, "field_gradient_local_t_per_m",
                   FieldGradientStack(samples));
    AddStructField(output, "field_sample_positions_m", SampleStack(
        count, [&](std::size_t index) -> const auto& {
            return samples[index].sample_positions_m;
        }));
    AddStructField(output, "field_samples_local_t", SampleStack(
        count, [&](std::size_t index) -> const auto& {
            return samples[index].sample_fields_local_t;
        }));
    AddStructField(output, "curvature_per_m", ScalarColumn(
        count, [&](std::size_t index) {
            return samples[index].curvature_per_m;
        }));
    AddStructField(output, "normal_gradient_per_m2", ScalarColumn(
        count, [&](std::size_t index) {
            return samples[index].normal_gradient_per_m2;
        }));
    AddStructField(output, "skew_gradient_per_m2", ScalarColumn(
        count, [&](std::size_t index) {
            return samples[index].skew_gradient_per_m2;
        }));
    AddStructField(output, "transverse_divergence_t_per_m", ScalarColumn(
        count, [&](std::size_t index) {
            return samples[index].transverse_divergence_t_per_m;
        }));
    AddStructField(output, "transverse_curl_mismatch_t_per_m", ScalarColumn(
        count, [&](std::size_t index) {
            return samples[index].transverse_curl_mismatch_t_per_m;
        }));
    AddStructField(output, "center_fit_bias_t", ScalarColumn(
        count, [&](std::size_t index) {
            return samples[index].center_fit_bias_t;
        }));
    AddStructField(output, "rms_fit_residual_t", ScalarColumn(
        count, [&](std::size_t index) {
            return samples[index].rms_fit_residual_t;
        }));
    AddStructField(output, "maximum_fit_residual_t", ScalarColumn(
        count, [&](std::size_t index) {
            return samples[index].maximum_fit_residual_t;
        }));
    AddStructField(output, "fit_rank", ScalarColumn(
        count, [&](std::size_t index) {
            return static_cast<double>(samples[index].fit_rank);
        }));
    AddStructField(output, "scaled_design_condition", ScalarColumn(
        count, [&](std::size_t index) {
            return samples[index].scaled_design_condition;
        }));
    AddStructField(output, "local_A_per_m", MatrixStack(
        count, [&](std::size_t index) -> const auto& {
            return samples[index].a_per_m;
        }));
    return output;
}

void Propagate(int nlhs, mxArray* plhs[], int nrhs,
               const mxArray* prhs[]) {
    if (nrhs != 2 || nlhs != 1)
        throw std::invalid_argument(
            "result = radia_mex('beam.transfer.propagate_variational', config)");
    const mxArray* config = prhs[1];
    RequireScalarStruct(config);
    const mxArray* schema = Field(config, "schema");
    if (schema && Text(schema, "schema") != "radia.beam.variational-map.v1")
        throw std::invalid_argument("unsupported beam variational-map schema");
    std::vector<radia::beam::DynamicsSegment6> segments =
        ParseSegments(config);
    radia::beam::VariationalOptions options;
    options.maximum_order = static_cast<unsigned>(PositiveInteger(
        Field(config, "maximum_order"), "maximum_order", 3, true));
    options.maximum_step_m = Scalar(
        Field(config, "maximum_step_m"), "maximum_step_m", 1.0e-3, true);
    options.maximum_steps = PositiveInteger(
        Field(config, "maximum_steps"), "maximum_steps", 1000000, true);
    options.maximum_region_pairs = PositiveInteger(
        Field(config, "maximum_region_pairs"), "maximum_region_pairs",
        100000, true);
    options.input_symmetry_tolerance = Scalar(
        Field(config, "input_symmetry_tolerance"),
        "input_symmetry_tolerance", 1.0e-12, true);
    plhs[0] = Report(radia::beam::PropagateVariationalMap(segments, options));
}

void PropagateGridFunction(std::shared_ptr<ngcomp::GridFunction> field,
                           int nlhs, mxArray* plhs[], int nrhs,
                           const mxArray* prhs[]) {
    if (nrhs != 3 || nlhs != 1)
        throw std::invalid_argument(
            "result = radia_mex('beam.transfer.from_grid_function', "
            "grid_function_handle, config)");
    GridFunctionInput input = ParseGridFunctionInput(prhs[2]);
    plhs[0] = GridFunctionReport(
        radia::beam::PropagateGridFunctionLinearMap(
            field, input.lengths_m, input.positions_m, input.tangents,
            input.names, input.options));
}

}  // namespace

bool DispatchBeamCommand(const std::string& command, int nlhs,
                         mxArray* plhs[], int nrhs,
                         const mxArray* prhs[]) {
    if (command != "beam.transfer.propagate_variational") return false;
    Propagate(nlhs, plhs, nrhs, prhs);
    return true;
}

void BeamTransferFromGridFunction(
        std::shared_ptr<ngcomp::GridFunction> field, int nlhs,
        mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    PropagateGridFunction(std::move(field), nlhs, plhs, nrhs, prhs);
}
