#include "radia_beam_mex_commands.h"

#include "rad_beam_ngsolve.h"
#include "rad_beam_dynamics.h"
#include "rad_beam_transfer.h"
#include "rad_lie_map_batch.h"
#include "rad_lie_map_kernel.h"
#include "rad_hdiv_field_evaluator.h"
#include "rad_orbit_tracker.h"

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
            "beam command requires a scalar configuration struct");
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

bool OptionalBoolean(const mxArray* value, const char* name,
                     bool fallback) {
    if (!value) return fallback;
    if (mxIsLogicalScalar(value)) return mxIsLogicalScalarTrue(value);
    const double number = Scalar(value, name, fallback ? 1.0 : 0.0, false);
    if (number != 0.0 && number != 1.0)
        throw std::invalid_argument(std::string(name) +
                                    " must be logical or 0/1");
    return number != 0.0;
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

radia::beam::GridFunctionMagneticInput ParseGridFunctionMagneticInput(
        const mxArray* value) {
    if (!value)
        return radia::beam::GridFunctionMagneticInput::HCurlVectorPotential;
    const std::string name = Text(value, "field_representation");
    if (name == "magnetic_flux_density")
        return radia::beam::GridFunctionMagneticInput::MagneticFluxDensity;
    if (name == "hcurl_vector_potential")
        return radia::beam::GridFunctionMagneticInput::HCurlVectorPotential;
    throw std::invalid_argument(
        "field_representation must be 'magnetic_flux_density' or "
        "'hcurl_vector_potential'");
}

const char* GridFunctionMagneticInputName(
        radia::beam::GridFunctionMagneticInput value) {
    return value ==
        radia::beam::GridFunctionMagneticInput::HCurlVectorPotential
        ? "hcurl_vector_potential"
        : "magnetic_flux_density";
}

GridFunctionInput ParseGridFunctionInput(const mxArray* config) {
    RequireScalarStruct(config);
    const mxArray* schema = Field(config, "schema");
    const std::string schema_name = schema ? Text(schema, "schema")
        : "radia.beam.grid-function-linear-map.v1";
    const bool multipole = schema_name ==
        "radia.beam.grid-function-multipole-map.v1";
    if (!multipole && schema_name !=
            "radia.beam.grid-function-linear-map.v1")
        throw std::invalid_argument(
            "unsupported beam GridFunction map schema");
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
    output.options.periodic_frame = OptionalBoolean(
        Field(config, "periodic_frame"), "periodic_frame", false);
    output.options.curvature_sign = Scalar(
        Field(config, "curvature_sign"), "curvature_sign", 1.0, true);
    output.options.gradient_sign = Scalar(
        Field(config, "gradient_sign"), "gradient_sign", 1.0, true);
    output.options.multipole_order = static_cast<unsigned>(PositiveInteger(
        Field(config, "multipole_order"), "multipole_order",
        multipole ? 3 : 1, true));
    output.options.maximum_map_order = static_cast<unsigned>(PositiveInteger(
        Field(config, "maximum_map_order"), "maximum_map_order",
        multipole ? 3 : 1, true));
    output.options.maximum_step_m = Scalar(
        Field(config, "maximum_step_m"), "maximum_step_m", 1.0e-3, true);
    output.options.maximum_steps = PositiveInteger(
        Field(config, "maximum_steps"), "maximum_steps", 1000000, true);
    output.options.magnetic_input = ParseGridFunctionMagneticInput(
        Field(config, "field_representation"));
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

mxArray* Tensor5Array(const radia::beam::Tensor5Map6& value) {
    const mwSize dimensions[] = {6, 6, 6, 6, 6};
    mxArray* output = mxCreateNumericArray(5, dimensions, mxDOUBLE_CLASS,
                                           mxREAL);
    double* data = mxGetPr(output);
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            for (std::size_t k = 0; k < 6; ++k)
                for (std::size_t l = 0; l < 6; ++l)
                    for (std::size_t m = 0; m < 6; ++m)
                        data[i + 6 * (j + 6 * (k + 6 * (l + 6 * m)))] =
                            value(i, j, k, l, m);
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

mxArray* CanonicalHamiltonianJetResult(
        const radia::beam::HamiltonianJet6& value) {
    const char* fields[] = {
        "schema", "backend", "coordinate_order", "poisson_pair_signs",
        "reference_beta", "reference_curvature_per_m",
        "field_curvature_per_m", "H2_per_m", "H3_per_m", "H4_per_m",
        "H5_per_m", "A_per_m", "F2_per_m", "F3_per_m", "F4_per_m"};
    mxArray* output = mxCreateStructMatrix(1, 1, 15, fields);
    mxSetField(output, 0, fields[0], mxCreateString(
        "radia.beam.canonical-hamiltonian-jet.result.v1"));
    mxSetField(output, 0, fields[1], mxCreateString("native-cpp-mex"));
    mxSetField(output, 0, fields[2], StringCell(
        {"x", "px_over_p0", "y", "py_over_p0", "ell", "delta"}));
    mxArray* signs = mxCreateDoubleMatrix(1, 3, mxREAL);
    mxGetPr(signs)[0] = 1.0;
    mxGetPr(signs)[1] = 1.0;
    mxGetPr(signs)[2] = -1.0;
    mxSetField(output, 0, fields[3], signs);
    mxSetField(output, 0, fields[4],
               mxCreateDoubleScalar(value.reference_beta));
    mxSetField(output, 0, fields[5],
               mxCreateDoubleScalar(value.reference_curvature_per_m));
    mxSetField(output, 0, fields[6],
               mxCreateDoubleScalar(value.field_curvature_per_m));
    mxSetField(output, 0, fields[7], MatrixArray(value.h2_per_m));
    mxSetField(output, 0, fields[8], Tensor3Array(value.h3_per_m));
    mxSetField(output, 0, fields[9], Tensor4Array(value.h4_per_m));
    mxSetField(output, 0, fields[10], Tensor5Array(value.h5_per_m));
    mxSetField(output, 0, fields[11], MatrixArray(value.dynamics.a_per_m));
    mxSetField(output, 0, fields[12],
               Tensor3Array(value.dynamics.f2_per_m));
    mxSetField(output, 0, fields[13],
               Tensor4Array(value.dynamics.f3_per_m));
    mxSetField(output, 0, fields[14],
               Tensor5Array(value.dynamics.f4_per_m));
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

template <typename Selector>
mxArray* MultipoleRows(std::size_t count, std::size_t columns,
                       Selector selector) {
    mxArray* output = mxCreateDoubleMatrix(count, columns, mxREAL);
    double* data = mxGetPr(output);
    for (std::size_t item = 0; item < count; ++item) {
        const auto& value = selector(item);
        for (std::size_t order = 0; order < columns; ++order)
            data[item + count * order] = value[order];
    }
    return output;
}

mxArray* GridFunctionReport(
        const radia::beam::GridFunctionTransferReport6& value) {
    mxArray* output = Report(value.transfer);
    const bool linear_schema = value.multipole_order == 1 &&
        value.transfer.maximum_order == 1;
    ReplaceStructField(output, "schema", mxCreateString(
        linear_schema
            ? "radia.beam.grid-function-linear-map.result.v1"
            : "radia.beam.grid-function-multipole-map.result.v1"));
    ReplaceStructField(output, "backend", mxCreateString(
        "native-cpp-ngsolve-gridfunction-mex"));
    AddStructField(output, "field_source",
                   mxCreateString("ngsolve.GridFunction"));
    AddStructField(output, "field_representation", mxCreateString(
        GridFunctionMagneticInputName(value.magnetic_input)));
    AddStructField(output, "magnetic_evaluation", mxCreateString(
        value.magnetic_input ==
            radia::beam::GridFunctionMagneticInput::HCurlVectorPotential
            ? "ngsolve-native-curl(A)"
            : "direct-B"));
    AddStructField(output, "grid_function_space_class", mxCreateString(
        value.grid_function_space_class.c_str()));
    AddStructField(output, "grid_function_space_order",
        mxCreateDoubleScalar(value.grid_function_space_order));
    AddStructField(output, "linearization_order", mxCreateDoubleScalar(
        static_cast<double>(value.multipole_order)));
    AddStructField(output, "magnetic_rigidity_t_m",
                   mxCreateDoubleScalar(value.magnetic_rigidity_t_m));
    AddStructField(output, "sample_radius_m",
                   mxCreateDoubleScalar(value.sample_radius_m));
    const std::string fit_model =
        "nine-point transverse harmonic multipole expansion through order " +
        std::to_string(value.multipole_order);
    AddStructField(output, "fit_model", mxCreateString(fit_model.c_str()));
    AddStructField(output, "frame_convention", mxCreateString(
        value.periodic_frame
            ? "right-handed periodic minimal-twist frame from Bishop/RMF "
              "holonomy correction"
            : "right-handed Bishop/RMF double reflection seeded by "
              "initial_horizontal"));
    AddStructField(output, "periodic_frame",
                   mxCreateLogicalScalar(value.periodic_frame));
    AddStructField(output, "frame_holonomy_correction_rad",
                   mxCreateDoubleScalar(
                       value.frame_holonomy_correction_rad));

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
    AddStructField(output, "multipole_normal_t_per_m_power", MultipoleRows(
        count, value.multipole_order + 1,
        [&](std::size_t index) -> const auto& {
            return samples[index].multipoles.normal_t_per_m_power;
        }));
    AddStructField(output, "multipole_skew_t_per_m_power", MultipoleRows(
        count, value.multipole_order + 1,
        [&](std::size_t index) -> const auto& {
            return samples[index].multipoles.skew_t_per_m_power;
        }));
    AddStructField(output, "multipole_rms_fit_residual_t", ScalarColumn(
        count, [&](std::size_t index) {
            return samples[index].multipole_rms_fit_residual_t;
        }));
    AddStructField(output, "multipole_maximum_fit_residual_t", ScalarColumn(
        count, [&](std::size_t index) {
            return samples[index].multipole_maximum_fit_residual_t;
        }));
    AddStructField(output, "multipole_fit_rank", ScalarColumn(
        count, [&](std::size_t index) {
            return static_cast<double>(samples[index].multipole_fit_rank);
        }));
    AddStructField(output, "multipole_scaled_design_condition", ScalarColumn(
        count, [&](std::size_t index) {
            return samples[index].multipole_scaled_design_condition;
        }));
    AddStructField(output, "local_A_per_m", MatrixStack(
        count, [&](std::size_t index) -> const auto& {
            return samples[index].a_per_m;
        }));
    AddStructField(output, "local_F2_per_m", Tensor3Stack(
        count, [&](std::size_t index) -> const auto& {
            return samples[index].dynamics_jet.f2_per_m;
        }));
    AddStructField(output, "local_F3_per_m", Tensor4Stack(
        count, [&](std::size_t index) -> const auto& {
            return samples[index].dynamics_jet.f3_per_m;
        }));
    return output;
}

radia::beam::Vec3 BeamVec3(const mxArray* value, const char* name) {
    const double* data = RealData(value, name);
    if (mxGetNumberOfElements(value) != 3)
        throw std::invalid_argument(std::string(name) +
                                    " must contain three entries");
    return {data[0], data[1], data[2]};
}

mxArray* BeamVec3Array(const radia::beam::Vec3& value) {
    mxArray* output = mxCreateDoubleMatrix(1, 3, mxREAL);
    double* data = mxGetPr(output);
    data[0] = value.x;
    data[1] = value.y;
    data[2] = value.z;
    return output;
}

radia::beam::ParticleSpecies ParseSpecies(const mxArray* value) {
    RequireScalarStruct(value);
    radia::beam::ParticleSpecies output;
    output.charge_c = Scalar(
        Field(value, "charge_c"), "species.charge_c", 0.0, false);
    output.rest_mass_kg = Scalar(
        Field(value, "rest_mass_kg"), "species.rest_mass_kg", 0.0,
        false);
    const mxArray* name = Field(value, "name");
    output.name = name ? Text(name, "species.name") : "custom";
    // Reuse the C++ factory as the canonical validation boundary.
    (void)radia::beam::ReferenceParticle::FromKineticEnergyEV(output, 0.0);
    return output;
}

radia::beam::CartesianState ParseState(const mxArray* value) {
    RequireScalarStruct(value);
    radia::beam::CartesianState output;
    output.position_m = BeamVec3(
        Field(value, "position_m"), "state.position_m");
    output.kinetic_momentum_kg_m_s = BeamVec3(
        Field(value, "kinetic_momentum_kg_m_s"),
        "state.kinetic_momentum_kg_m_s");
    output.time_s = Scalar(
        Field(value, "time_s"), "state.time_s", 0.0, true);
    output.path_length_m = Scalar(
        Field(value, "path_length_m"), "state.path_length_m", 0.0,
        true);
    return output;
}

radia::beam::IndependentVariable ParseIndependent(const mxArray* value) {
    const std::string name = value ? Text(value, "independent") : "time";
    if (name == "time") return radia::beam::IndependentVariable::time;
    if (name == "path_length")
        return radia::beam::IndependentVariable::path_length;
    if (name == "azimuth")
        return radia::beam::IndependentVariable::azimuth;
    throw std::invalid_argument(
        "independent must be 'time', 'path_length', or 'azimuth'");
}

std::shared_ptr<radia::beam::Field> ParseFieldObject(const mxArray* value) {
    RequireScalarStruct(value);
    const mxArray* type_value = Field(value, "type");
    const std::string type = type_value ? Text(type_value, "field.type")
                                        : "uniform";
    if (type == "zero") return std::make_shared<radia::beam::ZeroField>();
    if (type != "uniform")
        throw std::invalid_argument(
            "field.type must be 'zero' or 'uniform'");
    const mxArray* magnetic = Field(value, "magnetic_t");
    const mxArray* electric = Field(value, "electric_v_m");
    return std::make_shared<radia::beam::UniformField>(
        magnetic ? BeamVec3(magnetic, "field.magnetic_t")
                 : radia::beam::Vec3{},
        electric ? BeamVec3(electric, "field.electric_v_m")
                 : radia::beam::Vec3{});
}

mxArray* SpeciesStruct(const radia::beam::ParticleSpecies& value) {
    const char* fields[] = {"charge_c", "rest_mass_kg", "name"};
    mxArray* output = mxCreateStructMatrix(1, 1, 3, fields);
    mxSetField(output, 0, fields[0], mxCreateDoubleScalar(value.charge_c));
    mxSetField(output, 0, fields[1],
               mxCreateDoubleScalar(value.rest_mass_kg));
    mxSetField(output, 0, fields[2], mxCreateString(value.name.c_str()));
    return output;
}

mxArray* ReferenceParticleStruct(
        const radia::beam::ReferenceParticle& value) {
    const char* fields[] = {
        "species", "kinetic_energy_j", "momentum_kg_m_s",
        "magnetic_rigidity_t_m", "backend"};
    mxArray* output = mxCreateStructMatrix(1, 1, 5, fields);
    mxSetField(output, 0, fields[0], SpeciesStruct(value.species));
    mxSetField(output, 0, fields[1],
               mxCreateDoubleScalar(value.kinetic_energy_j));
    mxSetField(output, 0, fields[2],
               mxCreateDoubleScalar(value.momentum_kg_m_s));
    mxSetField(output, 0, fields[3],
               mxCreateDoubleScalar(value.magnetic_rigidity_t_m));
    mxSetField(output, 0, fields[4], mxCreateString("native-cpp-mex"));
    return output;
}

mxArray* StateStruct(const radia::beam::CartesianState& value) {
    const char* fields[] = {
        "position_m", "kinetic_momentum_kg_m_s", "time_s",
        "path_length_m"};
    mxArray* output = mxCreateStructMatrix(1, 1, 4, fields);
    mxSetField(output, 0, fields[0], BeamVec3Array(value.position_m));
    mxSetField(output, 0, fields[1],
               BeamVec3Array(value.kinetic_momentum_kg_m_s));
    mxSetField(output, 0, fields[2], mxCreateDoubleScalar(value.time_s));
    mxSetField(output, 0, fields[3],
               mxCreateDoubleScalar(value.path_length_m));
    return output;
}

const char* DomainStatusText(radia::beam::DomainStatus value) {
    switch (value) {
    case radia::beam::DomainStatus::inside:
        return "inside";
    case radia::beam::DomainStatus::outside:
        return "outside";
    case radia::beam::DomainStatus::boundary:
        return "boundary";
    case radia::beam::DomainStatus::invalid:
        return "invalid";
    }
    return "invalid";
}

mxArray* FieldSampleStruct(const radia::beam::FieldSample& value) {
    const char* fields[] = {"electric_v_m", "magnetic_t", "domain_status"};
    mxArray* output = mxCreateStructMatrix(1, 1, 3, fields);
    mxSetField(output, 0, fields[0], BeamVec3Array(value.electric_v_m));
    mxSetField(output, 0, fields[1], BeamVec3Array(value.magnetic_t));
    mxSetField(output, 0, fields[2],
               mxCreateString(DomainStatusText(value.domain_status)));
    return output;
}

mxArray* DerivativeStruct(const radia::beam::StateDerivative& value) {
    const char* fields[] = {
        "dposition_m", "dkinetic_momentum_kg_m_s", "dtime_s",
        "dpath_length_m", "field"};
    mxArray* output = mxCreateStructMatrix(1, 1, 5, fields);
    mxSetField(output, 0, fields[0], BeamVec3Array(value.dposition_m));
    mxSetField(output, 0, fields[1],
               BeamVec3Array(value.dkinetic_momentum_kg_m_s));
    mxSetField(output, 0, fields[2], mxCreateDoubleScalar(value.dtime_s));
    mxSetField(output, 0, fields[3],
               mxCreateDoubleScalar(value.dpath_length_m));
    mxSetField(output, 0, fields[4], FieldSampleStruct(value.field));
    return output;
}

mxArray* InvariantStruct(const radia::beam::InvariantReport& value) {
    const char* fields[] = {
        "momentum_kg_m_s", "relativistic_gamma", "kinetic_energy_j",
        "speed_m_s", "domain_status"};
    mxArray* output = mxCreateStructMatrix(1, 1, 5, fields);
    mxSetField(output, 0, fields[0],
               mxCreateDoubleScalar(value.momentum_kg_m_s));
    mxSetField(output, 0, fields[1],
               mxCreateDoubleScalar(value.relativistic_gamma));
    mxSetField(output, 0, fields[2],
               mxCreateDoubleScalar(value.kinetic_energy_j));
    mxSetField(output, 0, fields[3],
               mxCreateDoubleScalar(value.speed_m_s));
    mxSetField(output, 0, fields[4],
               mxCreateString(DomainStatusText(value.domain_status)));
    return output;
}

mxArray* StepStruct(const radia::beam::StepResult& value) {
    const char* fields[] = {
        "independent_before", "independent_after", "accepted_step",
        "state_before", "state_after", "rhs_before",
        "invariants_before", "invariants_after", "backend"};
    mxArray* output = mxCreateStructMatrix(1, 1, 9, fields);
    mxSetField(output, 0, fields[0],
               mxCreateDoubleScalar(value.independent_before));
    mxSetField(output, 0, fields[1],
               mxCreateDoubleScalar(value.independent_after));
    mxSetField(output, 0, fields[2],
               mxCreateDoubleScalar(value.accepted_step));
    mxSetField(output, 0, fields[3], StateStruct(value.state_before));
    mxSetField(output, 0, fields[4], StateStruct(value.state_after));
    mxSetField(output, 0, fields[5], DerivativeStruct(value.rhs_before));
    mxSetField(output, 0, fields[6],
               InvariantStruct(value.invariants_before));
    mxSetField(output, 0, fields[7],
               InvariantStruct(value.invariants_after));
    mxSetField(output, 0, fields[8], mxCreateString("native-cpp-mex"));
    return output;
}

mxArray* StepRecordStruct(const radia::beam::StepRecord& value) {
    const char* fields[] = {
        "independent_value", "attempted_step", "accepted_step",
        "accepted", "state_before", "state_after", "rhs_before",
        "invariants_before", "invariants_after"};
    mxArray* output = mxCreateStructMatrix(1, 1, 9, fields);
    mxSetField(output, 0, fields[0],
               mxCreateDoubleScalar(value.independent_value));
    mxSetField(output, 0, fields[1],
               mxCreateDoubleScalar(value.attempted_step));
    mxSetField(output, 0, fields[2],
               mxCreateDoubleScalar(value.accepted_step));
    mxSetField(output, 0, fields[3], mxCreateLogicalScalar(value.accepted));
    mxSetField(output, 0, fields[4], StateStruct(value.state_before));
    mxSetField(output, 0, fields[5], StateStruct(value.state_after));
    mxSetField(output, 0, fields[6], DerivativeStruct(value.rhs_before));
    mxSetField(output, 0, fields[7],
               InvariantStruct(value.invariants_before));
    mxSetField(output, 0, fields[8],
               InvariantStruct(value.invariants_after));
    return output;
}

mxArray* TrajectoryStruct(const radia::beam::Trajectory& value) {
    const char* fields[] = {"schema", "backend", "samples", "steps",
                            "summary"};
    mxArray* output = mxCreateStructMatrix(1, 1, 5, fields);
    mxSetField(output, 0, fields[0],
               mxCreateString("radia.beam.trajectory.result.v1"));
    mxSetField(output, 0, fields[1], mxCreateString("native-cpp-mex"));

    mxArray* samples = mxCreateCellMatrix(value.Samples().size(), 1);
    for (std::size_t index = 0; index < value.Samples().size(); ++index)
        mxSetCell(samples, index, StateStruct(value.Samples()[index]));
    mxSetField(output, 0, fields[2], samples);

    mxArray* steps = mxCreateCellMatrix(value.Steps().size(), 1);
    for (std::size_t index = 0; index < value.Steps().size(); ++index)
        mxSetCell(steps, index, StepRecordStruct(value.Steps()[index]));
    mxSetField(output, 0, fields[3], steps);

    const auto& summary = value.Summary();
    const char* summary_fields[] = {
        "accepted_steps", "independent_start", "independent_stop",
        "path_length_change_m", "momentum_conservation_applicable",
        "maximum_relative_momentum_error"};
    mxArray* summary_output = mxCreateStructMatrix(
        1, 1, 6, summary_fields);
    mxSetField(summary_output, 0, summary_fields[0],
               mxCreateDoubleScalar(
                   static_cast<double>(summary.accepted_steps)));
    mxSetField(summary_output, 0, summary_fields[1],
               mxCreateDoubleScalar(summary.independent_start));
    mxSetField(summary_output, 0, summary_fields[2],
               mxCreateDoubleScalar(summary.independent_stop));
    mxSetField(summary_output, 0, summary_fields[3],
               mxCreateDoubleScalar(summary.path_length_change_m));
    mxSetField(summary_output, 0, summary_fields[4],
               mxCreateLogicalScalar(
                   summary.momentum_conservation_applicable));
    mxSetField(summary_output, 0, summary_fields[5],
               mxCreateDoubleScalar(
                   summary.maximum_relative_momentum_error));
    mxSetField(output, 0, fields[4], summary_output);
    return output;
}

struct TrackingInput {
    radia::beam::ParticleSpecies species;
    radia::beam::CartesianState state;
    std::shared_ptr<radia::beam::Field> field;
    radia::beam::IndependentVariable independent;
};

TrackingInput ParseTrackingInput(const mxArray* config) {
    RequireScalarStruct(config);
    const mxArray* schema = Field(config, "schema");
    if (schema && Text(schema, "schema") != "radia.beam.tracking.v1")
        throw std::invalid_argument("unsupported beam tracking schema");
    TrackingInput output;
    output.species = ParseSpecies(Field(config, "species"));
    output.state = ParseState(Field(config, "state"));
    output.field = ParseFieldObject(Field(config, "field"));
    output.independent = ParseIndependent(Field(config, "independent"));
    return output;
}

std::shared_ptr<radia::beam::Stepper> ParseStepper(const mxArray* value) {
    const std::string name = value ? Text(value, "stepper") : "boris2";
    if (name == "boris2") return std::make_shared<radia::beam::Boris2>();
    if (name == "classical-rk4")
        return std::make_shared<radia::beam::ClassicalRK4>();
    throw std::invalid_argument(
        "stepper must be 'boris2' or 'classical-rk4'");
}

void ReferenceParticle(int nlhs, mxArray* plhs[], int nrhs,
                       const mxArray* prhs[]) {
    if (nrhs != 3 || nlhs != 1)
        throw std::invalid_argument(
            "reference = radia_mex('beam.reference_particle.from_kinetic_energy_ev', species, kinetic_energy_ev)");
    const auto species = ParseSpecies(prhs[1]);
    const double energy = Scalar(
        prhs[2], "kinetic_energy_ev", 0.0, false);
    plhs[0] = ReferenceParticleStruct(
        radia::beam::ReferenceParticle::FromKineticEnergyEV(
            species, energy));
}

void ParticleSpeciesPreset(const std::string& command, int nlhs,
                           mxArray* plhs[], int nrhs) {
    if (nrhs != 1 || nlhs != 1)
        throw std::invalid_argument(
            "species = radia_mex('beam.species.proton|electron')");
    plhs[0] = SpeciesStruct(
        command == "beam.species.proton"
            ? radia::beam::ParticleSpecies::Proton()
            : radia::beam::ParticleSpecies::Electron());
}

void BeamFieldSample(int nlhs, mxArray* plhs[], int nrhs,
                     const mxArray* prhs[]) {
    if (nrhs != 3 && nrhs != 4)
        throw std::invalid_argument(
            "sample = radia_mex('beam.field.sample', field, position_m [, time_s])");
    if (nlhs != 1)
        throw std::invalid_argument("beam.field.sample returns one result");
    auto field = ParseFieldObject(prhs[1]);
    const auto position = BeamVec3(prhs[2], "position_m");
    const double time = nrhs == 4
        ? Scalar(prhs[3], "time_s", 0.0, false) : 0.0;
    plhs[0] = FieldSampleStruct(field->Evaluate(position, time));
}

void LorentzRHS(int nlhs, mxArray* plhs[], int nrhs,
                const mxArray* prhs[]) {
    if (nrhs != 2 || nlhs != 1)
        throw std::invalid_argument(
            "rhs = radia_mex('beam.equation.rhs', config)");
    TrackingInput input = ParseTrackingInput(prhs[1]);
    radia::beam::LorentzEquation equation(
        input.species, input.field, input.independent);
    const double independent_value = Scalar(
        Field(prhs[1], "independent_value"), "independent_value", 0.0,
        true);
    plhs[0] = DerivativeStruct(
        equation.RHS(independent_value, input.state));
}

void BeamStep(int nlhs, mxArray* plhs[], int nrhs,
              const mxArray* prhs[]) {
    if (nrhs != 2 || nlhs != 1)
        throw std::invalid_argument(
            "result = radia_mex('beam.step', config)");
    TrackingInput input = ParseTrackingInput(prhs[1]);
    auto equation = std::make_shared<radia::beam::LorentzEquation>(
        input.species, input.field, input.independent);
    auto stepper = ParseStepper(Field(prhs[1], "stepper"));
    const double independent_value = Scalar(
        Field(prhs[1], "independent_value"), "independent_value", 0.0,
        true);
    const double step = Scalar(Field(prhs[1], "step"), "step", 0.0,
                               false);
    plhs[0] = StepStruct(
        stepper->Step(*equation, independent_value, input.state, step));
}

void BeamTrack(int nlhs, mxArray* plhs[], int nrhs,
               const mxArray* prhs[]) {
    if (nrhs != 2 || nlhs != 1)
        throw std::invalid_argument(
            "trajectory = radia_mex('beam.track', config)");
    TrackingInput input = ParseTrackingInput(prhs[1]);
    auto equation = std::make_shared<radia::beam::LorentzEquation>(
        input.species, input.field, input.independent);
    auto stepper = ParseStepper(Field(prhs[1], "stepper"));
    radia::beam::TrackPlan plan;
    plan.start = Scalar(Field(prhs[1], "start"), "start", 0.0, false);
    plan.stop = Scalar(Field(prhs[1], "stop"), "stop", 0.0, false);
    plan.maximum_step = Scalar(
        Field(prhs[1], "maximum_step"), "maximum_step", 0.0, false);
    plan.maximum_steps = PositiveInteger(
        Field(prhs[1], "maximum_steps"), "maximum_steps", 1000000,
        true);
    plhs[0] = TrajectoryStruct(
        radia::beam::Tracker(equation, stepper).Track(input.state, plan));
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

void CanonicalHamiltonianJet(int nlhs, mxArray* plhs[], int nrhs,
                             const mxArray* prhs[]) {
    if (nrhs != 2 || nlhs != 1)
        throw std::invalid_argument(
            "result = radia_mex('beam.hamiltonian.canonical_body_jet', config)");
    const mxArray* config = prhs[1];
    RequireScalarStruct(config);
    const mxArray* schema = Field(config, "schema");
    if (schema && Text(schema, "schema") !=
            "radia.beam.canonical-hamiltonian-jet.v1")
        throw std::invalid_argument(
            "unsupported beam canonical-Hamiltonian schema");
    const mxArray* coefficients_value = Field(config, "coefficients");
    const double* coefficients = RealData(
        coefficients_value, "coefficients");
    const std::size_t coefficient_count =
        mxGetNumberOfElements(coefficients_value);
    if (coefficient_count != 7 && coefficient_count != 9)
        throw std::invalid_argument(
            "coefficients must contain seven or nine entries");
    radia::beam::TransverseMagneticMultipoleExpansion expansion;
    expansion.order = coefficient_count == 9 ? 4 : 3;
    expansion.normal_t_per_m_power = {
        coefficients[0], coefficients[1], coefficients[3], coefficients[5],
        coefficient_count == 9 ? coefficients[7] : 0.0};
    expansion.skew_t_per_m_power = {
        0.0, coefficients[2], coefficients[4], coefficients[6],
        coefficient_count == 9 ? coefficients[8] : 0.0};
    const double rigidity = Scalar(
        Field(config, "magnetic_rigidity_t_m"),
        "magnetic_rigidity_t_m", 0.0, false);
    const double curvature_sign = Scalar(
        Field(config, "curvature_sign"), "curvature_sign", 1.0, true);
    const double gradient_sign = Scalar(
        Field(config, "gradient_sign"), "gradient_sign", 1.0, true);
    const double reference_beta = Scalar(
        Field(config, "reference_beta"), "reference_beta", 1.0, true);
    std::optional<double> reference_curvature_per_m;
    if (const mxArray* value = Field(config, "reference_curvature_per_m"))
        reference_curvature_per_m = Scalar(
            value, "reference_curvature_per_m", 0.0, false);
    plhs[0] = CanonicalHamiltonianJetResult(
        radia::beam::BuildCanonicalBodyHamiltonianJet(
            expansion, rigidity, curvature_sign, gradient_sign,
            reference_beta, reference_curvature_per_m));
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
        radia::beam::PropagateGridFunctionMultipoleMap(
            field, input.lengths_m, input.positions_m, input.tangents,
            input.names, input.options));
}

// ---------------------------------------------------------------------------
// Lie-map construction / factorization / tracking commands over the native
// rad_lie / rad_orbit kernels.  MATLAB stores column-major, the kernels are
// row-major; rank-k tensors are converted index by index, while (6,N) MATLAB
// blocks are byte-identical to (N,6) row-major buffers and pass through.
// ---------------------------------------------------------------------------

int Pow6(int rank) {
    int result = 1;
    for (int axis = 0; axis < rank; ++axis) result *= 6;
    return result;
}

void ReadRowMajorTensor(const mxArray* value, const char* name, int rank,
                        double* out) {
    RequireLeadingShape(value, name, static_cast<std::size_t>(rank), 1);
    const double* data = mxGetPr(value);
    const int size = Pow6(rank);
    int digits[5];
    for (int flat = 0; flat < size; ++flat) {
        int remainder = flat;
        for (int axis = rank-1; axis >= 0; --axis) {
            digits[axis] = remainder % 6;
            remainder /= 6;
        }
        int matlab_flat = 0;
        for (int axis = rank-1; axis >= 0; --axis)
            matlab_flat = matlab_flat * 6 + digits[axis];
        out[flat] = data[matlab_flat];
    }
}

mxArray* RowMajorTensorArray(const double* value, int rank) {
    const mwSize dimensions[5] = {6, 6, 6, 6, 6};
    mxArray* output = mxCreateNumericArray(
        static_cast<mwSize>(rank), dimensions, mxDOUBLE_CLASS, mxREAL);
    double* data = mxGetPr(output);
    const int size = Pow6(rank);
    int digits[5];
    for (int flat = 0; flat < size; ++flat) {
        int remainder = flat;
        for (int axis = rank-1; axis >= 0; --axis) {
            digits[axis] = remainder % 6;
            remainder /= 6;
        }
        int matlab_flat = 0;
        for (int axis = rank-1; axis >= 0; --axis)
            matlab_flat = matlab_flat * 6 + digits[axis];
        data[matlab_flat] = value[flat];
    }
    return output;
}

mxArray* RowVector4(const double values[4]) {
    mxArray* output = mxCreateDoubleMatrix(1, 4, mxREAL);
    double* data = mxGetPr(output);
    for (int index = 0; index < 4; ++index) data[index] = values[index];
    return output;
}

void LieMapTensorsSpoly(int nlhs, mxArray* plhs[], int nrhs,
                        const mxArray* prhs[]) {
    if (nrhs != 2 || nlhs != 1)
        throw std::invalid_argument(
            "result = radia_mex('beam.lie.map_tensors_spoly', config)");
    const mxArray* config = prhs[1];
    RequireScalarStruct(config);
    const mxArray* schema = Field(config, "schema");
    if (schema && Text(schema, "schema") != "radia.beam.lie-map-spoly.v1")
        throw std::invalid_argument("unsupported Lie-map spoly schema");
    const mxArray* Ay_value = Field(config, "Ay_t_m");
    const mxArray* As_value = Field(config, "As_t_m");
    if (!Ay_value || !As_value)
        throw std::invalid_argument(
            "Ay_t_m/As_t_m coefficient arrays are required");
    (void)RealData(Ay_value, "Ay_t_m");
    (void)RealData(As_value, "As_t_m");
    const mwSize rank = mxGetNumberOfDimensions(Ay_value);
    const mwSize* dimensions = mxGetDimensions(Ay_value);
    // MATLAB shape (d+1, d+1, k_s+1, n_segments); trailing singleton axes
    // are dropped by MATLAB, so rank 2..4 are all accepted.
    if (rank < 2 || rank > 4 || dimensions[0] != dimensions[1])
        throw std::invalid_argument(
            "Ay_t_m must have shape (d+1, d+1, k_s+1, n_segments)");
    const std::size_t transverse = dimensions[0];
    const std::size_t s_orders = rank >= 3 ? dimensions[2] : 1;
    const std::size_t segments = rank >= 4 ? dimensions[3] : 1;
    if (mxGetNumberOfDimensions(As_value) != rank
        || mxGetNumberOfElements(As_value)
            != mxGetNumberOfElements(Ay_value))
        throw std::invalid_argument("As_t_m must match the Ay_t_m shape");
    const std::size_t block = transverse * transverse;
    std::vector<double> Ay(segments * s_orders * block);
    std::vector<double> As(segments * s_orders * block);
    const double* Ay_data = mxGetPr(Ay_value);
    const double* As_data = mxGetPr(As_value);
    for (std::size_t segment = 0; segment < segments; ++segment)
        for (std::size_t order = 0; order < s_orders; ++order)
            for (std::size_t x_power = 0; x_power < transverse; ++x_power)
                for (std::size_t y_power = 0; y_power < transverse;
                     ++y_power) {
                    const std::size_t matlab_flat = x_power
                        + transverse * (y_power
                        + transverse * (order + s_orders * segment));
                    const std::size_t native_flat =
                        (segment * s_orders + order) * block
                        + x_power * transverse + y_power;
                    Ay[native_flat] = Ay_data[matlab_flat];
                    As[native_flat] = As_data[matlab_flat];
                }
    const mxArray* lengths_value = Field(config, "lengths_m");
    const double* lengths = RealData(lengths_value, "lengths_m");
    if (mxGetNumberOfElements(lengths_value) != segments)
        throw std::invalid_argument(
            "lengths_m must contain one entry per segment");
    const mxArray* curvature_value = Field(config, "curvature_per_m");
    const double* curvature_data = RealData(
        curvature_value, "curvature_per_m");
    std::vector<double> curvature;
    std::size_t curvature_columns = 1;
    if (mxGetNumberOfElements(curvature_value) == segments) {
        curvature.assign(curvature_data, curvature_data + segments);
    } else {
        const mwSize* curvature_dims = mxGetDimensions(curvature_value);
        if (mxGetNumberOfDimensions(curvature_value) != 2
            || static_cast<std::size_t>(curvature_dims[1]) != segments)
            throw std::invalid_argument(
                "curvature_per_m must be (n_segments) or "
                "(k_c+1, n_segments)");
        curvature_columns = curvature_dims[0];
        curvature.resize(segments * curvature_columns);
        for (std::size_t segment = 0; segment < segments; ++segment)
            for (std::size_t column = 0; column < curvature_columns;
                 ++column)
                curvature[segment*curvature_columns + column] =
                    curvature_data[column + curvature_columns*segment];
    }
    const double rigidity = Scalar(
        Field(config, "magnetic_rigidity_t_m"),
        "magnetic_rigidity_t_m", 0.0, false);
    const double curvature_sign = Scalar(
        Field(config, "curvature_sign"), "curvature_sign", 1.0, true);
    const double reference_beta = Scalar(
        Field(config, "reference_beta"), "reference_beta", 1.0, true);
    const bool covariant = OptionalBoolean(
        Field(config, "longitudinal_covariant"), "longitudinal_covariant",
        true);
    const double maximum_step_m = Scalar(
        Field(config, "maximum_step_m"), "maximum_step_m", 1.0e-3, true);
    const long long maximum_steps = static_cast<long long>(PositiveInteger(
        Field(config, "maximum_steps"), "maximum_steps", 1000000, true));
    const double orbit_tolerance = Scalar(
        Field(config, "reference_orbit_tolerance"),
        "reference_orbit_tolerance", 1.0e-8, true);
    double poisson[36];
    rad_lie::CanonicalPoissonMatrix6(poisson);
    std::vector<double> R(36), T(216), U(1296), V(7776);
    mxArray* linear = mxCreateDoubleMatrix(6, segments, mxREAL);
    double worst_linear = 0.0;
    rad_lie::LieMapTensorsFromSpolyArrays(
        Ay.data(), As.data(), segments, s_orders, transverse,
        lengths, curvature.data(), curvature_columns,
        rigidity, curvature_sign, reference_beta, covariant ? 1 : 0,
        poisson, maximum_step_m, maximum_steps, orbit_tolerance,
        R.data(), T.data(), U.data(), V.data(), mxGetPr(linear),
        &worst_linear);
    const char* names[] = {"R", "T", "U", "V", "hamiltonian_linear",
                           "worst_hamiltonian_linear"};
    mxArray* output = mxCreateStructMatrix(1, 1, 6, names);
    mxSetField(output, 0, "R", RowMajorTensorArray(R.data(), 2));
    mxSetField(output, 0, "T", RowMajorTensorArray(T.data(), 3));
    mxSetField(output, 0, "U", RowMajorTensorArray(U.data(), 4));
    mxSetField(output, 0, "V", RowMajorTensorArray(V.data(), 5));
    mxSetField(output, 0, "hamiltonian_linear", linear);
    mxSetField(output, 0, "worst_hamiltonian_linear",
               mxCreateDoubleScalar(worst_linear));
    plhs[0] = output;
}

void LieDragtFinnFactorize(int nlhs, mxArray* plhs[], int nrhs,
                           const mxArray* prhs[]) {
    if (nrhs != 2 || nlhs != 1)
        throw std::invalid_argument(
            "result = radia_mex('beam.lie.dragt_finn_factorize', config)");
    const mxArray* config = prhs[1];
    RequireScalarStruct(config);
    std::vector<double> R(36), T(216), U(1296), V(7776);
    ReadRowMajorTensor(Field(config, "R"), "R", 2, R.data());
    ReadRowMajorTensor(Field(config, "T"), "T", 3, T.data());
    ReadRowMajorTensor(Field(config, "U"), "U", 4, U.data());
    ReadRowMajorTensor(Field(config, "V"), "V", 5, V.data());
    double poisson[36];
    rad_lie::CanonicalPoissonMatrix6(poisson);
    std::vector<double> f3(216), f4(1296), f5(7776);
    std::vector<double> T_recon(216), U_recon(1296), V_recon(7776);
    rad_lie::DragtFinnDiagnostics diagnostics;
    rad_lie::DragtFinnFactorizeFourthOrder(
        R.data(), T.data(), U.data(), V.data(), poisson,
        f3.data(), f4.data(), f5.data(),
        T_recon.data(), U_recon.data(), V_recon.data(), &diagnostics);
    const char* names[] = {
        "R", "f3", "f4", "f5", "T", "U", "V",
        "f3_symmetry_defect", "f4_symmetry_defect", "f5_symmetry_defect",
        "relative_reconstruction_error", "raw_residual",
        "reconstructed_residual"};
    mxArray* output = mxCreateStructMatrix(1, 1, 13, names);
    mxSetField(output, 0, "R", RowMajorTensorArray(R.data(), 2));
    mxSetField(output, 0, "f3", RowMajorTensorArray(f3.data(), 3));
    mxSetField(output, 0, "f4", RowMajorTensorArray(f4.data(), 4));
    mxSetField(output, 0, "f5", RowMajorTensorArray(f5.data(), 5));
    mxSetField(output, 0, "T", RowMajorTensorArray(T_recon.data(), 3));
    mxSetField(output, 0, "U", RowMajorTensorArray(U_recon.data(), 4));
    mxSetField(output, 0, "V", RowMajorTensorArray(V_recon.data(), 5));
    mxSetField(output, 0, "f3_symmetry_defect",
               mxCreateDoubleScalar(diagnostics.f3_symmetry_defect));
    mxSetField(output, 0, "f4_symmetry_defect",
               mxCreateDoubleScalar(diagnostics.f4_symmetry_defect));
    mxSetField(output, 0, "f5_symmetry_defect",
               mxCreateDoubleScalar(diagnostics.f5_symmetry_defect));
    mxSetField(output, 0, "relative_reconstruction_error",
               mxCreateDoubleScalar(
                   diagnostics.relative_reconstruction_error));
    mxSetField(output, 0, "raw_residual",
               RowVector4(diagnostics.raw_residual));
    mxSetField(output, 0, "reconstructed_residual",
               RowVector4(diagnostics.reconstructed_residual));
    plhs[0] = output;
}

void LieApplyDragtFinnBatch(int nlhs, mxArray* plhs[], int nrhs,
                            const mxArray* prhs[]) {
    if (nrhs != 2 || nlhs != 1)
        throw std::invalid_argument(
            "states = radia_mex('beam.lie.apply_dragt_finn_batch', config)");
    const mxArray* config = prhs[1];
    RequireScalarStruct(config);
    std::vector<double> R(36), f3(216), f4(1296), f5(7776);
    ReadRowMajorTensor(Field(config, "R"), "R", 2, R.data());
    ReadRowMajorTensor(Field(config, "f3"), "f3", 3, f3.data());
    ReadRowMajorTensor(Field(config, "f4"), "f4", 4, f4.data());
    const mxArray* f5_value = Field(config, "f5");
    const bool has_f5 = f5_value != nullptr;
    if (has_f5)
        ReadRowMajorTensor(f5_value, "f5", 5, f5.data());
    const mxArray* states_value = Field(config, "states");
    (void)RealData(states_value, "states");
    if (mxGetNumberOfDimensions(states_value) != 2
        || mxGetDimensions(states_value)[0] != 6)
        throw std::invalid_argument("states must have shape (6, n_states)");
    const std::size_t count = mxGetDimensions(states_value)[1];
    const std::size_t substeps = PositiveInteger(
        Field(config, "generator_substeps"), "generator_substeps", 1, true);
    const double tolerance = Scalar(
        Field(config, "newton_tolerance"), "newton_tolerance", 1.0e-13,
        true);
    const std::size_t iteration_cap = PositiveInteger(
        Field(config, "maximum_newton_iterations"),
        "maximum_newton_iterations", 20, true);
    double poisson[36];
    rad_lie::CanonicalPoissonMatrix6(poisson);
    // A (6, N) MATLAB block is byte-identical to the (N, 6) row-major
    // ensemble the batch kernel consumes.
    mxArray* output = mxCreateDoubleMatrix(6, count, mxREAL);
    rad_lie::ApplyDragtFinnMapBatch(
        R.data(), f3.data(), f4.data(), has_f5 ? f5.data() : nullptr,
        poisson, mxGetPr(states_value), mxGetPr(output), count,
        static_cast<int>(substeps), tolerance,
        static_cast<int>(iteration_cap));
    plhs[0] = output;
}

mxArray* OrbitTrackResultStruct(
        mxArray* positions, mxArray* tangents, mxArray* stations,
        mxArray* curvature, const rad_orbit::OrbitTrackResult& tracked) {
    const char* names[] = {
        "positions_m", "tangents", "stations_m", "signed_curvature_per_m",
        "length_m", "out_of_plane_m", "out_of_plane_slope"};
    mxArray* output = mxCreateStructMatrix(1, 1, 7, names);
    mxSetField(output, 0, "positions_m", positions);
    mxSetField(output, 0, "tangents", tangents);
    mxSetField(output, 0, "stations_m", stations);
    mxSetField(output, 0, "signed_curvature_per_m", curvature);
    mxSetField(output, 0, "length_m", mxCreateDoubleScalar(tracked.length_m));
    mxSetField(output, 0, "out_of_plane_m",
               mxCreateDoubleScalar(tracked.out_of_plane_m));
    mxSetField(output, 0, "out_of_plane_slope",
               mxCreateDoubleScalar(tracked.out_of_plane_slope));
    return output;
}

void TrackReferenceOrbit3DCommand(int nlhs, mxArray* plhs[], int nrhs,
                                  const mxArray* prhs[]) {
    if (nrhs != 2 || nlhs != 1)
        throw std::invalid_argument(
            "result = radia_mex('beam.orbit.track_reference_3d', config)");
    const mxArray* config = prhs[1];
    RequireScalarStruct(config);
    const double radia_object = Scalar(
        Field(config, "radia_object"), "radia_object", 0.0, false);
    if (radia_object < 1.0 || radia_object != std::floor(radia_object))
        throw std::invalid_argument(
            "radia_object must be a positive Radia object key");
    const bool mirror_z = OptionalBoolean(
        Field(config, "mirror_z"), "mirror_z", false);
    const double rigidity = Scalar(
        Field(config, "magnetic_rigidity_t_m"),
        "magnetic_rigidity_t_m", 0.0, false);
    const mxArray* point_value = Field(config, "entrance_point_m");
    const mxArray* direction_value = Field(config, "entrance_direction");
    const double* point = RealData(point_value, "entrance_point_m");
    const double* direction = RealData(
        direction_value, "entrance_direction");
    if (mxGetNumberOfElements(point_value) != 3
        || mxGetNumberOfElements(direction_value) != 3)
        throw std::invalid_argument(
            "entrance_point_m/entrance_direction must have three entries");
    const double exit_x_m = Scalar(
        Field(config, "exit_x_m"), "exit_x_m", 0.0, false);
    const double step_m = Scalar(
        Field(config, "step_m"), "step_m", 1.0e-3, true);
    const double maximum_path_m = Scalar(
        Field(config, "maximum_path_m"), "maximum_path_m", 0.14, true);
    const double planarity_tolerance_m = Scalar(
        Field(config, "planarity_tolerance_m"), "planarity_tolerance_m",
        1.0e-6, true);
    const std::size_t station_count = PositiveInteger(
        Field(config, "station_count"), "station_count", 65, true);
    if (station_count < 2)
        throw std::invalid_argument("station_count must be at least 2");
    // (3, N) MATLAB blocks are byte-identical to (N, 3) row-major buffers.
    mxArray* positions = mxCreateDoubleMatrix(3, station_count, mxREAL);
    mxArray* tangents = mxCreateDoubleMatrix(3, station_count, mxREAL);
    mxArray* stations = mxCreateDoubleMatrix(station_count, 1, mxREAL);
    mxArray* curvature = mxCreateDoubleMatrix(station_count - 1, 1, mxREAL);
    const rad_orbit::OrbitTrackResult tracked =
        rad_orbit::TrackReferenceOrbit3D(
            nullptr, 0.0, static_cast<int>(radia_object), mirror_z ? 1 : 0,
            rigidity, point, direction, exit_x_m, step_m, maximum_path_m,
            planarity_tolerance_m, station_count, mxGetPr(positions),
            mxGetPr(tangents), mxGetPr(stations), mxGetPr(curvature));
    plhs[0] = OrbitTrackResultStruct(
        positions, tangents, stations, curvature, tracked);
}

}  // namespace

bool DispatchBeamCommand(const std::string& command, int nlhs,
                         mxArray* plhs[], int nrhs,
                         const mxArray* prhs[]) {
    if (command == "beam.species.proton" ||
        command == "beam.species.electron") {
        ParticleSpeciesPreset(command, nlhs, plhs, nrhs);
        return true;
    }
    if (command == "beam.reference_particle.from_kinetic_energy_ev") {
        ReferenceParticle(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "beam.field.sample") {
        BeamFieldSample(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "beam.equation.rhs") {
        LorentzRHS(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "beam.step") {
        BeamStep(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "beam.track") {
        BeamTrack(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "beam.transfer.propagate_variational") {
        Propagate(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "beam.hamiltonian.canonical_body_jet") {
        CanonicalHamiltonianJet(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "beam.lie.map_tensors_spoly") {
        LieMapTensorsSpoly(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "beam.lie.dragt_finn_factorize") {
        LieDragtFinnFactorize(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "beam.lie.apply_dragt_finn_batch") {
        LieApplyDragtFinnBatch(nlhs, plhs, nrhs, prhs);
        return true;
    }
    if (command == "beam.orbit.track_reference_3d") {
        TrackReferenceOrbit3DCommand(nlhs, plhs, nrhs, prhs);
        return true;
    }
    return false;
}

void BeamTransferFromGridFunction(
        std::shared_ptr<ngcomp::GridFunction> field, int nlhs,
        mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    PropagateGridFunction(std::move(field), nlhs, plhs, nrhs, prhs);
}

void BeamTrackGridFunction(
        std::shared_ptr<ngcomp::GridFunction> field, int nlhs,
        mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 3 || nlhs != 1)
        throw std::invalid_argument(
            "trajectory = radia_mex('beam.track.grid_function', "
            "grid_function_handle, config)");
    const mxArray* config = prhs[2];
    RequireScalarStruct(config);
    const mxArray* schema = Field(config, "schema");
    if (schema && Text(schema, "schema") != "radia.beam.tracking.v1")
        throw std::invalid_argument("unsupported beam tracking schema");
    const auto species = ParseSpecies(Field(config, "species"));
    const auto state = ParseState(Field(config, "state"));
    const auto independent = ParseIndependent(Field(config, "independent"));
    auto equation = std::make_shared<radia::beam::LorentzEquation>(
        species,
        std::make_shared<radia::beam::NGSolveGridFunctionField>(
            std::move(field), ParseGridFunctionMagneticInput(
                Field(config, "field_representation"))),
        independent);
    auto stepper = ParseStepper(Field(config, "stepper"));
    radia::beam::TrackPlan plan;
    plan.start = Scalar(Field(config, "start"), "start", 0.0, false);
    plan.stop = Scalar(Field(config, "stop"), "stop", 0.0, false);
    plan.maximum_step = Scalar(
        Field(config, "maximum_step"), "maximum_step", 0.0, false);
    plan.maximum_steps = PositiveInteger(
        Field(config, "maximum_steps"), "maximum_steps", 1000000, true);
    plhs[0] = TrajectoryStruct(
        radia::beam::Tracker(equation, stepper).Track(state, plan));
    ReplaceStructField(plhs[0], "backend", mxCreateString(
        "native-cpp-ngsolve-gridfunction-mex"));
}

void BeamTrackReferenceOrbitToPlane(
        std::shared_ptr<rad_hdiv::HDivFieldEvaluator> field, int nlhs,
        mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (!field || nrhs != 3 || nlhs != 1)
        throw std::invalid_argument(
            "result = radia_mex('beam.orbit.track_reference_to_plane', "
            "field_evaluator_handle, config)");
    const mxArray* config = prhs[2];
    RequireScalarStruct(config);
    const double iron_scale = Scalar(
        Field(config, "iron_scale"), "iron_scale", 1.0, true);
    const mxArray* algorithm_value = Field(config, "iron_algorithm");
    const std::string iron_algorithm = algorithm_value
        ? Text(algorithm_value, "iron_algorithm") : "auto";
    const int iron_algorithm_code = iron_algorithm == "auto" ? -1
        : (iron_algorithm == "direct" ? 0
        : (iron_algorithm == "tree" ? 1 : -2));
    if (iron_algorithm_code < -1)
        throw std::invalid_argument(
            "iron_algorithm must be 'auto', 'direct', or 'tree'");
    const double radia_object = Scalar(
        Field(config, "radia_object"), "radia_object", -1.0, true);
    if (!((radia_object == -1.0) ||
          (radia_object >= 1.0 && radia_object == std::floor(radia_object))))
        throw std::invalid_argument(
            "radia_object must be -1 (absent) or a positive object key");
    const bool mirror_z = OptionalBoolean(
        Field(config, "mirror_z"), "mirror_z", false);
    const double rigidity = Scalar(
        Field(config, "magnetic_rigidity_t_m"),
        "magnetic_rigidity_t_m", 0.0, false);
    const mxArray* constant_value = Field(config, "constant_field_t");
    const mxArray* point_value = Field(config, "entrance_point_m");
    const mxArray* direction_value = Field(config, "entrance_direction");
    const mxArray* normal_value = Field(config, "exit_plane_normal");
    const double* constant_field = RealData(
        constant_value, "constant_field_t");
    const double* point = RealData(point_value, "entrance_point_m");
    const double* direction = RealData(
        direction_value, "entrance_direction");
    const double* normal = RealData(normal_value, "exit_plane_normal");
    if (mxGetNumberOfElements(constant_value) != 3
        || mxGetNumberOfElements(point_value) != 3
        || mxGetNumberOfElements(direction_value) != 3
        || mxGetNumberOfElements(normal_value) != 3)
        throw std::invalid_argument(
            "constant field, entrance point/direction, and exit-plane "
            "normal must have three entries");
    const double exit_plane_offset_m = Scalar(
        Field(config, "exit_plane_offset_m"),
        "exit_plane_offset_m", 0.0, false);
    const double step_m = Scalar(
        Field(config, "step_m"), "step_m", 1.0e-3, true);
    const double maximum_path_m = Scalar(
        Field(config, "maximum_path_m"), "maximum_path_m", 0.14, true);
    const double planarity_tolerance_m = Scalar(
        Field(config, "planarity_tolerance_m"), "planarity_tolerance_m",
        1.0e-6, true);
    const std::size_t station_count = PositiveInteger(
        Field(config, "station_count"), "station_count", 65, true);
    if (station_count < 2)
        throw std::invalid_argument("station_count must be at least 2");
    // (3, N) MATLAB blocks are byte-identical to (N, 3) row-major buffers.
    mxArray* positions = mxCreateDoubleMatrix(3, station_count, mxREAL);
    mxArray* tangents = mxCreateDoubleMatrix(3, station_count, mxREAL);
    mxArray* stations = mxCreateDoubleMatrix(station_count, 1, mxREAL);
    mxArray* curvature = mxCreateDoubleMatrix(station_count - 1, 1, mxREAL);
    const rad_orbit::OrbitTrackResult tracked =
        rad_orbit::TrackReferenceOrbit3DToPlane(
            field.get(), iron_scale, iron_algorithm_code,
            static_cast<int>(radia_object), mirror_z ? 1 : 0,
            constant_field, rigidity, point, direction, normal,
            exit_plane_offset_m, step_m, maximum_path_m,
            planarity_tolerance_m, station_count, mxGetPr(positions),
            mxGetPr(tangents), mxGetPr(stations), mxGetPr(curvature));
    plhs[0] = OrbitTrackResultStruct(
        positions, tangents, stations, curvature, tracked);
}
