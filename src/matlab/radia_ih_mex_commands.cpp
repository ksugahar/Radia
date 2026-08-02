#include "radia_ih_mex_commands.h"

#include "radia_ih_runtime.h"

#include <algorithm>
#include <atomic>
#include <climits>
#include <cmath>
#include <complex>
#include <cstdint>
#include <initializer_list>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

template <class T>
using Registry = std::unordered_map<std::uint64_t, std::shared_ptr<T>>;

std::mutex registry_mutex;
Registry<radia::ih::EddyRuntime> eddy_registry;
Registry<radia::ih::ThermalRuntime> thermal_registry;
std::atomic<std::uint64_t> next_handle{UINT64_C(0x8000000000000000)};

const mxArray* field(const mxArray* value, const char* name) {
    return value && mxIsStruct(value) ? mxGetField(value, 0, name) : nullptr;
}

void require_scalar_struct(const mxArray* value, const char* context) {
    if (!value || !mxIsStruct(value) || mxGetNumberOfElements(value) != 1)
        throw std::invalid_argument(std::string(context) +
                                    " requires a scalar config struct");
}

std::vector<double> numbers(const mxArray* value, const char* name,
                            bool optional = false) {
    if (!value && optional) return {};
    if (!value || !mxIsDouble(value) || mxIsComplex(value))
        throw std::invalid_argument(std::string(name) +
                                    " must be a real double array");
    const std::size_t count = mxGetNumberOfElements(value);
    if (count == 0) return {};
    const double* data = mxGetPr(value);
    std::vector<double> result(data, data + count);
    if (!std::all_of(result.begin(), result.end(),
                     [](double item) { return std::isfinite(item); }))
        throw std::invalid_argument(std::string(name) +
                                    " must contain finite values");
    return result;
}

std::vector<double> row_major_numbers(const mxArray* value, const char* name,
                                      std::size_t expected_count,
                                      bool optional = false) {
    std::vector<double> result = numbers(value, name, optional);
    if (optional && !value) return result;
    if (result.size() != expected_count)
        throw std::invalid_argument(std::string(name) +
                                    " has the wrong number of values");
    if (mxGetM(value) > 1 && mxGetN(value) > 1)
        throw std::invalid_argument(std::string(name) +
            " must be supplied as an explicit row-major vector");
    return result;
}

std::vector<int> indices(const mxArray* value, const char* name) {
    const std::vector<double> source = numbers(value, name);
    std::vector<int> result;
    result.reserve(source.size());
    for (double item : source) {
        if (item < 0.0 || item > static_cast<double>(INT_MAX) ||
            std::floor(item) != item)
            throw std::invalid_argument(std::string(name) +
                                        " must contain nonnegative integers");
        result.push_back(static_cast<int>(item));
    }
    return result;
}

double scalar(const mxArray* value, const char* name, double fallback = 0.0,
              bool optional = false) {
    if (!value && optional) return fallback;
    if (!value || !mxIsNumeric(value) || mxIsComplex(value) ||
        mxGetNumberOfElements(value) != 1)
        throw std::invalid_argument(std::string(name) + " must be a scalar");
    const double result = mxGetScalar(value);
    if (!std::isfinite(result))
        throw std::invalid_argument(std::string(name) + " must be finite");
    return result;
}

int positive_integer(const mxArray* value, const char* name, int fallback = 0,
                     bool optional = false) {
    const double result = scalar(value, name, fallback, optional);
    if (result <= 0.0 || result > static_cast<double>(INT_MAX) ||
        std::floor(result) != result)
        throw std::invalid_argument(std::string(name) +
                                    " must be a positive integer");
    return static_cast<int>(result);
}

std::string text_value(const mxArray* value, const char* name,
                       const char* fallback = nullptr) {
    if (!value) {
        if (fallback) return fallback;
        throw std::invalid_argument(std::string(name) + " is required");
    }
    if (!mxIsChar(value))
        throw std::invalid_argument(std::string(name) + " must be text");
    char* buffer = mxArrayToUTF8String(value);
    if (!buffer)
        throw std::invalid_argument(std::string(name) + " is invalid");
    std::string result(buffer);
    mxFree(buffer);
    return result;
}

void require_choice(const std::string& value,
                    std::initializer_list<const char*> choices,
                    const char* name) {
    for (const char* choice : choices)
        if (value == choice) return;
    throw std::invalid_argument(std::string(name) + " has an unsupported value");
}

void validate_contract(const mxArray* config) {
    require_scalar_struct(config, "IH native runtime");
    if (text_value(field(config, "schema"), "schema") !=
        "radia.ih.simulink.native_sfunction.v1")
        throw std::invalid_argument("unsupported IH native configuration schema");
    if (const mxArray* backend = field(config, "backend")) {
        if (text_value(backend, "backend") !=
            "matlab-level2+radia-mex-handles")
            throw std::invalid_argument("unsupported IH native backend");
    }
    if (const mxArray* fallback = field(config, "python_fallback")) {
        if ((!mxIsLogical(fallback) && !mxIsNumeric(fallback)) ||
            mxIsComplex(fallback) || mxGetNumberOfElements(fallback) != 1 ||
            mxGetScalar(fallback) != 0.0)
            throw std::invalid_argument(
                "IH native runtime does not permit Python fallback");
    }
}

bool periodic_rotation(const mxArray* config) {
    const std::string mode =
        text_value(field(config, "rotation_mode"), "rotation_mode", "none");
    require_choice(mode, {"none", "periodic-uniform"}, "rotation_mode");
    return mode == "periodic-uniform";
}

std::vector<std::complex<double>> complex_values(
        const mxArray* config, const char* real_name, const char* imag_name,
        std::size_t expected_count, bool optional = false) {
    const mxArray* real_part = field(config, real_name);
    const mxArray* imag_part = field(config, imag_name);
    if (optional && !real_part && !imag_part) return {};
    if (!real_part || !imag_part)
        throw std::invalid_argument(std::string(real_name) + " and " +
                                    imag_name + " must be supplied together");
    const std::vector<double> real =
        row_major_numbers(real_part, real_name, expected_count);
    const std::vector<double> imag =
        row_major_numbers(imag_part, imag_name, expected_count);
    std::vector<std::complex<double>> result(expected_count);
    for (std::size_t index = 0; index < expected_count; ++index)
        result[index] = {real[index], imag[index]};
    return result;
}

radia::ih::CSRMatrix csr(const mxArray* config, const char* prefix, int n) {
    const std::string base(prefix);
    radia::ih::CSRMatrix result;
    result.n = n;
    result.row_ptr = indices(
        field(config, (base + "_row_ptr").c_str()),
        (base + "_row_ptr").c_str());
    result.col = indices(field(config, (base + "_col").c_str()),
                         (base + "_col").c_str());
    result.value = numbers(field(config, (base + "_value").c_str()),
                           (base + "_value").c_str());
    return result;
}

radia::ih::EddyConfig eddy_config(const mxArray* config) {
    validate_contract(config);
    if (text_value(field(config, "bh_mode"), "bh_mode") != "linear")
        throw std::invalid_argument(
            "IH native runtime currently supports only bh_mode='linear'");
    require_choice(text_value(field(config, "eddy_solver"), "eddy_solver"),
                   {"fem", "peec", "bem-a", "bim"}, "eddy_solver");

    radia::ih::EddyConfig result;
    result.n_unknown =
        positive_integer(field(config, "n_eddy_unknown"), "n_eddy_unknown");
    result.n_heat = positive_integer(field(config, "n_heat"), "n_heat");
    result.n_temperature =
        positive_integer(field(config, "n_temperature"), "n_temperature");
    const std::size_t matrix_size =
        static_cast<std::size_t>(result.n_unknown) * result.n_unknown;
    result.matrix = complex_values(config, "eddy_matrix_real",
                                   "eddy_matrix_imag", matrix_size);
    result.matrix_temperature_slope = complex_values(
        config, "eddy_matrix_temperature_slope_real",
        "eddy_matrix_temperature_slope_imag",
        static_cast<std::size_t>(result.n_temperature) * matrix_size, true);
    if (!result.matrix_temperature_slope.empty() &&
        !field(config, "bh_reference_temperature_K"))
        throw std::invalid_argument(
            "temperature-dependent Eddy data require bh_reference_temperature_K");
    result.reference_temperature_K = scalar(
        field(config, "bh_reference_temperature_K"),
        "bh_reference_temperature_K", 293.15, true);
    result.rhs_per_amp = complex_values(
        config, "eddy_rhs_real", "eddy_rhs_imag",
        static_cast<std::size_t>(result.n_unknown));
    result.heat_projection = row_major_numbers(
        field(config, "heat_projection"), "heat_projection",
        static_cast<std::size_t>(result.n_heat) * result.n_unknown);
    result.heat_weights = row_major_numbers(
        field(config, "heat_cell_weights"), "heat_cell_weights",
        static_cast<std::size_t>(result.n_heat));
    result.temperature_weights = row_major_numbers(
        field(config, "temperature_cell_weights"),
        "temperature_cell_weights",
        static_cast<std::size_t>(result.n_temperature));
    result.periodic_rotation = periodic_rotation(config);
    result.angle_origin_rad = scalar(field(config, "angle_origin_rad"),
                                     "angle_origin_rad", 0.0, true);
    return result;
}

radia::ih::ThermalConfig thermal_config(const mxArray* config) {
    validate_contract(config);
    require_choice(text_value(field(config, "thermal_solver"), "thermal_solver"),
                   {"fem"}, "thermal_solver");

    radia::ih::ThermalConfig result;
    const int n =
        positive_integer(field(config, "n_temperature"), "n_temperature");
    result.n_heat = positive_integer(field(config, "n_heat"), "n_heat");
    result.mass = csr(config, "mass", n);
    result.stiffness = csr(config, "stiffness", n);

    const bool convection_row = field(config, "convection_row_ptr") != nullptr;
    const bool convection_col = field(config, "convection_col") != nullptr;
    const bool convection_value = field(config, "convection_value") != nullptr;
    if ((convection_row || convection_col || convection_value) &&
        !(convection_row && convection_col && convection_value))
        throw std::invalid_argument(
            "convection CSR row, column, and value arrays are required together");
    result.has_convection = convection_row;
    if (result.has_convection) result.convection = csr(config, "convection", n);

    result.initial_temperature_K = row_major_numbers(
        field(config, "initial_temperature_K"), "initial_temperature_K",
        static_cast<std::size_t>(n));
    result.weights = row_major_numbers(
        field(config, "temperature_cell_weights"),
        "temperature_cell_weights", static_cast<std::size_t>(n));
    result.heat_to_temperature = row_major_numbers(
        field(config, "heat_to_temperature_projection"),
        "heat_to_temperature_projection",
        static_cast<std::size_t>(n) * result.n_heat);
    result.options.dt_s =
        scalar(field(config, "sample_time_s"), "sample_time_s");
    result.options.tolerance = scalar(field(config, "thermal_tolerance"),
                                      "thermal_tolerance", 1.0e-10, true);
    result.options.max_iterations = positive_integer(
        field(config, "thermal_max_iterations"),
        "thermal_max_iterations", 500, true);
    result.options.convection_W_per_m2K = scalar(
        field(config, "convection_W_per_m2K"),
        "convection_W_per_m2K", 0.0, true);
    result.periodic_rotation = periodic_rotation(config);
    result.angle_origin_rad = scalar(field(config, "angle_origin_rad"),
                                     "angle_origin_rad", 0.0, true);
    return result;
}

std::uint64_t input_handle(const mxArray* value) {
    if (!value || !mxIsUint64(value) || mxGetNumberOfElements(value) != 1)
        throw std::invalid_argument("handle must be a uint64 scalar");
    return *static_cast<const std::uint64_t*>(mxGetData(value));
}

mxArray* handle_output(std::uint64_t handle) {
    mxArray* result =
        mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    *static_cast<std::uint64_t*>(mxGetData(result)) = handle;
    return result;
}

mxArray* column(const std::vector<double>& values) {
    mxArray* result =
        mxCreateDoubleMatrix(static_cast<mwSize>(values.size()), 1, mxREAL);
    std::copy(values.begin(), values.end(), mxGetPr(result));
    return result;
}

void require_arity(int nrhs, int expected_rhs, int nlhs, int expected_lhs,
                   const char* usage) {
    if (nrhs != expected_rhs || nlhs != expected_lhs)
        throw std::invalid_argument(usage);
}

template <class T>
std::uint64_t insert(Registry<T>& registry, std::shared_ptr<T> value) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    std::uint64_t handle = next_handle.fetch_add(1);
    while (handle == 0 || eddy_registry.count(handle) != 0 ||
           thermal_registry.count(handle) != 0)
        handle = next_handle.fetch_add(1);
    registry.emplace(handle, std::move(value));
    mexLock();
    return handle;
}

template <class T>
std::shared_ptr<T> get(const Registry<T>& registry, std::uint64_t handle,
                       const char* message) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = registry.find(handle);
    if (found == registry.end()) throw std::invalid_argument(message);
    return found->second;
}

template <class T>
void erase(Registry<T>& registry, std::uint64_t handle, const char* message) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    if (registry.erase(handle) == 0) throw std::invalid_argument(message);
    mexUnlock();
}

}  // namespace

void CleanupIHHandles() {
    std::size_t count = 0;
    {
        std::lock_guard<std::mutex> guard(registry_mutex);
        count = eddy_registry.size() + thermal_registry.size();
        eddy_registry.clear();
        thermal_registry.clear();
    }
    for (std::size_t index = 0; index < count && mexIsLocked(); ++index)
        mexUnlock();
}

std::size_t IHHandleCount() {
    std::lock_guard<std::mutex> guard(registry_mutex);
    return eddy_registry.size() + thermal_registry.size();
}

bool DispatchIHCommand(const std::string& command, int nlhs, mxArray* plhs[],
                       int nrhs, const mxArray* prhs[]) {
    if (command == "ih.eddy.create") {
        require_arity(nrhs, 2, nlhs, 1,
                      "h = radia_mex('ih.eddy.create', config)");
        auto runtime = std::make_shared<radia::ih::EddyRuntime>(
            eddy_config(prhs[1]));
        plhs[0] = handle_output(insert(eddy_registry, std::move(runtime)));
        return true;
    }
    if (command == "ih.eddy.output") {
        require_arity(
            nrhs, 5, nlhs, 1,
            "q = radia_mex('ih.eddy.output', h, current_A, angle_rad, temperature_K)");
        auto runtime = get(eddy_registry, input_handle(prhs[1]),
                           "invalid, stale, or wrong-type IH Eddy handle");
        plhs[0] = column(runtime->output(
            scalar(prhs[2], "current_A"), scalar(prhs[3], "angle_rad"),
            numbers(prhs[4], "temperature_K")));
        return true;
    }
    if (command == "ih.eddy.destroy") {
        require_arity(nrhs, 2, nlhs, 0,
                      "radia_mex('ih.eddy.destroy', h)");
        erase(eddy_registry, input_handle(prhs[1]),
              "invalid, stale, or wrong-type IH Eddy handle");
        return true;
    }
    if (command == "ih.thermal.create") {
        require_arity(nrhs, 2, nlhs, 1,
                      "h = radia_mex('ih.thermal.create', config)");
        auto runtime = std::make_shared<radia::ih::ThermalRuntime>(
            thermal_config(prhs[1]));
        plhs[0] = handle_output(insert(thermal_registry, std::move(runtime)));
        return true;
    }
    if (command == "ih.thermal.output") {
        require_arity(nrhs, 2, nlhs, 1,
                      "T = radia_mex('ih.thermal.output', h)");
        auto runtime = get(thermal_registry, input_handle(prhs[1]),
                           "invalid, stale, or wrong-type IH Thermal handle");
        plhs[0] = column(runtime->output());
        return true;
    }
    if (command == "ih.thermal.update") {
        require_arity(
            nrhs, 5, nlhs, 0,
            "radia_mex('ih.thermal.update', h, heat, ambient_K, angle_rad)");
        auto runtime = get(thermal_registry, input_handle(prhs[1]),
                           "invalid, stale, or wrong-type IH Thermal handle");
        runtime->update(numbers(prhs[2], "heat_density_W_per_m3"),
                        scalar(prhs[3], "ambient_temperature_K"),
                        scalar(prhs[4], "angle_rad"));
        return true;
    }
    if (command == "ih.thermal.reset") {
        require_arity(nrhs, 2, nlhs, 0,
                      "radia_mex('ih.thermal.reset', h)");
        get(thermal_registry, input_handle(prhs[1]),
            "invalid, stale, or wrong-type IH Thermal handle")->reset();
        return true;
    }
    if (command == "ih.thermal.destroy") {
        require_arity(nrhs, 2, nlhs, 0,
                      "radia_mex('ih.thermal.destroy', h)");
        erase(thermal_registry, input_handle(prhs[1]),
              "invalid, stale, or wrong-type IH Thermal handle");
        return true;
    }
    return false;
}
