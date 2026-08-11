#include "radia_reactor_mex_commands.h"

#include "radia_nonlinear_reactor_runtime.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

std::mutex registry_mutex;
std::unordered_map<std::uint64_t,
                   std::shared_ptr<radia::reactor::Runtime>> registry;
std::atomic<std::uint64_t> next_handle{UINT64_C(0x9000000000000000)};

const mxArray* field(const mxArray* value, const char* name) {
    return value && mxIsStruct(value) ? mxGetField(value, 0, name) : nullptr;
}

std::string text_value(const mxArray* value, const char* name) {
    if (!value || !mxIsChar(value))
        throw std::invalid_argument(std::string(name) + " must be text");
    char* buffer = mxArrayToUTF8String(value);
    if (!buffer)
        throw std::invalid_argument(std::string(name) + " is invalid");
    std::string result(buffer);
    mxFree(buffer);
    return result;
}

std::vector<double> numbers(const mxArray* value, const char* name) {
    if (!value || !mxIsDouble(value) || mxIsComplex(value))
        throw std::invalid_argument(std::string(name) +
                                    " must be a real double array");
    const double* data = mxGetPr(value);
    std::vector<double> result(
        data, data + mxGetNumberOfElements(value));
    if (!std::all_of(result.begin(), result.end(),
                     [](double item) { return std::isfinite(item); }))
        throw std::invalid_argument(std::string(name) +
                                    " must contain finite values");
    return result;
}

std::vector<double> row_major(const mxArray* value, const char* name,
                              std::size_t count) {
    std::vector<double> result = numbers(value, name);
    if (result.size() != count || (mxGetM(value) > 1 && mxGetN(value) > 1))
        throw std::invalid_argument(std::string(name) +
            " must be an explicit row-major vector with the expected size");
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

int positive_integer(const mxArray* value, const char* name) {
    const double result = scalar(value, name);
    if (!(result >= 1.0) || result != std::floor(result) ||
        result > static_cast<double>(std::numeric_limits<int>::max()))
        throw std::invalid_argument(std::string(name) +
                                    " must be a positive integer");
    return static_cast<int>(result);
}

radia::reactor::Config parse_config(const mxArray* value) {
    if (!value || !mxIsStruct(value) || mxGetNumberOfElements(value) != 1)
        throw std::invalid_argument(
            "nonlinear reactor requires a scalar configuration struct");
    if (text_value(field(value, "schema"), "schema") !=
        "radia.simulink.nonlinear-hdiv-mmm-reactor.v1")
        throw std::invalid_argument(
            "unsupported nonlinear reactor configuration schema");
    if (text_value(field(value, "backend"), "backend") !=
        "matlab-level2+radia-mex-handle")
        throw std::invalid_argument("unsupported nonlinear reactor backend");
    if (const mxArray* fallback = field(value, "python_per_step")) {
        if (mxIsComplex(fallback) || mxGetNumberOfElements(fallback) != 1 ||
            mxGetScalar(fallback) != 0.0)
            throw std::invalid_argument(
                "nonlinear reactor does not permit Python per simulation step");
    }

    radia::reactor::Config config;
    config.n_modes = positive_integer(field(value, "n_modes"), "n_modes");
    config.n_samples =
        positive_integer(field(value, "n_samples"), "n_samples");
    const std::size_t n_modes = static_cast<std::size_t>(config.n_modes);
    const std::size_t n_samples = static_cast<std::size_t>(config.n_samples);
    config.demag = row_major(field(value, "demag_row_major"),
                             "demag_row_major", n_modes * n_modes);
    config.magnetization_modes = row_major(
        field(value, "magnetization_modes_row_major"),
        "magnetization_modes_row_major", n_modes * n_samples * 3);
    config.sample_weights = row_major(field(value, "sample_weights"),
                                      "sample_weights", n_samples);
    config.excitation_per_amp = row_major(
        field(value, "excitation_per_amp"), "excitation_per_amp", n_modes);
    config.magnetization_table_A_per_m = numbers(
        field(value, "magnetization_table_A_per_m"),
        "magnetization_table_A_per_m");
    config.field_table_A_per_m = numbers(
        field(value, "field_table_A_per_m"), "field_table_A_per_m");
    config.air_inductance_H = scalar(
        field(value, "air_inductance_H"), "air_inductance_H");
    config.winding_resistance_Ohm = scalar(
        field(value, "winding_resistance_Ohm"), "winding_resistance_Ohm");
    config.sample_time_s = scalar(
        field(value, "sample_time_s"), "sample_time_s");
    config.initial_current_A = scalar(
        field(value, "initial_current_A"), "initial_current_A", 0.0, true);
    config.residual_tolerance = scalar(
        field(value, "residual_tolerance"), "residual_tolerance", 1.0e-9, true);
    config.max_iterations = field(value, "max_iterations")
        ? positive_integer(field(value, "max_iterations"), "max_iterations")
        : 40;
    config.line_search_minimum = scalar(
        field(value, "line_search_minimum"), "line_search_minimum",
        0x1p-20, true);
    return config;
}

std::uint64_t input_handle(const mxArray* value) {
    if (!value || !mxIsUint64(value) || mxGetNumberOfElements(value) != 1)
        throw std::invalid_argument("reactor handle must be a uint64 scalar");
    return *static_cast<const std::uint64_t*>(mxGetData(value));
}

mxArray* handle_output(std::uint64_t handle) {
    mxArray* output = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    *static_cast<std::uint64_t*>(mxGetData(output)) = handle;
    return output;
}

mxArray* column(const std::vector<double>& values) {
    mxArray* output = mxCreateDoubleMatrix(
        static_cast<mwSize>(values.size()), 1, mxREAL);
    std::copy(values.begin(), values.end(), mxGetPr(output));
    return output;
}

mxArray* output_struct(const radia::reactor::Output& value) {
    const char* fields[] = {
        "voltage_V", "flux_linkage_Wb_turn", "differential_inductance_H",
        "peak_flux_density_T", "magnetic_energy_J",
        "residual_relative_norm", "nonlinear_iterations", "flux_density_T",
        "magnetization_coefficients"};
    mxArray* output = mxCreateStructMatrix(1, 1, 9, fields);
    mxSetField(output, 0, fields[0], mxCreateDoubleScalar(value.voltage_V));
    mxSetField(output, 0, fields[1],
               mxCreateDoubleScalar(value.flux_linkage_Wb_turn));
    mxSetField(output, 0, fields[2],
               mxCreateDoubleScalar(value.differential_inductance_H));
    mxSetField(output, 0, fields[3],
               mxCreateDoubleScalar(value.peak_flux_density_T));
    mxSetField(output, 0, fields[4],
               mxCreateDoubleScalar(value.magnetic_energy_J));
    mxSetField(output, 0, fields[5],
               mxCreateDoubleScalar(value.residual_relative_norm));
    mxSetField(output, 0, fields[6],
               mxCreateDoubleScalar(value.nonlinear_iterations));
    mxSetField(output, 0, fields[7], column(value.flux_density_T));
    mxSetField(output, 0, fields[8],
               column(value.magnetization_coefficients));
    return output;
}

mxArray* snapshot_struct(const radia::reactor::Snapshot& value) {
    const char* fields[] = {"previous_current_A",
                            "previous_flux_linkage_Wb_turn",
                            "accepted_steps",
                            "magnetization_coefficients"};
    mxArray* output = mxCreateStructMatrix(1, 1, 4, fields);
    mxSetField(output, 0, fields[0],
               mxCreateDoubleScalar(value.previous_current_A));
    mxSetField(output, 0, fields[1],
               mxCreateDoubleScalar(value.previous_flux_linkage_Wb_turn));
    mxSetField(output, 0, fields[2],
               mxCreateDoubleScalar(static_cast<double>(value.accepted_steps)));
    mxSetField(output, 0, fields[3],
               column(value.magnetization_coefficients));
    return output;
}

radia::reactor::Snapshot parse_snapshot(const mxArray* value, int modes) {
    if (!value || !mxIsStruct(value) || mxGetNumberOfElements(value) != 1)
        throw std::invalid_argument("reactor snapshot must be a scalar struct");
    radia::reactor::Snapshot snapshot;
    snapshot.previous_current_A = scalar(
        field(value, "previous_current_A"), "previous_current_A");
    snapshot.previous_flux_linkage_Wb_turn = scalar(
        field(value, "previous_flux_linkage_Wb_turn"),
        "previous_flux_linkage_Wb_turn");
    const double steps = scalar(field(value, "accepted_steps"), "accepted_steps");
    if (steps < 0.0 || steps != std::floor(steps))
        throw std::invalid_argument("accepted_steps must be a nonnegative integer");
    snapshot.accepted_steps = static_cast<std::uint64_t>(steps);
    snapshot.magnetization_coefficients = row_major(
        field(value, "magnetization_coefficients"),
        "magnetization_coefficients", static_cast<std::size_t>(modes));
    return snapshot;
}

void require_arity(int nrhs, int expected_rhs, int nlhs, int expected_lhs,
                   const char* usage) {
    if (nrhs != expected_rhs || nlhs != expected_lhs)
        throw std::invalid_argument(usage);
}

std::shared_ptr<radia::reactor::Runtime> get(std::uint64_t handle) {
    std::lock_guard<std::mutex> guard(registry_mutex);
    const auto found = registry.find(handle);
    if (found == registry.end())
        throw std::invalid_argument(
            "invalid, stale, or wrong-type nonlinear reactor handle");
    return found->second;
}

}  // namespace

void CleanupReactorHandles() {
    std::size_t count = 0;
    {
        std::lock_guard<std::mutex> guard(registry_mutex);
        count = registry.size();
        registry.clear();
    }
    for (std::size_t index = 0; index < count && mexIsLocked(); ++index)
        mexUnlock();
}

std::size_t ReactorHandleCount() {
    std::lock_guard<std::mutex> guard(registry_mutex);
    return registry.size();
}

bool DispatchReactorCommand(const std::string& command, int nlhs,
                            mxArray* plhs[], int nrhs,
                            const mxArray* prhs[]) {
    if (command == "reactor.create") {
        require_arity(nrhs, 2, nlhs, 1,
                      "h = radia_mex('reactor.create', config)");
        auto runtime = std::make_shared<radia::reactor::Runtime>(
            parse_config(prhs[1]));
        std::uint64_t handle = next_handle.fetch_add(1);
        {
            std::lock_guard<std::mutex> guard(registry_mutex);
            while (handle == 0 || registry.count(handle) != 0)
                handle = next_handle.fetch_add(1);
            registry.emplace(handle, std::move(runtime));
        }
        mexLock();
        plhs[0] = handle_output(handle);
        return true;
    }
    if (command == "reactor.output") {
        require_arity(nrhs, 3, nlhs, 1,
                      "y = radia_mex('reactor.output', h, current_A)");
        plhs[0] = output_struct(get(input_handle(prhs[1]))->output(
            scalar(prhs[2], "current_A")));
        return true;
    }
    if (command == "reactor.update") {
        require_arity(nrhs, 3, nlhs, 0,
                      "radia_mex('reactor.update', h, current_A)");
        get(input_handle(prhs[1]))->update(scalar(prhs[2], "current_A"));
        return true;
    }
    if (command == "reactor.snapshot") {
        require_arity(nrhs, 2, nlhs, 1,
                      "s = radia_mex('reactor.snapshot', h)");
        plhs[0] = snapshot_struct(get(input_handle(prhs[1]))->snapshot());
        return true;
    }
    if (command == "reactor.restore") {
        require_arity(nrhs, 3, nlhs, 0,
                      "radia_mex('reactor.restore', h, snapshot)");
        auto runtime = get(input_handle(prhs[1]));
        runtime->restore(parse_snapshot(prhs[2], runtime->mode_count()));
        return true;
    }
    if (command == "reactor.reset") {
        require_arity(nrhs, 2, nlhs, 0,
                      "radia_mex('reactor.reset', h)");
        get(input_handle(prhs[1]))->reset();
        return true;
    }
    if (command == "reactor.info") {
        require_arity(nrhs, 2, nlhs, 1,
                      "info = radia_mex('reactor.info', h)");
        auto runtime = get(input_handle(prhs[1]));
        const char* fields[] = {"n_modes", "n_samples", "accepted_steps"};
        plhs[0] = mxCreateStructMatrix(1, 1, 3, fields);
        mxSetField(plhs[0], 0, fields[0],
                   mxCreateDoubleScalar(runtime->mode_count()));
        mxSetField(plhs[0], 0, fields[1],
                   mxCreateDoubleScalar(runtime->sample_count()));
        mxSetField(plhs[0], 0, fields[2], mxCreateDoubleScalar(
            static_cast<double>(runtime->accepted_steps())));
        return true;
    }
    if (command == "reactor.destroy") {
        require_arity(nrhs, 2, nlhs, 0,
                      "radia_mex('reactor.destroy', h)");
        const std::uint64_t handle = input_handle(prhs[1]);
        {
            std::lock_guard<std::mutex> guard(registry_mutex);
            if (registry.erase(handle) == 0)
                throw std::invalid_argument(
                    "invalid, stale, or wrong-type nonlinear reactor handle");
        }
        mexUnlock();
        return true;
    }
    return false;
}
