#define S_FUNCTION_NAME radia_ih_thermal_sfun
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "mex.h"

#include "radia_ih_thermal.h"
#include "radia_ih_transport.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <climits>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

void set_error_status(SimStruct* S, const std::exception& error) {
    static thread_local char message[1024];
    std::snprintf(message, sizeof(message), "radia_ih_thermal_sfun: %s",
                  error.what());
    ssSetErrorStatus(S, message);
}

struct Context {
    int n_heat = 0;
    radia::ih::CSRMatrix mass;
    radia::ih::CSRMatrix stiffness;
    radia::ih::CSRMatrix convection;
    bool has_convection = false;
    radia::ih::ThermalState state;
    std::vector<double> initial_temperature_K;
    std::vector<double> weights;
    // Row-major [temperature][heat-region] source projection.
    std::vector<double> heat_to_temperature;
    radia::ih::ThermalStepOptions options;
    bool periodic_rotation = false;
    double angle_origin_rad = 0.0;
};

const mxArray* field(const mxArray* a, const char* name) {
    return a && mxIsStruct(a) ? mxGetField(a, 0, name) : nullptr;
}

std::vector<double> doubles(const mxArray* a, const char* name) {
    if (!a || !mxIsDouble(a) || mxIsComplex(a))
        throw std::invalid_argument(std::string("IH thermal config field '") + name + "' must be real double");
    const double* p = mxGetPr(a);
    std::vector<double> result(p, p + mxGetNumberOfElements(a));
    if (!std::all_of(result.begin(), result.end(),
                     [](double value) { return std::isfinite(value); }))
        throw std::invalid_argument(std::string("IH thermal config field '") +
                                    name + "' must contain finite values");
    return result;
}

std::vector<int> ints(const mxArray* a, const char* name) {
    if (!a || !mxIsDouble(a) || mxIsComplex(a))
        throw std::invalid_argument(std::string("IH thermal config field '") + name + "' must be numeric");
    std::vector<int> result;
    result.reserve(mxGetNumberOfElements(a));
    for (mwIndex i = 0; i < mxGetNumberOfElements(a); ++i) {
        const double value = mxGetDoubles(a)[i];
        if (!std::isfinite(value) || value < 0.0 ||
            value > static_cast<double>(INT_MAX) ||
            std::floor(value) != value)
            throw std::invalid_argument(std::string("IH thermal config field '") +
                                        name + "' must contain nonnegative integers");
        result.push_back(static_cast<int>(value));
    }
    return result;
}

double scalar(const mxArray* a, const char* name, double fallback) {
    if (!a) return fallback;
    if (!mxIsNumeric(a) || mxIsComplex(a) || mxGetNumberOfElements(a) != 1)
        throw std::invalid_argument(std::string("IH thermal config field '") + name + "' must be scalar");
    const double value = mxGetScalar(a);
    if (!std::isfinite(value))
        throw std::invalid_argument(std::string("IH thermal config field '") +
                                    name + "' must be finite");
    return value;
}

int positive_integer(const mxArray* a, const char* name, int fallback = 0) {
    const double value = scalar(a, name, static_cast<double>(fallback));
    if (value <= 0.0 || value > static_cast<double>(INT_MAX) ||
        std::floor(value) != value)
        throw std::invalid_argument(std::string("IH thermal config field '") +
                                    name + "' must be a positive integer");
    return static_cast<int>(value);
}

std::string string_value(const mxArray* a, const char* name,
                         const char* fallback) {
    if (!a) return fallback;
    if (!mxIsChar(a))
        throw std::invalid_argument(std::string("IH thermal config field '") +
                                    name + "' must be text");
    char* text = mxArrayToUTF8String(a);
    if (!text)
        throw std::invalid_argument(std::string("IH thermal config field '") +
                                    name + "' is invalid");
    const std::string value(text);
    mxFree(text);
    return value;
}

bool equivalent_periodic_angle(double angle, double origin) {
    const double period = 2.0 * std::acos(-1.0);
    return std::abs(std::remainder(angle - origin, period)) <= 1.0e-12;
}

radia::ih::CSRMatrix matrix(const mxArray* config, const char* prefix, int n) {
    const std::string base(prefix);
    const auto row = ints(field(config, (base + "_row_ptr").c_str()), (base + "_row_ptr").c_str());
    const auto col = ints(field(config, (base + "_col").c_str()), (base + "_col").c_str());
    const auto value = doubles(field(config, (base + "_value").c_str()), (base + "_value").c_str());
    radia::ih::CSRMatrix result;
    result.n = n; result.row_ptr = row; result.col = col; result.value = value;
    return result;
}

Context* make_context(const mxArray* config) {
    if (!config || !mxIsStruct(config))
        throw std::invalid_argument("IH thermal S-Function requires a configuration struct");
    const int n = positive_integer(
        field(config, "n_temperature"), "n_temperature");
    auto context = std::make_unique<Context>();
    context->n_heat = positive_integer(field(config, "n_heat"), "n_heat");
    context->mass = matrix(config, "mass", n);
    context->stiffness = matrix(config, "stiffness", n);
    if (field(config, "convection_row_ptr")) {
        context->convection = matrix(config, "convection", n);
        context->has_convection = true;
    }
    context->weights = doubles(field(config, "temperature_cell_weights"), "temperature_cell_weights");
    if (const mxArray* projection = field(config, "heat_to_temperature_projection")) {
        context->heat_to_temperature = doubles(projection, "heat_to_temperature_projection");
        if (context->heat_to_temperature.size() !=
            static_cast<std::size_t>(n * context->n_heat))
            throw std::invalid_argument("heat_to_temperature_projection has the wrong shape");
    } else if (context->n_heat == n) {
        context->heat_to_temperature.assign(static_cast<std::size_t>(n * n), 0.0);
        for (int i = 0; i < n; ++i)
            context->heat_to_temperature[static_cast<std::size_t>(i * n + i)] = 1.0;
    } else {
        throw std::invalid_argument(
            "heat_to_temperature_projection is required when n_heat differs from n_temperature");
    }
    context->state.temperature_K = doubles(field(config, "initial_temperature_K"), "initial_temperature_K");
    context->initial_temperature_K = context->state.temperature_K;
    if (context->weights.size() != static_cast<std::size_t>(n) ||
        context->state.temperature_K.size() != static_cast<std::size_t>(n))
        throw std::invalid_argument("IH thermal vectors do not match n_temperature");
    for (double weight : context->weights)
        if (!(weight > 0.0))
            throw std::invalid_argument(
                "temperature_cell_weights must be positive");
    for (double temperature : context->state.temperature_K)
        if (!(temperature > 0.0))
            throw std::invalid_argument(
                "initial_temperature_K must be positive");
    context->options.dt_s = scalar(field(config, "sample_time_s"), "sample_time_s", 0.0);
    context->options.tolerance = scalar(field(config, "thermal_tolerance"), "thermal_tolerance", 1.0e-10);
    context->options.convection_W_per_m2K = scalar(field(config, "convection_W_per_m2K"), "convection_W_per_m2K", 0.0);
    context->options.max_iterations = positive_integer(
        field(config, "thermal_max_iterations"), "thermal_max_iterations", 500);
    if (!(context->options.dt_s > 0.0))
        throw std::invalid_argument("sample_time_s must be positive");
    if (!(context->options.tolerance > 0.0))
        throw std::invalid_argument("thermal_tolerance must be positive");
    if (context->options.convection_W_per_m2K < 0.0)
        throw std::invalid_argument(
            "convection_W_per_m2K must be nonnegative");
    const std::string rotation_mode = string_value(
        field(config, "rotation_mode"), "rotation_mode", "none");
    if (rotation_mode == "periodic-uniform") context->periodic_rotation = true;
    else if (rotation_mode != "none")
        throw std::invalid_argument(
            "IH thermal rotation_mode must be 'none' or 'periodic-uniform'");
    if (const mxArray* origin = field(config, "angle_origin_rad"))
        context->angle_origin_rad = scalar(origin, "angle_origin_rad", 0.0);
    context->state.previous_angle_rad = context->angle_origin_rad;
    return context.release();
}

}  // namespace

static void mdlInitializeSizes(SimStruct* S) {
    ssSetNumSFcnParams(S, 1);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) return;
    const mxArray* config = ssGetSFcnParam(S, 0);
    if (!mxIsStruct(config)) { ssSetErrorStatus(S, "radia_ih_thermal_sfun: Parameters must be a struct."); return; }
    const mxArray* n_field = field(config, "n_temperature");
    if (!n_field || mxGetNumberOfElements(n_field) != 1 || !mxIsNumeric(n_field)) {
        ssSetErrorStatus(S, "radia_ih_thermal_sfun: missing scalar parameter n_temperature."); return;
    }
    const double n_value = mxGetScalar(n_field);
    if (!std::isfinite(n_value) || n_value <= 0.0 ||
        n_value > static_cast<double>(INT_MAX) ||
        std::floor(n_value) != n_value) {
        ssSetErrorStatus(S, "radia_ih_thermal_sfun: n_temperature must be positive."); return;
    }
    const int n = static_cast<int>(n_value);
    const mxArray* n_heat_field = field(config, "n_heat");
    if (!n_heat_field || mxGetNumberOfElements(n_heat_field) != 1 ||
        !mxIsNumeric(n_heat_field)) {
        ssSetErrorStatus(S, "radia_ih_thermal_sfun: missing scalar parameter n_heat."); return;
    }
    const double n_heat_value = mxGetScalar(n_heat_field);
    if (!std::isfinite(n_heat_value) || n_heat_value <= 0.0 ||
        n_heat_value > static_cast<double>(INT_MAX) ||
        std::floor(n_heat_value) != n_heat_value) {
        ssSetErrorStatus(S, "radia_ih_thermal_sfun: n_heat must be positive."); return;
    }
    const int n_heat = static_cast<int>(n_heat_value);
    if (!ssSetNumInputPorts(S, 3)) return;
    ssSetInputPortWidth(S, 0, n_heat);
    ssSetInputPortWidth(S, 1, 1);
    ssSetInputPortWidth(S, 2, 1);
    for (int port = 0; port < 3; ++port) {
        // Inputs are consumed exactly once from mdlUpdate. mdlOutputs only
        // publishes the accepted state, so an Eddy -> Thermal -> Eddy loop
        // has a well-defined one-step delay and no algebraic loop.
        ssSetInputPortDirectFeedThrough(S, port, 0);
        ssSetInputPortRequiredContiguous(S, port, 1);
        ssSetInputPortDataType(S, port, SS_DOUBLE);
    }
    if (!ssSetNumOutputPorts(S, 1)) return;
    ssSetOutputPortWidth(S, 0, n);
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);
    ssSetNumSampleTimes(S, 1);
    ssSetNumPWork(S, 1);
    if (!ssSetNumDWork(S, 1)) return;
    ssSetDWorkWidth(S, 0, n);
    ssSetDWorkDataType(S, 0, SS_DOUBLE);
    ssSetDWorkName(S, 0, "temperature_K");
    ssSetDWorkUsedAsDState(S, 0, 1);
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}

static void mdlInitializeSampleTimes(SimStruct* S) {
    const mxArray* config = ssGetSFcnParam(S, 0);
    const double dt = scalar(field(config, "sample_time_s"), "sample_time_s", 0.0);
    if (!(dt > 0.0) || !std::isfinite(dt)) {
        ssSetErrorStatus(S, "radia_ih_thermal_sfun: sample_time_s must be positive.");
        return;
    }
    ssSetSampleTime(S, 0, dt);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_START
static void mdlStart(SimStruct* S) {
    try {
        auto* context = make_context(ssGetSFcnParam(S, 0));
        ssSetInputPortWidth(S, 0, context->n_heat);
        ssSetOutputPortWidth(S, 0, static_cast<int>(context->state.temperature_K.size()));
        ssGetPWork(S)[0] = context;
        auto* state = static_cast<real_T*>(ssGetDWork(S, 0));
        std::copy(context->initial_temperature_K.begin(),
                  context->initial_temperature_K.end(), state);
    } catch (const std::exception& error) {
        set_error_status(S, error);
    }
}

#define MDL_INITIALIZE_CONDITIONS
static void mdlInitializeConditions(SimStruct* S) {
    auto* context = static_cast<Context*>(ssGetPWork(S)[0]);
    if (!context) return;
    context->state.temperature_K = context->initial_temperature_K;
    context->state.time_s = 0.0;
    context->state.previous_angle_rad = context->angle_origin_rad;
    auto* state = static_cast<real_T*>(ssGetDWork(S, 0));
    std::copy(context->initial_temperature_K.begin(),
              context->initial_temperature_K.end(), state);
}

static void mdlOutputs(SimStruct* S, int_T) {
    auto* context = static_cast<Context*>(ssGetPWork(S)[0]);
    if (!context) return;
    auto* output = static_cast<real_T*>(ssGetOutputPortSignal(S, 0));
    const auto* state = static_cast<const real_T*>(ssGetDWork(S, 0));
    std::copy(state, state + context->state.temperature_K.size(), output);
}

#define MDL_UPDATE
static void mdlUpdate(SimStruct* S, int_T) {
    auto* context = static_cast<Context*>(ssGetPWork(S)[0]);
    if (!context) return;
    try {
        const auto* heat = static_cast<const real_T*>(ssGetInputPortSignal(S, 0));
        const auto* ambient = static_cast<const real_T*>(ssGetInputPortSignal(S, 1));
        const auto* angle = static_cast<const real_T*>(ssGetInputPortSignal(S, 2));
        if (!std::isfinite(ambient[0]) || !std::isfinite(angle[0]))
            throw std::invalid_argument("IH thermal ambient temperature and angle must be finite");
        const auto* accepted_state = static_cast<const real_T*>(ssGetDWork(S, 0));
        context->state.temperature_K.assign(
            accepted_state,
            accepted_state + context->initial_temperature_K.size());
        std::vector<double> source(context->state.temperature_K.size(), 0.0);
        for (std::size_t i = 0; i < source.size(); ++i) {
            for (int j = 0; j < context->n_heat; ++j) {
                const double value = heat[j];
                if (!std::isfinite(value))
                    throw std::invalid_argument("IH thermal heat input must be finite");
                source[i] += context->heat_to_temperature[
                    i * static_cast<std::size_t>(context->n_heat) +
                    static_cast<std::size_t>(j)] * value;
            }
        }
        if (context->periodic_rotation) {
            std::vector<double> transported;
            radia::ih::transport_periodic(
                context->state.temperature_K, context->weights,
                angle[0] - context->state.previous_angle_rad, transported);
            context->state.temperature_K = std::move(transported);
        } else if (!equivalent_periodic_angle(
                       angle[0], context->state.previous_angle_rad)) {
            throw std::invalid_argument(
                "IH thermal received a changing angle but rotation_mode is 'none'");
        }
        context->options.ambient_temperature_K = ambient[0];
        radia::ih::advance_thermal(context->mass, context->stiffness,
                                   context->has_convection ? &context->convection : nullptr,
                                   source,
                                   context->weights, angle[0], context->options, context->state);
        auto* next_state = static_cast<real_T*>(ssGetDWork(S, 0));
        std::copy(context->state.temperature_K.begin(),
                  context->state.temperature_K.end(), next_state);
    } catch (const std::exception& error) {
        set_error_status(S, error);
    }
}

static void mdlTerminate(SimStruct* S) {
    delete static_cast<Context*>(ssGetPWork(S)[0]);
    ssGetPWork(S)[0] = nullptr;
}

#ifdef MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif
