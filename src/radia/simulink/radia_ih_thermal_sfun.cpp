#define S_FUNCTION_NAME radia_ih_thermal_sfun
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "mex.h"

#include "radia_ih_thermal.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <climits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Context {
    radia::ih::CSRMatrix mass;
    radia::ih::CSRMatrix stiffness;
    radia::ih::CSRMatrix convection;
    bool has_convection = false;
    radia::ih::ThermalState state;
    std::vector<double> weights;
    radia::ih::ThermalStepOptions options;
};

const mxArray* field(const mxArray* a, const char* name) {
    return a && mxIsStruct(a) ? mxGetField(a, 0, name) : nullptr;
}

std::vector<double> doubles(const mxArray* a, const char* name) {
    if (!a || !mxIsDouble(a) || mxIsComplex(a))
        throw std::invalid_argument(std::string("IH thermal config field '") + name + "' must be real double");
    const double* p = mxGetPr(a);
    return std::vector<double>(p, p + mxGetNumberOfElements(a));
}

std::vector<int> ints(const mxArray* a, const char* name) {
    if (!a || !mxIsDouble(a) || mxIsComplex(a))
        throw std::invalid_argument(std::string("IH thermal config field '") + name + "' must be numeric");
    std::vector<int> result;
    result.reserve(mxGetNumberOfElements(a));
    for (mwIndex i = 0; i < mxGetNumberOfElements(a); ++i) {
        result.push_back(static_cast<int>(mxGetDoubles(a)[i]));
    }
    return result;
}

double scalar(const mxArray* a, const char* name, double fallback) {
    if (!a) return fallback;
    if (!mxIsNumeric(a) || mxIsComplex(a) || mxGetNumberOfElements(a) != 1)
        throw std::invalid_argument(std::string("IH thermal config field '") + name + "' must be scalar");
    return mxGetScalar(a);
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
    const auto n_value = field(config, "n_temperature");
    const int n = static_cast<int>(scalar(n_value, "n_temperature", 0.0));
    if (n <= 0) throw std::invalid_argument("n_temperature must be positive");
    auto* context = new Context();
    context->mass = matrix(config, "mass", n);
    context->stiffness = matrix(config, "stiffness", n);
    if (field(config, "convection_row_ptr")) {
        context->convection = matrix(config, "convection", n);
        context->has_convection = true;
    }
    context->weights = doubles(field(config, "temperature_cell_weights"), "temperature_cell_weights");
    context->state.temperature_K = doubles(field(config, "initial_temperature_K"), "initial_temperature_K");
    if (context->weights.size() != static_cast<std::size_t>(n) ||
        context->state.temperature_K.size() != static_cast<std::size_t>(n))
        throw std::invalid_argument("IH thermal vectors do not match n_temperature");
    context->options.dt_s = scalar(field(config, "sample_time_s"), "sample_time_s", 0.0);
    context->options.tolerance = scalar(field(config, "thermal_tolerance"), "thermal_tolerance", 1.0e-10);
    context->options.convection_W_per_m2K = scalar(field(config, "convection_W_per_m2K"), "convection_W_per_m2K", 0.0);
    context->options.max_iterations = static_cast<int>(scalar(field(config, "thermal_max_iterations"), "thermal_max_iterations", 500));
    if (!(context->options.dt_s > 0.0)) throw std::invalid_argument("sample_time_s must be positive");
    return context;
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
    if (!std::isfinite(n_value) || n_value <= 0.0 || n_value > static_cast<double>(INT_MAX)) {
        ssSetErrorStatus(S, "radia_ih_thermal_sfun: n_temperature must be positive."); return;
    }
    const int n = static_cast<int>(n_value);
    if (!ssSetNumInputPorts(S, 3)) return;
    ssSetInputPortWidth(S, 0, n);
    ssSetInputPortWidth(S, 1, 1);
    ssSetInputPortWidth(S, 2, 1);
    for (int port = 0; port < 3; ++port) {
        ssSetInputPortDirectFeedThrough(S, port, 1);
        ssSetInputPortRequiredContiguous(S, port, 1);
        ssSetInputPortDataType(S, port, SS_DOUBLE);
    }
    if (!ssSetNumOutputPorts(S, 1)) return;
    ssSetOutputPortWidth(S, 0, n);
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);
    ssSetNumSampleTimes(S, 1);
    ssSetNumPWork(S, 1);
    ssSetNumDWork(S, 0);
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}

static void mdlInitializeSampleTimes(SimStruct* S) {
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_START
static void mdlStart(SimStruct* S) {
    try {
        auto* context = make_context(ssGetSFcnParam(S, 0));
        ssSetInputPortWidth(S, 0, static_cast<int>(context->state.temperature_K.size()));
        ssSetOutputPortWidth(S, 0, static_cast<int>(context->state.temperature_K.size()));
        ssGetPWork(S)[0] = context;
    } catch (const std::exception& error) {
        ssSetErrorStatus(S, error.what());
    }
}

static void mdlOutputs(SimStruct* S, int_T) {
    auto* context = static_cast<Context*>(ssGetPWork(S)[0]);
    if (!context) return;
    try {
        const auto* heat = static_cast<const real_T*>(ssGetInputPortSignal(S, 0));
        const auto* ambient = static_cast<const real_T*>(ssGetInputPortSignal(S, 1));
        const auto* angle = static_cast<const real_T*>(ssGetInputPortSignal(S, 2));
        context->options.ambient_temperature_K = ambient[0];
        radia::ih::advance_thermal(context->mass, context->stiffness,
                                   context->has_convection ? &context->convection : nullptr,
                                   std::vector<double>(heat, heat + context->state.temperature_K.size()),
                                   context->weights, angle[0], context->options, context->state);
        auto* output = static_cast<real_T*>(ssGetOutputPortSignal(S, 0));
        std::copy(context->state.temperature_K.begin(), context->state.temperature_K.end(), output);
    } catch (const std::exception& error) {
        ssSetErrorStatus(S, error.what());
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
