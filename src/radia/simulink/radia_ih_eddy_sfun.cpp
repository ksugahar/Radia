#define S_FUNCTION_NAME radia_ih_eddy_sfun
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "mex.h"

#include "radia_ih_transport.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <climits>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using Complex = std::complex<double>;

void set_error_status(SimStruct* S, const std::exception& error) {
    static thread_local char message[1024];
    std::snprintf(message, sizeof(message), "radia_ih_eddy_sfun: %s",
                  error.what());
    ssSetErrorStatus(S, message);
}

struct Context {
    int n_unknown = 0;
    int n_heat = 0;
    std::vector<Complex> matrix;
    // Optional temperature Jacobian, stored as [temperature][row][column].
    std::vector<Complex> matrix_temperature_slope;
    double reference_temperature_K = 293.15;
    std::vector<Complex> rhs_per_amp;
    std::vector<double> heat_projection;
    std::vector<double> heat_weights;
    std::vector<double> temperature_weights;
    std::vector<double> cached_heat;
    std::vector<double> reference_heat;
    std::vector<double> previous_temperature;
    bool have_cache = false;
    bool periodic_rotation = false;
    double angle_origin_rad = 0.0;
};

const mxArray* field(const mxArray* a, const char* name) {
    return a && mxIsStruct(a) ? mxGetField(a, 0, name) : nullptr;
}

std::vector<double> real_array(const mxArray* a, const char* name) {
    if (!a || !mxIsDouble(a) || mxIsComplex(a))
        throw std::invalid_argument(std::string("IH Eddy field '") + name + "' must be real double");
    std::vector<double> result(
        mxGetDoubles(a), mxGetDoubles(a) + mxGetNumberOfElements(a));
    if (!std::all_of(result.begin(), result.end(),
                     [](double value) { return std::isfinite(value); }))
        throw std::invalid_argument(std::string("IH Eddy field '") + name +
                                    "' must contain finite values");
    return result;
}

double scalar(const mxArray* a, const char* name) {
    if (!a || !mxIsNumeric(a) || mxIsComplex(a) || mxGetNumberOfElements(a) != 1)
        throw std::invalid_argument(std::string("IH Eddy field '") + name + "' must be a scalar");
    const double value = mxGetScalar(a);
    if (!std::isfinite(value))
        throw std::invalid_argument(std::string("IH Eddy field '") + name +
                                    "' must be finite");
    return value;
}

int positive_integer(const mxArray* a, const char* name) {
    const double value = scalar(a, name);
    if (value <= 0.0 || value > static_cast<double>(INT_MAX) ||
        std::floor(value) != value)
        throw std::invalid_argument(std::string("IH Eddy field '") + name +
                                    "' must be a positive integer");
    return static_cast<int>(value);
}

void validate_bh_mode(const mxArray* a) {
    if (!a) return;
    if (mxIsChar(a)) {
        char* text = mxArrayToUTF8String(a);
        if (!text) throw std::invalid_argument("IH Eddy field 'bh_mode' is invalid");
        const std::string value(text);
        mxFree(text);
        if (value == "linear") return;
        if (value == "nonlinear")
            throw std::invalid_argument(
                "IH native preview does not yet implement nonlinear BH iteration");
        throw std::invalid_argument("IH Eddy field 'bh_mode' must be 'linear'");
    }
    throw std::invalid_argument("IH Eddy field 'bh_mode' must be 'linear'");
}

std::string string_value(const mxArray* a, const char* name,
                         const char* fallback) {
    if (!a) return fallback;
    if (!mxIsChar(a))
        throw std::invalid_argument(std::string("IH Eddy field '") + name +
                                    "' must be text");
    char* text = mxArrayToUTF8String(a);
    if (!text)
        throw std::invalid_argument(std::string("IH Eddy field '") + name +
                                    "' is invalid");
    const std::string value(text);
    mxFree(text);
    return value;
}

bool equivalent_periodic_angle(double angle, double origin) {
    const double period = 2.0 * std::acos(-1.0);
    return std::abs(std::remainder(angle - origin, period)) <= 1.0e-12;
}

std::vector<Complex> complex_matrix(const mxArray* re, const mxArray* im,
                                    int rows, int cols, const char* name) {
    const auto r = real_array(re, name);
    const auto j = real_array(im, name);
    if (r.size() != static_cast<std::size_t>(rows * cols) || j.size() != r.size())
        throw std::invalid_argument(std::string("IH Eddy field '") + name + "' has the wrong shape");
    std::vector<Complex> result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) result[i] = Complex(r[i], j[i]);
    return result;
}

std::vector<Complex> solve(const std::vector<Complex>& a0,
                           const std::vector<Complex>& b0, int n) {
    auto a = a0; auto b = b0;
    for (int k = 0; k < n; ++k) {
        int pivot = k; double best = std::abs(a[static_cast<std::size_t>(k * n + k)]);
        for (int i = k + 1; i < n; ++i) {
            const double candidate = std::abs(a[static_cast<std::size_t>(i * n + k)]);
            if (candidate > best) { best = candidate; pivot = i; }
        }
        if (!(best > 0.0) || !std::isfinite(best))
            throw std::runtime_error("IH Eddy operator is singular");
        if (pivot != k) {
            for (int j = k; j < n; ++j)
                std::swap(a[static_cast<std::size_t>(k * n + j)], a[static_cast<std::size_t>(pivot * n + j)]);
            std::swap(b[static_cast<std::size_t>(k)], b[static_cast<std::size_t>(pivot)]);
        }
        for (int i = k + 1; i < n; ++i) {
            const Complex factor = a[static_cast<std::size_t>(i * n + k)] /
                                   a[static_cast<std::size_t>(k * n + k)];
            for (int j = k + 1; j < n; ++j)
                a[static_cast<std::size_t>(i * n + j)] -= factor * a[static_cast<std::size_t>(k * n + j)];
            b[static_cast<std::size_t>(i)] -= factor * b[static_cast<std::size_t>(k)];
        }
    }
    std::vector<Complex> x(static_cast<std::size_t>(n));
    for (int i = n - 1; i >= 0; --i) {
        Complex value = b[static_cast<std::size_t>(i)];
        for (int j = i + 1; j < n; ++j) value -= a[static_cast<std::size_t>(i * n + j)] * x[static_cast<std::size_t>(j)];
        x[static_cast<std::size_t>(i)] = value / a[static_cast<std::size_t>(i * n + i)];
    }
    return x;
}

Context* make_context(const mxArray* config) {
    if (!config || !mxIsStruct(config)) throw std::invalid_argument("IH Eddy requires a configuration struct");
    auto c = std::make_unique<Context>();
    c->n_unknown = positive_integer(
        field(config, "n_eddy_unknown"), "n_eddy_unknown");
    c->n_heat = positive_integer(field(config, "n_heat"), "n_heat");
    validate_bh_mode(field(config, "bh_mode"));
    if (c->n_unknown <= 0 || c->n_heat <= 0) throw std::invalid_argument("IH Eddy dimensions must be positive");
    c->matrix = complex_matrix(field(config, "eddy_matrix_real"), field(config, "eddy_matrix_imag"), c->n_unknown, c->n_unknown, "eddy_matrix");
    c->rhs_per_amp = complex_matrix(field(config, "eddy_rhs_real"), field(config, "eddy_rhs_imag"), c->n_unknown, 1, "eddy_rhs");
    c->heat_projection = real_array(field(config, "heat_projection"), "heat_projection");
    if (c->heat_projection.size() != static_cast<std::size_t>(c->n_heat * c->n_unknown))
        throw std::invalid_argument("heat_projection has the wrong shape");
    c->heat_weights = real_array(field(config, "heat_cell_weights"), "heat_cell_weights");
    if (c->heat_weights.size() != static_cast<std::size_t>(c->n_heat))
        throw std::invalid_argument("heat_cell_weights has the wrong shape");
    for (double weight : c->heat_weights)
        if (!(weight > 0.0) || !std::isfinite(weight))
            throw std::invalid_argument("heat_cell_weights must be finite and positive");
    c->cached_heat.assign(static_cast<std::size_t>(c->n_heat), 0.0);
    c->reference_heat.assign(static_cast<std::size_t>(c->n_heat), 0.0);
    const int n_temperature = positive_integer(
        field(config, "n_temperature"), "n_temperature");
    c->previous_temperature.assign(
        static_cast<std::size_t>(n_temperature), 0.0);
    c->temperature_weights = real_array(
        field(config, "temperature_cell_weights"), "temperature_cell_weights");
    if (c->temperature_weights.size() != c->previous_temperature.size())
        throw std::invalid_argument("temperature_cell_weights has the wrong shape");
    for (double weight : c->temperature_weights)
        if (!(weight > 0.0) || !std::isfinite(weight))
            throw std::invalid_argument(
                "temperature_cell_weights must be finite and positive");
    const std::string rotation_mode = string_value(
        field(config, "rotation_mode"), "rotation_mode", "none");
    if (rotation_mode == "periodic-uniform") c->periodic_rotation = true;
    else if (rotation_mode != "none")
        throw std::invalid_argument(
            "IH Eddy rotation_mode must be 'none' or 'periodic-uniform'");
    if (const mxArray* origin = field(config, "angle_origin_rad"))
        c->angle_origin_rad = scalar(origin, "angle_origin_rad");
    if (!std::isfinite(c->angle_origin_rad))
        throw std::invalid_argument("angle_origin_rad must be finite");
    if (const mxArray* reference = field(config, "bh_reference_temperature_K"))
        c->reference_temperature_K = scalar(reference, "bh_reference_temperature_K");
    const mxArray* slope_real = field(config, "eddy_matrix_temperature_slope_real");
    const mxArray* slope_imag = field(config, "eddy_matrix_temperature_slope_imag");
    if (slope_real || slope_imag) {
        c->matrix_temperature_slope = complex_matrix(
            slope_real, slope_imag, static_cast<int>(c->previous_temperature.size()),
            c->n_unknown * c->n_unknown, "eddy_matrix_temperature_slope");
    }
    return c.release();
}

std::vector<Complex> temperature_matrix(const Context& c, const double* temperature) {
    auto result = c.matrix;
    if (c.matrix_temperature_slope.empty()) return result;
    const std::size_t block = static_cast<std::size_t>(c.n_unknown * c.n_unknown);
    for (std::size_t t = 0; t < c.previous_temperature.size(); ++t) {
        const double delta = temperature[t] - c.reference_temperature_K;
        for (std::size_t k = 0; k < block; ++k)
            result[k] += delta * c.matrix_temperature_slope[t * block + k];
    }
    return result;
}

}  // namespace

static void mdlInitializeSizes(SimStruct* S) {
    ssSetNumSFcnParams(S, 1);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) return;
    const mxArray* config = ssGetSFcnParam(S, 0);
    if (!mxIsStruct(config)) { ssSetErrorStatus(S, "radia_ih_eddy_sfun: Parameters must be a struct."); return; }
    const auto positive_scalar = [config](const char* name) -> int {
        const mxArray* value = field(config, name);
        if (!value || mxGetNumberOfElements(value) != 1 || !mxIsNumeric(value))
            throw std::invalid_argument(std::string("missing scalar parameter: ") + name);
        const double number = mxGetScalar(value);
        if (!std::isfinite(number) || number <= 0.0 ||
            number > static_cast<double>(INT_MAX) ||
            std::floor(number) != number)
            throw std::invalid_argument(std::string("invalid positive parameter: ") + name);
        return static_cast<int>(number);
    };
    int n_temperature = 0;
    int n_heat = 0;
    try { n_temperature = positive_scalar("n_temperature"); n_heat = positive_scalar("n_heat"); }
    catch (const std::exception& error) { set_error_status(S, error); return; }
    if (!ssSetNumInputPorts(S, 3)) return;
    ssSetInputPortWidth(S, 0, 1);
    ssSetInputPortWidth(S, 1, 1);
    ssSetInputPortWidth(S, 2, n_temperature);
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortDirectFeedThrough(S, 1, 1);
    ssSetInputPortDirectFeedThrough(S, 2, 1);
    ssSetInputPortDataType(S, 0, SS_DOUBLE);
    ssSetInputPortDataType(S, 1, SS_DOUBLE);
    ssSetInputPortDataType(S, 2, SS_DOUBLE);
    ssSetInputPortRequiredContiguous(S, 0, 1);
    ssSetInputPortRequiredContiguous(S, 1, 1);
    ssSetInputPortRequiredContiguous(S, 2, 1);
    if (!ssSetNumOutputPorts(S, 1)) return;
    ssSetOutputPortWidth(S, 0, n_heat);
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);
    ssSetNumSampleTimes(S, 1);
    ssSetNumPWork(S, 1);
    ssSetNumDWork(S, 0);
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}

static void mdlInitializeSampleTimes(SimStruct* S) {
    const mxArray* config = ssGetSFcnParam(S, 0);
    const mxArray* sample_time = field(config, "sample_time_s");
    if (!sample_time || !mxIsNumeric(sample_time) || mxIsComplex(sample_time) ||
        mxGetNumberOfElements(sample_time) != 1 ||
        !(mxGetScalar(sample_time) > 0.0) || !std::isfinite(mxGetScalar(sample_time))) {
        ssSetErrorStatus(S, "radia_ih_eddy_sfun: sample_time_s must be positive.");
        return;
    }
    ssSetSampleTime(S, 0, mxGetScalar(sample_time));
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_START
static void mdlStart(SimStruct* S) {
    try {
        auto* c = make_context(ssGetSFcnParam(S, 0));
        ssSetOutputPortWidth(S, 0, c->n_heat);
        ssGetPWork(S)[0] = c;
    } catch (const std::exception& error) { set_error_status(S, error); }
}

static void mdlOutputs(SimStruct* S, int_T) {
    auto* c = static_cast<Context*>(ssGetPWork(S)[0]);
    if (!c) return;
    try {
        const double current = *static_cast<const real_T*>(ssGetInputPortSignal(S, 0));
        const double angle = *static_cast<const real_T*>(ssGetInputPortSignal(S, 1));
        const auto* temperature = static_cast<const real_T*>(ssGetInputPortSignal(S, 2));
        if (!std::isfinite(current) || !std::isfinite(angle)) throw std::invalid_argument("IH Eddy inputs must be finite");
        std::vector<double> temperature_at_angle;
        if (c->periodic_rotation) {
            radia::ih::transport_periodic(
                std::vector<double>(temperature,
                                    temperature + c->previous_temperature.size()),
                c->temperature_weights, angle - c->angle_origin_rad,
                temperature_at_angle);
        } else {
            if (!equivalent_periodic_angle(angle, c->angle_origin_rad))
                throw std::invalid_argument(
                    "IH Eddy received a changing angle but rotation_mode is 'none'");
            temperature_at_angle.assign(
                temperature, temperature + c->previous_temperature.size());
        }
        bool temperature_changed = !c->have_cache;
        for (std::size_t i = 0; i < c->previous_temperature.size(); ++i) {
            if (!std::isfinite(temperature_at_angle[i]))
                throw std::invalid_argument("IH Eddy temperature input must be finite");
            temperature_changed = temperature_changed ||
                temperature_at_angle[i] != c->previous_temperature[i];
        }
        const bool material_changed = !c->have_cache ||
            (!c->matrix_temperature_slope.empty() && temperature_changed);
        if (material_changed) {
            std::vector<Complex> rhs = c->rhs_per_amp;
            const auto matrix = temperature_matrix(*c, temperature_at_angle.data());
            const auto solution = solve(matrix, rhs, c->n_unknown);
            for (int i = 0; i < c->n_heat; ++i) {
                double q = 0.0;
                for (int j = 0; j < c->n_unknown; ++j) {
                    const double p = c->heat_projection[static_cast<std::size_t>(i * c->n_unknown + j)];
                    q += p * std::norm(solution[static_cast<std::size_t>(j)]);
                }
                c->reference_heat[static_cast<std::size_t>(i)] = q;
            }
            c->previous_temperature = temperature_at_angle;
            c->have_cache = true;
        }
        std::vector<double> heat_at_origin = c->reference_heat;
        const double scale = current * current;
        for (int i = 0; i < c->n_heat; ++i)
            heat_at_origin[static_cast<std::size_t>(i)] *= scale;
        if (c->periodic_rotation) {
            // The source is stationary. Express its heat distribution in
            // the rotating workpiece coordinate system.
            radia::ih::transport_periodic(
                heat_at_origin, c->heat_weights,
                -(angle - c->angle_origin_rad), c->cached_heat);
        } else {
            c->cached_heat = std::move(heat_at_origin);
        }
        auto* output = static_cast<real_T*>(ssGetOutputPortSignal(S, 0));
        std::copy(c->cached_heat.begin(), c->cached_heat.end(), output);
    } catch (const std::exception& error) { set_error_status(S, error); }
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
