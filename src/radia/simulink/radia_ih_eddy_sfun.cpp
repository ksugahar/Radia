#define S_FUNCTION_NAME radia_ih_eddy_sfun
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "mex.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <climits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Complex = std::complex<double>;

struct Context {
    int n_unknown = 0;
    int n_heat = 0;
    std::vector<Complex> matrix;
    // Optional temperature Jacobian, stored as [temperature][row][column].
    std::vector<Complex> matrix_temperature_slope;
    double reference_temperature_K = 293.15;
    std::vector<Complex> rhs_per_amp;
    std::vector<double> heat_projection;
    std::vector<double> cached_heat;
    std::vector<double> unit_current_heat;
    std::vector<double> previous_temperature;
    double previous_current = 0.0;
    bool have_cache = false;
    bool bh_linear = true;
    double previous_angle = 0.0;
};

const mxArray* field(const mxArray* a, const char* name) {
    return a && mxIsStruct(a) ? mxGetField(a, 0, name) : nullptr;
}

std::vector<double> real_array(const mxArray* a, const char* name) {
    if (!a || !mxIsDouble(a) || mxIsComplex(a))
        throw std::invalid_argument(std::string("IH Eddy field '") + name + "' must be real double");
    return std::vector<double>(mxGetDoubles(a), mxGetDoubles(a) + mxGetNumberOfElements(a));
}

double scalar(const mxArray* a, const char* name) {
    if (!a || !mxIsNumeric(a) || mxIsComplex(a) || mxGetNumberOfElements(a) != 1)
        throw std::invalid_argument(std::string("IH Eddy field '") + name + "' must be a scalar");
    return mxGetScalar(a);
}

bool string_is_linear(const mxArray* a) {
    if (!a) return true;
    if (mxIsChar(a)) {
        char* text = mxArrayToUTF8String(a);
        if (!text) throw std::invalid_argument("IH Eddy field 'bh_mode' is invalid");
        const bool linear = std::string(text) == "linear";
        mxFree(text);
        return linear;
    }
    throw std::invalid_argument("IH Eddy field 'bh_mode' must be 'linear' or 'nonlinear'");
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
    auto* c = new Context();
    c->n_unknown = static_cast<int>(scalar(field(config, "n_eddy_unknown"), "n_eddy_unknown"));
    c->n_heat = static_cast<int>(scalar(field(config, "n_heat"), "n_heat"));
    c->bh_linear = string_is_linear(field(config, "bh_mode"));
    if (c->n_unknown <= 0 || c->n_heat <= 0) throw std::invalid_argument("IH Eddy dimensions must be positive");
    c->matrix = complex_matrix(field(config, "eddy_matrix_real"), field(config, "eddy_matrix_imag"), c->n_unknown, c->n_unknown, "eddy_matrix");
    c->rhs_per_amp = complex_matrix(field(config, "eddy_rhs_real"), field(config, "eddy_rhs_imag"), c->n_unknown, 1, "eddy_rhs");
    c->heat_projection = real_array(field(config, "heat_projection"), "heat_projection");
    if (c->heat_projection.size() != static_cast<std::size_t>(c->n_heat * c->n_unknown))
        throw std::invalid_argument("heat_projection has the wrong shape");
    c->cached_heat.assign(static_cast<std::size_t>(c->n_heat), 0.0);
    c->unit_current_heat.assign(static_cast<std::size_t>(c->n_heat), 0.0);
    const mxArray* n_temperature = field(config, "n_temperature");
    if (!n_temperature || mxGetNumberOfElements(n_temperature) != 1)
        throw std::invalid_argument("IH Eddy requires n_temperature");
    c->previous_temperature.assign(static_cast<std::size_t>(mxGetScalar(n_temperature)), 0.0);
    if (const mxArray* reference = field(config, "bh_reference_temperature_K"))
        c->reference_temperature_K = scalar(reference, "bh_reference_temperature_K");
    const mxArray* slope_real = field(config, "eddy_matrix_temperature_slope_real");
    const mxArray* slope_imag = field(config, "eddy_matrix_temperature_slope_imag");
    if (slope_real || slope_imag) {
        c->matrix_temperature_slope = complex_matrix(
            slope_real, slope_imag, static_cast<int>(c->previous_temperature.size()),
            c->n_unknown * c->n_unknown, "eddy_matrix_temperature_slope");
    }
    return c;
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
        if (!std::isfinite(number) || number <= 0.0 || number > static_cast<double>(INT_MAX))
            throw std::invalid_argument(std::string("invalid positive parameter: ") + name);
        return static_cast<int>(number);
    };
    int n_temperature = 0;
    int n_heat = 0;
    try { n_temperature = positive_scalar("n_temperature"); n_heat = positive_scalar("n_heat"); }
    catch (const std::exception& error) { ssSetErrorStatus(S, error.what()); return; }
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
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_START
static void mdlStart(SimStruct* S) {
    try {
        auto* c = make_context(ssGetSFcnParam(S, 0));
        ssSetOutputPortWidth(S, 0, c->n_heat);
        ssGetPWork(S)[0] = c;
    } catch (const std::exception& error) { ssSetErrorStatus(S, error.what()); }
}

static void mdlOutputs(SimStruct* S, int_T) {
    auto* c = static_cast<Context*>(ssGetPWork(S)[0]);
    if (!c) return;
    try {
        const double current = *static_cast<const real_T*>(ssGetInputPortSignal(S, 0));
        const double angle = *static_cast<const real_T*>(ssGetInputPortSignal(S, 1));
        const auto* temperature = static_cast<const real_T*>(ssGetInputPortSignal(S, 2));
        if (!std::isfinite(current) || !std::isfinite(angle)) throw std::invalid_argument("IH Eddy inputs must be finite");
        bool temperature_changed = !c->have_cache;
        for (std::size_t i = 0; i < c->previous_temperature.size(); ++i) {
            if (!std::isfinite(temperature[i])) throw std::invalid_argument("IH Eddy temperature input must be finite");
            temperature_changed = temperature_changed || temperature[i] != c->previous_temperature[i];
        }
        const bool material_changed = !c->have_cache || temperature_changed;
        const bool current_requires_solve = !c->bh_linear &&
            (!c->have_cache || current != c->previous_current);
        if (material_changed || current_requires_solve) {
            std::vector<Complex> rhs = c->rhs_per_amp;
            const double solve_current = c->bh_linear ? 1.0 : current;
            for (auto& value : rhs) value *= solve_current;
            const auto matrix = temperature_matrix(*c, temperature);
            const auto solution = solve(matrix, rhs, c->n_unknown);
            for (int i = 0; i < c->n_heat; ++i) {
                double q = 0.0;
                for (int j = 0; j < c->n_unknown; ++j) {
                    const double p = c->heat_projection[static_cast<std::size_t>(i * c->n_unknown + j)];
                    q += p * std::norm(solution[static_cast<std::size_t>(j)]);
                }
                if (c->bh_linear) c->unit_current_heat[static_cast<std::size_t>(i)] = q;
                else c->cached_heat[static_cast<std::size_t>(i)] = q;
            }
            c->previous_current = current;
            std::copy(temperature, temperature + c->previous_temperature.size(), c->previous_temperature.begin());
            c->have_cache = true;
        }
        if (c->bh_linear) {
            const double scale = current * current;
            for (int i = 0; i < c->n_heat; ++i)
                c->cached_heat[static_cast<std::size_t>(i)] =
                    scale * c->unit_current_heat[static_cast<std::size_t>(i)];
        }
        c->previous_angle = angle;
        auto* output = static_cast<real_T*>(ssGetOutputPortSignal(S, 0));
        std::copy(c->cached_heat.begin(), c->cached_heat.end(), output);
    } catch (const std::exception& error) { ssSetErrorStatus(S, error.what()); }
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
