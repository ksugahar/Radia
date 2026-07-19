/*
 * Simulink C-MEX S-function adapter for the Radia angle-periodic motor ROM.
 *
 * The adapter deliberately owns a discrete internal time step.  This makes
 * the block a Co-Simulation-style plant block: Simulink supplies phase
 * voltages, load torque, and ambient temperature, while the Radia C ABI
 * advances the electromechanical and eddy-current state by one fixed step.
 * The numerical kernel remains rad_motor_rom_c.cpp; this file is only the
 * Simulink boundary and MATLAB-array marshaling layer.
 */
#define S_FUNCTION_NAME radia_motor_rom_sfun
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "mex.h"

#include "rad_motor_rom_c.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

enum ParameterIndex {
    kModelParam = 0,
    kSampleTimeParam = 1,
    kInitialStateParam = 2,
    kMaxIterationsParam = 3,
    kToleranceParam = 4,
    kParameterCount = 5,
};

struct MotorContext {
    RadMotorROMHandle* handle = nullptr;
    RadMotorROMModelData model{};
    RadMotorROMState state{};
    RadMotorROMInput input{};
    RadMotorROMStepOutput output{};

    std::size_t n_angle = 0;
    std::size_t n_phase = 0;
    std::size_t n_eddy = 0;
    std::size_t n_generalized = 0;
    double sample_time = 0.0;
    unsigned max_iterations = 30;
    double tolerance = 1.0e-11;
    double initial_time = 0.0;
    double initial_angle = 0.0;
    double initial_speed = 0.0;
    double initial_temperature = 293.15;
    std::vector<double> initial_currents;

    std::vector<double> inductance;
    std::vector<double> resistance;
    std::vector<double> pm_flux;
    std::vector<double> motion_flux;
    std::vector<double> cogging;
    std::vector<double> end_inductance;
    std::vector<double> end_resistance;
    std::vector<double> temperature_coefficient;
    std::vector<double> generalized_currents;
    std::vector<double> phase_voltage;
    std::vector<double> phase_flux;
    std::vector<double> y;
    std::string error;
};

const mxArray* field(const mxArray* object, const char* name) {
    if (object == nullptr || !mxIsStruct(object)) return nullptr;
    return mxGetField(object, 0, name);
}

bool numeric_scalar(const mxArray* value, double& result, std::string& error) {
    if (value == nullptr || !mxIsNumeric(value) || mxIsComplex(value) ||
        mxGetNumberOfElements(value) != 1) {
        error = "expected a real numeric scalar";
        return false;
    }
    result = mxGetScalar(value);
    if (!std::isfinite(result)) {
        error = "numeric scalar must be finite";
        return false;
    }
    return true;
}

bool numeric_scalar_allow_nan(const mxArray* value, double& result,
                              std::string& error) {
    if (value == nullptr || !mxIsNumeric(value) || mxIsComplex(value) ||
        mxGetNumberOfElements(value) != 1) {
        error = "expected a real numeric scalar";
        return false;
    }
    result = mxGetScalar(value);
    if (std::isinf(result)) {
        error = "numeric scalar must not be infinite";
        return false;
    }
    return true;
}

bool required_scalar(const mxArray* object, const char* name, double& result,
                     std::string& error) {
    const mxArray* value = field(object, name);
    if (value == nullptr) {
        error = std::string("missing field '") + name + "'";
        return false;
    }
    if (!numeric_scalar(value, result, error)) {
        error = std::string("field '") + name + "': " + error;
        return false;
    }
    return true;
}

bool required_scalar_allow_nan(const mxArray* object, const char* name,
                               double& result, std::string& error) {
    const mxArray* value = field(object, name);
    if (value == nullptr) {
        error = std::string("missing field '") + name + "'";
        return false;
    }
    if (!numeric_scalar_allow_nan(value, result, error)) {
        error = std::string("field '") + name + "': " + error;
        return false;
    }
    return true;
}

bool optional_scalar(const mxArray* object, const char* name, double fallback,
                     double& result, std::string& error) {
    const mxArray* value = field(object, name);
    if (value == nullptr) {
        result = fallback;
        return true;
    }
    if (!numeric_scalar(value, result, error)) {
        error = std::string("field '") + name + "': " + error;
        return false;
    }
    return true;
}

bool flag(const mxArray* object, const char* name, bool fallback,
          bool& result, std::string& error) {
    const mxArray* value = field(object, name);
    if (value == nullptr) {
        result = fallback;
        return true;
    }
    if (mxIsLogical(value) && mxGetNumberOfElements(value) == 1) {
        result = mxIsLogicalScalarTrue(value);
        return true;
    }
    double scalar = 0.0;
    if (!numeric_scalar(value, scalar, error)) {
        error = std::string("field '") + name + "': " + error;
        return false;
    }
    result = scalar != 0.0;
    return true;
}

bool real_double_array(const mxArray* value, const char* name,
                       std::string& error) {
    if (value == nullptr || !mxIsDouble(value) || mxIsComplex(value)) {
        error = std::string("field '") + name + "' must be a real double array";
        return false;
    }
    return true;
}

bool copy_vector(const mxArray* object, const char* name, std::size_t count,
                 std::vector<double>& target, std::string& error,
                 bool required = true) {
    const mxArray* value = field(object, name);
    if (value == nullptr) {
        if (!required) {
            target.clear();
            return true;
        }
        error = std::string("missing field '") + name + "'";
        return false;
    }
    if (!real_double_array(value, name, error) ||
        mxGetNumberOfElements(value) != count) {
        if (error.empty()) {
            std::ostringstream stream;
            stream << "field '" << name << "' must contain " << count
                   << " values";
            error = stream.str();
        }
        return false;
    }
    const double* source = mxGetPr(value);
    target.assign(source, source + count);
    return true;
}

bool copy_matrix(const mxArray* object, const char* name, std::size_t n,
                 std::vector<double>& target, std::string& error,
                 bool required = false) {
    const mxArray* value = field(object, name);
    if (value == nullptr) {
        if (!required) {
            target.clear();
            return true;
        }
        error = std::string("missing field '") + name + "'";
        return false;
    }
    if (!real_double_array(value, name, error) ||
        mxGetNumberOfElements(value) != n * n) {
        if (error.empty()) {
            std::ostringstream stream;
            stream << "field '" << name << "' must have shape [" << n << ", "
                   << n << "]";
            error = stream.str();
        }
        return false;
    }
    const double* source = mxGetPr(value);
    target.assign(n * n, 0.0);
    // MATLAB stores matrices column-major; the C ABI is row-major.
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            target[i * n + j] = source[i + n * j];
    return true;
}

bool copy_angle_matrix(const mxArray* object, const char* name,
                       std::size_t n_angle, std::size_t n,
                       std::vector<double>& target, std::string& error,
                       bool required = true) {
    const mxArray* value = field(object, name);
    if (value == nullptr) {
        if (!required) {
            target.clear();
            return true;
        }
        error = std::string("missing field '") + name + "'";
        return false;
    }
    if (!real_double_array(value, name, error) ||
        mxGetNumberOfElements(value) != n_angle * n * n) {
        if (error.empty()) {
            std::ostringstream stream;
            stream << "field '" << name << "' must have shape [" << n_angle
                   << ", " << n << ", " << n << "]";
            error = stream.str();
        }
        return false;
    }
    const double* source = mxGetPr(value);
    target.assign(n_angle * n * n, 0.0);
    // MATLAB index for [angle, row, column] is angle + na*(row + n*column).
    for (std::size_t a = 0; a < n_angle; ++a)
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                target[(a * n + i) * n + j] = source[a + n_angle * (i + n * j)];
    return true;
}

bool copy_angle_vector_field(const mxArray* object, const char* name,
                             std::size_t n_angle, std::size_t n,
                             std::vector<double>& target, std::string& error,
                             bool required = true) {
    const mxArray* value = field(object, name);
    if (value == nullptr) {
        if (!required) {
            target.clear();
            return true;
        }
        error = std::string("missing field '") + name + "'";
        return false;
    }
    if (!real_double_array(value, name, error) ||
        mxGetNumberOfElements(value) != n_angle * n) {
        if (error.empty()) {
            std::ostringstream stream;
            stream << "field '" << name << "' must have shape [" << n_angle
                   << ", " << n << "]";
            error = stream.str();
        }
        return false;
    }
    const double* source = mxGetPr(value);
    target.assign(n_angle * n, 0.0);
    for (std::size_t a = 0; a < n_angle; ++a)
        for (std::size_t i = 0; i < n; ++i)
            target[a * n + i] = source[a + n_angle * i];
    return true;
}

bool load_model(MotorContext& ctx, const mxArray* model, std::string& error) {
    if (model == nullptr || !mxIsStruct(model)) {
        error = "the first S-function parameter must be a RadiaMotorROM struct";
        return false;
    }

    double scalar = 0.0;
    if (!required_scalar(model, "n_phase", scalar, error) || scalar < 1.0 ||
        std::floor(scalar) != scalar)
        return false;
    ctx.n_phase = static_cast<std::size_t>(scalar);
    if (!required_scalar(model, "n_generalized", scalar, error) ||
        scalar < static_cast<double>(ctx.n_phase) || std::floor(scalar) != scalar)
        return false;
    ctx.n_generalized = static_cast<std::size_t>(scalar);
    ctx.n_eddy = ctx.n_generalized - ctx.n_phase;

    const mxArray* inductance = field(model, "inductance_H");
    if (!real_double_array(inductance, "inductance_H", error) ||
        mxGetNumberOfDimensions(inductance) != 3) {
        if (error.empty()) error = "inductance_H must be a real 3-D double array";
        return false;
    }
    const mwSize* dimensions = mxGetDimensions(inductance);
    ctx.n_angle = static_cast<std::size_t>(dimensions[0]);
    if (dimensions[1] != ctx.n_generalized || dimensions[2] != ctx.n_generalized ||
        ctx.n_angle < 3 || ctx.n_angle % 2 == 0)
        return false;

    bool has_motion = false;
    bool has_cogging = false;
    bool requires_hysteresis = false;
    if (!flag(model, "has_motional_v_cross_b", false, has_motion, error) ||
        !flag(model, "has_cogging_coenergy", false, has_cogging, error) ||
        !flag(model, "external_hysteresis_required", false, requires_hysteresis, error))
        return false;
    if (requires_hysteresis) {
        error = "external hysteresis callbacks are not available in the Simulink C-MEX adapter";
        return false;
    }

    if (!copy_angle_matrix(model, "inductance_H", ctx.n_angle, ctx.n_generalized,
                           ctx.inductance, error) ||
        !copy_angle_matrix(model, "resistance_ohm", ctx.n_angle, ctx.n_generalized,
                           ctx.resistance, error) ||
        !copy_angle_vector_field(model, "pm_flux_linkage_Wb", ctx.n_angle,
                                 ctx.n_generalized, ctx.pm_flux, error) ||
        !copy_matrix(model, "end_winding_inductance_H", ctx.n_generalized,
                     ctx.end_inductance, error) ||
        !copy_matrix(model, "end_winding_resistance_ohm", ctx.n_generalized,
                     ctx.end_resistance, error) ||
        !copy_vector(model, "resistance_temperature_coefficient_per_K",
                     ctx.n_generalized, ctx.temperature_coefficient, error))
        return false;

    if (has_motion &&
        !copy_angle_vector_field(model, "motion_flux_gradient_Wb_per_rad",
                                 ctx.n_angle, ctx.n_generalized, ctx.motion_flux,
                                 error))
        return false;
    if (has_cogging &&
        !copy_angle_vector_field(model, "cogging_coenergy_J", ctx.n_angle, 1,
                                 ctx.cogging, error))
        return false;

    double period = 0.0;
    double origin = 0.0;
    double skew = 0.0;
    double inertia = 0.0;
    double friction = 0.0;
    double reference_temperature = 0.0;
    double thermal_capacity = 0.0;
    double thermal_conductance = 0.0;
    if (!required_scalar(model, "angle_origin_rad", origin, error) ||
        !required_scalar(model, "period_rad", period, error) ||
        !required_scalar(model, "skew_span_rad", skew, error) ||
        !required_scalar(model, "inertia_kg_m2", inertia, error) ||
        !required_scalar(model, "viscous_friction_Nm_s", friction, error) ||
        !required_scalar(model, "reference_temperature_K", reference_temperature, error) ||
        !required_scalar_allow_nan(model, "thermal_capacity_J_per_K",
                                   thermal_capacity, error) ||
        !required_scalar(model, "thermal_conductance_W_per_K", thermal_conductance, error))
        return false;
    if (std::isnan(thermal_capacity)) thermal_capacity = 0.0;

    ctx.model.abi_version = RAD_MOTOR_ROM_ABI_VERSION;
    ctx.model.n_angle_samples = ctx.n_angle;
    ctx.model.n_phase = ctx.n_phase;
    ctx.model.n_generalized = ctx.n_generalized;
    ctx.model.angle_origin_rad = origin;
    ctx.model.period_rad = period;
    ctx.model.skew_span_rad = skew;
    ctx.model.inertia_kg_m2 = inertia;
    ctx.model.viscous_friction_Nm_s = friction;
    ctx.model.reference_temperature_K = reference_temperature;
    ctx.model.thermal_capacity_J_per_K = thermal_capacity;
    ctx.model.thermal_conductance_W_per_K = thermal_conductance;
    ctx.model.inductance_H = ctx.inductance.data();
    ctx.model.resistance_ohm = ctx.resistance.data();
    ctx.model.pm_flux_linkage_Wb = ctx.pm_flux.data();
    ctx.model.cogging_coenergy_J = has_cogging ? ctx.cogging.data() : nullptr;
    ctx.model.motion_flux_gradient_Wb_per_rad = has_motion ? ctx.motion_flux.data() : nullptr;
    ctx.model.end_winding_inductance_H = ctx.end_inductance.empty()
        ? nullptr : ctx.end_inductance.data();
    ctx.model.end_winding_resistance_ohm = ctx.end_resistance.empty()
        ? nullptr : ctx.end_resistance.data();
    ctx.model.resistance_temperature_coefficient_per_K =
        ctx.temperature_coefficient.data();
    ctx.model.hysteresis_trial = nullptr;
    ctx.model.hysteresis_user_data = nullptr;
    return true;
}

bool initialize_state(MotorContext& ctx, const mxArray* initial, std::string& error) {
    if (initial == nullptr || !mxIsStruct(initial)) {
        error = "the third S-function parameter must be an initial-state struct";
        return false;
    }
    double value = 0.0;
    if (!optional_scalar(initial, "time_s", 0.0, value, error)) return false;
    ctx.state.time_s = value;
    if (!optional_scalar(initial, "rotor_angle_rad", 0.0, value, error)) return false;
    ctx.state.rotor_angle_rad = value;
    if (!optional_scalar(initial, "rotor_speed_rad_s", 0.0, value, error)) return false;
    ctx.state.rotor_speed_rad_s = value;
    if (!optional_scalar(initial, "temperature_K", ctx.model.reference_temperature_K,
                        value, error))
        return false;
    ctx.state.temperature_K = value;
    if (!copy_vector(initial, "generalized_currents_A", ctx.n_generalized,
                     ctx.generalized_currents, error, false))
        return false;
    if (ctx.generalized_currents.empty())
        ctx.generalized_currents.assign(ctx.n_generalized, 0.0);
    ctx.initial_time = ctx.state.time_s;
    ctx.initial_angle = ctx.state.rotor_angle_rad;
    ctx.initial_speed = ctx.state.rotor_speed_rad_s;
    ctx.initial_temperature = ctx.state.temperature_K;
    ctx.initial_currents = ctx.generalized_currents;
    ctx.state.abi_version = RAD_MOTOR_ROM_ABI_VERSION;
    ctx.state.generalized_currents_A = ctx.generalized_currents.data();
    ctx.state.hysteresis_flux_linkage_Wb = nullptr;
    ctx.state.hysteresis_stored_energy_J = 0.0;
    return true;
}

void refresh_output(MotorContext& ctx) {
    std::fill(ctx.y.begin(), ctx.y.end(), 0.0);
    std::size_t offset = 0;
    for (std::size_t i = 0; i < ctx.n_phase; ++i)
        ctx.y[offset++] = ctx.generalized_currents[i];
    for (std::size_t i = ctx.n_phase; i < ctx.n_generalized; ++i)
        ctx.y[offset++] = ctx.generalized_currents[i];
    for (std::size_t i = 0; i < ctx.n_phase; ++i)
        ctx.y[offset++] = ctx.phase_flux[i];
    ctx.y[offset++] = ctx.state.rotor_angle_rad;
    ctx.y[offset++] = ctx.state.rotor_speed_rad_s;
    ctx.y[offset++] = ctx.output.electromagnetic_torque_Nm;
    ctx.y[offset++] = ctx.output.resistive_loss_W;
    ctx.y[offset++] = ctx.output.hysteresis_loss_W;
    ctx.y[offset++] = ctx.state.temperature_K;
    ctx.y[offset++] = ctx.output.energy_balance_residual_W;
    ctx.y[offset] = static_cast<double>(ctx.output.nonlinear_iterations);
}

void set_error(SimStruct* S, MotorContext* ctx, const std::string& message) {
    if (ctx != nullptr) ctx->error = message;
    ssSetErrorStatus(S, ctx == nullptr ? message.c_str() : ctx->error.c_str());
}

}  // namespace

static void mdlInitializeSizes(SimStruct* S) {
    ssSetNumSFcnParams(S, kParameterCount);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) return;

    const mxArray* model = ssGetSFcnParam(S, kModelParam);
    double n_phase = 0.0;
    double n_generalized = 0.0;
    std::string error;
    if (!required_scalar(model, "n_phase", n_phase, error) ||
        !required_scalar(model, "n_generalized", n_generalized, error) ||
        n_phase < 1.0 || n_generalized < n_phase ||
        std::floor(n_phase) != n_phase || std::floor(n_generalized) != n_generalized) {
        ssSetErrorStatus(S, "RadiaMotorROM must provide positive n_phase and n_generalized fields");
        return;
    }
    double sample_time = 0.0;
    double max_iterations = 0.0;
    double tolerance = 0.0;
    if (!numeric_scalar(ssGetSFcnParam(S, kSampleTimeParam), sample_time, error) ||
        !numeric_scalar(ssGetSFcnParam(S, kMaxIterationsParam), max_iterations, error) ||
        !numeric_scalar(ssGetSFcnParam(S, kToleranceParam), tolerance, error) ||
        sample_time <= 0.0 || max_iterations < 1.0 ||
        std::floor(max_iterations) != max_iterations || tolerance <= 0.0) {
        ssSetErrorStatus(S, "Invalid Radia motor sample time, iteration count, or tolerance");
        return;
    }

    ssSetSFcnParamTunable(S, kModelParam, 0);
    ssSetSFcnParamTunable(S, kSampleTimeParam, 0);
    ssSetSFcnParamTunable(S, kInitialStateParam, 0);
    ssSetSFcnParamTunable(S, kMaxIterationsParam, 0);
    ssSetSFcnParamTunable(S, kToleranceParam, 0);

    if (!ssSetNumInputPorts(S, 1)) return;
    ssSetInputPortWidth(S, 0, static_cast<int_T>(n_phase + 2.0));
    ssSetInputPortDataType(S, 0, SS_DOUBLE);
    // Outputs expose the previously accepted state; the current input is
    // consumed by mdlUpdate, so this block has no algebraic direct feedthrough.
    ssSetInputPortDirectFeedThrough(S, 0, 0);
    ssSetInputPortRequiredContiguous(S, 0, 1);

    if (!ssSetNumOutputPorts(S, 1)) return;
    const int_T output_width = static_cast<int_T>(2.0 * n_phase +
                                                   (n_generalized - n_phase) + 8.0);
    ssSetOutputPortWidth(S, 0, output_width);
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);

    ssSetNumSampleTimes(S, 1);
    ssSetNumPWork(S, 1);
    ssSetNumDWork(S, 0);
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}

static void mdlInitializeSampleTimes(SimStruct* S) {
    const double sample_time = mxGetScalar(ssGetSFcnParam(S, kSampleTimeParam));
    ssSetSampleTime(S, 0, sample_time);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_START
static void mdlStart(SimStruct* S) {
    MotorContext* ctx = nullptr;
    try {
        ctx = new MotorContext();
        ctx->sample_time = mxGetScalar(ssGetSFcnParam(S, kSampleTimeParam));
        ctx->max_iterations = static_cast<unsigned>(
            mxGetScalar(ssGetSFcnParam(S, kMaxIterationsParam)));
        ctx->tolerance = mxGetScalar(ssGetSFcnParam(S, kToleranceParam));
        std::string error;
        if (!load_model(*ctx, ssGetSFcnParam(S, kModelParam), error) ||
            !initialize_state(*ctx, ssGetSFcnParam(S, kInitialStateParam), error)) {
            set_error(S, ctx, error);
            delete ctx;
            return;
        }
        ctx->phase_voltage.assign(ctx->n_phase, 0.0);
        ctx->phase_flux.assign(ctx->n_phase, 0.0);
        ctx->y.assign(2 * ctx->n_phase + ctx->n_eddy + 8, 0.0);
        ctx->input.abi_version = RAD_MOTOR_ROM_ABI_VERSION;
        ctx->input.phase_voltages_V = ctx->phase_voltage.data();
        ctx->input.has_ambient_temperature = 1;
        ctx->output.abi_version = RAD_MOTOR_ROM_ABI_VERSION;
        ctx->output.phase_flux_linkage_Wb = ctx->phase_flux.data();
        ctx->output.speed_voltage_V = nullptr;
        ctx->handle = rad_motor_rom_create(&ctx->model);
        if (ctx->handle == nullptr) {
            set_error(S, ctx, rad_motor_rom_last_error(nullptr));
            delete ctx;
            return;
        }
        refresh_output(*ctx);
        ssGetPWork(S)[0] = ctx;
    } catch (const std::exception& exception) {
        set_error(S, ctx, exception.what());
        if (ctx != nullptr) {
            rad_motor_rom_destroy(ctx->handle);
            delete ctx;
        }
    } catch (...) {
        set_error(S, ctx, "unknown exception while initializing Radia motor S-function");
        if (ctx != nullptr) {
            rad_motor_rom_destroy(ctx->handle);
            delete ctx;
        }
    }
}

static void mdlInitializeConditions(SimStruct* S) {
    auto* ctx = static_cast<MotorContext*>(ssGetPWork(S)[0]);
    if (ctx == nullptr) return;
    ctx->state.time_s = ctx->initial_time;
    ctx->state.rotor_angle_rad = ctx->initial_angle;
    ctx->state.rotor_speed_rad_s = ctx->initial_speed;
    ctx->state.temperature_K = ctx->initial_temperature;
    ctx->generalized_currents = ctx->initial_currents;
    ctx->state.generalized_currents_A = ctx->generalized_currents.data();
    std::fill(ctx->phase_flux.begin(), ctx->phase_flux.end(), 0.0);
    ctx->output = RadMotorROMStepOutput{};
    ctx->output.abi_version = RAD_MOTOR_ROM_ABI_VERSION;
    ctx->output.phase_flux_linkage_Wb = ctx->phase_flux.data();
    refresh_output(*ctx);
}

static void mdlOutputs(SimStruct* S, int_T /*tid*/) {
    auto* ctx = static_cast<MotorContext*>(ssGetPWork(S)[0]);
    if (ctx == nullptr) return;
    real_T* y = static_cast<real_T*>(ssGetOutputPortSignal(S, 0));
    std::copy(ctx->y.begin(), ctx->y.end(), y);
}

#define MDL_UPDATE
static void mdlUpdate(SimStruct* S, int_T /*tid*/) {
    auto* ctx = static_cast<MotorContext*>(ssGetPWork(S)[0]);
    if (ctx == nullptr) return;
    const real_T* u = static_cast<const real_T*>(ssGetInputPortSignal(S, 0));
    std::copy(u, u + ctx->n_phase, ctx->phase_voltage.begin());
    ctx->input.load_torque_Nm = u[ctx->n_phase];
    ctx->input.ambient_temperature_K = u[ctx->n_phase + 1];
    const int status = rad_motor_rom_step(
        ctx->handle, &ctx->state, &ctx->input, ctx->sample_time,
        ctx->max_iterations, ctx->tolerance, &ctx->output);
    if (status != RAD_MOTOR_ROM_OK) {
        std::ostringstream stream;
        stream << "Radia motor ROM step failed (status " << status << "): "
               << rad_motor_rom_last_error(ctx->handle);
        set_error(S, ctx, stream.str());
        return;
    }
    refresh_output(*ctx);
}

static void mdlTerminate(SimStruct* S) {
    auto* ctx = static_cast<MotorContext*>(ssGetPWork(S)[0]);
    if (ctx == nullptr) return;
    rad_motor_rom_destroy(ctx->handle);
    delete ctx;
    ssGetPWork(S)[0] = nullptr;
}

#ifdef MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif
