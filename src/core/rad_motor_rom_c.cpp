#define RAD_MOTOR_ROM_BUILD
#include "rad_motor_rom_c.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double kTwoPi = 6.283185307179586476925286766559;
thread_local std::string g_last_error;

bool finite(double value) { return std::isfinite(value) != 0; }

bool all_finite(const double* values, size_t count) {
    if (values == nullptr) return false;
    for (size_t i = 0; i < count; ++i)
        if (!finite(values[i])) return false;
    return true;
}

double sinc(double value) {
    return std::abs(value) < 1.0e-14 ? 1.0 : std::sin(value) / value;
}

class FourierTable {
public:
    FourierTable() = default;

    FourierTable(const double* samples, size_t count, size_t width,
                 double origin, double period)
        : count_(count), width_(width), harmonics_((count - 1) / 2),
          origin_(origin), period_(period), constant_(width, 0.0),
          cosine_(harmonics_ * width, 0.0), sine_(harmonics_ * width, 0.0) {
        if (samples == nullptr) return;
        const double inv_count = 1.0 / static_cast<double>(count_);
        for (size_t j = 0; j < count_; ++j) {
            const double sample_phase = kTwoPi * static_cast<double>(j) * inv_count;
            const double* row = samples + j * width_;
            for (size_t c = 0; c < width_; ++c) constant_[c] += row[c] * inv_count;
            for (size_t h = 1; h <= harmonics_; ++h) {
                const double co = 2.0 * inv_count * std::cos(static_cast<double>(h) * sample_phase);
                const double si = 2.0 * inv_count * std::sin(static_cast<double>(h) * sample_phase);
                double* ac = cosine_.data() + (h - 1) * width_;
                double* bs = sine_.data() + (h - 1) * width_;
                for (size_t c = 0; c < width_; ++c) {
                    ac[c] += row[c] * co;
                    bs[c] += row[c] * si;
                }
            }
        }
    }

    void evaluate(double angle, int derivative, double skew, double* output) const {
        std::fill(output, output + width_, 0.0);
        if (derivative == 0) std::copy(constant_.begin(), constant_.end(), output);
        const double base_frequency = kTwoPi / period_;
        const double phase = base_frequency * (angle - origin_);
        for (size_t h = 1; h <= harmonics_; ++h) {
            const double frequency = base_frequency * static_cast<double>(h);
            const double co = std::cos(static_cast<double>(h) * phase);
            const double si = std::sin(static_cast<double>(h) * phase);
            const double skew_factor = sinc(0.5 * frequency * skew);
            const double* ac = cosine_.data() + (h - 1) * width_;
            const double* bs = sine_.data() + (h - 1) * width_;
            if (derivative == 0) {
                for (size_t c = 0; c < width_; ++c)
                    output[c] += skew_factor * (ac[c] * co + bs[c] * si);
            } else {
                for (size_t c = 0; c < width_; ++c)
                    output[c] += skew_factor * frequency * (-ac[c] * si + bs[c] * co);
            }
        }
    }

private:
    size_t count_ = 0;
    size_t width_ = 0;
    size_t harmonics_ = 0;
    double origin_ = 0.0;
    double period_ = kTwoPi;
    std::vector<double> constant_;
    std::vector<double> cosine_;
    std::vector<double> sine_;
};

double dot(const double* left, const double* right, size_t n) {
    double value = 0.0;
    for (size_t i = 0; i < n; ++i) value += left[i] * right[i];
    return value;
}

void matvec(const double* matrix, const double* vector, size_t n, double* output) {
    for (size_t i = 0; i < n; ++i) {
        double value = 0.0;
        for (size_t j = 0; j < n; ++j) value += matrix[i * n + j] * vector[j];
        output[i] = value;
    }
}

void symmetrize(double* matrix, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const double value = 0.5 * (matrix[i * n + j] + matrix[j * n + i]);
            matrix[i * n + j] = value;
            matrix[j * n + i] = value;
        }
    }
}

bool solve_dense(std::vector<double> matrix, std::vector<double> rhs,
                 size_t n, std::vector<double>& solution) {
    const double epsilon = 64.0 * std::numeric_limits<double>::epsilon();
    for (size_t k = 0; k < n; ++k) {
        size_t pivot = k;
        double pivot_abs = std::abs(matrix[k * n + k]);
        for (size_t i = k + 1; i < n; ++i) {
            const double candidate = std::abs(matrix[i * n + k]);
            if (candidate > pivot_abs) {
                pivot = i;
                pivot_abs = candidate;
            }
        }
        if (!(pivot_abs > epsilon) || !finite(pivot_abs)) return false;
        if (pivot != k) {
            for (size_t j = k; j < n; ++j)
                std::swap(matrix[k * n + j], matrix[pivot * n + j]);
            std::swap(rhs[k], rhs[pivot]);
        }
        for (size_t i = k + 1; i < n; ++i) {
            const double factor = matrix[i * n + k] / matrix[k * n + k];
            matrix[i * n + k] = 0.0;
            for (size_t j = k + 1; j < n; ++j)
                matrix[i * n + j] -= factor * matrix[k * n + j];
            rhs[i] -= factor * rhs[k];
        }
    }
    solution.assign(n, 0.0);
    for (size_t ii = n; ii-- > 0;) {
        double value = rhs[ii];
        for (size_t j = ii + 1; j < n; ++j)
            value -= matrix[ii * n + j] * solution[j];
        solution[ii] = value / matrix[ii * n + ii];
    }
    return true;
}

}  // namespace

struct RadMotorROMHandle {
    size_t na = 0;
    size_t np = 0;
    size_t n = 0;
    double period = kTwoPi;
    double skew = 0.0;
    double inertia = 0.0;
    double friction = 0.0;
    double reference_temperature = 293.15;
    double thermal_capacity = 0.0;
    double thermal_conductance = 0.0;
    FourierTable inductance;
    FourierTable resistance;
    FourierTable pm_flux;
    FourierTable cogging_coenergy;
    FourierTable motion_flux;
    bool has_cogging = false;
    bool has_motion = false;
    std::vector<double> end_inductance;
    std::vector<double> end_resistance;
    std::vector<double> temperature_coefficient;
    RadMotorROMHysteresisTrial hysteresis_trial = nullptr;
    void* hysteresis_user_data = nullptr;
    std::string last_error;

    void set_error(const std::string& value) {
        last_error = value;
        g_last_error = value;
    }

    void eval_inductance(double angle, int derivative, std::vector<double>& result) const {
        inductance.evaluate(angle, derivative, skew, result.data());
        if (derivative == 0)
            for (size_t i = 0; i < n * n; ++i) result[i] += end_inductance[i];
        symmetrize(result.data(), n);
    }

    void eval_resistance(double angle, double temperature, std::vector<double>& result) const {
        resistance.evaluate(angle, 0, skew, result.data());
        for (size_t i = 0; i < n * n; ++i) result[i] += end_resistance[i];
        std::vector<double> scale(n);
        for (size_t i = 0; i < n; ++i) {
            const double fi = 1.0 + temperature_coefficient[i] *
                (temperature - reference_temperature);
            if (!finite(fi) || fi <= 0.0)
                throw std::runtime_error(
                    "temperature-dependent resistance scale must remain positive");
            scale[i] = std::sqrt(fi);
        }
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                result[i * n + j] *= scale[i] * scale[j];
            }
        }
        symmetrize(result.data(), n);
    }

    void eval_pm(double angle, int derivative, std::vector<double>& result) const {
        pm_flux.evaluate(angle, derivative, skew, result.data());
    }

    void eval_motion(double angle, std::vector<double>& result) const {
        if (has_motion) motion_flux.evaluate(angle, 0, skew, result.data());
        else std::fill(result.begin(), result.end(), 0.0);
    }

    double eval_cogging(double angle, int derivative) const {
        if (!has_cogging) return 0.0;
        double value = 0.0;
        cogging_coenergy.evaluate(angle, derivative, skew, &value);
        return value;
    }
};

extern "C" {

uint32_t rad_motor_rom_abi_version(void) { return RAD_MOTOR_ROM_ABI_VERSION; }

RadMotorROMHandle* rad_motor_rom_create(const RadMotorROMModelData* model) {
    g_last_error.clear();
    if (model == nullptr || model->abi_version != RAD_MOTOR_ROM_ABI_VERSION) {
        g_last_error = "model is null or has an unsupported ABI version";
        return nullptr;
    }
    const size_t n = model->n_generalized;
    const size_t na = model->n_angle_samples;
    if (na < 3 || na % 2 == 0 ||
        model->n_phase == 0 || model->n_generalized < model->n_phase ||
        n > std::numeric_limits<size_t>::max() / n ||
        !finite(model->period_rad) || model->period_rad <= 0.0 ||
        !finite(model->angle_origin_rad) || !finite(model->skew_span_rad) ||
        !finite(model->inertia_kg_m2) || model->inertia_kg_m2 <= 0.0 ||
        !finite(model->viscous_friction_Nm_s) || model->viscous_friction_Nm_s < 0.0 ||
        !finite(model->reference_temperature_K) ||
        !finite(model->thermal_capacity_J_per_K) || model->thermal_capacity_J_per_K < 0.0 ||
        !finite(model->thermal_conductance_W_per_K) ||
        model->thermal_conductance_W_per_K < 0.0 ||
        model->inductance_H == nullptr || model->resistance_ohm == nullptr ||
        model->pm_flux_linkage_Wb == nullptr) {
        g_last_error = "invalid dimensions, physical constants, or required table pointers";
        return nullptr;
    }
    const size_t nn = n * n;
    if (na > std::numeric_limits<size_t>::max() / nn ||
        !all_finite(model->inductance_H, na * nn) ||
        !all_finite(model->resistance_ohm, na * nn) ||
        !all_finite(model->pm_flux_linkage_Wb, na * n) ||
        (model->cogging_coenergy_J != nullptr &&
         !all_finite(model->cogging_coenergy_J, na)) ||
        (model->motion_flux_gradient_Wb_per_rad != nullptr &&
         !all_finite(model->motion_flux_gradient_Wb_per_rad, na * n)) ||
        (model->end_winding_inductance_H != nullptr &&
         !all_finite(model->end_winding_inductance_H, nn)) ||
        (model->end_winding_resistance_ohm != nullptr &&
         !all_finite(model->end_winding_resistance_ohm, nn)) ||
        (model->resistance_temperature_coefficient_per_K != nullptr &&
         !all_finite(model->resistance_temperature_coefficient_per_K, n))) {
        g_last_error = "motor ROM tables must contain only finite values";
        return nullptr;
    }
    try {
        std::unique_ptr<RadMotorROMHandle> handle(new RadMotorROMHandle());
        handle->na = model->n_angle_samples;
        handle->np = model->n_phase;
        handle->n = model->n_generalized;
        handle->period = model->period_rad;
        handle->skew = model->skew_span_rad;
        handle->inertia = model->inertia_kg_m2;
        handle->friction = model->viscous_friction_Nm_s;
        handle->reference_temperature = model->reference_temperature_K;
        handle->thermal_capacity = model->thermal_capacity_J_per_K;
        handle->thermal_conductance = model->thermal_conductance_W_per_K;
        handle->inductance = FourierTable(
            model->inductance_H, handle->na, handle->n * handle->n,
            model->angle_origin_rad, handle->period);
        handle->resistance = FourierTable(
            model->resistance_ohm, handle->na, handle->n * handle->n,
            model->angle_origin_rad, handle->period);
        handle->pm_flux = FourierTable(
            model->pm_flux_linkage_Wb, handle->na, handle->n,
            model->angle_origin_rad, handle->period);
        handle->has_cogging = model->cogging_coenergy_J != nullptr;
        handle->cogging_coenergy = FourierTable(
            model->cogging_coenergy_J, handle->na, 1,
            model->angle_origin_rad, handle->period);
        handle->has_motion = model->motion_flux_gradient_Wb_per_rad != nullptr;
        handle->motion_flux = FourierTable(
            model->motion_flux_gradient_Wb_per_rad, handle->na, handle->n,
            model->angle_origin_rad, handle->period);
        handle->end_inductance.assign(handle->n * handle->n, 0.0);
        handle->end_resistance.assign(handle->n * handle->n, 0.0);
        handle->temperature_coefficient.assign(handle->n, 0.0);
        if (model->end_winding_inductance_H != nullptr)
            std::copy(model->end_winding_inductance_H,
                      model->end_winding_inductance_H + handle->n * handle->n,
                      handle->end_inductance.begin());
        if (model->end_winding_resistance_ohm != nullptr)
            std::copy(model->end_winding_resistance_ohm,
                      model->end_winding_resistance_ohm + handle->n * handle->n,
                      handle->end_resistance.begin());
        if (model->resistance_temperature_coefficient_per_K != nullptr)
            std::copy(model->resistance_temperature_coefficient_per_K,
                      model->resistance_temperature_coefficient_per_K + handle->n,
                      handle->temperature_coefficient.begin());
        handle->hysteresis_trial = model->hysteresis_trial;
        handle->hysteresis_user_data = model->hysteresis_user_data;
        return handle.release();
    } catch (const std::exception& error) {
        g_last_error = error.what();
    } catch (...) {
        g_last_error = "unknown error while constructing motor ROM";
    }
    return nullptr;
}

void rad_motor_rom_destroy(RadMotorROMHandle* handle) { delete handle; }

const char* rad_motor_rom_last_error(const RadMotorROMHandle* handle) {
    return handle == nullptr ? g_last_error.c_str() : handle->last_error.c_str();
}

int rad_motor_rom_step(RadMotorROMHandle* h, RadMotorROMState* state,
                       const RadMotorROMInput* input, double dt,
                       unsigned max_iterations, double tolerance,
                       RadMotorROMStepOutput* output) {
    if (h == nullptr || state == nullptr || input == nullptr || output == nullptr ||
        state->abi_version != RAD_MOTOR_ROM_ABI_VERSION ||
        input->abi_version != RAD_MOTOR_ROM_ABI_VERSION ||
        output->abi_version != RAD_MOTOR_ROM_ABI_VERSION ||
        state->generalized_currents_A == nullptr || input->phase_voltages_V == nullptr ||
        !finite(dt) || dt <= 0.0 || max_iterations == 0 ||
        !finite(tolerance) || tolerance <= 0.0) {
        if (h != nullptr) h->set_error("invalid step argument or ABI version");
        else g_last_error = "null motor ROM handle";
        return RAD_MOTOR_ROM_INVALID_ARGUMENT;
    }

    try {
        const size_t n = h->n;
        const double theta0 = state->rotor_angle_rad;
        const double omega0 = state->rotor_speed_rad_s;
        const double temperature0 = state->temperature_K;
        const double ambient = input->has_ambient_temperature
            ? input->ambient_temperature_K : temperature0;
        if (!finite(state->time_s) || !finite(theta0) || !finite(omega0) ||
            !finite(temperature0) || !finite(state->hysteresis_stored_energy_J) ||
            !finite(input->load_torque_Nm) || !finite(ambient) ||
            !all_finite(state->generalized_currents_A, n) ||
            !all_finite(input->phase_voltages_V, h->np) ||
            (state->hysteresis_flux_linkage_Wb != nullptr &&
             !all_finite(state->hysteresis_flux_linkage_Wb, n))) {
            h->set_error("motor ROM state and inputs must contain only finite values");
            return RAD_MOTOR_ROM_INVALID_ARGUMENT;
        }
        std::vector<double> q0(state->generalized_currents_A,
                               state->generalized_currents_A + n);
        std::vector<double> q1 = q0, qnext(n), qmid(n), voltage(n, 0.0);
        for (size_t i = 0; i < h->np; ++i) voltage[i] = input->phase_voltages_V[i];
        std::vector<double> L0(n * n), L1(n * n), Lmid(n * n), dL(n * n), R(n * n);
        std::vector<double> pm0(n), pm1(n), dpm(n), motion(n), lambda0(n), work(n),
            flux_linkage(n), rhs(n);
        std::vector<double> hflux0(n, 0.0), hflux1(n, 0.0);
        if (state->hysteresis_flux_linkage_Wb != nullptr)
            std::copy(state->hysteresis_flux_linkage_Wb,
                      state->hysteresis_flux_linkage_Wb + n, hflux0.begin());

        h->eval_inductance(theta0, 0, L0);
        h->eval_pm(theta0, 0, pm0);
        matvec(L0.data(), q0.data(), n, lambda0.data());
        for (size_t i = 0; i < n; ++i) lambda0[i] += pm0[i] + hflux0[i];

        double theta1 = theta0 + dt * omega0;
        double omega1 = omega0;
        double temperature1 = temperature0;
        double hysteresis_torque = 0.0;
        double hysteresis_energy1 = state->hysteresis_stored_energy_J;
        double hysteresis_dissipation = 0.0;
        unsigned iterations = 0;
        for (iterations = 1; iterations <= max_iterations; ++iterations) {
            const double theta_mid = 0.5 * (theta0 + theta1);
            const double omega_mid = 0.5 * (omega0 + omega1);
            const double temperature_mid = 0.5 * (temperature0 + temperature1);
            h->eval_resistance(theta_mid, temperature_mid, R);
            h->eval_motion(theta_mid, motion);
            hflux1 = hflux0;
            hysteresis_torque = 0.0;
            hysteresis_energy1 = state->hysteresis_stored_energy_J;
            hysteresis_dissipation = 0.0;
            if (h->hysteresis_trial != nullptr) {
                const int status = h->hysteresis_trial(
                    h->hysteresis_user_data, theta1, q1.data(), n,
                    hflux1.data(), &hysteresis_torque, &hysteresis_energy1,
                    &hysteresis_dissipation);
                if (status != 0 || !all_finite(hflux1.data(), n) ||
                    !finite(hysteresis_torque) || !finite(hysteresis_energy1) ||
                    !finite(hysteresis_dissipation) || hysteresis_dissipation < 0.0) {
                    h->set_error("hysteresis trial callback failed or returned invalid data");
                    return RAD_MOTOR_ROM_HYSTERESIS_ERROR;
                }
            }
            h->eval_pm(theta1, 0, pm1);
            for (size_t i = 0; i < n; ++i) {
                double rq0 = 0.0;
                for (size_t j = 0; j < n; ++j) rq0 += R[i * n + j] * q0[j];
                rhs[i] = lambda0[i] + dt * voltage[i] - 0.5 * dt * rq0
                    - dt * omega_mid * motion[i] - pm1[i] - hflux1[i];
            }
            h->eval_inductance(theta1, 0, L1);
            std::vector<double> system = L1;
            for (size_t i = 0; i < n * n; ++i) system[i] += 0.5 * dt * R[i];
            if (!solve_dense(std::move(system), rhs, n, qnext)) {
                h->set_error("implicit current system is singular");
                return RAD_MOTOR_ROM_SINGULAR_SYSTEM;
            }
            for (size_t i = 0; i < n; ++i) qmid[i] = 0.5 * (q0[i] + qnext[i]);
            h->eval_inductance(theta_mid, 1, dL);
            h->eval_pm(theta_mid, 1, dpm);
            matvec(dL.data(), qmid.data(), n, work.data());
            const double torque_reluctance = 0.5 * dot(qmid.data(), work.data(), n);
            const double torque_pm = dot(qmid.data(), dpm.data(), n);
            const double torque_cogging = h->eval_cogging(theta_mid, 1);
            const double torque_motion = dot(qmid.data(), motion.data(), n);
            const double torque = torque_reluctance + torque_pm + torque_cogging
                + torque_motion + hysteresis_torque;
            const double omega_next = omega0 + dt * (
                torque - input->load_torque_Nm - h->friction * omega_mid) / h->inertia;
            const double theta_next = theta0 + 0.5 * dt * (omega0 + omega_next);
            matvec(R.data(), qmid.data(), n, work.data());
            const double loss = dot(qmid.data(), work.data(), n);
            double temperature_next = temperature0;
            if (h->thermal_capacity > 0.0) {
                temperature_next = (
                    temperature0 + dt * (loss + h->thermal_conductance * ambient)
                    / h->thermal_capacity
                ) / (1.0 + dt * h->thermal_conductance / h->thermal_capacity);
            }
            double error = std::abs(theta_next - theta1);
            error = std::max(error, std::abs(omega_next - omega1));
            error = std::max(error, std::abs(temperature_next - temperature1));
            double qnorm = 0.0;
            for (size_t i = 0; i < n; ++i) {
                error = std::max(error, std::abs(qnext[i] - q1[i]));
                qnorm += qnext[i] * qnext[i];
            }
            q1 = qnext;
            theta1 = theta_next;
            omega1 = omega_next;
            temperature1 = temperature_next;
            if (error <= tolerance * std::max({1.0, std::sqrt(qnorm), std::abs(omega1)})) break;
        }
        if (iterations > max_iterations) {
            h->set_error("motor ROM fixed-point iteration did not converge");
            return RAD_MOTOR_ROM_NONCONVERGENCE;
        }

        if (h->hysteresis_trial != nullptr) {
            const int status = h->hysteresis_trial(
                h->hysteresis_user_data, theta1, q1.data(), n,
                hflux1.data(), &hysteresis_torque, &hysteresis_energy1,
                &hysteresis_dissipation);
            if (status != 0 || !all_finite(hflux1.data(), n) ||
                !finite(hysteresis_torque) || !finite(hysteresis_energy1) ||
                !finite(hysteresis_dissipation) || hysteresis_dissipation < 0.0) {
                h->set_error("final hysteresis trial callback failed or returned invalid data");
                return RAD_MOTOR_ROM_HYSTERESIS_ERROR;
            }
        }

        const double theta_mid = 0.5 * (theta0 + theta1);
        const double omega_mid = 0.5 * (omega0 + omega1);
        for (size_t i = 0; i < n; ++i) qmid[i] = 0.5 * (q0[i] + q1[i]);
        h->eval_resistance(theta_mid, 0.5 * (temperature0 + temperature1), R);
        h->eval_inductance(theta_mid, 1, dL);
        h->eval_pm(theta_mid, 1, dpm);
        h->eval_motion(theta_mid, motion);
        matvec(dL.data(), qmid.data(), n, work.data());
        const double reluctance = 0.5 * dot(qmid.data(), work.data(), n);
        const double permanent_magnet = dot(qmid.data(), dpm.data(), n);
        const double cogging = h->eval_cogging(theta_mid, 1);
        const double motional = dot(qmid.data(), motion.data(), n);
        const double torque = reluctance + permanent_magnet + cogging + motional
            + hysteresis_torque;
        matvec(R.data(), qmid.data(), n, work.data());
        const double resistive_loss = dot(qmid.data(), work.data(), n);
        const double electrical_power = dot(qmid.data(), voltage.data(), n);
        const double hysteresis_loss = hysteresis_dissipation / dt;

        h->eval_inductance(theta1, 0, L1);
        h->eval_pm(theta1, 0, pm1);
        matvec(L1.data(), q1.data(), n, work.data());
        const double magnetic0 = 0.5 * dot(q0.data(), lambda0.data(), n)
            - 0.5 * dot(q0.data(), pm0.data(), n)
            - 0.5 * dot(q0.data(), hflux0.data(), n);
        const double magnetic1 = 0.5 * dot(q1.data(), work.data(), n);
        const double stored0 = magnetic0 + 0.5 * h->inertia * omega0 * omega0
            - h->eval_cogging(theta0, 0)
            + state->hysteresis_stored_energy_J;
        const double stored1 = magnetic1 + 0.5 * h->inertia * omega1 * omega1
            - h->eval_cogging(theta1, 0)
            + hysteresis_energy1;
        const double load_power = input->load_torque_Nm * omega_mid;
        const double balance = electrical_power - resistive_loss - hysteresis_loss
            - load_power - h->friction * omega_mid * omega_mid
            - (stored1 - stored0) / dt;

        if (output->speed_voltage_V != nullptr) {
            for (size_t i = 0; i < n; ++i) {
                double value = dpm[i] + motion[i];
                for (size_t j = 0; j < n; ++j) value += dL[i * n + j] * qmid[j];
                output->speed_voltage_V[i] = omega_mid * value;
            }
        }
        if (output->phase_flux_linkage_Wb != nullptr) {
            matvec(L1.data(), q1.data(), n, flux_linkage.data());
            for (size_t i = 0; i < h->np; ++i)
                output->phase_flux_linkage_Wb[i] =
                    flux_linkage[i] + pm1[i] + hflux1[i];
        }
        std::copy(q1.begin(), q1.end(), state->generalized_currents_A);
        if (state->hysteresis_flux_linkage_Wb != nullptr)
            std::copy(hflux1.begin(), hflux1.end(), state->hysteresis_flux_linkage_Wb);
        state->time_s += dt;
        state->rotor_angle_rad = std::fmod(theta1, h->period);
        if (state->rotor_angle_rad < 0.0) state->rotor_angle_rad += h->period;
        state->rotor_speed_rad_s = omega1;
        state->temperature_K = temperature1;
        state->hysteresis_stored_energy_J = hysteresis_energy1;

        output->electromagnetic_torque_Nm = torque;
        output->reluctance_torque_Nm = reluctance;
        output->permanent_magnet_torque_Nm = permanent_magnet;
        output->cogging_torque_Nm = cogging;
        output->motional_lorentz_torque_Nm = motional;
        output->hysteresis_torque_Nm = hysteresis_torque;
        output->resistive_loss_W = resistive_loss;
        output->hysteresis_loss_W = hysteresis_loss;
        output->electrical_input_power_W = electrical_power;
        output->mechanical_load_power_W = load_power;
        output->stored_energy_J = stored1;
        output->energy_balance_residual_W = balance;
        output->nonlinear_iterations = iterations;
        h->last_error.clear();
        return RAD_MOTOR_ROM_OK;
    } catch (const std::exception& error) {
        h->set_error(error.what());
    } catch (...) {
        h->set_error("unknown motor ROM step error");
    }
    return RAD_MOTOR_ROM_INTERNAL_ERROR;
}

}  // extern "C"
