#ifndef RAD_MOTOR_ROM_C_H
#define RAD_MOTOR_ROM_C_H

#include <stddef.h>
#include <stdint.h>

#if defined(_WIN32)
#  if defined(RAD_MOTOR_ROM_BUILD)
#    define RAD_MOTOR_ROM_API __declspec(dllexport)
#  else
#    define RAD_MOTOR_ROM_API __declspec(dllimport)
#  endif
#else
#  define RAD_MOTOR_ROM_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define RAD_MOTOR_ROM_ABI_VERSION 1u

enum RadMotorROMStatus {
    RAD_MOTOR_ROM_OK = 0,
    RAD_MOTOR_ROM_INVALID_ARGUMENT = 1,
    RAD_MOTOR_ROM_INVALID_MODEL = 2,
    RAD_MOTOR_ROM_SINGULAR_SYSTEM = 3,
    RAD_MOTOR_ROM_NONCONVERGENCE = 4,
    RAD_MOTOR_ROM_HYSTERESIS_ERROR = 5,
    RAD_MOTOR_ROM_INTERNAL_ERROR = 6
};

typedef int (*RadMotorROMHysteresisTrial)(
    void* user_data,
    double rotor_angle_rad,
    const double* generalized_currents_A,
    size_t n_generalized,
    double* flux_linkage_Wb,
    double* torque_Nm,
    double* stored_energy_J,
    double* dissipated_energy_increment_J);

/* Every table is C-contiguous row-major with angle as the first index. */
typedef struct RadMotorROMModelData {
    uint32_t abi_version;
    size_t n_angle_samples;
    size_t n_phase;
    size_t n_generalized;
    double angle_origin_rad;
    double period_rad;
    double skew_span_rad;
    double inertia_kg_m2;
    double viscous_friction_Nm_s;
    double reference_temperature_K;
    double thermal_capacity_J_per_K;
    double thermal_conductance_W_per_K;
    const double* inductance_H;                 /* [na, n, n] */
    const double* resistance_ohm;               /* [na, n, n] */
    const double* pm_flux_linkage_Wb;            /* [na, n] */
    const double* cogging_coenergy_J;             /* [na], optional */
    const double* motion_flux_gradient_Wb_per_rad; /* [na, n], optional */
    const double* end_winding_inductance_H;      /* [n, n], optional */
    const double* end_winding_resistance_ohm;    /* [n, n], optional */
    const double* resistance_temperature_coefficient_per_K; /* [n], optional */
    RadMotorROMHysteresisTrial hysteresis_trial; /* optional, pure trial */
    void* hysteresis_user_data;
} RadMotorROMModelData;

typedef struct RadMotorROMState {
    uint32_t abi_version;
    double time_s;
    double rotor_angle_rad;
    double rotor_speed_rad_s;
    double temperature_K;
    double* generalized_currents_A;             /* [n], caller-owned */
    double* hysteresis_flux_linkage_Wb;          /* [n], optional */
    double hysteresis_stored_energy_J;
} RadMotorROMState;

typedef struct RadMotorROMInput {
    uint32_t abi_version;
    const double* phase_voltages_V;              /* [n_phase] */
    double load_torque_Nm;
    int has_ambient_temperature;
    double ambient_temperature_K;
} RadMotorROMInput;

typedef struct RadMotorROMStepOutput {
    uint32_t abi_version;
    double* phase_flux_linkage_Wb;                /* [n_phase], optional */
    double* speed_voltage_V;                      /* [n], optional */
    double electromagnetic_torque_Nm;
    double reluctance_torque_Nm;
    double permanent_magnet_torque_Nm;
    double cogging_torque_Nm;
    double motional_lorentz_torque_Nm;
    double hysteresis_torque_Nm;
    double resistive_loss_W;
    double hysteresis_loss_W;
    double electrical_input_power_W;
    double mechanical_load_power_W;
    double stored_energy_J;
    double energy_balance_residual_W;
    unsigned nonlinear_iterations;
} RadMotorROMStepOutput;

typedef struct RadMotorROMHandle RadMotorROMHandle;

RAD_MOTOR_ROM_API RadMotorROMHandle* rad_motor_rom_create(
    const RadMotorROMModelData* model);

RAD_MOTOR_ROM_API void rad_motor_rom_destroy(RadMotorROMHandle* handle);

RAD_MOTOR_ROM_API int rad_motor_rom_step(
    RadMotorROMHandle* handle,
    RadMotorROMState* state,
    const RadMotorROMInput* input,
    double dt_s,
    unsigned max_iterations,
    double tolerance,
    RadMotorROMStepOutput* output);

RAD_MOTOR_ROM_API const char* rad_motor_rom_last_error(
    const RadMotorROMHandle* handle);

RAD_MOTOR_ROM_API uint32_t rad_motor_rom_abi_version(void);

#ifdef __cplusplus
}
#endif

#endif
