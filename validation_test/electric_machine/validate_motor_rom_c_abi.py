"""Cross-language lock for the angle-periodic motor ROM C ABI.

Build ``src/core/rad_motor_rom_c.cpp`` as a shared library, then run this file
with ``--library``.  It advances the Python and C implementations through the
same 1000 implicit electromechanical steps and writes the maximum discrepancy.
"""
from __future__ import annotations

import argparse
import ctypes as ct
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import sys

import numpy as np

from radia.motor_rom import (
    AnglePeriodicMotorROM,
    MotorPortContract,
    MotorROMInput,
    PeriodicAngleTable,
)


ABI = 1
HERE = Path(__file__).resolve().parent
DEFAULT_RESULT = HERE / "motor_rom_c_abi_summary.json"
DoublePointer = ct.POINTER(ct.c_double)


class ModelData(ct.Structure):
    _fields_ = [
        ("abi_version", ct.c_uint32),
        ("n_angle_samples", ct.c_size_t),
        ("n_phase", ct.c_size_t),
        ("n_generalized", ct.c_size_t),
        ("angle_origin_rad", ct.c_double),
        ("period_rad", ct.c_double),
        ("skew_span_rad", ct.c_double),
        ("inertia_kg_m2", ct.c_double),
        ("viscous_friction_Nm_s", ct.c_double),
        ("reference_temperature_K", ct.c_double),
        ("thermal_capacity_J_per_K", ct.c_double),
        ("thermal_conductance_W_per_K", ct.c_double),
        ("inductance_H", DoublePointer),
        ("resistance_ohm", DoublePointer),
        ("pm_flux_linkage_Wb", DoublePointer),
        ("cogging_coenergy_J", DoublePointer),
        ("motion_flux_gradient_Wb_per_rad", DoublePointer),
        ("end_winding_inductance_H", DoublePointer),
        ("end_winding_resistance_ohm", DoublePointer),
        ("resistance_temperature_coefficient_per_K", DoublePointer),
        ("hysteresis_trial", ct.c_void_p),
        ("hysteresis_user_data", ct.c_void_p),
    ]


class State(ct.Structure):
    _fields_ = [
        ("abi_version", ct.c_uint32),
        ("time_s", ct.c_double),
        ("rotor_angle_rad", ct.c_double),
        ("rotor_speed_rad_s", ct.c_double),
        ("temperature_K", ct.c_double),
        ("generalized_currents_A", DoublePointer),
        ("hysteresis_flux_linkage_Wb", DoublePointer),
        ("hysteresis_stored_energy_J", ct.c_double),
    ]


class Input(ct.Structure):
    _fields_ = [
        ("abi_version", ct.c_uint32),
        ("phase_voltages_V", DoublePointer),
        ("load_torque_Nm", ct.c_double),
        ("has_ambient_temperature", ct.c_int),
        ("ambient_temperature_K", ct.c_double),
    ]


class Output(ct.Structure):
    _fields_ = [
        ("abi_version", ct.c_uint32),
        ("phase_flux_linkage_Wb", DoublePointer),
        ("speed_voltage_V", DoublePointer),
        ("electromagnetic_torque_Nm", ct.c_double),
        ("reluctance_torque_Nm", ct.c_double),
        ("permanent_magnet_torque_Nm", ct.c_double),
        ("cogging_torque_Nm", ct.c_double),
        ("motional_lorentz_torque_Nm", ct.c_double),
        ("hysteresis_torque_Nm", ct.c_double),
        ("resistive_loss_W", ct.c_double),
        ("hysteresis_loss_W", ct.c_double),
        ("electrical_input_power_W", ct.c_double),
        ("mechanical_load_power_W", ct.c_double),
        ("stored_energy_J", ct.c_double),
        ("energy_balance_residual_W", ct.c_double),
        ("nonlinear_iterations", ct.c_uint),
    ]


def pointer(array: np.ndarray) -> DoublePointer:
    return array.ctypes.data_as(DoublePointer)


def build_model():
    angles = np.linspace(0.0, 2.0 * np.pi, 33, endpoint=False)
    inductance = np.ascontiguousarray(
        [
            [
                [0.02 + 0.004 * np.cos(2.0 * angle), 0.001 * np.sin(angle)],
                [0.001 * np.sin(angle), 0.006],
            ]
            for angle in angles
        ],
        dtype=np.float64,
    )
    resistance = np.ascontiguousarray(
        np.repeat(np.diag([0.4, 2.0])[None, :, :], angles.size, axis=0),
        dtype=np.float64,
    )
    pm_flux = np.ascontiguousarray(
        [[0.03 * np.cos(angle), 0.01 * np.sin(angle)] for angle in angles],
        dtype=np.float64,
    )
    motion = np.ascontiguousarray(
        [[0.0, 0.003 * np.cos(angle)] for angle in angles], dtype=np.float64
    )
    cogging = np.ascontiguousarray(
        0.002 * np.cos(6.0 * angles), dtype=np.float64
    )
    alpha = np.ascontiguousarray([0.0039, 0.0039], dtype=np.float64)
    motor = AnglePeriodicMotorROM(
        MotorPortContract(("A",), ("eddy-bulk",)),
        PeriodicAngleTable(angles, inductance),
        PeriodicAngleTable(angles, resistance),
        PeriodicAngleTable(angles, pm_flux),
        inertia_kg_m2=0.02,
        viscous_friction_Nm_s=1.0e-3,
        motion_flux_gradient_Wb_per_rad=PeriodicAngleTable(angles, motion),
        cogging_coenergy_J=PeriodicAngleTable(angles, cogging),
        resistance_temperature_coefficient_per_K=alpha,
    )
    arrays = (angles, inductance, resistance, pm_flux, cogging, motion, alpha)
    model_data = ModelData(
        ABI,
        angles.size,
        1,
        2,
        angles[0],
        2.0 * np.pi,
        0.0,
        0.02,
        1.0e-3,
        293.15,
        0.0,
        0.0,
        pointer(inductance),
        pointer(resistance),
        pointer(pm_flux),
        pointer(cogging),
        pointer(motion),
        None,
        None,
        pointer(alpha),
        None,
        None,
    )
    return motor, model_data, arrays


def validate(library_path: Path) -> dict[str, object]:
    library = ct.CDLL(str(library_path))
    library.rad_motor_rom_abi_version.restype = ct.c_uint32
    library.rad_motor_rom_create.argtypes = [ct.POINTER(ModelData)]
    library.rad_motor_rom_create.restype = ct.c_void_p
    library.rad_motor_rom_destroy.argtypes = [ct.c_void_p]
    library.rad_motor_rom_step.argtypes = [
        ct.c_void_p,
        ct.POINTER(State),
        ct.POINTER(Input),
        ct.c_double,
        ct.c_uint,
        ct.c_double,
        ct.POINTER(Output),
    ]
    library.rad_motor_rom_step.restype = ct.c_int
    library.rad_motor_rom_last_error.argtypes = [ct.c_void_p]
    library.rad_motor_rom_last_error.restype = ct.c_char_p
    assert library.rad_motor_rom_abi_version() == ABI

    motor, model_data, arrays = build_model()
    handle = library.rad_motor_rom_create(ct.byref(model_data))
    if not handle:
        raise RuntimeError(library.rad_motor_rom_last_error(None).decode())

    currents = np.ascontiguousarray([1.0, 0.0], dtype=np.float64)
    hflux = np.ascontiguousarray([0.0, 0.0], dtype=np.float64)
    voltages = np.ascontiguousarray([2.0], dtype=np.float64)
    phase_flux = np.zeros(1, dtype=np.float64)
    speed_voltage = np.zeros(2, dtype=np.float64)
    c_state = State(ABI, 0.0, 0.2, 40.0, 293.15, pointer(currents), pointer(hflux), 0.0)
    c_input = Input(ABI, pointer(voltages), 0.1, 0, 293.15)
    c_output = Output(ABI, pointer(phase_flux), pointer(speed_voltage))
    py_state = motor.initial_state(
        rotor_angle_rad=0.2,
        rotor_speed_rad_s=40.0,
        phase_currents_A=(1.0,),
    )
    py_input = MotorROMInput(voltages, load_torque_Nm=0.1)

    maxima = {
        "rotor_angle_rad": 0.0,
        "rotor_speed_rad_s": 0.0,
        "generalized_currents_A": 0.0,
        "electromagnetic_torque_Nm": 0.0,
        "resistive_loss_W": 0.0,
        "energy_balance_residual_W_difference": 0.0,
        "python_energy_balance_residual_W": 0.0,
        "c_energy_balance_residual_W": 0.0,
    }
    nonpositive_scale_status = None
    nonpositive_scale_error = ""
    nonfinite_state_status = None
    nonfinite_state_error = ""
    dt = 2.0e-6
    try:
        for _ in range(1000):
            py_state, py_output = motor.step(py_state, py_input, dt)
            status = library.rad_motor_rom_step(
                handle, ct.byref(c_state), ct.byref(c_input), dt, 30, 1.0e-11,
                ct.byref(c_output),
            )
            if status != 0:
                raise RuntimeError(library.rad_motor_rom_last_error(handle).decode())
            maxima["rotor_angle_rad"] = max(
                maxima["rotor_angle_rad"], abs(c_state.rotor_angle_rad - py_state.rotor_angle_rad)
            )
            maxima["rotor_speed_rad_s"] = max(
                maxima["rotor_speed_rad_s"], abs(c_state.rotor_speed_rad_s - py_state.rotor_speed_rad_s)
            )
            maxima["generalized_currents_A"] = max(
                maxima["generalized_currents_A"],
                float(np.max(np.abs(currents - py_state.generalized_currents_A))),
            )
            maxima["electromagnetic_torque_Nm"] = max(
                maxima["electromagnetic_torque_Nm"],
                abs(c_output.electromagnetic_torque_Nm - py_output.electromagnetic_torque_Nm),
            )
            maxima["resistive_loss_W"] = max(
                maxima["resistive_loss_W"], abs(c_output.resistive_loss_W - py_output.resistive_loss_W)
            )
            maxima["energy_balance_residual_W_difference"] = max(
                maxima["energy_balance_residual_W_difference"],
                abs(c_output.energy_balance_residual_W - py_output.energy_balance_residual_W),
            )
            maxima["python_energy_balance_residual_W"] = max(
                maxima["python_energy_balance_residual_W"], abs(py_output.energy_balance_residual_W)
            )
            maxima["c_energy_balance_residual_W"] = max(
                maxima["c_energy_balance_residual_W"], abs(c_output.energy_balance_residual_W)
            )
        c_state.temperature_K = 0.0
        nonpositive_scale_status = library.rad_motor_rom_step(
            handle, ct.byref(c_state), ct.byref(c_input), dt, 30, 1.0e-11,
            ct.byref(c_output),
        )
        nonpositive_scale_error = library.rad_motor_rom_last_error(handle).decode()
        c_state.temperature_K = 293.15
        c_state.rotor_angle_rad = float("nan")
        nonfinite_state_status = library.rad_motor_rom_step(
            handle, ct.byref(c_state), ct.byref(c_input), dt, 30, 1.0e-11,
            ct.byref(c_output),
        )
        nonfinite_state_error = library.rad_motor_rom_last_error(handle).decode()
    finally:
        library.rad_motor_rom_destroy(handle)

    tolerances = {
        "state_and_output": 2.0e-10,
        "energy_balance_residual_W": 2.0e-7,
    }
    passed = (
        max(
            maxima["rotor_angle_rad"],
            maxima["rotor_speed_rad_s"],
            maxima["generalized_currents_A"],
            maxima["electromagnetic_torque_Nm"],
            maxima["resistive_loss_W"],
        ) < tolerances["state_and_output"]
        and maxima["python_energy_balance_residual_W"] < tolerances["energy_balance_residual_W"]
        and maxima["c_energy_balance_residual_W"] < tolerances["energy_balance_residual_W"]
        and nonpositive_scale_status == 6
        and "must remain positive" in nonpositive_scale_error
        and nonfinite_state_status == 1
        and "only finite values" in nonfinite_state_error
    )
    return {
        "schema": "radia.motor.c_abi_validation.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "library_name": library_path.name,
        "abi_version": ABI,
        "steps": 1000,
        "dt_s": dt,
        "max_abs_errors": maxima,
        "fail_loud_gate": {
            "nonpositive_temperature_scale_status": nonpositive_scale_status,
            "expected_status": 6,
            "error": nonpositive_scale_error,
            "nonfinite_state_status": nonfinite_state_status,
            "expected_nonfinite_state_status": 1,
            "nonfinite_state_error": nonfinite_state_error,
        },
        "tolerances": tolerances,
        "passed": passed,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=DEFAULT_RESULT)
    args = parser.parse_args()
    result = validate(args.library)
    args.out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
