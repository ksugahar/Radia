"""Portable export contracts for angle-periodic motor ROMs.

The bundle is intentionally solver-neutral: NPZ is the canonical numerical
payload, MAT is the Simulink loading surface, JSON fixes units and port order,
and an FMI 3.0.2 ``ModelVariables`` fragment maps the same scalar ports for an
FMU wrapper.  The actual stepping ABI is ``src/core/rad_motor_rom_c.h``.
"""
from __future__ import annotations

import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np

from .motor_rom import AnglePeriodicMotorROM, LoadAnglePeriodicMotorROM


def _port_variables(motor: AnglePeriodicMotorROM) -> list[dict[str, object]]:
    variables: list[dict[str, object]] = []

    def add(name, causality, unit, description, variability="continuous"):
        variables.append(
            {
                "name": name,
                "value_reference": len(variables),
                "causality": causality,
                "variability": variability,
                "unit": unit,
                "description": description,
            }
        )

    for name in motor.ports.phase_names:
        add(f"phase_voltage.{name}", "input", "V", f"Voltage of phase {name}")
    add("load_torque", "input", "N.m", "Mechanical load torque")
    add("ambient_temperature", "input", "K", "Ambient temperature")
    for name in motor.ports.phase_names:
        add(f"phase_current.{name}", "output", "A", f"Current of phase {name}")
        add(f"phase_flux_linkage.{name}", "output", "Wb", f"Flux linkage of phase {name}")
    for name in motor.ports.eddy_names:
        add(f"eddy_current.{name}", "output", "A", f"Internal eddy coordinate {name}")
    add("rotor_angle", "output", "rad", "Mechanical rotor angle")
    add("rotor_speed", "output", "rad/s", "Mechanical rotor angular speed")
    add("electromagnetic_torque", "output", "N.m", "Total electromagnetic torque")
    add("resistive_loss", "output", "W", "Ohmic loss")
    add("hysteresis_loss", "output", "W", "Hysteresis loss")
    add("temperature", "output", "K", "Lumped conductor temperature")
    add("energy_balance_residual", "output", "W", "Discrete power-balance residual")
    return variables


def MotorROMPortManifest(motor: AnglePeriodicMotorROM) -> dict[str, object]:
    """Return the versioned C/Simulink/FMI port and array contract."""

    return {
        "schema": "radia.motor.port_contract.v1",
        "model_schema": "radia.motor.angle_periodic_rom.v1",
        "array_order": "C-row-major",
        "generalized_current_order": list(motor.ports.generalized_names),
        "phase_names": list(motor.ports.phase_names),
        "eddy_names": list(motor.ports.eddy_names),
        "angle_interpolation": "odd-uniform-periodic-Fourier",
        "c_abi": {
            "abi_version": 1,
            "header": "src/core/rad_motor_rom_c.h",
            "library_base_name": "radia_motor_rom",
            "hysteresis_contract": "pure-trial-callback; commit only after accepted step",
        },
        "simulink": {
            "payload": "MAT mirrors canonical NPZ arrays",
            "step_function": "rad_motor_rom_step",
            "mex_s_function": "radia_motor_rom_sfun",
            "input_order": [
                "phase_voltages_V",
                "load_torque_Nm",
                "ambient_temperature_K",
            ],
            "output_order": [
                "phase_currents_A",
                "eddy_currents_A",
                "phase_flux_linkage_Wb",
                "rotor_angle_rad",
                "rotor_speed_rad_s",
                "electromagnetic_torque_Nm",
                "resistive_loss_W",
                "hysteresis_loss_W",
                "temperature_K",
                "energy_balance_residual_W",
                "nonlinear_iterations",
            ],
            "state_update": "internal-discrete-at-fixed-sample-time",
            "sample_time": "inherited/fixed communication step",
        },
        "fmi": {
            "version": "3.0.2",
            "interface": "Co-Simulation source boundary",
            "packaged_fmu": False,
            "note": "Compile the C ABI behind an FMI 3 wrapper; this bundle is not itself an FMU.",
        },
        "variables": _port_variables(motor),
        "features": motor.diagnostics(),
    }


def _write_fmi_variables_fragment(path: Path, manifest: dict[str, object]) -> None:
    root = ET.Element("ModelVariables")
    for variable in manifest["variables"]:
        attributes = {
            "name": str(variable["name"]),
            "valueReference": str(variable["value_reference"]),
            "causality": str(variable["causality"]),
            "variability": str(variable["variability"]),
            "unit": str(variable["unit"]),
            "description": str(variable["description"]),
        }
        ET.SubElement(root, "Float64", attributes)
    ET.indent(root, space="  ")
    path.write_text(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        + ET.tostring(root, encoding="unicode")
        + "\n",
        encoding="utf-8",
    )


def _write_matlab_loader(path: Path, mat_name: str, manifest_name: str) -> None:
    path.write_text(
        f"""% Load a Radia angle-periodic motor ROM for Simulink.
bundleDir = fileparts(mfilename('fullpath'));
RadiaMotorROM = load(fullfile(bundleDir, '{mat_name}'));
RadiaMotorROM.port_contract = jsondecode(fileread(fullfile(bundleDir, '{manifest_name}')));

assert(strcmp(RadiaMotorROM.schema, 'radia.motor.angle_periodic_rom.v1'));
assert(mod(numel(RadiaMotorROM.angles_rad), 2) == 1);
assert(size(RadiaMotorROM.inductance_H, 1) == numel(RadiaMotorROM.angles_rad));
assert(size(RadiaMotorROM.inductance_H, 2) == size(RadiaMotorROM.inductance_H, 3));

assignin('base', 'RadiaMotorROM', RadiaMotorROM);
fprintf('Radia motor ROM: %d angle samples, %d phases, %d internal eddy states\\n', ...
    numel(RadiaMotorROM.angles_rad), numel(RadiaMotorROM.phase_names), ...
    numel(RadiaMotorROM.eddy_names));
fprintf('C ABI: radia_motor_rom, version %d\\n', RadiaMotorROM.c_abi_version);
""",
        encoding="utf-8",
    )


def SaveMotorROMBundle(motor: AnglePeriodicMotorROM, path) -> dict[str, str]:
    """Write synchronized NPZ, MAT, JSON, MATLAB, and FMI-variable artifacts."""

    from scipy.io import savemat

    base = Path(path)
    if base.suffix:
        base = base.with_suffix("")
    base.parent.mkdir(parents=True, exist_ok=True)
    npz_path = base.with_suffix(".npz")
    mat_path = base.with_suffix(".mat")
    json_path = base.with_suffix(".json")
    matlab_path = base.parent / f"{base.name}_load.m"
    fmi_path = base.parent / f"{base.name}_fmi_model_variables.xml"

    motor.save_npz(npz_path)
    with np.load(npz_path, allow_pickle=False) as payload:
        mat_payload = {name: np.asarray(payload[name]) for name in payload.files}
    # Keep scalar contract fields as doubles: this is the least surprising
    # representation for MATLAB MEX/S-function parameter inspection.
    mat_payload["n_phase"] = np.asarray(motor.ports.n_phase, dtype=np.float64)
    mat_payload["n_eddy"] = np.asarray(motor.ports.n_eddy, dtype=np.float64)
    mat_payload["n_generalized"] = np.asarray(
        motor.ports.n_generalized, dtype=np.float64
    )
    if motor.thermal_capacity_J_per_K is None:
        mat_payload["thermal_capacity_J_per_K"] = np.asarray(0.0, dtype=np.float64)
    mat_payload["has_motional_v_cross_b"] = np.asarray(
        motor.motion_flux_gradient_Wb_per_rad is not None, dtype=bool
    )
    mat_payload["has_cogging_coenergy"] = np.asarray(
        motor.cogging_coenergy_J is not None, dtype=bool
    )
    mat_payload["external_hysteresis_required"] = np.asarray(
        motor.hysteresis_port is not None, dtype=bool
    )
    mat_payload["c_abi_version"] = np.asarray(1, dtype=np.uint32)
    savemat(mat_path, mat_payload, do_compression=True, long_field_names=True)

    manifest = MotorROMPortManifest(motor)
    json_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    _write_matlab_loader(matlab_path, mat_path.name, json_path.name)
    _write_fmi_variables_fragment(fmi_path, manifest)
    return {
        "npz": str(npz_path),
        "mat": str(mat_path),
        "json": str(json_path),
        "matlab_loader": str(matlab_path),
        "fmi_model_variables": str(fmi_path),
    }


def ValidateMotorROMBundle(path) -> dict[str, object]:
    """Read back a bundle and verify shape, port, and passivity contracts."""

    from scipy.io import loadmat

    base = Path(path)
    if base.suffix:
        base = base.with_suffix("")
    motor = LoadAnglePeriodicMotorROM(base.with_suffix(".npz"))
    manifest = json.loads(base.with_suffix(".json").read_text(encoding="utf-8"))
    mat = loadmat(base.with_suffix(".mat"), squeeze_me=True)
    if manifest.get("schema") != "radia.motor.port_contract.v1":
        raise ValueError("invalid motor port manifest schema")
    if int(np.asarray(mat["c_abi_version"]).item()) != 1:
        raise ValueError("MAT payload uses an unsupported C ABI")
    if list(motor.ports.generalized_names) != manifest["generalized_current_order"]:
        raise ValueError("generalized-current order differs between NPZ and manifest")
    diagnostics = motor.diagnostics()
    if not diagnostics["passive"]:
        raise ValueError("reloaded motor ROM is not passive")
    return {
        "passed": True,
        "angle_samples": int(motor.inductance_H.angles_rad.size),
        "n_phase": motor.ports.n_phase,
        "n_eddy": motor.ports.n_eddy,
        "c_abi_version": 1,
        "passive": True,
    }


__all__ = [
    "MotorROMPortManifest",
    "SaveMotorROMBundle",
    "ValidateMotorROMBundle",
]
