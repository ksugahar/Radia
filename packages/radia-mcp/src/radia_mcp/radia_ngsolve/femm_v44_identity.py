from __future__ import annotations

import math

from .femm_v45_identity import validate_public_v45_identity


def _sha(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, dict):
        return {}
    checks: dict[str, bool] = {}
    actuator = identity.get(
        "v44_axisymmetric_actuator_flux_force_coenergy_gap_derivative_axisfactor_mesh_mismatch"
    )
    if isinstance(actuator, dict):
        generation = str(actuator.get("generation", "")).strip()
        try:
            flux = float(actuator["flux_linkage_wb_turn"])
            result_flux = float(actuator["result_flux_linkage_wb_turn"])
            force = float(actuator["force_n"])
            result_force = float(actuator["result_force_n"])
            coenergy = float(actuator["coenergy_j"])
            result_coenergy = float(actuator["result_coenergy_j"])
            gap = float(actuator["gap_m"])
            result_gap = float(actuator["result_gap_m"])
            derivative = float(actuator["coenergy_force_derivative_n"])
            result_derivative = float(actuator["result_coenergy_force_derivative_n"])
            factor = float(actuator["axisymmetric_factor"])
            result_factor = float(actuator["result_axisymmetric_factor"])
        except (KeyError, TypeError, ValueError):
            return {"femm_v44_axisymmetric_actuator_identity": False}
        mirrored = all(
            actuator.get(field) == generation
            for field in ("flux_generation", "force_generation", "coenergy_generation", "gap_generation", "derivative_generation", "axis_factor_generation", "mesh_generation", "result_generation")
        )
        checks.update(
            {
                "femm_v44_axisymmetric_actuator_generation_closure": bool(generation) and mirrored,
                "femm_v44_axisymmetric_actuator_flux_force_coenergy_mirror": (
                    flux > 0.0 and result_flux == flux and coenergy > 0.0 and result_coenergy == coenergy
                    and math.isfinite(force) and result_force == force
                    and math.isclose(force, derivative, rel_tol=1.0e-12, abs_tol=1.0e-12)
                    and result_derivative == derivative
                ),
                "femm_v44_axisymmetric_actuator_gap_and_axis_factor": (
                    gap > 0.0 and result_gap == gap
                    and math.isclose(factor, 2.0 * math.pi, rel_tol=1.0e-12, abs_tol=1.0e-15)
                    and result_factor == factor
                ),
                "femm_v44_axisymmetric_actuator_mesh_and_result_owner": (
                    str(actuator.get("mesh_owner", "")).startswith("mesh:")
                    and actuator.get("result_mesh_owner") == actuator.get("mesh_owner")
                ),
                "femm_v44_axisymmetric_actuator_result_digest": (
                    _sha(actuator.get("result_sha256"))
                    and actuator.get("accepted_result_sha256") == actuator.get("result_sha256")
                ),
            }
        )
    fringe = identity.get(
        "v44_electrostatic_fringe_capacitance_charge_energy_interface_axisfactor_mesh_mismatch"
    )
    if isinstance(fringe, dict):
        generation = str(fringe.get("generation", "")).strip()
        try:
            voltage = float(fringe["voltage_v"])
            result_voltage = float(fringe["result_voltage_v"])
            capacitance = float(fringe["capacitance_f"])
            result_capacitance = float(fringe["result_capacitance_f"])
            charge = float(fringe["charge_c"])
            result_charge = float(fringe["result_charge_c"])
            energy = float(fringe["stored_energy_j"])
            result_energy = float(fringe["result_stored_energy_j"])
            flux = float(fringe["interface_flux_c"])
            result_flux = float(fringe["result_interface_flux_c"])
            factor = float(fringe["axisymmetric_factor"])
            result_factor = float(fringe["result_axisymmetric_factor"])
        except (KeyError, TypeError, ValueError):
            return {"femm_v44_electrostatic_fringe_identity": False}
        mirrored = all(
            fringe.get(field) == generation
            for field in ("capacitance_generation", "charge_generation", "energy_generation", "interface_generation", "axis_factor_generation", "mesh_generation", "result_generation")
        )
        checks.update(
            {
                "femm_v44_electrostatic_fringe_generation_closure": bool(generation) and mirrored,
                "femm_v44_electrostatic_fringe_charge_energy_closure": (
                    voltage > 0.0 and result_voltage == voltage and capacitance > 0.0
                    and result_capacitance == capacitance
                    and math.isclose(charge, capacitance * voltage, rel_tol=1.0e-12, abs_tol=1.0e-15)
                    and result_charge == charge
                    and math.isclose(energy, 0.5 * capacitance * voltage * voltage, rel_tol=1.0e-12, abs_tol=1.0e-15)
                    and result_energy == energy and result_flux == charge
                ),
                "femm_v44_electrostatic_fringe_axis_factor": (
                    math.isclose(factor, 2.0 * math.pi, rel_tol=1.0e-12, abs_tol=1.0e-15)
                    and result_factor == factor
                ),
                "femm_v44_electrostatic_fringe_mesh_and_result_owner": (
                    str(fringe.get("mesh_owner", "")).startswith("mesh:")
                    and fringe.get("result_mesh_owner") == fringe.get("mesh_owner")
                ),
                "femm_v44_electrostatic_fringe_result_digest": (
                    _sha(fringe.get("result_sha256"))
                    and fringe.get("accepted_result_sha256") == fringe.get("result_sha256")
                ),
            }
        )
    checks.update(validate_public_v45_identity(identity))
    return checks
