from __future__ import annotations

import math
from collections.abc import Mapping


def _sha(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _closed(row: Mapping[str, object], fields: tuple[str, ...]) -> bool:
    generation = str(row.get("generation", "")).strip()
    return bool(generation) and all(row.get(field) == generation for field in fields)


def validate_public_v45_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    force = identity.get("v45_public_axisymmetric_force_torque_coenergy_stress_contour_owner_mismatch")
    if isinstance(force, Mapping):
        checks["femm_v45_axisymmetric_force_contour_generation"] = _closed(
            force, ("force_generation", "torque_generation", "coenergy_generation", "stress_contour_generation", "axis_factor_generation", "mesh_generation", "result_generation")
        )
        checks["femm_v45_axisymmetric_force_contour_values"] = (
            force.get("force_method") == "weighted_stress_tensor"
            and force.get("result_force_method") == force.get("force_method")
            and force.get("torque_method") == "airgap_contour"
            and force.get("result_torque_method") == force.get("torque_method")
            and force.get("force_n") == force.get("result_force_n")
            and math.isclose(float(force.get("torque_nm")), float(force.get("result_torque_nm")), rel_tol=1e-12)
            and math.isclose(float(force.get("coenergy_j")), float(force.get("result_coenergy_j")), rel_tol=1e-12)
        )
        checks["femm_v45_axisymmetric_force_contour_owner"] = (
            math.isclose(float(force.get("axisymmetric_factor")), 2.0 * math.pi, rel_tol=1e-12)
            and force.get("result_axisymmetric_factor") == force.get("axisymmetric_factor")
            and str(force.get("contour_owner", "")).startswith("contour:")
            and force.get("result_contour_owner") == force.get("contour_owner")
            and str(force.get("mesh_owner", "")).startswith("mesh:")
            and force.get("result_mesh_owner") == force.get("mesh_owner")
            and _sha(force.get("result_sha256"))
            and force.get("accepted_result_sha256") == force.get("result_sha256")
        )
    fringe = identity.get("v45_public_electrostatic_fringe_charge_energy_capacitance_interface_flux_axisfactor_mismatch")
    if isinstance(fringe, Mapping):
        checks["femm_v45_electrostatic_fringe_generation"] = _closed(
            fringe, ("charge_generation", "energy_generation", "capacitance_generation", "interface_generation", "axis_factor_generation", "mesh_generation", "result_generation")
        )
        voltage = float(fringe.get("voltage_v"))
        capacitance = float(fringe.get("capacitance_f"))
        charge = float(fringe.get("charge_c"))
        energy = float(fringe.get("stored_energy_j"))
        checks["femm_v45_electrostatic_fringe_values"] = (
            voltage > 0.0 and fringe.get("result_voltage_v") == voltage
            and capacitance > 0.0 and fringe.get("result_capacitance_f") == capacitance
            and math.isclose(charge, capacitance * voltage, rel_tol=1e-12)
            and fringe.get("result_charge_c") == charge
            and math.isclose(energy, 0.5 * capacitance * voltage * voltage, rel_tol=1e-12)
            and fringe.get("result_stored_energy_j") == energy
            and fringe.get("result_interface_flux_c") == charge
        )
        checks["femm_v45_electrostatic_fringe_owner"] = (
            math.isclose(float(fringe.get("axisymmetric_factor")), 2.0 * math.pi, rel_tol=1e-12)
            and fringe.get("result_axisymmetric_factor") == fringe.get("axisymmetric_factor")
            and str(fringe.get("mesh_owner", "")).startswith("mesh:")
            and fringe.get("result_mesh_owner") == fringe.get("mesh_owner")
            and _sha(fringe.get("result_sha256"))
            and fringe.get("accepted_result_sha256") == fringe.get("result_sha256")
        )
    return checks
