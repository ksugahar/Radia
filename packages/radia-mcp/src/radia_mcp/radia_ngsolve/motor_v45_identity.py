from __future__ import annotations

import math
from collections.abc import Mapping


def _sha(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _closed(row: Mapping[str, object], fields: tuple[str, ...]) -> bool:
    generation = str(row.get("generation", "")).strip()
    return bool(generation) and all(row.get(field) == generation for field in fields)


def validate_public_v45_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    ipmsm = identity.get("v45_public_ipmsm_torque_ripple_radial_force_modal_power_efficiency_energy_mesh_owner_mismatch")
    if isinstance(ipmsm, Mapping):
        checks["motor_v45_ipmsm_generation"] = _closed(ipmsm, ("torque_generation", "radial_force_generation", "modal_generation", "power_generation", "efficiency_generation", "energy_generation", "mesh_generation", "result_generation"))
        electromagnetic = float(ipmsm.get("electromagnetic_power_w")); mechanical = float(ipmsm.get("mechanical_power_w")); efficiency = float(ipmsm.get("efficiency")); residual = float(ipmsm.get("energy_closure_residual"))
        checks["motor_v45_ipmsm_values"] = (float(ipmsm.get("torque_ripple_rms_nm")) >= 0.0 and ipmsm.get("result_torque_ripple_rms_nm") == ipmsm.get("torque_ripple_rms_nm") and ipmsm.get("radial_force_space_orders") == ipmsm.get("result_radial_force_space_orders") == [6, 12] and float(ipmsm.get("modal_excitation_n")) >= 0.0 and ipmsm.get("result_modal_excitation_n") == ipmsm.get("modal_excitation_n") and electromagnetic > 0.0 and ipmsm.get("result_electromagnetic_power_w") == electromagnetic and 0.0 <= mechanical <= electromagnetic and ipmsm.get("result_mechanical_power_w") == mechanical and math.isclose(efficiency, mechanical / electromagnetic, rel_tol=1e-12) and ipmsm.get("result_efficiency") == efficiency and 0.0 <= residual <= 1e-8 and ipmsm.get("result_energy_closure_residual") == residual)
        checks["motor_v45_ipmsm_owner"] = str(ipmsm.get("mesh_owner", "")).startswith("mesh:") and ipmsm.get("result_mesh_owner") == ipmsm.get("mesh_owner") and _sha(ipmsm.get("result_sha256")) and ipmsm.get("accepted_result_sha256") == ipmsm.get("result_sha256")
    induction = identity.get("v45_public_induction_machine_slip_rotor_loss_torque_heat_power_energy_result_mismatch")
    if isinstance(induction, Mapping):
        checks["motor_v45_induction_generation"] = _closed(induction, ("slip_generation", "rotor_loss_generation", "torque_generation", "heat_generation", "power_generation", "energy_generation", "result_generation"))
        slip = float(induction.get("slip")); rotor_loss = float(induction.get("rotor_copper_loss_w")); torque = float(induction.get("torque_nm")); heat = float(induction.get("heat_loss_w")); electrical = float(induction.get("electrical_power_w")); mechanical = float(induction.get("mechanical_power_w"))
        checks["motor_v45_induction_values"] = 0.0 <= slip < 1.0 and induction.get("result_slip") == slip and rotor_loss >= 0.0 and induction.get("result_rotor_copper_loss_w") == rotor_loss and torque >= 0.0 and induction.get("result_torque_nm") == torque and heat >= 0.0 and induction.get("result_heat_loss_w") == heat and electrical > 0.0 and 0.0 <= mechanical <= electrical and induction.get("result_mechanical_power_w") == mechanical and math.isclose(electrical, mechanical + rotor_loss + heat, rel_tol=0.0, abs_tol=1e-12)
        checks["motor_v45_induction_owner"] = _sha(induction.get("result_sha256")) and induction.get("accepted_result_sha256") == induction.get("result_sha256")
    return checks
