from __future__ import annotations

import math


def _sha(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, dict):
        return {}
    checks: dict[str, bool] = {}
    ipmsm = identity.get("v44_ipmsm_torque_ripple_radialforce_modal_power_efficiency_energy_mesh_mismatch")
    if isinstance(ipmsm, dict):
        generation = str(ipmsm.get("generation", "")).strip()
        try:
            ripple = float(ipmsm["torque_ripple_rms_nm"]); result_ripple = float(ipmsm["result_torque_ripple_rms_nm"])
            modal = float(ipmsm["modal_excitation_n"]); result_modal = float(ipmsm["result_modal_excitation_n"])
            electromagnetic = float(ipmsm["electromagnetic_power_w"]); result_electromagnetic = float(ipmsm["result_electromagnetic_power_w"])
            mechanical = float(ipmsm["mechanical_power_w"]); result_mechanical = float(ipmsm["result_mechanical_power_w"])
            efficiency = float(ipmsm["efficiency"]); result_efficiency = float(ipmsm["result_efficiency"])
            residual = float(ipmsm["energy_closure_residual"]); result_residual = float(ipmsm["result_energy_closure_residual"])
        except (KeyError, TypeError, ValueError):
            return {"motor_v44_ipmsm_identity": False}
        checks.update({
            "motor_v44_ipmsm_generation_closure": bool(generation) and all(ipmsm.get(k) == generation for k in ("torque_generation", "radial_force_generation", "modal_generation", "power_generation", "efficiency_generation", "energy_generation", "mesh_generation", "result_generation")),
            "motor_v44_ipmsm_force_modal_and_power_closure": ripple >= 0.0 and result_ripple == ripple and modal >= 0.0 and result_modal == modal and electromagnetic > 0.0 and result_electromagnetic == electromagnetic and 0.0 <= mechanical <= electromagnetic and result_mechanical == mechanical,
            "motor_v44_ipmsm_efficiency_and_energy_closure": math.isclose(efficiency, mechanical / electromagnetic, rel_tol=1.0e-12, abs_tol=1.0e-12) and result_efficiency == efficiency and 0.0 <= residual <= 1.0e-8 and result_residual == residual,
            "motor_v44_ipmsm_harmonics_mesh_owner": ipmsm.get("radial_force_space_orders") == [6, 12] and ipmsm.get("result_radial_force_space_orders") == [6, 12] and str(ipmsm.get("mesh_owner", "")).startswith("mesh:") and ipmsm.get("result_mesh_owner") == ipmsm.get("mesh_owner"),
            "motor_v44_ipmsm_result_digest": _sha(ipmsm.get("result_sha256")) and ipmsm.get("accepted_result_sha256") == ipmsm.get("result_sha256"),
        })
    induction = identity.get("v44_inductionmotor_slip_rotorloss_torque_current_power_heat_energy_mismatch")
    if isinstance(induction, dict):
        generation = str(induction.get("generation", "")).strip()
        try:
            slip = float(induction["slip"]); result_slip = float(induction["result_slip"])
            rotor_loss = float(induction["rotor_copper_loss_w"]); result_rotor_loss = float(induction["result_rotor_copper_loss_w"])
            torque = float(induction["torque_nm"]); result_torque = float(induction["result_torque_nm"])
            current = float(induction["stator_current_a"]); result_current = float(induction["result_stator_current_a"])
            electrical = float(induction["electrical_power_w"]); result_electrical = float(induction["result_electrical_power_w"])
            mechanical = float(induction["mechanical_power_w"]); result_mechanical = float(induction["result_mechanical_power_w"])
            heat = float(induction["heat_loss_w"]); result_heat = float(induction["result_heat_loss_w"])
            residual = float(induction["energy_closure_residual"]); result_residual = float(induction["result_energy_closure_residual"])
        except (KeyError, TypeError, ValueError):
            return {"motor_v44_induction_identity": False}
        checks.update({
            "motor_v44_induction_generation_closure": bool(generation) and all(induction.get(k) == generation for k in ("slip_generation", "rotor_loss_generation", "torque_generation", "current_generation", "power_generation", "heat_generation", "energy_generation", "mesh_generation", "result_generation")),
            "motor_v44_induction_slip_current_torque": 0.0 <= slip < 1.0 and result_slip == slip and rotor_loss >= 0.0 and result_rotor_loss == rotor_loss and torque >= 0.0 and result_torque == torque and current > 0.0 and result_current == current,
            "motor_v44_induction_power_heat_closure": electrical > 0.0 and result_electrical == electrical and 0.0 <= mechanical <= electrical and result_mechanical == mechanical and heat >= 0.0 and result_heat == heat and math.isclose(electrical, mechanical + rotor_loss + heat, rel_tol=0.0, abs_tol=1.0e-12),
            "motor_v44_induction_energy_mesh_owner": 0.0 <= residual <= 1.0e-8 and result_residual == residual and str(induction.get("mesh_owner", "")).startswith("mesh:") and induction.get("result_mesh_owner") == induction.get("mesh_owner"),
            "motor_v44_induction_result_digest": _sha(induction.get("result_sha256")) and induction.get("accepted_result_sha256") == induction.get("result_sha256"),
        })
    return checks
