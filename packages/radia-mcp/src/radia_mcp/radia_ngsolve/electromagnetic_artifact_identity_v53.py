"""Hysteresis-loss and electrostatic virtual-work identity checks for v53."""

from __future__ import annotations

import math
from collections.abc import Mapping

from .femm_artifact_identity_v54 import validate_public_identity as validate_public_v54_identity


HYSTERESIS = "hysteresis_complex_permeability_phasor_loss_material_owner_identity"
VIRTUAL_WORK = "electrostatic_virtualwork_voltage_charge_constraint_force_owner_identity"
_MU0 = 4.0e-7 * math.pi


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _hysteresis_ok(row: Mapping[str, object]) -> bool:
    permeability = row.get("relative_permeability")
    frequency = row.get("frequency_hz")
    h_rms = row.get("h_rms_a_per_m")
    loss = row.get("loss_density_w_m3")
    permeability_ok = (
        isinstance(permeability, Mapping)
        and set(permeability) == {"real", "imag"}
        and _finite(permeability["real"])
        and _finite(permeability["imag"])
        and float(permeability["real"]) > 0.0
        and float(permeability["imag"]) < 0.0
    )
    expected_loss = (
        2.0 * math.pi * float(frequency) * _MU0 * (-float(permeability["imag"])) * float(h_rms) ** 2
        if permeability_ok and _finite(frequency) and _finite(h_rms)
        else math.nan
    )
    return (
        _generations(row, "permeability_generation", "phasor_generation", "loss_generation", "material_generation", "owner_generation", "result_generation")
        and permeability_ok
        and row.get("result_relative_permeability") == permeability
        and row.get("phasor_convention") == "exp(+j_omega_t)"
        and row.get("result_phasor_convention") == row.get("phasor_convention")
        and _finite(frequency)
        and float(frequency) > 0.0
        and row.get("result_frequency_hz") == frequency
        and _finite(h_rms)
        and float(h_rms) > 0.0
        and row.get("result_h_rms_a_per_m") == h_rms
        and _finite(loss)
        and float(loss) > 0.0
        and math.isclose(float(loss), expected_loss, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and row.get("result_loss_density_w_m3") == loss
        and str(row.get("material_id") or "").startswith("material:")
        and row.get("result_material_id") == row.get("material_id")
        and str(row.get("material_owner") or "").startswith("material-owner:")
        and row.get("result_material_owner") == row.get("material_owner")
        and _result(row)
    )


def _virtual_work_ok(row: Mapping[str, object]) -> bool:
    displacement = row.get("virtual_displacement_m")
    before = row.get("coenergy_before_j")
    after = row.get("coenergy_after_j")
    force = row.get("force_n")
    expected_force = (float(after) - float(before)) / float(displacement) if all(_finite(value) for value in (displacement, before, after)) and float(displacement) != 0.0 else math.nan
    return (
        _generations(row, "path_generation", "constraint_generation", "energy_generation", "force_generation", "owner_generation", "result_generation")
        and row.get("virtual_work_path") == "constant_voltage_coenergy"
        and row.get("result_virtual_work_path") == row.get("virtual_work_path")
        and row.get("constraint_mode") == "fixed_voltage"
        and row.get("result_constraint_mode") == row.get("constraint_mode")
        and _finite(row.get("voltage_v"))
        and float(row["voltage_v"]) != 0.0
        and row.get("result_voltage_v") == row.get("voltage_v")
        and _finite(row.get("charge_c"))
        and row.get("result_charge_c") == row.get("charge_c")
        and _finite(displacement)
        and float(displacement) != 0.0
        and row.get("result_virtual_displacement_m") == displacement
        and _finite(before)
        and _finite(after)
        and row.get("result_coenergy_before_j") == before
        and row.get("result_coenergy_after_j") == after
        and _finite(force)
        and math.isclose(float(force), expected_force, rel_tol=1.0e-10, abs_tol=1.0e-12)
        and row.get("result_force_n") == force
        and str(row.get("force_owner") or "").startswith("force:")
        and row.get("result_force_owner") == row.get("force_owner")
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks = validate_public_v54_identity(identity)
    hysteresis = identity.get(HYSTERESIS)
    virtual_work = identity.get(VIRTUAL_WORK)
    if hysteresis is not None:
        checks["v53_hysteresis_permeability_phasor_loss_material_owner"] = isinstance(hysteresis, Mapping) and _hysteresis_ok(hysteresis)
    if virtual_work is not None:
        checks["v53_electrostatic_virtual_work_constraint_force_owner"] = isinstance(virtual_work, Mapping) and _virtual_work_ok(virtual_work)
    return checks
