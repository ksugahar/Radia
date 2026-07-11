"""Material-parameter gate for variable permanent-magnet models."""

from __future__ import annotations

import math
from typing import Any, Mapping


def variable_magnet_material_parameter_gate(
    parameters: Mapping[str, Any],
    *,
    parameter_authority: str,
    study_label_is_parameter_authority: bool,
) -> dict[str, Any]:
    """Validate sign, ordering, and provenance of a variable-magnet parameter set."""

    try:
        values = {name: float(parameters[name]) for name in ("iHc", "Br", "mur", "mug", "Br0")}
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("parameters must contain numeric iHc, Br, mur, mug, and Br0") from exc
    checks = {
        "parameters_finite": all(math.isfinite(value) for value in values.values()),
        "coercive_field_uses_negative_internal_sign": values["iHc"] < 0.0,
        "remanence_positive": values["Br"] > 0.0,
        "relative_permeability_physical": values["mur"] >= 1.0,
        "gain_not_below_relative_permeability": values["mug"] >= values["mur"],
        "initial_remanence_within_full_remanence": 0.0 <= values["Br0"] <= values["Br"],
        "solver_message_is_parameter_authority": parameter_authority == "saved_solver_report_message",
        "study_label_not_used_as_parameter_value": study_label_is_parameter_authority is False,
    }
    return {
        "policy": "variable_magnet_material_parameter_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "parameters": values,
        "metrics": {
            "initial_to_full_remanence_ratio": values["Br0"] / values["Br"] if values["Br"] else math.inf,
            "gain_to_relative_permeability_ratio": values["mug"] / values["mur"] if values["mur"] else math.inf,
        },
        "lesson": "Use the saved solver message as material-parameter authority; a study label is an identifier and may not encode the effective Br0 value.",
    }
