import copy
import json

from radia_mcp.radia_ngsolve.server import (
    harmonic_magnetic_force_triplet_closure_gate,
)


def summary() -> dict:
    return {
        "frequency_hz": 1000.0,
        "quantity_dimension": "3d_total",
        "force_unit": "N",
        "component_frame": "global",
        "methods": [
            {
                "method": "material_surface",
                "role": "body",
                "force": [0.0, 0.0, -1.0],
            },
            {
                "method": "maxwell_stress",
                "role": "body",
                "force": [0.0, 0.0, -1.02],
            },
            {
                "method": "coil_lorentz",
                "role": "source",
                "force": [0.0, 0.0, 1.021],
            },
        ],
    }


def call(payload: dict) -> dict:
    return json.loads(harmonic_magnetic_force_triplet_closure_gate(json.dumps(payload)))


def test_accepts_independent_body_methods_and_action_reaction_closure():
    result = call(summary())
    assert result["status"] == "ok"
    assert result["metrics"]["dominant_axis"] == "z"
    assert result["metrics"]["action_reaction_relative_residual"] < 0.01


def test_rejects_source_force_with_body_force_sign():
    payload = copy.deepcopy(summary())
    payload["methods"][2]["force"][2] = -1.021
    result = call(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["source_force_opposes_body_force"] is False


def test_rejects_body_method_disagreement():
    payload = copy.deepcopy(summary())
    payload["methods"][0]["force"][2] = -0.8
    result = call(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["body_method_difference_within_tolerance"] is False


def test_rejects_total_force_labeled_as_per_length():
    payload = summary()
    payload["force_unit"] = "N/m"
    result = call(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["force_dimension_and_unit_consistent"] is False
