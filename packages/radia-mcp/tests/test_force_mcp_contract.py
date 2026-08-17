from __future__ import annotations

import json
from types import SimpleNamespace

from radia_mcp.force import knowledge, server


def test_force_server_exposes_common_knowledge_and_calculation_tools():
    tools = server.mcp._tool_manager._tools
    assert {
        "force",
        "force_methods",
        "force_recipe",
        "force_extras",
        "force_validation_guide",
        "force_result",
        "force_lorentz",
        "force_maxwell_surface",
        "force_time_average_lorentz",
        "force_time_average_maxwell_surface",
        "force_virtual_work",
        "force_coenergy_torque",
        "force_air_gap_torque",
        "force_air_gap_torque_samples",
        "force_time_average_air_gap_torque_samples",
        "force_method_selection_gate",
        "force_method_agreement_gate",
        "force_action_reaction_gate",
        "force_weight_equilibrium_gate",
        "force_status",
        "force_topics",
    }.issubset(tools)


def test_force_is_the_canonical_router_for_compatibility_front_doors():
    from radia_mcp.differential_forms.server import (
        differential_forms_em_force_extras,
        differential_forms_em_force_recipe,
        differential_forms_forces,
    )
    from radia_mcp.maglev.knowledge import get_knowledge as get_maglev_knowledge
    from radia_mcp.motor.server import motor_em_force_extras, motor_em_force_recipe
    from radia_mcp.radia_ngsolve.server import force_validation

    assert motor_em_force_recipe("method_choice") == knowledge.get_force_recipe(
        "method_choice"
    )
    assert motor_em_force_extras("lorentz_canonical") == knowledge.get_force_extras(
        "lorentz_canonical"
    )
    assert differential_forms_forces("overview") == knowledge.get_force_methods(
        "overview"
    )
    assert differential_forms_em_force_recipe("common_pitfalls") == (
        knowledge.get_force_recipe("common_pitfalls")
    )
    assert differential_forms_em_force_extras("meissner_force") == (
        knowledge.get_force_extras("meissner_force")
    )
    assert force_validation("method_map") == knowledge.get_force_validation(
        "method_map"
    )
    maglev_force = get_maglev_knowledge("force_computation")
    assert knowledge.get_force_knowledge("maglev") in maglev_force
    assert knowledge.get_force_recipe("method_choice") in maglev_force


class _Result(list):
    def tolist(self):
        return list(self)


def _build_result(
    force_n,
    torque_nm,
    *,
    method,
    frame,
    pivot_m,
    field_convention,
    amplitude,
    dimensionality="3d",
    per_unit_depth=False,
):
    return {
        "schema": "radia.force-result/v1",
        "method": method,
        "frame": frame,
        "pivot_m": pivot_m,
        "field_convention": field_convention,
        "phasor_amplitude": amplitude,
        "dimensionality": dimensionality,
        "per_unit_depth": per_unit_depth,
        "force_N": (
            None
            if force_n is None
            else force_n.tolist() if hasattr(force_n, "tolist") else list(force_n)
        ),
        "torque_Nm": (
            None
            if torque_nm is None
            else torque_nm.tolist() if hasattr(torque_nm, "tolist") else list(torque_nm)
        ),
    }


def test_force_result_normalizes_application_owned_resultants(monkeypatch):
    monkeypatch.setattr(
        server,
        "_load_force_api",
        lambda: SimpleNamespace(force_torque_result=_build_result),
    )
    payload = json.loads(
        server.force_result(
            None,
            [0.0, 0.0, 4.5],
            "coenergy_virtual_work",
            frame="rotor",
            pivot_m=[0.0, 0.0, 0.1],
            dimensionality="3d",
        )
    )
    assert payload["status"] == "ok"
    assert payload["force_N"] is None
    assert payload["torque_Nm"] == [0.0, 0.0, 4.5]
    assert payload["frame"] == "rotor"


def test_force_lorentz_delegates_to_radia_force(monkeypatch):
    calls = {}

    def integrate(current_density, magnetic_flux_density, weights):
        calls["args"] = (current_density, magnetic_flux_density, weights)
        return _Result([1.0, 2.0, 3.0])

    monkeypatch.setattr(
        server,
        "_load_force_api",
        lambda: SimpleNamespace(
            integrate_lorentz_force=integrate,
            force_torque_result=_build_result,
        ),
    )
    payload = json.loads(
        server.force_lorentz(
            [[1.0, 0.0, 0.0]],
            [[0.0, 1.0, 0.0]],
            [0.5],
        )
    )

    assert payload["status"] == "ok"
    assert payload["method"] == "lorentz_body_force"
    assert payload["force_N"] == [1.0, 2.0, 3.0]
    assert calls["args"][2] == [0.5]


def test_force_maxwell_surface_delegates_to_radia_force(monkeypatch):
    calls = {}

    def integrate(field, normals, weights, *, permeability_H_per_m):
        calls["permeability"] = permeability_H_per_m
        return _Result([0.0, 0.0, 4.0])

    monkeypatch.setattr(
        server,
        "_load_force_api",
        lambda: SimpleNamespace(
            integrate_maxwell_surface_force=integrate,
            force_torque_result=_build_result,
        ),
    )
    payload = json.loads(
        server.force_maxwell_surface(
            [[0.0, 0.0, 1.0]],
            [[0.0, 0.0, 1.0]],
            [2.0],
        )
    )

    assert payload["status"] == "ok"
    assert payload["method"] == "maxwell_surface_stress_air"
    assert payload["force_N"] == [0.0, 0.0, 4.0]
    assert calls["permeability"] == payload["permeability_H_per_m"]


def test_force_lorentz_returns_torque_when_sample_points_are_supplied(monkeypatch):
    def integrate(*args, **kwargs):
        return _Result([1.0, 0.0, 0.0]), _Result([0.0, 0.0, -2.0])

    monkeypatch.setattr(
        server,
        "_load_force_api",
        lambda: SimpleNamespace(
            integrate_lorentz_force_and_torque=integrate,
            force_torque_result=_build_result,
        ),
    )
    payload = json.loads(
        server.force_lorentz(
            [[0.0, 0.0, 1.0]],
            [[0.0, 1.0, 0.0]],
            [1.0],
            sample_points_m=[[0.0, 2.0, 0.0]],
            pivot_m=[0.0, 0.0, 0.0],
            frame="rotor",
        )
    )
    assert payload["force_N"] == [1.0, 0.0, 0.0]
    assert payload["torque_Nm"] == [0.0, 0.0, -2.0]
    assert payload["frame"] == "rotor"


def test_common_force_gates_cover_motor_and_maglev_acceptance():
    method = json.loads(
        server.force_method_selection_gate(
            "coil",
            "lorentz_body_force",
        )
    )
    assert method["status"] == "ok"

    first = {
        "method": "lorentz_body_force",
        "frame": "global_cartesian",
        "pivot_m": [0.0, 0.0, 0.0],
        "dimensionality": "3d",
        "per_unit_depth": False,
        "field_convention": "static",
        "force_N": [10.0, 0.0, 0.0],
        "torque_Nm": [0.0, 0.0, 2.0],
    }
    second = {
        **first,
        "method": "maxwell_surface_stress_air",
        "force_N": [9.9, 0.0, 0.0],
        "torque_Nm": [0.0, 0.0, 2.01],
    }
    agreement = json.loads(server.force_method_agreement_gate(first, second))
    assert agreement["status"] == "ok"

    reaction = json.loads(
        server.force_action_reaction_gate([3.0, 0.0, 0.0], [-3.0, 0.0, 0.0])
    )
    assert reaction["status"] == "ok"

    equilibrium = json.loads(
        server.force_weight_equilibrium_gate([0.0, 0.0, 9.80665], 1.0)
    )
    assert equilibrium["status"] == "ok"


def test_motor_and_maglev_force_gates_forward_to_common_layer():
    from radia_mcp.maglev.server import (
        maglev_force_torque_method_agreement_gate,
        maglev_force_weight_equilibrium_gate,
    )
    from radia_mcp.motor.server import motor_force_torque_method_agreement_gate

    primary = {
        "method": "lorentz_body_force",
        "frame": "global_cartesian",
        "pivot_m": [0.0, 0.0, 0.0],
        "dimensionality": "3d",
        "per_unit_depth": False,
        "field_convention": "static",
        "force_N": [10.0, 0.0, 0.0],
        "torque_Nm": [0.0, 0.0, 2.0],
    }
    independent = {
        **primary,
        "method": "maxwell_surface_stress_air",
        "force_N": [9.9, 0.0, 0.0],
        "torque_Nm": [0.0, 0.0, 2.01],
    }
    common = json.loads(server.force_method_agreement_gate(primary, independent))
    motor = json.loads(
        motor_force_torque_method_agreement_gate(primary, independent)
    )
    maglev = json.loads(
        maglev_force_torque_method_agreement_gate(primary, independent)
    )
    assert motor == common
    assert maglev == common

    common_equilibrium = json.loads(
        server.force_weight_equilibrium_gate([0.0, 0.0, 9.80665], 1.0)
    )
    maglev_equilibrium = json.loads(
        maglev_force_weight_equilibrium_gate([0.0, 0.0, 9.80665], 1.0)
    )
    assert maglev_equilibrium == common_equilibrium


def test_common_force_gates_reject_vacuous_and_nonfinite_evidence():
    vacuous = json.loads(server.force_method_agreement_gate({}, {}))
    assert vacuous["status"] == "needs_attention"
    assert not vacuous["checks"]["at_least_one_comparable_resultant"]

    nonfinite = json.loads(
        server.force_method_agreement_gate(
            {"method": "a", "force_N": [1.0, 0.0, 0.0]},
            {"method": "b", "force_N": [1.0, 0.0, 0.0]},
            maximum_force_relative_difference=float("nan"),
        )
    )
    assert nonfinite["status"] == "invalid_input"

    incomplete_torque_pair = json.loads(
        server.force_action_reaction_gate(
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            torque_a_nm=[0.0, 0.0, 1.0],
        )
    )
    assert incomplete_torque_pair["status"] == "invalid_input"


def test_force_calculators_fail_loudly_when_radia_is_unavailable(monkeypatch):
    def unavailable():
        raise RuntimeError("install radia")

    monkeypatch.setattr(server, "_load_force_api", unavailable)
    payload = json.loads(server.force_lorentz([], [], []))
    assert payload["status"] == "unavailable"
    assert "install radia" in payload["error"]
