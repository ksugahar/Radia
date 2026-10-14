from __future__ import annotations

import cmath
import json

import pytest

from radia_mcp.motor.server import motor_circuit_age_application_plan
from radia_mcp.radia_ngsolve.airgap_element import (
    annular_rotation_phase,
    planar_translation_phase,
)
from radia_mcp.radia_ngsolve.circuit_excitation import (
    compile_circuit_age_application,
)


def _series_payload(problem_kind: str, element_family: str = "P2_curved") -> dict:
    return {
        "problem_kind": problem_kind,
        "element_family": element_family,
        "circuits": [
            {
                "name": "phase_a",
                "connection": "series",
                "current_a": [3.0, -4.0],
                "frequency_hz": 50.0,
            }
        ],
        "regions": [
            {
                "name": "coil_plus",
                "circuit": "phase_a",
                "turns": 20,
                "area_m2": 2.0e-4,
            },
            {
                "name": "coil_minus",
                "circuit": "phase_a",
                "turns": -20,
                "area_m2": 2.0e-4,
            },
        ],
    }


@pytest.mark.parametrize("problem_kind", ["planar_2d", "axisymmetric"])
def test_series_current_density_is_shared_by_planar_and_axisymmetric(problem_kind):
    result = compile_circuit_age_application(_series_payload(problem_kind))

    assert result["status"] == "compiled"
    assert result["problem_kind"] == problem_kind
    assert result["element_family"] == "P2_curved"
    assert result["assembly_contract"]["kelvin_open_boundary_compatible"] is True
    assert result["assembly_contract"]["eddy_transient_compatible"] is True
    assert result["circuits"][0]["additional_circuit_unknowns"] == []
    expected = [20.0 * 3.0 / 2.0e-4, 20.0 * -4.0 / 2.0e-4]
    assert result["regions"][0]["impressed_current_density_a_per_m2"] == expected
    assert result["regions"][1]["impressed_current_density_a_per_m2"] == [
        -expected[0],
        -expected[1],
    ]


def test_parallel_current_is_a_constraint_not_an_equal_split():
    payload = {
        "problem_kind": "planar_2d",
        "element_family": "Q2",
        "circuits": [
            {"name": "bars", "connection": "parallel", "current_a": 12.0}
        ],
        "regions": [
            {
                "name": "bar_wide",
                "circuit": "bars",
                "turns": 1,
                "area_m2": 4.0e-4,
                "conductivity_s_per_m": 5.8e7,
            },
            {
                "name": "bar_narrow",
                "circuit": "bars",
                "turns": 1,
                "area_m2": 1.0e-4,
                "conductivity_s_per_m": 5.8e7,
            },
        ],
    }
    result = compile_circuit_age_application(payload)
    circuit = result["circuits"][0]

    assert circuit["equal_current_split_assumed"] is False
    assert circuit["current_constraint"] == "sum(I_branch(region)) = I_circuit"
    assert circuit["additional_circuit_unknowns"] == [
        "V_common:bars",
        "I_branch:bar_wide",
        "I_branch:bar_narrow",
    ]
    assert all(
        region["impressed_current_density_a_per_m2"] is None
        for region in result["regions"]
    )


def test_rotary_age_uses_harmonic_phase_without_mesh_rebuild():
    payload = _series_payload("planar_2d", "P1")
    payload["motion"] = {
        "kind": "annular_age",
        "position_rad": 0.25,
        "harmonics": [1, 3],
    }
    motion = compile_circuit_age_application(payload)["motion"]

    assert motion["mesh_rebuild"] is False
    assert motion["mechanical_observable"] == "mesh-independent harmonic torque"
    factor = complex(*motion["phase_factors"]["3"])
    assert factor == pytest.approx(annular_rotation_phase(3, 0.25))
    assert factor == pytest.approx(cmath.exp(-1j * 3.0 * 0.25))


def test_planar_age_translation_uses_wavenumber_phase():
    payload = _series_payload("planar_2d", "Q1")
    payload["motion"] = {
        "kind": "planar_age",
        "position_m": 0.002,
        "wavenumbers_per_m": [100.0, 200.0],
    }
    motion = compile_circuit_age_application(payload)["motion"]

    assert motion["mesh_rebuild"] is False
    assert motion["mechanical_observable"] == "mesh-independent harmonic thrust"
    factor = complex(*motion["phase_factors"]["100"])
    assert factor == pytest.approx(planar_translation_phase(100.0, 0.002))
    assert factor == pytest.approx(cmath.exp(-1j * 100.0 * 0.002))


def test_age_phase_helpers_reject_invalid_modes_without_touching_a_mesh():
    with pytest.raises(ValueError, match="positive integer"):
        annular_rotation_phase(0, 0.0)
    with pytest.raises(ValueError, match="positive and finite"):
        planar_translation_phase(-1.0, 0.0)


def test_orphan_circuit_and_axisymmetric_age_fail_loudly():
    orphan = _series_payload("planar_2d")
    orphan["regions"][0]["circuit"] = "undefined"
    with pytest.raises(ValueError, match="undefined circuit"):
        compile_circuit_age_application(orphan)

    axis_age = _series_payload("axisymmetric")
    axis_age["motion"] = {
        "kind": "annular_age",
        "position_rad": 0.0,
        "harmonics": [1],
    }
    with pytest.raises(ValueError, match="planar_2d"):
        compile_circuit_age_application(axis_age)


def test_motor_mcp_tool_exposes_compiler_and_rejects_bad_json():
    compiled = json.loads(
        motor_circuit_age_application_plan(json.dumps(_series_payload("planar_2d")))
    )
    assert compiled["status"] == "compiled"
    assert compiled["schema"] == "radia.circuit-age-application.v1"

    invalid = json.loads(motor_circuit_age_application_plan("{"))
    assert invalid["status"] == "invalid_input"
