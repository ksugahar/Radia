"""Regression tests for neutral magnetic-pole and CoilBuilder exchange."""

from __future__ import annotations

import copy
import json

from radia.coil_builder import CoilBuilder
from radia_mcp.electromagnet.coil_exchange import build_shared_pole_coil_exchange
from radia_mcp.electromagnet.server import electromagnet_shared_pole_coil_exchange


def _contract() -> dict:
    return {
        "schema": "cae-ai-lab.shared-pole-coil.v1",
        "component_label": "test_dipole",
        "coordinate_system": "right-handed Cartesian",
        "length_unit": "m",
        "current_unit": "A",
        "bodies": [
            {
                "name": "upper_pole",
                "role": "magnetic_pole",
                "kind": "axis_aligned_box",
                "lower_m": [-0.04, -0.02, 0.01],
                "upper_m": [0.04, 0.02, 0.03],
                "material_label": "pole_steel",
            }
        ],
        "coils": [
            {
                "name": "main_winding",
                "kind": "rounded_rectangle_racetrack",
                "corner_centerline_radius_m": 0.012,
                "straight_x_m": 0.08,
                "straight_y_m": 0.12,
                "cross_section_width_m": 0.004,
                "cross_section_height_m": 0.006,
                "centre_m": [0.0, 0.0, 0.0],
                "orientation_rows": [
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
                "current_A": 250.0,
                "turns": 8,
                "current_direction": "segment_path_start_to_end",
                "conductor_model": "solid_uniform_current_density",
                "region_label": "main_winding",
                "current_density_A_per_m2": 250.0 / (0.004 * 0.006),
            }
        ],
    }


def test_shared_contract_builds_closed_coilbuilder_and_modelir_fragment():
    result = build_shared_pole_coil_exchange(_contract())

    assert result["status"] == "accepted"
    assert result["accepted"] is True
    assert len(result["identity_sha256"]) == 64
    assert result["shared_contract"]["bodies"][0]["role"] == "magnetic_pole"
    identity = result["radia_adapters"][0]["coil_identity"]
    assert identity["closed"] is True
    assert len(identity["segments"]) == 8
    assert identity["ampere_turns_A"] == 2000.0
    coil = _contract()["coils"][0]
    native = CoilBuilder.rounded_rectangle_racetrack(
        coil["current_A"],
        corner_centerline_radius=coil["corner_centerline_radius_m"],
        straight_x=coil["straight_x_m"],
        straight_y=coil["straight_y_m"],
        width=coil["cross_section_width_m"],
        height=coil["cross_section_height_m"],
        centre=coil["centre_m"],
        orientation=coil["orientation_rows"],
    ).to_identity_manifest(turns=coil["turns"], region_label=coil["region_label"])
    assert identity == native
    fragment = result["modelir_fragment"]
    assert fragment["geometry"]["bodies"][0]["name"] == "upper_pole"
    assert fragment["sources"][0]["kind"] == "coil"
    assert result["constitutive_policy"]["changed_by_exchange"] is False


def test_stale_identity_rejects_current_reversal():
    positive = build_shared_pole_coil_exchange(_contract())
    reversed_contract = copy.deepcopy(_contract())
    reversed_contract["identity_sha256"] = positive["identity_sha256"]
    reversed_contract["coils"][0]["current_A"] *= -1.0
    reversed_contract["coils"][0]["current_density_A_per_m2"] *= -1.0

    result = build_shared_pole_coil_exchange(reversed_contract)

    assert result["status"] == "rejected"
    assert "identity_sha256" in result["issues"][0]


def test_left_handed_frame_and_unsupported_geometry_are_rejected():
    left_handed = _contract()
    left_handed["coils"][0]["orientation_rows"][2][2] = -1.0
    result = build_shared_pole_coil_exchange(left_handed)
    assert result["status"] == "rejected"
    assert "right-handed" in result["issues"][0]

    unsupported = _contract()
    unsupported["bodies"][0]["kind"] = "freeform_surface"
    result = build_shared_pole_coil_exchange(unsupported)
    assert result["status"] == "rejected"
    assert "unsupported body kind" in result["issues"][0]


def test_mcp_wrapper_returns_the_same_exchange_contract():
    result = electromagnet_shared_pole_coil_exchange(json.dumps(_contract()))

    assert result["status"] == "accepted"
    assert result["policy"] == "shared_pole_coil_exchange_v1"


def test_current_density_must_match_equivalent_current():
    payload = _contract()
    payload["coils"][0]["current_density_A_per_m2"] *= 1.01

    result = build_shared_pole_coil_exchange(payload)

    assert result["status"] == "rejected"
    assert "current density" in result["issues"][0]
