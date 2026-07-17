# -*- coding: utf-8 -*-
"""Fast tests for the public-safe radia-motor dual-lane catalog."""

import os
import sys


_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.motor.dual_lane_training_catalog import (  # noqa: E402
    FORBIDDEN_PUBLIC_MARKERS,
    format_motor_dual_lane_training_catalog,
    motor_dual_lane_training_catalog,
    motor_dual_lane_training_catalog_gate,
    route_dual_lane_training_case,
)


def test_dual_lane_training_catalog_has_wide_scrubbed_case_set():
    cases = motor_dual_lane_training_catalog()
    assert len(cases) >= 50
    gate = motor_dual_lane_training_catalog_gate()
    assert gate["status"] == "PASS"
    assert gate["count"] >= 50
    assert gate["forbidden_hits"] == []
    assert gate["checks"]["all_cases_have_age_lane"] is True
    assert gate["checks"]["all_cases_have_mmm_eddy_lane"] is True
    assert gate["checks"]["covers_wide_machine_families"] is True


def test_each_training_case_routes_to_age_and_mmm_eddy():
    for case in motor_dual_lane_training_catalog():
        assert case["radia_motor_age"]["lane_id"] == "radia-motor-age"
        assert case["radia_motor_age"]["validation_lane"] == "ngsolve_age"
        assert case["radia_motor_age"]["targets"]
        assert case["radia_motor_mmm_eddy"]["lane_id"] == "radia-motor-mmm-eddy"
        assert (
            case["radia_motor_mmm_eddy"]["validation_lane"]
            == "hdiv_mmm_hcurl_eddy_bubble"
        )
        assert case["radia_motor_mmm_eddy"]["targets"]
        assert case["teaching_gate"]


def test_catalog_keeps_source_names_out_of_public_text():
    text = format_motor_dual_lane_training_catalog("all")
    for marker in FORBIDDEN_PUBLIC_MARKERS:
        assert marker not in text
    assert "source-native provenance is private" in text


def test_catalog_closes_external_readable_reference_gap():
    gate = motor_dual_lane_training_catalog_gate()
    assert gate["source_seed_classes"]["external_readable_machine_fem_reference"] >= 4
    assert gate["source_seed_classes"]["external_current_machine_example_library"] >= 20
    titles = "\n".join(case["title"] for case in motor_dual_lane_training_catalog())
    assert "Nonlinear C-core" in titles
    assert "SRM static torque" in titles
    assert "Induction-machine radial air-gap" in titles
    assert "Surface-PM air-gap" in titles
    assert "Axial-flux PM" in titles
    assert "Outer-rotor BLDC" in titles
    assert "Doubly-fed induction generator" in titles
    assert "Distributed winding-function" in titles


def test_training_route_selects_matching_case_and_both_lanes():
    route = route_dual_lane_training_case("IPM saliency MTPA")
    assert route["selected_case"]["case_id"] == "ipm_saliency_mtpa"
    calls = "\n".join(route["next_public_calls"])
    assert "motor_age_validation_plan" in calls
    assert "hdiv_mmm_hcurl_eddy_bubble" in calls


def test_training_route_covers_wide_machine_families():
    route = route_dual_lane_training_case("BLDC outer rotor polarity")
    assert route["selected_case"]["case_id"] == "bldc_outer_rotor_polarity"
    route = route_dual_lane_training_case("DFIG slip power")
    assert route["selected_case"]["case_id"] == "dfig_slip_power_coupling"
    route = route_dual_lane_training_case("distributed winding function")
    assert route["selected_case"]["case_id"] == "winding_function_distributed_layout"
