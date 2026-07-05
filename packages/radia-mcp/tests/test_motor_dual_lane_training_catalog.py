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


def test_dual_lane_training_catalog_has_exactly_30_scrubbed_cases():
    cases = motor_dual_lane_training_catalog()
    assert len(cases) == 30
    gate = motor_dual_lane_training_catalog_gate()
    assert gate["status"] == "PASS"
    assert gate["count"] == 30
    assert gate["forbidden_hits"] == []
    assert gate["checks"]["all_cases_have_age_lane"] is True
    assert gate["checks"]["all_cases_have_vim_lane"] is True


def test_each_training_case_routes_to_age_and_vim():
    for case in motor_dual_lane_training_catalog():
        assert case["radia_motor_age"]["lane_id"] == "radia-motor-age"
        assert case["radia_motor_age"]["validation_lane"] == "ngsolve_age"
        assert case["radia_motor_age"]["targets"]
        assert case["radia_motor_vim"]["lane_id"] == "radia-motor-vim"
        assert case["radia_motor_vim"]["validation_lane"] == "hdiv_vim_reduced_fem"
        assert case["radia_motor_vim"]["targets"]
        assert case["teaching_gate"]


def test_catalog_keeps_source_names_out_of_public_text():
    text = format_motor_dual_lane_training_catalog("all")
    for marker in FORBIDDEN_PUBLIC_MARKERS:
        assert marker not in text
    assert "source-native provenance is private" in text


def test_catalog_closes_external_readable_reference_gap():
    gate = motor_dual_lane_training_catalog_gate()
    assert gate["source_seed_classes"]["external_readable_machine_fem_reference"] >= 4
    titles = "\n".join(case["title"] for case in motor_dual_lane_training_catalog())
    assert "Nonlinear C-core" in titles
    assert "SRM static torque" in titles
    assert "Induction-machine radial air-gap" in titles
    assert "Surface-PM air-gap" in titles


def test_training_route_selects_matching_case_and_both_lanes():
    route = route_dual_lane_training_case("IPM saliency MTPA")
    assert route["selected_case"]["case_id"] == "ipm_saliency_mtpa"
    calls = "\n".join(route["next_public_calls"])
    assert "motor_age_validation_plan" in calls
    assert "hdiv_vim_reduced_fem" in calls
