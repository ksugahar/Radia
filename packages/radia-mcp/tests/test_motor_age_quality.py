# -*- coding: utf-8 -*-
"""Fast tests for the radia-motor AGE quality routing layer."""

import os
import sys


_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.motor.age_quality_knowledge import (  # noqa: E402
    format_age_validation_plan,
    get_age_quality_report,
    route_age_validation_plan,
)


def test_age_quality_report_names_public_gates():
    report = get_age_quality_report("gate_matrix")
    assert "age_analytic_dtn" in report
    assert "age_eddy_machine" in report
    assert "validation_test/radia_mcp/test_airgap_eddy_machine.py" in report
    assert "gold_age_invariant" in get_age_quality_report("publication_policy")


def test_ipm_plan_routes_to_dq_and_field_weakening_gates():
    plan = route_age_validation_plan("IPM hairpin motor MTPA field weakening")
    assert plan["family"] == "ipm"
    assert "dq_control_layer" in plan["required_gate_ids"]
    assert "age_ipm_synchronous_torque" in plan["required_gate_ids"]
    assert "validation_test/radia_mcp/test_field_weakening.py" in plan["pytest_targets"]
    text = format_age_validation_plan(plan)
    assert "Ld/Lq saliency" in text
    assert "MTPA current angle" in text


def test_induction_plan_routes_to_eddy_and_slip_gates():
    plan = route_age_validation_plan("induction cage rotor slip loss")
    assert plan["family"] == "induction"
    assert "age_eddy_machine" in plan["required_gate_ids"]
    assert "induction_slip_layer" in plan["required_gate_ids"]
    assert "validation_test/radia_mcp/test_motor_induction_coupling.py" in plan["pytest_targets"]
    text = format_age_validation_plan(plan)
    assert "rotor eddy loss" in text
    assert "fully meshed reference" in text


def test_runbook_does_not_embed_private_absolute_paths():
    runbook = get_age_quality_report("runbook")
    for drive in ("S", "W", "C"):
        assert f"{drive}:" + "\\" not in runbook
