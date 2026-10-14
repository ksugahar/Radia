import copy
import json

from radia_mcp.build123d.path_sweep_gate import (
    build123d_path_sweep_handoff_gate,
    build123d_path_sweep_source_contract_gate,
)
from radia_mcp.build123d.server import (
    build123d_path_sweep_handoff_gate as mcp_handoff_gate,
    build123d_path_sweep_source_contract_gate as mcp_source_gate,
)


STEP_SHA = "4" * 64


def _handoff():
    return {
        "length_unit": "mm",
        "step_sha256": STEP_SHA,
        "volume_oracle_policy": "cross_kernel_mass_properties",
        "analytic": {
            "path_length_mm": 245.06410171031138,
            "section_area_mm2": 400.0,
        },
        "native": {
            "path_length_mm": 245.06410171031138,
            "volume_mm3": 91398.2236861551,
            "area_mm2": 19742.386437027966,
            "bbox_size_mm": [170.0, 100.0, 20.0],
            "solid_count": 1,
            "is_valid": True,
        },
        "external": {
            "step_sha256": STEP_SHA,
            "volume_mm3": 91395.74122778537,
            "area_mm2": 19742.386437027955,
            "bbox_size_mm": [170.0, 100.0, 20.0],
            "volume_count": 1,
            "euler_characteristic": 2,
        },
    }


def _source():
    return {
        "ecosystem": "build123d",
        "example_id": "introductory-ex14",
        "api_mode": "algebra",
        "path_segments": ["JernArc", "JernArc", "Line"],
        "explicit_path_keyword": True,
        "profile_plane": "XZ",
        "profile_type": "Rectangle",
        "build123d_version": "0.10.0",
        "validity_attribute_type": "bool",
        "validity_access": "property",
        "run_error": "",
    }


def test_path_sweep_handoff_accepts_live_shape_and_does_not_promote_naive_volume():
    result = build123d_path_sweep_handoff_gate(_handoff())
    assert result["status"] == "ok"
    assert result["naive_tube_volume_is_oracle"] is False
    assert result["metrics"]["naive_tube_volume_relative_gap"] > 0.07
    assert json.loads(mcp_handoff_gate(json.dumps(_handoff())))["status"] == "ok"


def test_path_sweep_handoff_rejects_external_volume_loss_and_naive_oracle():
    bad = copy.deepcopy(_handoff())
    bad["external"]["volume_mm3"] *= 0.98
    bad["volume_oracle_policy"] = "section_area_times_path_length"
    result = build123d_path_sweep_handoff_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["external_volume_matches"] is False
    assert result["checks"]["finished_volume_uses_cross_kernel_oracle"] is False


def test_path_sweep_source_contract_accepts_build123d_010_property_access():
    result = build123d_path_sweep_source_contract_gate(_source())
    assert result["status"] == "ok"
    assert result["diagnosis"] == "source_contract_ok"
    assert json.loads(mcp_source_gate(json.dumps(_source())))["status"] == "ok"


def test_path_sweep_source_contract_diagnoses_bool_property_called_as_method():
    bad = _source()
    bad["validity_access"] = "method"
    bad["run_error"] = "TypeError: 'bool' object is not callable"
    result = build123d_path_sweep_source_contract_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["diagnosis"] == "is_valid_property_called_as_method"
