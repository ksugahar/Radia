import copy
import json

from radia_mcp.build123d.external_cad_gate import (
    build123d_upstream_example_roundtrip_gate,
    external_cad_mass_topology_crosscheck_gate,
    step_portability_diagnosis_gate,
)
from radia_mcp.build123d.server import (
    build123d_external_cad_mass_topology_gate as mcp_external_gate,
    build123d_step_portability_diagnosis_gate as mcp_portability_gate,
    build123d_upstream_example_roundtrip_gate as mcp_upstream_gate,
)


STEP_SHA = "a" * 64
TOPOLOGY = {"vertices": 36, "edges": 54, "faces": 25, "solids": 1, "euler_characteristic": 7}


def _upstream():
    native = {
        "volume": 44436.46039467075,
        "area": 13163.451211657966,
        "centroid": [0.0, 0.0, 4.902011830393274],
        "bbox_size": [80.0, 60.0, 10.0],
        **TOPOLOGY,
    }
    imported = copy.deepcopy(native)
    imported["volume"] *= 1.0 + 4.0e-15
    imported["area"] *= 1.0 + 4.0e-15
    return {
        "source_kind": "upstream_native_example",
        "source_sha256": "b" * 64,
        "upstream_commit": "c" * 40,
        "build123d_version": "0.10.0",
        "step_sha256": STEP_SHA,
        "native": native,
        "roundtrip": imported,
        "timings_s": {"source_build": 0.03, "step_export": 0.02, "step_reimport": 0.04},
    }


def _reference():
    return {
        "source": "build123d",
        "step_sha256": STEP_SHA,
        "volume": 44436.46039467075,
        "area": 13163.451211657966,
        "center_of_mass": [0.0, 0.0, 4.902011830393274],
        "bbox_size": [80.0, 60.0, 10.0],
        **TOPOLOGY,
    }


def _external():
    return {
        "source": "external_cad",
        "step_sha256": STEP_SHA,
        "volume": 44436.51846333852,
        "area": 13163.451211657968,
        "representative_center": [0.0, 0.0, 5.0],
        "center_semantics": "entity_center_excluded",
        "bbox_size": [80.0, 60.0, 10.0],
        **TOPOLOGY,
    }


def test_upstream_roundtrip_live_shape_passes_and_dispatches():
    result = build123d_upstream_example_roundtrip_gate(_upstream())
    assert result["status"] == "ok"
    assert result["metrics"]["euler_characteristic"] == 7
    assert json.loads(mcp_upstream_gate(json.dumps(_upstream())))["status"] == "ok"


def test_upstream_roundtrip_rejects_topology_drift():
    bad = _upstream()
    bad["roundtrip"]["faces"] = 24
    result = build123d_upstream_example_roundtrip_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["brep_topology_roundtrip_matches"] is False


def test_external_kernel_passes_with_entity_center_explicitly_excluded():
    result = external_cad_mass_topology_crosscheck_gate(_reference(), _external())
    assert result["status"] == "ok"
    assert result["metrics"]["center_comparison_performed"] is False
    assert result["metrics"]["volume_relative_error"] < 2.0e-6
    assert json.loads(mcp_external_gate(json.dumps(_reference()), json.dumps(_external())))["status"] == "ok"


def test_external_kernel_rejects_ambiguous_center_and_topology_drift():
    bad = _external()
    bad["center_semantics"] = "centroid"
    bad["faces"] = 24
    result = json.loads(mcp_external_gate(json.dumps(_reference()), json.dumps(bad)))
    assert result["status"] == "invalid_input"

    bad["center_semantics"] = "entity_center_excluded"
    result = external_cad_mass_topology_crosscheck_gate(_reference(), bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["brep_topology_matches"] is False


def _portability(external_volume=5767.508755516695):
    return {
        "source_kind": "upstream_native_example",
        "upstream_commit": "c" * 40,
        "source_sha256": "b" * 64,
        "step_sha256": "a" * 64,
        "native_volume": 5767.5000452165295,
        "self_roundtrip_volume": 5767.50004521653,
        "external_imports": [
            {"mode": "heal", "volume": external_volume, "volume_count": 1},
            {"mode": "noheal", "volume": external_volume, "volume_count": 1},
        ],
    }


def test_step_portability_diagnosis_accepts_cross_kernel_control_and_dispatches():
    payload = _portability()
    result = step_portability_diagnosis_gate(payload)
    assert result["status"] == "ok"
    assert result["diagnosis"] == "portable"
    assert json.loads(mcp_portability_gate(json.dumps(payload)))["status"] == "ok"


def test_step_portability_diagnosis_locates_external_translation_loss():
    payload = _portability(external_volume=282701.24890846136)
    payload["native_volume"] = 363795.07369811134
    payload["self_roundtrip_volume"] = 363795.073698225
    result = step_portability_diagnosis_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["diagnosis"] == "external_kernel_translation_loss"
    assert result["healing_not_root_cause"] is True
    assert result["checks"]["self_roundtrip_matches"] is True
    assert result["checks"]["external_volume_matches"] is False
