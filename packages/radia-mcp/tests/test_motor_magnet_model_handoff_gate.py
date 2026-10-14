import json

from radia_mcp.motor.magnet_model_handoff_gate import magnet_model_handoff_gate
from radia_mcp.motor.server import motor_magnet_model_handoff_gate


def _args():
    return dict(
        nonlinear_tolerance=0.01,
        source_result_artifact_id="excitation-field-result",
        source_result_digest="a" * 64,
        magnet_control_artifact_id="magnet-control",
        magnet_control_digest="b" * 64,
        magnet_geometry_artifact_id="magnet-geometry",
        magnet_geometry_digest="c" * 64,
        numbering_policy="preserve",
        element_id_offset=0,
        node_id_offset=0,
        material_mapping_count=4,
        geometry_transform="identity",
    )


def test_magnet_model_handoff_accepts_converged_two_file_package():
    result = magnet_model_handoff_gate([[0.31, 0.082, 0.031], [0.040, 0.0061]], **_args())
    assert result["status"] == "ok"
    assert result["handoff_ready"] is True


def test_magnet_model_handoff_rejects_stale_digest_and_numbering_drift():
    args = _args()
    args.update(source_result_digest="stale", element_id_offset=7, geometry_transform="unknown")
    result = magnet_model_handoff_gate([[0.31, 0.082, 0.031], [0.040, 0.061]], **args)
    assert result["status"] == "needs_attention"
    assert result["checks"]["terminal_residual_meets_tolerance"] is False
    assert result["checks"]["numbering_policy_and_offsets_consistent"] is False


def test_magnet_model_handoff_mcp_dispatches_json():
    result = json.loads(motor_magnet_model_handoff_gate(
        [[0.31, 0.082, 0.031], [0.040, 0.0061]],
        0.01, "excitation-field-result", "a" * 64,
        "magnet-control", "b" * 64, "magnet-geometry", "c" * 64,
        "preserve", 0, 0, 4, "identity",
    ))
    assert result["status"] == "ok"
