from __future__ import annotations

import copy
import json
from pathlib import Path

from radia_mcp.maglev.knowledge import get_knowledge
from radia_mcp.maglev.server import team28_cycle_averaged_motion_gate as mcp_gate
from radia_mcp.maglev.team28_dynamic_gate import team28_cycle_averaged_motion_gate


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "validation_test"
    / "maglev"
    / "team28_coilbuilder_dynamic_simulink_results.json"
)


def _artifact() -> dict:
    return json.loads(RESULT.read_text(encoding="utf-8"))


def test_gate_accepts_validated_cycle_averaged_mechanical_motion() -> None:
    result = team28_cycle_averaged_motion_gate(_artifact())

    assert result["status"] == "ok"
    assert result["validated_scope"] == "cycle_averaged_mechanical_motion"
    assert result["metrics"]["excitation_frequency_hz"] == 50.0
    assert result["metrics"]["force_family_snapshot_count"] == 25.0
    assert result["metrics"]["eddy_state_order"] == 3.0


def test_mcp_gate_rejects_full_electromagnetic_transient_overclaim() -> None:
    result = json.loads(
        mcp_gate(_artifact(), claim_scope="full_electromagnetic_transient")
    )

    assert result["status"] == "needs_attention"
    assert result["checks"]["claimed_scope_matches_artifact"] is False
    assert "motion_induced_emf" in result["unsupported_claims"]


def test_gate_rejects_frequency_motion_and_damping_lineage_mismatch() -> None:
    artifact = copy.deepcopy(_artifact())
    contract = artifact["model_contract"]
    contract["excitation_frequency_hz"] = 60.0
    contract["motional_emf_included"] = True
    contract["damping_identified_from_measurement"] = True
    artifact["source_artifact"] = r"C:\private\team28.slx"

    result = team28_cycle_averaged_motion_gate(artifact)

    assert result["status"] == "needs_attention"
    assert result["checks"]["excitation_frequency_matches"] is False
    assert result["checks"]["motional_emf_is_explicitly_excluded"] is False
    assert result["checks"]["damping_is_not_claimed_as_identified_data"] is False
    assert result["checks"]["artifact_paths_are_portable_and_relative"] is False


def test_dispatcher_teaches_the_scope_boundary() -> None:
    lesson = get_knowledge("team28_dynamic_scope")

    assert "cycle-averaged lift" in lesson
    assert "not** a full electromagnetic transient" in lesson
    assert "team28_cycle_averaged_motion_gate" in lesson
