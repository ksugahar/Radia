"""Balanced MCP learning-method profile shared within this repository."""
from __future__ import annotations

_STAGE_SPECS = (
    ("baseline_gap", "Prove the pre-change capability gap.", "A missing behavior is reproduced.", "A claimed gap without a baseline observation is rejected."),
    ("source_controls", "Bind source-native positive and negative controls.", "The owning ecosystem example is accepted.", "A synthetic-only seed is rejected for source learning."),
    ("structured_output", "Publish a stable structured result contract.", "Schema, columns, and units are present.", "A value-only response is rejected."),
    ("input_validation", "Validate required inputs at the MCP boundary.", "A complete typed request is accepted.", "Missing, malformed, or unsafe input is rejected."),
    ("security_boundary", "Keep path, provenance, and execution boundaries explicit.", "An allowed resource stays within policy.", "An arbitrary path or private payload is denied."),
    ("timeout_cancel_progress", "Bound long work and expose progress/cancel semantics.", "A bounded task reports completion or progress.", "An unbounded or unowned cancellation request is rejected."),
    ("source_provenance", "Tie learning to a source-native artifact and version.", "Source identity, version, and run date are recorded.", "Unattributed numbers are rejected."),
    ("artifact_feedback", "Feed verified result artifacts back into MCP behavior.", "A replayable artifact promotes a rule or gate.", "Artifact-only storage without behavior change is rejected."),
    ("protocol_smoke", "Exercise the actual MCP protocol surface.", "tools/list and the changed tool call succeed.", "Unit-test-only evidence without a protocol probe is incomplete."),
    ("balance_audit", "Require equal public and source capability gains.", "Both lanes have probes, tests, and commits.", "A radia-only or source-only change cannot close the round."),
)


def build_balanced_learning_profile(server: str, public_owner: str, source_owner: str) -> dict:
    stages = [
        {
            "round": index,
            "capability_id": capability_id,
            "objective": objective,
            "positive_control": positive,
            "negative_control": negative,
        }
        for index, (capability_id, objective, positive, negative)
        in enumerate(_STAGE_SPECS, start=1)
    ]
    return {
        "schema": "cae-ai-lab.balanced-mcp-learning-profile.v1",
        "policy": "equal_capability_gain_v1",
        "server": server,
        "public_owner": public_owner,
        "source_owner": source_owner,
        "stage_count": len(stages),
        "stages": stages,
        "workflow_roles": {
            "detect": "Observe the baseline and runtime/source availability.",
            "check": "Run typed validation, positive controls, and negative controls.",
            "run": "Execute only through the owning bounded workflow.",
            "test": "Replay protocol, focused tests, artifact gate, and commit evidence.",
        },
        "protocol_policy": {
            "inspector_cli": "Prefer tools/list plus tools/call through MCP Inspector CLI.",
            "conformance": "Rotate official MCP conformance scenarios across one full loop.",
            "fallback": "Document why a direct FastMCP/protocol probe was used instead.",
        },
        "completion_rule": (
            "A stage is verified only when public and source lanes each gain behavior, "
            "pass positive and negative probes, pass focused verification, and record "
            "a non-pending commit."
        ),
    }
