"""Balanced MCP learning-method profile for the public SPICE package."""
from __future__ import annotations

_STAGE_IDS = (
    "baseline_gap", "source_controls", "structured_output", "input_validation",
    "security_boundary", "timeout_cancel_progress", "source_provenance",
    "artifact_feedback", "protocol_smoke", "balance_audit",
)


def build_balanced_learning_profile() -> dict:
    return {
        "schema": "cae-ai-lab.balanced-mcp-learning-profile.v1",
        "policy": "equal_capability_gain_v1",
        "server": "radia-spice-lab",
        "public_owner": "radia-spice-lab",
        "source_owner": "private LTSpice source MCP",
        "stage_count": len(_STAGE_IDS),
        "stages": [
            {
                "round": index,
                "capability_id": capability_id,
                "objective": f"Verify {capability_id.replace('_', ' ')} learning behavior.",
                "positive_control": "Complete source-native evidence is accepted.",
                "negative_control": "Incomplete or one-sided evidence is rejected.",
            }
            for index, capability_id in enumerate(_STAGE_IDS, start=1)
        ],
        "workflow_roles": {
            "detect": "Observe baseline and available parser/runtime capabilities.",
            "check": "Run typed positive and negative controls.",
            "run": "Use bounded public-safe conversion and parsing tools.",
            "test": "Run protocol, focused tests, artifact gate, and commit evidence.",
        },
        "protocol_policy": {
            "inspector_cli": "Prefer tools/list plus tools/call through MCP Inspector CLI.",
            "conformance": "Rotate official MCP conformance scenarios across a full loop.",
            "fallback": "Record why a direct FastMCP/protocol probe was used instead.",
        },
        "completion_rule": "Both public and source lanes need behavior, controls, tests, and non-pending commits.",
    }
