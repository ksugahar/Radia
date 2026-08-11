"""Balanced MCP learning-method profile for the public SPICE package."""
from __future__ import annotations

_STAGE_IDS = (
    "baseline_gap", "source_controls", "structured_output", "input_validation",
    "security_boundary", "timeout_cancel_progress", "source_provenance",
    "artifact_feedback", "protocol_smoke", "balance_audit",
)


def build_balanced_learning_profile() -> dict:
    profile = {
        "schema": "cae-ai-lab.balanced-mcp-learning-profile.v1",
        "policy": "equal_capability_gain_v1",
        "server": "radia.ltspice",
        "public_owner": "radia",
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
    profile["self_check"] = validate_balanced_learning_profile(profile)
    return profile


def validate_balanced_learning_profile(profile: dict) -> dict:
    stages = profile.get("stages") if isinstance(profile, dict) else None
    if not isinstance(stages, list):
        stages = []
    ids = [str(row.get("capability_id") or "") for row in stages if isinstance(row, dict)]
    controls_complete = len(stages) == 10 and all(
        isinstance(row, dict)
        and str(row.get("positive_control") or "").strip()
        and str(row.get("negative_control") or "").strip()
        for row in stages
    )
    checks = {
        "stage_ids_match": ids == list(_STAGE_IDS),
        "rounds_ordered": [row.get("round") for row in stages] == list(range(1, 11)),
        "controls_complete": controls_complete,
        "workflow_roles_complete": set(profile.get("workflow_roles") or {})
        == {"detect", "check", "run", "test"},
    }
    return {
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
    }
