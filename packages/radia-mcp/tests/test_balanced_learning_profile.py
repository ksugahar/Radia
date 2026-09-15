from radia_mcp.radia_ngsolve.server import balanced_mcp_learning_profile
from radia_mcp.common.learning_quality import validate_balanced_learning_profile


def test_radia_mcp_exposes_balanced_learning_profile_tool():
    profile = balanced_mcp_learning_profile()
    assert profile["policy"] == "equal_capability_gain_v1"
    assert profile["stage_count"] == 10
    assert len({row["capability_id"] for row in profile["stages"]}) == 10
    assert set(profile["workflow_roles"]) == {"detect", "check", "run", "test"}
    assert profile["self_check"]["status"] == "ok"
    bad = {**profile, "stages": [dict(row) for row in profile["stages"]]}
    bad["stages"][6]["negative_control"] = ""
    assert validate_balanced_learning_profile(bad)["status"] == "needs_attention"
