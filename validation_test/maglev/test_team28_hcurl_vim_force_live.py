"""Live numerical regression, separate from stored TEAM28 evidence checks.

Run explicitly on an idle compute host; this is not ordinary fast CI.
"""
import json
from validation_test.maglev.team28_hcurl_vim_force import run


def test_team28_hcurl_vim_force_live(tmp_path):
    result = run([0.025], outer_quad=4)
    output = tmp_path / "team28_hcurl_vim_force_live.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    assert len(result["cases"]) == 1
    assert result["hcurl_vim_force_acceptance_complete"], result["checks"]
    assert all(result["checks"].values())
