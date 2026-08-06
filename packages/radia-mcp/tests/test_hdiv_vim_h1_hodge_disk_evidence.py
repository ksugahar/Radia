import json
from pathlib import Path

from radia_mcp.radia_ngsolve.knowledge.hdiv_vim import (
    get_hdiv_vim_documentation,
)


REPO = Path(__file__).resolve().parents[3]
RESULT = (
    REPO
    / "validation_test"
    / "vim_coupled"
    / "results_h1_hodge_bdm2_disk.json"
)


def test_mcp_h1_hodge_disk_claim_matches_the_heavy_artifact():
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    coupled = " ".join(get_hdiv_vim_documentation("eddy_bubble").split())
    rows = result["h1_hodge_3d_cases"]

    assert result["pass"] is True
    assert result["checks"]["bdm2_h_and_h1_p_ladder_is_strict"] is True
    assert result["checks"]["fine_bdm2_h1p4_error_below_one_percent"] is True
    assert rows[-1]["reference_relative_error"] < 0.01
    assert rows[-1]["max_snapshot_relative_residual"] < 1.0e-8
    assert "1.0916977441" in coupled
    assert "4.65%" in coupled
    assert "2.24%" in coupled
    assert "1.44%" in coupled
    assert "0.98%" in coupled
    assert "not a universal open-boundary" in coupled
