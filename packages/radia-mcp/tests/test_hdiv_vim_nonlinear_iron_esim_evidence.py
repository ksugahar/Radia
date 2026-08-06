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
    / "results_nonlinear_iron_esim_coupling.json"
)


def test_mcp_local_esim_claim_matches_same_region_iron_artifact():
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    coupled = " ".join(get_hdiv_vim_documentation("eddy_bubble").split())
    rows = result["amplitude_ladder"]

    assert result["pass"] is True
    assert all(result["checks"].values())
    assert result["identity"]["magnetic_region"] == "iron"
    assert result["identity"]["conductive_region"] == "iron"
    assert all(row["local_esim"]["converged"] for row in rows)
    assert all(row["local_surface_impedance"]["passive"] for row in rows)
    assert max(row["mixed"]["residual_relative_norm"] for row in rows) < 1.0e-10
    assert max(
        row["mixed"]["fixed_gram_replay_relative_difference"] for row in rows
    ) < 1.0e-12
    assert "50/1000/5000 A/m" in coupled
    assert "bulk nonlinear B-H operator" in coupled
    assert "is not yet implemented" in coupled
