import json
from pathlib import Path


SUMMARY = (
    Path(__file__).resolve().parents[1]
    / "validation_test"
    / "maglev"
    / "team28_hcurl_vim_force_summary.json"
)


def test_team28_p6_hcurl_vim_force_acceptance_record():
    result = json.loads(SUMMARY.read_text(encoding="utf-8"))

    assert result["schema"] == "radia.team28.hcurl-vim-force.v1"
    assert result["runtime"]["hostname"] == "mdx"
    assert result["hcurl_vim_force_acceptance_complete"] is True
    assert all(result["checks"].values())
    assert len(result["cases"]) == 3
    assert result["maximum_force_relative_error"] < 0.01
    assert result["outer_quadrature_relative_force_change"] < 0.001
    assert all(case["parent_ndof"] > 20_000 for case in result["cases"])
    assert all(case["evrs_rank"] == 3 for case in result["cases"])
    assert all(
        case["cln_handoff"]["state_order"] == 3
        and case["cln_handoff"]["port_count"] == 1
        and case["cln_handoff"]["passive"] is True
        for case in result["cases"]
    )
    assert all(
        case["interaction"]["kernel_epsilon_m"] is None
        for case in result["cases"]
    )
