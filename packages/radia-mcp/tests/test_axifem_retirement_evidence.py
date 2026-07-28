from __future__ import annotations

import copy

from radia_mcp.fem.axifem_retirement import FAMILIES, _sha, validate_axifem_element_evidence


def _evidence() -> dict:
    rows = []
    for index, (family, (cell_type, order, curved)) in enumerate(FAMILIES.items()):
        rows.append(
            {
                "family": family,
                "cell_type": cell_type,
                "order": order,
                "curved_geometry": curved,
                "vol_roundtrip": True,
                "vol_sha256": format(index + 1, "064x"),
                "mesh_contract_sha256": format(index + 11, "064x"),
                "identity_error_l2_sq": {
                    "interpolation": 1.0e-20,
                    "gradient": 2.0e-20,
                    "field": 3.0e-20,
                },
            }
        )
    evidence = {
        "schema": "radia.axifem-element-evidence.v2",
        "executed_at_utc": "2026-07-28T17:00:00Z",
        "execution_version": {"radia": "4.95.24", "ngsolve": "6.2.2604", "python": "3.12.10"},
        "git_head": "1" * 40,
        "git_dirty": False,
        "mesh_route": "Netgen .vol -> ngsolve.Mesh(path)",
        "test_summary": {"passed": 45, "failed": 0},
        "elements": rows,
        "p2_curved_metrics": {
            "straight_volume_error_percent": 3.897,
            "curved_volume_error_percent": 0.0248,
            "straight_total_flux_error_percent": 5.0744,
            "curved_total_flux_error_percent": 1.074,
        },
        "q2_curved_metrics": {
            "maximum_straight_equivalence_relative_error": 9.61e-12,
            "annular_successive_changes": [0.055019, 0.025835, 0.00606],
        },
    }
    evidence["evidence_payload_sha256"] = _sha(evidence)
    return evidence


def test_accepts_all_six_current_vol_backed_element_paths() -> None:
    result = validate_axifem_element_evidence(_evidence())
    assert result["status"] == "accepted"
    assert result["pass"] is True
    assert set(result["mesh_contract_sha256_by_family"]) == set(FAMILIES)


def test_rejects_missing_family_dirty_tree_and_tampered_payload() -> None:
    evidence = _evidence()
    evidence["elements"].pop()
    evidence["git_dirty"] = True
    result = validate_axifem_element_evidence(evidence)
    assert result["status"] == "rejected"
    assert result["checks"]["families_exact"] is False
    assert result["checks"]["worktree_clean"] is False
    assert result["checks"]["evidence_payload_sha256"] is False


def test_rejects_curved_geometry_and_convergence_overclaims() -> None:
    evidence = _evidence()
    bad = next(row for row in evidence["elements"] if row["family"] == "Q2_curved")
    bad["curved_geometry"] = False
    evidence["p2_curved_metrics"]["curved_volume_error_percent"] = 1.0
    evidence["q2_curved_metrics"]["annular_successive_changes"] = [0.01, 0.02, 0.03]
    evidence["evidence_payload_sha256"] = _sha({k: v for k, v in evidence.items() if k != "evidence_payload_sha256"})
    result = validate_axifem_element_evidence(evidence)
    assert result["status"] == "rejected"
    assert result["checks"]["Q2_curved_contract"] is False
    assert result["checks"]["p2_curved_geometry_improves_tenfold"] is False
    assert result["checks"]["q2_curved_annular_converges"] is False
