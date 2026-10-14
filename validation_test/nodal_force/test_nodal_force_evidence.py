import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / "nodal_force_results.json"


def test_nodal_force_evidence_matches_declared_accuracy():
    data = json.loads(EVIDENCE.read_text(encoding="utf-8"))
    finest = data["results"][-1]

    derived_checks = {
        "finest_mst_error_lt_0p02_pct": abs(finest["err_MST_pct"]) < 0.02,
        "finest_nodal_error_lt_0p02_pct": abs(finest["err_nodal_pct"]) < 0.02,
        "finest_mst_nodal_diff_lt_0p01_pct": (
            data["mst_nodal_diff_pct_finest"] < 0.01
        ),
    }

    assert data["schema"] == "radia.validation.nodal_force.v1"
    assert data["source_notebook"] == "docs/nodal_force/nodal_force.ipynb"
    assert data["checks"] == derived_checks
    assert all(derived_checks.values())
    assert data["F_analytical_N_per_m"] < 0.0
    assert data["F_dipole_N_per_m"] < 0.0
    assert [row["maxh"] for row in data["results"]] == sorted(
        (row["maxh"] for row in data["results"]), reverse=True
    )
