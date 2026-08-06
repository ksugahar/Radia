import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[3]
RESULT = (
    REPO
    / "validation_test"
    / "vim_coupled"
    / "results_hcurl_eddy_bubble_disk.json"
)


def _load_result() -> dict:
    return json.loads(RESULT.read_text(encoding="utf-8"))


def test_disk_validation_evidence_records_versioned_passing_run():
    result = _load_result()

    assert result["schema"] == "cae-ai-lab.solver-run.v1"
    assert result["pass"] is True
    assert result["created_at_utc"]
    assert result["tool_versions"]["radia"] != "unknown"
    assert result["tool_versions"]["ngsolve"] != "unknown"
    assert result["checks"]["validation_passed"] is True


def test_disk_validation_evidence_freezes_h_and_p_positive_controls():
    result = _load_result()
    h_rows = result["positive_h_refinement"]
    p_rows = result["positive_p_refinement"]

    assert [row["maxh_m"] for row in h_rows] == [0.003, 0.002, 0.0015, 0.001]
    assert [row["parent"]["order"] for row in p_rows] == [1, 2, 3]
    assert max(row["modal"]["port_dominant_abs_error_pct"] for row in h_rows) < 2.0
    assert max(row["modal"]["port_dominant_abs_error_pct"] for row in p_rows) < 2.0
    assert all(
        row["modal"]["port_dominant_residue_fraction"] > 0.85
        for row in h_rows + p_rows
    )
    assert all(
        row["training"]["vector_potentials"]
        == ["A_uniform", "r2_A", "r4_A", "z2_A"]
        for row in h_rows + p_rows
    )


def test_disk_validation_evidence_promotes_single_port_counterfactual():
    result = _load_result()
    negative = result["single_port_negative_control"]

    assert [row["training"]["krylov_steps"] for row in negative] == [8, 16]
    assert all(row["training"]["port_count"] == 1 for row in negative)
    assert min(
        row["modal"]["port_dominant_abs_error_pct"] for row in negative
    ) > 4.0
    assert (
        result["errors"]["best_single_port_error_pct"]
        - result["errors"]["best_polynomial_error_pct"]
    ) > 3.0


def test_disk_validation_evidence_requires_projection_and_passivity():
    result = _load_result()
    rows = (
        result["positive_h_refinement"]
        + result["positive_p_refinement"]
        + result["single_port_negative_control"]
    )

    assert max(
        row["interaction"]["projection_relative_residual"] for row in rows
    ) <= 1.0e-9
    assert min(
        row["modal"]["minimum_resistance_eigenvalue"] for row in rows
    ) >= -1.0e-12
    assert min(
        row["modal"]["minimum_inductance_eigenvalue"] for row in rows
    ) >= -1.0e-12
