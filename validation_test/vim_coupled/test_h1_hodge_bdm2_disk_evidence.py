import json
from pathlib import Path


RESULT = Path(__file__).with_name("results_h1_hodge_bdm2_disk.json")


def _result() -> dict:
    return json.loads(RESULT.read_text(encoding="utf-8"))


def test_h1_hodge_disk_artifact_is_versioned_and_reproducible():
    result = _result()

    assert result["schema"] == "radia.validation.h1-hodge-bdm2-axisymmetric-disk.v1"
    assert result["created_at_utc"]
    assert result["tool_versions"]["radia"] != "unknown"
    assert len(result["tool_versions"]["radia_source_head"]) == 40
    assert result["tool_versions"]["radia_source_head"] == result["tool_versions"][
        "radia_source_head_end"
    ]
    assert all(len(value) == 64 for value in result["source_fingerprints"].values())
    assert result["identity"]["meshes_tracked"] is False
    assert RESULT.with_name("validate_h1_hodge_bdm2_disk.py").is_file()


def test_axisymmetric_q2_static_reference_is_fine_and_real():
    reference = _result()["axisymmetric_q2_static_reference"]

    assert reference["frequency_hz"] == 0.0
    assert reference["elements"] == 18432
    assert reference["dofs"] == 74305
    assert reference["normalized_Bz"][0] > 1.0
    assert abs(reference["normalized_Bz"][1]) < 1.0e-14


def test_bdm2_h_and_h1_p_ladder_reaches_one_percent():
    result = _result()
    rows = result["h1_hodge_3d_cases"]
    errors = [row["reference_relative_error"] for row in rows]

    assert [(row["maxh_m"], row["hdiv_order"], row["h1_order"]) for row in rows] == [
        (0.002, 1, 2),
        (0.002, 2, 3),
        (0.001, 2, 3),
        (0.001, 2, 4),
    ]
    assert errors[1:] == sorted(errors[1:], reverse=True)
    assert errors[1] < errors[0]
    assert errors[-1] < 0.01
    assert max(row["max_snapshot_relative_residual"] for row in rows) < 1.0e-8
    assert all(
        -1.0e-8 <= row["reduced_demag_generalized_eigenvalue"] <= 1.0 + 1.0e-5
        for row in rows
    )
    assert all(row["operator"]["unit_stiffness"] is True for row in rows)
    assert result["pass"] is True


def test_h1_hodge_disk_claim_is_bounded():
    claims = _result()["claim_boundary"]

    assert "below one-percent error" in claims["established"]
    assert "exact open-boundary H1 formulation" in claims["not_established"]
    assert "transient mixed-solver accuracy" in claims["not_established"]
