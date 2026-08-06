import json
from pathlib import Path


RESULT = Path(__file__).with_name(
    "results_mapped_hex_bdm2_hodge_reference.json"
)


def _result() -> dict:
    return json.loads(RESULT.read_text(encoding="utf-8"))


def test_hodge_reference_is_versioned_and_reproducible():
    result = _result()

    assert result["schema"] == "radia.validation.mapped-hex-bdm2-hodge-reference.v2"
    assert result["created_at_utc"]
    assert result["tool_versions"]["radia"] != "unknown"
    assert len(result["tool_versions"]["radia_source_head"]) == 40
    assert result["tool_versions"]["radia_source_head"] == result["tool_versions"][
        "radia_source_head_end"
    ]
    assert all(len(digest) == 64 for digest in result["source_fingerprints"].values())
    assert result["identity"]["body_active_dofs"] == 207
    assert result["identity"]["meshes_tracked"] is False
    assert RESULT.with_name("validate_mapped_hex_bdm2_hodge_reference.py").is_file()


def test_h1_hodge_projection_bounds_the_same_mapped_bdm2_space():
    result = _result()
    rows = result["h1_hodge_reference"]

    assert [row["h1_order"] for row in rows] == [2, 3]
    assert {row["hdiv_active_dofs"] for row in rows} == {207}
    assert all(row["spectrum"]["minimum"] >= -1.0e-8 for row in rows)
    assert all(row["spectrum"]["maximum"] <= 1.0 + 1.0e-5 for row in rows)
    assert all(row["spectrum"]["outside_count"] == 0 for row in rows)
    assert result["checks"]["h1_hodge_spectra_are_contractions"] is True
    assert all(row["backend"]["operator"] == "C.T @ K^-1 @ C" for row in rows)
    assert all(row["backend"]["unit_stiffness"] is True for row in rows)
    assert all(
        row["backend"]["contraction_contract"] == "standard-unit-H1 metric"
        for row in rows
    )
    assert result["checks"]["h1_hodge_public_operator_path_verified"] is True


def test_charge_bem_violation_is_localized_without_overclaiming_h1_accuracy():
    result = _result()
    charge = result["charge_bem_diagnostic"]

    assert charge["affinity"]["nonaffine_cell_count"] == 2
    assert charge["hdiv_dofs"] == 207
    assert charge["spectrum"]["maximum"] > 1.1
    assert charge["spectrum"]["outside_count"] > 0
    assert result["checks"]["charge_bem_spectrum_violation_reproduced"] is True
    assert result["pass"] is True
    assert "not_established" in result["claim_boundary"]
    assert "open-boundary H1 accuracy reference" in result["claim_boundary"][
        "not_established"
    ]


def test_mapped_bdm2_h1_response_feeds_the_shared_hcurl_mixed_solver():
    result = _result()
    mixed = result["h1_hcurl_mixed_reference"]

    assert mixed["mesh_elements"] == 54
    assert mixed["hdiv_active_dofs"] == 207
    assert mixed["hdiv_modes"] == 2
    assert mixed["demag_backend"] == "H1HodgeDemagOperator"
    assert mixed["snapshot_backend"] == "ngsolve-mass-preconditioned-cg"
    assert mixed["max_snapshot_relative_residual"] < 1.0e-8
    assert mixed["reduced_demag_generalized_spectrum"]["minimum"] >= -1.0e-8
    assert mixed["reduced_demag_generalized_spectrum"]["maximum"] <= 1.0 + 1.0e-5
    assert mixed["hcurl_interaction_backend"] == "hacapk-sampled-laplace"
    assert mixed["mixed_relative_residual"] < 1.0e-10
    assert mixed["average_joule_loss"] > 0.0
    assert mixed["magnetization_coefficient_norm"] > 0.0
    assert mixed["eddy_coefficient_norm"] > 0.0
    assert result["checks"]["mapped_bdm2_response_uses_generic_ngsolve_cg"] is True
    assert result["checks"][
        "mapped_bdm2_hcurl_mixed_solve_is_finite_and_converged"
    ] is True
