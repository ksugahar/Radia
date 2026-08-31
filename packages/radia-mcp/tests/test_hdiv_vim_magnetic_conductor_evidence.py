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
    / "results_magnetic_conductor_disk_adjudication.json"
)
HODGE_RESULT = RESULT.with_name("results_mapped_hex_bdm2_hodge_reference.json")


def test_mcp_teaches_only_the_verified_magnetic_conductor_claims():
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    hodge = json.loads(HODGE_RESULT.read_text(encoding="utf-8"))
    coupled = " ".join(get_hdiv_vim_documentation("eddy_bubble").split())

    assert result["checks"]["bdm1_hex_fine_error_below_2pct"] is True
    assert result["checks"]["bdm2_hex_nonmonotone_gap_detected"] is True
    assert result["checks"]["coupled_sampled_h_fine_error_below_2pct"] is True
    assert result["checks"]["coupled_h_ladder_fine_reference_required"] is True
    assert result["checks"]["production_direct_q2_hex_execution_verified"] is True
    assert result["checks"]["production_direct_q2_hex_backend_observed"] is True
    assert result["checks"]["production_direct_q2_hex_h_convergence_verified"] is True
    assert result["checks"]["mapped_hex_bdm2_spectrum_violation_reproduced"] is True
    assert result["checks"]["mapped_hex_bdm2_dense_aca_parity_verified"] is True
    assert result["checks"]["mapped_hex_bdm2_material_gate_verified"] is True
    assert result["full_profile_replay"]["full_profile_passed"] is True
    assert result["full_profile_replay"]["source_head_stable_during_run"] is True
    assert hodge["checks"]["h1_hodge_spectra_are_contractions"] is True
    assert hodge["checks"]["charge_hacapk_spectrum_is_a_contraction"] is True
    assert hodge["checks"][
        "charge_hacapk_material_solve_and_field_are_production_finite"
    ] is True
    assert "mapped-HEX BDM1" in coupled
    assert "BDM2" in coupled
    assert "not an accuracy oracle" in coupled
    assert "direct-Q2" in coupled
    assert "hex_geometry_backend" in coupled
    assert "8.03e-7" in coupled
    assert "strict 32/96/384-HEX h ladder" in coupled
    assert "Nsample x Nsample" in coupled
    assert "1819.7 s to 339.8 s" in coupled
    assert "mapped/non-affine pure-HEX BDM2 primal material and field lane" in coupled
    assert "ACA compression was not the cause" in coupled
    assert "production C++ composite operator" in coupled
    assert "Shape derivatives remain fail-loud" in coupled
    assert "same 207 active mapped-body BDM2 DoFs" in coupled
    assert "0.9961" in coupled
    assert "q9/q12" in coupled
    assert "Cubit 2025.12" in coupled
    assert "1.06e-10" in coupled
    assert "not an open-boundary accuracy oracle" in coupled
    assert "universal solver-superiority claim" in coupled
