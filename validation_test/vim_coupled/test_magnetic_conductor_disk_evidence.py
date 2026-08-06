import importlib.util
import json
from pathlib import Path

import pytest


RESULT = Path(__file__).with_name(
    "results_magnetic_conductor_disk_adjudication.json"
)


def _load_result() -> dict:
    return json.loads(RESULT.read_text(encoding="utf-8"))


def _load_validation_driver():
    path = RESULT.with_name("validate_magnetic_conductor_disk.py")
    spec = importlib.util.spec_from_file_location(
        "validate_magnetic_conductor_disk_contract", path
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_magnetic_conductor_disk_identity_and_reference_are_frozen():
    result = _load_result()
    identity = result["identity"]

    assert result["schema"] == "radia.validation.vim-coupled.magnetic-conductor-disk.v1"
    assert result["executed_at_utc"]
    assert result["execution_version"]["radia"] != "unknown"
    assert identity["geometry"] == "solid circular disk"
    assert identity["radius_m"] == 0.01
    assert identity["thickness_m"] == 0.0005
    assert identity["material"]["mu_r"] == 100.0
    assert identity["material"]["sigma_S_per_m"] == 1.0e7
    assert [row["frequency_hz"] for row in result["axisymmetric_q2_reference"]] == [
        100.0,
        10_000.0,
    ]
    assert all(
        row["mesh_parameters"]
        == {
            "nr_disk": 160,
            "nz_disk": 32,
            "nr_air": 32,
            "nz_air": 32,
            "outer_radius_m": 0.5,
            "outer_half_height_m": 0.5,
        }
        for row in result["axisymmetric_q2_reference"]
    )
    assert (RESULT.parent / "validate_magnetic_conductor_disk.py").is_file()
    assert (RESULT.parent / "build_magnetic_conductor_disk_hex.py").is_file()


def test_execution_identity_covers_direct_python_and_native_operators():
    fingerprints = _load_validation_driver()._source_fingerprints()

    assert {
        "src/radia/axifem.pyd",
        "src/radia/_radia_pybind.pyd",
        "src/radia/vim/_eddy_hybrid.py",
        "src/radia/vim/_hcurl_tet_interaction.py",
        "validation_test/vim_coupled/validate_magnetic_conductor_disk.py",
        "validation_test/vim_coupled/build_magnetic_conductor_disk_hex.py",
    } <= fingerprints.keys()
    assert all(len(digest) == 64 for digest in fingerprints.values())


def test_full_3d_hcurl_p_and_h_refinement_reach_two_percent():
    result = _load_result()
    rows = result["full_3d_hcurl_A_form"]
    p_rows = [
        row
        for row in rows
        if row["frequency_hz"] == 10_000.0 and row["disk_maxh_m"] == 0.002
    ]
    fine = next(
        row
        for row in rows
        if row["frequency_hz"] == 10_000.0
        and row["order"] == 3
        and row["disk_maxh_m"] == 0.001
    )

    assert [row["order"] for row in p_rows] == [1, 2, 3]
    assert [row["reference_relative_error"] for row in p_rows] == sorted(
        (row["reference_relative_error"] for row in p_rows), reverse=True
    )
    assert fine["reference_relative_error"] < 0.02
    assert max(row["relative_residual"] for row in rows) < 1.0e-8


def test_bdm1_hex_h_refines_but_bdm2_is_not_promoted():
    result = _load_result()
    static = result["hdiv_mmm_static_hex"]
    bdm1 = static["bdm1_h_ladder"]
    bdm2 = static["bdm2_negative_control"]
    axisymmetric = result["axisymmetric_q2_reference"][0]["normalized_Bz"]
    reference = complex(*axisymmetric)

    assert [row["size_mm"] for row in bdm1] == [2.0, 1.0, 0.5]
    assert [row["reference_relative_error"] for row in bdm1] == sorted(
        (row["reference_relative_error"] for row in bdm1), reverse=True
    )
    assert bdm1[-1]["reference_relative_error"] < 0.02
    assert all(
        row["reference_relative_error"]
        == pytest.approx(
            abs(complex(row["normalized_Bz"], 0.0) - reference) / abs(reference)
        )
        for row in bdm1
    )
    assert bdm2[1]["normalized_Bz"] > bdm2[0]["normalized_Bz"]
    assert static["bdm2_status"].startswith("not promoted")


def test_mapped_hex_bdm2_spectrum_gate_is_versioned_and_root_caused():
    result = _load_result()
    gate = result["mapped_hex_bdm2_spectrum_gate"]
    cases = {row["case"]: row for row in gate["minimal_generalized_spectra"]}
    disk = gate["full_disk_32_hex"]

    assert gate["executed_at_utc"]
    assert len(gate["radia_source_head"]) == 40
    assert gate["status"] == "material_solve_rejected_for_nonaffine_hex_bdm2"
    assert cases["two_affine_hex"]["outside_count"] == 0
    assert cases["two_warped_hex"]["maximum"] > 1.1
    assert cases["four_hex_annular_sector"]["minimum"] < -0.01
    assert cases["four_hex_annular_sector"]["maximum"] > 1.7
    assert cases["two_warped_hex_high_near"]["outside_count"] == 0
    assert cases["four_hex_annular_sector_high_near"]["outside_count"] > 0
    assert disk["aca_hmatrix"]["minimum"] < -5.0
    assert disk["all_dense_gram"]["minimum"] < -5.0
    assert disk["aca_dense_maximum_difference"] < 3.0e-4
    assert disk["aca_dense_minimum_difference"] < 1.0e-5
    assert disk["volume_only"]["minimum"] >= 0.0
    assert disk["boundary_only"]["minimum"] >= 0.0
    assert result["checks"]["mapped_hex_bdm2_material_gate_verified"] is True
    assert result["checks"]["mapped_hex_bdm2_space_hodge_contraction_verified"] is True


def test_validation_driver_records_expected_mapped_hex_bdm2_rejection(monkeypatch):
    driver = _load_validation_driver()

    def reject(_mesh_path, *, order):
        assert order == 2
        raise NotImplementedError(
            "vim.Solve: mapped/non-affine HEX BDM2 material solve is gated because "
            "the current separate volume/surface charge quadrature is unsafe; use "
            "mapped HEX BDM1 instead"
        )

    monkeypatch.setattr(driver, "solve_hdiv_static", reject)
    row = driver.probe_mapped_hex_bdm2_material_gate(Path("unused.vol"))

    assert row["status"] == "rejected_as_expected"
    assert row["expected_gate"] is True
    assert row["error_type"] == "NotImplementedError"
    assert "mapped/non-affine HEX BDM2" in row["message"]


def test_live_runner_replays_independent_references_and_bdm2_gate():
    result = _load_result()
    replay = result["mapped_hex_bdm2_gate_replay"]

    assert replay["radia"] == "4.95.42"
    assert len(replay["radia_source_head"]) == 40
    assert replay["source_fingerprints_stable_during_run"] is True
    assert replay["full_3d_hcurl"]["reference_relative_error"] < 0.02
    assert replay["mapped_hex_bdm1"]["role"].startswith("coarse positive")
    assert replay["coupled_hdiv_mmm_hcurl_eddy_bubble"][
        "reference_relative_error"
    ] > 0.02
    assert replay["coupled_hdiv_mmm_hcurl_eddy_bubble"][
        "relative_residual"
    ] < 1.0e-10
    assert replay["direct_sampled_parity"]["observable_relative_difference"] < 1.0e-4
    assert replay["direct_sampled_parity"]["joule_relative_difference"] < 1.0e-4
    assert {row["status"] for row in replay["mapped_hex_bdm2_gate"]} == {
        "rejected_as_expected"
    }
    assert replay["pass"] is True
    assert result["checks"]["mapped_hex_bdm2_live_gate_replay_verified"] is True


def test_coupled_smoke_is_bounded_to_mechanics_and_routing():
    result = _load_result()
    coupled = result["coupled_hdiv_mmm_hcurl_eddy_bubble_smoke"]
    checks = result["checks"]

    assert coupled["interaction"].endswith("not the accuracy oracle")
    assert coupled["relative_residual"] < 1.0e-10
    assert coupled["joule_loss"] >= 0.0
    assert coupled["stale_10khz_route_rejected"] is True
    assert checks["production_direct_q2_hex_execution_verified"] is True
    assert checks["production_direct_q2_hex_h_convergence_verified"] is True
    assert "universal speed or accuracy superiority over another solver" in result[
        "adjudication"
    ]["not_established"]


def test_sampled_coupled_h_ladder_reaches_two_percent():
    result = _load_result()
    rows = result["coupled_sampled_h_ladder"]["rows"]
    errors = [row["reference_relative_error"] for row in rows]

    assert [row["size_mm"] for row in rows] == [4.0, 2.0, 1.0]
    assert [row["hexes"] for row in rows] == [32, 96, 384]
    assert errors == sorted(errors, reverse=True)
    assert errors[-1] < 0.02
    assert max(row["relative_residual"] for row in rows) < 1.0e-10
    assert min(row["joule_loss"] for row in rows) >= 0.0
    assert result["checks"]["coupled_sampled_h_fine_error_below_2pct"] is True
    assert result["checks"]["production_direct_q2_hex_h_convergence_verified"] is True
    fidelity = result["coupled_sampled_h_ladder"]["reference_fidelity_gate"]
    assert fidelity["accepted_fine_q2_elements"] == 18432
    assert fidelity["accepted_fine_error"] < 0.02
    assert fidelity["rejected_coarse_q2_elements"] == 6656
    assert fidelity["rejected_coarse_error"] > 0.02


def test_production_direct_q2_hex_replay_records_backend_and_parity():
    result = _load_result()
    replay = result["production_direct_q2_hex_replay"]
    checks = result["checks"]

    assert replay["radia"] == "4.95.42"
    assert replay["source_head_stable_during_run"] is True
    assert replay["source_fingerprints_stable_during_run"] is True
    assert replay["solver_run_contract_validated"] is True
    assert replay["full_profile_passed"] is True
    assert replay["hcurl_diagonal_backend"] == "direct-q2-hex-reference-density"
    assert replay["reference_density"] == "curl(T)*abs(det(dX/dxi))"
    assert replay["tensor_degree"] == 2
    assert replay["hexes"] == 32
    assert replay["charge_count"] == 32 * 27
    assert replay["affine_geometry_residual"] > 0.4
    assert replay["relative_residual"] < 1.0e-10
    assert replay["joule_loss"] >= 0.0
    assert replay["sampled_observable_relative_difference"] < 1.0e-4
    assert replay["sampled_joule_relative_difference"] < 1.0e-4
    assert checks["production_direct_q2_hex_backend_observed"] is True
    assert checks["production_direct_q2_hex_sampled_observable_parity_below_1e-4"]
    assert checks["production_direct_q2_hex_sampled_joule_parity_below_1e-4"]


def test_production_direct_q2_hex_h_ladder_is_strict_and_convergent():
    result = _load_result()
    ladder = result["production_direct_q2_hex_h_ladder"]
    rows = ladder["rows"]
    errors = [row["reference_relative_error"] for row in rows]

    assert ladder["solver_run_contract_validated"] is True
    assert ladder["direct_h_ladder_passed"] is True
    assert ladder["source_head_stable_during_run"] is True
    assert ladder["source_fingerprints_stable_during_run"] is True
    assert "src/radia/vim/_hcurl_tet_interaction.py" in ladder[
        "source_fingerprints"
    ]
    assert "src/radia/_radia_pybind.pyd" in ladder["source_fingerprints"]
    assert [row["hexes"] for row in rows] == [32, 96, 384]
    assert errors == sorted(errors, reverse=True)
    assert errors[-1] < 0.02
    assert max(row["sampled_observable_relative_difference"] for row in rows) < 1.0e-4
    assert max(row["sampled_joule_relative_difference"] for row in rows) < 1.0e-4
    assert max(row["relative_residual"] for row in rows) < 1.0e-10
    assert min(row["joule_loss"] for row in rows) >= 0.0
    assert all(row["charge_count"] < row["unpruned_charge_count"] for row in rows)
    assert ladder["timing"]["observed_speedup"] > 5.0
    assert result["checks"]["production_direct_q2_hex_h_convergence_verified"]


def test_full_profile_replay_is_versioned_and_source_stable():
    replay = _load_result()["full_profile_replay"]

    assert replay["executed_at_utc"]
    assert replay["radia"] == "4.95.41"
    assert len(replay["radia_source_head"]) == 40
    assert replay["source_head_stable_during_run"] is True
    assert replay["source_fingerprints_stable_during_run"] is True
    assert replay["solver_run_contract_validated"] is True
    assert replay["full_profile_passed"] is True
    assert replay["max_axisymmetric_component_difference_from_stored"] == 0.0
    assert replay["max_hcurl_reference_error_difference_from_stored"] == 0.0
    assert replay["max_bdm1_normalized_Bz_difference_from_stored"] < 1.0e-14
