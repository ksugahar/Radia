"""Unit tests for the ESRF #3 Cubit three-engine mesh contract."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
BUILDER_PATH = (
    ROOT / "validation_test" / "esrf_three_engine" / "build_hybrid_undulator_mesh.py"
)
RUNNER_PATH = (
    ROOT / "validation_test" / "esrf_three_engine" / "run_hybrid_undulator_three_engine.py"
)
PROFILE_PATH = (
    ROOT / "validation_test" / "esrf_three_engine" / "profile_hybrid_undulator_hdiv_setup.py"
)
COIL_YOKE_BUILDER_PATH = (
    ROOT / "validation_test" / "esrf_three_engine" / "build_esrf_coil_yoke_kelvin_mesh.py"
)
COIL_YOKE_RUNNER_PATH = (
    ROOT / "validation_test" / "esrf_three_engine" / "run_coil_yoke_three_engine.py"
)
COIL_YOKE_CONTRACT_PATH = (
    ROOT / "validation_test" / "esrf_three_engine" / "esrf_coil_yoke.py"
)


def _builder_module():
    spec = importlib.util.spec_from_file_location("_esrf_three_engine_builder", BUILDER_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _runner_module():
    spec = importlib.util.spec_from_file_location("_esrf_three_engine_runner", RUNNER_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _coil_yoke_builder_module():
    spec = importlib.util.spec_from_file_location(
        "_esrf_coil_yoke_builder", COIL_YOKE_BUILDER_PATH
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _coil_yoke_contract_module():
    spec = importlib.util.spec_from_file_location(
        "_esrf_coil_yoke_contract", COIL_YOKE_CONTRACT_PATH
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_hybrid_undulator_three_engine_mesh_keeps_pm_as_given_source():
    builder = _builder_module()
    contract = builder.three_engine_material_contract()
    assert contract["response_materials"] == ("iron", "air", "kelvin")
    assert contract["pm_fem_material"] == "air"
    assert contract["pm_source"] == "fixed-given MagnetizationSource"
    assert contract["mixed_h1_source_potential"] == "global_physical"
    assert set(contract["required_boundaries"]) == {
        "iron_air_interface", "kelvin_int", "kelvin_ext"
    }


def test_hybrid_reflection_signature_preserves_material_and_parity():
    builder = _builder_module()
    upper = np.asarray(
        [
            [-0.01, -0.02, 0.03],
            [0.01, -0.02, 0.03],
            [0.00, 0.02, 0.03],
            [0.00, 0.00, 0.06],
        ]
    )
    lower = upper.copy()
    lower[:, 2] *= -1.0
    direct = builder._reflection_element_signature("iron", upper, reflected=False)
    reflected = builder._reflection_element_signature("iron", lower, reflected=True)
    assert direct == reflected
    assert len(direct) == len("iron") + 1 + upper.size * 8
    assert direct != builder._reflection_element_signature("air", upper, reflected=False)


def test_hybrid_undulator_journal_preserves_cubit_and_kelvin_contract(tmp_path):
    builder = _builder_module()
    assets = tmp_path / "assets"
    assets.mkdir()
    (assets / "iron.step").write_text("placeholder", encoding="utf-8")
    (assets / "magnet.step").write_text("placeholder", encoding="utf-8")
    journal, vol = builder.write_journal(
        tmp_path / "mesh",
        assets_dir=assets,
        kelvin_radius_m=0.18,
        iron_size_m=0.004,
        air_size_m=0.006,
        kelvin_size_m=0.012,
        curve_order=2,
    )
    assert journal.is_file()
    assert vol.name == "hybrid_undulator_kelvin.vol"
    text = journal.read_text(encoding="utf-8")
    assert "build_hybrid_undulator_inside_cubit.py" in text
    assert "export netgen" in text
    assert "order 2" in text
    inside = journal.with_name("build_hybrid_undulator_inside_cubit.py")
    inside_text = inside.read_text(encoding="utf-8")
    assert "build_inside_cubit" in inside_text
    assert "iron.step" in inside_text
    assert "magnet.step" in inside_text
    assert "kelvin_radius_m=0.17999999999999999" in inside_text


def test_hybrid_undulator_runner_keeps_fixed_magnetization_out_of_material_unknowns():
    runner = _runner_module()
    assert runner.MIXED_DOMAIN.reduced_materials == ("air",)
    assert runner.MIXED_DOMAIN.total_materials == ("iron", "kelvin")
    assert runner.MIXED_DOMAIN.nonlinear_materials == ("iron",)
    assert runner.MIXED_DOMAIN_LABEL == "H1 TOSCA mixed total/reduced Omega"
    points = runner.observation_points()
    assert points.shape == (27, 3)
    assert (points[:, 0] == 0.0).all()
    assert (points[:, 2] == 0.0).all()
    source_points = runner.source_accuracy_points()
    assert source_points.shape == (79, 3)
    assert np.max(np.abs(source_points[:, 0])) <= 0.015
    assert np.max(np.abs(source_points[:, 2])) <= 0.003
    # A source mesh is deliberately independent of the iron response mesh.
    # Any override is valid only with a retained direct-field reference; this
    # prevents a performance experiment from silently changing the PM source.
    text = RUNNER_PATH.read_text(encoding="utf-8")
    assert "--source-mesh" in text
    assert "--source-reference-mesh" in text
    assert "_source_mesh_direct_accuracy" in text
    assert "source_mesh_overridden" in text
    assert 'default="picard-mass-riesz"' in text
    assert "nonlinear_solver=nonlinear_solver" in text
    assert "--hdiv-gram-backend" in text
    assert "--hdiv-exact-dense-memory-mb" in text
    assert "gram_backend=gram_backend" in text


def test_hybrid_undulator_hdiv_profile_separates_source_from_gram_setup():
    text = PROFILE_PATH.read_text(encoding="utf-8")
    assert "source_present\": False" in text
    assert "source_present\": True" in text
    assert "source-projection" in text
    assert "source-field-compare" in text
    assert "--source-field-algorithm" in text
    assert "charge_gram_wall_s" in text
    assert "vim.Solve(" in text
    assert "--gram-backend" in text
    assert "--exact-dense-memory-mb" in text


def test_coil_yoke_journal_uses_kelvin_not_a_finite_outer_air_box(tmp_path):
    builder = _coil_yoke_builder_module()
    assets = tmp_path / "assets"
    assets.mkdir()
    (assets / "iron.step").write_text("placeholder", encoding="utf-8")
    journal, vol = builder.write_journal(
        tmp_path / "mesh",
        assets_dir=assets,
        kelvin_radius_m=0.16,
        iron_size_m=0.006,
        air_size_m=0.016,
        kelvin_size_m=0.035,
        gap_size_m=0.003,
        gap_refinement_radius_m=0.026,
        gap_refinement_half_length_m=0.032,
        beam_axis=0,
        curve_order=2,
    )
    contract = builder.three_engine_material_contract()
    assert contract["response_materials"] == ("iron", "air", "kelvin")
    assert contract["coil_source"] == "shared mesh-free CoilBuilder solid current"
    assert contract["mixed_h1_source_potential"] == (
        "total_volume_hodge_with_harmonic_cut"
    )
    assert contract["finite_outer_air_box_forbidden"] is True
    assert vol.name == "coil_yoke_kelvin.vol"
    assert 'export netgen' in journal.read_text(encoding="utf-8")
    inside = journal.with_name("build_esrf_coil_yoke_inside_cubit.py")
    text = inside.read_text(encoding="utf-8")
    assert "build_inside_cubit" in text
    assert 'symmetry=["z"]' in COIL_YOKE_BUILDER_PATH.read_text(encoding="utf-8")
    assert "iron.step" in text


def test_coil_yoke_runner_keeps_one_coil_source_and_three_required_engines():
    text = COIL_YOKE_RUNNER_PATH.read_text(encoding="utf-8")
    source_text = (
        ROOT / "validation_test" / "esrf_three_engine" / "esrf_coil_yoke.py"
    ).read_text(encoding="utf-8")
    assert "build_radia_coil_source" in text
    assert "source_is_meshed\": False" in source_text
    assert "solve_hdiv" in text
    assert "solve_reduced_a" in text
    assert "solve_omega" in text
    assert 'source_potential_contract="total_hodge"' in (
        ROOT / "validation_test" / "c_type_three_engine" / "run_three_engine.py"
    ).read_text(encoding="utf-8")
    assert "require_static_electromagnet_three_engine_contract" in text
    assert "finite_outer_air_box_forbidden\": True" in text
    assert "fem_kelvin_mesh_shared\": True" in text
    assert "fem_periodic_kelvin_mesh_shared" not in text
    assert "--resume" in text
    assert "--mixed-nonlinear-maximum-iterations" in text
    assert "nonlinear_maximum_iterations=mixed_nonlinear_maximum_iterations" in text
    assert "--hdiv-image" in text
    assert "image=options.hdiv_image" in text
    assert "hdiv_image=options.hdiv_image" in text
    assert "observation_volume_quadrature" in text
    assert "average_observation_field" in text
    assert '"tensor_gauss_2x2x2"' in text
    assert '"source_potential_contract": "total_hodge"' in text
    assert '"linear_solver": options.reduced_a_solver' in text
    shared = (
        ROOT / "validation_test" / "c_type_three_engine" / "run_three_engine.py"
    ).read_text(encoding="utf-8")
    assert '"hmat_stats": dict(result.get("hmat_stats") or {})' in shared
    assert '"charge_gram_wall_s"' in shared


def test_coil_yoke_observation_volume_average_preserves_affine_fields():
    helper = _coil_yoke_contract_module()
    centres, samples = helper.observation_volume_quadrature(7, half_width_m=2.0e-5)
    assert samples.shape == (8 * len(centres), 3)
    affine = np.column_stack(
        (2.0 * samples[:, 0] - samples[:, 1], samples[:, 2], 1.0 + samples[:, 0])
    )
    averaged = helper.average_observation_field(affine, len(centres))
    expected = np.column_stack(
        (2.0 * centres[:, 0] - centres[:, 1], centres[:, 2], 1.0 + centres[:, 0])
    )
    np.testing.assert_allclose(averaged, expected, atol=1.0e-15, rtol=0.0)


def _coil_yoke_runner_module(monkeypatch):
    monkeypatch.syspath_prepend(str(COIL_YOKE_RUNNER_PATH.parent))
    spec = importlib.util.spec_from_file_location(
        "_esrf_coil_yoke_runner", COIL_YOKE_RUNNER_PATH
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _converged_diagnostics(converged=True):
    return {
        "formulation": "test",
        "nonlinear": True,
        "nonlinear_stats": {"converged": converged, "iterations": 7},
    }


def test_coil_yoke_runner_declares_picard_state_and_anderson_controls():
    text = COIL_YOKE_RUNNER_PATH.read_text(encoding="utf-8")
    assert "--mixed-anderson-depth" in text
    assert "--reduced-a-anderson-depth" in text
    assert "--mixed-relaxation" in text
    assert "MixedOmegaPicardNotConverged" in text
    assert 'CHECKPOINT_SCHEMA = "radia.validation.esrf-coil-yoke-checkpoint.v3"' in text
    assert 'STATE_SCHEMA = "radia.validation.esrf-coil-yoke-picard-state.v1"' in text
    assert "_write_state(" in text
    assert "_read_state(" in text
    # The iteration cap is provenance, never checkpoint identity.
    common = text.split("common = _checkpoint_contract(", 1)[1].split(")", 1)[0]
    assert "nonlinear_maximum_iterations" not in common
    assert '"nonlinear_maximum_iterations": iteration_caps' in text
    shared = (
        ROOT / "validation_test" / "c_type_three_engine" / "run_three_engine.py"
    ).read_text(encoding="utf-8")
    assert "anderson_depth=int(anderson_depth)" in shared
    assert "nonlinear_mu_r_initial=mu_r_initial" in shared
    assert "nu_initial=nu_initial" in shared


def test_coil_yoke_checkpoint_reader_accepts_only_converged_matching_results(tmp_path, monkeypatch):
    runner = _coil_yoke_runner_module(monkeypatch)
    contract = runner._checkpoint_contract(
        case=7, nonlinear_tolerance=2.0e-5, engine="reduced_a",
        engine_settings={"linear_solver": "direct", "relaxation": 0.1})
    field = np.arange(6.0).reshape(2, 3)
    path = tmp_path / "run.reduced_a.checkpoint.json"

    runner._write_checkpoint(path, contract, field, _converged_diagnostics(),
                             {"nonlinear_maximum_iterations": 80})
    stored = json.loads(path.read_text(encoding="utf-8"))
    assert stored["schema"] == runner.CHECKPOINT_SCHEMA
    assert stored["provenance"]["nonlinear_maximum_iterations"] == 80
    resumed_field, resumed_diagnostics = runner._read_checkpoint(path, contract)
    np.testing.assert_array_equal(resumed_field, field)
    assert resumed_diagnostics["nonlinear_stats"]["converged"] is True

    with pytest.raises(RuntimeError, match="contract changed"):
        runner._read_checkpoint(path, dict(contract, nonlinear_tolerance=1.0e-6))
    with pytest.raises(RuntimeError, match="non-converged"):
        runner._write_checkpoint(path, contract, field, _converged_diagnostics(False), {})

    # A v2 checkpoint carried the cap in its identity; a converged one still resumes.
    legacy = tmp_path / "legacy.reduced_a.checkpoint.json"
    legacy_contract = {**contract, "nonlinear_maximum_iterations": 80,
                       "engine_settings": {**contract["engine_settings"],
                                           "nonlinear_maximum_iterations": 80}}
    legacy.write_text(json.dumps({
        "schema": "radia.validation.esrf-coil-yoke-checkpoint.v2",
        "contract": legacy_contract,
        "field_T": field.tolist(),
        "diagnostics": _converged_diagnostics(),
    }), encoding="utf-8")
    resumed_field, _ = runner._read_checkpoint(legacy, contract)
    np.testing.assert_array_equal(resumed_field, field)
    # ... but a non-converged v2 result is refused instead of reused.
    legacy.write_text(json.dumps({
        "schema": "radia.validation.esrf-coil-yoke-checkpoint.v2",
        "contract": legacy_contract,
        "field_T": field.tolist(),
        "diagnostics": _converged_diagnostics(False),
    }), encoding="utf-8")
    with pytest.raises(RuntimeError, match="non-converged"):
        runner._read_checkpoint(legacy, contract)
    legacy.write_text(json.dumps({"schema": "unknown", "contract": contract,
                                  "field_T": [], "diagnostics": {}}), encoding="utf-8")
    with pytest.raises(RuntimeError, match="unsupported checkpoint schema"):
        runner._read_checkpoint(legacy, contract)
    assert runner._read_checkpoint(tmp_path / "missing.json", contract) is None


def test_coil_yoke_picard_state_roundtrip_is_explicitly_partial(tmp_path, monkeypatch):
    runner = _coil_yoke_runner_module(monkeypatch)
    contract = runner._checkpoint_contract(
        case=6, engine="mixed_total_reduced_omega",
        engine_settings={"relaxation": 0.3, "anderson_depth": 2})
    stats = {
        "iterations": 80, "converged": False, "relative_B_change": 9.03e-5,
        "contraction_rate_estimate": 0.97,
        "element_numbers": [3, 5, 8], "mu_r_elements": [900.0, 12.5, 1.0],
        "history": [{"iteration": 1}, {"iteration": 2, "relative_B_change": 0.5}],
    }
    state = runner._picard_state("mixed_total_reduced_omega", stats, "mu_r_elements")
    assert state["mu_r_elements"] == [900.0, 12.5, 1.0]
    assert state["element_numbers"] == [3, 5, 8]
    assert state["iterations"] == 80
    assert state["resumed_iterations"] == 0
    path = tmp_path / "run.mixed_total_reduced_omega.state.json"
    runner._write_state(path, contract, state, {"nonlinear_maximum_iterations": 80})
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["schema"] == runner.STATE_SCHEMA
    assert payload["converged"] is False
    assert runner._read_state(path, contract) == state
    with pytest.raises(RuntimeError, match="state contract changed"):
        runner._read_state(path, dict(contract, engine="reduced_a"))
    assert runner._read_state(tmp_path / "missing.json", contract) is None
