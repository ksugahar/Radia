import importlib.util
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs" / "stream_function"
EVIDENCE = Path(__file__).resolve().parent / "demos"
EXPECTED = {
    "demo_active_shield.json",
    "demo_fmm_biot_savart.json",
    "demo_pareto_results.json",
    "demo_reg_hyperparam_results.json",
    "demo_regcoil_fusion.json",
    "demo_regcoil_fusion_advanced.json",
    "demo_regularized_aca_results.json",
    "demo_shim_coil_purity.json",
    "pareto_cylinder.json",
    "pareto_cylinder_deform.json",
    "pareto_deform.json",
    "pareto_geometry_nsga.json",
    "pareto_tikhonov_aca.json",
}


def _load(name):
    return json.loads((EVIDENCE / name).read_text(encoding="utf-8"))


def _numbers(value):
    if isinstance(value, dict):
        for child in value.values():
            yield from _numbers(child)
    elif isinstance(value, list):
        for child in value:
            yield from _numbers(child)
    elif value is None:
        return
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        yield float(value)


def test_stream_function_demo_evidence_inventory_is_complete_and_finite():
    assert {path.name for path in EVIDENCE.glob("*.json")} == EXPECTED
    for name in EXPECTED:
        assert all(math.isfinite(value) for value in _numbers(_load(name))), name


def test_stream_function_demo_evidence_preserves_physical_claims():
    shield = _load("demo_active_shield.json")
    assert shield["stray_shielded"] < shield["stray_unshielded"]
    assert shield["shielding_factor"] > 10.0

    fmm = _load("demo_fmm_biot_savart.json")
    assert fmm["rel_err_ext_A_vs_B"] < 0.1

    pareto = _load("demo_pareto_results.json")
    assert pareto["n_trials"] == len(pareto["all_trials"])
    assert pareto["n_pareto"] == len(pareto["pareto"]) > 1

    hyper = _load("demo_reg_hyperparam_results.json")
    assert hyper["best_chain_rms"] < hyper["baseline_uniform_sigma_chain_rms"]
    assert hyper["aca_rank"] > 0

    regularized = _load("demo_regularized_aca_results.json")
    assert regularized["aca_rank"] > 0
    assert len(regularized["results"]) >= 5

    fusion = _load("demo_regcoil_fusion.json")
    assert fusion["net_current"]["betti1_winding_surface"] == 2
    assert all(
        target["bn_residual_rel"] < 1.0e-6
        for target in fusion["producible_targets"].values()
    )

    advanced = _load("demo_regcoil_fusion_advanced.json")
    assert 0.9 < advanced["force_stress"]["ratio_mean"] < 1.1
    assert advanced["focus_standoff"]["complexity_monotonic_in_gap"] is True
    assert advanced["focus_shape"]["complexity_reduction_vs_circular"] > 0.2


def test_stream_function_pareto_records_have_nonempty_fronts():
    keys = {
        "pareto_cylinder.json": ("envelope",),
        "pareto_cylinder_deform.json": ("flat_front", "opt_front"),
        "pareto_deform.json": ("flat_front", "opt_front"),
        "pareto_geometry_nsga.json": ("nsga_front", "nsga_lower_envelope"),
        "pareto_tikhonov_aca.json": ("nondominated_h1", "nondominated_linf"),
    }
    for name, front_keys in keys.items():
        data = _load(name)
        assert all(data[key] for key in front_keys), name


def test_docs_default_json_routes_to_validation_and_custom_output_is_preserved(tmp_path):
    spec = importlib.util.spec_from_file_location(
        "stream_function_validation_output", DOCS / "_validation_output.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    assert module.validation_output("record.json", DOCS) == EVIDENCE / "record.json"
    assert module.validation_json_for_basename(
        DOCS / "record", "record.json"
    ) == EVIDENCE / "record.json"
    assert module.validation_output("record.json", tmp_path) == tmp_path / "record.json"
    assert module.validation_json_for_basename(
        tmp_path / "custom", "record.json"
    ) == tmp_path / "custom.json"
