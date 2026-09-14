import copy
import importlib.util
import json
import hashlib
from pathlib import Path


PATH = Path(__file__).resolve().parents[1] / "validation_test/omega_quadrature/assess_plateau.py"
SPEC = importlib.util.spec_from_file_location("omega_plateau", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_ci_candidate_bonus12_16_replays_recorded_evidence():
    root = PATH.parents[1] / "esrf_three_engine/results/candidate_59b094d8"
    raw = root / "omega_algebraic_bonus12_16.json"
    report = root / "omega_plateau_bonus12_16.json"
    assert hashlib.sha256(raw.read_bytes()).hexdigest() == (
        "7cc8befa4eeb1f7ae056c816e2c36da5e07ce58bc36e99a0af7a65c8b343ca15")
    assert hashlib.sha256(report.read_bytes()).hexdigest() == (
        "fa1577ebf844f56e81ae68362756ae84e88301a00146fbe5b2b7506df6a99a18")
    payload = json.loads(raw.read_text(encoding="utf-8"))
    wheel = payload["runtime"]["direct_url"]["archive_info"]["hashes"]["sha256"]
    assert wheel == "46e3aa1ddd89010e21419e1d28f6e44e403cf95014d7e5486109eca4232145b9"
    result = MODULE.assess(payload, 12, 16)
    assert result == json.loads(report.read_text(encoding="utf-8"))
    assert result["passed"] is True


def payload(delta=1e-5):
    def row(bonus, value):
        return {
            "bonus": bonus,
            "order": 2,
            "source_hodge": {"bonus_intorder": bonus, "relative_harmonic_norm": 0.0025},
            "linear_residual": {"free_dofs": {"relative": 1e-12}},
            "block_action_residual": {"blocks": {
                name: {"action_relative": 1e-12}
                for name in ("phi_reduced", "phi_total", "interface_constraint")}},
            "field_observations": {"samples_m": [[0, 0, 0], [1, 0, 0]],
                                   "B_samples_T": [[1, 0, 0], [0, value, 0]]},
        }
    return {"completed": True, "source_unchanged": True, "mesh_unchanged": True,
            "runtime": {"mode": "wheel", "direct_url": {"archive_info": {
                "hashes": {"sha256": "a" * 64}}}},
            "rows": [row(8, delta), row(12, 0.0)]}


def test_formal_plateau_passes_only_with_matched_high_rule_evidence():
    result = MODULE.assess(payload())
    assert result["passed"] is True
    assert result["scope"].startswith("ESRF6 nominal order-2")
    assert json.loads(json.dumps(result, allow_nan=False))["passed"] is True


def test_formal_plateau_rejects_field_drift_and_hodge_mismatch():
    drifting = payload(delta=0.1)
    assert MODULE.assess(drifting)["passed"] is False
    mismatched = copy.deepcopy(payload())
    mismatched["rows"][1]["source_hodge"]["bonus_intorder"] = 8
    result = MODULE.assess(mismatched)
    assert result["passed"] is False
    assert result["checks"]["hodge_bonus_12"] is False
