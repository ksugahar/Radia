import copy
import importlib.util
import json
from pathlib import Path


PATH = Path(__file__).resolve().parents[1] / "validation_test/omega_quadrature/assess_plateau.py"
SPEC = importlib.util.spec_from_file_location("omega_plateau", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


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
