import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
RESULT = HERE / "team28_axisym_fem_height_sweep.json"


def test_full_fem_height_sweep_is_fresh_and_passes():
    payload = json.loads(RESULT.read_text(encoding="utf-8"))

    assert payload["schema"] == "radia.team28-axisym-fem-height-sweep.v1"
    assert payload["pass"] is True
    assert all(payload["checks"].values())
    for relative, digest in payload["source_sha256"].items():
        data = (REPO_ROOT / relative).read_bytes().replace(b"\r\n", b"\n")
        if hashlib.sha256(data).hexdigest() != digest:
            spec = importlib.util.spec_from_file_location("privacy_source", REPO_ROOT / "validation_test/radia_mcp/privacy_source_contract.py")
            contract = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(contract)
            assert contract.verifies_redaction(relative, digest, data.decode(), payload.get("source_privacy_redactions", {}).get(relative, {})), relative
    assert payload["dZ_mm"] == list(range(-7, 18))
    assert payload["max_abs_legacy_minus_lab_N"] < 1.0e-3
    assert 0.2 < payload["equilibrium_dZ_mm"] < 0.3
