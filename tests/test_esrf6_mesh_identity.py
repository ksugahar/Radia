"""Fast input-selection regression; no Radia native import or solver."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "case6_inputs", ROOT / "validation_test/esrf_three_engine/check_case6_inputs.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


@pytest.fixture
def inputs(tmp_path):
    iron, fem = tmp_path / "iron.vol", tmp_path / "fem.vol"
    iron.write_bytes(b"selected iron")
    fem.write_bytes(b"selected FEM")
    contract = tmp_path / "contract.json"
    contract.write_text(json.dumps({
        "schema": "radia.validation.case6-mesh-contract.v1", "case": 6,
        "iron": {"sha256": hashlib.sha256(iron.read_bytes()).hexdigest()},
        "fem": {"sha256": hashlib.sha256(fem.read_bytes()).hexdigest()},
    }))
    return iron, fem, contract


def test_matching_input_remains_numerically_unaccepted(inputs):
    result = MODULE.verify_inputs(*inputs)
    assert result["identity_passed"] is True
    assert result["numerical_acceptance"].startswith("HOLD")


@pytest.mark.parametrize("index", [0, 1])
def test_old_or_modified_mesh_is_rejected(inputs, index):
    inputs[index].write_bytes(b"obsolete unmerged mesh")
    with pytest.raises(ValueError, match="identity mismatch"):
        MODULE.verify_inputs(*inputs)


def test_swapped_roles_are_rejected(inputs):
    iron, fem, contract = inputs
    with pytest.raises(ValueError, match="iron mesh identity mismatch"):
        MODULE.verify_inputs(fem, iron, contract)


def test_missing_file_is_not_a_fallback(inputs):
    inputs[0].unlink()
    with pytest.raises(FileNotFoundError):
        MODULE.verify_inputs(*inputs)


def test_bad_contract_is_rejected(inputs):
    inputs[2].write_text('{"schema":"unknown","case":6}')
    with pytest.raises(ValueError, match="not an ESRF6"):
        MODULE.verify_inputs(*inputs)


def test_recorded_contract_is_the_bytes_actually_validated(inputs, monkeypatch):
    original_bytes = inputs[2].read_bytes()
    original_sha256 = MODULE.sha256

    def change_contract_after_read(path):
        inputs[2].write_text('{"schema":"replacement"}')
        return original_sha256(path)

    monkeypatch.setattr(MODULE, "sha256", change_contract_after_read)
    result = MODULE.verify_inputs(*inputs)
    assert result["contract_sha256"] == hashlib.sha256(original_bytes).hexdigest()
    assert result["contract_sha256"] != original_sha256(inputs[2])
