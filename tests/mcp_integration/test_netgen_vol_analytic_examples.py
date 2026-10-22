"""Lock calculation outputs to frozen prior records; analytic acceptance is separate."""
import hashlib
import json
from pathlib import Path
import pytest
from validation_test.radia_mcp.netgen_vol_examples import CASES, example_module

FIXTURES = Path(__file__).with_name("fixtures") / "netgen_vol_analytic"

@pytest.mark.parametrize("module_name", CASES, ids=CASES)
def test_frozen_summary(module_name, tmp_path, monkeypatch):
    manifest = json.loads((FIXTURES / "manifest.json").read_text(encoding="utf-8"))
    reference = manifest["cases"][module_name]
    frozen = (FIXTURES / reference["file"]).read_bytes()
    assert hashlib.sha256(frozen).hexdigest() == reference["sha256"]
    with example_module(module_name) as module:
        def forbidden_acceptance(*args, **kwargs):
            raise AssertionError("Analytic acceptance belongs to validation_test")
        monkeypatch.setattr(module, "validate_summary", forbidden_acceptance)
        summary = module.build_summary()
    output = tmp_path / reference["file"]
    output.write_text(json.dumps(summary, sort_keys=True), encoding="utf-8")
    assert json.loads(output.read_text(encoding="utf-8")) == json.loads(frozen)
