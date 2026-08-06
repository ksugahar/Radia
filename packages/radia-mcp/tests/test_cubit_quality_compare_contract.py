"""Fast contracts for cubit_netgen_quality_compare (no Cubit, no
netgen needed -- error paths and input validation only; the real
three-way comparison lives in validation_test/radia_mcp)."""

import json

from radia_mcp.cubit.server import cubit_netgen_quality_compare


def test_missing_step_is_input_error(tmp_path):
    out = json.loads(cubit_netgen_quality_compare(
        str(tmp_path / "nope.step")))
    assert out["status"] == "error"
    assert out["kind"] == "input"


def test_unknown_scheme_is_input_error(tmp_path):
    step = tmp_path / "x.step"
    step.write_text("dummy", encoding="utf-8")
    out = json.loads(cubit_netgen_quality_compare(
        str(step), schemes=["netgen", "voxels"]))
    assert out["status"] == "error"
    assert out["kind"] == "input"
    assert "voxels" in out["error"]
