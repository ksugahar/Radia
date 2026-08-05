"""Fast error-contract test for the canonical cubit_check_vol gate."""

import json

from radia_mcp.cubit.server import cubit_check_vol


def test_missing_vol_returns_structured_error(tmp_path):
    out = json.loads(cubit_check_vol(str(tmp_path / "does_not_exist.vol")))
    assert out["status"] == "error"
    assert out["stage"] in {"import", "input", "check"}
    assert "error" in out
