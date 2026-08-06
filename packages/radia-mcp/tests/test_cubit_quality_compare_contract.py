"""Fast contracts for cubit_netgen_quality_compare (no Cubit, no
netgen needed -- error paths and input validation only; the real
three-way comparison lives in validation_test/radia_mcp)."""

import json

import pytest

from radia_mcp.cubit import server


cubit_netgen_quality_compare = server.cubit_netgen_quality_compare


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


@pytest.mark.parametrize("kwargs", [
    {"order": "bad"},
    {"netgen_maxh": -1.0, "cubit_size": 1.0},
    {"netgen_maxh": 1.0, "cubit_size": float("nan")},
    {"threshold": float("inf")},
    {"schemes": []},
    {"schemes": ["netgen", "netgen"]},
])
def test_invalid_options_fail_before_meshing(tmp_path, kwargs):
    step = tmp_path / "x.step"
    step.write_text("dummy", encoding="ascii")
    out = json.loads(cubit_netgen_quality_compare(str(step), **kwargs))
    assert out["status"] == "error"
    assert out["kind"] == "input"


def test_cubit_referee_exception_stays_in_route_row(tmp_path, monkeypatch):
    step = tmp_path / "x.step"
    step.write_text("dummy", encoding="ascii")
    monkeypatch.setattr(
        server, "_run_batch",
        lambda *_args, **_kwargs: {"status": "ok", "summary": {}})

    from radia_mcp.gmsh import msh_inspect

    def fail_referee(*_args, **_kwargs):
        raise RuntimeError("broken quality file")

    monkeypatch.setattr(msh_inspect, "mesh_quality", fail_referee)
    out = json.loads(cubit_netgen_quality_compare(
        str(step), netgen_maxh=1.0, cubit_size=1.0,
        schemes=["cubit_tet"]))
    assert out["status"] == "ok"
    assert out["all_routes_completed"] is False
    assert out["rows"][0]["status"] == "error"
    assert out["rows"][0]["kind"] == "referee"
    assert "broken quality file" in out["rows"][0]["error"]
