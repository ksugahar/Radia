"""Solver-free contracts of the motor_periodic_step_reuse_run tool and its worker.

Guards: request bounds (only the qualified synthetic ring, amg/sparsecholesky,
no unknown keys); missing-runtime reporting without an artifact; one owned
artifact file per write; the MCP tool's registration and invalid-input path.
The solver-backed checks (dense oracle, true residual, setup reuse, artifact
hash, concurrent workers) run in the solver-numerics lane:
validation_test/radia_mcp/test_motor_periodic_step_reuse.py.
"""

from __future__ import annotations

import asyncio
import json
import sys

import pytest

from radia_mcp.motor import periodic_step_reuse_worker as worker


def _request(tmp_path, **overrides):
    request = {
        "case": "synthetic_periodic_ring", "n": 12, "steps": 8, "dt": 0.01, "pitch": 1,
        "signed": True, "method": "sparsecholesky", "rhs": "random", "seed": 42,
        "output_dir": str(tmp_path),
    }
    request.update(overrides)
    return request


@pytest.mark.parametrize("overrides", [
    {"case": "motor"},
    {"n": 3},
    {"n": worker.MAX_DOFS + 1},
    {"n": True},
    {"steps": 0},
    {"steps": worker.MAX_STEPS + 1},
    {"pitch": 0},
    {"pitch": 12},
    {"dt": 0.0},
    {"dt": float("nan")},
    {"signed": 1},
    {"method": "ams"},
    {"method": "lu"},
    {"rhs": "ones"},
    {"seed": -1},
    {"output_dir": "relative/dir"},
    {"extra": 1},
])
def test_rejects_requests_outside_the_qualified_case(tmp_path, overrides):
    with pytest.raises(ValueError):
        worker.parse_request(_request(tmp_path, **overrides))


def test_missing_output_dir_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="does not exist"):
        worker.parse_request(_request(tmp_path, output_dir=str(tmp_path / "absent")))


def test_missing_runtime_is_reported_not_crashed(tmp_path, monkeypatch, capsys):
    monkeypatch.setitem(sys.modules, "radia.periodic_step_reuse", None)
    with pytest.raises(worker.UnsupportedRuntime):
        worker.run(_request(tmp_path))
    monkeypatch.setattr(sys, "stdin", _Stdin(json.dumps(_request(tmp_path))))
    assert worker.main() == 0
    result = json.loads(capsys.readouterr().out)
    assert result["status"] == "unsupported_runtime", f"actual={result}"
    assert list(tmp_path.iterdir()) == [], "no artifact may be written"


class _Stdin:
    def __init__(self, text):
        self._text = text

    def read(self):
        return self._text


def test_tool_is_registered_and_rejects_invalid_json():
    from radia_mcp.motor import server

    tools = asyncio.run(server.mcp.list_tools())
    assert "motor_periodic_step_reuse_run" in {tool.name for tool in tools}
    result = json.loads(asyncio.run(server.motor_periodic_step_reuse_run("{not json")))
    assert result["status"] == "invalid_input", f"actual={result}"
    assert result["schema"] == worker.SCHEMA


def test_artifact_writer_never_reuses_a_file(tmp_path):
    params = worker.parse_request(_request(tmp_path))
    first = worker.write_artifact(tmp_path, params, b"one\n")
    second = worker.write_artifact(tmp_path, params, b"two\n")
    assert first != second
    assert first.read_bytes() == b"one\n" and second.read_bytes() == b"two\n"
    assert not list(tmp_path.glob("*.partial"))
