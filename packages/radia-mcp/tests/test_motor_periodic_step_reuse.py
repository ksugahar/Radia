"""Surrogate tests for the motor_periodic_step_reuse_run tool and its owned worker.

Guards: request bounds (only the qualified synthetic ring, amg/sparsecholesky,
no unknown keys); missing-runtime reporting; on a runtime that provides
radia.periodic_step_reuse, signed and unsigned pitch-mapped states against an
independent dense per-step oracle, the per-step true residual against b_eff,
explicit zero right-hand side, one observed setup reused by every step, the
artifact hash over its final bytes, and a worker gate that rejects a perturbed
solution; the MCP tool's registration and invalid-input path.
Assumes: the ring is the SPD scalar system of Radia's own surrogate tests.
"""

from __future__ import annotations

import asyncio
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
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


def _solver_module():
    pytest.importorskip("ngsolve")
    return pytest.importorskip("radia.periodic_step_reuse")


def _dense_oracle(request):
    # Independent reference: dense A, explicit signed permutation matrix with
    # P[i, perm[i]] = sign[i], numpy.linalg.solve per step, same seeded loads.
    n, dt = request["n"], request["dt"]
    K, M = worker.ring_system(n)
    A = K.toarray() + M.toarray() / dt
    P = np.zeros((n, n))
    rows = np.arange(n)
    signs = np.where(rows % 2 == 0, 1.0, -1.0) if request["signed"] else np.ones(n)
    P[rows, (rows + request["pitch"]) % n] = signs
    rng = np.random.default_rng(request["seed"])
    x = np.zeros(n)
    for _ in range(request["steps"]):
        rhs = rng.standard_normal(n) if request["rhs"] == "random" else np.zeros(n)
        x = np.linalg.solve(A, rhs + M.toarray() @ (P @ x) / dt)
    return x


@pytest.mark.parametrize("method", ["sparsecholesky", "amg"])
@pytest.mark.parametrize("signed", [False, True])
def test_matches_dense_oracle_with_one_observed_setup(tmp_path, method, signed):
    _solver_module()
    from radia._residual_gate import RELATIVE_LIMIT

    request = _request(tmp_path, method=method, signed=signed, pitch=3)
    result = worker.run(request)
    expected = _dense_oracle(request)

    final = np.array(json.loads(Path(result["artifact"]["path"]).read_bytes())["final_state"])
    error = np.linalg.norm(final - expected) / np.linalg.norm(expected)
    assert error <= 1e-8, f"actual relative error={error}, expected<=1e-8"
    assert result["status"] == "ok"
    assert len(result["steps"]) == request["steps"]
    for record in result["steps"]:
        assert record["relative_residual"] <= RELATIVE_LIMIT, f"actual={record}"
        assert record["solver_relative_residual"] <= RELATIVE_LIMIT, f"actual={record}"
        assert record["b_eff_norm"] > 0.0 and record["zero_rhs"] is False
    assert result["residual_max"] == max(r["relative_residual"] for r in result["steps"])
    observed = result["observed"]
    assert observed["setup_builds_observed"] == 1, f"actual={observed}"
    assert observed["setup_builds_during_steps"] == 0, f"actual={observed}"
    assert observed["reuse_count"] == request["steps"], f"actual={observed}"
    assert observed["solver_steps_solved"] == request["steps"], f"actual={observed}"
    assert result["scope"]["motor_validation"] is False and result["scope"]["ams"] is False


def test_artifact_hash_covers_final_bytes(tmp_path):
    _solver_module()
    result = worker.run(_request(tmp_path))
    artifact = result["artifact"]
    data = Path(artifact["path"]).read_bytes()
    assert hashlib.sha256(data).hexdigest() == artifact["sha256"]
    assert len(data) == artifact["bytes"]
    body = json.loads(data)
    assert body["schema"] == worker.SCHEMA
    assert body["provenance"]["entry"] == "radia.periodic_step_reuse.PeriodicStepReuseSolver"
    assert body["steps"] == result["steps"]


def test_zero_rhs_is_explicit(tmp_path):
    _solver_module()
    result = worker.run(_request(tmp_path, rhs="zero", steps=3))
    for record in result["steps"]:
        assert record["zero_rhs"] is True and record["b_eff_norm"] == 0.0, f"actual={record}"
        assert record["relative_residual"] == 0.0, f"actual={record}"
    final = json.loads(Path(result["artifact"]["path"]).read_bytes())["final_state"]
    assert final == [0.0] * 12


def test_worker_gate_rejects_a_perturbed_solution(tmp_path, monkeypatch):
    module = _solver_module()
    original = module.PeriodicStepReuseSolver.step

    def perturbed(self, rhs, x_prev=None):
        return original(self, rhs, x_prev=x_prev) + 1e-3

    monkeypatch.setattr(module.PeriodicStepReuseSolver, "step", perturbed)
    with pytest.raises(RuntimeError, match="independent true relative residual"):
        worker.run(_request(tmp_path))
    assert list(tmp_path.iterdir()) == [], "no artifact may be written"


def test_tool_is_registered_and_rejects_invalid_json():
    from radia_mcp.motor import server

    tools = asyncio.run(server.mcp.list_tools())
    assert "motor_periodic_step_reuse_run" in {tool.name for tool in tools}
    result = json.loads(asyncio.run(server.motor_periodic_step_reuse_run("{not json")))
    assert result["status"] == "invalid_input", f"actual={result}"
    assert result["schema"] == worker.SCHEMA


def test_identical_concurrent_calls_own_separate_artifacts(tmp_path, monkeypatch):
    # Two identical requests through the real tool run in two worker processes
    # at once; each must report its own finished file whose final bytes match
    # the reported SHA-256, and no partial file may remain.
    _solver_module()
    from radia_mcp.motor import server

    monkeypatch.setenv("PYTHONPATH", os.pathsep.join(p for p in sys.path if p))
    request_json = json.dumps(_request(tmp_path))

    async def both():
        return await asyncio.gather(server.motor_periodic_step_reuse_run(request_json),
                                    server.motor_periodic_step_reuse_run(request_json))

    results = [json.loads(text) for text in asyncio.run(both())]
    assert [r["status"] for r in results] == ["ok", "ok"], f"actual={results}"
    paths = [Path(r["artifact"]["path"]) for r in results]
    assert paths[0] != paths[1], f"actual paths={paths}"
    for result, path in zip(results, paths):
        data = path.read_bytes()
        assert hashlib.sha256(data).hexdigest() == result["artifact"]["sha256"]
        assert len(data) == result["artifact"]["bytes"]
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted(p.name for p in paths)


def test_artifact_writer_never_reuses_a_file(tmp_path):
    params = worker.parse_request(_request(tmp_path))
    first = worker.write_artifact(tmp_path, params, b"one\n")
    second = worker.write_artifact(tmp_path, params, b"two\n")
    assert first != second
    assert first.read_bytes() == b"one\n" and second.read_bytes() == b"two\n"
    assert not list(tmp_path.glob("*.partial"))
