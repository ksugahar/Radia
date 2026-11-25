"""Solver-backed checks of the motor_periodic_step_reuse_run worker.

Runs Radia's PeriodicStepReuseSolver through the radia-mcp motor worker on the
qualified synthetic periodic ring and checks: signed and unsigned pitch-mapped
states against an independent dense per-step oracle, each step's true relative
residual against b_eff, explicit zero right-hand side, one observed setup reused
by every step, the artifact SHA-256 over its final bytes, the worker gate on a
perturbed solution, and two identical concurrent worker processes writing
separate artifacts.  Needs the CPU Radia build and NGSolve; the package tests
keep the solver-free contracts.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import radia.periodic_step_reuse as solver_module
from radia._residual_gate import RELATIVE_LIMIT
from radia_mcp.motor import periodic_step_reuse_worker as worker


def _request(tmp_path, **overrides):
    request = {
        "case": "synthetic_periodic_ring", "n": 12, "steps": 8, "dt": 0.01, "pitch": 1,
        "signed": True, "method": "sparsecholesky", "rhs": "random", "seed": 42,
        "output_dir": str(tmp_path),
    }
    request.update(overrides)
    return request


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


def _artifact(result):
    return json.loads(Path(result["artifact"]["path"]).read_bytes())


@pytest.mark.parametrize("method", ["sparsecholesky", "amg"])
@pytest.mark.parametrize("signed", [False, True])
def test_matches_dense_oracle_with_one_observed_setup(tmp_path, method, signed):
    request = _request(tmp_path, method=method, signed=signed, pitch=3)
    result = worker.run(request)
    expected = _dense_oracle(request)

    final = np.array(_artifact(result)["final_state"])
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
    result = worker.run(_request(tmp_path, rhs="zero", steps=3))
    for record in result["steps"]:
        assert record["zero_rhs"] is True and record["b_eff_norm"] == 0.0, f"actual={record}"
        assert record["relative_residual"] == 0.0, f"actual={record}"
    assert _artifact(result)["final_state"] == [0.0] * 12


def test_worker_gate_rejects_a_perturbed_solution(tmp_path, monkeypatch):
    original = solver_module.PeriodicStepReuseSolver.step

    def perturbed(self, rhs, x_prev=None):
        return original(self, rhs, x_prev=x_prev) + 1e-3

    monkeypatch.setattr(solver_module.PeriodicStepReuseSolver, "step", perturbed)
    with pytest.raises(RuntimeError, match="independent true relative residual"):
        worker.run(_request(tmp_path))
    assert list(tmp_path.iterdir()) == [], "no artifact may be written"


def test_identical_concurrent_workers_own_separate_artifacts(tmp_path):
    # Two identical requests in two worker processes at once; each must report
    # its own finished file whose final bytes match the reported SHA-256.
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(p for p in sys.path if p))
    payload = json.dumps(_request(tmp_path)).encode("utf-8")
    processes = [subprocess.Popen([sys.executable, "-B", worker.__file__], env=env,
                                  stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
                 for _ in range(2)]
    outputs = [process.communicate(payload, timeout=120) for process in processes]
    results = []
    for process, (stdout, stderr) in zip(processes, outputs):
        assert process.returncode == 0, stderr.decode("utf-8", errors="replace")[-1000:]
        results.append(json.loads(stdout.decode("utf-8").strip().splitlines()[-1]))
    assert [r["status"] for r in results] == ["ok", "ok"], f"actual={results}"
    paths = [Path(r["artifact"]["path"]) for r in results]
    assert paths[0] != paths[1], f"actual paths={paths}"
    for result, path in zip(results, paths):
        data = path.read_bytes()
        assert hashlib.sha256(data).hexdigest() == result["artifact"]["sha256"]
        assert len(data) == result["artifact"]["bytes"]
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted(p.name for p in paths)
