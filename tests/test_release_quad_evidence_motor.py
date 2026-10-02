"""release_quad evidence-motor rebinds the proof, then runs the acceptance tests unchanged."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

_TOOL = Path(__file__).resolve().parents[1] / "tools" / "release_quad.py"
_SPEC = importlib.util.spec_from_file_location("radia_release_quad_motor", _TOOL)
assert _SPEC is not None and _SPEC.loader is not None
release_quad = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(release_quad)


def _snapshot():
    return {
        rel: (release_quad.REPO / rel).read_bytes() for rel in release_quad._MOTOR_ACCEPTANCE_TESTS
    }


def _recording_run(returncodes):
    calls = []

    def fake_run(cmd, **kw):
        calls.append(list(cmd))
        return SimpleNamespace(returncode=returncodes[len(calls) - 1])

    return calls, fake_run


def test_motor_paths_exist_in_the_repository():
    paths = (
        release_quad._MOTOR_ARTIFACT,
        release_quad._MOTOR_PROOF,
        release_quad._MOTOR_PROOF_GENERATOR,
        *release_quad._MOTOR_ACCEPTANCE_TESTS,
    )
    assert [rel for rel in paths if not (release_quad.REPO / rel).is_file()] == []
    assert "tests/mcp_integration/test_matlab_source_contract.py" in (
        release_quad._MOTOR_ACCEPTANCE_TESTS
    )


def test_proof_is_regenerated_before_acceptance():
    calls, fake_run = _recording_run([0, 0])

    assert release_quad._motor_post_regeneration(fake_run) == 0
    assert calls[0] == [sys.executable, release_quad._MOTOR_PROOF_GENERATOR]
    assert calls[1][1:3] == ["-m", "pytest"]
    assert calls[1][-len(release_quad._MOTOR_ACCEPTANCE_TESTS) :] == list(
        release_quad._MOTOR_ACCEPTANCE_TESTS
    )
    assert len(calls) == 2


def test_proof_failure_stops_before_acceptance():
    calls, fake_run = _recording_run([1])

    assert release_quad._motor_post_regeneration(fake_run) == 4
    assert len(calls) == 1


def test_acceptance_failure_propagates_without_editing_tests():
    before = _snapshot()
    calls, fake_run = _recording_run([0, 1])

    assert release_quad._motor_post_regeneration(fake_run) == 4
    assert len(calls) == 2
    assert _snapshot() == before


def test_acceptance_fails_on_a_missing_test_without_running(monkeypatch):
    monkeypatch.setattr(release_quad, "_MOTOR_ACCEPTANCE_TESTS", ("tests/does_not_exist_motor.py",))

    def fake_run(cmd, **kw):
        raise AssertionError("pytest must not run when a test file is missing")

    assert release_quad._motor_acceptance(fake_run) == 4
