"""release_quad evidence-motor runs the existing acceptance tests and never edits them."""
from __future__ import annotations

import importlib.util
from pathlib import Path
from types import SimpleNamespace

_TOOL = Path(__file__).resolve().parents[1] / "tools" / "release_quad.py"
_SPEC = importlib.util.spec_from_file_location("radia_release_quad_motor", _TOOL)
assert _SPEC is not None and _SPEC.loader is not None
release_quad = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(release_quad)


def _snapshot():
    return {
        rel: (release_quad.REPO / rel).read_bytes()
        for rel in release_quad._MOTOR_ACCEPTANCE_TESTS
    }


def test_acceptance_paths_exist_in_the_repository():
    assert release_quad._MOTOR_ACCEPTANCE_TESTS
    missing = [
        rel for rel in release_quad._MOTOR_ACCEPTANCE_TESTS
        if not (release_quad.REPO / rel).is_file()
    ]
    assert missing == []


def test_acceptance_runs_the_tests_unchanged_and_propagates_failure():
    before = _snapshot()
    calls = []

    def fake_run(cmd, **kw):
        calls.append(cmd)
        return SimpleNamespace(returncode=1)

    assert release_quad._motor_acceptance(fake_run) == 4
    assert len(calls) == 1
    assert calls[0][1:3] == ["-m", "pytest"]
    assert list(calls[0][-len(release_quad._MOTOR_ACCEPTANCE_TESTS):]) == list(
        release_quad._MOTOR_ACCEPTANCE_TESTS
    )
    assert _snapshot() == before


def test_acceptance_passes_only_when_pytest_passes():
    assert release_quad._motor_acceptance(lambda cmd, **kw: SimpleNamespace(returncode=0)) == 0


def test_acceptance_fails_on_a_missing_test_without_running(monkeypatch):
    monkeypatch.setattr(
        release_quad, "_MOTOR_ACCEPTANCE_TESTS", ("tests/does_not_exist_motor.py",)
    )

    def fake_run(cmd, **kw):
        raise AssertionError("pytest must not run when a test file is missing")

    assert release_quad._motor_acceptance(fake_run) == 4
