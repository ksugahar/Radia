from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace


TOOL = Path(__file__).resolve().parents[1] / "tools" / "release_quad.py"
SPEC = importlib.util.spec_from_file_location("radia_release_quad_optuna", TOOL)
assert SPEC is not None and SPEC.loader is not None
release_quad = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(release_quad)


def test_release_quad_tracks_the_independent_optuna_version():
    versions = release_quad._read_repo_versions()
    assert versions["radia-optuna"] == "0.1.4"
    assert versions["optuna.__version__"] == "0.1.4"


def test_optuna_candidate_records_every_machine_for_one_exact_wheel(
    monkeypatch, tmp_path
):
    wheel = tmp_path / "radia_optuna-0.1.1-py3-none-win_amd64.whl"
    wheel.write_bytes(b"exact-ci-wheel")
    digest = hashlib.sha256(wheel.read_bytes()).hexdigest()
    monkeypatch.setattr(release_quad, "OPTUNA_GATE_ROOT", tmp_path / "state")
    monkeypatch.setattr(
        release_quad,
        "_download_verified_optuna_ci_wheel",
        lambda _run_id: (
            {
                "ci_run_id": "12345",
                "ci_url": "https://example.invalid/run/12345",
                "commit": "a" * 40,
                "wheel": str(wheel),
                "wheel_sha256": digest,
                "version": "0.1.1",
            },
            "verified",
        ),
    )
    calls = []

    def pass_target(key, candidate_wheel, candidate_hash):
        calls.append((key, candidate_wheel, candidate_hash))
        return True, release_quad.OPTUNA_SUCCESS_MARKER

    monkeypatch.setattr(release_quad, "_run_optuna_candidate_target", pass_target)
    args = argparse.Namespace(ci_run_id="12345", target="all")
    assert release_quad.cmd_optuna_candidate(args) == 0
    assert [key for key, _, _ in calls] == list(release_quad.SIMULINK_TARGETS)
    state = json.loads(
        release_quad._optuna_state_path(digest).read_text(encoding="utf-8")
    )
    assert state["wheel_sha256"] == digest
    assert set(state["targets"]) == set(release_quad.SIMULINK_TARGETS)
    assert {row["status"] for row in state["targets"].values()} == {"passed"}


def test_optuna_done_requires_exact_head_hash_version_and_four_targets(
    monkeypatch, tmp_path
):
    wheel = tmp_path / "radia_optuna-0.1.1-py3-none-win_amd64.whl"
    wheel.write_bytes(b"four-machine-candidate")
    digest = hashlib.sha256(wheel.read_bytes()).hexdigest()
    head = "b" * 40
    monkeypatch.setattr(release_quad, "OPTUNA_GATE_ROOT", tmp_path / "state")
    monkeypatch.setattr(
        release_quad,
        "_verify_optuna_wheel",
        lambda _wheel: ({"ok": True, "version": "0.1.1"}, "verified"),
    )
    monkeypatch.setattr(
        release_quad, "_optuna_release_source_ready", lambda: (True, head)
    )
    state = {
        "schema": "radia.release-quad.optuna-candidate.v1",
        "ci_run_id": "12345",
        "commit": head,
        "wheel": str(wheel),
        "wheel_sha256": digest,
        "version": "0.1.1",
        "targets": {
            key: {"status": "passed"} for key in release_quad.SIMULINK_TARGETS
        },
    }
    path = release_quad._optuna_state_path(digest)
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(state), encoding="utf-8")

    rc, checked = release_quad._verify_optuna_candidate_state(str(wheel))
    assert rc == 0
    assert checked == state

    state["targets"]["mdx"]["status"] = "failed"
    path.write_text(json.dumps(state), encoding="utf-8")
    rc, checked = release_quad._verify_optuna_candidate_state(str(wheel))
    assert rc == 4
    assert checked is None


def test_installed_wheel_runner_emits_quad_success_marker_and_checks_notices():
    root = Path(__file__).resolve().parents[1]
    runner = (
        root / "packages/radia-optuna/tests/run_installed_wheel_simulink.ps1"
    ).read_text(encoding="utf-8")
    doctor = (
        root / "packages/radia-optuna/src/radia_optuna/cli.py"
    ).read_text(encoding="utf-8")
    assert "RADIA_OPTUNA_WHEEL_SIMULINK_OK" in runner
    assert "PythonExecutable" in runner
    assert "PreverifiedWheelSha256" in runner
    assert "Get-FileHash" in runner
    assert "upstream_notices_complete" in doctor
    assert "notices_complete" in doctor


def test_local_candidate_decodes_matlab_output_as_utf8(monkeypatch, tmp_path):
    wheel = tmp_path / "radia_optuna-0.1.1-py3-none-win_amd64.whl"
    wheel.write_bytes(b"platform-wheel")
    observed = {}

    def completed(_command, **kwargs):
        observed["command"] = _command
        observed.update(kwargs)
        return SimpleNamespace(
            returncode=0,
            stdout=f"・ MATLAB ready\n{release_quad.OPTUNA_SUCCESS_MARKER}\n",
            stderr="",
        )

    monkeypatch.setattr(release_quad.subprocess, "run", completed)
    passed, output = release_quad._run_optuna_candidate_target(
        "lab", wheel, hashlib.sha256(wheel.read_bytes()).hexdigest()
    )

    assert passed is True
    assert "・ MATLAB ready" in output
    assert observed["encoding"] == "utf-8"
    assert observed["errors"] == "replace"
    assert "-PreverifiedWheelSha256" in observed["command"]
    assert hashlib.sha256(wheel.read_bytes()).hexdigest() in observed["command"]
