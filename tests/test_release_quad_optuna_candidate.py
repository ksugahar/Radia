from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import tomllib
from pathlib import Path
from types import SimpleNamespace


TOOL = Path(__file__).resolve().parents[1] / "tools" / "release_quad.py"
SPEC = importlib.util.spec_from_file_location("radia_release_quad_optuna", TOOL)
assert SPEC is not None and SPEC.loader is not None
release_quad = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(release_quad)


def test_release_quad_tracks_the_independent_optuna_version():
    versions = release_quad._read_repo_versions()
    package = TOOL.parents[1] / "packages/radia-optuna/pyproject.toml"
    expected = tomllib.loads(package.read_text(encoding="utf-8"))["project"]["version"]
    assert versions["radia-optuna"] == expected
    assert versions["optuna.__version__"] == expected


def test_optuna_candidate_requires_the_independent_distribution_workflow():
    source = TOOL.read_text(encoding="utf-8")
    assert '"workflowName": "radia-optuna"' in source
    assert '"workflowName": "CI"' not in source


def test_optuna_workflow_covers_its_release_quad_acceptance_contract():
    """The independent distribution CI runs when its acceptance code changes."""
    import yaml

    root = Path(__file__).resolve().parents[1]
    workflow = yaml.safe_load(
        (root / ".github/workflows/radia-optuna.yml").read_text(encoding="utf-8")
    )
    triggers = workflow[True] if True in workflow else workflow["on"]
    owned = {"tools/release_quad.py", "tests/test_release_quad_optuna_candidate.py"}
    for event in ("push", "pull_request"):
        assert owned <= set(triggers[event]["paths"]), event
    build = next(
        step["run"] for step in workflow["jobs"]["build-test"]["steps"]
        if step.get("name") == "Test, build, and verify radia-optuna wheel"
    )
    assert "-m pytest tests\\test_release_quad_optuna_candidate.py" in build
    # The exact successful push-event binding stays in force.
    source = TOOL.read_text(encoding="utf-8")
    assert '"event": "push"' in source and '"headSha": head' in source


def test_optuna_workflow_installs_pinned_pybind_before_native_build():
    root = Path(__file__).resolve().parents[1]
    workflow = (root / ".github/workflows/radia-optuna.yml").read_text(
        encoding="utf-8"
    )
    dependency = workflow.index(
        "& $env:RADIA_OPTUNA_CI_PYTHON -m pip install pybind11==3.0.2 ninja"
    )
    build = workflow.index("& .\\Build.ps1 -OptunaMexOnly")
    assert dependency < build
    assert (
        "& $env:RADIA_OPTUNA_CI_PYTHON -m pip install "
        "build wheel pytest optuna==5.0.0"
    ) in workflow


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

    state["targets"]["mdx2"]["status"] = "failed"
    path.write_text(json.dumps(state), encoding="utf-8")
    rc, checked = release_quad._verify_optuna_candidate_state(str(wheel))
    assert rc == 4
    assert checked is None


def _fake_git(head, origin_main, dirty=""):
    def run_git(*argv, **_kwargs):
        if argv[:2] == ("rev-parse", "HEAD"):
            return SimpleNamespace(returncode=0, stdout=head + "\n", stderr="")
        if argv[:2] == ("rev-parse", "origin/main"):
            return SimpleNamespace(returncode=0, stdout=origin_main + "\n", stderr="")
        if argv[0] == "status":
            return SimpleNamespace(returncode=0, stdout=dirty, stderr="")
        if argv[0] == "fetch":
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        raise AssertionError(f"unexpected git call: {argv}")

    return run_git


def test_optuna_source_identity_is_the_checkout_not_the_solver_tag(monkeypatch):
    main = "d" * 40

    def solver_tag():
        raise AssertionError("radia-optuna must not resolve the solver release tag")

    monkeypatch.setattr(release_quad, "_release_commit", solver_tag)
    monkeypatch.setattr(release_quad, "_release_head", solver_tag)
    monkeypatch.setattr(release_quad, "_git", _fake_git(main, main))
    assert release_quad._optuna_release_source_ready() == (True, main)


def test_optuna_source_rejects_a_checkout_behind_main_or_dirty(monkeypatch):
    main = "d" * 40
    behind = "1" * 40
    monkeypatch.setattr(release_quad, "_git", _fake_git(behind, main))
    ready, message = release_quad._optuna_release_source_ready()
    assert ready is False
    assert behind in message and main in message

    monkeypatch.setattr(release_quad, "_git", _fake_git(main, main, dirty=" M x.py\n"))
    assert release_quad._optuna_release_source_ready() == (
        False,
        "tracked release source is dirty",
    )


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
