from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
import tomllib
import types
from pathlib import Path
from types import SimpleNamespace

import pytest


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


def test_genuine_optuna_build_runs_its_release_quad_candidate_contract():
    """A triggered radia-optuna build checks the candidate contract at its own SHA.

    Whether release tooling triggers that workflow is owned by
    test_test_tier_policy (it must not; fast-contracts covers it).
    """
    import yaml

    root = Path(__file__).resolve().parents[1]
    workflow = yaml.safe_load(
        (root / ".github/workflows/radia-optuna.yml").read_text(encoding="utf-8")
    )
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

    def pass_target(key, candidate_wheel, candidate_hash, engine_session=None):
        calls.append((key, candidate_wheel, candidate_hash, engine_session))
        return True, release_quad.OPTUNA_SUCCESS_MARKER

    monkeypatch.setattr(release_quad, "_run_optuna_candidate_target", pass_target)
    args = argparse.Namespace(
        ci_run_id="12345", target="all", engine_session=["mdx2=radia_shared"]
    )
    assert release_quad.cmd_optuna_candidate(args) == 0
    assert [call[0] for call in calls] == list(release_quad.SIMULINK_TARGETS)
    assert {call[0]: call[3] for call in calls if call[3]} == {"mdx2": "radia_shared"}
    state = json.loads(
        release_quad._optuna_state_path(digest).read_text(encoding="utf-8")
    )
    assert state["wheel_sha256"] == digest
    assert set(state["targets"]) == set(release_quad.SIMULINK_TARGETS)
    assert {row["status"] for row in state["targets"].values()} == {"passed"}
    assert state["targets"]["mdx2"]["engine_session"] == "radia_shared"
    assert state["targets"]["lab"]["engine_session"] is None


def test_optuna_candidate_rejects_ambiguous_engine_sessions_before_any_target(
    monkeypatch, tmp_path
):
    wheel = tmp_path / "radia_optuna-0.1.1-py3-none-win_amd64.whl"
    wheel.write_bytes(b"exact-ci-wheel")
    monkeypatch.setattr(release_quad, "OPTUNA_GATE_ROOT", tmp_path / "state")
    monkeypatch.setattr(
        release_quad,
        "_download_verified_optuna_ci_wheel",
        lambda _run_id: (
            {
                "ci_run_id": "12345",
                "commit": "a" * 40,
                "wheel": str(wheel),
                "wheel_sha256": hashlib.sha256(wheel.read_bytes()).hexdigest(),
                "version": "0.1.1",
            },
            "verified",
        ),
    )

    def no_target(*_args, **_kwargs):
        raise AssertionError("no MATLAB target may run after a bad session map")

    monkeypatch.setattr(release_quad, "_run_optuna_candidate_target", no_target)
    for target, sessions in (
        ("mdx1", ["mdx2=radia_shared"]),            # unselected target
        ("mdx2", ["mdx2=a", "mdx2=b"]),             # duplicate target
        ("mdx2", ["mdx2="]),                        # empty name
        ("mdx2", ["radia_shared"]),                 # missing HOST=
    ):
        args = argparse.Namespace(
            ci_run_id="12345", target=target, engine_session=sessions
        )
        assert release_quad.cmd_optuna_candidate(args) == 2


def test_optuna_candidate_parser_accepts_repeatable_engine_sessions(monkeypatch):
    captured = {}

    def record(args):
        captured["args"] = args
        return 0

    monkeypatch.setattr(release_quad, "cmd_optuna_candidate", record)
    monkeypatch.setattr(
        release_quad.sys, "argv",
        ["release_quad.py", "optuna-candidate", "--ci-run-id", "1",
         "--target", "lab,mdx2", "--engine-session", "lab=a",
         "--engine-session", "mdx2=b"],
    )
    try:
        release_quad.main()
    except SystemExit as stop:
        assert stop.code in (0, None)
    assert captured["args"].engine_session == ["lab=a", "mdx2=b"]


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
    worker = observed["command"][observed["command"].index("-EngineWorker") + 1]
    assert Path(worker) == release_quad.REPO / "tools/verify_simulink_release.py"
    # No selected session: the runner may only own a MATLAB it starts itself.
    assert "-EngineSession" not in observed["command"]

    release_quad._run_optuna_candidate_target(
        "lab", wheel, hashlib.sha256(wheel.read_bytes()).hexdigest(), "radia_lab"
    )
    command = observed["command"]
    assert command[command.index("-EngineSession") + 1] == "radia_lab"


def test_remote_candidate_copies_the_engine_worker_and_selects_the_session(
    monkeypatch, tmp_path
):
    wheel = tmp_path / "radia_optuna-0.1.1-py3-none-win_amd64.whl"
    wheel.write_bytes(b"platform-wheel")
    digest = hashlib.sha256(wheel.read_bytes()).hexdigest()
    copies, scripts = [], []

    def completed(command, **kwargs):
        if command[0] == "scp":
            copies.append((Path(command[1]).name, command[2]))
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        assert command[0] == "ssh" and command[1] == "mdx2"
        scripts.append(kwargs["input"])
        return SimpleNamespace(
            returncode=0, stdout=release_quad.OPTUNA_SUCCESS_MARKER, stderr=""
        )

    monkeypatch.setattr(release_quad.subprocess, "run", completed)
    passed, _ = release_quad._run_optuna_candidate_target(
        "mdx2", wheel, digest, "o'shared"
    )
    assert passed is True
    remote_root = f"mdx2:C:/temp/radia-release-quad/optuna-{digest[:16]}"
    assert ("verify_simulink_release.py",
            f"{remote_root}/verify_simulink_release.py") in copies
    invocation = scripts[-1]
    worker = "C:\\temp\\radia-release-quad\\optuna-" + digest[:16] + \
        "\\verify_simulink_release.py"
    assert f"-EngineWorker '{worker}'" in invocation
    # The PowerShell literal doubles an embedded quote.
    assert "-EngineSession 'o''shared'" in invocation

    scripts.clear()
    release_quad._run_optuna_candidate_target("mdx2", wheel, digest)
    assert "-EngineWorker" in scripts[-1]
    assert "-EngineSession" not in scripts[-1]


RUNNER = (
    Path(__file__).resolve().parents[1]
    / "packages/radia-optuna/tests/run_installed_wheel_simulink.ps1"
)


def _runner_expression(session: str) -> str:
    """Evaluate the runner's own expression lines with fixed inputs."""
    lines = RUNNER.read_text(encoding="utf-8").splitlines()
    source = [line.strip() for line in lines
              if line.strip().startswith(("$pathReset = ", "$batch = "))]
    assert len(source) == 2
    script = "\n".join([
        f"$EngineSession = '{session}'",
        "$matlabPathLiteral = 'C:\\venv\\matlab'",
        "$testDirectoryLiteral = 'C:\\tests\\matlab'",
        "$simulinkEvidenceLiteral = 'C:\\temp\\run\\simulink-evidence.json'",
        *source,
        "[Console]::Out.Write($batch)",
    ])
    result = subprocess.run(
        ["pwsh", "-NoProfile", "-NonInteractive", "-Command", "-"],
        input=script, capture_output=True, text=True, encoding="utf-8",
        check=True, timeout=60,
    )
    return result.stdout


def _top_level_statements(expression: str) -> list[str]:
    statements, depth, quoted, current = [], 0, False, ""
    for character in expression:
        if character == "'":
            quoted = not quoted
        elif not quoted and character in "([{":
            depth += 1
        elif not quoted and character in ")]}":
            depth -= 1
        if character == ";" and depth == 0 and not quoted:
            statements.append(current.strip())
            current = ""
        else:
            current += character
    return [statement for statement in statements + [current.strip()] if statement]


@pytest.mark.parametrize("blocked", [None, "model", "mex"])
def test_runner_expression_leaves_a_borrowed_session_as_found(
    monkeypatch, tmp_path, blocked
):
    """Drive the shared worker with the runner's real expression."""
    expression = _runner_expression("owned_by_user")
    statements = _top_level_statements(expression)
    # A borrowed session keeps its path: no reset, candidates prepended.
    assert not any("restoredefaultpath" in statement for statement in statements)
    assert statements[-1].startswith("run_installed_wheel_acceptance(")

    worker_spec = importlib.util.spec_from_file_location(
        "optuna_engine_worker", TOOL.parent / "verify_simulink_release.py"
    )
    worker = importlib.util.module_from_spec(worker_spec)
    worker_spec.loader.exec_module(worker)
    base = {"result": "user result", "fileId": 7, "cleanupFile": "user"}
    state = {"path": "user path", "pwd": "user folder",
             "env": {"PATH": "user PATH"}}
    evaluated = []

    class Engine:
        def feature(self, _name): return 123
        def find_system(self, *_args): return ["user_model"] if blocked == "model" else []
        def inmem(self, **_kwargs): return [], ["optuna_mex"] if blocked == "mex" else []
        def path(self, *args, **_kwargs):
            if args:
                state["path"] = args[0]
            return state["path"]
        def pwd(self): return state["pwd"]
        def cd(self, value, **_kwargs): state["pwd"] = value
        def getenv(self, name): return state["env"].get(name, "")
        def setenv(self, name, value, **_kwargs): state["env"][name] = value
        def matlabroot(self): return str(tmp_path)
        def eval(self, text, **_kwargs):
            evaluated.append(text)
            # MATLAB assigns base variables only for top-level "name =" forms.
            for statement in _top_level_statements(text):
                name, separator, _ = statement.partition("=")
                if separator and name.strip().isidentifier():
                    base[name.strip()] = "overwritten"
                if statement.startswith("addpath("):
                    state["path"] = "candidate;" + state["path"]
            if text == expression:
                state["pwd"] = "C:\\temp\\scratch"
                raise RuntimeError("acceptance failed")
        def quit(self): pytest.fail("a borrowed MATLAB must stay alive")

    api = types.ModuleType("matlab.engine")
    api.find_matlab = lambda: ("owned_by_user",)
    api.connect_matlab = lambda name: Engine()
    api.start_matlab = lambda *_: pytest.fail("must reuse the selected Engine")
    parent = types.ModuleType("matlab")
    parent.engine = api
    monkeypatch.setitem(sys.modules, "matlab", parent)
    monkeypatch.setitem(sys.modules, "matlab.engine", api)
    monkeypatch.setattr(worker, "_matlab_process_ids", lambda: {123})

    with pytest.raises(RuntimeError):
        worker._engine_worker(str(tmp_path), expression, "owned_by_user")
    if blocked:
        # Refused before the acceptance touched the session.
        assert evaluated == []
    else:
        assert evaluated == [expression, "clear radia_mex optuna_mex"]
    assert base == {"result": "user result", "fileId": 7, "cleanupFile": "user"}
    assert state["path"] == "user path" and state["pwd"] == "user folder"
    assert state["env"]["PATH"] == "user PATH"


def test_runner_expression_resets_the_path_only_for_an_owned_session():
    statements = _top_level_statements(_runner_expression(""))
    assert statements[0] == "restoredefaultpath"
    assert statements[-1].startswith("run_installed_wheel_acceptance(")


def test_installed_wheel_runner_uses_the_isolated_engine_worker():
    runner = RUNNER.read_text(encoding="utf-8")
    # No unconditional MATLAB launch: the worker owns every start or reuse,
    # through the Engine installed into the isolated venv from this MATLAB.
    assert "-batch" not in runner
    assert "& $venvPython @workerArguments" in runner
    assert "extern\\engines\\python" in runner
    assert "matlab_engine = $engineEvidence" in runner
    # A borrowed session gets no license retry.
    assert "$maxMatlabAttempts = if ($EngineSession) { 1 } else { 3 }" in runner
