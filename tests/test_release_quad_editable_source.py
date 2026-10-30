from __future__ import annotations

import base64
import hashlib
import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


_TOOL = Path(__file__).resolve().parents[1] / "tools" / "release_quad.py"
_SPEC = importlib.util.spec_from_file_location("radia_release_quad_tool", _TOOL)
assert _SPEC is not None and _SPEC.loader is not None
release_quad = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(release_quad)
_GIT = shutil.which("git")
assert _GIT is not None, "Git is required by the release-quad contract tests"




@pytest.fixture(autouse=True)
def explicit_role_bindings(monkeypatch):
    for role in ("100", "VALIDATION"):
        host="build-test" if role=="100" else "validation-test"
        for key,value in {"SSH":host,"HOSTNAME":host,"PYTHON":f"C:/Runtime/{host}/python.exe",
                          "QA_PYTHON":f"C:/QA/{host}/Scripts/python.exe",
                          "EDITABLE_REPO":f"C:/Sources/{host}","EDITABLE_COMMIT":"a"*40,
                          "EDITABLE_PATH_ALIASES":"[]"}.items():
            monkeypatch.setenv(f"RADIA_RELEASE_{key}_{role}",value)


def _git(repo, *args):
    return subprocess.run(
        [_GIT, "-C", str(repo), *args],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _fake_wheel_site(root, *, shadow=False):
    """Create a tiny Radia wheel install with a valid hash-checked RECORD."""
    site = root / "site"
    package = site / "radia"
    dist_info = site / "radia-5.1.0.dist-info"
    package.mkdir(parents=True)
    dist_info.mkdir()
    files = {
        "radia/__init__.py": b"__version__ = '5.1.0'\n",
        "radia-5.1.0.dist-info/METADATA": (
            b"Metadata-Version: 2.1\nName: radia\nVersion: 5.1.0\n"
        ),
    }
    for relative, content in files.items():
        (site / relative).write_bytes(content)
    rows = []
    for relative, content in files.items():
        digest = base64.urlsafe_b64encode(hashlib.sha256(content).digest())
        rows.append(f"{relative},sha256={digest.rstrip(b'=').decode()},{len(content)}")
    rows.append("radia-5.1.0.dist-info/RECORD,,")
    (dist_info / "RECORD").write_text("\n".join(rows) + "\n", encoding="utf-8")
    paths = [site]
    if shadow:
        shadow_root = root / "shadow"
        shadow_package = shadow_root / "radia"
        shadow_package.mkdir(parents=True)
        (shadow_package / "__init__.py").write_bytes(files["radia/__init__.py"])
        paths.insert(0, shadow_root)
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(str(path) for path in paths)
    return env


def test_release_quad_git_helper_trusts_only_its_active_worktree(monkeypatch):
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(release_quad.subprocess, "run", fake_run)
    release_quad._git("status", "--short")

    command, kwargs = calls[0]
    assert command[:3] == [
        release_quad.GIT_EXE,
        "-c",
        f"safe.directory={release_quad.REPO.resolve().as_posix()}",
    ]
    assert command[3:] == [
        "-C", str(release_quad.REPO), "status", "--short"
    ]
    assert kwargs == {"capture_output": True, "text": True, "check": True}


def test_release_head_uses_version_tag_commit(monkeypatch):
    monkeypatch.setattr(release_quad, "_radia_version", lambda: "5.1.0")
    calls = []

    def fake_git(*args, **kwargs):
        calls.append((args, kwargs))
        return subprocess.CompletedProcess(
            args, 0, stdout="a" * 40 + "\n", stderr="")

    monkeypatch.setattr(release_quad, "_git", fake_git)

    assert release_quad._release_head() == "a" * 40
    assert calls == [(('rev-list', '-n', '1', 'v5.1.0'), {'check': False})]





def test_editable_runtime_gate_requires_exact_native_manifest():
    gate = release_quad.EDITABLE_RELEASE_VERIFY
    assert "release_native_payloads.json" in gate
    assert 'radia.release-native-payloads.v1' in gate
    assert '"peec_matrices.pyd"' in gate and '"radia_motor_rom.dll"' in gate
    assert "native payload set differs" in gate
    # The retired native extension name is spelled in two pieces in the gate.
    assert 'name == "c" + "ln_core.pyd"' in gate
    assert '"locked-old" in lowered' in gate



















def test_restore_editable_is_a_tombstone_that_restores_nothing(monkeypatch, capsys):
    """The command stays only to refuse and name its replacement."""
    events = []
    monkeypatch.setattr(release_quad, "run",
                        lambda *a, **k: events.append(a) or None)
    assert release_quad.cmd_restore_editable(None) == 2
    out = capsys.readouterr().out
    assert "repoint" in out
    assert not events




def test_unc_normalization_covers_synthetic_canonical_and_release_worktrees(monkeypatch):
    monkeypatch.setattr(release_quad, "_PATH_ALIASES", (
        ("//192.0.2.10/share/repo/deep/", "x:/deep-root/"),
        ("//192.0.2.10/share/repo/", "x:/synthetic-root/"),
        ("//192.0.2.11/share/repo/", "x:/synthetic-root/"),
        ("w:/synthetic-repo/", "x:/synthetic-root/"),
    ))
    assert release_quad._norm_path(r"\\192.0.2.10\share\repo\deep\Candidate\\") == "x:/deep-root/candidate"
    assert release_quad._norm_path(r"\\192.0.2.10\share\repo\Release\src\radia") == "x:/synthetic-root/release/src/radia"
    assert release_quad._norm_path(r"\\192.0.2.11\SHARE\REPO\main\packages") == "x:/synthetic-root/main/packages"
    assert release_quad._norm_path(r"W:\synthetic-repo\main\\") == "x:/synthetic-root/main"
    assert release_quad._norm_path(r"\\unrelated-server\share\repo\main") == "//unrelated-server/share/repo/main"


def test_remote_default_aliases_match_local_normalization_contract():
    for alias, _ in release_quad._PATH_ALIASES[1:]:
        assert alias in release_quad.REMOTE_EDITABLE_VERIFY




def test_simulink_candidate_accepts_its_exact_tag_when_controller_is_newer(
        tmp_path, monkeypatch):
    repo = tmp_path / "release-controller"
    repo.mkdir()
    _git(repo, "init")
    _git(repo, "config", "user.name", "Radia Test")
    _git(repo, "config", "user.email", "radia-test@example.invalid")
    tracked = repo / "tracked.txt"
    tracked.write_text("candidate\n", encoding="ascii")
    _git(repo, "add", "tracked.txt")
    _git(repo, "commit", "-m", "candidate source")
    candidate = _git(repo, "rev-parse", "HEAD")
    _git(repo, "tag", "-a", "v4.95.75", "-m", "release", candidate)
    tracked.write_text("controller repair\n", encoding="ascii")
    _git(repo, "commit", "-am", "controller repair")
    head = _git(repo, "rev-parse", "HEAD")
    monkeypatch.setattr(release_quad, "REPO", repo)

    valid, message = release_quad._simulink_candidate_commit_is_release_anchored(
        {"commit": candidate, "radia_version": "4.95.75"}, head)

    assert valid
    assert "v4.95.75" in message

    valid, message = release_quad._simulink_candidate_commit_is_release_anchored(
        {"commit": candidate, "radia_version": "4.95.76"}, head)
    assert not valid
    assert "release tag" in message


def test_local_release_source_requires_exact_sha_and_tracked_clean(
        tmp_path, monkeypatch):
    repo = tmp_path / "release-source"
    repo.mkdir()
    _git(repo, "init")
    (repo / "tracked.txt").write_text("release\n", encoding="ascii")
    _git(repo, "add", "tracked.txt")
    _git(
        repo,
        "-c",
        "user.name=Radia Test",
        "-c",
        "user.email=radia-test@example.invalid",
        "commit",
        "-m",
        "release source",
    )
    head = _git(repo, "rev-parse", "HEAD")

    assert Path(release_quad.GIT_EXE).is_absolute()
    monkeypatch.setenv("PATH", "")
    assert release_quad._verify_local_release_source(str(repo), head) == 0
    assert release_quad._verify_local_release_source(str(repo), "0" * 40) == 4

    (repo / "tracked.txt").write_text("parallel WIP\n", encoding="ascii")
    assert release_quad._verify_local_release_source(str(repo), head) == 4


def _build_release_repo(root):
    """A tagged v5.2.1 release commit followed by a tooling repair commit."""
    repo = root / "controller"
    repo.mkdir()
    _git(repo, "init")
    _git(repo, "config", "user.name", "Radia Test")
    _git(repo, "config", "user.email", "radia-test@example.invalid")
    (repo / "pyproject.toml").write_text('# Release metadata — 日本語\n[project]\nname = "radia"\nversion = "5.2.1"\n',
                                         encoding="utf-8")
    (repo / "tool.py").write_text("release\n", encoding="ascii")
    _git(repo, "add", "pyproject.toml", "tool.py")
    _git(repo, "commit", "-m", "release")
    release = _git(repo, "rev-parse", "HEAD")
    _git(repo, "tag", "-a", "v5.2.1", "-m", "release", release)
    (repo / "tool.py").write_text("repaired release tooling\n", encoding="ascii")
    _git(repo, "commit", "-am", "tooling repair")
    return repo, release


@pytest.fixture(scope="module")
def release_template(tmp_path_factory):
    """Build the release controller and a clone of it once per module.

    Every test works on its own copy (a directory copy keeps the commits,
    tag and SHAs), so the per-test Git process cost is only what the test
    itself exercises.
    """
    root = tmp_path_factory.mktemp("release-template")
    repo, release = _build_release_repo(root)
    source = root / "release-source"
    _git(root, "clone", "-q", str(repo), str(source))
    return repo, source, release


def _release_repo(tmp_path, template):
    repo, _source, release = template
    copy = tmp_path / "controller"
    shutil.copytree(repo, copy)
    return copy, release


def _release_source(tmp_path, template):
    _repo, source, _release = template
    copy = tmp_path / "release-source"
    shutil.copytree(source, copy)
    return copy


def test_later_tooling_controller_must_descend_cleanly_from_the_release(
        tmp_path, release_template, monkeypatch):
    repo, release = _release_repo(tmp_path, release_template)
    monkeypatch.setattr(subprocess, "_text_encoding", lambda: "cp932")
    assert release_quad._verify_release_controller(repo, release, "5.2.1") == 0
    # Same version declared at the release commit is part of the identity.
    assert release_quad._verify_release_controller(repo, release, "5.2.2") == 4
    # Tracked changes in the controller are refused.
    (repo / "tool.py").write_text("unreviewed\n", encoding="ascii")
    assert release_quad._verify_release_controller(repo, release, "5.2.1") == 4
    _git(repo, "checkout", "--", "tool.py")
    # A controller that does not contain the release commit is refused.
    _git(repo, "checkout", "--orphan", "unrelated")
    _git(repo, "commit", "-m", "unrelated history")
    assert release_quad._verify_release_controller(repo, release, "5.2.1") == 4


def _done_until_source_gates(monkeypatch, controller, release):
    monkeypatch.setattr(release_quad, "REPO", controller)
    monkeypatch.setattr(release_quad, "cmd_preflight", lambda _args: 0)
    monkeypatch.setattr(release_quad, "_release_head", lambda: release)
    monkeypatch.setattr(release_quad, "_radia_version", lambda: "5.2.1")
    passed = AssertionError("source gates passed")

    def next_gate():
        raise passed

    monkeypatch.setattr(release_quad, "_verify_head_release_tag", next_gate)
    return passed


def test_done_accepts_a_later_controller_with_a_separate_exact_tag_source(
        tmp_path, monkeypatch, release_template):
    from argparse import Namespace
    controller, release = _release_repo(tmp_path, release_template)
    source = _release_source(tmp_path, release_template)
    _git(source, "checkout", "-q", "--detach", "v5.2.1")
    passed = _done_until_source_gates(monkeypatch, controller, release)
    with pytest.raises(AssertionError) as raised:
        release_quad.cmd_done(Namespace(simulink_package=None, release_source=str(source)))
    assert raised.value is passed
    # Without a separate source the later controller is itself the source: refused.
    assert release_quad.cmd_done(Namespace(simulink_package=None, release_source=None)) == 4


@pytest.mark.parametrize("case", ["stale-source", "dirty-source", "dirty-controller",
                                  "unrelated-controller"])
def test_done_rejects_a_bad_source_or_controller(tmp_path, monkeypatch, case, release_template):
    from argparse import Namespace
    controller, release = _release_repo(tmp_path, release_template)
    source = _release_source(tmp_path, release_template)
    if case != "stale-source":
        _git(source, "checkout", "-q", "--detach", "v5.2.1")
    if case == "dirty-source":
        (source / "tool.py").write_text("parallel WIP\n", encoding="ascii")
    elif case == "dirty-controller":
        (controller / "tool.py").write_text("unreviewed\n", encoding="ascii")
    elif case == "unrelated-controller":
        _git(controller, "checkout", "-q", "--orphan", "unrelated")
        _git(controller, "commit", "-q", "-m", "unrelated history")
    _done_until_source_gates(monkeypatch, controller, release)
    assert release_quad.cmd_done(
        Namespace(simulink_package=None, release_source=str(source))) == 4




def test_release_tag_gate_requires_declared_version_at_exact_head(monkeypatch):
    monkeypatch.setattr(
        release_quad, "_read_repo_versions", lambda: {"radia": "4.95.80"}
    )
    monkeypatch.setattr(release_quad, "_release_head", lambda: "a" * 40)
    monkeypatch.setattr(release_quad, "_release_tag_commit", lambda _version: "a" * 40)
    assert release_quad._verify_head_release_tag() == 0

    monkeypatch.setattr(release_quad, "_release_tag_commit", lambda _version: None)
    assert release_quad._verify_head_release_tag() == 4

    monkeypatch.setattr(release_quad, "_release_tag_commit", lambda _version: "b" * 40)
    assert release_quad._verify_head_release_tag() == 4










def test_ci_check_runs_use_latest_attempt_for_each_name(monkeypatch):
    runs = {
        "check_runs": [
            {
                "id": 1,
                "name": "fast-contracts",
                "status": "completed",
                "conclusion": "cancelled",
            },
            {
                "id": 2,
                "name": "fast-contracts",
                "status": "completed",
                "conclusion": "success",
            },
            {
                "id": 3,
                "name": "policy-checks",
                "status": "completed",
                "conclusion": "success",
            },
        ]
    }
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", lambda _path: (runs, {}))

    ok, message = release_quad._check_github_hosted_workflows(
        "a" * 40, required_names=None, timeout_sec=1, poll_sec=0
    )

    assert ok is True
    assert "2 latest SHA-bound check-runs GREEN" in message


def _green(name, run_id):
    return {"id": run_id, "name": name,
            "status": "completed", "conclusion": "success"}


def test_absent_release_check_is_not_passing_evidence(monkeypatch):
    """A SHA the release workflow never ran on must not verify.

    Radia Native Release triggers only on a v* tag or a manual dispatch, so a
    commit can be green on every check that did run while the release build
    has never seen it. Treating that silence as success is what spent tags
    4.95.86..4.95.89.
    """
    runs = {"check_runs": [_green("policy-checks", 1), _green("mcp-matrix", 2)]}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", lambda _path: (runs, {}))

    ok, message = release_quad._check_github_hosted_workflows(
        "a" * 40, required_names=None,
        require_present=[release_quad.RELEASE_CHECK_RUN],
        timeout_sec=1, poll_sec=0, registration_grace_sec=0,
    )

    assert ok is False
    assert release_quad.RELEASE_CHECK_RUN in message


def test_present_release_check_verifies_with_the_rest(monkeypatch):
    """require_present adds a requirement; it does not narrow the check set."""
    runs = {"check_runs": [
        _green("policy-checks", 1),
        _green(release_quad.RELEASE_CHECK_RUN, 2),
    ]}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", lambda _path: (runs, {}))

    ok, message = release_quad._check_github_hosted_workflows(
        "a" * 40, required_names=None,
        require_present=[release_quad.RELEASE_CHECK_RUN],
        timeout_sec=1, poll_sec=0,
    )

    assert ok is True
    assert "2 latest SHA-bound check-runs GREEN" in message


def test_require_present_still_fails_on_an_unrelated_red(monkeypatch):
    """Requiring the release build by name must not excuse any other check."""
    runs = {"check_runs": [
        _green(release_quad.RELEASE_CHECK_RUN, 1),
        {"id": 2, "name": "policy-checks",
         "status": "completed", "conclusion": "failure"},
    ]}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", lambda _path: (runs, {}))

    ok, _message = release_quad._check_github_hosted_workflows(
        "a" * 40, required_names=None,
        require_present=[release_quad.RELEASE_CHECK_RUN],
        timeout_sec=1, poll_sec=0,
    )

    assert ok is False


def test_release_check_run_name_matches_the_workflow_job():
    """RELEASE_CHECK_RUN must name a job that build-test.yml actually defines."""
    workflow = (Path(__file__).resolve().parents[1]
                / ".github" / "workflows" / "build-test.yml").read_text(encoding="utf-8")
    assert f"\n    name: {release_quad.RELEASE_CHECK_RUN}\n" in workflow
    assert release_quad.RELEASE_CHECK_RUN != "build-test"


def test_other_distribution_build_does_not_certify_native_release(monkeypatch):
    runs = {"check_runs": [_green("build-test", 1), _green("policy-checks", 2)]}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", lambda _path: (runs, {}))
    ok, message = release_quad._check_github_hosted_workflows(
        "a" * 40, require_present=[release_quad.RELEASE_CHECK_RUN],
        timeout_sec=1, poll_sec=0, registration_grace_sec=0,
    )
    assert not ok
    assert release_quad.RELEASE_CHECK_RUN in message


@pytest.mark.parametrize("conclusion", ["skipped", "neutral", "failure", "cancelled"])
@pytest.mark.parametrize("required_names", [None, ["policy-checks"]])
def test_native_evidence_must_succeed_even_when_other_checks_are_selected(
    monkeypatch, conclusion, required_names
):
    native = _green(release_quad.RELEASE_CHECK_RUN, 2)
    native["conclusion"] = conclusion
    runs = {"check_runs": [_green("policy-checks", 1), native]}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", lambda _path: (runs, {}))
    ok, message = release_quad._check_github_hosted_workflows(
        "a" * 40, required_names=required_names,
        require_present=[release_quad.RELEASE_CHECK_RUN],
        timeout_sec=1, poll_sec=0, registration_grace_sec=0,
    )
    assert not ok
    assert f"{release_quad.RELEASE_CHECK_RUN}: {conclusion}" in message


@pytest.mark.parametrize("workflow_path", sorted(release_quad.INDEPENDENT_DISTRIBUTION_WORKFLOWS))
@pytest.mark.parametrize("status,conclusion", [("queued", None), ("completed", "failure")])
def test_cpu_release_does_not_wait_for_independent_distribution_ci(
    monkeypatch, workflow_path, status, conclusion
):
    sha = "a" * 40
    native = _green(release_quad.RELEASE_CHECK_RUN, 1)
    cpu = _green("fast-contracts", 2)
    independent = dict(id=3, name="fast-contracts", status=status, conclusion=conclusion,
        details_url="https://github.com/owner/repo/actions/runs/50/job/60")
    def github(path):
        if "/check-runs?" in path:
            return {"check_runs": [native, cpu, independent]}, {}
        assert path == "repos/owner/repo/actions/runs/50"
        return dict(head_sha=sha, path=workflow_path), {}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", github)
    monkeypatch.setattr(release_quad, "_release_head", lambda: sha)
    assert release_quad.cmd_ci_verify(None) == 0


@pytest.mark.parametrize("workflow_path", [".github/workflows/radia.yml", "unknown-shared-check.yml"])
def test_cpu_release_still_rejects_cpu_and_unclassified_failures(monkeypatch, workflow_path):
    sha = "a" * 40
    failed = dict(id=2, name="shared-contracts", status="completed", conclusion="failure",
        details_url="https://github.com/owner/repo/actions/runs/50/job/60")
    def github(path):
        if "/check-runs?" in path:
            return {"check_runs": [_green(release_quad.RELEASE_CHECK_RUN, 1), failed]}, {}
        return dict(head_sha=sha, path=workflow_path), {}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", github)
    ok, message = release_quad._check_github_hosted_workflows(sha,
        require_present=[release_quad.RELEASE_CHECK_RUN],
        excluded_workflows=release_quad.INDEPENDENT_DISTRIBUTION_WORKFLOWS,
        timeout_sec=0, poll_sec=0)
    assert not ok and "shared-contracts: failure" in message


@pytest.mark.parametrize("ownership", ["lookup-error", "wrong-sha", "missing-path"])
def test_cpu_release_cannot_ignore_unverified_workflow_ownership(monkeypatch, ownership):
    sha = "a" * 40
    gpu = dict(id=2, name="cuda-mex", status="completed", conclusion="failure",
        details_url="https://github.com/owner/repo/actions/runs/50/job/60")
    def github(path):
        if "/check-runs?" in path:
            return {"check_runs": [_green(release_quad.RELEASE_CHECK_RUN, 1), gpu]}, {}
        if ownership == "lookup-error":
            raise RuntimeError("ownership API unavailable")
        return dict(head_sha="b" * 40 if ownership == "wrong-sha" else sha,
            path=None if ownership == "missing-path" else ".github/workflows/radia-simulink.yml"), {}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", github)
    ok, _ = release_quad._check_github_hosted_workflows(sha,
        require_present=[release_quad.RELEASE_CHECK_RUN],
        excluded_workflows=release_quad.INDEPENDENT_DISTRIBUTION_WORKFLOWS,
        timeout_sec=0, poll_sec=0)
    assert not ok


def test_gpu_named_native_check_cannot_replace_cpu_native_acceptance(monkeypatch):
    sha = "a" * 40
    gpu = _green(release_quad.RELEASE_CHECK_RUN, 1)
    gpu["details_url"] = "https://github.com/owner/repo/actions/runs/50/job/60"
    def github(path):
        if "/check-runs?" in path:
            return {"check_runs": [gpu]}, {}
        return dict(head_sha=sha, path=".github/workflows/radia-simulink.yml"), {}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", github)
    ok, message = release_quad._check_github_hosted_workflows(sha,
        require_present=[release_quad.RELEASE_CHECK_RUN],
        excluded_workflows=release_quad.INDEPENDENT_DISTRIBUTION_WORKFLOWS,
        timeout_sec=0, poll_sec=0, registration_grace_sec=0)
    assert not ok and release_quad.RELEASE_CHECK_RUN in message


@pytest.mark.parametrize("host", ["100"])
def test_explicit_editable_root_does_not_resolve_controller_version_tag(monkeypatch, host):
    env = getattr(release_quad, "EDITABLE_REPO_" + host.upper() + "_ENV")
    monkeypatch.setenv(env, "C:/release-quad/v5.1.0-ccccccccc/")
    monkeypatch.setattr(release_quad, "_release_commit",
                        lambda: pytest.fail("explicit root must not resolve a controller tag"))
    assert getattr(release_quad, "_editable_repo_" + host)() == "C:/release-quad/v5.1.0-ccccccccc/"




@pytest.mark.parametrize(
    "expected_version,shadow,wrong_root,expected_rc,expected_message",
    [
        ("5.1.0", False, False, 0, "validation wheel verified"),
        ("5.1.1", False, False, 1, "Radia version mismatch"),
        ("5.1.0", True, False, 1, "Shadowed Radia import"),
        ("5.1.0", False, True, 1, "Radia wheel must be under site-packages"),
    ],
)
def test_wheel_runtime_probe_checks_installed_version_and_import_origin(
    tmp_path, expected_version, shadow, wrong_root, expected_rc, expected_message
):
    env = _fake_wheel_site(tmp_path, shadow=shadow)
    code = release_quad._wheel_runtime_probe_code(expected_version, "validation")
    # Model the interpreter's install scheme without installing a test package
    # into the real environment; the production probe itself runs unchanged.
    install_root = str(tmp_path / ("other-site" if wrong_root else "site"))
    setup = (
        "import sysconfig\n"
        f"sysconfig.get_paths = lambda: {{'purelib': {install_root!r}, 'platlib': {install_root!r}}}\n"
    )

    result = subprocess.run(
        [sys.executable, "-c", setup + code],
        capture_output=True,
        text=True,
        env=env,
        timeout=30,
    )

    assert result.returncode == expected_rc
    assert expected_message in result.stdout + result.stderr


@pytest.mark.parametrize("status,conclusion", [("queued", None), ("completed", "failure")])
def test_cpu_release_excludes_a100_nightly(monkeypatch, status, conclusion):
    sha = "a" * 40
    gpu = dict(id=2, name="validate", status=status, conclusion=conclusion,
        details_url="https://github.com/owner/repo/actions/runs/50/job/60")
    def github(path):
        if "/check-runs?" in path:
            return {"check_runs": [_green(release_quad.RELEASE_CHECK_RUN, 1), gpu]}, {}
        return dict(head_sha=sha, path=".github/workflows/gpu-validation.yml"), {}
    monkeypatch.setattr(release_quad, "_git_repo_owner_name", lambda: "owner/repo")
    monkeypatch.setattr(release_quad, "gh_get", github)
    ok, message = release_quad._check_github_hosted_workflows(sha,
        require_present=[release_quad.RELEASE_CHECK_RUN],
        excluded_workflows=release_quad.INDEPENDENT_DISTRIBUTION_WORKFLOWS,
        timeout_sec=0, poll_sec=0)
    assert ok, message


@pytest.mark.parametrize("missing", ["RADIA_RELEASE_SSH_VALIDATION", "RADIA_RELEASE_PYTHON_VALIDATION"])
def test_validation_config_fails_before_any_phase8_action(monkeypatch, missing):
    """Guards: no fallback or partial deployment. Assumes: remote actions are observed mocks."""
    from argparse import Namespace
    monkeypatch.delenv(missing)
    monkeypatch.setattr(release_quad, "run", lambda *a, **kw: pytest.fail("No remote action"))
    monkeypatch.setattr(release_quad, "cmd_preflight", lambda *a: pytest.fail("No preflight before config"))
    with pytest.raises(ValueError, match=missing):
        release_quad.cmd_phase8(Namespace(target="100,validation"))


@pytest.mark.parametrize("name,value", [
    ("RADIA_RELEASE_SSH_VALIDATION", "-oProxyCommand=bad"),
    ("RADIA_RELEASE_SSH_VALIDATION", "target extra"),
    ("RADIA_RELEASE_PYTHON_VALIDATION", "python"),
    ("RADIA_RELEASE_PYTHON_VALIDATION", "C:relative/python.exe"),
    ("RADIA_RELEASE_PYTHON_VALIDATION", "C:/WindowsApps/python.exe"),
    ("RADIA_RELEASE_PYTHON_VALIDATION", 'C:/x/python.exe";bad'),
])
def test_validation_rejects_ambiguous_runtime_configuration(monkeypatch,name,value):
    """Guards: explicit executable/alias only. Assumes: Windows target semantics."""
    monkeypatch.setenv(name,value)
    with pytest.raises(ValueError,match=name):release_quad._solver_target("validation")



@pytest.mark.parametrize('role',['100','validation'])
def test_explicit_source_is_required_without_historic_defaults(monkeypatch,role):
    """Guards: no implicit source choice. Assumes: synthetic private bindings."""
    setting='RADIA_RELEASE_EDITABLE_REPO_'+role.upper()
    monkeypatch.delenv(setting)
    with pytest.raises(ValueError,match=setting):release_quad._editable_target(role)


def test_deployment_plan_names_normal_editable_and_separate_qa(monkeypatch,capsys):
    """Guards: dry planning never installs. Assumes: no actual host lookup."""
    import json
    monkeypatch.setattr(release_quad,'run',lambda *a,**k:pytest.fail('planning mutated a host'))
    assert release_quad.cmd_deployment_plan(None)==0
    plans=json.loads(capsys.readouterr().out)
    assert {p['host'] for p in plans}=={'100','validation'}
    for plan in plans:
        bound=release_quad._editable_target(plan['host'])
        assert plan['solver_runtime']=='editable' and plan['python']==bound['python']
        assert plan['source']==bound['source'] and plan['source_commit']=='a'*40
        assert plan['isolated_wheel_qa_required'] is True


@pytest.mark.parametrize('targets',['validation,100','all'])
def test_phase8_checks_all_bindings_before_first_mutation(monkeypatch,targets):
    """Guards: a later invalid target cannot leave partial installs. Assumes: mocked actions."""
    from argparse import Namespace
    seen=[]
    monkeypatch.setattr(release_quad,'_release_head',lambda:'a'*40)
    monkeypatch.setattr(release_quad,'cmd_preflight',lambda a:seen.append('preflight') or 0)
    monkeypatch.setattr(release_quad,'_deploy_validation',lambda:seen.append('validation') or 0)
    monkeypatch.setattr(release_quad,'_deploy_100',lambda:seen.append('100') or 0)
    assert release_quad.cmd_phase8(Namespace(target=targets))==0
    assert seen==['preflight','validation','100']
    seen.clear();monkeypatch.setenv('RADIA_RELEASE_EDITABLE_COMMIT_100','b'*40)
    with pytest.raises(ValueError,match='release tag'):release_quad.cmd_phase8(Namespace(target=targets))
    assert seen==[]


@pytest.mark.parametrize('role',['100','validation'])
@pytest.mark.parametrize('installer_rc',[0,3])
def test_role_deployment_checks_identity_and_installer_before_intent(monkeypatch,role,installer_rc):
    """Guards: install failure never records success. Assumes: observed remote transport."""
    calls=[]
    monkeypatch.setattr(release_quad,'_check_remote_target',lambda target:calls.append(('target',target['role'])) or 0)
    monkeypatch.setattr(release_quad,'run',lambda cmd,**kw:calls.append(('install',cmd)) or subprocess.CompletedProcess(cmd,installer_rc))
    monkeypatch.setattr(release_quad,'_record_release_intent_remote',lambda actual,reason:calls.append(('record',actual)) or 0)
    monkeypatch.setattr(release_quad,'_verify_role_editable',lambda actual:calls.append(('verify',actual)) or 0)
    assert release_quad._deploy_editable_role(role)==(3 if installer_rc else 0)
    assert calls[0]==('target',role)
    assert calls[1][1][:2]==['ssh',release_quad._editable_target(role)['ssh']]
    assert [c[0] for c in calls]==(['target','install'] if installer_rc else ['target','install','record','verify'])


@pytest.mark.parametrize('role',['100','validation'])
def test_installer_program_checks_before_pip_and_changes_only_bound_radia(monkeypatch,role):
    """Guards: identity/native checks precede a single editable install. Assumes: Python control flow is executed with observed subprocesses."""
    target=release_quad._editable_target(role);seen=[]
    monkeypatch.setattr(release_quad,'_target_check_code',lambda t:"seen.append('identity')")
    monkeypatch.setattr(release_quad.subprocess,'_contract_event',seen.append,raising=False)
    monkeypatch.setattr(release_quad,'EDITABLE_RELEASE_VERIFY',"import subprocess,sys; subprocess._contract_event('native-pre' if 'preinstall' in sys.argv else 'native-post')")
    monkeypatch.setattr(release_quad,'SOLVER_INSTALL_GUARD',"seen.append('lock')")
    monkeypatch.setattr(release_quad.subprocess,'run',lambda cmd,**kw:seen.append(cmd))
    scope={'seen':seen};exec(release_quad._editable_install_program(target),scope)
    assert seen[:3]==['identity','native-pre','lock']
    assert seen[3][:5]==[target['python'],'-I','-m','pip','install']
    assert seen[3][-2:]==['-e',target['source']+'[gpu]']
    assert seen[4:6]==['identity','native-post']
    assert seen[-1]==[target['python'],'-I','-m','pip','check']


@pytest.mark.parametrize('failure',[None,0,1,2,3])
def test_final_pip_checks_cover_four_explicit_interpreters_and_fail_fast(monkeypatch,failure):
    """Guards: all normal and QA dependencies checked. Assumes: exit codes are synthetic."""
    calls=[]
    def run(cmd,**kw):
        rc=int(len(calls)==failure);calls.append(cmd)
        return subprocess.CompletedProcess(cmd,rc)
    monkeypatch.setattr(release_quad,'run',run)
    assert release_quad._verify_final_pip_checks()==(0 if failure is None else 4)
    expected=[]
    for role in release_quad.SOLVER_TARGETS:
        target=release_quad._editable_target(role);expected.append((target['ssh'],target['python']))
    expected += [(host,python) for _,host,python in (release_quad._solver_qa_target(role) for role in release_quad.SOLVER_TARGETS)]
    assert len(calls)==(4 if failure is None else failure+1)
    for cmd,(host,python) in zip(calls,expected):
        assert cmd[:2]==['ssh',host]
        script=base64.b64decode(cmd[-1]).decode('utf-16le')
        assert release_quad._runtime_guard_powershell(python) in script
        assert '& '+release_quad._powershell_literal(python)+' -I -m pip check' in script
        assert 'exit $LASTEXITCODE' in script


@pytest.mark.parametrize('role',['100','validation'])
@pytest.mark.parametrize('rc',[0,1])
def test_isolated_qa_verifier_propagates_payload_check_exit(monkeypatch,role,rc):
    """Guards: QA is selected separately and failures propagate. Assumes: no SSH executed."""
    calls=[]
    monkeypatch.setattr(release_quad,'_check_remote_target',lambda target,**kw:calls.append(('identity',target,kw)) or 0)
    monkeypatch.setattr(release_quad,'_verify_wheel_runtime',lambda host,label,**kw:calls.append(('wheel',host,label,kw)) or rc)
    assert release_quad._verify_qa_wheel(role)==rc
    label,host,python=release_quad._solver_qa_target(role)
    assert calls[0][2]==dict(python=python,source=False,qa=True)
    assert calls[1]==('wheel',host,label+' QA',dict(python_path=python,isolated=True))


def test_done_requires_both_normal_and_both_qa_views_before_phase9(monkeypatch):
    """Guards: no missing role/QA completion. Assumes: earlier gates are observed mocks."""
    from argparse import Namespace
    calls=[]
    monkeypatch.setattr(release_quad,'cmd_preflight',lambda a:0)
    monkeypatch.setattr(release_quad,'_release_head',lambda:'a'*40)
    monkeypatch.setattr(release_quad,'_verify_local_release_source',lambda *a:0)
    monkeypatch.setattr(release_quad,'_verify_head_release_tag',lambda:0)
    monkeypatch.setattr(release_quad,'_verify_role_editable',lambda role:calls.append(('normal',role)) or 0)
    monkeypatch.setattr(release_quad,'_verify_qa_wheel',lambda role:calls.append(('qa',role)) or 0)
    monkeypatch.setattr(release_quad,'cmd_phase9',lambda a:4)
    assert release_quad.cmd_done(Namespace(simulink_package=None))==4
    assert set(calls)=={('normal','100'),('normal','validation'),('qa','100'),('qa','validation')} and len(calls)==4


@pytest.mark.parametrize('defect',['hostname','python','commit','dirty'])
def test_real_target_identity_guard_rejects_before_installer(tmp_path,release_template,defect):
    """Guards: target identity and tracked source are real preconditions. Assumes: disposable local Git source and current Python emulate a target."""
    import platform
    repo,sha=_release_repo(tmp_path,release_template)
    _git(repo,'checkout','-q','--detach',sha)
    target={'binding':{'hostname':platform.node()},'python':sys.executable,'source':str(repo),'commit':sha}
    ok=subprocess.run([sys.executable,'-I','-c',release_quad._target_check_code(target)],capture_output=True,text=True,timeout=30)
    assert ok.returncode==0,ok.stderr
    expected={'hostname':'Selected hostname mismatch','python':'Selected Python mismatch','commit':'Selected source commit mismatch','dirty':'Selected source has tracked changes'}[defect]
    if defect=='hostname':target['binding']['hostname']='wrong-host-for-test'
    elif defect=='python':target['python']=str(tmp_path/'missing-python.exe')
    elif defect=='commit':target['commit']='0'*40
    else:(repo/'tool.py').write_text('dirty\n',encoding='utf-8')
    bad=subprocess.run([sys.executable,'-I','-c',release_quad._target_check_code(target)],capture_output=True,text=True,timeout=30)
    assert bad.returncode!=0 and expected in bad.stderr,(bad.returncode,bad.stdout,bad.stderr)

@pytest.mark.parametrize('defect',[None,'missing-peec','missing-motor','wrong-source','dirty','wrong-schema','wrong-binary-hash'])
def test_editable_native_gate_checks_complete_cpu_payload_and_source_receipts(tmp_path,defect):
    """Guards: hash-only manifests cannot conceal missing/stale CPU binaries.

    Assumes: synthetic payload bytes; the independently tested GPU verifier is
    replaced with an empty module list to isolate the CPU provenance contract.
    """
    import json
    repo=tmp_path/'source';package=repo/'src/radia';package.mkdir(parents=True)
    tools=repo/'tools';tools.mkdir()
    (tools/'python_gpu_payload.py').write_text('MODULES=()\n',encoding='utf-8')
    sha='a'*40
    names=['_radia_pybind.pyd','axifem.pyd','peec_matrices.pyd','sparsesolv_ngsolve.pyd','radia_motor_rom.dll']
    hashes={}
    for name in names:
        payload=(name+' synthetic bytes').encode();(package/name).write_bytes(payload)
        hashes[name]=hashlib.sha256(payload).hexdigest()
        if name.endswith('.pyd'):
            receipt=dict(schema='radia.native-build-provenance.v1',source_commit=sha,source_dirty=False,binary_name=name,binary_sha256=hashes[name],binary_bytes=len(payload))
            if name=='_radia_pybind.pyd':
                if defect=='wrong-source':receipt['source_commit']='b'*40
                elif defect=='dirty':receipt['source_dirty']=True
                elif defect=='wrong-schema':receipt['schema']='other'
                elif defect=='wrong-binary-hash':receipt['binary_sha256']='0'*64
            (package/(name+'.build.json')).write_text(json.dumps(receipt),encoding='utf-8')
    if defect in ('missing-peec','missing-motor'):
        name='peec_matrices.pyd' if defect=='missing-peec' else 'radia_motor_rom.dll'
        (package/name).unlink();hashes.pop(name)
    (package/'release_native_payloads.json').write_text(json.dumps(dict(schema='radia.release-native-payloads.v1',source_commit=sha,files=hashes)),encoding='utf-8')
    result=subprocess.run([sys.executable,'-I','-c',release_quad.EDITABLE_RELEASE_VERIFY,str(repo),sha,'test-ng','test-netgen','preinstall'],capture_output=True,text=True,timeout=30)
    assert (result.returncode==0)==(defect is None),(defect,result.returncode,result.stdout,result.stderr)
    if defect is not None:
        assert ('omits required' if defect.startswith('missing-') else 'native source provenance differs') in result.stderr
