from __future__ import annotations

import base64
import hashlib
import importlib.util
import inspect
import json
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
    for role in ("VALIDATION",):
        host="validation-test"
        for key,value in {"SSH":host,"HOSTNAME":host,"PYTHON":f"C:/Runtime/{host}/python.exe"}.items():
            monkeypatch.setenv(f"RADIA_RELEASE_{key}_{role}",value)


def _git(repo, *args):
    return subprocess.run(
        [_GIT, "-C", str(repo), *args],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


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


def test_deployment_plan_names_four_lanes_and_validation_runtime_only(monkeypatch,capsys):
    """Guards: dry planning never installs; GPU lanes never run on the CPU build role. Assumes: no host lookup."""
    import json
    monkeypatch.setattr(release_quad,'run',lambda *a,**k:pytest.fail('planning mutated a host'))
    monkeypatch.setattr(release_quad,'_wheel_target',lambda *a:pytest.fail('planning resolved a host'))
    assert release_quad.cmd_deployment_plan(None)==0
    plan=json.loads(capsys.readouterr().out)
    lanes={row['lane']:(row['build_role'],row['validation_role']) for row in plan['lanes']}
    assert lanes==release_quad.RELEASE_LANES
    assert set(lanes)=={'cpu','linux-python-gpu','windows-python-gpu','matlab'}
    assert all(lanes[name][0]!='cpu-build' for name in lanes if name!='cpu'),lanes
    assert plan['solver_runtime_roles']==['validation'] and plan['acceptance_roles']==['validation']
    assert plan['published_wheel_lane']=='windows-python-gpu'


# ---------------------------------------------------------------------------
# Verified-wheel deployment on the validation role (no editable route)
# ---------------------------------------------------------------------------

def _record_rows(files):
    rows = []
    for relative, content in files.items():
        digest = base64.urlsafe_b64encode(hashlib.sha256(content).digest())
        rows.append(f"{relative},sha256={digest.rstrip(b'=').decode()},{len(content)}")
    return rows


def _fake_wheel_site(root, *, shadow=False, editable=False):
    """Create a tiny Radia install with a valid hash-checked RECORD and its source wheel."""
    import zipfile
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
    wheel = root / "radia-5.1.0-py3-none-any.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        for relative, content in files.items():
            archive.writestr(relative, content)
        archive.writestr("radia-5.1.0.dist-info/RECORD",
                         "\n".join(_record_rows(files) + ["radia-5.1.0.dist-info/RECORD,,"]) + "\n")
    if editable:
        files["radia-5.1.0.dist-info/direct_url.json"] = json.dumps(
            {"url": "file:///C:/Sources/radia", "dir_info": {"editable": True}}).encode()
    for relative, content in files.items():
        (site / relative).write_bytes(content)
    rows = _record_rows(files) + ["radia-5.1.0.dist-info/RECORD,,"]
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
    return env, wheel


@pytest.mark.parametrize(
    "expected_version,shadow,wrong_root,editable,other_wheel,expected_rc,expected_message",
    [
        ("5.1.0", False, False, False, False, 0, "validation wheel verified"),
        ("5.1.1", False, False, False, False, 1, "Radia version mismatch"),
        ("5.1.0", True, False, False, False, 1, "Shadowed Radia import"),
        ("5.1.0", False, True, False, False, 1, "Radia wheel must be under site-packages"),
        ("5.1.0", False, False, True, False, 1, "must use a verified wheel, not an editable install"),
        ("5.1.0", False, False, False, True, 1, "Installed files differ from the candidate wheel"),
    ],
    ids=["exact", "version", "shadow", "root", "editable", "other-wheel"],
)
def test_wheel_runtime_probe_binds_the_install_to_the_candidate_wheel(
    tmp_path, expected_version, shadow, wrong_root, editable, other_wheel, expected_rc, expected_message
):
    """Guards: the validation runtime is the exact non-editable candidate wheel.

    Assumes: a synthetic site and wheel; the production probe runs unchanged.
    """
    import zipfile
    env, wheel = _fake_wheel_site(tmp_path, shadow=shadow, editable=editable)
    if other_wheel:
        changed = {"radia/__init__.py": b"__version__ = '5.1.0'\n# other build\n"}
        wheel = tmp_path / "other" / wheel.name
        wheel.parent.mkdir()
        with zipfile.ZipFile(wheel, "w") as archive:
            archive.writestr("radia/__init__.py", changed["radia/__init__.py"])
            archive.writestr("radia-5.1.0.dist-info/RECORD", "\n".join(_record_rows(changed)) + "\n")
    code = release_quad._wheel_runtime_probe_code(expected_version, "validation", str(wheel))
    install_root = str(tmp_path / ("other-site" if wrong_root else "site"))
    setup = (
        "import sysconfig\n"
        f"sysconfig.get_paths = lambda: {{'purelib': {install_root!r}, 'platlib': {install_root!r}}}\n"
    )
    result = subprocess.run([sys.executable, "-c", setup + code],
                            capture_output=True, text=True, env=env, timeout=30)
    assert result.returncode == expected_rc, (result.returncode, result.stdout, result.stderr)
    assert expected_message in result.stdout + result.stderr


def _phase8_args(**updates):
    from argparse import Namespace
    values = dict(target="validation", source_sha="a" * 40, manifest="manifest.json")
    values.update(updates)
    return Namespace(**values)


@pytest.mark.parametrize("missing", ["RADIA_RELEASE_SSH_VALIDATION", "RADIA_RELEASE_PYTHON_VALIDATION",
                                     "RADIA_RELEASE_HOSTNAME_VALIDATION"])
def test_validation_config_fails_before_any_phase8_action(monkeypatch, missing):
    """Guards: no fallback or partial deployment. Assumes: remote actions are observed mocks."""
    monkeypatch.delenv(missing)
    monkeypatch.setattr(release_quad, "run", lambda *a, **kw: pytest.fail("No remote action"))
    monkeypatch.setattr(release_quad, "cmd_preflight", lambda *a: pytest.fail("No preflight before config"))
    with pytest.raises(ValueError, match=missing):
        release_quad.cmd_phase8(_phase8_args())


def test_phase8_rejects_the_build_host_and_unbound_candidates_before_mutation(monkeypatch, tmp_path):
    """Guards: only the exact tagged candidate wheel reaches the validation role. Assumes: mocked actions."""
    seen = []
    wheel = tmp_path / "radia-9.9.9-cp312-cp312-win_amd64.whl"
    monkeypatch.setattr(release_quad, "_release_head", lambda: "a" * 40)
    monkeypatch.setattr(release_quad, "cmd_preflight", lambda a: seen.append("preflight") or 0)
    monkeypatch.setattr(release_quad, "_release_wheel",
                        lambda manifest, sha: seen.append(("wheel", str(manifest), sha)) or (wheel, "f" * 64))
    monkeypatch.setattr(release_quad, "_fork_wheel_lock", lambda sha: {"sha": sha})
    monkeypatch.setattr(release_quad, "_deploy_wheel",
                        lambda role, path, digest, lock: seen.append(("deploy", role, path, digest, lock)) or 0)
    for target in ("100", "all", "validation,100"):
        assert release_quad.cmd_phase8(_phase8_args(target=target)) == 2
    for sha in ("b" * 40, "A" * 40, "a" * 39, None):
        assert release_quad.cmd_phase8(_phase8_args(source_sha=sha)) == 4
    assert seen == []
    assert release_quad.cmd_phase8(_phase8_args()) == 0
    assert seen == ["preflight", ("wheel", "manifest.json", "a" * 40),
                    ("deploy", "validation", wheel, "f" * 64, {"sha": "a" * 40})]
    seen.clear()
    monkeypatch.setattr(release_quad, "_release_wheel",
                        lambda *a: (_ for _ in ()).throw(ValueError("SHA-256 mismatch")))
    assert release_quad.cmd_phase8(_phase8_args()) == 4
    assert seen == ["preflight"]


@pytest.mark.parametrize("defect", [None, "changed", "none", "two", "wrong-version"])
def test_release_wheel_is_the_one_hash_verified_published_lane_wheel(monkeypatch, tmp_path, defect):
    """Guards: phase8 installs only the receipt-hashed published wheel. Assumes: receipts are mocked."""
    monkeypatch.setattr(release_quad, "_radia_version", lambda: "9.9.9")
    lane = tmp_path / "windows-python-gpu"
    lane.mkdir()
    names = {"two": ["radia-9.9.9-cp312-cp312-win_amd64.whl", "radia-9.9.9-cp313-cp313-win_amd64.whl"],
             "none": [], "wrong-version": ["radia-9.9.8-cp312-cp312-win_amd64.whl"]}.get(
        defect, ["radia-9.9.9-cp312-cp312-win_amd64.whl"])
    validated = {"windows-python-gpu/receipt.json": hashlib.sha256(b"{}").hexdigest()}
    (lane / "receipt.json").write_bytes(b"{}")
    for name in names:
        (lane / name).write_bytes(name.encode())
        validated["windows-python-gpu/" + name] = hashlib.sha256(name.encode()).hexdigest()
    if defect == "changed":
        (lane / names[0]).write_bytes(b"replaced after receipt verification")
    seen = []
    monkeypatch.setattr(release_quad, "_verify_candidate_receipts",
                        lambda manifest, sha: seen.append(sha) or {"windows-python-gpu": validated})
    manifest = tmp_path / "manifest.json"
    if defect is None:
        path, digest = release_quad._release_wheel(manifest, "a" * 40)
        assert path == (lane / names[0]).resolve() and digest == validated["windows-python-gpu/" + names[0]]
    else:
        with pytest.raises(ValueError, match="SHA-256 mismatch" if defect == "changed" else "exactly one"):
            release_quad._release_wheel(manifest, "a" * 40)
    assert seen == ["a" * 40]


def _lock():
    return json.loads((_TOOL.parent / "ci_fork_wheels.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize("staged_ok", [True, False])
def test_installer_hash_checks_staged_wheels_and_never_installs_editable(monkeypatch, tmp_path, staged_ok):
    """Guards: hash and identity checks precede pip; the exact Radia wheel and locked
    fork pair are installed, never -e. Assumes: observed subprocesses.
    """
    monkeypatch.setattr(release_quad, "_radia_version", lambda: "9.9.9")
    wheel = tmp_path / "radia-9.9.9-cp312-cp312-win_amd64.whl"
    fork = tmp_path / "ngsolve-fork.whl"
    wheel.write_bytes(b"candidate wheel")
    fork.write_bytes(b"fork wheel")
    staged = {str(wheel): hashlib.sha256(b"candidate wheel").hexdigest(),
              str(fork): hashlib.sha256(b"fork wheel" if staged_ok else b"other").hexdigest()}
    target = release_quad._wheel_target("validation")
    seen = []
    monkeypatch.setattr(release_quad, "_target_check_code", lambda t: "seen.append('identity')\n")
    monkeypatch.setattr(release_quad, "SOLVER_INSTALL_GUARD", "seen.append('lock')")

    def observe(cmd, **kwargs):
        assert kwargs.get("check") is True
        seen.append(("fresh", cmd[0], "verify_installed" in cmd[3]) if cmd[1:3] == ["-I", "-c"] else cmd)
    monkeypatch.setattr(release_quad.subprocess, "run", observe)
    program = release_quad._wheel_install_program(target, staged, _lock())
    if not staged_ok:
        with pytest.raises(SystemExit, match="SHA-256 differs"):
            exec(program, {"seen": seen})
        assert seen == ["identity"]
        return
    exec(program, {"seen": seen})
    python = target["python"]
    assert seen == [
        "identity", "lock",
        [python, "-I", "-m", "pip", "install", "--no-cache-dir", "--only-binary=:all:",
         str(wheel) + "[gpu]", str(fork)],
        [python, "-I", "-m", "pip", "install", "--no-cache-dir", "--only-binary=:all:",
         "--force-reinstall", "--no-deps", str(fork), str(wheel)],
        ("fresh", python, True),
        [python, "-I", "-m", "pip", "check"],
    ]
    installs = [cmd for cmd in seen if isinstance(cmd, list) and "install" in cmd]
    assert all("-e" not in cmd and "--editable" not in cmd for cmd in installs), installs


@pytest.mark.parametrize("defect", [None, "stock", "wrong-base", "missing-dlls"])
def test_fork_lock_comes_from_the_deployment_source_and_names_a_radia_build(monkeypatch, defect):
    """Guards: stock or mismatched NGSolve/Netgen pins cannot enter the lock. Assumes: git show is observed."""
    lock = _lock()
    if defect == "stock":
        lock["wheels"]["ngsolve"]["version"] = "6.2.2608"
    elif defect == "wrong-base":
        lock["wheels"]["netgen-mesher"]["version"] = "6.2.2607+radia.b2894d9f"
    elif defect == "missing-dlls":
        lock["dlls"] = {}
    calls = []
    monkeypatch.setattr(release_quad, "_git", lambda *a, **k: calls.append(a) or
                        subprocess.CompletedProcess(a, 0, json.dumps(lock), ""))
    if defect is None:
        assert release_quad._fork_wheel_lock("c" * 40) == lock
    else:
        with pytest.raises(ValueError):
            release_quad._fork_wheel_lock("c" * 40)
    assert calls == [("show", "c" * 40 + ":tools/ci_fork_wheels.json")]


@pytest.mark.parametrize("defect", [None, "stock", "wrong-build", "dll"])
def test_installed_fork_check_accepts_only_the_locked_pair(tmp_path, defect):
    """Guards: the target-side check rejects stock or wrong-build NGSolve/Netgen and
    foreign DLLs. Assumes: a synthetic prefix with fake packages and metadata.
    """
    lock = _lock()
    site = tmp_path / "site"
    for package in ("netgen", "ngsolve"):
        (site / package).mkdir(parents=True)
        (site / package / "__init__.py").write_text("", encoding="utf-8")
    dll_bytes = {name: name.encode() for name in lock["dlls"]}
    lock["dlls"] = {name: hashlib.sha256(data).hexdigest() for name, data in dll_bytes.items()}
    for name, data in dll_bytes.items():
        (site / "netgen" / name).write_bytes(data + (b"x" if defect == "dll" else b""))
    (site / "radia").mkdir()
    (site / "radia" / "__init__.py").write_text("", encoding="utf-8")
    (site / "radia" / "_cuda_runtime.py").write_text("def _pins():\n    return {}\n", encoding="utf-8")
    versions = {name: wheel["version"] for name, wheel in lock["wheels"].items()}
    if defect == "stock":
        versions["ngsolve"] = "6.2.2608"
    elif defect == "wrong-build":
        versions["netgen-mesher"] = "6.2.2608+radia.00000000"
    setup = (f"import sys, importlib.metadata\nsys.prefix = {str(tmp_path)!r}\nsys.path.insert(0, {str(site)!r})\n"
             f"importlib.metadata.version = {versions!r}.__getitem__\n")
    result = subprocess.run([sys.executable, "-I", "-c", setup + release_quad._solver_pin_check_code(lock)],
                            capture_output=True, text=True, timeout=30)
    assert (result.returncode == 0) == (defect is None), (defect, result.stdout, result.stderr)


def _fake_lock(tmp_path):
    lock = _lock()
    paths = []
    for wheel in lock["wheels"].values():
        path = tmp_path / wheel["filename"]
        path.write_bytes(b"fork")
        paths.append(path)
    return lock, paths


@pytest.mark.parametrize("failure", [None, "identity", "mkdir", "scp", "install", "cleanup"])
def test_deploy_wheel_owns_its_staging_and_reports_cleanup_failure(monkeypatch, tmp_path, failure):
    """Guards: exclusive job-owned staging, per-file cleanup only, and failures
    (including cleanup) are nonzero. Assumes: observed transport.
    """
    wheel = tmp_path / "radia-9.9.9-cp312-cp312-win_amd64.whl"
    wheel.write_bytes(b"w")
    monkeypatch.setattr(release_quad, "_radia_version", lambda: "9.9.9")
    lock, forks = _fake_lock(tmp_path)
    monkeypatch.setattr(release_quad, "_fetch_fork_wheels", lambda lk, scratch: forks)
    calls, scripts = [], []
    monkeypatch.setattr(release_quad, "_check_remote_target",
                        lambda target: calls.append("identity") or int(failure == "identity"))

    def ssh(host, script):
        kind = ("mkdir" if "New-Item" in script else "cleanup" if "Remove-Item" in script else "install")
        calls.append(kind)
        scripts.append((kind, script))
        return int(failure == kind)
    monkeypatch.setattr(release_quad, "_ssh_pwsh", ssh)
    monkeypatch.setattr(release_quad, "run", lambda cmd, **kw: calls.append("scp") or
                        subprocess.CompletedProcess(cmd, int(failure == "scp")))
    rc = release_quad._deploy_wheel("validation", wheel, "f" * 64, lock)
    assert rc == (0 if failure is None else 3)
    expected = {"identity": ["identity"], "mkdir": ["identity", "mkdir"],
                "scp": ["identity", "mkdir", "scp", "cleanup"]}.get(
        failure, ["identity", "mkdir", "scp", "scp", "scp", "install", "cleanup"])
    assert calls == expected
    mkdir = dict(scripts).get("mkdir")
    if mkdir:
        assert "-Force" not in mkdir and "exit 1" in mkdir
    cleanup = dict(scripts).get("cleanup")
    if cleanup:
        assert "-Recurse" not in cleanup
        assert cleanup.count("Remove-Item") == (1 if failure == "scp" else 3)
        assert "[IO.Directory]::Delete(" in cleanup and ", $false)" in cleanup


def test_same_wheel_retries_never_share_a_staging_directory(monkeypatch, tmp_path):
    """Guards: two deployments of one wheel use distinct directories. Assumes: observed transport."""
    wheel = tmp_path / "radia-9.9.9-cp312-cp312-win_amd64.whl"
    wheel.write_bytes(b"w")
    monkeypatch.setattr(release_quad, "_radia_version", lambda: "9.9.9")
    lock, forks = _fake_lock(tmp_path)
    monkeypatch.setattr(release_quad, "_fetch_fork_wheels", lambda lk, scratch: forks)
    monkeypatch.setattr(release_quad, "_check_remote_target", lambda target: 0)
    created = []
    monkeypatch.setattr(release_quad, "_ssh_pwsh",
                        lambda host, script: created.append(script) if "New-Item" in script else 0)
    monkeypatch.setattr(release_quad, "run", lambda cmd, **kw: subprocess.CompletedProcess(cmd, 0))
    for _ in range(2):
        assert release_quad._deploy_wheel("validation", wheel, "f" * 64, lock) == 0
    assert len(created) == 2 and created[0] != created[1]


def test_validation_role_has_no_editable_route():
    """Guards: the editable validation route stays deleted (latest-only policy)."""
    for name in ("cmd_verify_editable", "cmd_repoint", "cmd_restore_editable", "_editable_target",
                 "_deploy_editable_role", "_editable_install_program", "_record_release_intent_remote",
                 "_solver_qa_target", "editable_intent", "EDITABLE_RELEASE_VERIFY"):
        assert not hasattr(release_quad, name), name
    assert not (_TOOL.parent / "editable_intent.py").exists()
    result = subprocess.run([sys.executable, str(_TOOL), "--help"], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    for command in ("verify-editable", "repoint", "restore-editable"):
        assert command not in result.stdout, result.stdout
    assert "deployment-plan" in result.stdout


@pytest.mark.parametrize("defect", ["hostname", "python"])
def test_real_target_identity_guard_rejects_before_installer(tmp_path, defect):
    """Guards: target identity is a real precondition. Assumes: the current Python emulates a target."""
    import platform
    target = {"hostname": platform.node(), "python": sys.executable}
    good = subprocess.run([sys.executable, "-I", "-c", release_quad._target_check_code(target)],
                          capture_output=True, text=True, timeout=30)
    assert good.returncode == 0, good.stderr
    if defect == "hostname":
        target["hostname"] = "wrong-host-for-test"
    else:
        target["python"] = str(tmp_path / "missing-python.exe")
    bad = subprocess.run([sys.executable, "-I", "-c", release_quad._target_check_code(target)],
                         capture_output=True, text=True, timeout=30)
    expected = {"hostname": "Selected hostname mismatch", "python": "Selected Python mismatch"}[defect]
    assert bad.returncode != 0 and expected in bad.stderr, (bad.returncode, bad.stdout, bad.stderr)

@pytest.mark.skipif(shutil.which("pwsh") is None, reason="PowerShell 7 is required to execute the cleanup script")
def test_cleanup_keeps_unexpected_files_and_fails_without_prompting(tmp_path):
    """Guards: staging cleanup deletes only placed files, never prompts, fails on leftovers. Assumes: local pwsh."""
    staging = tmp_path / "stage"
    staging.mkdir()
    placed = staging / "radia.whl"
    placed.write_bytes(b"w")
    sentinel = staging / "other-owner.txt"
    sentinel.write_bytes(b"keep")
    literal = release_quad._powershell_literal(str(staging))
    script = (f"if (Test-Path -LiteralPath {release_quad._powershell_literal(str(placed))}) "
              f"{{ Remove-Item -LiteralPath {release_quad._powershell_literal(str(placed))} -Force }}; "
              f"try {{ [IO.Directory]::Delete({literal}, $false) }} catch {{ exit 1 }}; "
              f"if (Test-Path -LiteralPath {literal}) {{ exit 1 }}; exit 0")
    source = inspect.getsource(release_quad._deploy_wheel)
    assert "[IO.Directory]::Delete({literal}, $false)" in source, "production cleanup changed; update this replica"
    proc = subprocess.run(["pwsh", "-NoProfile", "-Command", script], stdin=subprocess.DEVNULL,
                          capture_output=True, text=True, timeout=60)
    assert proc.returncode == 1, (proc.returncode, proc.stdout, proc.stderr)
    assert not placed.exists() and sentinel.read_bytes() == b"keep"
    sentinel.unlink()
    placed.write_bytes(b"w")
    proc = subprocess.run(["pwsh", "-NoProfile", "-Command", script], stdin=subprocess.DEVNULL,
                          capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0 and not staging.exists(), (proc.returncode, proc.stderr)
