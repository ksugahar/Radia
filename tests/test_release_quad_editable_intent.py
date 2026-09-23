"""release_quad follows the recorded editable intent: verify, done, repoint, Phase 8."""
from __future__ import annotations

import importlib.util
import inspect
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

_TOOL = Path(__file__).resolve().parents[1] / "tools" / "release_quad.py"
_SPEC = importlib.util.spec_from_file_location("radia_release_quad_intent", _TOOL)
assert _SPEC is not None and _SPEC.loader is not None
release_quad = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(release_quad)
intent = release_quad.editable_intent


@pytest.fixture(autouse=True)
def isolated_record(tmp_path, monkeypatch):
    """Never read the machine's record or the release overrides from the environment."""
    path = tmp_path / "editable-intent.json"
    monkeypatch.setenv(intent.INTENT_FILE_ENV, str(path))
    monkeypatch.delenv(release_quad.EDITABLE_REPO_LAB_ENV, raising=False)
    monkeypatch.delenv(release_quad.EDITABLE_REPO_100_ENV, raising=False)
    return path


def _record(path, package, source):
    data = intent.load_intent(path)
    intent.set_entry(data, package, {"source": source, "commit": None, "tracked_clean": None,
                                     "recorded_at": "2026-09-12T00:00:00Z", "recorded_by": "test",
                                     "recorded_via": "test", "reason": "fixture",
                                     "pushed_refs": None, "previous": None})
    intent.save_intent(data, path)


def test_remote_verify_ships_the_intent_tool_as_the_script():
    assert release_quad.REMOTE_EDITABLE_VERIFY == (
        _TOOL.with_name("editable_intent.py").read_text(encoding="utf-8"))
    assert "//192.168.121.100/work/00_cae/radia/" in release_quad.REMOTE_EDITABLE_VERIFY


def test_expected_lab_packages_come_from_the_record_unless_overridden(isolated_record, monkeypatch):
    assert release_quad._expected_lab_editable_packages() == [
        ("radia", None), ("cubit-mesh-export", None), ("radia-mcp", None),
        ("mcp-server-document", None),
    ]
    _record(isolated_record, "radia-mcp", "S:/Radia/release-quad/x/packages/radia-mcp")
    assert dict(release_quad._expected_lab_editable_packages())["radia-mcp"] == (
        "S:/Radia/release-quad/x/packages/radia-mcp")

    monkeypatch.setenv(release_quad.EDITABLE_REPO_LAB_ENV, "S:/Radia/release-quad/override")
    assert dict(release_quad._expected_lab_editable_packages())["radia"] == (
        "S:/Radia/release-quad/override")


def test_unrecorded_package_is_unverified_not_drift_and_gets_no_target(monkeypatch, capsys):
    monkeypatch.setattr(release_quad, "_pip_show", lambda name: {
        "version": "1.4.53", "location": "site-packages",
        "editable project location": "S:/Radia/release-quad/x/packages/radia-mcp"})
    monkeypatch.setattr(release_quad, "_fresh_import_origin",
                        lambda name: "S:/Radia/release-quad/x/packages/radia-mcp/src/radia_mcp/__init__.py")
    report = {}
    assert release_quad._verify_lab_editable([("radia-mcp", None)], report) == 1
    assert report == {"ok": 0, "drift": 0, "unverified": 1, "missing": 0}
    out = capsys.readouterr().out
    assert "UNVERIFIED" in out
    assert "record-current" in out
    for forbidden in ("01_GitHub", "pip install -e", "uninstall", "Stop-Process", "2026-05-27"):
        assert forbidden not in out


def test_drift_hint_names_repoint_and_no_default_target(monkeypatch, capsys):
    monkeypatch.setattr(release_quad, "_pip_show", lambda name: {
        "version": "1.4.53", "location": "site-packages",
        "editable project location": "S:/Radia/release-quad/x/packages/radia-mcp"})
    monkeypatch.setattr(release_quad, "_fresh_import_origin", lambda name: None)
    report = {}
    assert release_quad._verify_lab_editable(
        [("radia-mcp", "S:/Radia/release-quad/y/packages/radia-mcp")], report) == 1
    assert report["drift"] == 1 and report["unverified"] == 0
    out = capsys.readouterr().out
    assert "DRIFT" in out
    assert "repoint --package <pkg> --source <intended>" in out
    assert "No default target" in out
    for forbidden in ("01_GitHub", "pip uninstall", "uninstall -y", "Stop-Process"):
        assert forbidden not in out


def test_verify_editable_exit_codes_distinguish_unverified_from_drift(monkeypatch):
    def lab(counts):
        def fake(packages=None, report=None):
            if report is not None:
                report.update(counts)
            return counts["drift"] + counts["unverified"]
        return fake

    monkeypatch.setattr(release_quad, "_verify_lab_editable", lab({"drift": 0, "unverified": 1}))
    monkeypatch.setattr(release_quad, "_verify_100_editable", lab({"drift": 0, "unverified": 0}))
    assert release_quad.cmd_verify_editable(None) == 5

    monkeypatch.setattr(release_quad, "_verify_lab_editable", lab({"drift": 0, "unverified": 0}))
    monkeypatch.setattr(release_quad, "_verify_100_editable", lab({"drift": 1, "unverified": 0}))
    assert release_quad.cmd_verify_editable(None) == 1

    monkeypatch.setattr(release_quad, "_verify_100_editable", lab({"drift": 0, "unverified": 0}))
    assert release_quad.cmd_verify_editable(None) == 0


def _remote_ok(counts):
    return (0, {"schema": intent.SCHEMA, "host": "100", "interpreter": "python",
                "record_file": "C:/ProgramData/Radia/editable-intent.json", "rows": [],
                "counts": counts, "exit_code": 0}, "")


def test_remote_verify_passes_explicit_expectations_only_with_the_release_override(monkeypatch):
    seen = []
    monkeypatch.setattr(release_quad, "_remote_editable_intent",
                        lambda host, argv, **kw: seen.append((host, list(argv))) or
                        _remote_ok({"match": 3, "drift": 0, "unverified": 0}))
    assert release_quad._verify_100_editable() == 0
    assert seen[-1] == (release_quad.SSH_100, ["--json", "verify"])

    monkeypatch.setenv(release_quad.EDITABLE_REPO_100_ENV, r"W:\00_CAE\Radia\release-quad\x")
    assert release_quad._verify_100_editable() == 0
    _host, argv = seen[-1]
    assert argv[:2] == ["--json", "verify"]
    assert "--expect" in argv
    assert r"radia=W:\00_CAE\Radia\release-quad\x" in argv
    # The coupled release was retired on 2026-09-16: the solver release-quad
    # states an expectation for radia alone and neither repoints nor
    # version-gates the independently released packages.
    assert not any(arg.startswith(("radia-mcp=", "cubit-mesh-export="))
                   for arg in argv)


def test_remote_verify_counts_unverified_and_drift(monkeypatch):
    report = {}
    monkeypatch.setattr(release_quad, "_remote_editable_intent",
                        lambda host, argv, **kw: _remote_ok({"match": 2, "drift": 0, "unverified": 1}))
    assert release_quad._verify_100_editable(None, report) == 1
    assert report["unverified"] == 1 and report["drift"] == 0

    monkeypatch.setattr(release_quad, "_remote_editable_intent",
                        lambda host, argv, **kw: (1, None, "ssh: connect failed"))
    assert release_quad._verify_100_editable(None, report) == 1
    assert report["drift"] == 1


def test_done_reports_unverified_without_override_or_record(monkeypatch, capsys):
    monkeypatch.setattr(release_quad, "cmd_preflight", lambda _args: 0)

    def must_not_run(*_args, **_kwargs):
        raise AssertionError("no source can be verified before an expectation exists")

    monkeypatch.setattr(release_quad, "_verify_local_release_source", must_not_run)
    monkeypatch.setattr(release_quad, "_verify_lab_editable", must_not_run)
    args = SimpleNamespace(simulink_package=None)
    assert release_quad.cmd_done(args) == 5
    out = capsys.readouterr().out
    assert "UNVERIFIED" in out
    assert "record-current" in out
    assert "01_GitHub" not in out


def test_done_verifies_the_recorded_source_when_no_override_is_set(isolated_record, monkeypatch):
    _record(isolated_record, "radia", "S:/Radia/release-quad/recorded/")
    monkeypatch.setattr(release_quad, "cmd_preflight", lambda _args: 0)
    monkeypatch.setattr(release_quad, "_release_head", lambda: "a" * 40)
    seen = {}

    def capture(repo, sha):
        seen["repo"] = repo
        return 4

    monkeypatch.setattr(release_quad, "_verify_local_release_source", capture)
    assert release_quad.cmd_done(SimpleNamespace(simulink_package=None)) == 4
    assert seen["repo"] == "S:/Radia/release-quad/recorded"


def test_done_prefers_the_release_override_over_the_record(isolated_record, monkeypatch):
    _record(isolated_record, "radia", "S:/Radia/release-quad/recorded")
    monkeypatch.setenv(release_quad.EDITABLE_REPO_LAB_ENV, "S:/Radia/release-quad/override/")
    assert release_quad._done_active_lab_source() == "S:/Radia/release-quad/override"


def test_lab_deploy_records_the_installed_source_after_pip_install(monkeypatch):
    """The record is written after the install is verified, and only then.

    Phase 8 installs the numerical solver alone. It does not stop a process and
    does not uninstall cubit-mesh-export or radia-mcp, which release on their
    own, so this asserts their absence as well as the ordering.
    """
    events = []
    monkeypatch.setenv(release_quad.EDITABLE_REPO_LAB_ENV, "S:/Radia/release-quad/x")
    monkeypatch.setattr(release_quad, "_release_head", lambda: "a" * 40)
    monkeypatch.setattr(release_quad, "_verify_local_release_source", lambda repo, sha: 0)
    monkeypatch.setattr(release_quad, "_verify_lab_editable",
                        lambda packages: events.append(("verify", packages)) or 0)

    def fake_run(command, **kwargs):
        events.append(" ".join(str(c) for c in command))
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(release_quad, "run", fake_run)
    monkeypatch.setattr(release_quad, "_record_release_intent_lab",
                        lambda repo: events.append(("record", repo)) or 0)

    assert release_quad._deploy_lab() == 0
    commands = [e for e in events if isinstance(e, str)]
    install = next(i for i, e in enumerate(events)
                   if isinstance(e, str) and "pip install" in e)
    verify = next(i for i, e in enumerate(events)
                  if isinstance(e, tuple) and e[0] == "verify")
    record = events.index(("record", "S:/Radia/release-quad/x"))
    assert install < verify < record
    assert not any("uninstall" in c for c in commands)
    assert not any("cubit-mesh-export" in c or "radia-mcp" in c for c in commands)


def test_lab_record_helper_adopts_only_solver_without_pip(monkeypatch):
    seen = {}
    monkeypatch.setattr(release_quad, "_release_head", lambda: "b" * 40)

    def fake_record_current(packages, reason, **kwargs):
        seen.update(packages=packages, reason=reason, kwargs=kwargs)
        return {"results": [], "record_file": "r", "exit_code": 0}

    monkeypatch.setattr(intent, "record_current", fake_record_current)
    assert release_quad._record_release_intent_lab("S:/Radia/release-quad/x") == 0
    assert seen["packages"] == ["radia"]
    assert "phase8" in seen["reason"] and "bbbbbbbbbbbb" in seen["reason"]
    assert seen["kwargs"] == {"via": "release-quad phase8"}


def test_remote_record_helper_uses_record_current_with_the_host_tool(monkeypatch):
    seen = {}
    monkeypatch.setattr(release_quad, "_release_head", lambda: "c" * 40)
    monkeypatch.setattr(release_quad, "_remote_editable_intent",
                        lambda host, argv, **kw: seen.update(host=host, argv=list(argv)) or
                        (0, {"results": [], "record_file": "r", "exit_code": 0}, ""))
    assert release_quad._record_release_intent_remote("100", "100号機", r"W:\x") == 0
    assert seen["host"] == "100"
    assert seen["argv"][:3] == ["--json", "repoint", "--record-current"]
    assert seen["argv"].count("--package") == 1
    assert seen["argv"][-2:] == ["--package", "radia"]
    assert "--source" not in seen["argv"]


def test_repoint_on_100_forwards_host_namespace_arguments(monkeypatch):
    seen = {}
    monkeypatch.setattr(release_quad, "_remote_editable_intent",
                        lambda host, argv, **kw: seen.update(host=host, argv=list(argv)) or
                        (0, {"results": [], "record_file": "r", "exit_code": 0, "dry_run": True}, ""))
    args = SimpleNamespace(host="100", package=["radia-mcp"],
                           source=[r"W:\00_CAE\Radia\release-quad\x\packages\radia-mcp"],
                           reason="try the reviewed tree", via=None, record_current=False,
                           rollback=False, dry_run=True, require_pushed=False)
    assert release_quad.cmd_repoint(args) == 0
    assert seen["host"] == release_quad.SSH_100
    argv = seen["argv"]
    assert argv[:2] == ["--json", "repoint"]
    assert r"W:\00_CAE\Radia\release-quad\x\packages\radia-mcp" in argv
    assert "--dry-run" in argv and "--reason" in argv and "try the reviewed tree" in argv
    assert "--record-current" not in argv and "--rollback" not in argv


def test_repoint_on_lab_delegates_to_the_intent_tool(monkeypatch):
    seen = {}
    monkeypatch.setattr(intent, "main", lambda argv: seen.update(argv=list(argv)) or 0)
    args = SimpleNamespace(host="lab", package=["radia"], source=["S:/Radia/release-quad/x"],
                           reason="release worktree", via=None, record_current=False,
                           rollback=False, dry_run=False, require_pushed=True)
    assert release_quad.cmd_repoint(args) == 0
    assert seen["argv"][0] == "repoint"
    assert "--require-pushed" in seen["argv"]
    assert "--json" not in seen["argv"]

    undo = SimpleNamespace(host="lab", package=["radia"], source=[], reason="", via=None,
                           record_current=False, rollback=True, dry_run=False, require_pushed=False)
    assert release_quad.cmd_repoint(undo) == 0
    assert "--rollback" in seen["argv"]


def test_repoint_refuses_incomplete_requests():
    base = {"host": "lab", "via": None, "record_current": False, "rollback": False,
            "dry_run": False, "require_pushed": False}
    assert release_quad.cmd_repoint(SimpleNamespace(package=["radia"], source=[], reason="r", **base)) == 2
    assert release_quad.cmd_repoint(SimpleNamespace(package=["radia"], source=["x"], reason="", **base)) == 2
    adopt = {**base, "record_current": True}
    assert release_quad.cmd_repoint(SimpleNamespace(package=[], source=["x"], reason="r", **adopt)) == 2
    both = {**base, "record_current": True, "rollback": True}
    assert release_quad.cmd_repoint(SimpleNamespace(package=[], source=[], reason="r", **both)) == 2
    other = {**base, "host": "hibino"}
    assert release_quad.cmd_repoint(SimpleNamespace(package=["radia"], source=["x"], reason="r", **other)) == 2


def test_restore_editable_is_removed_and_names_no_target(capsys):
    assert not hasattr(release_quad, "_restore_lab_canonical_editable")
    assert not hasattr(release_quad, "_restore_100_canonical_editable")
    assert release_quad.cmd_restore_editable(None) == 2
    out = capsys.readouterr().out
    assert "repoint" in out
    for forbidden in ("01_GitHub", "uninstall -y", "Stop-Process"):
        assert forbidden not in out


def test_repoint_paths_never_uninstall_or_stop_processes():
    sources = [inspect.getsource(getattr(release_quad, name)) for name in (
        "cmd_repoint", "cmd_restore_editable", "cmd_verify_editable",
        "_record_release_intent_lab", "_record_release_intent_remote",
        "_verify_remote_editable", "_remote_editable_intent",
    )]
    sources.append(release_quad.REMOTE_EDITABLE_VERIFY)
    for text in sources:
        for forbidden in ("pip uninstall", "uninstall -y", "Stop-Process", "taskkill",
                          "_kill_mcp_local", "_kill_cubit_local"):
            assert forbidden not in text


def test_cli_registers_repoint_and_keeps_restore_as_a_refusal():
    text = _TOOL.read_text(encoding="utf-8")
    assert '"repoint":          cmd_repoint,' in text
    assert 'sub.add_parser(\n        "repoint",' in text
    assert "removed: no default development source" in text
