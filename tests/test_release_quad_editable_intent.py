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
    assert release_quad._PATH_ALIASES[-1][0] in release_quad.REMOTE_EDITABLE_VERIFY








def test_verify_editable_checks_lab_wheel_and_100_editable(monkeypatch):
    calls = []
    monkeypatch.setattr(release_quad, "_verify_lab_wheel", lambda: calls.append("lab-wheel") or 0)
    monkeypatch.setattr(release_quad, "_verify_100_editable", lambda: calls.append("100-editable") or 0)
    assert release_quad.cmd_verify_editable(None) == 0
    assert calls == ["lab-wheel", "100-editable"]

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


@pytest.mark.parametrize("explicit", [False, True])
def test_solver_remote_verify_does_not_expand_to_other_recorded_products(monkeypatch, explicit):
    """Exercise the remote CLI parser: --expect alone does not limit packages."""
    source = "W:/release/radia"
    monkeypatch.setenv(release_quad.EDITABLE_REPO_100_ENV, source)
    selected = []
    result = _remote_ok({"match": 1, "drift": 0, "unverified": 0})[1]

    def verify(packages, expectations, **kwargs):
        selected.extend(packages)
        assert expectations == {"radia": source}
        return result

    def remote(host, argv):
        assert intent.main(argv) == 0
        return 0, result, ""

    monkeypatch.setattr(intent, "verify", verify)
    monkeypatch.setattr(release_quad, "_remote_editable_intent", remote)
    assert release_quad._verify_100_editable([("radia", source)] if explicit else None) == 0
    assert selected == ["radia"]


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


@pytest.mark.parametrize("record", [False, True])
@pytest.mark.parametrize("source_rc", [0, 4])
def test_done_default_source_is_controller_and_ignores_legacy_lab_intent(
        isolated_record, monkeypatch, record, source_rc):
    # Legacy LAB editable records no longer select or gate the runtime source.
    if record:
        _record(isolated_record, "radia", "S:/Radia/old-lab-editable/")
    assert release_quad._done_release_source(SimpleNamespace()) == str(release_quad.REPO)
    monkeypatch.setattr(release_quad, "cmd_preflight", lambda _args: 0)
    monkeypatch.setattr(release_quad, "_release_head", lambda: "a" * 40)
    seen = []
    monkeypatch.setattr(release_quad, "_verify_local_release_source",
                        lambda repo, sha: seen.append((repo, sha)) or source_rc)
    stop = AssertionError("stop after the source gate")
    def later_gate(_args):
        raise stop
    monkeypatch.setattr(release_quad, "cmd_temp_shadows", later_gate)
    args = SimpleNamespace(simulink_package=None, release_source=None)
    if source_rc:
        assert release_quad.cmd_done(args) == source_rc
    else:
        with pytest.raises(AssertionError) as raised:
            release_quad.cmd_done(args)
        assert raised.value is stop
    assert seen == [(str(release_quad.REPO), "a" * 40)]

@pytest.mark.parametrize("rc", [0, 2, 3])
def test_lab_deploy_uses_wheel_route_and_propagates_failure(monkeypatch, rc):
    calls = []
    monkeypatch.setattr(release_quad, "_deploy_pypi",
                        lambda host, label, **kw: calls.append((host, label, kw)) or rc)
    assert release_quad._deploy_lab() == rc
    assert calls == [("102", "LAB", {})]



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


def test_repoint_on_lab_is_refused_before_action(monkeypatch, capsys):
    actions = []
    monkeypatch.setattr(intent, "main", lambda argv: actions.append(argv) or 0)
    monkeypatch.setattr(release_quad, "_remote_editable_intent",
                        lambda *args, **kwargs: actions.append(args))
    args = SimpleNamespace(host="lab", package=["radia"], source=["S:/old"],
                           reason="test", via=None, record_current=False,
                           rollback=False, dry_run=False, require_pushed=False)
    assert release_quad.cmd_repoint(args) == 2
    assert actions == []
    assert "only on 100号機" in capsys.readouterr().out

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
        "_record_release_intent_remote",
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
