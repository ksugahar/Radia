"""Recorded editable intent: record, verify, repoint, rollback (tools/editable_intent.py)."""
from __future__ import annotations

import base64
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

_TOOL = Path(__file__).resolve().parents[1] / "tools" / "editable_intent.py"
_SPEC = importlib.util.spec_from_file_location("radia_editable_intent_under_test", _TOOL)
assert _SPEC is not None and _SPEC.loader is not None
ei = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(ei)


@pytest.fixture
def intent_file(tmp_path, monkeypatch):
    path = tmp_path / "state" / "editable-intent.json"
    monkeypatch.setenv(ei.INTENT_FILE_ENV, str(path))
    return path


@pytest.fixture
def sources(tmp_path):
    old = tmp_path / "old-tree" / "packages" / "radia-mcp"
    new = tmp_path / "new-tree" / "packages" / "radia-mcp"
    for tree in (old, new):
        tree.mkdir(parents=True)
        (tree / "pyproject.toml").write_text("[project]\nname='radia-mcp'\n", encoding="utf-8")
    return str(old), str(new)


class FakeInterpreter:
    """Installed registration and fresh-process import as seen by the tool."""

    def __init__(self, source, installed=True, editable=True):
        self.source = source
        self.installed = installed
        self.editable = editable
        self.pip_calls = []

    def observe(self, package):
        if not self.installed:
            return {"installed": False, "version": None, "editable": False, "source": None, "pth": []}
        return {"installed": True, "version": "1.4.53", "editable": self.editable,
                "source": self.source if self.editable else None, "pth": []}

    def origin(self, package, executable=None):
        if not self.installed or not self.editable:
            return None
        return str(Path(self.source) / "src" / "radia_mcp" / "__init__.py")

    def runner(self, returncode=0, move_to=None):
        def run(command, **kwargs):
            self.pip_calls.append(command)
            if returncode == 0 and move_to is not None:
                self.source = move_to
            return subprocess.CompletedProcess(command, returncode, stdout="pip output", stderr="")
        return run


@pytest.fixture
def fake(monkeypatch, sources):
    interpreter = FakeInterpreter(sources[0])
    monkeypatch.setattr(ei, "observe_registration", interpreter.observe)
    monkeypatch.setattr(ei, "fresh_import_origin", interpreter.origin)
    return interpreter


def _log_events(intent_file):
    log = ei.log_file_path(intent_file)
    if not log.exists():
        return []
    return [json.loads(line)["event"] for line in log.read_text(encoding="utf-8").splitlines()]


def test_record_is_keyed_by_interpreter_and_package(intent_file):
    data = ei.load_intent(intent_file)
    ei.set_entry(data, "radia-mcp", {"source": "S:/Radia/release-quad/x/packages/radia-mcp"})
    ei.save_intent(data, intent_file)

    loaded = ei.load_intent(intent_file)
    assert loaded["schema"] == ei.SCHEMA
    assert ei.recorded_entry(loaded, "radia-mcp")["source"].endswith("radia-mcp")
    assert ei.recorded_entry(loaded, "radia") is None
    assert ei.recorded_entry(loaded, "radia-mcp", executable="C:/other/python.exe") is None
    assert ei.recorded_packages(loaded) == ["radia-mcp"]


def test_norm_path_unifies_nas_spellings():
    unc = r"\\192.168.121.100\work\00_CAE\Radia\release-quad\x\packages\radia-mcp"
    assert ei.norm_path(unc) == "s:/radia/release-quad/x/packages/radia-mcp"
    assert ei.norm_path(r"W:\00_CAE\Radia\release-quad\x\\") == "s:/radia/release-quad/x"
    assert ei.norm_path("S:/Radia/release-quad/X/") == ei.norm_path(unc[: -len(r"\packages\radia-mcp")])


def test_verify_without_record_is_unverified_and_names_no_target(intent_file, fake, sources):
    report = ei.verify(["radia-mcp"])
    row = report["rows"][0]
    assert row["status"] == "unverified"
    assert row["expectation_basis"] == "none"
    assert row["registration"] == sources[0]
    assert report["exit_code"] == ei.EXIT_UNVERIFIED == 5

    text = ei.format_verify(report)
    assert "UNVERIFIED" in text
    assert "record-current" in text
    for forbidden in ("01_GitHub", "pip uninstall", "Stop-Process", "canonical"):
        assert forbidden not in text


def test_verify_compares_against_the_record_or_an_explicit_expectation(intent_file, fake, sources):
    old, new = sources
    data = ei.load_intent(intent_file)
    ei.set_entry(data, "radia-mcp", {"source": old})
    ei.save_intent(data, intent_file)

    assert ei.verify(["radia-mcp"])["rows"][0]["status"] == "match"

    data = ei.load_intent(intent_file)
    ei.set_entry(data, "radia-mcp", {"source": new})
    ei.save_intent(data, intent_file)
    drifted = ei.verify(["radia-mcp"])
    assert drifted["rows"][0]["status"] == "drift"
    assert drifted["exit_code"] == ei.EXIT_DRIFT
    assert "No default target" in ei.format_verify(drifted)

    explicit = ei.verify(["radia-mcp"], {"radia-mcp": old})
    assert explicit["rows"][0]["status"] == "match"
    assert explicit["rows"][0]["expectation_basis"] == "explicit"


def test_verify_reports_not_installed_and_not_editable(intent_file, monkeypatch):
    absent = FakeInterpreter("x", installed=False)
    monkeypatch.setattr(ei, "observe_registration", absent.observe)
    monkeypatch.setattr(ei, "fresh_import_origin", absent.origin)
    assert ei.verify(["radia"])["rows"][0]["status"] == "not_installed"

    wheel = FakeInterpreter("x", editable=False)
    monkeypatch.setattr(ei, "observe_registration", wheel.observe)
    monkeypatch.setattr(ei, "fresh_import_origin", wheel.origin)
    assert ei.verify(["radia"])["rows"][0]["status"] == "not_editable"


def test_repoint_records_previous_before_pip_and_intent_after_verification(intent_file, fake, sources):
    old, new = sources
    report = ei.repoint([("radia-mcp", new)], "try the reviewed tree",
                        runner=fake.runner(move_to=new))

    assert report["exit_code"] == 0
    assert report["results"][0]["status"] == "repointed"
    command = fake.pip_calls[0]
    assert command == [sys.executable, "-m", "pip", "install", "--no-deps", "--no-cache-dir",
                       "--no-build-isolation", "-e", new]
    assert "uninstall" not in " ".join(command)

    entry = ei.recorded_entry(ei.load_intent(intent_file), "radia-mcp")
    assert entry["source"] == new
    assert entry["reason"] == "try the reviewed tree"
    assert entry["recorded_via"] == "repoint"
    assert entry["previous"]["source"] == old
    assert _log_events(intent_file) == ["repoint-pre", "repoint-post"]
    assert "reconnected" in report["note"]


def test_repoint_blocked_pip_leaves_the_record_untouched(intent_file, fake, sources):
    old, new = sources
    report = ei.repoint([("radia-mcp", new)], "try", runner=fake.runner(returncode=1))

    assert report["exit_code"] == ei.EXIT_ACTION == 3
    result = report["results"][0]
    assert result["status"] == "blocked"
    assert "rollback" in result["hint"]
    assert "uninstall" not in result["hint"]
    assert ei.recorded_entry(ei.load_intent(intent_file), "radia-mcp") is None
    assert _log_events(intent_file) == ["repoint-pre", "repoint-blocked"]
    assert fake.source == old


def test_repoint_verification_failure_leaves_the_record_untouched(intent_file, fake, sources):
    _old, new = sources
    report = ei.repoint([("radia-mcp", new)], "try", runner=fake.runner(returncode=0))

    assert report["exit_code"] == ei.EXIT_VERIFY == 4
    assert report["results"][0]["status"] == "verify_failed"
    assert ei.recorded_entry(ei.load_intent(intent_file), "radia-mcp") is None
    assert _log_events(intent_file) == ["repoint-pre", "repoint-verify-failed"]


def test_repoint_stops_at_the_first_failed_package(intent_file, fake, sources, tmp_path):
    _old, new = sources
    missing = str(tmp_path / "does-not-exist")
    report = ei.repoint([("radia-mcp", missing), ("radia", new)], "try",
                        runner=fake.runner(move_to=new))
    assert report["exit_code"] == ei.EXIT_PRECONDITION
    assert [r["status"] for r in report["results"]] == ["refused"]
    assert fake.pip_calls == []


def test_rollback_reinstalls_the_previous_recorded_pointer(intent_file, fake, sources):
    old, new = sources
    assert ei.repoint([("radia-mcp", new)], "try", runner=fake.runner(move_to=new))["exit_code"] == 0

    report = ei.rollback(["radia-mcp"], "undo", runner=fake.runner(move_to=old))
    assert report["exit_code"] == 0
    assert fake.pip_calls[-1][-1] == old
    entry = ei.recorded_entry(ei.load_intent(intent_file), "radia-mcp")
    assert entry["source"] == old
    assert entry["recorded_via"] == "rollback"
    assert entry["previous"]["source"] == new


def test_rollback_without_a_previous_pointer_is_refused(intent_file, fake):
    report = ei.rollback(["radia-mcp"], "undo", runner=fake.runner())
    assert report["exit_code"] == ei.EXIT_PRECONDITION
    assert report["results"][0]["status"] == "refused"
    assert fake.pip_calls == []


def test_record_current_adopts_only_an_editable_installation(intent_file, fake, sources, monkeypatch):
    report = ei.record_current(["radia-mcp"], "adopt the pointer codex verified")
    assert report["exit_code"] == 0
    entry = ei.recorded_entry(ei.load_intent(intent_file), "radia-mcp")
    assert entry["source"] == sources[0]
    assert entry["recorded_via"] == "record-current"
    assert _log_events(intent_file) == ["record-current"]

    fake.editable = False
    refused = ei.record_current(["radia-mcp"], "adopt")
    assert refused["exit_code"] == ei.EXIT_PRECONDITION
    assert refused["results"][0]["status"] == "refused"


def test_require_pushed_refuses_a_commit_no_remote_ref_contains(intent_file, fake, sources, monkeypatch):
    _old, new = sources
    monkeypatch.setattr(ei, "describe_source", lambda source: {
        "source": str(source), "exists": True, "project": True,
        "commit": "f" * 40, "tracked_clean": True})
    monkeypatch.setattr(ei, "pushed_refs", lambda source, commit: [])

    report = ei.repoint([("radia-mcp", new)], "handoff", require_pushed=True,
                        runner=fake.runner(move_to=new))
    assert report["exit_code"] == ei.EXIT_PRECONDITION
    assert "--require-pushed" in report["results"][0]["detail"]
    assert fake.pip_calls == []

    monkeypatch.setattr(ei, "pushed_refs", lambda source, commit: ["origin/main"])
    report = ei.repoint([("radia-mcp", new)], "handoff", require_pushed=True,
                        runner=fake.runner(move_to=new))
    assert report["exit_code"] == 0
    assert report["results"][0]["entry"]["pushed_refs"] == ["origin/main"]


def test_dry_run_installs_and_records_nothing(intent_file, fake, sources):
    _old, new = sources

    def must_not_run(command, **kwargs):
        raise AssertionError("dry run must not call pip")

    report = ei.repoint([("radia-mcp", new)], "plan", dry_run=True, runner=must_not_run)
    assert report["exit_code"] == 0
    assert report["results"][0]["status"] == "planned"
    assert report["results"][0]["command"][-2:] == ["-e", new]
    assert not intent_file.exists()
    assert "dry run" in ei.format_action(report)


def test_describe_source_reports_missing_and_non_project_trees(tmp_path):
    missing = ei.describe_source(tmp_path / "nowhere")
    assert missing["exists"] is False and missing["commit"] is None
    plain = tmp_path / "plain"
    plain.mkdir()
    described = ei.describe_source(plain)
    assert described["exists"] is True and described["project"] is False
    assert described["repository"] is None and described["commit"] is None


def test_describe_source_finds_the_commit_from_a_package_subdirectory(tmp_path):
    """packages/radia-mcp lives below the repository root; git runs at that root."""
    repo = tmp_path / "tree"
    package = repo / "packages" / "radia-mcp"
    package.mkdir(parents=True)
    (package / "pyproject.toml").write_text("[project]\nname='radia-mcp'\n", encoding="utf-8")

    def git(*args):
        return subprocess.run(["git", "-C", str(repo), *args], check=True,
                              capture_output=True, text=True).stdout.strip()

    git("init")
    git("add", "packages/radia-mcp/pyproject.toml")
    git("-c", "user.name=Radia Test", "-c", "user.email=radia-test@example.invalid",
        "commit", "-m", "package")
    head = git("rev-parse", "HEAD").lower()

    described = ei.describe_source(package)
    assert described["repository"] == str(repo)
    assert described["commit"] == head
    assert described["tracked_clean"] is True
    assert ei.pushed_refs(package, head) == []

    (package / "pyproject.toml").write_text("[project]\nname='radia-mcp'\n# wip\n", encoding="utf-8")
    assert ei.describe_source(package)["tracked_clean"] is False


def test_cli_accepts_a_base64_argv_token_and_prints_json(intent_file, capsys):
    token = base64.b64encode(json.dumps(["--json", "show"]).encode("utf-8")).decode("ascii")
    assert ei.main(["--argv-b64", token]) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["schema"] == ei.SCHEMA
    assert report["packages"] == {}


def test_cli_repoint_requires_pairs_and_a_reason(intent_file):
    with pytest.raises(SystemExit):
        ei.main(["repoint", "--package", "radia-mcp"])
    with pytest.raises(SystemExit):
        ei.main(["repoint", "--package", "radia-mcp", "--source", "x"])
    with pytest.raises(SystemExit):
        ei.main(["repoint", "--record-current", "--source", "x", "--reason", "r"])


def test_tool_source_never_uninstalls_or_stops_processes():
    source = _TOOL.read_text(encoding="utf-8")
    for forbidden in ("pip uninstall", "uninstall -y", "Stop-Process", "taskkill", "os.kill"):
        assert forbidden not in source
