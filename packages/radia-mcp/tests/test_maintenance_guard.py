"""No real installation, configuration or live client is touched by these tests."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
from radia_mcp import _editable_activation as a
from radia_mcp import _maintenance_guard as g
from radia_mcp import maintenance as m

AUTH = dict(owner="test-owner", reason="approved test", clients_idle=True)


@pytest.fixture(autouse=True)
def isolated(monkeypatch, tmp_path):
    monkeypatch.setattr(g, "STATE_ROOT", tmp_path / "guard")


def operation(**kwargs):
    return g.change(
        **{**AUTH, **kwargs},
        scope={"operation": "test"},
        expected={"source": "old"},
        proposed={"source": "new"},
    )


@pytest.mark.parametrize("kwargs", [{"owner": ""}, {"reason": ""}, {"clients_idle": False}])
def test_unknown_owner_or_busy_clients_never_begin(kwargs):
    with pytest.raises(ValueError):
        with operation(**kwargs):
            pytest.fail("must not mutate")
    assert not g.STATE_ROOT.exists()


def test_completed_receipts_remain_and_never_claim_live_verification():
    with operation() as first:
        first["observed"] = {"source": "new"}
    with operation() as second:
        pass
    receipts = list((g.STATE_ROOT / "receipts").glob("*.json"))
    assert len(receipts) == 2
    assert second["previous_change_id"] == first["change_id"]
    assert second["status"] == "completed"
    assert second["live_status"] == "unverified"
    assert second["client_reconnect"] == "not-performed"


def test_failed_change_requires_explicit_reconciliation_and_redacts_errors():
    with pytest.raises(RuntimeError):
        with operation() as failed:
            raise RuntimeError("secret-token-do-not-log")
    raw = (g.STATE_ROOT / "state.json").read_text()
    assert "secret-token" not in raw
    assert json.loads(raw)["status"] == "failed"
    with pytest.raises(ValueError, match="Unresolved"):
        with operation():
            pytest.fail("must not retry")
    with pytest.raises(ValueError, match="matching"):
        g.reconcile("bad", **AUTH, observe=lambda _: {})
    result = g.reconcile(failed["change_id"], **AUTH, observe=lambda _: {"checked": True})
    assert result["status"] == "reconciled"
    assert result["observed"]["checked"]
    old = json.loads((g.STATE_ROOT / "receipts" / (failed["change_id"] + ".json")).read_text())
    assert old["status"] == "failed"
    with operation():
        pass


def _child(script):
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(Path(g.__file__).parents[1])
    return subprocess.run(
        [sys.executable, "-c", script, str(g.STATE_ROOT)],
        env=environment,
        capture_output=True,
        text=True,
        timeout=20,
    )


def test_another_process_cannot_acquire_or_reconcile_an_active_owner():
    with operation() as record:
        result = _child(
            "import sys;from pathlib import Path;from radia_mcp import _maintenance_guard as g;"
            "g.STATE_ROOT=Path(sys.argv[1]);\n"
            "try:\n"
            f" g.reconcile('{record['change_id']}',owner='other',reason='test',clients_idle=True,observe=lambda _:{{}})\n"
            "except ValueError as e:\n print(str(e));sys.exit(7)\n"
        )
        assert result.returncode == 7, result.stderr
        assert "deployment lock" in result.stdout


def test_process_death_releases_os_lock_but_not_unresolved_change():
    result = _child(
        "import os,sys;from pathlib import Path;from radia_mcp import _maintenance_guard as g;"
        "g.STATE_ROOT=Path(sys.argv[1]);"
        "ctx=g.change(owner='dead',reason='test',clients_idle=True,scope={},expected={},proposed={});"
        "ctx.__enter__();os._exit(42)"
    )
    assert result.returncode == 42, result.stderr
    with pytest.raises(ValueError, match="Unresolved"):
        with operation():
            pytest.fail("must not steal ownership")
    assert json.loads((g.STATE_ROOT / "state.json").read_text())["status"] == "running"


def test_atomic_publish_failure_preserves_previous_complete_json(monkeypatch, tmp_path):
    target = tmp_path / "state.json"
    g._atomic(target, {"before": 1})
    monkeypatch.setattr(g.os, "replace", lambda *_: (_ for _ in ()).throw(OSError("test")))
    with pytest.raises(OSError):
        g._atomic(target, {"after": 2})
    assert json.loads(target.read_text()) == {"before": 1}
    assert not list(tmp_path.glob("*.tmp"))


def test_config_apply_requires_scope_and_records_only_digests(tmp_path):
    path = tmp_path / "config.json"
    before = b'{"secret":"private"}'
    after = b'{"secret":"retained"}'
    path.write_bytes(before)
    with pytest.raises(ValueError):
        m.apply_config(path, before, after)
    assert path.read_bytes() == before
    m.apply_config(path, before, after, **AUTH, targets="test client")
    raw = (g.STATE_ROOT / "state.json").read_text()
    assert "private" not in raw and "retained" not in raw
    assert json.loads(raw)["observed"]["sha256"] == m._digest(after)


def test_config_compare_happens_inside_guard_and_preserves_concurrent_edit(tmp_path):
    path = tmp_path / "config.json"
    path.write_bytes(b"other owner's data")
    with pytest.raises(ValueError, match="changed"):
        m.apply_config(path, b"old", b"new", **AUTH, targets="test client")
    assert path.read_bytes() == b"other owner's data"
    assert json.loads((g.STATE_ROOT / "state.json").read_text())["status"] == "failed"


@pytest.fixture
def activation(monkeypatch, tmp_path):
    old, new = tmp_path / "old", tmp_path / "new"
    old.mkdir()
    new.mkdir()
    snapshots = {
        old: {"source": g.identity(old), "commit": "a" * 40, "version": "1.4.53"},
        new: {"source": g.identity(new), "commit": "b" * 40, "version": "1.4.53"},
    }
    monkeypatch.setattr(a, "_snapshot", lambda path, commit, clean: snapshots[path])
    state = {**snapshots[old], "version": "1.4.53"}
    monkeypatch.setattr(a, "installed", lambda: dict(state))
    calls = []

    def pip(path):
        calls.append(path)
        state.update(snapshots[path])

    monkeypatch.setattr(a, "_pip", pip)
    return old, new, state, calls


def test_same_version_different_source_is_activated_and_not_called_live(activation):
    old, new, state, calls = activation
    result = a.activate_editable(new, "b" * 40, old, "a" * 40, **AUTH, targets="test clients")
    assert calls == [new]
    assert result["observed"]["version"] == "1.4.53"
    assert result["live_status"] == "unverified"
    assert result["observed"]["source"] == g.identity(new)


def test_unexpected_registration_blocks_install(activation):
    old, new, state, calls = activation
    state["source"] = "a third session changed this"
    with pytest.raises(ValueError, match="since approval"):
        a.activate_editable(new, "b" * 40, old, "a" * 40, **AUTH, targets="clients")
    assert calls == []


def test_install_failure_never_rolls_back_or_retries(activation, monkeypatch):
    old, new, _, calls = activation

    def failing(path):
        calls.append(path)
        raise subprocess.CalledProcessError(1, "pip", stderr=b"private-output")

    monkeypatch.setattr(a, "_pip", failing)
    with pytest.raises(subprocess.CalledProcessError):
        a.activate_editable(new, "b" * 40, old, "a" * 40, **AUTH, targets="clients")
    assert calls == [new]
    assert "private-output" not in (g.STATE_ROOT / "state.json").read_text()


def test_alias_normalization_uses_resolved_paths(tmp_path):
    folder = tmp_path / "source"
    folder.mkdir()
    assert g.identity(folder) == g.identity(folder / ".." / "source")


def test_pip_is_scoped_and_ignores_task_overrides(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setenv("PYTHONPATH", "unrelated-task")
    monkeypatch.setenv("PYTHONHOME", "unrelated-runtime")
    monkeypatch.setattr(a.subprocess, "run", lambda args, **kw: calls.append((args, kw)))
    a._pip(tmp_path)
    args, kw = calls[0]
    assert args[:5] == [sys.executable, "-I", "-m", "pip", "--isolated"]
    assert args[-4:] == ["--prefix", sys.prefix, "-e", str(tmp_path)]
    assert "--no-deps" in args and "--no-build-isolation" in args
    assert "PYTHONPATH" not in kw["env"] and "PYTHONHOME" not in kw["env"]
    assert kw["check"] and kw["timeout"] == 300


def test_registration_and_fresh_resolution_mismatch_rejected(monkeypatch, tmp_path):
    source = tmp_path / "source"
    module = source / "src/radia_mcp/__init__.py"
    module.parent.mkdir(parents=True)
    module.touch()
    other = tmp_path / "other.py"
    other.touch()
    payload = {
        "version": "same",
        "direct": {"url": source.as_uri(), "dir_info": {"editable": True}},
        "module": str(other),
    }
    monkeypatch.setattr(
        a.subprocess,
        "run",
        lambda *args, **kw: subprocess.CompletedProcess(args, 0, stdout=json.dumps(payload)),
    )
    with pytest.raises(ValueError, match="differs from editable"):
        a.installed()


def test_cli_config_mutation_denied_without_authority(tmp_path, capsys):
    path = tmp_path / "client.json"
    path.write_text("{}")
    assert m.main(["config", str(path), "--apply"]) == 2
    assert path.read_text() == "{}"
    assert "Maintenance failed" in capsys.readouterr().err


def test_activation_receipt_records_each_pending_client(activation):
    old, new, _, _ = activation
    result = a.activate_editable(
        new,
        "b" * 40,
        old,
        "a" * 40,
        **AUTH,
        targets="host/user/codex/server,host/user/claude/server",
    )
    receipt = json.loads(Path(result["receipt"]).read_text())
    assert len(receipt["client_results"]) == 2
    assert all(
        c["status"] == "unverified" and c["observed_at"] is None for c in receipt["client_results"]
    )


def test_real_git_snapshot_checks_commit_and_tracked_changes(tmp_path):
    repo = tmp_path / "repo"
    package = repo / "packages/radia-mcp"
    code = package / "src/radia_mcp/__init__.py"
    code.parent.mkdir(parents=True)
    code.write_text('__version__ = "test"')
    (package / "pyproject.toml").write_text('[project]\nname="radia-mcp"\nversion="test"\n')

    def git(*args):
        return subprocess.run(
            ["git", "-C", str(repo), *args], check=True, capture_output=True, text=True
        ).stdout.strip()

    git("init")
    git("add", ".")
    git("-c", "user.name=Test", "-c", "user.email=test@example.invalid", "commit", "-m", "fixture")
    sha = git("rev-parse", "HEAD")
    assert a._snapshot(package, sha, clean=True)["commit"] == sha
    with pytest.raises(ValueError, match="full Git SHA"):
        a._snapshot(package, sha[:9], clean=True)
    with pytest.raises(ValueError, match="differs"):
        a._snapshot(package, "0" * 40, clean=True)
    code.write_text("changed")
    with pytest.raises(ValueError, match="uncommitted changes"):
        a._snapshot(package, sha, clean=True)
