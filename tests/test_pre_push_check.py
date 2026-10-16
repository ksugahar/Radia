"""Outgoing ref coverage and fail-closed dispatch, without SSH or live hooks."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

spec = importlib.util.spec_from_file_location(
    "pre_push_check", Path(__file__).parents[1] / "tools/pre_push_check.py")
check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check)
A, B, ZERO = "a" * 40, "b" * 40, "0" * 40
REAL_LOCAL = check.run_local_candidate


@pytest.fixture(autouse=True)
def isolate_local_checks(monkeypatch):
    monkeypatch.setattr(check, "run_local_candidate", lambda base, head: 0)


def record(ref, head=B, base=A):
    return f"refs/heads/local {head} {ref} {base}\n"


def test_all_branches_and_new_branch_exact_shas():
    assert check.candidates([record("refs/heads/main"),
                             record("refs/heads/codex/new", base=ZERO)]) == [(A, B), (ZERO, B)]


def test_deletes_tags_and_duplicate_pairs():
    assert check.candidates([record("refs/heads/old", head=ZERO),
                             record("refs/tags/v1"), record("refs/heads/one"),
                             record("refs/heads/two")]) == [(A, B)]


@pytest.mark.parametrize("bad", ["malformed", record("refs/heads/main", head="HEAD")])
def test_invalid_later_record_prevents_all_remote_work(monkeypatch, bad):
    monkeypatch.setattr(check.subprocess, "run", lambda *a, **kw: pytest.fail("remote started"))
    assert check.main([record("refs/heads/valid"), bad]) == 1


def test_failed_gate_stops_multi_ref_push(monkeypatch):
    calls = []
    def run(command, **kwargs):
        calls.append(command)
        return SimpleNamespace(returncode=7)
    monkeypatch.setattr(check.subprocess, "run", run)
    assert check.main([record("refs/heads/one"), record("refs/heads/two", base=ZERO)]) == 7
    assert len(calls) == 1
    assert calls[0][-4:] == ["--base", A, "--head", B]


def test_transport_failure_is_not_success(monkeypatch):
    def run(*args, **kwargs):
        raise OSError("unavailable")
    monkeypatch.setattr(check.subprocess, "run", run)
    assert check.main([record("refs/heads/new", base=ZERO)]) == 1


def test_success_dispatches_every_distinct_candidate(monkeypatch):
    calls = []
    monkeypatch.setattr(check.subprocess, "run", lambda cmd, **kw:
                        calls.append(cmd) or SimpleNamespace(returncode=0))
    assert check.main([record("refs/heads/main"), record("refs/heads/new", base=ZERO)]) == 0
    assert [cmd[-4:] for cmd in calls] == [["--base", A, "--head", B], ["--base", ZERO, "--head", B]]


def test_missing_helper_blocks_push(monkeypatch, tmp_path):
    monkeypatch.setattr(check, "ROOT", tmp_path)
    assert check.main([record("refs/heads/new")]) == 1


def test_local_failure_blocks_remote_transport(monkeypatch):
    monkeypatch.setattr(check, "run_local_candidate", lambda base, head: 5)
    monkeypatch.setattr(check.subprocess, "run", lambda *a, **kw: pytest.fail("SSH launched"))
    assert check.main([record("refs/heads/new")]) == 5


def test_local_candidate_precedes_each_remote_preflight(monkeypatch):
    events = []
    monkeypatch.setattr(check, "run_local_candidate", lambda base, head:
                        events.append(("local", base, head)) or 0)
    monkeypatch.setattr(check.subprocess, "run", lambda cmd, **kw:
                        events.append(("remote", cmd[-3], cmd[-1])) or SimpleNamespace(returncode=0))
    assert check.main([record("refs/heads/main"), record("refs/heads/new", base=ZERO)]) == 0
    assert events == [("local", A, B), ("remote", A, B),
                      ("local", ZERO, B), ("remote", ZERO, B)]


@pytest.mark.parametrize("code", [0, 7])
def test_local_checker_uses_detached_candidate_and_cleans_only_its_tree(monkeypatch, tmp_path, code):
    parent = tmp_path / "owned"
    parent.mkdir()
    monkeypatch.setattr(check.tempfile, "mkdtemp", lambda **kw: str(parent))
    calls = []
    def run(command, **kwargs):
        calls.append((command, kwargs))
        if "worktree" in command and "add" in command:
            (parent / "candidate").mkdir()
        elif "--local-push-candidate" in command:
            assert kwargs["cwd"] == parent / "candidate"
            assert command[-4:] == ["--local-push-candidate", B, "--local-push-base", A]
            return SimpleNamespace(returncode=code)
        elif "worktree" in command and "remove" in command:
            (parent / "candidate").rmdir()
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(check.subprocess, "run", run)
    assert REAL_LOCAL(A, B) == code
    assert not parent.exists()
    assert any("--detach" in command for command, _ in calls)


def test_local_cleanup_failure_blocks_and_retains_owned_path(monkeypatch, tmp_path):
    parent = tmp_path / "owned"
    parent.mkdir()
    monkeypatch.setattr(check.tempfile, "mkdtemp", lambda **kw: str(parent))
    def run(command, **kwargs):
        if "add" in command:
            (parent / "candidate").mkdir()
        return SimpleNamespace(returncode=2 if "remove" in command else 0)
    monkeypatch.setattr(check.subprocess, "run", run)
    assert REAL_LOCAL(A, B) == 1
    assert (parent / "candidate").is_dir()

def test_cleanup_exception_preserves_original_failure(monkeypatch, tmp_path, capsys):
    parent = tmp_path/'owned'
    parent.mkdir()
    monkeypatch.setattr(check.tempfile, 'mkdtemp', lambda **kw: str(parent))
    def run(command, **kw):
        if 'add' in command:
            (parent/'candidate').mkdir()
        if '--local-push-candidate' in command:
            raise OSError('original check failure')
        if 'remove' in command:
            raise ValueError('cleanup refusal')
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(check.subprocess, 'run', run)
    assert REAL_LOCAL(A, B) == 1
    error = capsys.readouterr().err
    assert 'original check failure' in error and 'cleanup refusal' in error
    assert (parent/'candidate').exists()


def test_missing_base_suggests_fetch_without_checkout(monkeypatch, capsys):
    import subprocess
    def run(command, **kw):
        assert 'cat-file' in command
        raise subprocess.CalledProcessError(128, command)
    monkeypatch.setattr(check.subprocess, 'run', run)
    assert REAL_LOCAL(A, B) == 1
    assert 'fetch the remote base commit/history' in capsys.readouterr().err
