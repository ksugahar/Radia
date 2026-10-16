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
