"""Outgoing ref coverage and fail-closed dispatch, without SSH or live hooks."""
import importlib.util
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

spec = importlib.util.spec_from_file_location(
    "pre_push_check", Path(__file__).parents[1] / "tools/pre_push_check.py")
check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check)
A, B, ZERO = "a" * 40, "b" * 40, "0" * 40
REAL_LOCAL = check.run_local_candidate
REAL_EFFECTIVE_BASE = check.effective_base


@pytest.fixture(autouse=True)
def isolate_local_checks(monkeypatch):
    monkeypatch.setattr(check, "run_local_candidate", lambda base, head: 0)
    monkeypatch.setattr(check, "effective_base", lambda base, head: base)


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


def _git(cwd, *args):
    return subprocess.run(["git", "-c", f"safe.directory={cwd}", *args], cwd=cwd, check=True,
                          capture_output=True, text=True, encoding="utf-8").stdout.strip()


STUB = """import json, os, subprocess, sys, time
head = subprocess.run(['git', '-c', 'safe.directory=' + os.getcwd(), 'rev-parse', 'HEAD'],
                      capture_output=True, text=True).stdout.strip()
marker = open('marker.txt').read().strip()
json.dump({'argv': sys.argv[1:], 'head': head, 'marker': marker, 'cwd': os.getcwd(),
           'git_env': sorted(k for k in os.environ if k.startswith('GIT_'))},
          open(os.environ['STUB_RECORD'], 'w'))
time.sleep(float(os.environ.get('STUB_SLEEP', '0')))
raise SystemExit(int(os.environ.get('STUB_CODE', '0')))
"""


@pytest.fixture
def source_repo(tmp_path, monkeypatch):
    """A throwaway clone: base and candidate commits, checkout left dirty at base."""
    repo = tmp_path / "source"
    repo.mkdir()
    _git(repo, "init", "-q")
    _git(repo, "config", "user.email", "t@example.invalid")
    _git(repo, "config", "user.name", "t")
    (repo / "tools").mkdir()
    (repo / "tools" / "ci_preflight.py").write_text(STUB, encoding="utf-8")
    (repo / "marker.txt").write_text("base\n", encoding="utf-8")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "base")
    base = _git(repo, "rev-parse", "HEAD")
    (repo / "marker.txt").write_text("candidate\n", encoding="utf-8")
    _git(repo, "commit", "-q", "-am", "candidate")
    head = _git(repo, "rev-parse", "HEAD")
    _git(repo, "checkout", "-q", base)  # a different checkout ...
    (repo / "marker.txt").write_text("dirty WIP\n", encoding="utf-8")  # ... with WIP
    monkeypatch.setattr(check, "ROOT", repo)
    owned = tmp_path / "owned"
    owned.mkdir()
    monkeypatch.setattr(check.tempfile, "mkdtemp", lambda **kw: str(owned))
    record = tmp_path / "record.json"
    monkeypatch.setenv("STUB_RECORD", str(record))
    monkeypatch.setenv("GIT_DIR", str(tmp_path / "leaked-hook-git-dir"))  # hook-local state
    return SimpleNamespace(repo=repo, base=base, head=head, owned=owned, record=record)


def _real_git_calls(monkeypatch):
    calls = []
    real = check.subprocess.run

    def spy(command, **kw):
        calls.append(list(command))
        return real(command, **kw)

    monkeypatch.setattr(check.subprocess, "run", spy)
    return calls


@pytest.mark.parametrize("code", [0, 7])
def test_exact_candidate_export_from_dirty_other_checkout(monkeypatch, source_repo, code):
    """Guards: the exact outgoing SHA is checked (not checkout HEAD or WIP), the
    child gets no GIT_* state, its exit code propagates, no worktree is used
    and the owned export is removed."""
    import json
    monkeypatch.setenv("STUB_CODE", str(code))
    monkeypatch.delenv("GIT_DIR", raising=False)
    calls = _real_git_calls(monkeypatch)
    monkeypatch.setenv("GIT_DIR", str(source_repo.owned.parent / "leaked-hook-git-dir"))
    assert REAL_LOCAL(source_repo.base, source_repo.head) == code
    seen = json.loads(source_repo.record.read_text(encoding="utf-8"))
    assert (seen["head"], seen["marker"]) == (source_repo.head, "candidate"), seen
    assert seen["argv"] == ["--local-push-candidate", source_repo.head,
                            "--local-push-base", source_repo.base], seen
    assert seen["git_env"] == [], seen
    assert not any("worktree" in c for c in calls), calls
    monkeypatch.delenv("GIT_DIR")
    assert _git(source_repo.repo, "worktree", "list", "--porcelain").count("worktree ") == 1
    assert (source_repo.repo / "marker.txt").read_text(encoding="utf-8") == "dirty WIP\n"
    assert not source_repo.owned.exists()


def test_run_longer_than_the_old_limit_passes_and_overrun_fails(monkeypatch, source_repo):
    """Guards: a check is bounded by the configured budget, not a fixed 180 s.
    Scaled: a 2 s run passes under a 6 s budget and fails under a 1 s budget."""
    monkeypatch.setenv("STUB_SLEEP", "2")
    monkeypatch.setattr(check, "LOCAL_CANDIDATE_BUDGET_S", 6)
    assert REAL_LOCAL(source_repo.base, source_repo.head) == 0
    source_repo.owned.mkdir()
    monkeypatch.setattr(check, "LOCAL_CANDIDATE_BUDGET_S", 1)
    assert REAL_LOCAL(source_repo.base, source_repo.head) == 1
    assert not source_repo.owned.exists()


def test_budgets_cover_the_measured_policy_test_time():
    """Guards: the real budgets exceed the measured 255.76 s policy-test run."""
    spec = importlib.util.spec_from_file_location(
        "local_push_checks", Path(__file__).parents[1] / "tools/local_push_checks.py")
    local = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(local)
    assert local.COMMAND_BUDGETS_S["compact_tests"] >= 2 * 255.76, local.COMMAND_BUDGETS_S
    assert check.LOCAL_CANDIDATE_BUDGET_S >= sum(local.COMMAND_BUDGETS_S.values()) + 60


def test_base_must_be_an_ancestor_of_the_candidate(monkeypatch, source_repo, capsys):
    assert REAL_LOCAL(source_repo.head, source_repo.base) == 1
    assert "ancestor" in capsys.readouterr().err
    assert source_repo.owned.exists() and not any(source_repo.owned.iterdir())


def test_cleanup_refuses_links_and_retains_export(monkeypatch, tmp_path):
    parent = tmp_path / "owned"
    checkout = parent / "candidate"
    checkout.mkdir(parents=True)
    target = tmp_path / "outside"
    target.mkdir()
    try:
        (checkout / "link").symlink_to(target, target_is_directory=True)
    except OSError:
        pytest.skip("symlink creation not permitted")
    with pytest.raises(ValueError, match="link"):
        check._remove_owned(checkout, parent)
    assert checkout.exists() and target.exists()


def test_missing_base_suggests_fetch_without_checkout(monkeypatch, capsys):
    def run(command, **kw):
        assert 'cat-file' in command
        raise subprocess.CalledProcessError(128, command)
    monkeypatch.setattr(check.subprocess, 'run', run)
    assert REAL_LOCAL(A, B) == 1
    assert 'fetch the remote base commit/history' in capsys.readouterr().err


def test_local_and_remote_checks_use_the_effective_base(monkeypatch):
    events = []
    main_base = "d" * 40
    monkeypatch.setattr(check, "effective_base", lambda base, head: main_base)
    monkeypatch.setattr(check, "run_local_candidate", lambda base, head:
                        events.append(("local", base, head)) or 0)
    monkeypatch.setattr(check.subprocess, "run", lambda cmd, **kw:
                        events.append(("remote", cmd[-3], cmd[-1])) or SimpleNamespace(returncode=0))
    assert check.main([record("refs/heads/branch")]) == 0
    assert events == [("local", main_base, B), ("remote", main_base, B)]


def test_update_that_merges_main_is_compared_with_main(monkeypatch, tmp_path):
    """Guards: merging origin/main into a pushed branch does not count main's changes."""
    repo = tmp_path / "merge-repo"
    repo.mkdir()
    _git(repo, "init", "-q", "-b", "main")
    _git(repo, "config", "user.email", "t@example.invalid")
    _git(repo, "config", "user.name", "t")

    def commit(name):
        (repo / name).write_text(name, encoding="utf-8")
        _git(repo, "add", name)
        _git(repo, "commit", "-q", "-m", name)
        return _git(repo, "rev-parse", "HEAD")

    commit("root")
    _git(repo, "checkout", "-q", "-b", "feature")
    pushed = commit("feature-1")
    plain_update = commit("feature-2")
    _git(repo, "checkout", "-q", "main")
    main_tip = commit("main-2")
    _git(repo, "update-ref", "refs/remotes/origin/main", main_tip)
    _git(repo, "checkout", "-q", "feature")
    _git(repo, "merge", "-q", "--no-ff", "-m", "merge main", "main")
    merged = _git(repo, "rev-parse", "HEAD")
    monkeypatch.setattr(check, "ROOT", repo)
    assert REAL_EFFECTIVE_BASE(plain_update, merged) == main_tip
    assert REAL_EFFECTIVE_BASE(pushed, plain_update) == pushed
    assert REAL_EFFECTIVE_BASE(ZERO, merged) == ZERO
    assert REAL_EFFECTIVE_BASE("c" * 40, merged) == "c" * 40
