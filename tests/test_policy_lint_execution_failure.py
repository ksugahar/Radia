"""A lint that could not run has not passed.

`tools/policy_lint.py` shells out to git for most of its checks. It used to
map a git failure -- exit 128, git missing, a bad pathspec -- onto "no
matches", which every caller then read as "no violation". Simulating exit
code 128 made every policy PASS. Whatever prevents a check from running must
make the run fail, loudly, and this pins that in both directions: the failure
path fails, and the normal path still passes.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load():
    spec = importlib.util.spec_from_file_location(
        "policy_lint_under_test", ROOT / "tools" / "policy_lint.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("rc", [128, 2, 255])
def test_a_failing_git_fails_the_lint(monkeypatch, capsys, rc):
    """Exit code 128 from git must not read as a clean repository."""
    lint = _load()
    monkeypatch.setattr(lint, "_sh", lambda cmd: (rc, ""))
    assert lint.main([]) == 1
    captured = capsys.readouterr()
    assert "could not run" in captured.err
    assert "PASS" not in captured.out


def test_a_missing_git_fails_the_lint(monkeypatch, capsys):
    """The subprocess itself failing is the same situation, one layer down."""
    lint = _load()

    def missing(cmd):
        raise FileNotFoundError("git")

    monkeypatch.setattr(lint, "_sh", missing)
    assert lint.main([]) == 1
    assert "could not run" in capsys.readouterr().err


def test_grep_and_ls_files_distinguish_no_match_from_failure():
    """git grep exits 1 for no matches; that is a result, not an error."""
    lint = _load()
    with pytest.raises(lint.LintExecutionError):
        lint._sh = lambda cmd: (128, "")
        lint._git_grep("anything", ())
    lint._sh = lambda cmd: (1, "")
    assert lint._git_grep("anything", ()) == []
    with pytest.raises(lint.LintExecutionError):
        lint._sh = lambda cmd: (2, "")
        lint._git_ls_files("*.pyd")
    lint._sh = lambda cmd: (0, "a\n\nb\n")
    assert lint._git_ls_files("*.pyd") == ["a", "b"]


def test_the_real_lint_still_passes_on_this_tree():
    """The fix must not have turned an honest pass into a failure."""
    lint = _load()
    assert lint.main(["--quiet"]) == 0


def test_git_trust_is_scoped_to_the_active_checkout(monkeypatch):
    from types import SimpleNamespace
    lint = _load()
    calls = []
    def run(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(returncode=128, stdout="", stderr="git diagnostic")
    monkeypatch.setattr(lint.subprocess, "run", run)
    assert lint._sh(["git", "grep", "needle"]) == (128, "git diagnostic")
    assert calls[0][0] == ["git", "-c", f"safe.directory={lint.REPO}", "grep", "needle"]
    assert calls[0][1]["cwd"] == lint.REPO

@pytest.mark.parametrize('private_text', [
    '10.' + '99.88.77',
    'C:' + '/Users/SyntheticOperator/scratch',
    chr(92) * 2 + 'synthetic-server' + chr(92) + 'share',
])
def test_exporter_release_skill_is_scanned_for_runtime_identity(monkeypatch, tmp_path, private_text):
    """The release guide must not regain a blanket operational exemption."""
    lint = _load()
    relative = '.agents/skills/release-cubit-mesh-export/SKILL.md'
    source = tmp_path / relative
    source.parent.mkdir(parents=True)
    source.write_text('Run on ' + private_text, encoding='utf-8')
    monkeypatch.setattr(lint, 'REPO', str(tmp_path))
    monkeypatch.setattr(lint, '_git_grep', lambda *a, **k: [relative + ':1:synthetic identity'])
    monkeypatch.setattr(lint, '_git_ls_files', lambda *a: [relative])
    monkeypatch.setattr(lint, 'check_bare_mdx_growth', lambda: [])
    assert lint.check_public_runtime_privacy()
    source.write_text('Use the editable development host and wheel acceptance host. '
                      'Resolve SSH aliases and absolute paths from private configuration.',
                      encoding='utf-8')
    assert lint.check_public_runtime_privacy() == []
