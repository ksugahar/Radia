"""Focused tests for the mdx pre-push candidate boundary."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _load_module():
    path = ROOT / "tools" / "ci_preflight_mdx.py"
    spec = importlib.util.spec_from_file_location("ci_preflight_mdx", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _git(repo: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=repo, text=True, encoding="utf-8"
    ).strip()


def test_new_remote_branch_uses_main_merge_base(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-b", "main")
    _git(repo, "config", "user.email", "ci@example.invalid")
    _git(repo, "config", "user.name", "CI Test")
    (repo / "base.txt").write_text("base\n", encoding="utf-8")
    _git(repo, "add", "base.txt")
    _git(repo, "commit", "-m", "base")
    base = _git(repo, "rev-parse", "HEAD")
    _git(repo, "update-ref", "refs/remotes/origin/main", base)
    _git(repo, "switch", "-c", "feature")
    (repo / "feature.txt").write_text("feature\n", encoding="utf-8")
    _git(repo, "add", "feature.txt")
    _git(repo, "commit", "-m", "feature")
    head = _git(repo, "rev-parse", "HEAD")

    module = _load_module()
    module.ROOT = repo

    assert module.resolve_candidate_base(module.ZERO, head) == base


def _runner(host, *, busy=False, status='online'):
    return {'labels': [{'name': 'mdx'}, {'name': host}], 'busy': busy, 'status': status}


@pytest.mark.parametrize('runner', [_runner('mdx1', busy=True), _runner('mdx1', status='offline')])
def test_explicit_host_cannot_bypass_availability(runner):
    module = _load_module()
    with pytest.raises(RuntimeError, match='No idle'):
        module.select_idle_host([runner, _runner('mdx2')], 'mdx1')


def test_auto_selects_idle_host():
    module = _load_module()
    assert module.select_idle_host([_runner('mdx1', busy=True), _runner('mdx2')]) == 'mdx2'


def _mock_candidate(monkeypatch, module, tmp_path):
    monkeypatch.setattr(module.subprocess, 'check_output', lambda *a, **kw:
        json.dumps({'runners': [_runner('mdx1')]}))
    monkeypatch.setattr(module, 'create_bundle', lambda *a: 'a' * 40)
    real_tempdir = module.tempfile.TemporaryDirectory
    monkeypatch.setattr(module.tempfile, 'TemporaryDirectory', lambda **kw:
        real_tempdir(dir=tmp_path))


def test_partial_upload_is_cleaned_on_failure(monkeypatch, tmp_path):
    module = _load_module()
    _mock_candidate(monkeypatch, module, tmp_path)
    scripts = []
    monkeypatch.setattr(module, 'remote_command', lambda script, host: scripts.append(script))
    def fail_upload(*args, **kwargs):
        raise RuntimeError('upload failed')
    monkeypatch.setattr(module, 'run', fail_upload)
    assert module.main(['--base', 'a' * 40, '--head', 'b' * 40, '--host', 'mdx1']) == 1
    assert scripts[-1].startswith('Remove-Item -LiteralPath ')
    assert '.bundle' in scripts[-1]


def test_generated_script_locks_before_setup_and_always_cleans(monkeypatch, tmp_path):
    module = _load_module()
    _mock_candidate(monkeypatch, module, tmp_path)
    scripts = []
    monkeypatch.setattr(module, 'remote_command', lambda script, host: scripts.append(script))
    monkeypatch.setattr(module, 'run', lambda *a, **kw: None)
    assert module.main(['--base', 'a' * 40, '--head', 'b' * 40, '--host', 'mdx1']) == 0
    script = next(s for s in scripts if '$lock = $null' in s)
    assert script.index('try {') < script.index('[IO.File]::Open') < script.index('& $git clone')
    assert script.index('Runner.Worker') < script.index('& $git clone')
    assert '[IO.FileShare]::None' in script
    assert script.index('} finally {') > script.index('worktree add')
    assert '$lock.Dispose()' in script
    # Parse the generated PowerShell without executing any remote command.
    import shutil
    pwsh = shutil.which('pwsh')
    if pwsh:
        file = tmp_path / 'probe.ps1'
        file.write_text(script, encoding='utf-8')
        command = '$e=$null; $t=$null; [void][System.Management.Automation.Language.Parser]::ParseFile($args[0],[ref]$t,[ref]$e); if($e.Count){$e | Out-String | Write-Error; exit 1}'
        parser = tmp_path / 'parse.ps1'
        parser.write_text(command, encoding='utf-8')
        subprocess.run([pwsh, '-NoProfile', '-File', str(parser), str(file)], check=True, timeout=30)
        # Exercise setup failure with only this test's scratch paths and no Git/network work.
        isolated = script.replace(r'C:\temp\radia-preflight', str(tmp_path))
        isolated = isolated.replace(r'C:\actions-runner\tools\PortableGit\bin\git.exe',
                                    str(tmp_path / 'missing-git.exe'))
        bundle_line = next(line for line in isolated.splitlines() if line.startswith('$bundle = '))
        bundle = Path(bundle_line.split("'", 2)[1])
        bundle.parent.mkdir(parents=True, exist_ok=True)
        bundle.write_bytes(b'partial upload')
        file.write_text(isolated, encoding='utf-8')
        result = subprocess.run([pwsh, '-NoProfile', '-File', str(file)],
                                capture_output=True, text=True, timeout=30)
        assert result.returncode != 0
        assert 'Git is unavailable' in result.stderr
        assert not bundle.exists()
