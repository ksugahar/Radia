"""Focused tests for the 100 pre-push candidate boundary."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _load_module():
    path = ROOT / "tools" / "ci_preflight_remote.py"
    spec = importlib.util.spec_from_file_location("ci_preflight_remote", path)
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


@pytest.mark.parametrize('mapping,host,expected', [
    ('', '100', '100'),
    ('100=compute', '100', 'compute'),
    (' 100 = compute ', '100', 'compute'),
])
def test_ssh_destination_mapping(mapping, host, expected):
    assert _load_module().resolve_ssh_destination(host, mapping) == expected


def test_unset_ssh_mapping_keeps_runner_label(monkeypatch):
    monkeypatch.delenv('RADIA_PREFLIGHT_SSH_ALIAS', raising=False)
    assert _load_module().resolve_ssh_destination('100') == '100'


@pytest.mark.parametrize('mapping', [
    '100', '1001=', '=100', '1001=100,', '1001=100,1001=1002',
    'other=100', '1001=-F', '1001=100 value', '1001=100:22', '1001=100=more',
])
def test_invalid_ssh_mapping_fails_loudly(mapping):
    with pytest.raises(RuntimeError, match='Invalid RADIA_PREFLIGHT_SSH_ALIAS'):
        _load_module().resolve_ssh_destination('100', mapping)


def _mock_candidate(monkeypatch, module, tmp_path):
    monkeypatch.delenv('RADIA_PREFLIGHT_SSH_ALIAS', raising=False)
    monkeypatch.setattr(module.subprocess, 'check_output', lambda *a, **kw: pytest.fail('Preflight must not query GitHub runner inventory'))
    monkeypatch.setattr(module, 'create_bundle', lambda *a: 'a' * 40)
    real_tempdir = module.tempfile.TemporaryDirectory
    monkeypatch.setattr(module.tempfile, 'TemporaryDirectory', lambda **kw:
        real_tempdir(dir=tmp_path))


def test_retired_host_mapping_fails_before_creating_bundle(monkeypatch, tmp_path, capsys):
    module = _load_module()
    _mock_candidate(monkeypatch, module, tmp_path)
    monkeypatch.setenv('RADIA_PREFLIGHT_SSH_ALIAS', '1002=1002')
    def forbidden(*args, **kwargs):
        pytest.fail('Incomplete mapping must fail before bundle or transport work')
    monkeypatch.setattr(module, 'create_bundle', forbidden)
    monkeypatch.setattr(module, 'run', forbidden)
    assert module.main(['--base', 'a' * 40, '--head', 'b' * 40]) == 1
    assert 'Invalid RADIA_PREFLIGHT_SSH_ALIAS mapping' in capsys.readouterr().err


def test_alias_used_for_setup_upload_execution_and_cleanup(monkeypatch, tmp_path, capsys):
    module = _load_module()
    _mock_candidate(monkeypatch, module, tmp_path)
    monkeypatch.setenv('RADIA_PREFLIGHT_SSH_ALIAS', '100=compute')
    commands = []
    monkeypatch.setattr(module, 'run', lambda command, **kw: commands.append(command))
    assert module.main(['--base', 'a' * 40, '--head', 'b' * 40, '--host', '100']) == 0
    assert [command[1] for command in commands if command[0] == 'ssh'] == ['compute'] * 3
    upload = next(command for command in commands if command[0] == 'scp')
    assert upload[2].startswith('compute:')
    assert 'Preflight host: 100' in capsys.readouterr().out


def test_unreachable_alias_aborts_without_label_fallback(monkeypatch, tmp_path, capsys):
    module = _load_module()
    _mock_candidate(monkeypatch, module, tmp_path)
    monkeypatch.setenv('RADIA_PREFLIGHT_SSH_ALIAS', '100=compute')
    commands = []
    def failed_connection(command, **kwargs):
        commands.append(command)
        raise RuntimeError('SSH connection failed')
    monkeypatch.setattr(module, 'run', failed_connection)
    assert module.main(['--base', 'a' * 40, '--head', 'b' * 40]) == 1
    assert commands and all(command[:2] == ['ssh', 'compute'] for command in commands)
    assert '100 preflight failed: SSH connection failed' in capsys.readouterr().err


def test_partial_upload_is_cleaned_on_failure(monkeypatch, tmp_path):
    module = _load_module()
    _mock_candidate(monkeypatch, module, tmp_path)
    scripts = []
    monkeypatch.setattr(module, 'remote_command', lambda script, host: scripts.append(script))
    def fail_upload(*args, **kwargs):
        raise RuntimeError('upload failed')
    monkeypatch.setattr(module, 'run', fail_upload)
    assert module.main(['--base', 'a' * 40, '--head', 'b' * 40, '--host', '100']) == 1
    assert scripts[-1].startswith('if (Test-Path -LiteralPath ')
    assert scripts[-1].endswith('exit 0')
    assert '.bundle' in scripts[-1]


@pytest.mark.parametrize('busy', [False, True])
def test_generated_script_locks_before_setup_and_always_cleans(monkeypatch, tmp_path, busy):
    module = _load_module()
    _mock_candidate(monkeypatch, module, tmp_path)
    scripts = []
    monkeypatch.setattr(module, 'remote_command', lambda script, host: scripts.append(script))
    monkeypatch.setattr(module, 'run', lambda *a, **kw: None)
    assert module.main(['--base', 'a' * 40, '--head', 'b' * 40, '--host', '100']) == 0
    script = next(s for s in scripts if '$lock = $null' in s)
    assert script.index('try {') < script.index('[IO.File]::Open') < script.index('& $git clone')
    assert script.index('Runner.Worker') < script.index('& $git clone')
    assert '[IO.FileShare]::None' in script
    assert script.index('} finally {') > script.index('Export-Candidate $git $repo $worktree $head')
    assert script.rindex('Remove-CandidateExport $root $worktree') > script.index('Export-Candidate $git')
    assert 'worktree add' not in script and 'worktree remove' not in script
    assert 'Set-Content -LiteralPath $alternates' not in script
    assert 'if ($exported -ne $head)' in script
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
        isolated = isolated.replace(r'C:\Program Files\Git\cmd\git.exe',
                                    str(tmp_path / 'missing-git.exe'))
        bundle_line = next(line for line in isolated.splitlines() if line.startswith('$bundle = '))
        bundle = Path(bundle_line.split("'", 2)[1])
        bundle.parent.mkdir(parents=True, exist_ok=True)
        bundle.write_bytes(b'partial upload')
        # Never depend on the host runner (including the one executing this test).
        process_result = "[pscustomobject]@{ Id = 123 }" if busy else "$null"
        isolated = (
            "function Get-Process { param($Name, $ErrorAction) "
            f"return {process_result} " + "}\n" + isolated
        )
        file.write_text(isolated, encoding='utf-8')
        # pwsh writes its failure text in the console code page (cp932 on a
        # Japanese Windows), so a strict utf-8 decode raised inside the reader
        # thread and left result.stderr as None.
        result = subprocess.run([pwsh, '-NoProfile', '-File', str(file)],
                                capture_output=True, text=True, timeout=30,
                                encoding='utf-8', errors='replace')
        assert result.returncode != 0
        expected = 'build runner is occupied' if busy else 'Git is unavai100le'
        assert expected in result.stderr
        assert not bundle.exists()


def test_preflight_refuses_other_roles_without_transport():
    with pytest.raises(RuntimeError, match="build target"):
        _load_module().resolve_ssh_destination("validation", "")


def _pwsh_export_fixture(tmp_path):
    """A cached repository with one commit and the owned work root, all scratch."""
    import shutil
    pwsh, git = shutil.which('pwsh'), shutil.which('git')
    if not (pwsh and git):
        pytest.skip('pwsh and git are required for the generated-PowerShell export test')
    repo = tmp_path / 'root' / 'repo'
    repo.mkdir(parents=True)
    env = {k: v for k, v in __import__('os').environ.items() if not k.startswith('GIT_')}
    def g(*args):
        return subprocess.run([git, '-c', f'safe.directory={repo}', *args], cwd=repo, env=env,
                              check=True, capture_output=True, text=True).stdout.strip()
    g('init', '-q'); g('config', 'user.email', 't@example.invalid'); g('config', 'user.name', 't')
    (repo / 'a.txt').write_text('candidate\n', encoding='utf-8')
    g('add', '-A'); g('commit', '-q', '-m', 'c')
    (tmp_path / 'root' / 'work').mkdir()
    return pwsh, git, env, tmp_path / 'root', repo, g('rev-parse', 'HEAD')


def _run_export(tmp_path, pwsh, git, env, root, repo, head, export, probe=''):
    module = _load_module()
    body = (module.EXPORT_FUNCTIONS + "\n$ErrorActionPreference = 'Stop'\n$script:exportOwned = $false\n"
            f"try {{\n  Export-Candidate '{git}' '{repo}' '{export}' '{head}'\n{probe}\n}} finally {{\n"
            f"  Remove-CandidateExport '{root}' '{export}'\n}}\n")
    file = tmp_path / 'export.ps1'
    file.write_text(body, encoding='utf-8')
    return subprocess.run([pwsh, '-NoProfile', '-File', str(file)], capture_output=True, text=True,
                          encoding='utf-8', errors='replace', timeout=120, env=env)


def test_generated_export_runs_on_real_git_and_cleans_up(tmp_path):
    """Guards: the build-host export reads the candidate through LF alternates,
    verifies HEAD, is not a registered worktree, and is removed afterwards."""
    pwsh, git, env, root, repo, head = _pwsh_export_fixture(tmp_path)
    export = root / 'work' / 'run1'
    record = tmp_path / 'probe.txt'
    probe = ("  $alt = Join-Path $args_export '.git/objects/info/alternates'\n"
             .replace('$args_export', f"'{export}'")
             + "  $bytes = [IO.File]::ReadAllBytes($alt)\n"
             + f"  $headNow = (& '{git}' -c 'safe.directory={export}' -C '{export}' rev-parse HEAD).Trim()\n"
             + f"  Set-Content -LiteralPath '{record}' -Value ('{{0}} {{1}} {{2}}' -f $bytes[-1], ($bytes -contains 13), $headNow)\n")
    result = _run_export(tmp_path, pwsh, git, env, root, repo, head, export, probe)
    assert result.returncode == 0, result.stderr
    last, has_cr, seen = record.read_text(encoding='utf-8').split()
    assert (last, has_cr, seen) == ('10', 'False', head)
    assert not export.exists()
    listed = subprocess.run([git, '-c', f'safe.directory={repo}', 'worktree', 'list', '--porcelain'],
                            cwd=repo, env=env, capture_output=True, text=True, check=True).stdout
    assert listed.count('worktree ') == 1, listed


def test_generated_export_never_removes_a_preexisting_path(tmp_path):
    """Guards: a path that existed before this run is refused and kept, sentinel intact."""
    pwsh, git, env, root, repo, head = _pwsh_export_fixture(tmp_path)
    export = root / 'work' / 'run1'
    export.mkdir()
    sentinel = export / 'sentinel.txt'
    sentinel.write_text('keep', encoding='utf-8')
    result = _run_export(tmp_path, pwsh, git, env, root, repo, head, export)
    assert result.returncode != 0
    assert 'already exists' in result.stderr, result.stderr
    assert sentinel.read_text(encoding='utf-8') == 'keep'


def test_forwarded_env_lines_pass_locations_only():
    import importlib.util
    from pathlib import Path
    spec = importlib.util.spec_from_file_location(
        "ci_preflight_remote_fwd", Path(__file__).parents[1] / "tools/ci_preflight_remote.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    lines = module.forwarded_env_lines({
        "RADIA_PUBLIC_NAME_RULES_FILE": "D:/rules/o'k.json",
        "RADIA_OPTUNA_SOURCE": "  ",
        "UNRELATED": "x",
    })
    assert lines == ["$env:RADIA_PUBLIC_NAME_RULES_FILE = 'D:/rules/o''k.json'"], lines


def test_failed_remote_launch_does_not_disclose_encoded_script(monkeypatch):
    """Guards: a failing ssh launch reports status only, never the payload."""
    import base64
    import subprocess as sp
    module = _load_module()
    marker = "SYNTHETIC-PRIVATE-LOCATION-7f3c"
    script = "\n".join(module.forwarded_env_lines(
        {"RADIA_PUBLIC_NAME_RULES_FILE": f"D:/{marker}/rules.json"}))
    encoded = base64.b64encode(script.encode("utf-16le")).decode("ascii")
    monkeypatch.setattr(module.subprocess, "run",
                        lambda *a, **k: sp.CompletedProcess(a[0], 1))
    with pytest.raises(RuntimeError) as excinfo:
        module.remote_command(script, "100")
    message = str(excinfo.value)
    assert "failed (1)" in message, message
    assert marker not in message and encoded not in message, message
