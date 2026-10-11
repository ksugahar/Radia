"""Keep the formerly omitted compiled module in the parity inventory."""
import json
import yaml
import subprocess
import shutil
import os
import runpy
from types import SimpleNamespace
import pytest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def matlab_impact_step():
    """The checked workflow step, parsed once for the module."""
    workflow = yaml.safe_load((ROOT / ".github/workflows/sparsesolv.yml").read_text())
    return next(s for s in workflow["jobs"]["ams-regression"]["steps"]
                if s.get("id") == "matlab-impact")


@pytest.fixture(scope="module")
def base_checkout(tmp_path_factory):
    """An initialised repository with an empty base commit, built once.

    Each case copies it and adds its own change commit, so the base
    repository is not rebuilt per case.
    """
    git = shutil.which("git")
    if not git:
        pytest.skip("PowerShell/Git runner contract")
    repo = tmp_path_factory.mktemp("impact-base") / "checkout with spaces"
    repo.mkdir()
    env = {k: v for k, v in os.environ.items()
           if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    for args in (("init",), ("commit", "--allow-empty", "-m", "base")):
        subprocess.run([git, "-C", str(repo), "-c", "user.name=CI test",
                        "-c", "user.email=ci@example.invalid", *args],
                       env=env, check=True, capture_output=True)
    return repo


def test_owned_engine_does_not_inspect_or_attach_to_caller_sessions(monkeypatch):
    runner = runpy.run_path(str(ROOT/"validation_test/ngsolve_matlab_parity/run_sparsesolv_parity.py"))
    def forbidden(*args, **kwargs):
        raise AssertionError("Caller-owned MATLAB sessions must remain untouched")
    monkeypatch.setattr(subprocess, "run", forbidden)
    started = []
    api = SimpleNamespace(find_matlab=forbidden, connect_matlab=forbidden,
                          start_matlab=lambda options: started.append(options) or "owned")
    assert runner["start_owned_engine"](api) == "owned"
    assert started == ["-nodesktop -nosplash -singleCompThread"]

def test_engine_timeout_retains_evidence_when_cleanup_fails(tmp_path, monkeypatch):
    runner = runpy.run_path(str(ROOT/"validation_test/ngsolve_matlab_parity/run_sparsesolv_parity.py"))
    output = tmp_path/"nested/failure.json"

    def fail_cleanup(command, **kwargs):
        assert command == ["taskkill", "/PID", "12345", "/T", "/F"]
        assert json.loads(output.read_text())["passed"] is False
        raise OSError("access denied")

    monkeypatch.setattr(subprocess, "run", fail_cleanup)
    runner["record_timeout"](SimpleNamespace(pid=12345), output)
    record = json.loads(output.read_text())
    assert record["passed"] is False
    assert record["cleanup_error"] == "access denied"


@pytest.mark.parametrize("argv,expected", [([], 900), (["--timeout", "1200"], 1200)])
def test_parity_runner_waits_for_the_worker_with_the_configured_timeout(
        tmp_path, monkeypatch, argv, expected):
    runner = runpy.run_path(str(ROOT/"validation_test/ngsolve_matlab_parity/run_sparsesolv_parity.py"))
    calls = {}

    class Child:
        pid = 4242

        def wait(self, timeout=None):
            calls["timeout"] = timeout
            return 0

    def popen(command, **kwargs):
        calls["command"] = command
        return Child()

    monkeypatch.setattr(subprocess, "Popen", popen)
    monkeypatch.setattr("sys.argv", ["run_sparsesolv_parity.py", "--output",
                                     str(tmp_path/"r.json"), *argv])
    with pytest.raises(SystemExit) as stop:
        runner["main"]()
    assert stop.value.code == 0
    assert calls["timeout"] == expected
    assert "--worker" in calls["command"] and "--timeout" not in calls["command"]


def test_sparsesolv_has_native_and_explicit_fallback_owners():
    manifest = json.loads((ROOT / "matlab/python_api_parity_manifest.json").read_text())
    entry = next(e for e in manifest["binary_extensions"] if e["python"] == "sparsesolv_ngsolve.pyd")
    assert entry["status"] == "focused-native-commands"
    assert (ROOT / entry["fallback"]).is_file()
    source = (ROOT / "src/matlab/radia_mex.cpp").read_text()
    for command in entry["native_commands"]:
        assert f'command == "{command}"' in source
    for owner in entry["matlab"].split("; "):
        assert (ROOT / owner).is_file()


def test_ams_setup_uses_ngsolve_gradient_without_taskmanager():
    source = (ROOT / "src/matlab/radia_mex.cpp").read_text()
    body = source.split("void SparseSolvAMS(")[1].split("void SparseSolvIC(")[0]
    assert "hc->CreateGradient()" in body
    assert "RegionTaskManager" not in body
    assert "a.fespace != space.fespace" in body


def test_matlab_lane_runs_engine_and_retains_json_on_mdx():
    workflow = yaml.safe_load((ROOT/".github/workflows/sparsesolv.yml").read_text())
    job = workflow["jobs"]["ams-regression"]
    assert "100" in job["runs-on"]
    step = next(s for s in job["steps"] if s.get("name") == "Build and verify native MATLAB parity")
    assert step["if"] == (
        "steps.native-impact.outputs.required == 'true' && "
        "steps.matlab-impact.outputs.required == 'true'"
    )
    assert "-MatlabMexOnly" in step["run"]
    assert "run_sparsesolv_parity.py --output" in step["run"]
    assert "sparsesolv-matlab.json" in job["steps"][-1]["with"]["path"]


def test_mex_runtime_boundary_has_native_ci_coverage():
    workflow = yaml.safe_load((ROOT/".github/workflows/sparsesolv.yml").read_text())
    triggers = workflow.get("on", workflow.get(True))
    for event in ("push", "pull_request"):
        for path in ("matlab/+radia/+internal/callMex.m", "matlab/+radia/setup.m",
                     "tests/matlab/test_mex_runtime_setup.m", "tests/matlab/test_beam_transfer_mex.m"):
            assert path in triggers[event]["paths"]
    runner = (ROOT/"validation_test/ngsolve_matlab_parity/run_sparsesolv_parity.py").read_text()
    assert 'root/"tests/matlab/test_mex_runtime_setup.m"' in runner
    assert 'root/"tests/matlab/test_beam_transfer_mex.m"' in runner
    assert 'eng.setenv("RADIA_PYTHON_EXECUTABLE", sys.executable' in runner
    assert "runtests(testfiles)" in runner


def test_missing_diff_base_selects_matlab_without_failing_step(tmp_path, matlab_impact_step):
    pwsh = shutil.which("pwsh")
    if not pwsh:
        pytest.skip("PowerShell runner contract")
    script = matlab_impact_step["run"].replace("${{ github.event_name }}", "pull_request")
    output = tmp_path/"output"
    result = subprocess.run([pwsh, "-NoProfile", "-Command",
        "function git { $global:LASTEXITCODE=1 }; " + script],
        env={**os.environ, "GITHUB_OUTPUT": str(output),
             "GITHUB_EVENT_NAME": "pull_request"}, capture_output=True,
        text=True, encoding="utf-8", errors="replace")
    assert result.returncode == 0, result.stderr
    assert output.read_text().strip() == "required=true"


IMPACT_CASES = [
    ('tools/run_test_tier.py', 'pull_request', 'false'),
    ('matlab/+radia/+sparsesolv/AMS.m', 'pull_request', 'true'),
    ('matlab/+radia/+python/sparsesolv.m', 'push', 'true'),
    ('matlab/+radia/+internal/callMex.m', 'pull_request', 'true'),
    ('matlab/+radia/setup.m', 'push', 'true'),
    ('tests/matlab/test_mex_runtime_setup.m', 'pull_request', 'true'),
    ('docs/intro.md', 'workflow_dispatch', 'true'),
    ('src/core/rad_ngsolve_radia_field.h', 'push', 'true'),
    ('matlab/+radia/RadiaField.m', 'push', 'true'),
    ('tests/matlab/test_radiafield_mex.m', 'push', 'true'),
]

_DRIVER = r'''
$cases = Get-Content -LiteralPath $args[0] -Raw | ConvertFrom-Json
foreach ($case in $cases) {
  $env:GIT_TEST_ASSUME_DIFFERENT_OWNER = '1'
  $env:GIT_CONFIG_NOSYSTEM = '1'
  $env:GIT_CONFIG_GLOBAL = $case.gitconfig
  $env:GITHUB_WORKSPACE = $case.workspace
  $env:GITHUB_OUTPUT = $case.output
  $env:GITHUB_EVENT_NAME = $case.event
  $env:GITHUB_EVENT_PATH = $case.event_path
  Set-Location -LiteralPath $case.cwd
  $global:LASTEXITCODE = 0
  try {
    & $case.script *> $case.log
    $code = $LASTEXITCODE
  } catch {
    $code = 99
    $_ | Out-File -LiteralPath $case.log -Append -Encoding utf8
  }
  Set-Content -LiteralPath $case.code -Value $code -Encoding ascii
}
'''


@pytest.fixture(scope="module")
def impact_runs(tmp_path_factory, matlab_impact_step, base_checkout):
    """Run the workflow step for every case in one PowerShell process.

    Each case gets its own copied checkout, change commit, event file and
    GITHUB_OUTPUT, and the step runs as a script file, as the runner does
    (`pwsh -command ". '{0}'"`), so its final `exit 0` ends that case only.
    """
    pwsh, git = shutil.which('pwsh'), shutil.which('git')
    if not pwsh or not git:
        pytest.skip('PowerShell/Git runner contract')
    env = {k: v for k, v in os.environ.items()
           if k not in ('GIT_DIR', 'GIT_WORK_TREE', 'GIT_INDEX_FILE')}
    root = tmp_path_factory.mktemp('impact-cases')
    cases, runs = [], {}
    for index, (changed, event, _expected) in enumerate(IMPACT_CASES):
        case_dir = root / f'case-{index}'
        repo = case_dir / 'checkout with spaces'
        shutil.copytree(base_checkout, repo)

        def command(*args, repo=repo):
            subprocess.run([git, '-C', str(repo), '-c', 'user.name=CI test',
                            '-c', 'user.email=ci@example.invalid', *args],
                           env=env, check=True, capture_output=True)

        source = repo / changed
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_text('changed\n')
        command('add', changed)
        command('commit', '-m', 'change')
        base = subprocess.check_output([git, '-C', str(repo), 'rev-parse', 'HEAD^'],
                                       env=env, text=True).strip()
        event_path = case_dir / 'event.json'
        event_path.write_text(json.dumps({'before': base}), encoding='utf-8')
        script = case_dir / 'step.ps1'
        script.write_text(matlab_impact_step['run'].replace('${{ github.event_name }}', event),
                          encoding='utf-8')
        record = {'cwd': str(case_dir), 'workspace': str(repo), 'event': event,
                  'output': str(case_dir / 'github-output'),
                  'event_path': str(event_path), 'gitconfig': str(case_dir / 'empty-gitconfig'),
                  'script': str(script), 'log': str(case_dir / 'step.log'),
                  'code': str(case_dir / 'step.code')}
        untrusted = subprocess.run(
            [git, '-C', str(repo), 'rev-parse', '--show-toplevel'],
            env={**env, 'GIT_TEST_ASSUME_DIFFERENT_OWNER': '1', 'GIT_CONFIG_NOSYSTEM': '1',
                 'GIT_CONFIG_GLOBAL': record['gitconfig']},
            capture_output=True, text=True, encoding='utf-8', errors='replace')
        cases.append(record)
        runs[(changed, event)] = (record, untrusted)
    plan = root / 'cases.json'
    plan.write_text(json.dumps(cases), encoding='utf-8')
    driver = root / 'driver.ps1'
    driver.write_text(_DRIVER, encoding='utf-8')
    subprocess.run([pwsh, '-NoProfile', '-File', str(driver), str(plan)], cwd=root, env=env,
                   capture_output=True, text=True, timeout=120,
                   encoding='utf-8', errors='replace')
    return runs


@pytest.mark.parametrize('changed,event,expected', IMPACT_CASES)
def test_impact_uses_checkout_even_outside_repository(changed, event, expected, impact_runs):
    record, untrusted = impact_runs[(changed, event)]
    assert untrusted.returncode != 0 and 'dubious ownership' in untrusted.stderr
    log = Path(record['log']).read_text(encoding='utf-8', errors='replace') \
        if Path(record['log']).exists() else ''
    assert Path(record['code']).read_text().strip() == '0', log
    assert Path(record['output']).read_text().strip() == f'required={expected}'
