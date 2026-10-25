"""Guards: native/CTest failure propagation and scoped build cleanup; signed Cubit
payload bytes and mtimes cannot be replaced by an ordinary build.

Assumes: production Build.ps1 AST/guards execute on isolated fixture paths;
Windows PowerShell 7/cmd.exe, native calls and deletion stubbed, no compilation.
"""

from pathlib import Path
import re
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
TARGETS = (
    "sparsesolv_ngsolve", "axifem",
)


@pytest.mark.parametrize("target", TARGETS)
@pytest.mark.parametrize("exit_code", [0, 7])
def test_native_build_guard_propagates_failure(tmp_path, target, exit_code):
    source = (ROOT / "Build.ps1").read_text(encoding="utf-8-sig")
    match = re.search(
        rf'"\$CMAKE_EXE" --build [^\n]*--target {target} -j \$BuildJobs\n'
        r'(?P<guard>\s*if errorlevel 1 \([^)]*\))', source,
    )
    assert match, f"missing immediate failure guard for {target}"
    guard = match.group("guard")
    assert f"ERROR: {target} build failed" in guard
    assert "exit /b 1" in guard
    _exercise_guard(tmp_path, guard, exit_code)


@pytest.mark.parametrize("target", ["cubit_mesh_curver", "cubit_callback_contract", "cubit_mesh_export_ccm", "ctest"])
@pytest.mark.parametrize("exit_code", [0, 7])
def test_package_native_guards_propagate_failure(tmp_path, target, exit_code):
    script = ROOT / "packages/cubit-mesh-export/tools/build_native.ps1"
    source = script.read_text(encoding="utf-8")
    command = '& $ctest --test-dir' if target == 'ctest' else f'--target {target} -j $BuildJobs'
    line = next(i for i, text in enumerate(source.splitlines()) if command in text)
    guard = source.splitlines()[line + 1]
    assert re.fullmatch(r'\s*if \(\$LASTEXITCODE -ne 0\) \{ throw "[^"]+" \}', guard), guard
    harness = tmp_path / 'package-guard.ps1'
    harness.write_text(f'$global:LASTEXITCODE={exit_code}\n{guard}\n"CONTINUED"\n', encoding='utf-8')
    result = subprocess.run(['pwsh', '-NoProfile', '-File', str(harness)],
                            capture_output=True, encoding='utf-8', errors='replace', timeout=20)
    assert (result.returncode == 0) is (exit_code == 0), result.stderr
    assert ('CONTINUED' in result.stdout) is (exit_code == 0)


def _exercise_guard(tmp_path, guard, exit_code):
    if sys.platform != "win32":
        pytest.skip("cmd.exe execution requires Windows")
    batch = tmp_path / "guard.bat"
    batch.write_text(
        f'@echo off\n"{sys.executable}" -c "import sys; sys.exit({exit_code})"\n'
        f'{guard}\necho CONTINUED\nexit /b 0\n', encoding="ascii",
    )
    result = subprocess.run(
        ["cmd.exe", "/d", "/c", str(batch)], capture_output=True,
        text=True, timeout=15,
    )
    assert result.returncode == (1 if exit_code else 0)
    assert ("CONTINUED" in result.stdout) == (exit_code == 0)


def test_root_cpu_and_package_native_ownership():
    root = (ROOT / "Build.ps1").read_text(encoding="ascii")
    assert "BuildCubitPlugin" not in root and "src\\cubit_plugin" not in root
    for target in ("sparsesolv_ngsolve", "axifem", "peec_matrices", "_radia_pybind"):
        assert target in root, target
    package = (ROOT / "packages/cubit-mesh-export/tools/build_native.ps1").read_text(encoding="utf-8")
    assert "-DCUBIT_CALLBACK_CONTRACT_TESTS=ON" in package
    configs = package.split('& $cmake -G Ninja')[1:]
    assert len(configs) == 2
    for config in configs:
        assert re.search(r'if \(\$LASTEXITCODE -ne 0\) \{ throw "cmake configure', config)
    assert "record" not in package and "Copy-Item" not in package


@pytest.mark.parametrize("radia_only,axi_only,sdk", [
    (False, False, True), (False, False, False),
    (True, False, True), (False, True, True),
])
def test_rebuild_cleans_only_selected_build_directories(tmp_path, radia_only, axi_only, sdk):
    result = _execute_build_scope(tmp_path, "cleanup", radia_only, axi_only, sdk)
    expected = [str((tmp_path / "project/build-msvc").resolve())]
    assert result["removed"] == expected, (result, expected)


def test_ordinary_build_preserves_signed_cubit_payload_bytes_and_mtimes(tmp_path):
    radia_only = False
    import os
    import shutil
    project = tmp_path / "project"
    package = project / "packages/cubit-mesh-export/src/cubit_mesh_export"
    raw = project / "src/cubit_plugin"
    package.mkdir(parents=True)
    # Satisfy the unchanged required Radia-output checks; the ABI invocation
    # is stubbed below, so these sentinels are never loaded as native code.
    (project / "build-msvc").mkdir()
    (project / "src/radia").mkdir(parents=True)
    for name in ("_radia_pybind.cp312-win_amd64.pyd", "peec_matrices.cp312-win_amd64.pyd"):
        (project / "build-msvc" / name).write_bytes(b"required Radia build-output sentinel")
    for directory in (raw / "build-pyd", raw / "build-ccm"):
        directory.mkdir(parents=True)
    files = [package / name for name in ("cubit_mesh_export.ccm", "cubit_mesh_curver.pyd", "native_payloads.json")]
    for path in files:
        path.write_bytes(b"signed release sentinel -- never load as native")
        os.utime(path, (1_600_000_000, 1_600_000_000))
    (raw / "build-pyd/cubit_mesh_curver.cp312-win_amd64.pyd").write_bytes(b"raw candidate")
    (raw / "build-ccm/cubit_mesh_export.ccm").write_bytes(b"raw candidate")
    (raw / "geometry.cpp").write_bytes(b"// source mtime sentinel\n")
    os.utime(raw / "geometry.cpp", (1_700_000_000, 1_700_000_000))
    before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in files}
    # A tiny independent index makes tracked-byte preservation observable.
    git = shutil.which("git")
    assert git, "Git is required for tracked-payload contracts"
    def git_run(*args):
        return subprocess.run([git, *args], cwd=project, capture_output=True,
                              encoding="utf-8", errors="replace", timeout=15)
    assert git_run("init", "-q").returncode == 0
    assert git_run("add", "-f", "--", *[str(p.relative_to(project)) for p in files]).returncode == 0
    result = _execute_build_scope(tmp_path, "payload", radia_only, False, True)
    after = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in files}
    assert after == before, (result, before, after)
    diff = git_run("diff", "--exit-code", "--")
    assert diff.returncode == 0, diff.stdout + diff.stderr


def _execute_build_scope(tmp_path, action, radia_only, axi_only, sdk):
    """Execute the production AST block; replace only native calls/deletion."""
    import json
    import shutil
    if sys.platform != "win32":
        pytest.skip("PowerShell build-scope execution requires Windows")
    pwsh = shutil.which("pwsh")
    assert pwsh, "PowerShell 7 is required for Windows build contracts"
    harness = tmp_path / "build-scope.ps1"
    harness.write_text(r'''
param([string]$BuildScript,[string]$Project,[string]$Action,[int]$Radia,[int]$Axi,[int]$Sdk)
$ErrorActionPreference='Stop'
[Console]::OutputEncoding=[Text.UTF8Encoding]::new($false)
$tokens=$null; $errors=$null
$ast=[System.Management.Automation.Language.Parser]::ParseFile($BuildScript,[ref]$tokens,[ref]$errors)
if ($errors.Count) { throw "Build.ps1 parse errors: $errors" }
$PROJECT_DIR=$Project; $BUILD_DIR=Join-Path $Project 'build-msvc'
$RadiaOnly=[bool]$Radia; $AxiFemOnly=[bool]$Axi; $BuildCubitPlugin=[bool]$Sdk; $Rebuild=$true
$removed=[Collections.Generic.List[string]]::new()
if ($Action -eq 'cleanup') {
    $nodes=@($ast.FindAll({param($n)
        $n -is [System.Management.Automation.Language.IfStatementAst] -and
        $n.Extent.Text -match '^if\s*\(\$Rebuild\)'
    },$true))
    if ($nodes.Count -ne 1) { throw "Expected one cleanup block, got $($nodes.Count)" }
    # No filesystem deletion: observe the actual production deletion requests.
    function Test-Path { param([string]$LiteralPath) return $true }
    function Remove-Item { param([string]$LiteralPath,[switch]$Recurse,[switch]$Force)
        $removed.Add([IO.Path]::GetFullPath($LiteralPath))
    }
} elseif ($Action -eq 'payload') {
    $nodes=@($ast.FindAll({param($n)
        $n -is [System.Management.Automation.Language.IfStatementAst] -and
        $n.Extent.Text -match '^if\s*\(-not\s+\$AxiFemOnly\)' -and
        $n.Extent.Text.Contains('$modules = @(') -and
        $n.Extent.Text.Contains('check_native_module_abi.py')
    },$true))
    if ($nodes.Count -ne 1) { throw "Expected one payload block, got $($nodes.Count)" }
    # The production copy/touch logic executes on isolated fake files. Never
    # invoke a compiler, Python ABI loader, Cubit, installation or deployment.
    $PythonExecutable='Invoke-PythonStub'
    function Invoke-PythonStub { $global:LASTEXITCODE=0 }
    function Write-NativeBuildProvenance { }
} else { throw 'Unknown action' }
& ([scriptblock]::Create($nodes[0].Extent.Text))
$record=@{action=$Action; removed=@($removed.ToArray()); radia=$RadiaOnly; axi=$AxiFemOnly; sdk=$BuildCubitPlugin}
Write-Output ('BUILD_SCOPE_RESULT=' + ($record | ConvertTo-Json -Compress))
''', encoding="utf-8")
    result = subprocess.run(
        [pwsh, "-NoProfile", "-NonInteractive", "-File", str(harness),
         "-BuildScript", str(ROOT / "Build.ps1"), "-Project", str(tmp_path / "project"),
         "-Action", action, "-Radia", str(int(radia_only)), "-Axi", str(int(axi_only)), "-Sdk", str(int(sdk))],
        capture_output=True, encoding="utf-8", errors="replace", timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    lines = [s.removeprefix("BUILD_SCOPE_RESULT=") for s in result.stdout.splitlines()
             if s.startswith("BUILD_SCOPE_RESULT=")]
    assert len(lines) == 1, result.stdout
    return json.loads(lines[0])


def test_builders_write_only_to_declared_roots():
    import importlib.util
    spec = importlib.util.spec_from_file_location('build_write_policy', ROOT / 'tests/build_write_policy.py')
    policy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(policy)
    for relative in ('Build.ps1', 'packages/cubit-mesh-export/tools/build_native.ps1'):
        rows = policy.audit(ROOT / relative)
        assert rows, relative


def test_cme_native_workflow_and_cli_scope(tmp_path):
    import json
    import os
    import yaml
    workflow = yaml.safe_load((ROOT / '.github/workflows/cubit-mesh-export.yml').read_text(encoding='utf-8'))
    jobs = workflow['jobs']
    native = jobs['native-100']
    assert native['runs-on'] == ['self-hosted', 'Windows', 'X64', '100']
    assert native['timeout-minutes'] == 45
    assert native['concurrency'] == {'group': 'cme-native-100', 'cancel-in-progress': False}
    run = '\n'.join(s.get('run','') for s in native['steps'])
    assert 'build_native.ps1 -Rebuild -BuildJobs 4' in run
    assert 'exceeds_pr_budget' in run
    assert 'native-100' in jobs['publish']['needs']
    # Real CLI and real temporary Git history; no live checkout mutation.
    import shutil
    repo = tmp_path / 'repo'
    (repo / 'tools').mkdir(parents=True)
    for name in ('sparsesolv_ci_impact.py','ci_distribution_contract.py','ci_distribution_manifest.json'):
        shutil.copyfile(ROOT/'tools'/name,repo/'tools'/name)
    def git(*args):
        return subprocess.check_output(['git','-C',str(repo),*args],encoding='utf-8',errors='replace').strip()
    git('init','-q')
    git('config','user.name','Contract')
    git('config','user.email','contract@example.invalid')
    git('add','.')
    git('commit','-qm','initial')
    for path, required in [('Build.ps1',False),('CMakeLists.txt',False),
                           ('packages/cubit-mesh-export/tools/build_native.ps1',True),
                           ('src/cubit_plugin/NetgenCurver.cpp',True)]:
        before = git('rev-parse','HEAD')
        file = repo/path
        file.parent.mkdir(parents=True,exist_ok=True)
        file.write_bytes(b'# selector input fixture\n')
        git('add','--',path)
        git('commit','-qm','input')
        out = tmp_path / 'out.txt'
        out.write_bytes(b'')
        event.write_text(json.dumps({'before':before,'ref':'refs/heads/main'}),encoding='utf-8')
        env = os.environ.copy()
        env.update(GITHUB_EVENT_NAME='push',GITHUB_EVENT_PATH=str(event),GITHUB_OUTPUT=str(out))
        result = subprocess.run([sys.executable, str(repo/'tools/sparsesolv_ci_impact.py'), '--lane', 'cubit-mesh-export'],
                                cwd=repo,env=env,capture_output=True,encoding='utf-8',errors='replace',timeout=30)
        assert result.returncode == 0, result.stderr
        assert f'required={str(required).lower()}' in out.read_text(encoding='utf-8'), (path,result.stdout)
