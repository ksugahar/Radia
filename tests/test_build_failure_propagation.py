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
    "cubit_mesh_curver", "cubit_mesh_export_ccm", "cubit_callback_contract",
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


@pytest.mark.parametrize("exit_code", [0, 7])
def test_native_ctest_guard_propagates_failure(tmp_path, exit_code):
    source = (ROOT / "Build.ps1").read_text(encoding="utf-8-sig")
    match = re.search(
        r'"\$CTEST_EXE" --test-dir "%CUBIT_PLUGIN_BUILD%" --output-on-failure\n'
        r'(?P<guard>\s*if errorlevel 1 \([^)]*\))', source,
    )
    assert match, "missing immediate CTest failure guard"
    guard = match.group("guard")
    assert "ERROR: native Cubit callback contract failed" in guard
    assert "exit /b 1" in guard
    _exercise_guard(tmp_path, guard, exit_code)


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


def test_plugin_configuration_failure_cannot_use_stale_build_tree():
    source = (ROOT / "Build.ps1").read_text(encoding="utf-8-sig")
    commands = re.findall(
        r'"\$CMAKE_EXE" -G Ninja [^\n]*"\$PROJECT_DIR\\src\\cubit_plugin"\n'
        r'(?P<guard>\s*if errorlevel 1 \([^)]*\))', source,
    )
    assert len(commands) == 2
    assert all("configuration failed" in guard and "exit /b 1" in guard
               for guard in commands)


def test_optional_cubit_build_does_not_reuse_stale_runner_environment():
    source = (ROOT / "Build.ps1").read_text(encoding="utf-8-sig")
    plugin_block = source[source.index('set "CUBIT_DIR='):source.index(
        'echo Build completed.', source.index('set "CUBIT_DIR='))]
    assert '$DetectedCubitCmakeDir = if ($CubitInstallDir)' in source
    assert '$BuildCubitPlugin = [bool]($CubitInstallDir -and' in source
    assert ('Test-Path (Join-Path $CubitInstallDir '
            '"cmake\\CubitConfig.cmake")') in source
    assert 'set "CUBIT_DIR=$DetectedCubitCmakeDir"' in plugin_block
    assert 'set "CUBIT_DIR=$CubitCmakeDir"' not in plugin_block
    assert 'if /I "$BuildCubitPlugin"=="True" (' in plugin_block
    assert 'if defined CUBIT_DIR' not in plugin_block
    assert plugin_block.count('-DCubit_DIR="$DetectedCubitCmakeDir"') == 2
    assert '-DCubit_DIR="$CubitCmakeDir"' not in plugin_block
    assert '!CUBIT_DIR!' not in plugin_block
    assert '%CUBIT_DIR%' not in plugin_block
    guarded_block = plugin_block[
        plugin_block.index('if /I "$BuildCubitPlugin"=="True" ('):
        plugin_block.index(') else (')
    ]
    for line in guarded_block.splitlines():
        if re.match(r"\s*(?:echo|rem)\b", line, flags=re.IGNORECASE):
            assert not re.search(r"(?<!\^)[()]", line), line


@pytest.mark.parametrize("radia_only,axi_only,sdk", [
    (False, False, True), (False, False, False),
    (True, False, True), (False, True, True),
])
def test_rebuild_cleans_only_selected_build_directories(tmp_path, radia_only, axi_only, sdk):
    result = _execute_build_scope(tmp_path, "cleanup", radia_only, axi_only, sdk)
    expected = [str((tmp_path / "project/build-msvc").resolve())]
    if sdk and not (radia_only or axi_only):
        expected += [str((tmp_path / "project/src/cubit_plugin" / name).resolve())
                     for name in ("build-pyd", "build-ccm")]
    assert result["removed"] == expected, (result, expected)


@pytest.mark.parametrize("radia_only", [False, True])
def test_ordinary_build_preserves_signed_cubit_payload_bytes_and_mtimes(tmp_path, radia_only):
    import os
    import shutil
    project = tmp_path / "project"
    package = project / "packages/cubit-mesh-export/src/cubit_mesh_export"
    raw = project / "src/cubit_plugin"
    package.mkdir(parents=True)
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
