"""Guards: native failure propagation and scoped build cleanup.

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


@pytest.mark.parametrize("radia_only,axi_only,sdk", [
    (False, False, True), (False, False, False),
    (True, False, True), (False, True, True),
])
def test_rebuild_cleans_only_selected_build_directories(tmp_path, radia_only, axi_only, sdk):
    result = _execute_build_scope(tmp_path, "cleanup", radia_only, axi_only, sdk)
    expected = [str((tmp_path / "project/build-msvc").resolve())]
    assert result["removed"] == expected, (result, expected)


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
