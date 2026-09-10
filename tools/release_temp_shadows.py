"""Retire the explicitly named Omega override after its main integration."""
import json
import subprocess

HOSTS = ("mdx1", "mdx2", "hibino")
SCRIPT = r'''
$ErrorActionPreference = 'Stop'
$target = 'C:\temp\radia-omega-test'
$apply = __APPLY__
$result = [ordered]@{target=$target; passed=$false; removed=$false; blockers=@()}
try {
    $expected = [IO.Path]::GetFullPath($target)
    if ($expected -ne 'C:\temp\radia-omega-test') { throw 'Unexpected cleanup target' }
    $parent = Get-Item -LiteralPath 'C:\temp' -Force
    if ($parent.Attributes -band [IO.FileAttributes]::ReparsePoint) {
        throw 'Cleanup parent is a reparse point'
    }
    foreach ($scope in @('Process','User','Machine')) {
        $value = [Environment]::GetEnvironmentVariable('PYTHONPATH', $scope)
        if ($value -and $value.Replace('/','\').ToLowerInvariant().Contains($target.ToLowerInvariant())) {
            $result.blockers += "PYTHONPATH still references override ($scope)"
        }
    }
    if (Test-Path -LiteralPath $target) {
        $items = @(Get-Item -LiteralPath $target -Force)
        $items += @(Get-ChildItem -LiteralPath $target -Force -Recurse)
        if (@($items | Where-Object { $_.Attributes -band [IO.FileAttributes]::ReparsePoint }).Count) {
            $result.blockers += 'Reparse point in cleanup tree'
        }
        $active = @(Get-CimInstance Win32_Process | Where-Object {
            $_.Name -match '^(python(w)?|MATLAB)\.exe$'
        } | Select-Object ProcessId, Name)
        if ($active.Count) {
            $result.blockers += 'Python/MATLAB processes active; cannot prove override unused'
            $result.active_processes = $active
        }
        if ($apply -and -not $result.blockers.Count) {
            Remove-Item -LiteralPath 'C:\temp\radia-omega-test' -Recurse -Force
            $result.removed = $true
        }
    }
    $result.exists = Test-Path -LiteralPath $target
    $result.passed = -not $result.exists -and -not $result.blockers.Count
} catch { $result.error = $_.Exception.Message }

$result | ConvertTo-Json -Depth 6 -Compress

'''


def inspect_shadows(apply=False):
    results = {}
    for host in HOSTS:
        try:
            result = subprocess.run(
                ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", host,
                 "pwsh -NoProfile -NonInteractive -ExecutionPolicy Bypass -Command -"],
                input=SCRIPT.replace("__APPLY__", "$true" if apply else "$false"),
                text=True, capture_output=True, timeout=90, check=True)
            lines = [line for line in result.stdout.splitlines() if line.startswith("{")]
            if not lines:
                raise ValueError("Remote cleanup returned no JSON report")
            report = json.loads(lines[-1])
            if not isinstance(report, dict) or not isinstance(report.get("passed"), bool):
                raise ValueError("Remote cleanup returned an invalid report")
            results[host] = report
        except (OSError, subprocess.SubprocessError, ValueError, IndexError) as exc:
            results[host] = {"passed": False, "error": str(exc)}
    return results
