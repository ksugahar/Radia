# A running Windows executable cannot be replaced in place.  The laboratory
# build policy explicitly permits terminating repository-local Eqnedit64
# instances before a rebuild, so target only these two resolved paths and
# leave deployment copies and every unrelated process alone.

$ErrorActionPreference = 'Stop'
$root = [IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$targets = @(
    [IO.Path]::GetFullPath((Join-Path $root 'build\Eqnedit64.exe')),
    [IO.Path]::GetFullPath((Join-Path $root 'dist\Eqnedit64.exe'))
)

$stopped = 0
foreach ($process in @(Get-Process -Name 'Eqnedit64' -ErrorAction SilentlyContinue)) {
    try {
        $path = [IO.Path]::GetFullPath($process.Path)
    } catch {
        continue
    }
    if (-not ($targets -ccontains $path)) { continue }
    Write-Host "[INFO] stopping repository Eqnedit64 process $($process.Id): $path"
    Stop-Process -Id $process.Id -Force
    Wait-Process -Id $process.Id -Timeout 10 -ErrorAction SilentlyContinue
    $stopped++
}

if ($stopped) {
    Write-Host "[OK] stopped $stopped repository-local Eqnedit64 instance(s)"
}
