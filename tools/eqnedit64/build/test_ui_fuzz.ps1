# Seeded random operations through the real window procedures.
#
# All seeds run in ONE hidden process. Each seed resets the editor model and is
# logged before execution, preserving deterministic replay while avoiding a
# private-font register/unregister cycle for every seed.

param(
    [string]$AppPath,
    [ValidateRange(1, 10000)]
    [int]$Seeds = 24,
    [ValidateRange(1, 1000000)]
    [int]$Ops = 3000,
    [ValidateRange(1, 3600)]
    [int]$TimeoutSeconds = 120
)

$ErrorActionPreference = 'Stop'

$root = Split-Path -Parent $PSScriptRoot
$exe = if ($AppPath) {
    (Resolve-Path -LiteralPath $AppPath -ErrorAction Stop).Path
} else {
    Join-Path $root 'build\Eqnedit64.exe'
}
if (-not $AppPath -and -not (Test-Path $exe)) {
    $exe = Join-Path $root 'dist\Eqnedit64.exe'
}
if (-not (Test-Path $exe)) {
    Write-Output 'skip  Eqnedit64.exe has not been built'
    return
}

$tempRoot = 'C:\temp'
New-Item -ItemType Directory -Path $tempRoot -Force | Out-Null
$log = Join-Path $tempRoot ("eqn-ui-fuzz-" + [guid]::NewGuid().ToString('N') + '.log')
$sessionId = (Get-Process -Id $PID).SessionId
function Get-SessionFontDriverHostIds {
    return (@(Get-Process -Name fontdrvhost -ErrorAction SilentlyContinue |
        Where-Object { $_.SessionId -eq $sessionId } |
        Sort-Object Id | ForEach-Object { $_.Id }) -join ',')
}

$fontHostsBefore = Get-SessionFontDriverHostIds
$process = Start-Process -FilePath $exe `
    -ArgumentList @('--ui-fuzz-batch', $Seeds, $Ops) `
    -PassThru -WindowStyle Hidden -RedirectStandardError $log
$finished = $process.WaitForExit($TimeoutSeconds * 1000)
$progress = @(Get-Content -LiteralPath $log -ErrorAction SilentlyContinue)
if (-not $finished) {
    try { $process.Kill() } catch { }
    $process.WaitForExit(5000) | Out-Null
    $last = $progress | Select-Object -Last 1
    Remove-Item -LiteralPath $log -ErrorAction SilentlyContinue
    throw "GUI fuzz timed out after ${TimeoutSeconds}s; last progress: $last"
}

$fontHostsAfter = Get-SessionFontDriverHostIds
$last = $progress | Select-Object -Last 1
Remove-Item -LiteralPath $log -ErrorAction SilentlyContinue
if ($process.ExitCode -ne 0) {
    throw ("GUI fuzz exit 0x{0:X8}; last progress: {1}" -f `
           $process.ExitCode, $last)
}
if ($fontHostsBefore -and $fontHostsAfter -ne $fontHostsBefore) {
    throw ("Windows fontdrvhost restarted during structural fuzz: " +
           "$fontHostsBefore -> $fontHostsAfter")
}

Write-Output ("ok    {0} seeds x {1} operations through one real-window process" -f `
              $Seeds, $Ops)
