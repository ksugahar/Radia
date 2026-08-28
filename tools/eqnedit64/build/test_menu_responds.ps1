# Every menu item must return, through the real WM_COMMAND path.
#
# The application exercises the complete menu in ONE hidden process and
# dismisses its own modal dialogs. Progress is flushed before each item, so a
# timeout names the last command without creating hundreds of short-lived font
# registrations in the user's Windows session.

param(
    [string]$AppPath,
    [ValidateRange(1, 3600)]
    [int]$TimeoutSeconds = 180
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
    exit 0
}

$tempRoot = 'C:\temp'
New-Item -ItemType Directory -Path $tempRoot -Force | Out-Null
$log = Join-Path $tempRoot ("eqn-menu-" + [guid]::NewGuid().ToString('N') + '.log')
$sessionId = (Get-Process -Id $PID).SessionId
function Get-SessionFontDriverHostIds {
    return (@(Get-Process -Name fontdrvhost -ErrorAction SilentlyContinue |
        Where-Object { $_.SessionId -eq $sessionId } |
        Sort-Object Id | ForEach-Object { $_.Id }) -join ',')
}

$fontHostsBefore = Get-SessionFontDriverHostIds
$process = Start-Process -FilePath $exe -ArgumentList '--menu-responds-all' `
    -PassThru -WindowStyle Hidden -RedirectStandardOutput $log
$finished = $process.WaitForExit($TimeoutSeconds * 1000)
$lines = @(Get-Content -LiteralPath $log -Encoding Unicode -ErrorAction SilentlyContinue)
if (-not $finished) {
    try { $process.Kill() } catch { }
    $process.WaitForExit(5000) | Out-Null
    $last = $lines | Where-Object { $_ -like 'BEGIN*' } | Select-Object -Last 1
    Remove-Item -LiteralPath $log -ErrorAction SilentlyContinue
    Write-Output ("FAIL  menu command did not return within {0}s: {1}" -f `
                  $TimeoutSeconds, $last)
    exit 1
}

$fontHostsAfter = Get-SessionFontDriverHostIds
$begun = @($lines | Where-Object { $_ -like 'BEGIN*' }).Count
$done = @($lines | Where-Object { $_ -like 'DONE*' }).Count
$lastBegin = $lines | Where-Object { $_ -like 'BEGIN*' } | Select-Object -Last 1
Remove-Item -LiteralPath $log -ErrorAction SilentlyContinue

if ($process.ExitCode -ne 0) {
    Write-Output ("FAIL  menu suite exit {0}; last command: {1}" -f `
                  $process.ExitCode, $lastBegin)
    exit 1
}
if ($begun -eq 0 -or $done -ne $begun) {
    Write-Output ("FAIL  menu progress incomplete: {0} begun, {1} completed" -f `
                  $begun, $done)
    exit 1
}
if ($fontHostsBefore -and $fontHostsAfter -ne $fontHostsBefore) {
    Write-Output ("FAIL  Windows fontdrvhost restarted during menu test: " +
                  "$fontHostsBefore -> $fontHostsAfter")
    exit 1
}

Write-Output ("ok    all {0} menu items returned in one process" -f $done)
exit 0
