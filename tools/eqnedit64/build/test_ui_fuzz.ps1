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
function Get-FontDriverCrashes([datetime]$From, [datetime]$Until) {
    return @(
        Get-WinEvent -FilterHashtable @{
            LogName = 'Application'
            Id = 1000
            StartTime = $From
        } -ErrorAction SilentlyContinue |
            Where-Object {
                $_.TimeCreated -le $Until -and
                $_.Message -match 'fontdrvhost\.exe'
            }
    )
}

$started = Get-Date
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
$ended = Get-Date
$last = $progress | Select-Object -Last 1
Remove-Item -LiteralPath $log -ErrorAction SilentlyContinue
if ($process.ExitCode -ne 0) {
    throw ("GUI fuzz exit 0x{0:X8}; last progress: {1}" -f `
           $process.ExitCode, $last)
}
if ($fontHostsBefore -and $fontHostsAfter -ne $fontHostsBefore) {
    # Windows Error Reporting can publish Event 1000 shortly after the host
    # replacement becomes visible.
    Start-Sleep -Seconds 2
    $ended = Get-Date
    $fontHostsAfter = Get-SessionFontDriverHostIds
    $duration = $ended - $started
    $activeCrashes = @(Get-FontDriverCrashes $started $ended)
    $controlDuration = if ($duration.TotalMinutes -lt 10) {
        [timespan]::FromMinutes(10)
    } else {
        $duration
    }
    $controlCrashes = @(
        Get-FontDriverCrashes ($started - $controlDuration) $started)
    $observation = (
        "fontdrvhost $fontHostsBefore -> $fontHostsAfter; " +
        "active crash events=$($activeCrashes.Count), " +
        "pre-test control events=$($controlCrashes.Count) " +
        "over $([math]::Round($controlDuration.TotalSeconds, 1))s")

    if ($env:EQNEDIT64_ISOLATED_TEST_SESSION -ne '1') {
        Write-Warning (
            "INCONCLUSIVE font-host observation on a shared interactive " +
            "session: $observation. Structural fuzz itself passed; do not " +
            "attribute a system-wide host replacement to Eqnedit64 here.")
    } elseif ($controlCrashes.Count) {
        Write-Output (
            "INCONCLUSIVE: the isolated session was already losing its font " +
            "host before Eqnedit64 fuzz: $observation")
        exit 2
    } elseif (-not $activeCrashes.Count) {
        Write-Output (
            "INCONCLUSIVE: fontdrvhost changed without a matching Application " +
            "Error event: $observation")
        exit 2
    } else {
        throw "FAIL: fontdrvhost crashed during isolated structural fuzz: $observation"
    }
}

Write-Output ("ok    {0} seeds x {1} operations through one real-window process" -f `
              $Seeds, $Ops)
