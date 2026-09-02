[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$AppPath,
    [ValidateRange(2, 256)]
    [int]$Iterations = 32
)

$ErrorActionPreference = 'Stop'
if ($env:EQNEDIT64_ISOLATED_TEST_SESSION -ne '1') {
    throw @'
The font-lifecycle stress test repeatedly starts a private-font process.
Run it only in disposable CI, a VM, or a dedicated Windows user session with
EQNEDIT64_ISOLATED_TEST_SESSION=1. Never run it on the interactive lab desktop.
'@
}

$app = (Resolve-Path -LiteralPath $AppPath).Path
$sessionId = (Get-Process -Id $PID).SessionId
function Get-SessionFontHostIds {
    return @(
        Get-Process -Name fontdrvhost -ErrorAction SilentlyContinue |
            Where-Object SessionId -eq $sessionId |
            Sort-Object Id |
            Select-Object -ExpandProperty Id
    )
}
function Get-FontHostCrashes([datetime]$From, [datetime]$Until) {
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
$before = @(Get-SessionFontHostIds)
$changedAt = $null
for ($index = 0; $index -lt $Iterations; $index++) {
    $process = Start-Process -FilePath $app `
        -ArgumentList '--status-layout-test' -WindowStyle Hidden `
        -Wait -PassThru
    if ($process.ExitCode -ne 0) {
        throw "Eqnedit64 font lifecycle iteration $index failed: $($process.ExitCode)"
    }
    $current = @(Get-SessionFontHostIds)
    if ($before.Count -and
        (($current -join ',') -cne ($before -join ','))) {
        $changedAt = $index
        break
    }
}

# Windows Error Reporting can arrive just after the crashing process exits.
Start-Sleep -Seconds 2
$after = @(Get-SessionFontHostIds)
$ended = Get-Date
$duration = $ended - $started
$crashes = @(Get-FontHostCrashes $started $ended)
$controlDuration = if ($duration.TotalMinutes -lt 10) {
    [timespan]::FromMinutes(10)
} else {
    $duration
}
$controlCrashes = @(
    Get-FontHostCrashes ($started - $controlDuration) $started)
$changed = $before.Count -and
    (($after -join ',') -cne ($before -join ','))
if ($changed -or $crashes.Count) {
    $observation = (
        "fontdrvhost $($before -join ',') -> $($after -join ','); " +
        "changed at iteration=$changedAt; active crash events=$($crashes.Count); " +
        "pre-test control events=$($controlCrashes.Count) " +
        "over $([math]::Round($controlDuration.TotalSeconds, 1))s")
    if ($controlCrashes.Count) {
        Write-Output (
            "INCONCLUSIVE: the supposedly isolated session was already " +
            "losing fontdrvhost before this lifecycle probe: $observation")
        exit 2
    }
    if ($changed -and -not $crashes.Count) {
        Write-Output (
            "INCONCLUSIVE: fontdrvhost changed without a matching Application " +
            "Error event: $observation")
        exit 2
    }
    throw "FAIL: fontdrvhost crashed during the isolated Eqnedit64 lifecycle probe: $observation"
}

Write-Output ("PASS: {0} private-font lifecycles; fontdrvhost remained {1}" -f
    $Iterations, ($after -join ','))
