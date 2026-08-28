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

$started = Get-Date
$before = @(Get-SessionFontHostIds)
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
        throw ("fontdrvhost changed during Eqnedit64 iteration ${index}: " +
            "$($before -join ',') -> $($current -join ',')")
    }
}

# Windows Error Reporting can arrive just after the crashing process exits.
Start-Sleep -Seconds 2
$after = @(Get-SessionFontHostIds)
if ($before.Count -and (($after -join ',') -cne ($before -join ','))) {
    throw "fontdrvhost changed after Eqnedit64 stress: $($before -join ',') -> $($after -join ',')"
}
$crashes = @(
    Get-WinEvent -FilterHashtable @{
        LogName = 'Application'
        Id = 1000
        StartTime = $started
    } -ErrorAction SilentlyContinue |
        Where-Object Message -match 'fontdrvhost\.exe'
)
if ($crashes.Count) {
    throw "Eqnedit64 stress caused $($crashes.Count) fontdrvhost crash event(s)."
}

Write-Output ("PASS: {0} private-font lifecycles; fontdrvhost remained {1}" -f
    $Iterations, ($after -join ','))
