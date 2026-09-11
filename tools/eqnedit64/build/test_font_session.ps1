[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)] [string]$AppPath,
    [ValidateRange(2, 256)] [int]$Iterations = 32,
    [string]$OutputPath = 'C:\temp\Eqnedit64-font-session.json'
)
$ErrorActionPreference = 'Stop'
if ($env:EQNEDIT64_ISOLATED_TEST_SESSION -ne '1') {
    throw 'Run only in disposable CI or an isolated VM. Never run on the interactive lab desktop.'
}
function Get-SessionFontHostIds {
    return @(Get-Process -Name fontdrvhost -ErrorAction SilentlyContinue |
        Where-Object SessionId -eq $sessionId | Sort-Object Id |
        Select-Object -ExpandProperty Id)
}
function Get-FontHostCrashes([datetime]$From, [datetime]$Until) {
    try {
        $events = @(Get-WinEvent -FilterHashtable @{
            LogName = 'Application'; Id = 1000; StartTime = $From
        } -ErrorAction Stop)
    } catch {
        # An empty query is normal; access/provider failures are not evidence.
        if ($_.FullyQualifiedErrorId -like 'NoMatchingEventsFound*') { return @() }
        throw
    }
    return @($events | Where-Object {
        $_.TimeCreated -le $Until -and $_.Message -match 'fontdrvhost\.exe'
    } | ForEach-Object {
        [ordered]@{ time = $_.TimeCreated.ToUniversalTime().ToString('o')
            record_id = $_.RecordId; message = $_.Message }
    })
}
$result = [ordered]@{
    schema = 'eqnedit64.font-session.v1'
    status = 'INCONCLUSIVE'; reason = ''; host = $env:COMPUTERNAME
    os_version = [Environment]::OSVersion.VersionString
    started_utc = (Get-Date).ToUniversalTime().ToString('o')
    app_path = $AppPath; app_sha256 = $null; session_id = $null
    before = @(); after = @(); control_events = @(); active_events = @()
    iterations_requested = $Iterations; iterations = @()
}
$exitCode = 2
try {
    $app = (Resolve-Path -LiteralPath $AppPath).Path
    $result.app_path = $app
    $result.app_sha256 = (Get-FileHash -LiteralPath $app -Algorithm SHA256).Hash
    $sessionId = (Get-Process -Id $PID).SessionId
    $result.session_id = $sessionId
    $started = Get-Date
    $before = @(Get-SessionFontHostIds)
    $result.before = $before
    $result.control_events = @(Get-FontHostCrashes ($started - [timespan]::FromMinutes(10)) $started)
    if (-not $before.Count) {
        $result.reason = 'No fontdrvhost baseline in the test session.'
    } elseif ($result.control_events.Count) {
        $result.reason = 'Font-host crashes in the pre-test control window; no executable was started.'
    } else {
        $iterationFailed = $false
        for ($index = 0; $index -lt $Iterations; $index++) {
            $process = Start-Process -FilePath $app -ArgumentList '--status-layout-test' `
                -WindowStyle Hidden -Wait -PassThru
            $current = @(Get-SessionFontHostIds)
            $result.iterations += [ordered]@{
                index = $index; exit_code = $process.ExitCode
                font_host_ids = $current; observed_utc = (Get-Date).ToUniversalTime().ToString('o')
            }
            if ($process.ExitCode -ne 0) { $iterationFailed = $true }
            if ($iterationFailed -or (($current -join ',') -cne ($before -join ','))) { break }
        }
        # Windows Error Reporting can arrive just after process exit.
        Start-Sleep -Seconds 2
        $after = @(Get-SessionFontHostIds)
        $result.after = $after
        $result.active_events = @(Get-FontHostCrashes $started (Get-Date))
        $changed = (($after -join ',') -cne ($before -join ',')) -or
            @($result.iterations | Where-Object {
                ($_.font_host_ids -join ',') -cne ($before -join ',')
            }).Count -gt 0
        if ($result.active_events.Count -or $iterationFailed) {
            $result.status = 'FAIL'; $exitCode = 1
            $result.reason = 'Executable failure or font-host crash during lifecycle probe.'
        } elseif ($changed) {
            $result.reason = 'Font-host identity changed without a matching crash event.'
        } else {
            $result.status = 'PASS'; $exitCode = 0
            $result.reason = 'All private-font lifecycles passed with stable font-host identity.'
        }
    }
} catch {
    $result.reason = 'Observation failed: ' + $_.Exception.Message
} finally {
    $result.ended_utc = (Get-Date).ToUniversalTime().ToString('o')
    $destination = [IO.Path]::GetFullPath($OutputPath)
    [IO.Directory]::CreateDirectory([IO.Path]::GetDirectoryName($destination)) | Out-Null
    $result | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $destination -Encoding utf8
}
Write-Output ("{0}: {1} Evidence: {2}" -f $result.status, $result.reason, $destination)
exit $exitCode
