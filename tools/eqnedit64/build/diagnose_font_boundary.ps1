[CmdletBinding()]
param(
    [Parameter(Mandatory)] [ValidateSet('self', 'model-idle', 'model-self')] [string]$Arm,
    [Parameter(Mandatory)] [string]$OutputDirectory
)
$ErrorActionPreference = 'Stop'
if ($env:EQNEDIT64_ISOLATED_TEST_SESSION -ne '1') { throw 'Disposable CI only.' }
$directory = [IO.Path]::GetFullPath($OutputDirectory)
if (Test-Path -LiteralPath $directory) { throw 'Use a fresh evidence directory.' }
New-Item -ItemType Directory -Path $directory | Out-Null
$env:EQNEDIT64_FONT_TRACE_DIR = $directory
$app = (Resolve-Path 'build_cmake/Release/Eqnedit64.exe').Path
$session = (Get-Process -Id $PID).SessionId
$started = Get-Date
$result = [ordered]@{
    arm = $Arm; status = 'INCONCLUSIVE'; reason = ''; commands = @()
    source_sha = $env:GITHUB_SHA; host = $env:COMPUTERNAME; session = $session
    app_sha256 = (Get-FileHash $app -Algorithm SHA256).Hash
    started_utc = $started.ToUniversalTime().ToString('o'); events = @()
}
function CrashEvents([datetime]$from) {
    try { $events = @(Get-WinEvent -FilterHashtable @{LogName='Application'; Id=1000; StartTime=$from} -ErrorAction Stop) }
    catch { if ($_.FullyQualifiedErrorId -like 'NoMatchingEventsFound*') { return @() }; throw }
    @($events | Where-Object Message -Match 'fontdrvhost\.exe' | ForEach-Object {
        [ordered]@{utc=$_.TimeCreated.ToUniversalTime().ToString('o'); record_id=$_.RecordId; message=$_.Message}
    })
}
function InvokeObserved([string]$file, [string]$argument) {
    $begin = (Get-Date).ToUniversalTime().ToString('o')
    $process = Start-Process -FilePath $file -ArgumentList $argument -WindowStyle Hidden -PassThru -Wait
    $result.commands += [ordered]@{argument=$argument; pid=$process.Id; exit_code=$process.ExitCode
        started_utc=$begin; ended_utc=(Get-Date).ToUniversalTime().ToString('o')}
    if ($process.ExitCode -ne 0) { throw "Command failed: $argument ($($process.ExitCode))" }
}
$job = $null
try {
    $result.control_events = @(CrashEvents ($started.AddMinutes(-10)))
    if ($result.control_events.Count) { throw 'Unhealthy pre-test control window.' }
    $job = Start-Job -ArgumentList $session, $directory -ScriptBlock {
        param($session, $directory)
        $ErrorActionPreference = 'Stop'
        $known = @{}
        $writer = [IO.StreamWriter]::new((Join-Path $directory 'font-host.jsonl'), $false)
        try {
            while (-not (Test-Path (Join-Path $directory 'stop'))) {
                $current = @(Get-Process -Name fontdrvhost -ErrorAction SilentlyContinue | Where-Object SessionId -EQ $session)
                $currentIds = @($current | Select-Object -ExpandProperty Id | Sort-Object)
                foreach ($process in $current) {
                    if (-not $known.ContainsKey($process.Id)) {
                        $queryError = $null
                        try { $null = $process.Handle } # Retain exit information if Windows permits it.
                        catch { $queryError = $_.Exception.Message }
                        $known[$process.Id] = @{process=$process; query_error=$queryError}
                    } else { $process.Dispose() }
                }
                $states = @(foreach ($entry in $known.Values) {
                    $process = $entry.process
                    $exited = $null; $exitUtc = $null
                    if (-not $entry.query_error) {
                        try {
                            $exited = $process.HasExited
                            if ($exited) { $exitUtc = $process.ExitTime.ToUniversalTime().ToString('o') }
                        } catch { $entry.query_error = $_.Exception.Message }
                    }
                    [ordered]@{pid=$process.Id; exited=$exited
                        exit_utc=$exitUtc; query_error=$entry.query_error}
                })
                $writer.WriteLine(([ordered]@{utc=(Get-Date).ToUniversalTime().ToString('o')
                    current_ids=$currentIds; processes=$states} | ConvertTo-Json -Compress -Depth 5))
                $writer.Flush()
                if (-not (Test-Path (Join-Path $directory 'ready'))) {
                    if (-not $current.Count) { throw 'No session font-host baseline.' }
                    [IO.File]::WriteAllText((Join-Path $directory 'ready'), 'ready')
                }
                Start-Sleep -Milliseconds 100
            }
        } finally {
            $writer.Dispose()
            foreach ($entry in $known.Values) { $entry.process.Dispose() }
        }
    }
    $deadline = (Get-Date).AddSeconds(20)
    while (-not (Test-Path (Join-Path $directory 'ready'))) {
        if ($job.State -eq 'Failed' -or (Get-Date) -gt $deadline) { throw 'Font-host monitor did not become ready.' }
        Start-Sleep -Milliseconds 100
    }
    if ($Arm -ne 'self') { InvokeObserved (Get-Command python).Source 'tests/run_model_tests.py' }
    if ($Arm -ne 'model-idle') { InvokeObserved $app '--self-test' }
    Start-Sleep -Seconds 3
    $result.events = @(CrashEvents $started)
    $result.status = 'OBSERVED'
} catch { $result.reason = $_.Exception.Message }
finally {
    try { $result.events = @(CrashEvents $started) }
    catch { $result.status = 'INCONCLUSIVE'; $result.reason += ' Event observation failed: ' + $_.Exception.Message }
    if ($job) {
        [IO.File]::WriteAllText((Join-Path $directory 'stop'), 'stop')
        $done = Wait-Job $job -Timeout 10
        if (-not $done -or $job.State -ne 'Completed') {
            $result.status = 'INCONCLUSIVE'; $result.reason += ' Monitor failed or timed out.'
            Stop-Job $job
        }
        Receive-Job $job -ErrorAction Continue | Out-Null
        Remove-Job $job
    }
    $timeline = Join-Path $directory 'font-host.jsonl'
    if (Test-Path $timeline) {
        $samples = @(Get-Content $timeline | ConvertFrom-Json)
        $result.samples = $samples.Count
        if (-not $samples.Count) { $result.status = 'INCONCLUSIVE'; $result.reason += ' Empty timeline.' }
        $result.font_host_exited = @($samples.processes | Where-Object exited -EQ $true).Count -gt 0
        $result.font_host_changed = @($samples | Where-Object {
            ($_.current_ids -join ',') -ne ($samples[0].current_ids -join ',')
        }).Count -gt 0
    }
    $result.ended_utc = (Get-Date).ToUniversalTime().ToString('o')
    $result | ConvertTo-Json -Depth 8 | Set-Content (Join-Path $directory 'result.json') -Encoding utf8
}
if ($result.status -ne 'OBSERVED') { exit 2 }
if ($result.events.Count -or $result.font_host_exited -or $result.font_host_changed) { exit 1 }
# A diagnostic non-reproduction is not product acceptance.
exit 0
