[CmdletBinding()]
param(
    [Parameter(Mandatory)] [ValidateSet('self', 'model-idle', 'model-self', 'probe')] [string]$Arm,
    [Parameter(Mandatory)] [string]$OutputDirectory,
    [ValidateSet('exit', 'remove')] [string]$ReleaseMode = 'exit',
    [ValidateSet('none', 'measure')] [string]$Measurement = 'none',
    [ValidateRange(1, 64)] [int]$MaxLaunches = 64,
    [ValidateRange(-1, 10)] [int]$ModelPrefix = -1
)
$ErrorActionPreference = 'Stop'
if ($env:EQNEDIT64_ISOLATED_TEST_SESSION -ne '1') { throw 'Disposable CI only.' }
$directory = [IO.Path]::GetFullPath($OutputDirectory)
if (Test-Path -LiteralPath $directory) { throw 'Use a fresh evidence directory.' }
New-Item -ItemType Directory -Path $directory | Out-Null
$env:EQNEDIT64_FONT_TRACE_DIR = $directory
$app = (Resolve-Path $(if ($Arm -eq 'probe') { 'build_cmake/Release/font_lifecycle_probe.exe' } else { 'build_cmake/Release/Eqnedit64.exe' })).Path
$session = (Get-Process -Id $PID).SessionId
$started = Get-Date
$result = [ordered]@{
    arm = $Arm; status = 'INCONCLUSIVE'; reason = ''; commands = @()
    source_sha = $env:GITHUB_SHA; host = $env:COMPUTERNAME; session = $session
    payload_source_sha = $env:EQNEDIT64_PAYLOAD_SOURCE_SHA; model_prefix = $ModelPrefix
    os_version = [Environment]::OSVersion.VersionString
    app_sha256 = (Get-FileHash $app -Algorithm SHA256).Hash
    started_utc = $started.ToUniversalTime().ToString('o'); events = @()
}
if ($Arm -eq 'probe') {
    $font = (Resolve-Path 'assets/latinmodern-math.otf').Path
    $result.font_sha256 = (Get-FileHash $font -Algorithm SHA256).Hash
    $result.release_mode = $ReleaseMode; $result.measurement = $Measurement
    $result.max_launches = $MaxLaunches; $result.baseline_idle_s = 3; $result.tail_idle_s = 3
    $result.first_failure_k = $null
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
    $index = $result.commands.Count
    $process = Start-Process -FilePath $file -ArgumentList $argument -WindowStyle Hidden -PassThru -Wait `
        -RedirectStandardOutput (Join-Path $directory "command-$index.stdout") `
        -RedirectStandardError (Join-Path $directory "command-$index.stderr")
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
                        $startUtc = $null
                        try {
                            $null = $process.Handle # Retain exit information if Windows permits it.
                            $startUtc = $process.StartTime.ToUniversalTime().ToString('o')
                        }
                        catch { $queryError = $_.Exception.Message }
                        $known[$process.Id] = @{process=$process; query_error=$queryError; start_utc=$startUtc}
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
                        start_utc=$entry.start_utc; exit_utc=$exitUtc; query_error=$entry.query_error}
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
    if ($Arm -eq 'probe') {
        $baseline = @((Get-Content (Join-Path $directory 'font-host.jsonl') -TotalCount 1 | ConvertFrom-Json).current_ids)
        Start-Sleep -Seconds 3
        $current = @(Get-Process -Name fontdrvhost -ErrorAction SilentlyContinue | Where-Object SessionId -EQ $session | Select-Object -ExpandProperty Id | Sort-Object)
        if (($current -join ',') -ne ($baseline -join ',')) { throw 'Font-host changed during baseline idle.' }
        $result.status = 'NO_REPRODUCTION_WITHIN_BOUND'
        for ($k = 1; $k -le $MaxLaunches; $k++) {
            if ($job.State -ne 'Running') { throw 'Font-host monitor stopped during the trial.' }
            $current = @(Get-Process -Name fontdrvhost -ErrorAction SilentlyContinue | Where-Object SessionId -EQ $session | Select-Object -ExpandProperty Id | Sort-Object)
            if (($current -join ',') -ne ($baseline -join ',')) {
                $result.status = 'HOST_CHANGED'; $result.first_failure_k = $k - 1
                $result.reason = 'Host change observed before the next launch.'
                break
            }
            $begin = (Get-Date).ToUniversalTime().ToString('o')
            $process = Start-Process -FilePath $app -ArgumentList ('"' + $font + '" ' + $ReleaseMode + ' ' + $Measurement + ' 1') `
                -WindowStyle Hidden -Wait -PassThru `
                -RedirectStandardOutput (Join-Path $directory ('probe-{0:d3}.jsonl' -f $k)) `
                -RedirectStandardError (Join-Path $directory ('probe-{0:d3}.stderr' -f $k))
            $current = @(Get-Process -Name fontdrvhost -ErrorAction SilentlyContinue | Where-Object SessionId -EQ $session | Select-Object -ExpandProperty Id | Sort-Object)
            $result.commands += [ordered]@{k=$k; pid=$process.Id; exit_code=$process.ExitCode
                started_utc=$begin; ended_utc=(Get-Date).ToUniversalTime().ToString('o'); font_host_ids=$current}
            $changed = ($current -join ',') -ne ($baseline -join ',')
            if ($process.ExitCode -ne 0 -or $changed) {
                $result.status = $(if ($changed) {'HOST_CHANGED'} else {'CHILD_FAILURE'})
                $result.first_failure_k = $k
                break
            }
        }
    } else {
        if ($Arm -ne 'self') {
            $arguments = 'tests/run_model_tests.py'
            if ($ModelPrefix -ge 0) { $arguments += " --diagnostic-prefix $ModelPrefix" }
            InvokeObserved (Get-Command python).Source $arguments
        }
        if ($Arm -ne 'model-idle') { InvokeObserved $app '--self-test' }
        $result.status = 'OBSERVED'
    }
    Start-Sleep -Seconds 3
    $result.events = @(CrashEvents $started)
} catch { $result.status = 'INCONCLUSIVE'; $result.reason = $_.Exception.Message }
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
        $result.query_errors = @($samples.processes | Where-Object query_error | Select-Object -ExpandProperty query_error -Unique)
        if ($result.query_errors.Count) { $result.status = 'INCONCLUSIVE'; $result.reason += ' Process observation errors.' }
        if (-not $samples.Count) { $result.status = 'INCONCLUSIVE'; $result.reason += ' Empty timeline.' }
        $result.font_host_exited = @($samples.processes | Where-Object exited -EQ $true).Count -gt 0
        $result.font_host_changed = @($samples | Where-Object {
            ($_.current_ids -join ',') -ne ($samples[0].current_ids -join ',')
        }).Count -gt 0
    }
    if ($Arm -eq 'probe' -and $result.status -eq 'NO_REPRODUCTION_WITHIN_BOUND' -and
        ($result.events.Count -or $result.font_host_exited -or $result.font_host_changed)) {
        $result.status = 'HOST_FAILURE_IN_OBSERVATION_WINDOW'
        $result.reason = 'Event or timeline detects host failure; exact causative launch is not inferred.'
    }
    $result.ended_utc = (Get-Date).ToUniversalTime().ToString('o')
    $result | ConvertTo-Json -Depth 8 | Set-Content (Join-Path $directory 'result.json') -Encoding utf8
}
if ($result.status -eq 'INCONCLUSIVE') { exit 2 }
if ($result.status -in @('HOST_CHANGED', 'CHILD_FAILURE', 'HOST_FAILURE_IN_OBSERVATION_WINDOW')) { exit 1 }
if ($result.events.Count -or $result.font_host_exited -or $result.font_host_changed) { exit 1 }
# A diagnostic non-reproduction is not product acceptance.
exit 0
