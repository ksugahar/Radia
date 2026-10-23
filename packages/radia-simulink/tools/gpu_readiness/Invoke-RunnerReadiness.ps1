# Manager-dispatched, read-only probe. No install, activation, build or process kill.
param(
    [Parameter(Mandatory)][string]$RunnerRoot,
    [Parameter(Mandatory)][string]$ExpectedSid,
    [Parameter(Mandatory)][string]$ExpectedRunnerName,
    [Parameter(Mandatory)][string]$MatlabRoot,
    [Parameter(Mandatory)][string]$ProbePython,
    [Parameter(Mandatory)][string]$ExpectedGPU,
    [Parameter(Mandatory)][string]$ProbeDir,
    [Parameter(Mandatory)][string]$EvidenceDir
)
$ErrorActionPreference = 'Stop'
function Fingerprint([string]$Value) {
    return [Convert]::ToHexString([Security.Cryptography.SHA256]::HashData([Text.Encoding]::UTF8.GetBytes($Value))).ToLowerInvariant()
}
$summary = @{ passed=$false; state='started'; schema='radia.runner-readiness.v1'; context=$null }
if (-not (Test-Path -LiteralPath $EvidenceDir -PathType Container)) { throw 'Manager must create the owned evidence directory first.' }
$summaryFile = Join-Path $EvidenceDir 'summary.json'
if (Test-Path -LiteralPath $summaryFile) { throw 'Fresh evidence directory required.' }
try {
    foreach ($value in @($RunnerRoot,$ExpectedSid,$ExpectedRunnerName,$MatlabRoot,$ProbePython,$ExpectedGPU)) {
        if ([string]::IsNullOrWhiteSpace($value)) { throw 'Every expected identity must be explicit.' }
    }
    if ($env:GITHUB_EVENT_NAME -ne 'workflow_dispatch') { throw 'Manual dispatch required.' }
    if ($env:GITHUB_ACTIONS -ne 'true' -or $env:RUNNER_NAME -ne $ExpectedRunnerName) { throw 'Actual approved runner job required.' }
    $chain = @(); $cursor = $PID
    for ($index=0; $index -lt 64 -and $cursor -gt 0; $index++) {
        $process = Get-CimInstance Win32_Process -Filter "ProcessId = $cursor"
        if (-not $process) { break }
        $chain += $process
        if ($process.ParentProcessId -eq $cursor) { break }
        $cursor = $process.ParentProcessId
    }
    $workerPath = [IO.Path]::GetFullPath((Join-Path $RunnerRoot 'bin/Runner.Worker.exe'))
    $listenerPath = [IO.Path]::GetFullPath((Join-Path $RunnerRoot 'bin/Runner.Listener.exe'))
    $worker = @($chain | Where-Object { $_.ExecutablePath -eq $workerPath })
    $listener = @($chain | Where-Object { $_.ExecutablePath -eq $listenerPath })
    if ($worker.Count -ne 1 -or $listener.Count -ne 1) { throw 'Expected Worker and Listener must be ancestors, not just running elsewhere.' }
    $current = $chain[0]; $identities = @()
    if ($current.SessionId -eq 0) { throw 'Interactive authenticated account required; service session refused.' }
    foreach ($process in @($current,$worker[0],$listener[0])) {
        $sid = Invoke-CimMethod -InputObject $process -MethodName GetOwnerSid
        if ($sid.ReturnValue -ne 0 -or $sid.Sid -ne $ExpectedSid -or $process.SessionId -ne $current.SessionId) { throw 'Runner account/session mismatch.' }
        $identities += @{ pid=$process.ProcessId; sid_sha256=(Fingerprint $sid.Sid); session=$process.SessionId; executable_sha256=(Fingerprint $process.ExecutablePath) }
    }
    $summary.context = @{ identities=$identities; runner_sha256=(Fingerprint $env:RUNNER_NAME); run_id=$env:GITHUB_RUN_ID; attempt=$env:GITHUB_RUN_ATTEMPT; source=$env:GITHUB_SHA }
    # Manager must have admitted an idle GPU window before this job; don't kill jobs.
    $compute = & nvidia-smi.exe --query-compute-apps=pid --format=csv,noheader,nounits
    if ($LASTEXITCODE -ne 0) { throw 'GPU-idle query failed.' }
    if (@($compute | Where-Object { $_.Trim() }).Count -gt 0) { throw 'GPU compute processes exist or idle status is ambiguous; manager must resolve.' }
    if (-not (Test-Path -LiteralPath (Join-Path $MatlabRoot 'bin/matlab.exe') -PathType Leaf) -or -not (Test-Path -LiteralPath $ProbePython -PathType Leaf)) { throw 'Existing MATLAB and Python executables required; nothing will be installed.' }
    $vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio/Installer/vswhere.exe'
    if (-not (Test-Path -LiteralPath $vswhere -PathType Leaf)) { throw 'VS discovery missing.' }
    $vs = & $vswhere -latest -products '*' -version '[17.0,18.0)' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
    if ($LASTEXITCODE -ne 0) { throw 'VS2022 discovery failed.' }
    if (@($vs).Count -ne 1 -or -not $vs) { throw 'VS2022 C++ instance missing.' }
    $msvc = (Get-Content -LiteralPath (Join-Path $vs 'VC/Auxiliary/Build/Microsoft.VCToolsVersion.default.txt') -Raw).Trim()
    if ($msvc -ne '14.44.35207') { throw 'Approved MSVC 14.44.35207 required; VS2026 is excluded.' }
    $summary.toolchain = @{ vs_major=17; msvc=$msvc; root_sha256=(Fingerprint $vs); scope='discovery only, no compiler or mexcuda execution' }
    $env:OMP_NUM_THREADS='2'; $env:MKL_NUM_THREADS='2'; $env:OPENBLAS_NUM_THREADS='2'
    $env:PYTHONUTF8='1'; $env:PYTHONIOENCODING='utf-8'; $env:PYTHONDONTWRITEBYTECODE='1'
    $batchLog = Join-Path $EvidenceDir 'matlab-batch.log'
    $engineFile = Join-Path $EvidenceDir 'engine.json'
    if ((Test-Path -LiteralPath $batchLog) -or (Test-Path -LiteralPath $engineFile)) { throw 'Fresh outputs required.' }
    $escapedDir=$ProbeDir.Replace("'","''"); $escapedRoot=$MatlabRoot.Replace("'","''"); $escapedGPU=$ExpectedGPU.Replace("'","''")
    $expression="addpath('$escapedDir'); r=probe_matlab_readiness('$escapedRoot','$escapedGPU'); value=jsondecode(r); assert(value.passed,'Readiness failed');"
    & (Join-Path $MatlabRoot 'bin/matlab.exe') -wait -batch $expression *> $batchLog
    if ($LASTEXITCODE -ne 0) { throw 'MATLAB batch readiness failed; inspect log.' }
    $lines=@(Get-Content -LiteralPath $batchLog | Where-Object { $_.StartsWith('RADIA_READINESS:') })
    if ($lines.Count -ne 1) { throw 'Expected one MATLAB readiness record.' }
    $batch=$lines[0].Substring('RADIA_READINESS:'.Length) | ConvertFrom-Json
    if ($batch.passed -ne $true) { throw 'Batch record did not pass.' }
    & $ProbePython (Join-Path $ProbeDir 'probe_engine.py') --matlab-root $MatlabRoot --gpu $ExpectedGPU --probe-dir $ProbeDir --result $engineFile *> (Join-Path $EvidenceDir 'engine.log')
    if ($LASTEXITCODE -ne 0) { throw 'Python Engine readiness failed; inspect record.' }
    $engine=Get-Content -LiteralPath $engineFile -Raw | ConvertFrom-Json
    if ($engine.passed -ne $true -or $engine.matlab.gpu.uuid -ne $batch.gpu.uuid -or $engine.matlab.matlabroot -ne $batch.matlabroot) { throw 'Engine/batch identity mismatch.' }
    # Upload only portable readiness facts; raw logs/paths/account identity stay in RUNNER_TEMP.
    $summary.batch = @{ passed=$batch.passed; release=$batch.release; version=$batch.matlab_version; licenses=@{}; gpu=@{ name=$batch.gpu.name; uuid_sha256=(Fingerprint $batch.gpu.uuid); compute_capability=$batch.gpu.compute_capability; toolkit=$batch.gpu.toolkit; total_memory_bytes=$batch.gpu.total_memory_bytes } }
    foreach ($feature in @('MATLAB','Simulink','Distrib_Computing_Toolbox')) {
        $summary.batch.licenses[$feature] = @{ exists=$batch.licenses.$feature.exists; checkout=$batch.licenses.$feature.checkout }
    }
    $summary.engine = @{ passed=$engine.passed; engine_version=$engine.engine_version; python_version=$engine.python_version.Split([Environment]::NewLine)[0]; shutdown=$engine.shutdown; same_installation=$true; same_gpu=$true }
    if ($engine.shutdown -ne 'complete') { throw 'Owned Engine did not report shutdown.' }
    $summary.passed=$true
} catch {
    $summary.error_type=$_.Exception.GetType().Name; $summary.error_code='readiness_failed'; Write-Host 'Readiness failed; sanitized evidence records failure.'
} finally {
    $summary.state='finished'
    $summary | ConvertTo-Json -Depth 20 | Out-File -LiteralPath $summaryFile -Encoding utf8 -NoClobber
}

if (-not $summary.passed) { exit 1 }
