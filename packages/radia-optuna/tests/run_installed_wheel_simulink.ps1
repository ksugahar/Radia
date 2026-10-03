param(
    [Parameter(Mandatory = $true)]
    [string]$Wheel,
    [string]$MatlabExecutable = 'matlab',
    [string]$PythonExecutable = 'python',
    [string]$EvidenceOutput = '',
    [string]$PreverifiedWheelSha256 = '',
    # An explicitly shared MATLAB Engine to reuse (never quit). Without it a
    # dedicated session is started only when no MATLAB process or shared
    # Engine exists; an existing MATLAB without a selected session refuses.
    [string]$EngineSession = '',
    # The solver release Engine worker that owns this policy.
    [string]$EngineWorker = ''
)

$ErrorActionPreference = 'Stop'
$wheelPath = (Resolve-Path -LiteralPath $Wheel).ProviderPath
$testDirectory = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot 'matlab')).ProviderPath
$packageRoot = Split-Path -Parent $PSScriptRoot
$wheelVerifier = Join-Path $packageRoot 'verify_wheel.py'
$runRoot = Join-Path 'C:\temp' ('radia-optuna-wheel-simulink-' + [guid]::NewGuid().ToString('N'))
$resolvedTempRoot = [IO.Path]::GetFullPath('C:\temp') + [IO.Path]::DirectorySeparatorChar
$resolvedRunRoot = [IO.Path]::GetFullPath($runRoot)
if (-not $resolvedRunRoot.StartsWith($resolvedTempRoot, [StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to create the isolated test outside C:\temp: $resolvedRunRoot"
}

function ConvertTo-MatlabLiteral([string]$Value) {
    return $Value.Replace("'", "''")
}

New-Item -ItemType Directory -Path $resolvedRunRoot | Out-Null
try {
    if ($PreverifiedWheelSha256) {
        if ($PreverifiedWheelSha256 -notmatch '^[0-9a-fA-F]{64}$') {
            throw 'PreverifiedWheelSha256 must contain exactly 64 hexadecimal characters'
        }
        $actualWheelSha256 = (Get-FileHash -LiteralPath $wheelPath -Algorithm SHA256).Hash.ToLowerInvariant()
        if ($actualWheelSha256 -ne $PreverifiedWheelSha256.ToLowerInvariant()) {
            throw "Preverified wheel SHA256 mismatch: expected $PreverifiedWheelSha256, got $actualWheelSha256"
        }
        $wheelVerification = [ordered]@{
            schema = 'radia-optuna.preverified-wheel.v1'
            ok = $true
            verification_mode = 'release-quad-exact-sha256'
            sha256 = $actualWheelSha256
        }
    } else {
        if (-not (Test-Path -LiteralPath $wheelVerifier -PathType Leaf)) {
            throw "Repository wheel verifier is missing: $wheelVerifier"
        }
        $wheelVerificationJson = (& $PythonExecutable $wheelVerifier $wheelPath --json | Out-String).Trim()
        if ($LASTEXITCODE -ne 0) {
            throw "Wheel verification failed with exit code $LASTEXITCODE"
        }
        $wheelVerification = $wheelVerificationJson | ConvertFrom-Json
    }

    $venv = Join-Path $resolvedRunRoot 'venv'
    & $PythonExecutable -m venv $venv
    if ($LASTEXITCODE -ne 0) { throw "Python venv creation failed with exit code $LASTEXITCODE" }

    $venvPython = Join-Path $venv 'Scripts\python.exe'
    & $venvPython -m pip install --disable-pip-version-check --no-deps $wheelPath
    if ($LASTEXITCODE -ne 0) { throw "Installing the wheel failed with exit code $LASTEXITCODE" }

    $doctor = Join-Path $venv 'Scripts\radia-optuna-doctor.exe'
    $doctorJson = (& $doctor --json | Out-String).Trim()
    if ($LASTEXITCODE -ne 0) { throw "radia-optuna-doctor failed with exit code $LASTEXITCODE" }
    $doctorEvidence = $doctorJson | ConvertFrom-Json

    $installedMatlabPath = (& $venvPython -c 'from radia_optuna import matlab_path; print(matlab_path())').Trim()
    if (-not (Test-Path -LiteralPath $installedMatlabPath -PathType Container)) {
        throw "Installed MATLAB directory is missing: $installedMatlabPath"
    }

    $matlabPathLiteral = ConvertTo-MatlabLiteral $installedMatlabPath
    $testDirectoryLiteral = ConvertTo-MatlabLiteral $testDirectory
    $simulinkEvidencePath = Join-Path $resolvedRunRoot 'simulink-evidence.json'
    $simulinkEvidenceLiteral = ConvertTo-MatlabLiteral $simulinkEvidencePath
    # A dedicated session starts from the default path. A borrowed session
    # keeps its caller's path (the worker restores it afterwards) because
    # restoredefaultpath would cut that session's own MCP/tool functions;
    # the candidate directories are prepended so they take precedence. The
    # acceptance runs in a function workspace, leaving base variables alone.
    $pathReset = if ($EngineSession) { '' } else { 'restoredefaultpath; ' }
    $batch = "${pathReset}addpath('$matlabPathLiteral','-begin'); addpath('$testDirectoryLiteral','-begin'); run_installed_wheel_acceptance('$matlabPathLiteral','$simulinkEvidenceLiteral');"

    # MATLAB runs through the solver release Engine worker, which reuses only
    # the selected shared session (PID-checked, with path, folder and
    # environment restored and no loaded diagrams or Radia/Optuna MEX), or
    # owns one new session when no MATLAB process or shared Engine exists.
    if (-not $EngineWorker) {
        $EngineWorker = Join-Path $PSScriptRoot '..\..\..\tools\verify_simulink_release.py'
    }
    if (-not (Test-Path -LiteralPath $EngineWorker -PathType Leaf)) {
        throw "MATLAB Engine worker is missing: $EngineWorker"
    }
    $matlabCommand = (Get-Command -Name $MatlabExecutable -CommandType Application -ErrorAction Stop |
        Select-Object -First 1).Source
    $matlabRoot = Split-Path -Parent (Split-Path -Parent $matlabCommand)

    # The worker runs in the isolated venv with MathWorks' matlabengine
    # release matched to this MATLAB, never a system interpreter's Engine.
    # Its in-tree setup.py derives the MATLAB root from the build folder, so
    # it is neither copied nor built inside the MATLAB installation.
    [xml]$versionInfo = Get-Content -LiteralPath (Join-Path $matlabRoot 'VersionInfo.xml') -Raw
    $matlabVersion = [string]$versionInfo.MathWorks_version_info.version
    if ($matlabVersion -notmatch '^(\d+)\.(\d+)\.') {
        throw "Cannot read the MATLAB version from $matlabRoot"
    }
    $engineRequirement = "matlabengine==$($Matches[1]).$($Matches[2]).*"
    & $venvPython -m pip install --disable-pip-version-check $engineRequirement
    if ($LASTEXITCODE -ne 0) { throw "Installing $engineRequirement failed with exit code $LASTEXITCODE" }
    # The Engine records the MATLAB it binds to; it must be this MATLAB.
    $engineProbe = (& $venvPython -c 'import importlib.metadata as m, json, os, matlab.engine; lines = open(os.path.join(os.path.dirname(matlab.engine.__file__), "_arch.txt")).read().splitlines(); print(json.dumps({"version": m.version("matlabengine"), "bin": lines[1]}))' | Out-String).Trim()
    if ($LASTEXITCODE -ne 0 -or -not $engineProbe) {
        throw "The isolated venv cannot import matlab.engine"
    }
    $engineBinding = $engineProbe | ConvertFrom-Json
    $expectedBin = [IO.Path]::GetFullPath((Join-Path $matlabRoot 'bin')).TrimEnd('\') + '\'
    if (-not ([IO.Path]::GetFullPath($engineBinding.bin) + '\').StartsWith($expectedBin, [StringComparison]::OrdinalIgnoreCase)) {
        throw "matlabengine binds to $($engineBinding.bin), not $matlabRoot"
    }
    $engineEvidence = [ordered]@{
        distribution = 'matlabengine'
        requirement = $engineRequirement
        version = $engineBinding.version
        matlab_version = $matlabVersion
        matlab_bin = $engineBinding.bin
        interpreter = $venvPython
        session = if ($EngineSession) { $EngineSession } else { $null }
        ownership = if ($EngineSession) { 'borrowed' } else { 'owned' }
    }
    $workerArguments = @('-X', 'utf8', '-s', $EngineWorker, '--engine-worker', $matlabRoot, $batch)
    if ($EngineSession) { $workerArguments += $EngineSession }

    # MathWorks' online license service can transiently reject an otherwise
    # valid start with error 5202. Retry only that failure, and only for an
    # owned session; a borrowed session, a MATLAB assertion, numerical
    # mismatch, or any other error still fails on the first attempt.
    $matlabLog = Join-Path $resolvedRunRoot 'matlab-engine.log'
    $maxMatlabAttempts = if ($EngineSession) { 1 } else { 3 }
    for ($attempt = 1; $attempt -le $maxMatlabAttempts; $attempt++) {
        & $venvPython @workerArguments 2>&1 | Tee-Object -FilePath $matlabLog
        $matlabExitCode = $LASTEXITCODE
        if ($matlabExitCode -eq 0) { break }

        $matlabOutput = Get-Content -LiteralPath $matlabLog -Raw -ErrorAction SilentlyContinue
        $isLicenseService5202 = $matlabOutput -match '(?<!\d)5202(?!\d)'
        if (-not $isLicenseService5202 -or $attempt -eq $maxMatlabAttempts) {
            throw "Installed-wheel Simulink test failed with exit code $matlabExitCode"
        }

        $delaySeconds = 10 * $attempt
        Write-Warning "MathWorks license service returned 5202; retrying the same MATLAB Engine command in $delaySeconds seconds (attempt $($attempt + 1)/$maxMatlabAttempts)."
        Start-Sleep -Seconds $delaySeconds
    }
    $simulinkEvidence = Get-Content -LiteralPath $simulinkEvidencePath -Raw |
        ConvertFrom-Json
    $evidence = [ordered]@{
        schema = 'radia-optuna.installed-wheel-evidence.v1'
        ok = $true
        wheel_verification = $wheelVerification
        doctor = $doctorEvidence
        matlab_engine = $engineEvidence
        simulink_e2e = $simulinkEvidence
        table_resume = $simulinkEvidence.table_resume
    }
    $encodedEvidence = $evidence | ConvertTo-Json -Depth 12
    if ($EvidenceOutput) {
        $resolvedEvidenceOutput = [IO.Path]::GetFullPath($EvidenceOutput)
        $evidenceParent = Split-Path -Parent $resolvedEvidenceOutput
        if ($evidenceParent -and -not (Test-Path -LiteralPath $evidenceParent)) {
            New-Item -ItemType Directory -Path $evidenceParent -Force | Out-Null
        }
        [IO.File]::WriteAllText(
            $resolvedEvidenceOutput,
            $encodedEvidence + [Environment]::NewLine,
            [Text.UTF8Encoding]::new($false))
    }
    Write-Output $encodedEvidence
    Write-Output 'RADIA_OPTUNA_WHEEL_SIMULINK_OK'
} finally {
    $checkedRunRoot = [IO.Path]::GetFullPath($resolvedRunRoot)
    if ($checkedRunRoot.StartsWith($resolvedTempRoot, [StringComparison]::OrdinalIgnoreCase) -and
            (Test-Path -LiteralPath $checkedRunRoot)) {
        Remove-Item -LiteralPath $checkedRunRoot -Recurse -Force
    }
}
