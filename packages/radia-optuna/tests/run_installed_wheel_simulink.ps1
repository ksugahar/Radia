param(
    [Parameter(Mandatory = $true)]
    [string]$Wheel,
    [string]$MatlabExecutable = 'matlab',
    [string]$PythonExecutable = 'python',
    [string]$EvidenceOutput = '',
    [string]$PreverifiedWheelSha256 = ''
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
    $batch = "restoredefaultpath; addpath('$matlabPathLiteral'); addpath('$testDirectoryLiteral'); result=test_standalone_simulink('$matlabPathLiteral'); assert(result.ok); fileId=fopen('$simulinkEvidenceLiteral','w'); assert(fileId>=0); cleanupFile=onCleanup(@()fclose(fileId)); fprintf(fileId,'%s\n',jsonencode(result,PrettyPrint=true)); clear cleanupFile;"

    # MathWorks' online license service can transiently reject an otherwise
    # valid batch start with error 5202. Retry only that startup failure; a
    # MATLAB assertion, numerical mismatch, or any other error still fails on
    # the first attempt.
    $matlabLog = Join-Path $resolvedRunRoot 'matlab-batch.log'
    $maxMatlabAttempts = 3
    for ($attempt = 1; $attempt -le $maxMatlabAttempts; $attempt++) {
        & $MatlabExecutable -batch $batch 2>&1 | Tee-Object -FilePath $matlabLog
        $matlabExitCode = $LASTEXITCODE
        if ($matlabExitCode -eq 0) { break }

        $matlabOutput = Get-Content -LiteralPath $matlabLog -Raw -ErrorAction SilentlyContinue
        $isLicenseService5202 = $matlabOutput -match '(?<!\d)5202(?!\d)'
        if (-not $isLicenseService5202 -or $attempt -eq $maxMatlabAttempts) {
            throw "Installed-wheel Simulink test failed with exit code $matlabExitCode"
        }

        $delaySeconds = 10 * $attempt
        Write-Warning "MathWorks license service returned 5202; retrying the same MATLAB batch command in $delaySeconds seconds (attempt $($attempt + 1)/$maxMatlabAttempts)."
        Start-Sleep -Seconds $delaySeconds
    }
    $simulinkEvidence = Get-Content -LiteralPath $simulinkEvidencePath -Raw |
        ConvertFrom-Json
    $evidence = [ordered]@{
        schema = 'radia-optuna.installed-wheel-evidence.v1'
        ok = $true
        wheel_verification = $wheelVerification
        doctor = $doctorEvidence
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
