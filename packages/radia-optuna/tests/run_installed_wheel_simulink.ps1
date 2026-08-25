param(
    [Parameter(Mandatory = $true)]
    [string]$Wheel,
    [string]$MatlabExecutable = 'matlab',
    [string]$PythonExecutable = 'python'
)

$ErrorActionPreference = 'Stop'
$wheelPath = (Resolve-Path -LiteralPath $Wheel).ProviderPath
$testDirectory = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot 'matlab')).ProviderPath
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
    $venv = Join-Path $resolvedRunRoot 'venv'
    & $PythonExecutable -m venv $venv
    if ($LASTEXITCODE -ne 0) { throw "Python venv creation failed with exit code $LASTEXITCODE" }

    $venvPython = Join-Path $venv 'Scripts\python.exe'
    & $venvPython -m pip install --disable-pip-version-check --no-deps $wheelPath
    if ($LASTEXITCODE -ne 0) { throw "Installing the wheel failed with exit code $LASTEXITCODE" }

    $doctor = Join-Path $venv 'Scripts\radia-optuna-doctor.exe'
    & $doctor --json
    if ($LASTEXITCODE -ne 0) { throw "radia-optuna-doctor failed with exit code $LASTEXITCODE" }

    $installedMatlabPath = (& $venvPython -c 'from radia_optuna import matlab_path; print(matlab_path())').Trim()
    if (-not (Test-Path -LiteralPath $installedMatlabPath -PathType Container)) {
        throw "Installed MATLAB directory is missing: $installedMatlabPath"
    }

    $matlabPathLiteral = ConvertTo-MatlabLiteral $installedMatlabPath
    $testDirectoryLiteral = ConvertTo-MatlabLiteral $testDirectory
    $batch = "restoredefaultpath; addpath('$matlabPathLiteral'); addpath('$testDirectoryLiteral'); result=test_standalone_simulink('$matlabPathLiteral'); assert(result.ok);"

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
    Write-Output 'RADIA_OPTUNA_WHEEL_SIMULINK_OK'
} finally {
    $checkedRunRoot = [IO.Path]::GetFullPath($resolvedRunRoot)
    if ($checkedRunRoot.StartsWith($resolvedTempRoot, [StringComparison]::OrdinalIgnoreCase) -and
            (Test-Path -LiteralPath $checkedRunRoot)) {
        Remove-Item -LiteralPath $checkedRunRoot -Recurse -Force
    }
}
