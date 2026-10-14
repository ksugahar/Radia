param(
    [ValidateRange(1, 101)]
    [int]$Repeats = 7,
    [string]$Output = ""
)

$ErrorActionPreference = "Stop"
$RepoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$ProbeRoot = $RepoRoot
foreach ($Drive in Get-PSDrive -PSProvider FileSystem) {
    if ($Drive.DisplayRoot -and
            $RepoRoot.StartsWith($Drive.DisplayRoot, [StringComparison]::OrdinalIgnoreCase)) {
        $Suffix = $RepoRoot.Substring($Drive.DisplayRoot.Length).TrimStart("\")
        $ProbeRoot = Join-Path $Drive.Root $Suffix
        break
    }
}
$MatlabDirPath = Join-Path $ProbeRoot "matlab"
$MatlabDir = $MatlabDirPath.Replace("\", "/").Replace("'", "''")
$ManifestPath = Join-Path $ProbeRoot "packages\radia-optuna\src\radia_optuna\manifest.json"
if (-not (Test-Path -LiteralPath $ManifestPath)) {
    throw "Missing radia-optuna manifest: $ManifestPath"
}
$ExpectedCommandCount = [int](
    Get-Content -LiteralPath $ManifestPath -Raw | ConvertFrom-Json
).native_command_count
$Batch = @"
addpath('$MatlabDir'); clear optuna_mex; started=tic; info=optuna_mex('api.info'); elapsed=toc(started); assert(info.command_count==$ExpectedCommandCount); fprintf('OPTUNA_MEX_FIRST_CALL_S=%.9f\n',elapsed)
"@

$Measurements = @()
foreach ($Repeat in 1..$Repeats) {
    $Transcript = & matlab -batch $Batch 2>&1
    if ($LASTEXITCODE -ne 0) {
        throw "MATLAB cold-start probe $Repeat failed:`n$($Transcript -join [Environment]::NewLine)"
    }
    $Sentinel = $Transcript | Where-Object {
        $_ -match '^OPTUNA_MEX_FIRST_CALL_S=([0-9.]+)$'
    } | Select-Object -Last 1
    if (-not $Sentinel -or $Sentinel -notmatch '=([0-9.]+)$') {
        throw "MATLAB cold-start probe $Repeat did not report its sentinel."
    }
    $Measurements += [double]::Parse(
        $Matches[1], [Globalization.CultureInfo]::InvariantCulture)
}

$Binary = Join-Path $RepoRoot "matlab\optuna_mex.mexw64"
if (-not (Test-Path -LiteralPath $Binary)) {
    throw "Missing $Binary. Run Build.ps1 -OptunaMexOnly first."
}
$Sorted = @($Measurements | Sort-Object)
$Median = if ($Sorted.Count % 2 -eq 1) {
    $Sorted[[int][Math]::Floor($Sorted.Count / 2)]
} else {
    ($Sorted[$Sorted.Count / 2 - 1] + $Sorted[$Sorted.Count / 2]) / 2
}
$Result = [ordered]@{
    schema = "radia.validation.optuna-mex-first-call.v1"
    generated_at = (Get-Date).ToUniversalTime().ToString("o")
    host = [Environment]::MachineName
    repeats = $Repeats
    first_call_seconds = $Measurements
    median_first_call_seconds = $Median
    binary_size_bytes = (Get-Item -LiteralPath $Binary).Length
    binary_sha256 = (Get-FileHash -LiteralPath $Binary -Algorithm SHA256).Hash.ToLowerInvariant()
    gateway_directory = $MatlabDirPath
    boundary = "First optuna_mex api.info call in each fresh MATLAB process; excludes MATLAB executable startup."
}
$Encoded = $Result | ConvertTo-Json -Depth 4
if ($Output) {
    [IO.File]::WriteAllText(
        [IO.Path]::GetFullPath($Output), $Encoded + [Environment]::NewLine,
        [Text.UTF8Encoding]::new($false))
}
$Encoded
