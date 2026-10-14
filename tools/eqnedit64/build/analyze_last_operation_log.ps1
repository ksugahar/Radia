param(
    [string]$LogPath,
    [ValidateSet('structure', 'full')]
    [string]$Privacy = 'structure',
    [string]$OutputPath
)

$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot

if (-not $LogPath) {
    $logDirectory = Join-Path $env:LOCALAPPDATA 'Eqnedit64\logs'
    $latest = Get-ChildItem -LiteralPath $logDirectory -Filter 'operation-*.log' -File |
        Sort-Object LastWriteTime -Descending |
        Select-Object -First 1
    if (-not $latest) {
        throw "No Eqnedit64 operation log was found in $logDirectory"
    }
    $LogPath = $latest.FullName
}

if (-not (Test-Path -LiteralPath $LogPath)) {
    throw "Operation log does not exist: $LogPath"
}

if (-not $OutputPath) {
    $stamp = Get-Date -Format 'yyyyMMdd-HHmmss'
    $OutputPath = "C:\temp\Eqnedit64-usability-$stamp.json"
}

$analyzer = Join-Path $projectRoot 'tools\analyze_usability_trace.py'
python $analyzer $LogPath --output $OutputPath --privacy $Privacy
if ($LASTEXITCODE -ne 0) {
    throw "Usability trace analysis failed with exit code $LASTEXITCODE"
}

Write-Host "LLM usability-review bundle: $OutputPath"
