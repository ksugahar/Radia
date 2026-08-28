# Final, non-interactive Eqnedit64 release gate.
#
# This is intentionally stricter than build_tests.bat: it adds ASan, the
# extended real-window-procedure endurance run and release identity checks.
# Every UI process is started hidden by the underlying tests;
# no SendInput, cursor movement, or foreground activation is used.

[CmdletBinding()]
param(
    [ValidateRange(1, 10000)]
    [int]$EnduranceSeeds = 100,
    [ValidateRange(1, 1000000)]
    [int]$EnduranceOps = 5000,
    [switch]$SkipBuild,
    [string]$ReportPath = 'C:\temp\Eqnedit64-final-acceptance.json'
)

$ErrorActionPreference = 'Stop'
if ($env:EQNEDIT64_ISOLATED_TEST_SESSION -ne '1') {
    throw @'
The full release gate starts and exits several private-font processes.
Run it only in a disposable CI, VM, or dedicated Windows user session, then
set EQNEDIT64_ISOLATED_TEST_SESSION=1 in that isolated session.
'@
}
$root = Split-Path -Parent $PSScriptRoot
$buildExe = Join-Path $root 'build\Eqnedit64.exe'
$distExe = Join-Path $root 'dist\Eqnedit64.exe'
$started = Get-Date

function Invoke-Batch([string]$RelativePath) {
    Push-Location $root
    try {
        & cmd.exe /d /c $RelativePath
        if ($LASTEXITCODE -ne 0) {
            throw "$RelativePath failed with exit code $LASTEXITCODE"
        }
    } finally {
        Pop-Location
    }
}

function Invoke-PowerShellScript([string]$RelativePath, [object[]]$Arguments) {
    $script = Join-Path $root $RelativePath
    & pwsh -NoProfile -ExecutionPolicy Bypass -File $script @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "$RelativePath failed with exit code $LASTEXITCODE"
    }
}

if (-not $SkipBuild) {
    Invoke-Batch 'build\build_tests.bat'
    Invoke-Batch 'build\build_asan.bat'
    Invoke-PowerShellScript 'build\test_asan.ps1' @()
}

foreach ($required in @($buildExe, $distExe)) {
    if (-not (Test-Path -LiteralPath $required -PathType Leaf)) {
        throw "Required release input is missing: $required"
    }
}

# The long normal-build run is separate from ASan: it exercises the signed,
# portable executable that will actually be deployed.
Invoke-PowerShellScript 'build\test_ui_fuzz.ps1' @(
    '-AppPath', $distExe, '-Seeds', $EnduranceSeeds, '-Ops', $EnduranceOps,
    '-TimeoutSeconds', 120)

$paintReport = 'C:\temp\Eqnedit64-paint-benchmark.txt'
$paintError = 'C:\temp\Eqnedit64-paint-benchmark.err.txt'
$paint = Start-Process -FilePath $distExe -ArgumentList '--paint-bench' `
    -WindowStyle Hidden -Wait -PassThru `
    -RedirectStandardOutput $paintReport -RedirectStandardError $paintError
if ($paint.ExitCode -ne 0) {
    $details = Get-Content -LiteralPath $paintReport -Raw -ErrorAction SilentlyContinue
    throw "Offscreen paint performance gate failed ($($paint.ExitCode)):`n$details"
}
$paintMeasurements = Get-Content -LiteralPath $paintReport -Raw

$releaseFiles = @(Get-ChildItem -LiteralPath (Split-Path -Parent $distExe) -File)
if ($releaseFiles.Count -ne 1 -or $releaseFiles[0].Name -cne 'Eqnedit64.exe') {
    throw "dist must contain Eqnedit64.exe only: $($releaseFiles.Name -join ', ')"
}
$buildHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $buildExe).Hash
$distHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $distExe).Hash
if ($buildHash -ne $distHash) {
    throw 'dist\Eqnedit64.exe is not the exact executable exercised by the gate.'
}

$signature = Get-AuthenticodeSignature -LiteralPath $distExe
if ($signature.Status -ne 'Valid' -or
    $signature.SignerCertificate.Subject -cne 'CN=ksugahar' -or
    -not ($signature.SignerCertificate.EnhancedKeyUsageList.ObjectId -contains
        '1.3.6.1.5.5.7.3.3')) {
    throw "Release signature is invalid: $($signature.Status) / $($signature.SignerCertificate.Subject)"
}

$stamp = Join-Path $root 'src\build_stamp.h'
$buildCommit = 'unknown'
if (Test-Path -LiteralPath $stamp) {
    $match = Select-String -LiteralPath $stamp -Pattern 'BUILD_COMMIT' |
        Select-Object -First 1
    if ($match) { $buildCommit = $match.Line -replace '.*"(.*)".*', '$1' }
}
$gitCommit = (& git -C $root rev-parse HEAD).Trim()
$report = [ordered]@{
    schema = 'eqnedit64.final-acceptance.v1'
    accepted_at = (Get-Date).ToString('o')
    elapsed_seconds = [math]::Round(((Get-Date) - $started).TotalSeconds, 3)
    source_commit = $gitCommit
    executable_build_commit = $buildCommit
    executable = $distExe
    sha256 = $distHash
    signer = $signature.SignerCertificate.Subject
    distribution_files = @($releaseFiles.Name)
    gui_endurance = [ordered]@{
        seeds = $EnduranceSeeds
        operations_per_seed = $EnduranceOps
        total_operations = $EnduranceSeeds * $EnduranceOps
    }
    paint_benchmark = $paintMeasurements.Trim()
    background_suite = if ($SkipBuild) { 'skipped by caller' } else { 'passed' }
    asan_suite = if ($SkipBuild) { 'skipped by caller' } else { 'passed' }
}
$reportFolder = Split-Path -Parent $ReportPath
if (-not (Test-Path -LiteralPath $reportFolder)) {
    New-Item -ItemType Directory -Path $reportFolder | Out-Null
}
$report | ConvertTo-Json -Depth 5 |
    Set-Content -LiteralPath $ReportPath -Encoding utf8

Write-Host 'PASS: Eqnedit64 final acceptance gate'
Write-Host "  commit: $buildCommit"
Write-Host "  sha256: $distHash"
Write-Host "  signer: $($signature.SignerCertificate.Subject)"
Write-Host "  endurance: $EnduranceSeeds x $EnduranceOps = $($EnduranceSeeds * $EnduranceOps) operations"
Write-Host "  paint benchmark: all cached cases below 5 ms"
Write-Host "  report: $ReportPath"
