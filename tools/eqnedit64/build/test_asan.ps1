param(
    [string]$AppPath,
    [ValidateRange(1, 10000)]
    [int]$Seeds = 24,
    [ValidateRange(1, 1000000)]
    [int]$Ops = 3000,
    [ValidateRange(1, 3600)]
    [int]$TimeoutSeconds = 600
)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
$exe = if ($AppPath) {
    (Resolve-Path -LiteralPath $AppPath -ErrorAction Stop).Path
} else {
    (Resolve-Path -LiteralPath (Join-Path $root `
        'build\Eqnedit64_asan.exe') -ErrorAction Stop).Path
}
$previousOptions = $env:ASAN_OPTIONS
$runtime = Join-Path (Split-Path -Parent $exe) 'clang_rt.asan_dynamic-x86_64.dll'
if (-not (Test-Path -LiteralPath $runtime -PathType Leaf)) {
    throw "MSVC ASan runtime missing beside diagnostic executable: $runtime. Rebuild with build_asan.bat."
}
$runtimeInfo = Get-Item -LiteralPath $runtime
$runtimeHash = (Get-FileHash -LiteralPath $runtime -Algorithm SHA256).Hash
Write-Host "ASan runtime: $runtime; version=$($runtimeInfo.VersionInfo.FileVersion); SHA256=$runtimeHash"

try {
    $env:ASAN_OPTIONS = 'abort_on_error=1:halt_on_error=1:detect_leaks=0'
    foreach ($mode in @('--self-test', '--status-layout-test',
            '--visual-scale-test', '--ui-interaction-test')) {
        $process = Start-Process -FilePath $exe -ArgumentList $mode `
            -PassThru -WindowStyle Hidden
        if (-not $process.WaitForExit(60000)) {
            $process.Kill()
            throw "ASan $mode timed out."
        }
        if ($process.ExitCode -ne 0) {
            throw "ASan $mode failed with exit code $($process.ExitCode)."
        }
        Write-Host "PASS: ASan $mode"
    }

    & (Join-Path $PSScriptRoot 'test_ui_fuzz.ps1') -AppPath $exe `
        -Seeds $Seeds -Ops $Ops -TimeoutSeconds $TimeoutSeconds
    Write-Host 'PASS: ASan hidden tests and GUI fuzz'
} finally {
    $env:ASAN_OPTIONS = $previousOptions
}
