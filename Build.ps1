#==============================================================================
# Build - Unified build script for Radia project
# Usage: .\Build.ps1 [-BuildType Release|Debug] [-Clean] [-Test]
#==============================================================================

param(
    [ValidateSet("Release", "Debug", "RelWithDebInfo")]
    [string]$BuildType = "Release",

    [switch]$Clean,
    [switch]$Test
)

Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "  Radia Project - Build" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Build Type: $BuildType" -ForegroundColor White
Write-Host ""

# Build radia.pyd
Write-Host "[1/2] Building radia.pyd..." -ForegroundColor Cyan
$radiaArgs = @("-BuildType", $BuildType)
if ($Clean) { $radiaArgs += "-Clean" }

& powershell -ExecutionPolicy Bypass -File .\BuildRadiaInternal.ps1 @radiaArgs

if ($LASTEXITCODE -ne 0) {
    Write-Host "[ERROR] Failed to build radia.pyd" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "[2/2] Building rad_ngsolve.pyd..." -ForegroundColor Cyan

# Load VS environment and build rad_ngsolve
& powershell -ExecutionPolicy Bypass -Command {
    Import-Module 'C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\Tools\Microsoft.VisualStudio.DevShell.dll'
    Enter-VsDevShell -VsInstallPath 'C:\Program Files\Microsoft Visual Studio\2022\Community' -SkipAutomaticLocation
    cmake --build build --config $using:BuildType --target rad_ngsolve
}

if ($LASTEXITCODE -ne 0) {
    Write-Host "[ERROR] Failed to build rad_ngsolve.pyd" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "========================================" -ForegroundColor Green
Write-Host "  Build Complete" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
Write-Host ""

# Show outputs
Write-Host "Built modules:" -ForegroundColor Yellow

$radiaPath = "build\$BuildType\radia.pyd"
if (Test-Path $radiaPath) {
    $size = [math]::Round((Get-Item $radiaPath).Length / 1KB, 1)
    Write-Host "  [OK] radia.pyd        ($size KB)" -ForegroundColor Green
} else {
    Write-Host "  [ ] radia.pyd        (not found)" -ForegroundColor Gray
}

$ngSolvePath = "build\$BuildType\rad_ngsolve.pyd"
if (Test-Path $ngSolvePath) {
    $size = [math]::Round((Get-Item $ngSolvePath).Length / 1KB, 1)
    Write-Host "  [OK] rad_ngsolve.pyd  ($size KB)" -ForegroundColor Green
} else {
    Write-Host "  [ ] rad_ngsolve.pyd  (not found)" -ForegroundColor Gray
}

Write-Host ""

# Run tests if requested
if ($Test) {
    Write-Host "Running tests..." -ForegroundColor Cyan
    
    Write-Host "  Testing radia..." -ForegroundColor Gray
    python -c "import sys; sys.path.insert(0, r'build\$BuildType'); import radia as rad; print(f'  [OK] radia {rad.UtiVer()}')"
    
    Write-Host "  Testing rad_ngsolve..." -ForegroundColor Gray
    python -c "import sys; sys.path.insert(0, r'build\$BuildType'); import ngsolve; import rad_ngsolve; print('  [OK] rad_ngsolve')"
    
    Write-Host ""
}

Write-Host "Usage:" -ForegroundColor Yellow
Write-Host "  python -c `"import sys; sys.path.insert(0, r'build\$BuildType'); import radia as rad; print(rad.UtiVer())`"" -ForegroundColor Gray
Write-Host ""
