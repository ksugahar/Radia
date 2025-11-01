#==============================================================================
# rad_ngsolve Build Script
# Builds rad_ngsolve.pyd following build_rad_ngsolve.cmd pattern
#==============================================================================

param(
    [ValidateSet("Release", "Debug")]
    [string]$BuildType = "Release"
)

Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "  rad_ngsolve Build Script" -ForegroundColor Cyan  
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""

# Clean build directory
Write-Host "[1/4] Cleaning build directory..." -ForegroundColor Cyan
if (Test-Path build) {
    Remove-Item -Recurse -Force build
}
New-Item -ItemType Directory -Path build | Out-Null
Write-Host "      Build directory cleaned" -ForegroundColor Green

# Configure with CMake
Write-Host ""
Write-Host "[2/4] Configuring CMake..." -ForegroundColor Cyan
cmake -S . -B build -G "Visual Studio 17 2022"

if ($LASTEXITCODE -ne 0) {
    Write-Host "[ERROR] CMake configuration failed" -ForegroundColor Red
    exit 1
}
Write-Host "      Configuration complete" -ForegroundColor Green

# Build
Write-Host ""
Write-Host "[3/4] Building rad_ngsolve ($BuildType)..." -ForegroundColor Cyan
cmake --build build --config $BuildType --target rad_ngsolve

if ($LASTEXITCODE -ne 0) {
    Write-Host "[ERROR] Build failed" -ForegroundColor Red
    exit 1
}
Write-Host "      Build complete" -ForegroundColor Green

# Test module
Write-Host ""
Write-Host "[4/4] Testing module..." -ForegroundColor Cyan

$modulePath = "build/$BuildType/rad_ngsolve.pyd"
if (Test-Path $modulePath) {
    Write-Host "      Module found: $modulePath" -ForegroundColor Green
    
    # Test import
    $testScript = @"
import sys
sys.path.insert(0, r'build/$BuildType')
import ngsolve
import rad_ngsolve
print('[OK] rad_ngsolve imported successfully')
print('[OK] Available:', [x for x in dir(rad_ngsolve) if not x.startswith('_')])
"@
    
    $testScript | python -
    
    if ($LASTEXITCODE -eq 0) {
        Write-Host ""
        Write-Host "========================================" -ForegroundColor Green
        Write-Host "  BUILD SUCCESS" -ForegroundColor Green
        Write-Host "========================================" -ForegroundColor Green
        Write-Host ""
        Write-Host "Module: $modulePath" -ForegroundColor White
        Write-Host ""
        Write-Host "Usage:" -ForegroundColor Yellow
        Write-Host "  import ngsolve" -ForegroundColor Gray
        Write-Host "  import rad_ngsolve" -ForegroundColor Gray
        Write-Host ""
    } else {
        Write-Host "[WARN] Module built but import test failed" -ForegroundColor Yellow
    }
} else {
    Write-Host "[ERROR] Module not found at: $modulePath" -ForegroundColor Red
    exit 1
}
