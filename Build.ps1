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
Write-Host "[2/2] Building radia_ngsolve.pyd..." -ForegroundColor Cyan

# Load VS environment and build radia_ngsolve
& powershell -ExecutionPolicy Bypass -Command {
    Import-Module 'C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\Tools\Microsoft.VisualStudio.DevShell.dll'
    Enter-VsDevShell -VsInstallPath 'C:\Program Files\Microsoft Visual Studio\2022\Community' -SkipAutomaticLocation
    cmake --build build --config $using:BuildType --target radia_ngsolve
}

if ($LASTEXITCODE -ne 0) {
    Write-Host "[ERROR] Failed to build radia_ngsolve.pyd" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "========================================" -ForegroundColor Green
Write-Host "  Build Complete" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
Write-Host ""

# Copy .pyd files to src/radia for PyPI packaging
Write-Host "Copying to src/radia for PyPI packaging..." -ForegroundColor Cyan

# Create src/radia directory if it doesn't exist
if (-not (Test-Path "src\radia")) {
    New-Item -ItemType Directory -Path "src\radia" -Force | Out-Null
}

# Copy radia.pyd - check multiple possible locations
$radiaSrcPaths = @(
    "build\$BuildType\radia.cp312-win_amd64.pyd",
    "build\$BuildType\radia.pyd",
    "build\radia.cp312-win_amd64.pyd",
    "build\radia.pyd"
)
$radiaDst = "src\radia\radia.pyd"
$radiaCopied = $false

foreach ($radiaSrc in $radiaSrcPaths) {
    if (Test-Path $radiaSrc) {
        Copy-Item $radiaSrc $radiaDst -Force
        Write-Host "  [OK] Copied radia.pyd to src/radia/ (from $radiaSrc)" -ForegroundColor Green
        $radiaCopied = $true
        break
    }
}

if (-not $radiaCopied) {
    Write-Host "  [WARN] radia.pyd not found in build directory" -ForegroundColor Yellow
}

# Copy radia_ngsolve.pyd - check multiple possible locations
$ngSolveSrcPaths = @(
    "build\$BuildType\radia_ngsolve.pyd",
    "build\radia_ngsolve.pyd"
)
$ngSolveDst = "src\radia\radia_ngsolve.pyd"
$ngSolveCopied = $false

foreach ($ngSolveSrc in $ngSolveSrcPaths) {
    if (Test-Path $ngSolveSrc) {
        Copy-Item $ngSolveSrc $ngSolveDst -Force
        Write-Host "  [OK] Copied radia_ngsolve.pyd to src/radia/ (from $ngSolveSrc)" -ForegroundColor Green
        $ngSolveCopied = $true
        break
    }
}

if (-not $ngSolveCopied) {
    Write-Host "  [WARN] radia_ngsolve.pyd not found in build directory" -ForegroundColor Yellow
}

Write-Host ""

# Show outputs
Write-Host "Built modules:" -ForegroundColor Yellow

# Check for radia.pyd in build directory
$radiaFound = $false
foreach ($radiaSrc in $radiaSrcPaths) {
    if (Test-Path $radiaSrc) {
        $size = [math]::Round((Get-Item $radiaSrc).Length / 1KB, 1)
        Write-Host "  [OK] radia.pyd            ($size KB) - $radiaSrc" -ForegroundColor Green
        $radiaFound = $true
        break
    }
}
if (-not $radiaFound) {
    Write-Host "  [ ] radia.pyd            (not found)" -ForegroundColor Gray
}

# Check for radia_ngsolve.pyd in build directory
$ngSolveFound = $false
foreach ($ngSolveSrc in $ngSolveSrcPaths) {
    if (Test-Path $ngSolveSrc) {
        $size = [math]::Round((Get-Item $ngSolveSrc).Length / 1KB, 1)
        Write-Host "  [OK] radia_ngsolve.pyd    ($size KB) - $ngSolveSrc" -ForegroundColor Green
        $ngSolveFound = $true
        break
    }
}
if (-not $ngSolveFound) {
    Write-Host "  [ ] radia_ngsolve.pyd    (not found)" -ForegroundColor Gray
}

Write-Host ""

# Run tests if requested
if ($Test) {
    Write-Host "Running tests..." -ForegroundColor Cyan

    Write-Host "  Testing radia..." -ForegroundColor Gray
    python -c "import sys; sys.path.insert(0, r'build\$BuildType'); import radia as rad; print(f'  [OK] radia {rad.UtiVer()}')"

    Write-Host "  Testing radia_ngsolve..." -ForegroundColor Gray
    python -c "import sys; sys.path.insert(0, r'build\$BuildType'); import ngsolve; import radia_ngsolve; print('  [OK] radia_ngsolve')"

    Write-Host ""
}

Write-Host "Usage:" -ForegroundColor Yellow
Write-Host "  python -c `"import sys; sys.path.insert(0, r'build\$BuildType'); import radia as rad; print(rad.UtiVer())`"" -ForegroundColor Gray
Write-Host ""
