#==============================================================================
# Build.ps1 - Build Radia with Intel oneAPI compiler
#
# IMPORTANT: Intel oneAPI compiler is REQUIRED for building Radia.
# MSVC is NOT supported.
#
# Usage:
#   powershell.exe -ExecutionPolicy Bypass -File Build.ps1
#   powershell.exe -ExecutionPolicy Bypass -File Build.ps1 -Rebuild
#   powershell.exe -ExecutionPolicy Bypass -File Build.ps1 -Test
#
# Options:
#   -Rebuild  Clean build directory before building
#   -Test     Run import test after build
#==============================================================================

param(
    [switch]$Rebuild,
    [switch]$Test,
    [switch]$NoOpenMP   # Disable OpenMP (for debugging)
)

$ErrorActionPreference = "Stop"

# Intel oneAPI paths
$INTEL_ONEAPI = "C:\Program Files (x86)\Intel\oneAPI"
$INTEL_COMPILER = "$INTEL_ONEAPI\compiler\latest"
$INTEL_MKL = "$INTEL_ONEAPI\mkl\latest"

# Verify Intel oneAPI is installed
if (-not (Test-Path "$INTEL_COMPILER\bin\icx-cl.exe")) {
    Write-Host ""
    Write-Host "========================================" -ForegroundColor Red
    Write-Host "  ERROR: Intel oneAPI Not Found" -ForegroundColor Red
    Write-Host "========================================" -ForegroundColor Red
    Write-Host ""
    Write-Host "Intel oneAPI compiler is REQUIRED for building Radia." -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Expected location:" -ForegroundColor Gray
    Write-Host "  $INTEL_COMPILER\bin\icx-cl.exe" -ForegroundColor Gray
    Write-Host ""
    Write-Host "Please install:" -ForegroundColor Cyan
    Write-Host "  1. Intel oneAPI Base Toolkit" -ForegroundColor White
    Write-Host "  2. Intel oneAPI HPC Toolkit" -ForegroundColor White
    Write-Host ""
    Write-Host "Download from: https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html" -ForegroundColor Gray
    Write-Host ""
    exit 1
}

Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "  Radia Build (Intel oneAPI)" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Compiler: $INTEL_COMPILER\bin\icx-cl.exe" -ForegroundColor Gray
Write-Host "MKL:      $INTEL_MKL" -ForegroundColor Gray
Write-Host ""

# Project directory
$PROJECT_DIR = Split-Path -Parent $MyInvocation.MyCommand.Path
$BUILD_DIR = "$PROJECT_DIR\build-intel"

# Clean build if requested
if ($Rebuild -and (Test-Path $BUILD_DIR)) {
    Write-Host "Cleaning build directory..." -ForegroundColor Yellow
    Remove-Item -Recurse -Force $BUILD_DIR
}

# Create build directory
if (-not (Test-Path $BUILD_DIR)) {
    New-Item -ItemType Directory -Path $BUILD_DIR | Out-Null
}

# CMake executable
$CMAKE_EXE = "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
if (-not (Test-Path $CMAKE_EXE)) {
    Write-Host "ERROR: CMake not found at $CMAKE_EXE" -ForegroundColor Red
    Write-Host "Please install Visual Studio 2022 with CMake component" -ForegroundColor Yellow
    exit 1
}

# Intel compiler paths
$ICX_CL = "$INTEL_COMPILER\bin\icx-cl.exe"
$INTEL_LIB = "$INTEL_COMPILER\lib"

# OpenMP flag
if ($NoOpenMP) {
    $OPENMP_FLAG = "OFF"
    Write-Host "OpenMP: DISABLED (debug mode)" -ForegroundColor Yellow
} else {
    $OPENMP_FLAG = "ON"
    Write-Host "OpenMP: ENABLED" -ForegroundColor Green
}

# Create batch file to run with Intel environment
$BatchContent = @"
@echo off
REM Set up Visual Studio environment (suppress vswhere.exe warning)
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1

REM Add Intel compiler paths
set PATH=$INTEL_COMPILER\bin;%PATH%
set LIB=$INTEL_LIB;%LIB%
set INCLUDE=$INTEL_COMPILER\include;%INCLUDE%

REM Add MKL paths
set LIB=$INTEL_MKL\lib;%LIB%
set INCLUDE=$INTEL_MKL\include;%INCLUDE%
set MKLROOT=$INTEL_MKL

cd /d "$BUILD_DIR"

echo Configuring CMake with Intel compiler...
"$CMAKE_EXE" "$PROJECT_DIR" ^
    -G "Ninja" ^
    -DCMAKE_C_COMPILER="$ICX_CL" ^
    -DCMAKE_CXX_COMPILER="$ICX_CL" ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DRADIA_ENABLE_OPENMP=$OPENMP_FLAG

if errorlevel 1 exit /b 1

echo Building radia...
"$CMAKE_EXE" --build . --config Release --target radia -j

if errorlevel 1 exit /b 1

echo Build completed.
"@

$BatchFile = "$PROJECT_DIR\build_temp.bat"
$BatchContent | Out-File -FilePath $BatchFile -Encoding ascii

try {
    # Run the batch file
    Write-Host "Building..." -ForegroundColor Cyan
    cmd /c $BatchFile
    $BuildResult = $LASTEXITCODE

    if ($BuildResult -ne 0) {
        throw "Build failed with exit code $BuildResult"
    }

    # Copy .pyd to package directory
    $PYD_SOURCE = "$BUILD_DIR\radia.cp312-win_amd64.pyd"
    $PYD_DEST = "$PROJECT_DIR\src\radia\radia.pyd"

    if (Test-Path $PYD_SOURCE) {
        Copy-Item $PYD_SOURCE $PYD_DEST -Force
        Write-Host "Copied radia.pyd to src/radia/" -ForegroundColor Green
    } else {
        Write-Host "WARNING: Could not find $PYD_SOURCE" -ForegroundColor Yellow
    }

    # Copy Intel DLLs for distribution
    $INTEL_DLLS = @(
        "$INTEL_MKL\bin\mkl_rt.2.dll",        # Intel MKL runtime
        "$INTEL_MKL\bin\mkl_vml_avx2.2.dll",  # MKL Vector Math Library AVX2
        "$INTEL_MKL\bin\mkl_vml_def.2.dll",   # MKL Vector Math Library default
        "$INTEL_COMPILER\bin\libiomp5md.dll", # Intel OpenMP
        "$INTEL_COMPILER\bin\libmmd.dll",     # Intel Math library
        "$INTEL_COMPILER\bin\svml_dispmd.dll" # Intel Short Vector Math Library
    )

    foreach ($dll in $INTEL_DLLS) {
        if (Test-Path $dll) {
            $dllName = Split-Path -Leaf $dll
            Copy-Item $dll "$PROJECT_DIR\src\radia\" -Force
            Write-Host "Copied $dllName to src/radia/" -ForegroundColor Green
        } else {
            Write-Host "WARNING: DLL not found: $dll" -ForegroundColor Yellow
        }
    }

    Write-Host ""
    Write-Host "========================================" -ForegroundColor Green
    Write-Host "  Build Completed Successfully" -ForegroundColor Green
    Write-Host "========================================" -ForegroundColor Green
    Write-Host ""
}
catch {
    Write-Host ""
    Write-Host "========================================" -ForegroundColor Red
    Write-Host "  Build Failed" -ForegroundColor Red
    Write-Host "========================================" -ForegroundColor Red
    Write-Host ""
    Write-Host "Error: $_" -ForegroundColor Red
    Write-Host ""
    exit 1
}
finally {
    # Clean up temp batch file
    if (Test-Path $BatchFile) {
        Remove-Item $BatchFile -Force
    }
}

# Test import if requested
if ($Test) {
    Write-Host "Testing import..." -ForegroundColor Cyan
    python -c "import radia; print(f'radia version: {radia.__version__}')"

    if ($LASTEXITCODE -eq 0) {
        Write-Host "Import test passed!" -ForegroundColor Green
    } else {
        Write-Host "Import test failed!" -ForegroundColor Red
    }
    Write-Host ""
}

# Show usage
Write-Host "Usage:" -ForegroundColor Yellow
Write-Host "  import radia" -ForegroundColor Gray
Write-Host "  print(radia.__version__)" -ForegroundColor Gray
Write-Host ""
