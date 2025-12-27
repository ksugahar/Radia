#==============================================================================
# BuildMSVC.ps1 - Build Radia with MSVC + Intel MKL
#
# This script builds both radia.pyd and radia_ngsolve.pyd using MSVC compiler
# with Intel MKL for BLAS/LAPACK operations.
#
# Use this script when building radia_ngsolve, which requires MSVC for
# compatibility with MSVC-compiled NGSolve libraries.
#
# Usage:
#   powershell.exe -ExecutionPolicy Bypass -File BuildMSVC.ps1
#   powershell.exe -ExecutionPolicy Bypass -File BuildMSVC.ps1 -Rebuild
#   powershell.exe -ExecutionPolicy Bypass -File BuildMSVC.ps1 -Test
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

# Intel MKL path (required for BLAS/LAPACK)
$INTEL_ONEAPI = "C:\Program Files (x86)\Intel\oneAPI"
$INTEL_MKL = "$INTEL_ONEAPI\mkl\latest"

# Verify Intel MKL is installed
if (-not (Test-Path "$INTEL_MKL\lib\mkl_rt.lib")) {
    Write-Host ""
    Write-Host "========================================" -ForegroundColor Red
    Write-Host "  ERROR: Intel MKL Not Found" -ForegroundColor Red
    Write-Host "========================================" -ForegroundColor Red
    Write-Host ""
    Write-Host "Intel MKL is REQUIRED for BLAS/LAPACK operations." -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Expected location:" -ForegroundColor Gray
    Write-Host "  $INTEL_MKL\lib\mkl_rt.lib" -ForegroundColor Gray
    Write-Host ""
    Write-Host "Please install Intel oneAPI Base Toolkit" -ForegroundColor Cyan
    Write-Host "Download from: https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html" -ForegroundColor Gray
    Write-Host ""
    exit 1
}

Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "  Radia Build (MSVC + Intel MKL)" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Compiler: MSVC (Visual Studio 2022)" -ForegroundColor Gray
Write-Host "MKL:      $INTEL_MKL" -ForegroundColor Gray
Write-Host ""

# Project directory
$PROJECT_DIR = Split-Path -Parent $MyInvocation.MyCommand.Path
$BUILD_DIR = "$PROJECT_DIR\build-msvc"

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

# OpenMP flag
if ($NoOpenMP) {
    $OPENMP_FLAG = "OFF"
    Write-Host "OpenMP: DISABLED (debug mode)" -ForegroundColor Yellow
} else {
    $OPENMP_FLAG = "ON"
    Write-Host "OpenMP: ENABLED" -ForegroundColor Green
}

# Create batch file to run with Visual Studio environment
$BatchContent = @"
@echo off
REM Set up Visual Studio environment
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1

REM Add MKL paths
set LIB=$INTEL_MKL\lib;%LIB%
set INCLUDE=$INTEL_MKL\include;%INCLUDE%
set MKLROOT=$INTEL_MKL

cd /d "$BUILD_DIR"

echo Configuring CMake with MSVC...
"$CMAKE_EXE" "$PROJECT_DIR" ^
    -G "Ninja" ^
    -DCMAKE_C_COMPILER=cl ^
    -DCMAKE_CXX_COMPILER=cl ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DRADIA_ENABLE_OPENMP=$OPENMP_FLAG

if errorlevel 1 exit /b 1

echo Building radia...
"$CMAKE_EXE" --build . --config Release --target radia -j

if errorlevel 1 exit /b 1

echo Building radia_ngsolve...
"$CMAKE_EXE" --build . --config Release --target radia_ngsolve -j

if errorlevel 1 exit /b 1

echo Build completed.
"@

$BatchFile = "$PROJECT_DIR\build_temp_msvc.bat"
$BatchContent | Out-File -FilePath $BatchFile -Encoding ascii

try {
    # Run the batch file
    Write-Host "Building..." -ForegroundColor Cyan
    cmd /c $BatchFile
    $BuildResult = $LASTEXITCODE

    if ($BuildResult -ne 0) {
        throw "Build failed with exit code $BuildResult"
    }

    # Copy radia.pyd to package directory
    $PYD_SOURCE = "$BUILD_DIR\radia.cp312-win_amd64.pyd"
    $PYD_DEST = "$PROJECT_DIR\src\radia\radia.pyd"

    if (Test-Path $PYD_SOURCE) {
        Copy-Item $PYD_SOURCE $PYD_DEST -Force
        Write-Host "Copied radia.pyd to src/radia/" -ForegroundColor Green
    } else {
        Write-Host "WARNING: Could not find $PYD_SOURCE" -ForegroundColor Yellow
    }

    # Copy radia_ngsolve.pyd to package directory
    # NGSolve uses add_ngsolve_python_module which may not add version tag
    $NGSOLVE_PYD_SOURCE = "$BUILD_DIR\radia_ngsolve.pyd"
    $NGSOLVE_PYD_DEST = "$PROJECT_DIR\src\radia\radia_ngsolve.pyd"

    if (Test-Path $NGSOLVE_PYD_SOURCE) {
        Copy-Item $NGSOLVE_PYD_SOURCE $NGSOLVE_PYD_DEST -Force
        Write-Host "Copied radia_ngsolve.pyd to src/radia/" -ForegroundColor Green
    } else {
        # Try with version tag
        $NGSOLVE_PYD_SOURCE_VERSIONED = "$BUILD_DIR\radia_ngsolve.cp312-win_amd64.pyd"
        if (Test-Path $NGSOLVE_PYD_SOURCE_VERSIONED) {
            Copy-Item $NGSOLVE_PYD_SOURCE_VERSIONED $NGSOLVE_PYD_DEST -Force
            Write-Host "Copied radia_ngsolve.pyd to src/radia/" -ForegroundColor Green
        } else {
            Write-Host "WARNING: Could not find radia_ngsolve.pyd" -ForegroundColor Yellow
        }
    }

    # Copy Intel MKL DLLs for distribution
    $INTEL_DLLS = @(
        "$INTEL_MKL\bin\mkl_rt.2.dll"         # Intel MKL runtime
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
    python -c "import sys; sys.path.insert(0, r'$PROJECT_DIR\src\radia'); import radia; print(f'radia version: {radia.__version__}')"

    if ($LASTEXITCODE -eq 0) {
        Write-Host "radia import test passed!" -ForegroundColor Green
    } else {
        Write-Host "radia import test failed!" -ForegroundColor Red
    }

    python -c "import sys; sys.path.insert(0, r'$PROJECT_DIR\src\radia'); import radia_ngsolve; print('radia_ngsolve import OK')"

    if ($LASTEXITCODE -eq 0) {
        Write-Host "radia_ngsolve import test passed!" -ForegroundColor Green
    } else {
        Write-Host "radia_ngsolve import test failed!" -ForegroundColor Red
    }
    Write-Host ""
}

# Show usage
Write-Host "Usage:" -ForegroundColor Yellow
Write-Host "  import radia" -ForegroundColor Gray
Write-Host "  import radia_ngsolve" -ForegroundColor Gray
Write-Host "  print(radia.__version__)" -ForegroundColor Gray
Write-Host ""
