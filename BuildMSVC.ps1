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
    # mkl_rt.X.dll is SDL (Single Dynamic Library) that loads other DLLs at runtime
    # Use wildcard patterns to be version-agnostic (e.g., mkl_rt.2.dll, mkl_rt.3.dll)

    Write-Host "Copying Intel MKL DLLs..." -ForegroundColor Cyan

    # MKL DLL patterns (version-agnostic)
    $MKL_DLL_PATTERNS = @(
        "mkl_rt.*.dll",            # MKL SDL runtime
        "mkl_core.*.dll",          # MKL core
        "mkl_intel_thread.*.dll",  # MKL threading (Intel OpenMP)
        "mkl_def.*.dll",           # Default CPU kernels
        "mkl_avx2.*.dll",          # AVX2 optimized kernels
        "mkl_vml_def.*.dll",       # Vector math library (default)
        "mkl_vml_avx2.*.dll"       # Vector math library (AVX2)
    )

    foreach ($pattern in $MKL_DLL_PATTERNS) {
        $dlls = Get-ChildItem -Path "$INTEL_MKL\bin" -Filter $pattern -ErrorAction SilentlyContinue
        if ($dlls) {
            foreach ($dll in $dlls) {
                Copy-Item $dll.FullName "$PROJECT_DIR\src\radia\" -Force
                Write-Host "  Copied $($dll.Name)" -ForegroundColor Green
            }
        } else {
            Write-Host "  WARNING: No DLL matching $pattern found" -ForegroundColor Yellow
        }
    }

    # Intel OpenMP and compiler runtime DLLs
    $INTEL_COMPILER = "$INTEL_ONEAPI\compiler\latest"
    $RUNTIME_DLL_PATTERNS = @(
        "libiomp5md.dll",     # Intel OpenMP runtime
        "libmmd.dll",         # Intel math library
        "svml_dispmd.dll"     # Intel short vector math library
    )

    Write-Host "Copying Intel compiler runtime DLLs..." -ForegroundColor Cyan
    foreach ($pattern in $RUNTIME_DLL_PATTERNS) {
        $dlls = Get-ChildItem -Path "$INTEL_COMPILER\bin" -Filter $pattern -ErrorAction SilentlyContinue
        if ($dlls) {
            foreach ($dll in $dlls) {
                Copy-Item $dll.FullName "$PROJECT_DIR\src\radia\" -Force
                Write-Host "  Copied $($dll.Name)" -ForegroundColor Green
            }
        } else {
            Write-Host "  WARNING: No DLL matching $pattern found" -ForegroundColor Yellow
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
    Write-Host ""
    Write-Host "========================================" -ForegroundColor Cyan
    Write-Host "  Running Tests" -ForegroundColor Cyan
    Write-Host "========================================" -ForegroundColor Cyan
    Write-Host ""

    # First, test basic import
    Write-Host "Testing radia import..." -ForegroundColor Cyan
    python -c "import sys; sys.path.insert(0, r'$PROJECT_DIR\src\radia'); import radia; print(f'radia version: {radia.__version__}')"

    if ($LASTEXITCODE -eq 0) {
        Write-Host "radia import test passed!" -ForegroundColor Green
    } else {
        Write-Host "radia import test failed!" -ForegroundColor Red
        exit 1
    }

    Write-Host "Testing radia_ngsolve import..." -ForegroundColor Cyan
    python -c "import sys; sys.path.insert(0, r'$PROJECT_DIR\src\radia'); import radia_ngsolve; print('radia_ngsolve import OK')"

    if ($LASTEXITCODE -eq 0) {
        Write-Host "radia_ngsolve import test passed!" -ForegroundColor Green
    } else {
        Write-Host "radia_ngsolve import test failed!" -ForegroundColor Yellow
        Write-Host "(NGSolve may not be installed - continuing)" -ForegroundColor Gray
    }

    Write-Host ""

    # Run pytest on basic tests (fast tests only by default)
    Write-Host "Running pytest (basic tests)..." -ForegroundColor Cyan
    Write-Host ""

    # Check if pytest is available
    $pytestCheck = python -c "import pytest; print('ok')" 2>&1
    if ($LASTEXITCODE -ne 0) {
        Write-Host "pytest not installed. Install with: pip install pytest" -ForegroundColor Yellow
        Write-Host "Skipping pytest..." -ForegroundColor Gray
    } else {
        # Run pytest with basic marker (fast tests only)
        # Use -m "basic" to run only basic tests, or -m "not slow and not benchmark" for quick tests
        Push-Location $PROJECT_DIR
        try {
            python -m pytest tests/ -m "basic" -v --tb=short
            $TestResult = $LASTEXITCODE

            if ($TestResult -eq 0) {
                Write-Host ""
                Write-Host "========================================" -ForegroundColor Green
                Write-Host "  All Tests Passed!" -ForegroundColor Green
                Write-Host "========================================" -ForegroundColor Green
            } elseif ($TestResult -eq 5) {
                # Exit code 5 means no tests were collected (no tests matched the marker)
                Write-Host ""
                Write-Host "No tests with 'basic' marker found. Running quick tests..." -ForegroundColor Yellow
                python -m pytest tests/test_radia.py -v --tb=short
                $TestResult = $LASTEXITCODE
                if ($TestResult -eq 0) {
                    Write-Host ""
                    Write-Host "========================================" -ForegroundColor Green
                    Write-Host "  All Tests Passed!" -ForegroundColor Green
                    Write-Host "========================================" -ForegroundColor Green
                } else {
                    Write-Host ""
                    Write-Host "========================================" -ForegroundColor Red
                    Write-Host "  Some Tests Failed!" -ForegroundColor Red
                    Write-Host "========================================" -ForegroundColor Red
                    exit $TestResult
                }
            } else {
                Write-Host ""
                Write-Host "========================================" -ForegroundColor Red
                Write-Host "  Some Tests Failed!" -ForegroundColor Red
                Write-Host "========================================" -ForegroundColor Red
                exit $TestResult
            }
        }
        finally {
            Pop-Location
        }
    }
    Write-Host ""
}

# Show usage
Write-Host "Usage:" -ForegroundColor Yellow
Write-Host "  import radia" -ForegroundColor Gray
Write-Host "  import radia_ngsolve" -ForegroundColor Gray
Write-Host "  print(radia.__version__)" -ForegroundColor Gray
Write-Host ""
