#==============================================================================
# Build.ps1 - Build Radia with MSVC + Intel MKL
#
# This script builds both radia.pyd and radia_ngsolve.pyd using MSVC compiler
# with Intel MKL for BLAS/LAPACK operations.
#
# Usage:
#   .\Build.ps1
#   .\Build.ps1 -Rebuild
#   .\Build.ps1 -Test
#   .\Build.ps1 -RadiaOnly
#
# Options:
#   -Rebuild    Clean build directory before building
#   -Test       Run import test after build
#   -NoOpenMP   Disable OpenMP (for debugging)
#   -NoExaFMM   Disable ExaFMM (for debugging, ExaFMM is ON by default)
#   -RadiaOnly  Build only radia.pyd (skip radia_ngsolve)
#   -Verbose    Show detailed build output
#==============================================================================

param(
    [switch]$Rebuild,
    [switch]$Test,
    [switch]$NoOpenMP,     # Disable OpenMP (for debugging)
    [switch]$NoExaFMM,     # Disable ExaFMM (for debugging, ExaFMM is enabled by default)
    [switch]$RadiaOnly,    # Build only radia.pyd
    [switch]$Verbose       # Show detailed build output
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

# ExaFMM flag (enabled by default for FMM-accelerated conductor field computation)
if ($NoExaFMM) {
    $EXAFMM_FLAG = "OFF"
    Write-Host "ExaFMM: DISABLED (debug mode)" -ForegroundColor Yellow
} else {
    $EXAFMM_FLAG = "ON"
    Write-Host "ExaFMM: ENABLED" -ForegroundColor Green
}

# Determine build targets and NGSolve flag
if ($RadiaOnly) {
    $BUILD_TARGETS = "_radia"
    $NGSOLVE_FLAG = "OFF"
    Write-Host "Build target: _radia only (no NGSolve)" -ForegroundColor Gray
} else {
    $BUILD_TARGETS = "_radia radia_ngsolve"
    $NGSOLVE_FLAG = "ON"
    Write-Host "Build targets: _radia, radia_ngsolve" -ForegroundColor Gray
}

# Create batch file to run with Visual Studio environment
$BatchContent = @"
@echo off
setlocal enabledelayedexpansion

REM Set up Visual Studio environment
echo Setting up Visual Studio environment...
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
if errorlevel 1 (
    echo ERROR: Failed to set up Visual Studio environment
    exit /b 1
)

REM Add MKL paths
set LIB=$INTEL_MKL\lib;%LIB%
set INCLUDE=$INTEL_MKL\include;%INCLUDE%
set MKLROOT=$INTEL_MKL

cd /d "$BUILD_DIR"

echo.
echo ========================================
echo   Configuring CMake with MSVC...
echo ========================================
echo.
"$CMAKE_EXE" "$PROJECT_DIR" ^
    -G "Ninja" ^
    -DCMAKE_C_COMPILER=cl ^
    -DCMAKE_CXX_COMPILER=cl ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DRADIA_ENABLE_OPENMP=$OPENMP_FLAG ^
    -DRADIA_ENABLE_EXAFMM=$EXAFMM_FLAG ^
    -DRADIA_BUILD_NGSOLVE=$NGSOLVE_FLAG

if errorlevel 1 (
    echo ERROR: CMake configuration failed
    exit /b 1
)

REM Legacy _radia module is no longer built
REM pybind11 module (_radia_pybind) is now the default
REM See __init__.py which imports from _radia_pybind

"@

# Add radia_ngsolve build if not RadiaOnly
if (-not $RadiaOnly) {
    $BatchContent += @"

echo.
echo ========================================
echo   Building radia_ngsolve...
echo ========================================
echo.
"$CMAKE_EXE" --build . --config Release --target radia_ngsolve -j

if errorlevel 1 (
    echo WARNING: radia_ngsolve build failed (NGSolve may not be installed)
    REM Continue anyway - radia.pyd is the main target
)

echo.
echo ========================================
echo   Building peec_matrices...
echo ========================================
echo.
"$CMAKE_EXE" --build . --config Release --target peec_matrices -j

if errorlevel 1 (
    echo WARNING: peec_matrices build failed
)

echo.
echo ========================================
echo   Building cln_core...
echo ========================================
echo.
"$CMAKE_EXE" --build . --config Release --target cln_core -j

if errorlevel 1 (
    echo WARNING: cln_core build failed
)

echo.
echo ========================================
echo   Building mmm_core...
echo ========================================
echo.
"$CMAKE_EXE" --build . --config Release --target mmm_core -j

if errorlevel 1 (
    echo WARNING: mmm_core build failed
)

echo.
echo ========================================
echo   Building _radia_pybind...
echo ========================================
echo.
"$CMAKE_EXE" --build . --config Release --target _radia_pybind -j

if errorlevel 1 (
    echo ERROR: _radia_pybind build failed
    exit /b 1
)

"@
} else {
    # RadiaOnly mode: Just build _radia_pybind
    $BatchContent += @"

echo.
echo ========================================
echo   Building _radia_pybind (RadiaOnly)...
echo ========================================
echo.
"$CMAKE_EXE" --build . --config Release --target _radia_pybind -j

if errorlevel 1 (
    echo ERROR: _radia_pybind build failed
    exit /b 1
)

"@
}

$BatchContent += @"

echo.
echo Build completed successfully.
exit /b 0
"@

$BatchFile = "$PROJECT_DIR\build_temp_msvc.bat"
$BatchContent | Out-File -FilePath $BatchFile -Encoding ascii

# Build log file
$BuildLog = "$PROJECT_DIR\build_log.txt"

try {
    # Run the batch file with output capture
    Write-Host ""
    Write-Host "Building..." -ForegroundColor Cyan
    Write-Host ""

    if ($Verbose) {
        # Verbose mode: show output in real-time
        $process = Start-Process -FilePath "cmd.exe" -ArgumentList "/c", $BatchFile -NoNewWindow -PassThru -Wait
        $BuildResult = $process.ExitCode
    } else {
        # Normal mode: capture output to log file and show summary
        $process = Start-Process -FilePath "cmd.exe" -ArgumentList "/c", "$BatchFile > `"$BuildLog`" 2>&1" -NoNewWindow -PassThru -Wait
        $BuildResult = $process.ExitCode

        # Show last 50 lines of build log
        if (Test-Path $BuildLog) {
            $logContent = Get-Content $BuildLog -Tail 50
            foreach ($line in $logContent) {
                if ($line -match "error|ERROR") {
                    Write-Host $line -ForegroundColor Red
                } elseif ($line -match "warning|WARNING") {
                    Write-Host $line -ForegroundColor Yellow
                } elseif ($line -match "Building|Linking|Compiling") {
                    Write-Host $line -ForegroundColor Cyan
                } else {
                    Write-Host $line
                }
            }
        }
    }

    if ($BuildResult -ne 0) {
        throw "Build failed with exit code $BuildResult"
    }

    # Legacy _radia.pyd is no longer built or copied
    # pybind11 module (_radia_pybind.pyd) is now the default
    # __init__.py imports from _radia_pybind

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

    # Copy peec_matrices.pyd to package directory (pybind11 module)
    $PEEC_PYD_SOURCE = "$BUILD_DIR\peec_matrices.cp312-win_amd64.pyd"
    $PEEC_PYD_DEST = "$PROJECT_DIR\src\radia\peec_matrices.pyd"

    if (Test-Path $PEEC_PYD_SOURCE) {
        Copy-Item $PEEC_PYD_SOURCE $PEEC_PYD_DEST -Force
        Write-Host "Copied peec_matrices.pyd to src/radia/" -ForegroundColor Green
    } else {
        Write-Host "WARNING: Could not find $PEEC_PYD_SOURCE" -ForegroundColor Yellow
    }

    # Copy cln_core.pyd to package directory (pybind11 module)
    $CLN_PYD_SOURCE = "$BUILD_DIR\cln_core.cp312-win_amd64.pyd"
    $CLN_PYD_DEST = "$PROJECT_DIR\src\radia\cln_core.pyd"

    if (Test-Path $CLN_PYD_SOURCE) {
        Copy-Item $CLN_PYD_SOURCE $CLN_PYD_DEST -Force
        Write-Host "Copied cln_core.pyd to src/radia/" -ForegroundColor Green
    } else {
        Write-Host "WARNING: Could not find $CLN_PYD_SOURCE" -ForegroundColor Yellow
    }

    # Copy mmm_core.pyd to package directory (pybind11 module)
    $MMM_PYD_SOURCE = "$BUILD_DIR\mmm_core.cp312-win_amd64.pyd"
    $MMM_PYD_DEST = "$PROJECT_DIR\src\radia\mmm_core.pyd"

    if (Test-Path $MMM_PYD_SOURCE) {
        Copy-Item $MMM_PYD_SOURCE $MMM_PYD_DEST -Force
        Write-Host "Copied mmm_core.pyd to src/radia/" -ForegroundColor Green
    } else {
        Write-Host "WARNING: Could not find $MMM_PYD_SOURCE" -ForegroundColor Yellow
    }

    # Copy _radia_pybind.pyd to package directory (main Python module)
    # This is the primary pybind11 module that __init__.py imports
    $RADIA_PYBIND_SOURCE = "$BUILD_DIR\_radia_pybind.cp312-win_amd64.pyd"
    $RADIA_PYBIND_DEST = "$PROJECT_DIR\src\radia\_radia_pybind.pyd"

    if (Test-Path $RADIA_PYBIND_SOURCE) {
        Copy-Item $RADIA_PYBIND_SOURCE $RADIA_PYBIND_DEST -Force
        Write-Host "Copied _radia_pybind.pyd to src/radia/" -ForegroundColor Green
    } else {
        Write-Host "ERROR: _radia_pybind.pyd not found - build failed?" -ForegroundColor Red
    }

    # MKL DLLs are NOT bundled in the wheel (v2.0.0+)
    # MKL is installed via pip dependency: mkl>=2024.2.0, intel-cmplr-lib-rt>=2024.2.0
    # DLLs are found at runtime via {sys.prefix}/Library/bin/ (see __init__.py)
    Write-Host "MKL DLLs: NOT bundled (pip dependency)" -ForegroundColor Cyan

    # Show build summary
    Write-Host ""
    Write-Host "========================================" -ForegroundColor Green
    Write-Host "  Build Completed Successfully" -ForegroundColor Green
    Write-Host "========================================" -ForegroundColor Green
    Write-Host ""

    # Show PYD file info
    Write-Host "Output:" -ForegroundColor Yellow
    if (Test-Path $RADIA_PYBIND_DEST) {
        $pybindInfo = Get-Item $RADIA_PYBIND_DEST
        Write-Host "  _radia_pybind.pyd: $([math]::Round($pybindInfo.Length / 1MB, 2)) MB" -ForegroundColor Gray
        Write-Host "  Modified:  $($pybindInfo.LastWriteTime)" -ForegroundColor Gray
    }
    if (Test-Path $NGSOLVE_PYD_DEST) {
        $ngInfo = Get-Item $NGSOLVE_PYD_DEST
        Write-Host "  radia_ngsolve.pyd: $([math]::Round($ngInfo.Length / 1MB, 2)) MB" -ForegroundColor Gray
    }
    if (Test-Path $PEEC_PYD_DEST) {
        $peecInfo = Get-Item $PEEC_PYD_DEST
        Write-Host "  peec_matrices.pyd: $([math]::Round($peecInfo.Length / 1MB, 2)) MB" -ForegroundColor Gray
    }
    if (Test-Path $CLN_PYD_DEST) {
        $clnInfo = Get-Item $CLN_PYD_DEST
        Write-Host "  cln_core.pyd: $([math]::Round($clnInfo.Length / 1MB, 2)) MB" -ForegroundColor Gray
    }
    if (Test-Path $MMM_PYD_DEST) {
        $mmmInfo = Get-Item $MMM_PYD_DEST
        Write-Host "  mmm_core.pyd: $([math]::Round($mmmInfo.Length / 1MB, 2)) MB" -ForegroundColor Gray
    }
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

    # Show error details from build log
    if (Test-Path $BuildLog) {
        Write-Host "Last 30 lines of build log:" -ForegroundColor Yellow
        Get-Content $BuildLog -Tail 30 | ForEach-Object {
            if ($_ -match "error|ERROR") {
                Write-Host $_ -ForegroundColor Red
            } else {
                Write-Host $_
            }
        }
        Write-Host ""
        Write-Host "Full log: $BuildLog" -ForegroundColor Gray
    }
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
    python -c "import pytest" 2>&1 | Out-Null
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
