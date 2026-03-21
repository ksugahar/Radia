#==============================================================================
# Build.ps1 - Build Radia with MSVC + Intel MKL + NGSolve
#
# Builds _radia_pybind.pyd and auxiliary modules using MSVC compiler
# with Intel MKL for BLAS/LAPACK and NGSolve for TaskManager parallelization.
#
# Usage:
#   powershell.exe -ExecutionPolicy Bypass -File Build.ps1
#   powershell.exe -ExecutionPolicy Bypass -File Build.ps1 -Rebuild
#   powershell.exe -ExecutionPolicy Bypass -File Build.ps1 -Test
#
# Options:
#   -Rebuild    Clean build directory before building
#   -Test       Run import test + pytest after build
#   -Verbose    Show detailed build output
#
# Requirements:
#   - Visual Studio 2022 (MSVC compiler)
#   - Intel oneAPI Base Toolkit (MKL only, NOT the compiler)
#   - NGSolve (pip install or source build)
#   - Python 3.12 with pybind11
#==============================================================================

param(
    [switch]$Rebuild,
    [switch]$Test,
    [switch]$Verbose,
    [switch]$RadiaOnly
)

$ErrorActionPreference = "Stop"

# ============================================================================
# Paths
# ============================================================================

$PROJECT_DIR = Split-Path -Parent $MyInvocation.MyCommand.Path
$BUILD_DIR = "$PROJECT_DIR\build-msvc"

# Intel MKL (required for BLAS/LAPACK)
if ($env:MKLROOT -and (Test-Path $env:MKLROOT)) {
    $INTEL_MKL = $env:MKLROOT
} else {
    $INTEL_MKL = "C:\Program Files (x86)\Intel\oneAPI\mkl\latest"
}
if (-not (Test-Path "$INTEL_MKL\lib\mkl_rt.lib")) {
    Write-Host "ERROR: Intel MKL not found at $INTEL_MKL" -ForegroundColor Red
    Write-Host "Install Intel oneAPI Base Toolkit (MKL component)" -ForegroundColor Yellow
    exit 1
}

# NGSolve (optional override via NGSOLVE_DIR environment variable)
$NGSolveCMakeArgs = ""
if ($env:NGSOLVE_DIR -and (Test-Path "$env:NGSOLVE_DIR\NGSolveConfig.cmake")) {
    $NGSolveCMakeArgs = " ^`n    -DNGSolve_DIR=`"$env:NGSOLVE_DIR`" ^`n    -DNetgen_DIR=`"$env:NGSOLVE_DIR`""
    Write-Host "NGSolve: $env:NGSOLVE_DIR (from env)" -ForegroundColor Gray
}

# CMake (from Visual Studio 2022)
$CMAKE_EXE = $null
$VSWHERE = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
if (Test-Path $VSWHERE) {
    $VS_PATH = & $VSWHERE -latest -property installationPath 2>$null
    if ($VS_PATH) {
        $CMAKE_CANDIDATE = "$VS_PATH\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
        if (Test-Path $CMAKE_CANDIDATE) {
            $CMAKE_EXE = $CMAKE_CANDIDATE
        }
    }
}
if (-not $CMAKE_EXE) {
    # Try BuildTools edition
    $CMAKE_EXE = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\2022\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
}
if (-not (Test-Path $CMAKE_EXE)) {
    # Try Community edition
    $CMAKE_EXE = "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
}
if (-not (Test-Path $CMAKE_EXE)) {
    # Fallback: pip-installed cmake
    $PIP_CMAKE = & python -c "import shutil; print(shutil.which('cmake') or '')" 2>$null
    if ($PIP_CMAKE -and (Test-Path $PIP_CMAKE)) {
        $CMAKE_EXE = $PIP_CMAKE
    }
}
if (-not (Test-Path $CMAKE_EXE)) {
    Write-Host "ERROR: CMake not found. Install Visual Studio 2022 with CMake or pip install cmake." -ForegroundColor Red
    exit 1
}

# ============================================================================
# Setup
# ============================================================================

Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "  Radia Build (MSVC + MKL + NGSolve)" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "MKL:   $INTEL_MKL" -ForegroundColor Gray
Write-Host "Build: $BUILD_DIR" -ForegroundColor Gray
Write-Host ""

if ($Rebuild -and (Test-Path $BUILD_DIR)) {
    Write-Host "Cleaning build directory..." -ForegroundColor Yellow
    Remove-Item -Recurse -Force $BUILD_DIR
}
if (-not (Test-Path $BUILD_DIR)) {
    New-Item -ItemType Directory -Path $BUILD_DIR | Out-Null
}

# ============================================================================
# Build (via batch file for vcvars64 environment)
# ============================================================================

$BatchContent = @"
@echo off
setlocal enabledelayedexpansion

call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" > nul 2>&1

set LIB=$INTEL_MKL\lib;%LIB%
set INCLUDE=$INTEL_MKL\include;%INCLUDE%
set MKLROOT=$INTEL_MKL

cd /d "$BUILD_DIR"

echo.
echo ========================================
echo   CMake Configure
echo ========================================
"$CMAKE_EXE" "$PROJECT_DIR" ^
    -G "Ninja" ^
    -DCMAKE_C_COMPILER=cl ^
    -DCMAKE_CXX_COMPILER=cl ^
    -DCMAKE_BUILD_TYPE=Release$NGSolveCMakeArgs

if errorlevel 1 (
    echo ERROR: CMake configuration failed
    exit /b 1
)

echo.
echo ========================================
echo   Building _radia_pybind
echo ========================================
"$CMAKE_EXE" --build . --config Release --target _radia_pybind -j
if errorlevel 1 (
    echo ERROR: _radia_pybind build failed
    exit /b 1
)

echo.
echo ========================================
echo   Building peec_matrices
echo ========================================
"$CMAKE_EXE" --build . --config Release --target peec_matrices -j
if errorlevel 1 ( echo WARNING: peec_matrices build failed )

echo.
echo ========================================
echo   Building cln_core
echo ========================================
"$CMAKE_EXE" --build . --config Release --target cln_core -j
if errorlevel 1 ( echo WARNING: cln_core build failed )

echo.
echo ========================================
echo   Building mmm_core
echo ========================================
"$CMAKE_EXE" --build . --config Release --target mmm_core -j
if errorlevel 1 ( echo WARNING: mmm_core build failed )

echo.
echo Build completed.
exit /b 0
"@

$BatchFile = "$PROJECT_DIR\build_temp_msvc.bat"
$BuildLog = "$PROJECT_DIR\build_log.txt"
$BatchContent | Out-File -FilePath $BatchFile -Encoding ascii

try {
    Write-Host "Building..." -ForegroundColor Cyan

    if ($Verbose) {
        $process = Start-Process -FilePath "cmd.exe" -ArgumentList "/c", $BatchFile -NoNewWindow -PassThru -Wait
        $BuildResult = $process.ExitCode
    } else {
        $process = Start-Process -FilePath "cmd.exe" -ArgumentList "/c", "$BatchFile > `"$BuildLog`" 2>&1" -NoNewWindow -PassThru -Wait
        $BuildResult = $process.ExitCode
        if (Test-Path $BuildLog) {
            Get-Content $BuildLog -Tail 30 | ForEach-Object {
                if ($_ -match "error|ERROR") { Write-Host $_ -ForegroundColor Red }
                elseif ($_ -match "warning|WARNING") { Write-Host $_ -ForegroundColor Yellow }
                else { Write-Host $_ }
            }
        }
    }

    if ($BuildResult -ne 0) { throw "Build failed with exit code $BuildResult" }

    # ========================================================================
    # Copy .pyd files to src/radia/
    # ========================================================================

    # Delete stale radia_ngsolve.pyd (merged into _radia_pybind, now pure Python)
    $STALE_NGSOLVE_PYD = "$PROJECT_DIR\src\radia\radia_ngsolve.pyd"
    if (Test-Path $STALE_NGSOLVE_PYD) {
        Remove-Item $STALE_NGSOLVE_PYD -Force
        Write-Host "  Removed stale radia_ngsolve.pyd" -ForegroundColor Yellow
    }

    $modules = @(
        @{ src = "_radia_pybind.cp312-win_amd64.pyd"; dst = "_radia_pybind.pyd"; required = $true },
        @{ src = "peec_matrices.cp312-win_amd64.pyd"; dst = "peec_matrices.pyd"; required = $false },
        @{ src = "cln_core.cp312-win_amd64.pyd";      dst = "cln_core.pyd";      required = $false },
        @{ src = "mmm_core.cp312-win_amd64.pyd";       dst = "mmm_core.pyd";      required = $false }
    )

    foreach ($mod in $modules) {
        $srcPath = "$BUILD_DIR\$($mod.src)"
        $dstPath = "$PROJECT_DIR\src\radia\$($mod.dst)"
        if (Test-Path $srcPath) {
            Copy-Item $srcPath $dstPath -Force
            $info = Get-Item $dstPath
            Write-Host "  $($mod.dst): $([math]::Round($info.Length / 1MB, 2)) MB" -ForegroundColor Green
        } elseif ($mod.required) {
            Write-Host "ERROR: $($mod.src) not found" -ForegroundColor Red
            exit 1
        } else {
            Write-Host "  $($mod.dst): skipped" -ForegroundColor Yellow
        }
    }

    Write-Host ""
    Write-Host "========================================" -ForegroundColor Green
    Write-Host "  Build Completed" -ForegroundColor Green
    Write-Host "========================================" -ForegroundColor Green
}
catch {
    Write-Host ""
    Write-Host "========================================" -ForegroundColor Red
    Write-Host "  Build Failed: $_" -ForegroundColor Red
    Write-Host "========================================" -ForegroundColor Red
    if (Test-Path $BuildLog) { Write-Host "Log: $BuildLog" -ForegroundColor Gray }
    exit 1
}
finally {
    if (Test-Path $BatchFile) { Remove-Item $BatchFile -Force }
}

# ============================================================================
# Test (optional)
# ============================================================================

if ($Test) {
    Write-Host ""
    Write-Host "Running tests..." -ForegroundColor Cyan

    python -c "import radia; print(f'radia {radia.__version__} OK')"
    if ($LASTEXITCODE -ne 0) { Write-Host "Import failed!" -ForegroundColor Red; exit 1 }

    Push-Location $PROJECT_DIR
    try { python -m pytest tests/ -v --tb=short }
    finally { Pop-Location }
}
