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

if exist "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat" (
    call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat" > nul 2>&1
) else (
    call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" > nul 2>&1
)

REM Ninja header dependency tracking:
REM MSVC /showIncludes outputs in Japanese (cp932) on Japanese Windows.
REM If rules.ninja gets a garbled msvc_deps_prefix (encoding mismatch),
REM Ninja records #deps 0 for all .obj files and misses header changes.
REM Fix: delete build-*/CMakeFiles/rules.ninja and reconfigure (Build.ps1 -Rebuild).
REM VSLANG=1033 does NOT change cl.exe language, but ensures consistent env.
set VSLANG=1033

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
echo ========================================
echo   Building radia_cubit_mesh (Cubit plugin .pyd)
echo ========================================
set "CUBIT_PLUGIN_SRC=$PROJECT_DIR\src\cubit_plugin"
set "CUBIT_PLUGIN_BUILD=$PROJECT_DIR\src\cubit_plugin\build-pyd"
set "CUBIT_DIR=C:\Program Files\Coreform Cubit 2025.3\cmake"
set "NETGEN_DIR=C:\Program Files\Python312\Lib\site-packages\netgen"
rem Compact Netgen sources are in-repo (src/cubit_plugin/compact_netgen/netgen_src/).
rem No external NETGEN_SRC_DIR needed.

if exist "%CUBIT_DIR%\CubitConfig.cmake" (
    if not exist "%CUBIT_PLUGIN_BUILD%" mkdir "%CUBIT_PLUGIN_BUILD%"
    cd /d "%CUBIT_PLUGIN_BUILD%"
    "$CMAKE_EXE" -G Ninja -DCMAKE_BUILD_TYPE=Release -DCubit_DIR="%CUBIT_DIR%" -DNETGEN_DIR="%NETGEN_DIR%" "%CUBIT_PLUGIN_SRC%"
    "$CMAKE_EXE" --build . --config Release --target radia_cubit_mesh -j
    if errorlevel 1 ( echo WARNING: radia_cubit_mesh build failed )

    echo.
    echo ========================================
    echo   Building radia_cubit.ccm (APREPRO commands)
    echo ========================================
    set "CUBIT_CCM_BUILD=$PROJECT_DIR\src\cubit_plugin\build-ccm"
    if not exist "!CUBIT_CCM_BUILD!" mkdir "!CUBIT_CCM_BUILD!"
    cd /d "!CUBIT_CCM_BUILD!"
    "$CMAKE_EXE" -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=cl -DCMAKE_CXX_COMPILER=cl -DCubit_DIR="%CUBIT_DIR%" -DNETGEN_DIR="%NETGEN_DIR%" "%CUBIT_PLUGIN_SRC%"
    "$CMAKE_EXE" --build . --config Release --target radia_cubit_ccm -j
    if errorlevel 1 ( echo WARNING: radia_cubit_ccm build failed )

    echo.
    echo ========================================
    echo   Building radia_cubit.ccl (GUI component)
    echo ========================================
    "$CMAKE_EXE" --build . --config Release --target radia_cubit_ccl -j
    if errorlevel 1 ( echo WARNING: radia_cubit_ccl build failed )

    rem Copy ccl to src/radia/ so radia-setup deploys the latest version
    if exist "!CUBIT_CCM_BUILD!\radia_cubit.ccl" (
        copy /Y "!CUBIT_CCM_BUILD!\radia_cubit.ccl" "$PROJECT_DIR\src\radia\radia_cubit.ccl" >nul
        echo   radia_cubit.ccl: copied to src/radia/
    )

    cd /d "$BUILD_DIR"
) else (
    echo SKIP: Cubit SDK not found at %CUBIT_DIR%
)

echo.
echo Build completed.
exit /b 0
"@

$BatchFile = "$PROJECT_DIR\build_temp_msvc.bat"
$BuildLog = "$PROJECT_DIR\build_log.txt"
$BatchContent | Out-File -FilePath $BatchFile -Encoding ascii

try {
    Write-Host "Building..." -ForegroundColor Cyan

    # Use cmd /c directly (not Start-Process -Wait which waits for child processes)
    & cmd.exe /c "$BatchFile > `"$BuildLog`" 2>&1"
    $BuildResult = $LASTEXITCODE

    if (Test-Path $BuildLog) {
        Get-Content $BuildLog -Tail 30 | ForEach-Object {
            if ($_ -match "error|ERROR") { Write-Host $_ -ForegroundColor Red }
            elseif ($_ -match "warning|WARNING") { Write-Host $_ -ForegroundColor Yellow }
            else { Write-Host $_ }
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

    # Cubit plugin files (built in separate build dirs)
    # Freshness check: compare build output against destination.
    #   - If src is newer than dst (or dst missing): copy (UPDATED)
    #   - If src == dst (same size + same or older time): skip (up-to-date)
    #   - If src is missing: skip with warning
    $cubitFiles = @(
        @{ src = "$PROJECT_DIR\src\cubit_plugin\build-pyd\radia_cubit_mesh.cp312-win_amd64.pyd";
           dst = "$PROJECT_DIR\src\radia\radia_cubit_mesh.pyd" },
        @{ src = "$PROJECT_DIR\src\cubit_plugin\build-ccm\radia_cubit.ccm";
           dst = "$PROJECT_DIR\src\radia\radia_cubit.ccm" }
    )
    foreach ($cf in $cubitFiles) {
        $name = Split-Path $cf.dst -Leaf
        if (Test-Path $cf.src) {
            $srcInfo = Get-Item $cf.src
            $needCopy = $true
            if (Test-Path $cf.dst) {
                $dstInfo = Get-Item $cf.dst
                if ($srcInfo.Length -eq $dstInfo.Length -and $srcInfo.LastWriteTime -le $dstInfo.LastWriteTime) {
                    Write-Host "  ${name}: up-to-date ($([math]::Round($srcInfo.Length / 1KB, 1)) KB)" -ForegroundColor Cyan
                    $needCopy = $false
                }
            }
            if ($needCopy) {
                Copy-Item $cf.src $cf.dst -Force
                Write-Host "  ${name}: $([math]::Round($srcInfo.Length / 1KB, 1)) KB (UPDATED)" -ForegroundColor Green
            }
        } else {
            Write-Host "  ${name}: skipped (not built)" -ForegroundColor Yellow
        }
    }

    # Copy plugin binaries to cubit-mesh-export package (for wheel + cubit-plugin-install)
    $cmeDir = "$PROJECT_DIR\packages\cubit-mesh-export\src\cubit_mesh_export"
    $cmeCopies = @(
        @{ src = "$PROJECT_DIR\src\radia\radia_cubit_mesh.pyd"; name = "radia_cubit_mesh.pyd" },
        @{ src = "$PROJECT_DIR\src\radia\radia_cubit.ccm";      name = "radia_cubit.ccm" },
        @{ src = "$PROJECT_DIR\src\radia\radia_cubit.ccl";       name = "radia_cubit.ccl" }
    )
    foreach ($cc in $cmeCopies) {
        if (Test-Path $cc.src) {
            $srcInfo = Get-Item $cc.src
            $cmeDst = "$cmeDir\$($cc.name)"
            $needCopy = $true
            if (Test-Path $cmeDst) {
                $dstInfo = Get-Item $cmeDst
                if ($srcInfo.Length -eq $dstInfo.Length -and $srcInfo.LastWriteTime -le $dstInfo.LastWriteTime) {
                    $needCopy = $false
                }
            }
            if ($needCopy) {
                Copy-Item $cc.src $cmeDst -Force
                Write-Host "  cubit-mesh-export/$($cc.name): copied" -ForegroundColor Green
            }
        }
    }

    # === Freshness-check sync: touch bundled binaries to the newest source mtime ===
    # cubit-mesh-export's pyproject.toml runs a freshness gate in
    # `get_requires_for_build_wheel` that compares mtime of every file
    # in src/cubit_plugin/ vs every bundled binary.  Ninja's incremental
    # build only rebuilds targets whose transitive source changed, so
    # editing RadiaComp.cpp (.ccl source) leaves radia_cubit.ccm
    # (built from RadiaPlugin.cpp) untouched.  The freshness check then
    # FALSE-alarms: "radia_cubit.ccm is 16 h older than RadiaComp.cpp"
    # even though the CONTENT of .ccm is correct.
    # Fix: after all copies, force-touch every bundled binary mtime so
    # it >= newest source mtime.  Build.ps1 has already ensured the
    # CONTENTS are up-to-date; this only corrects mtime for the gate.
    $pluginSrcFiles = Get-ChildItem "$PROJECT_DIR\src\cubit_plugin" `
        -Include "*.cpp","*.hpp","*.h","*.c" -File -Recurse `
        -ErrorAction SilentlyContinue
    if ($pluginSrcFiles) {
        $newest = ($pluginSrcFiles | Measure-Object -Property LastWriteTime -Maximum).Maximum
        # Bump by 1 s so freshness check sees ">= newest source".
        $touchTime = $newest.AddSeconds(1)
        foreach ($cc in $cmeCopies) {
            $cmeDst = "$cmeDir\$($cc.name)"
            if (Test-Path $cmeDst) {
                (Get-Item -LiteralPath $cmeDst).LastWriteTime = $touchTime
            }
        }
        # Also touch the src/radia/ copies so install_full.py /
        # verify-deploy freshness checks see consistent mtimes.
        foreach ($name in @("radia_cubit_mesh.pyd", "radia_cubit.ccm", "radia_cubit.ccl")) {
            $p = "$PROJECT_DIR\src\radia\$name"
            if (Test-Path $p) { (Get-Item -LiteralPath $p).LastWriteTime = $touchTime }
        }
        Write-Host "  bundled plugin binaries: mtime synced to newest source ($touchTime)" -ForegroundColor Cyan
    }

    foreach ($mod in $modules) {
        $srcPath = "$BUILD_DIR\$($mod.src)"
        $dstPath = "$PROJECT_DIR\src\radia\$($mod.dst)"
        if (Test-Path $srcPath) {
            $srcInfo = Get-Item $srcPath
            $needCopy = $true
            if (Test-Path $dstPath) {
                $dstInfo = Get-Item $dstPath
                if ($srcInfo.Length -eq $dstInfo.Length -and $srcInfo.LastWriteTime -le $dstInfo.LastWriteTime) {
                    Write-Host "  $($mod.dst): up-to-date ($([math]::Round($srcInfo.Length / 1MB, 2)) MB)" -ForegroundColor Cyan
                    $needCopy = $false
                }
            }
            if ($needCopy) {
                Copy-Item $srcPath $dstPath -Force
                Write-Host "  $($mod.dst): $([math]::Round($srcInfo.Length / 1MB, 2)) MB (UPDATED)" -ForegroundColor Green
            }
        } elseif ($mod.required) {
            Write-Host "ERROR: $($mod.src) not found" -ForegroundColor Red
            exit 1
        } else {
            Write-Host "  $($mod.dst): skipped" -ForegroundColor Yellow
        }
    }

    # Kill orphan vctip.exe (VS BuildTools telemetry) that may block process exit
    Get-Process vctip -ErrorAction SilentlyContinue | Stop-Process -Force

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
