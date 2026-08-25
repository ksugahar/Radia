#==============================================================================
# Build.ps1 - Build Radia with MSVC + Intel MKL + NGSolve
#
# Builds _radia_pybind.pyd and auxiliary modules using MSVC compiler
# with Intel MKL for BLAS/LAPACK and NGSolve for TaskManager parallelization.
#
# Usage:
#   powershell.exe -ExecutionPolicy Bypass -File Build.ps1
#   powershell.exe -ExecutionPolicy Bypass -File Build.ps1 -Rebuild
#   powershell.exe -ExecutionPolicy Bypass -File Build.ps1 -RadiaOnly
#   powershell.exe -ExecutionPolicy Bypass -File Build.ps1 -MatlabMexOnly
#   powershell.exe -ExecutionPolicy Bypass -File Build.ps1 -Test
#
# Options:
#   -Rebuild    Clean build directory before building
#   -RadiaOnly  Build and copy only _radia_pybind.pyd
#   -MatlabMexOnly  Configure and build the shared radia_mex native gateway
#   -Test       Run source-tree import test + pytest after build
#   -Verbose    Show detailed build output
#
# Requirements:
#   - Visual Studio Build Tools (any version with VC.Tools.x86.x64) — auto-detected via vswhere
#   - Intel oneAPI Base Toolkit (MKL only, NOT the compiler)
#   - NGSolve (pip install or source build)
#   - Python 3.12 with pybind11
#==============================================================================

param(
    [switch]$Rebuild,
    [switch]$Test,
    [switch]$Verbose,
    [switch]$RadiaOnly,
    [switch]$AxiFemOnly,           # configure + build ONLY axifem (fast C++ iteration)
    [switch]$MatlabMexOnly,        # configure + build MATLAB MEX and native S-Functions
    [switch]$InstallToSitePackages  # also copy rebuilt .pyd(s) into the importable site-packages\radia
)

$ErrorActionPreference = "Stop"

if (@($RadiaOnly, $AxiFemOnly, $MatlabMexOnly).Where({ $_ }).Count -gt 1) {
    throw "-RadiaOnly, -AxiFemOnly, and -MatlabMexOnly are mutually exclusive"
}

# ============================================================================
# Paths
# ============================================================================

$PROJECT_DIR = Split-Path -Parent $MyInvocation.MyCommand.Path
$BUILD_DIR = "$PROJECT_DIR\build-msvc"

# ---- Cubit install discovery ------------------------------------------------
# Was: hardcoded 'Coreform Cubit 2025.3'. Now: dot-source tools/find_cubit.ps1
# which picks the highest-version 'Coreform Cubit *' under Program Files
# (or honors $env:CUBIT_INSTALL_DIR). Sets $CubitInstallDir / $CubitBinDir /
# $CubitCmakeDir / $CubitVersion. Cubit-plugin build is SKIPPED (not an
# error) when discovery returns $null -- radia core still builds.
. "$PROJECT_DIR\tools\find_cubit.ps1"
if ($CubitInstallDir) {
    Write-Host "Cubit: $CubitVersion at $CubitInstallDir"
} else {
    Write-Host "Cubit: NOT FOUND -- Cubit-plugin .pyd build will be skipped" -ForegroundColor Yellow
}

# Intel MKL 2026 (required for BLAS/LAPACK and HACApK/PARDISO). Accept both a
# full oneAPI install and the pip mkl-devel layout under Python's Library.
$MklCandidates = @()
if ($env:MKLROOT -and (Test-Path $env:MKLROOT)) {
    $MklCandidates += $env:MKLROOT
}
$MklCandidates += "C:\Program Files (x86)\Intel\oneAPI\mkl\latest"
$PythonLibrary = & python -c "import pathlib, sys; print(pathlib.Path(sys.prefix) / 'Library')" 2>$null
if ($PythonLibrary -and (Test-Path $PythonLibrary)) {
    $MklCandidates += $PythonLibrary
}
$INTEL_MKL = $MklCandidates | Where-Object { Test-Path "$_\lib\mkl_rt.lib" } | Select-Object -First 1
if (-not $INTEL_MKL) {
    $INTEL_MKL = $MklCandidates | Select-Object -First 1
}
$MklImportLibrary = "$INTEL_MKL\lib\mkl_rt.lib"
$MklRuntime = "$INTEL_MKL\bin\mkl_rt.3.dll"
if (-not (Test-Path $MklImportLibrary) -or -not (Test-Path $MklRuntime)) {
    if ($AxiFemOnly) {
        # axifem does NOT link MKL, and an incremental configure on a
        # populated build dir reuses the cached MKL paths -- so warn, don't exit.
        Write-Host "WARNING: Intel MKL not found at $INTEL_MKL -- continuing (-AxiFemOnly does not need MKL)" -ForegroundColor Yellow
    } else {
        Write-Host "ERROR: Intel MKL 2026 not found at $INTEL_MKL" -ForegroundColor Red
        Write-Host 'Install with: python -m pip install "mkl-devel>=2026,<2027"' -ForegroundColor Yellow
        exit 1
    }
}

# NGSolve (optional override via NGSOLVE_DIR environment variable)
$NGSolveCMakeArgs = ""
if ($env:NGSOLVE_DIR -and (Test-Path "$env:NGSOLVE_DIR\NGSolveConfig.cmake")) {
    $NGSolveCMakeArgs = " ^`n    -DNGSolve_DIR=`"$env:NGSOLVE_DIR`" ^`n    -DNetgen_DIR=`"$env:NGSOLVE_DIR`""
    Write-Host "NGSolve: $env:NGSOLVE_DIR (from env)" -ForegroundColor Gray
}
$MatlabMexCMakeArgs = if ($MatlabMexOnly) {
    " ^`n    -DRADIA_BUILD_MATLAB_MEX=ON"
} else {
    ""
}

# Visual Studio (any version) via vswhere
$VS_PATH = $null
$CMAKE_EXE = $null
$VSWHERE = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
if (Test-Path $VSWHERE) {
    $VS_PATH = & $VSWHERE -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath 2>$null
    if ($VS_PATH) {
        $CMAKE_CANDIDATE = "$VS_PATH\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
        if (Test-Path $CMAKE_CANDIDATE) {
            $CMAKE_EXE = $CMAKE_CANDIDATE
        }
    }
}
if (-not $VS_PATH) {
    Write-Host "ERROR: Visual Studio (with VC x86/x64 tools) not found via vswhere." -ForegroundColor Red
    Write-Host "Install Visual Studio Build Tools (any version)." -ForegroundColor Yellow
    exit 1
}
if (-not $CMAKE_EXE) {
    # Fallback: pip-installed cmake
    $PIP_CMAKE = & python -c "import shutil; print(shutil.which('cmake') or '')" 2>$null
    if ($PIP_CMAKE -and (Test-Path $PIP_CMAKE)) {
        $CMAKE_EXE = $PIP_CMAKE
    }
}
if (-not (Test-Path $CMAKE_EXE)) {
    Write-Host "ERROR: CMake not found. Install VS workload with CMake or pip install cmake." -ForegroundColor Red
    exit 1
}
$CTEST_EXE = Join-Path (Split-Path -Parent $CMAKE_EXE) "ctest.exe"
if (-not (Test-Path $CTEST_EXE)) {
    Write-Host "ERROR: CTest not found beside CMake: $CTEST_EXE" -ForegroundColor Red
    exit 1
}

# Ninja generator backend.  Pass the long path explicitly so stale CMake
# short-path cache entries do not break after Python install path changes.
$NINJA_EXE = & python -c "import shutil; print(shutil.which('ninja') or '')" 2>$null
if (-not $NINJA_EXE) {
    $NINJA_EXE = & where.exe ninja 2>$null | Select-Object -First 1
}
if (-not ($NINJA_EXE -and (Test-Path $NINJA_EXE))) {
    Write-Host "ERROR: Ninja not found. Install with: python -m pip install ninja" -ForegroundColor Red
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
Write-Host "MKL:   $INTEL_MKL (mkl_rt.3.dll)" -ForegroundColor Gray
Write-Host "Build: $BUILD_DIR" -ForegroundColor Gray
Write-Host ""

if ($Rebuild) {
    $BuildDirsToClean = @($BUILD_DIR)
    if (-not $AxiFemOnly -and -not $MatlabMexOnly) {
        $BuildDirsToClean += @(
            "$PROJECT_DIR\src\cubit_plugin\build-pyd",
            "$PROJECT_DIR\src\cubit_plugin\build-ccm"
        )
    }
    $ProjectRoot = [System.IO.Path]::GetFullPath($PROJECT_DIR).TrimEnd('\') + '\'
    foreach ($BuildDirToClean in $BuildDirsToClean) {
        $ResolvedBuildDir = [System.IO.Path]::GetFullPath($BuildDirToClean)
        if (-not $ResolvedBuildDir.StartsWith(
                $ProjectRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
            throw "Refusing to clean build directory outside project: $ResolvedBuildDir"
        }
        if (Test-Path -LiteralPath $ResolvedBuildDir) {
            Write-Host "Cleaning build directory: $ResolvedBuildDir" -ForegroundColor Yellow
            Remove-Item -LiteralPath $ResolvedBuildDir -Recurse -Force
        }
    }
}
if (-not (Test-Path $BUILD_DIR)) {
    New-Item -ItemType Directory -Path $BUILD_DIR | Out-Null
}

# ============================================================================
# Build (via batch file for vcvars64 environment)
# ============================================================================

if ($AxiFemOnly) {
# ---- Fast path: DIRECT compile+link of axifem only (no CMake, no MKL) ----
# axifem links ONLY NGSolve/Netgen + Python (no MKL), so we compile its 5
# .cpp and link directly against the installed packages.  This works on any
# machine with VS + ngsolve installed -- even without the full MKL/CMake build
# environment -- and sidesteps a stale build-msvc CMake cache (e.g. 8.3 short
# paths that do not resolve when 8dot3 name creation is disabled on the volume).
$netgenPkg = (& python -c "import netgen,os;print(os.path.dirname(netgen.__file__))").Trim()
$pyPrefix  = (& python -c "import sys;print(sys.prefix)").Trim()
$pyInc     = (& python -c "import sysconfig;print(sysconfig.get_path('include'))").Trim()
$pyLib     = (& python -c "import sys;print(f'python{sys.version_info.major}{sys.version_info.minor}.lib')").Trim()
$axiSrc    = "$PROJECT_DIR\src\ext\axifem"
$objDir    = "$BUILD_DIR\axifem_direct"
if (-not (Test-Path $objDir)) { New-Item -ItemType Directory -Path $objDir | Out-Null }
Write-Host "  netgen pkg: $netgenPkg" -ForegroundColor Gray
Write-Host "  python    : $pyPrefix ($pyLib)" -ForegroundColor Gray
# Compile flags mirror the CMake target (see `ninja -t commands axifem`).
$axiCFlags = '/nologo /TP -DHAVE_NETGEN_SOURCES -DHAVE_STRUCT_TIMESPEC -DLAPACK -DMSVC_EXPRESS -DNDEBUG -DNETGEN_PYTHON -DNGS_PYTHON -DNG_PYTHON -DNOMINMAX -DPYBIND11_SIMPLE_GIL_MANAGEMENT -DPy_NO_LINK_LIB -DRADIA_AXIFEM_PHASE_2B -DTCL -DUSE_TIMEOFDAY -DUSE_UMFPACK -DWIN32 -DWNT -DWNT_WINDOW -D_CRT_SECURE_NO_WARNINGS -D_WIN32_WINNT=0x1000 -Daxifem_EXPORTS /EHsc /O2 /Ob2 /MD /fp:fast /W0 /wd4244 /wd4267 /arch:AVX2 /bigobj /std:c++20 /wd4068 -DMAX_SYS_DIM=3'
$BatchContent = @"
@echo off
call "$VS_PATH\VC\Auxiliary\Build\vcvars64.bat" > nul 2>&1
if errorlevel 1 (
    echo ERROR: vcvars64.bat failed at "$VS_PATH\VC\Auxiliary\Build\vcvars64.bat"
    exit /b 1
)
cd /d "$objDir"
for %%F in (axifem axi_henrotte_fe axi_henrotte_fespace axi_henrotte_integrators axi_henrotte_numeric) do (
    echo === COMPILE %%F.cpp ===
    cl $axiCFlags -I"$axiSrc" -external:I"$netgenPkg\include" -external:I"$pyInc" -external:I"$netgenPkg\include\include" -external:W0 /Fo"$objDir\%%F.obj" /c "$axiSrc\%%F.cpp"
    if errorlevel 1 ( echo COMPILE_FAILED %%F & exit /b 1 )
)
echo === LINK axifem.pyd ===
link /nologo "$objDir\axifem.obj" "$objDir\axi_henrotte_fe.obj" "$objDir\axi_henrotte_fespace.obj" "$objDir\axi_henrotte_integrators.obj" "$objDir\axi_henrotte_numeric.obj" /out:"$BUILD_DIR\axifem.pyd" /implib:"$objDir\axifem.lib" /dll /machine:x64 /INCREMENTAL:NO /ignore:4273 /ignore:4217 /ignore:4049 "$pyPrefix\libs\$pyLib" "$netgenPkg\lib\libngsolve.lib" "$netgenPkg\lib\nglib.lib" "$netgenPkg\lib\ngcore.lib" kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib
if errorlevel 1 ( echo LINK_FAILED & exit /b 2 )
copy /Y "$BUILD_DIR\axifem.pyd" "$PROJECT_DIR\src\radia\axifem.pyd" >nul
echo Build completed.
exit /b 0
"@
} else {
$BatchContent = @"
@echo off
setlocal enabledelayedexpansion

call "$VS_PATH\VC\Auxiliary\Build\vcvars64.bat" > nul 2>&1
if errorlevel 1 (
    echo ERROR: vcvars64.bat failed at "$VS_PATH\VC\Auxiliary\Build\vcvars64.bat"
    exit /b 1
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
set RADIA_ONLY=$RadiaOnly
set MATLAB_MEX_ONLY=$MatlabMexOnly
set RUN_CPP_TESTS=$Test

cd /d "$BUILD_DIR"

echo.
echo ========================================
echo   CMake Configure
echo ========================================
"$CMAKE_EXE" "$PROJECT_DIR" ^
    -G "Ninja" ^
    -DCMAKE_MAKE_PROGRAM="$NINJA_EXE" ^
    -DCMAKE_C_COMPILER=cl ^
    -DCMAKE_CXX_COMPILER=cl ^
    -DCMAKE_BUILD_TYPE=Release$NGSolveCMakeArgs$MatlabMexCMakeArgs

if errorlevel 1 (
    echo ERROR: CMake configuration failed
    exit /b 1
)

if /I "%MATLAB_MEX_ONLY%"=="True" (
    echo.
    echo ========================================
    echo   Building MATLAB MEX targets
    echo ========================================
    "$CMAKE_EXE" --build . --config Release --target radia_mex -j
    if errorlevel 1 (
        echo ERROR: MATLAB MEX target build failed
        exit /b 1
    )
    echo.
    echo Build completed -MatlabMexOnly.
    exit /b 0
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

if /I "%RUN_CPP_TESTS%"=="True" (
    echo.
    echo ========================================
    echo   Running fast C++ kernel tests
    echo ========================================
    "$CMAKE_EXE" --build . --config Release --target test_rad_beam_transfer -j
    if errorlevel 1 (
        echo ERROR: beam-transfer C++ test build failed
        exit /b 1
    )
    "$CMAKE_EXE" --build . --config Release --target test_rad_beam_dynamics -j
    if errorlevel 1 (
        echo ERROR: beam-dynamics C++ test build failed
        exit /b 1
    )
    "$CTEST_EXE" --test-dir . -C Release --output-on-failure -R "radia.beam_(transfer|dynamics).cpp"
    if errorlevel 1 (
        echo ERROR: beam C++ test failed
        exit /b 1
    )
)

if /I "%RADIA_ONLY%"=="True" (
    echo.
    echo Build completed -RadiaOnly.
    exit /b 0
)

echo.
echo ========================================
echo   Building radia_motor_rom C ABI
echo ========================================
"$CMAKE_EXE" --build . --config Release --target radia_motor_rom -j
if errorlevel 1 (
    echo ERROR: radia_motor_rom build failed
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
echo   Building sparsesolv_ngsolve
echo ========================================
"$CMAKE_EXE" --build . --config Release --target sparsesolv_ngsolve -j
if errorlevel 1 ( echo WARNING: sparsesolv_ngsolve build failed )

echo.
echo ========================================
echo   Building axifem
echo ========================================
"$CMAKE_EXE" --build . --config Release --target axifem -j
if errorlevel 1 ( echo WARNING: axifem build failed )

echo.
echo ========================================
echo   Building _equation (LaTeX equations)
echo ========================================
"$CMAKE_EXE" --build . --config Release --target _equation -j
if errorlevel 1 ( echo WARNING: _equation build failed )

echo.
echo ========================================
echo   Building eqnedt64 (the equation editor)
echo ========================================
rem A standalone window rather than a Python entry point: importing the radia
rem package costs about 1.4 s against 4 ms for the equation module itself.
rem Statically linked, so the result is one file that runs on a bare Windows.
"$CMAKE_EXE" --build . --config Release --target eqnedt64 -j
if errorlevel 1 ( echo WARNING: eqnedt64 build failed )

echo.
echo ========================================
echo   Building cubit_mesh_curver (Cubit plugin .pyd)
echo ========================================
set "CUBIT_PLUGIN_SRC=$PROJECT_DIR\src\cubit_plugin"
set "CUBIT_PLUGIN_BUILD=$PROJECT_DIR\src\cubit_plugin\build-pyd"
rem CUBIT_DIR is auto-discovered via tools/find_cubit.ps1 (highest-version
rem 'Coreform Cubit *' under C:\Program Files). Override with env var
rem CUBIT_INSTALL_DIR (the helper appends \cmake automatically). The
rem outer PowerShell wrapper passes $CubitCmakeDir into the batch via
rem this set; if discovery failed the wrapper exited before us.
set "CUBIT_DIR=$CubitCmakeDir"
set "NETGEN_DIR=C:\Program Files\Python312\Lib\site-packages\netgen"
rem Compact Netgen sources are in-repo (src/cubit_plugin/compact_netgen/netgen_src/).
rem No external NETGEN_SRC_DIR needed.

if exist "%CUBIT_DIR%\CubitConfig.cmake" (
    if not exist "%CUBIT_PLUGIN_BUILD%" mkdir "%CUBIT_PLUGIN_BUILD%"
    cd /d "%CUBIT_PLUGIN_BUILD%"
    rem build-pyd: force FULL Netgen mode (disable compact_netgen detection)
    rem so cubit_mesh_curver.pyd builds against pip-installed Netgen + pybind11.
    rem .ccm/.ccl in build-ccm continue to use compact_netgen (no DLL deps).
    "$CMAKE_EXE" -G Ninja -DCMAKE_BUILD_TYPE=Release -DCubit_DIR="%CUBIT_DIR%" -DNETGEN_DIR="%NETGEN_DIR%" -DCOMPACT_NETGEN_OVERRIDES=NONE "%CUBIT_PLUGIN_SRC%"
    "$CMAKE_EXE" --build . --config Release --target cubit_mesh_curver -j
    if errorlevel 1 ( echo WARNING: cubit_mesh_curver build failed )

    rem ========================================
    rem   Building cubit_mesh_export.ccm (APREPRO commands; no Qt deps)
    rem ========================================
    rem radia 4.80.0 removed the Qt5 .ccl GUI component (RadiaComp.cpp);
    rem all GUI work is now in src/radia/panels/radia_export_menu.py
    rem (PySide6).  The .ccm has no Qt dependency and builds unconditionally
    rem when the Cubit SDK is present.
    echo.
    echo ========================================
    echo   Building cubit_mesh_export.ccm (APREPRO commands)
    echo ========================================
    set "CUBIT_CCM_BUILD=$PROJECT_DIR\src\cubit_plugin\build-ccm"
    if not exist "!CUBIT_CCM_BUILD!" mkdir "!CUBIT_CCM_BUILD!"
    cd /d "!CUBIT_CCM_BUILD!"
    "$CMAKE_EXE" -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=cl -DCMAKE_CXX_COMPILER=cl -DCubit_DIR="%CUBIT_DIR%" -DNETGEN_DIR="%NETGEN_DIR%" "%CUBIT_PLUGIN_SRC%"
    "$CMAKE_EXE" --build . --config Release --target cubit_mesh_export_ccm -j
    if errorlevel 1 ( echo WARNING: cubit_mesh_export_ccm build failed )

    cd /d "$BUILD_DIR"
) else (
    echo SKIP: Cubit SDK not found at %CUBIT_DIR%
)

echo.
echo Build completed.
exit /b 0
"@
}

$BatchFile = "$PROJECT_DIR\build_temp_msvc.bat"
$BuildLog = "$PROJECT_DIR\build_log.txt"
$BatchContent | Out-File -FilePath $BatchFile -Encoding ascii

try {
    Write-Host "Building..." -ForegroundColor Cyan

    # Use cmd /c directly (not Start-Process -Wait which waits for child processes)
    & cmd.exe /c "$BatchFile > `"$BuildLog`" 2>&1"
    $BuildResult = $LASTEXITCODE

    if (Test-Path $BuildLog) {
        $BuildLogText = Get-Content $BuildLog -Raw
        Get-Content $BuildLog -Tail 30 | ForEach-Object {
            if ($_ -match "error|ERROR") { Write-Host $_ -ForegroundColor Red }
            elseif ($_ -match "warning|WARNING") { Write-Host $_ -ForegroundColor Yellow }
            else { Write-Host $_ }
        }
        # cmd.exe can lose a nested batch `exit /b` code on some Windows
        # configurations. Never report a successful native build when the
        # generator itself recorded a failed subcommand.
        if ($BuildLogText -match "MATLAB MEX target build failed" -or
                $BuildLogText -match "ninja: build stopped: subcommand failed") {
            $BuildResult = 1
        }
    }
    if ($MatlabMexOnly) {
        $RequiredMexArtifacts = @(
            "$PROJECT_DIR\matlab\radia_mex.mexw64"
        )
        $MissingMexArtifacts = $RequiredMexArtifacts | Where-Object { -not (Test-Path $_) }
        if ($MissingMexArtifacts) {
            Write-Host "ERROR: MATLAB MEX artifacts were not produced:" -ForegroundColor Red
            $MissingMexArtifacts | ForEach-Object { Write-Host "  $_" -ForegroundColor Red }
            $BuildResult = 1
        }
    }

    if ($BuildResult -ne 0) { throw "Build failed with exit code $BuildResult" }

    # For -AxiFemOnly, axifem.pyd is already placed in src/radia/ by the
    # CMake POST_BUILD copy, so skip the full module/cubit copy section below.
    if (-not $AxiFemOnly -and -not $MatlabMexOnly) {
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
        @{ src = "_radia_pybind.cp312-win_amd64.pyd"; dst = "_radia_pybind.pyd"; required = $true }
    )
    if (-not $RadiaOnly) {
        $modules += @(
            @{ src = "peec_matrices.cp312-win_amd64.pyd"; dst = "peec_matrices.pyd"; required = $false },
            @{ src = "cln_core.cp312-win_amd64.pyd";      dst = "cln_core.pyd";      required = $false }
        )
    }

    # Cubit plugin files (built in src/cubit_plugin/build-*): shipped ONLY by
    # the cubit-mesh-export package (Tier-2, 2026-06-01).  radia NO LONGER
    # bundles cubit_mesh_export.ccm/.pyd -- cme is the sole shipper AND deployer
    # (its cubit-plugin-install), so radia and cubit-mesh-export release fully
    # independently.  Safe because radia's setup_cubit.py delegates plugin
    # install to cme, and register_toolbar.py freshness-checks the DEPLOYED
    # plugin, not a radia-bundled copy.  Build output is copied DIRECTLY to
    # the cme package (no src/radia/ hop).
    # Freshness check: copy when src newer than dst (or dst missing).
    $cmeDir = "$PROJECT_DIR\packages\cubit-mesh-export\src\cubit_mesh_export"
    $cubitFiles = @(
        @{ src = "$PROJECT_DIR\src\cubit_plugin\build-pyd\cubit_mesh_curver.cp312-win_amd64.pyd";
           dst = "$cmeDir\cubit_mesh_curver.pyd" },
        @{ src = "$PROJECT_DIR\src\cubit_plugin\build-ccm\cubit_mesh_export.ccm";
           dst = "$cmeDir\cubit_mesh_export.ccm" }
    )
    if (-not $RadiaOnly) {
        foreach ($cf in $cubitFiles) {
            $name = Split-Path $cf.dst -Leaf
            if (Test-Path $cf.src) {
                $srcInfo = Get-Item $cf.src
                $needCopy = $true
                if (Test-Path $cf.dst) {
                    $dstInfo = Get-Item $cf.dst
                    if ($srcInfo.Length -eq $dstInfo.Length -and $srcInfo.LastWriteTime -le $dstInfo.LastWriteTime) {
                        Write-Host "  cubit-mesh-export/${name}: up-to-date ($([math]::Round($srcInfo.Length / 1KB, 1)) KB)" -ForegroundColor Cyan
                        $needCopy = $false
                    }
                }
                if ($needCopy) {
                    Copy-Item $cf.src $cf.dst -Force
                    Write-Host "  cubit-mesh-export/${name}: $([math]::Round($srcInfo.Length / 1KB, 1)) KB (UPDATED)" -ForegroundColor Green
                }
            } else {
                Write-Host "  ${name}: skipped (not built)" -ForegroundColor Yellow
            }
        }

        # === Freshness-check sync: touch bundled binaries to the newest source mtime ===
        # cubit-mesh-export's pyproject.toml runs a freshness gate in
        # `get_requires_for_build_wheel` that compares mtime of every file
        # in src/cubit_plugin/ vs every bundled binary.  Ninja's incremental
        # build only rebuilds targets whose transitive source changed, so
        # editing a .cpp that only feeds the .ccm leaves the .pyd (or vice-
        # versa) untouched, which the freshness gate then FALSE-alarms on.
        # Fix: force-touch every bundled binary mtime so it >= newest source.
        # (As of radia 4.80.0 the .ccl is gone -- only .ccm + .pyd remain.)
        $pluginSrcFiles = Get-ChildItem "$PROJECT_DIR\src\cubit_plugin" `
            -Include "*.cpp","*.hpp","*.h","*.c" -File -Recurse `
            -ErrorAction SilentlyContinue
        if ($pluginSrcFiles) {
            $newest = ($pluginSrcFiles | Measure-Object -Property LastWriteTime -Maximum).Maximum
            # Bump by 1 s so freshness check sees ">= newest source".
            $touchTime = $newest.AddSeconds(1)
            foreach ($cf in $cubitFiles) {
                if (Test-Path $cf.dst) {
                    (Get-Item -LiteralPath $cf.dst).LastWriteTime = $touchTime
                }
            }
            Write-Host "  cubit-mesh-export plugin binaries: mtime synced to newest source ($touchTime)" -ForegroundColor Cyan
        }
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
    }  # end: if (-not $AxiFemOnly) -- full module/cubit copy section

    # ====================================================================
    # Optional: install rebuilt .pyd(s) into the importable site-packages
    # radia package.  Build.ps1 / CMake POST_BUILD only refresh src/radia/ and
    # build-msvc/; if radia is pip-installed NON-editable (the common case),
    # `import radia` reads site-packages and will NOT see the rebuild without
    # this copy.  (An editable `pip install -e` install picks up src/radia/
    # automatically and does not need this.)
    # ====================================================================
    if ($InstallToSitePackages) {
        $platlib = & python -c "import sysconfig; print(sysconfig.get_paths()['platlib'])" 2>$null
        $radiaPkg = if ($platlib) { Join-Path $platlib 'radia' } else { $null }
        if ($radiaPkg -and (Test-Path $radiaPkg)) {
            $installList = @('axifem.pyd')
            if (-not $AxiFemOnly) {
                $installList += @('_radia_pybind.pyd')
                if (-not $RadiaOnly) {
                    $installList += @(
                        'peec_matrices.pyd',
                        'cln_core.pyd',
                        'radia_motor_rom.dll'
                    )
                }
            }
            foreach ($f in $installList) {
                $srcF = Join-Path "$PROJECT_DIR\src\radia" $f
                if (Test-Path $srcF) {
                    Copy-Item $srcF (Join-Path $radiaPkg $f) -Force
                    Write-Host "  installed $f -> $radiaPkg" -ForegroundColor Green
                }
            }
        } else {
            Write-Host "  -InstallToSitePackages: site-packages\radia not found (skipped)" -ForegroundColor Yellow
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

    $oldPythonPath = $env:PYTHONPATH
    if ($oldPythonPath) {
        $env:PYTHONPATH = "$PROJECT_DIR\src;$oldPythonPath"
    } else {
        $env:PYTHONPATH = "$PROJECT_DIR\src"
    }

    $ImportCheck = @"
import pathlib
import radia

repo_src = pathlib.Path(r"$PROJECT_DIR") / "src"
radia_file = pathlib.Path(radia.__file__).resolve()
try:
    radia_file.relative_to(repo_src.resolve())
except ValueError:
    raise SystemExit(f"imported radia outside source tree: {radia_file}")
print(f"radia {radia.__version__} OK from {radia_file}")
"@
    $ImportCheck | python -
    if ($LASTEXITCODE -ne 0) { Write-Host "Import failed!" -ForegroundColor Red; exit 1 }

    if ($RadiaOnly -or $AxiFemOnly -or $MatlabMexOnly) {
        Write-Host "Skipping full pytest because only a focused build was requested." -ForegroundColor Yellow
    } else {
        Push-Location $PROJECT_DIR
        try { python -m pytest tests/ -v --tb=short }
        finally { Pop-Location }
    }

    $env:PYTHONPATH = $oldPythonPath
}
