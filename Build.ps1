#==============================================================================
# Build.ps1 - Radia CPU pybind11 / GPU mexcuda build entry points
#
# Builds _radia_pybind.pyd and auxiliary modules using MSVC compiler
# with Intel MKL for BLAS/LAPACK. Python/NGSolve CPU execution owns TaskManager.
# CUDA MEX compilation and MATLAB/GPU execution do not use TaskManager.
#
# Usage:
#   pwsh -NoProfile -ExecutionPolicy Bypass -File Build.ps1
#   pwsh -NoProfile -ExecutionPolicy Bypass -File Build.ps1 -Rebuild
#   pwsh -NoProfile -ExecutionPolicy Bypass -File Build.ps1 -RadiaOnly
#   pwsh -NoProfile -ExecutionPolicy Bypass -File Build.ps1 -RadiaOnly -RequireNativeProvenance
#   pwsh -NoProfile -ExecutionPolicy Bypass -File Build.ps1 -CudaMexOnly
#   pwsh -NoProfile -ExecutionPolicy Bypass -File Build.ps1 -AmsCudaOnly
#   pwsh -NoProfile -ExecutionPolicy Bypass -File Build.ps1 -OptunaMexOnly
#   pwsh -NoProfile -ExecutionPolicy Bypass -File Build.ps1 -Test
#
# Options:
#   -Rebuild    Clean build directory before building
#   -RadiaOnly  Build and copy only _radia_pybind.pyd
#   -CudaMexOnly  Stage the CUDA-only NGSolve assembled-operator bridge
#   -AmsCudaOnly  Stage the experimental AMS CUDA gateway and snapshot adapter
#   -OptunaMexOnly  Build the Optuna CUDA kernels and CPU host control using mexcuda
#   -RequireNativeProvenance  Reject dirty source before building; otherwise a
#                 dirty development build succeeds without a provenance sidecar
#   -Test       Run source-tree import test + pytest after build
#   -Verbose    Show detailed build output
#
# Requirements:
#   - Visual Studio Build Tools (any version with VC.Tools.x86.x64) — auto-detected via vswhere
#   - MATLAB with Parallel Computing Toolbox and a supported mexcuda compiler
#   - Intel MKL from `python -m pip install mkl-devel`
#   - NGSolve (pip install or source build)
#   - Python 3.12 with pybind11
#==============================================================================

param(
    [ValidateRange(1, 64)]
    [int]$BuildJobs = 4,           # bound compiler memory on shared CI hosts
    [switch]$Rebuild,
    [switch]$Test,
    [switch]$Verbose,
    [switch]$RadiaOnly,
    [switch]$HMatrixCudaOnly,    # isolated resident HACApK leaf application
    [string]$HMatrixMklRoot,     # explicit isolated MKL for complete C HACApK assembly
    [switch]$HMatrixFaultInjection, # isolated diagnostic hooks, never a release default
    [switch]$HDivPythonBindings,  # MEX-free shared resident operator component
    [string]$GpuBuildAdmission,  # private reviewed HDiv GPU worker/task/source admission
    [string]$GpuBuildAdmissionSha256,
    [switch]$ProjectionCudaOnly, # experimental NGSolve target replay bridge
    [switch]$ProjectionLifecycleDiagnostic, # isolated weak-owner lifetime probe, not a release default
    [switch]$IhCudaOnly,          # isolated persistent coupled IH CUDA gateway
    [switch]$FieldCudaOnly,       # isolated retained-cloud RadFld CUDA gateway
    [switch]$FmmCudaOnly,         # isolated numeric Laplace FMM CUDA gateway
    [switch]$AmsCudaOnly,         # experimental AMS snapshot adapter + CUDA MEX
    [switch]$AxiFemOnly,           # configure + build ONLY axifem (fast C++ iteration)
    [switch]$CsrPythonBindings,     # optional direct Python CSR CUDA binding
    [switch]$ResidentAmgPython,    # explicit staged resident algebraic GPU validation module
    [switch]$CudaMexOnly,          # isolated CUDA-only assembled-operator gateway
    [switch]$OptunaMexOnly,        # mexcuda kernels + CPU host control
    [switch]$RequireNativeProvenance, # reject dirty source before a provenance-bearing build
    [switch]$InstallToSitePackages  # also copy rebuilt .pyd(s) into the importable site-packages\radia
)

$ErrorActionPreference = "Stop"
if ($CsrPythonBindings -and -not $CudaMexOnly) { throw '-CsrPythonBindings requires -CudaMexOnly' }
if ($HDivPythonBindings -and -not $HMatrixCudaOnly) { throw '-HDivPythonBindings requires -HMatrixCudaOnly' }
if ($ResidentAmgPython -and -not $CudaMexOnly) { throw "Resident AMG Python requires -CudaMexOnly" }
if ($ProjectionLifecycleDiagnostic -and -not $ProjectionCudaOnly) { throw 'Projection lifecycle diagnostic requires the isolated Projection target' }
if ([bool]$GpuBuildAdmission -ne [bool]$GpuBuildAdmissionSha256) { throw 'GPU build admission and hash must be supplied together' }
if ($GpuBuildAdmission -and -not ($HMatrixCudaOnly -or $ProjectionCudaOnly -or $FmmCudaOnly -or $FieldCudaOnly -or $IhCudaOnly -or $AmsCudaOnly -or $CudaMexOnly)) { throw 'GPU build admission requires an isolated native CUDA target' }
if ($FmmCudaOnly) {
    if ($RadiaOnly -or $HMatrixCudaOnly -or $ProjectionCudaOnly -or $IhCudaOnly -or $FieldCudaOnly -or $AmsCudaOnly -or $AxiFemOnly -or $CudaMexOnly -or $OptunaMexOnly -or $InstallToSitePackages -or $Rebuild -or $Test -or $HMatrixMklRoot -or $args.Count) {
        throw '-FmmCudaOnly is an isolated experimental build'
    }
    $fmmBuildArgs = @()
    if ($GpuBuildAdmission) { $fmmBuildArgs += @('--build-admission', $GpuBuildAdmission, '--build-admission-sha256', $GpuBuildAdmissionSha256) }
    if ($RequireNativeProvenance) { $fmmBuildArgs += '--require-clean-source' }
    & python (Join-Path $PSScriptRoot 'packages/radia-simulink/tools/build_fmm_cuda.py') @fmmBuildArgs
    if ($LASTEXITCODE -ne 0) { throw 'FMM CUDA build failed' }
    return
}
if ($HMatrixMklRoot -and -not $HMatrixCudaOnly) { throw '-HMatrixMklRoot requires -HMatrixCudaOnly' }
if ($HMatrixFaultInjection -and -not ($HMatrixCudaOnly -and $HMatrixMklRoot)) { throw '-HMatrixFaultInjection requires the full isolated HMatrix CUDA target' }
if ($args.Count -gt 0) {
    throw "Unsupported build arguments: $($args -join ' '). Radia CPU uses -RadiaOnly; GPU uses -CudaMexOnly."
}

if ($HMatrixCudaOnly) {
    if ($RadiaOnly -or $ProjectionCudaOnly -or $IhCudaOnly -or $FieldCudaOnly -or $AxiFemOnly -or $AmsCudaOnly -or $CudaMexOnly -or $OptunaMexOnly -or $InstallToSitePackages -or $Rebuild -or $Test) {
        throw '-HMatrixCudaOnly is an isolated experimental build'
    }
    $hMatrixBuildArgs = @()
    if ($HMatrixMklRoot) { $hMatrixBuildArgs += @('--mkl-root', $HMatrixMklRoot) }
    if ($HMatrixFaultInjection) { $hMatrixBuildArgs += '--fault-injection' }
    if ($HDivPythonBindings) { $hMatrixBuildArgs += '--python-bindings' }
    if ($RequireNativeProvenance) { $hMatrixBuildArgs += '--require-clean-source' }
    if ($GpuBuildAdmission) { $hMatrixBuildArgs += @('--build-admission', $GpuBuildAdmission, '--build-admission-sha256', $GpuBuildAdmissionSha256) }
    & python (Join-Path $PSScriptRoot 'packages/radia-simulink/tools/build_hmatrix_cuda.py') @hMatrixBuildArgs
    if ($LASTEXITCODE -ne 0) { throw 'H-matrix CUDA experimental build failed' }
    return
}

if ($ProjectionCudaOnly) {
    if ($RadiaOnly -or $IhCudaOnly -or $FieldCudaOnly -or $AxiFemOnly -or $AmsCudaOnly -or $CudaMexOnly -or $OptunaMexOnly -or $InstallToSitePackages -or $Rebuild -or $Test) {
        throw '-ProjectionCudaOnly is an isolated experimental build'
    }
    $projectionBuildArgs = @()
    if ($ProjectionLifecycleDiagnostic) { $projectionBuildArgs += '--lifecycle-diagnostic' }
    if ($RequireNativeProvenance) { $projectionBuildArgs += '--require-clean-source' }
    if ($GpuBuildAdmission) { $projectionBuildArgs += @('--build-admission', $GpuBuildAdmission, '--build-admission-sha256', $GpuBuildAdmissionSha256) }
    & python (Join-Path $PSScriptRoot 'packages/radia-simulink/tools/build_projection_cuda.py') @projectionBuildArgs
    if ($LASTEXITCODE -ne 0) { throw 'Projection CUDA experimental build failed' }
    return
}

if ($IhCudaOnly) {
    if ($RadiaOnly -or $FieldCudaOnly -or $AxiFemOnly -or $AmsCudaOnly -or $CudaMexOnly -or $OptunaMexOnly -or $InstallToSitePackages -or $Rebuild -or $Test) {
        throw '-IhCudaOnly is an isolated experimental build; incompatible options supplied'
    }
    $ihBuildArgs = @()
    if ($GpuBuildAdmission) { $ihBuildArgs += @('--build-admission', $GpuBuildAdmission, '--build-admission-sha256', $GpuBuildAdmissionSha256) }
    if ($RequireNativeProvenance) { $ihBuildArgs += '--require-clean-source' }
    & python (Join-Path $PSScriptRoot 'packages/radia-simulink/tools/build_ih_cuda.py') @ihBuildArgs
    if ($LASTEXITCODE -ne 0) { throw 'IH CUDA experimental build failed' }
    return
}

if ($FieldCudaOnly) {
    if ($RadiaOnly -or $AxiFemOnly -or $AmsCudaOnly -or $CudaMexOnly -or $OptunaMexOnly -or $InstallToSitePackages -or $Rebuild -or $Test) {
        throw '-FieldCudaOnly is an isolated experimental build; incompatible options supplied'
    }
    $fieldBuildArgs = @()
    if ($GpuBuildAdmission) { $fieldBuildArgs += @('--build-admission', $GpuBuildAdmission, '--build-admission-sha256', $GpuBuildAdmissionSha256) }
    if ($RequireNativeProvenance) { $fieldBuildArgs += '--require-clean-source' }
    & python (Join-Path $PSScriptRoot 'packages/radia-simulink/tools/build_field_cuda.py') @fieldBuildArgs
    if ($LASTEXITCODE -ne 0) { throw 'RadFld CUDA experimental build failed' }
    return
}

if ($AmsCudaOnly) {
    if ($RadiaOnly -or $AxiFemOnly -or $CudaMexOnly -or $OptunaMexOnly -or $InstallToSitePackages -or $Rebuild) {
        throw '-AmsCudaOnly is an isolated experimental build; incompatible options supplied'
    }
    $amsBuildArgs = @()
    if ($GpuBuildAdmission) { $amsBuildArgs += @('--build-admission', $GpuBuildAdmission, '--build-admission-sha256', $GpuBuildAdmissionSha256) }
    if ($RequireNativeProvenance) { $amsBuildArgs += '--require-clean-source' }
    & python (Join-Path $PSScriptRoot 'packages/radia-simulink/tools/build_ams_cuda.py') @amsBuildArgs
    if ($LASTEXITCODE -ne 0) { throw 'AMS CUDA experimental build failed' }
    return
}

function Get-NativeBuildResult {
    param(
        [int]$ExitCode,
        [string]$LogText,
        [string[]]$RequiredArtifacts = @()
    )
    if ($ExitCode -ne 0) { return $ExitCode }
    # Nested cmd.exe batches can lose exit /b codes. Existing files do not
    # make a failed build valid, so inspect the log before checking artifacts.
    if ($LogText -match '(?im)^(?:ERROR: .*\b(?:build|configuration) failed|ninja: (?:error:|build stopped:)|CMake Error)') {
        return 1
    }
    foreach ($artifact in $RequiredArtifacts) {
        if (-not (Test-Path -LiteralPath $artifact -PathType Leaf)) {
            Write-Host "ERROR: required native artifact was not produced: $artifact" -ForegroundColor Red
            return 1
        }
    }
    return 0
}

if (@($RadiaOnly, $AxiFemOnly, $CudaMexOnly, $OptunaMexOnly).Where({ $_ }).Count -gt 1) {
    throw "-RadiaOnly, -AxiFemOnly, -CudaMexOnly, and -OptunaMexOnly are mutually exclusive"
}
if ($CudaMexOnly) {
    if ($InstallToSitePackages) { throw 'CUDA MEX builds remain staged; no runtime replacement' }
    $cudaRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
    if ($RequireNativeProvenance) {
        $cudaStatus = & git -C $cudaRoot status --porcelain
        if ($LASTEXITCODE -ne 0 -or $cudaStatus) { throw 'CUDA native provenance requires a clean source checkout' }
    }
    $cudaMexBuildArgs = @()
    if ($CsrPythonBindings) { $cudaMexBuildArgs += '--python-bindings' }
    if ($ResidentAmgPython) { $cudaMexBuildArgs += "--resident-amg-python" }
    if ($GpuBuildAdmission) { $cudaMexBuildArgs += @('--build-admission', $GpuBuildAdmission, '--build-admission-sha256', $GpuBuildAdmissionSha256) }
    & python (Join-Path $cudaRoot 'packages/radia-simulink/tools/build_cuda_mex.py') --out (Join-Path $cudaRoot 'build-cuda-mex') @cudaMexBuildArgs
    if ($LASTEXITCODE -ne 0) { throw 'CUDA MEX build failed' }
    return
}
if ($OptunaMexOnly) {
    if ($InstallToSitePackages -or $Rebuild) {
        throw 'Optuna mexcuda build is isolated; install/rebuild switches are unsupported'
    }
    $optunaRoot = $PSScriptRoot.Replace("'", "''")
    if ($RequireNativeProvenance) {
        $optunaStatus = & git -C $PSScriptRoot status --porcelain
        if ($LASTEXITCODE -ne 0 -or $optunaStatus) { throw 'Optuna provenance requires clean source' }
    }
    & matlab -batch "addpath('$optunaRoot/packages/radia-optuna'); buildOptunaMex;"
    if ($LASTEXITCODE -ne 0) { throw 'Optuna mexcuda build failed' }
    if ($Test) {
        & matlab -batch "addpath('$optunaRoot/validation_test/optimization'); validate_optuna_cuda('$optunaRoot'); addpath('$optunaRoot/validation_test/radia_optuna'); validate_cuda('$optunaRoot','$optunaRoot/packages/radia-optuna/matlab','$optunaRoot/build-optuna-cuda/validation.json');"
        if ($LASTEXITCODE -ne 0) { throw 'Optuna CUDA validation failed' }
    }
    return
}

# ============================================================================
# Paths
# ============================================================================

$PROJECT_DIR = Split-Path -Parent $MyInvocation.MyCommand.Path
$BUILD_DIR = "$PROJECT_DIR\build-msvc"
$NativeBuildProvenanceScript = "$PROJECT_DIR\tools\native_build_provenance.ps1"
. $NativeBuildProvenanceScript
$NativeBuildSourceIdentity = $null
$NativeProvenanceBinaries = @()
if (-not $AxiFemOnly) {
    $NativeProvenanceBinaries += "$PROJECT_DIR\src\radia\_radia_pybind.pyd"
}
if ($NativeProvenanceBinaries.Count -gt 0) {
    $NativeBuildSourceIdentity = Get-NativeBuildSourceIdentity -RepoRoot $PROJECT_DIR
    foreach ($NativeBinary in $NativeProvenanceBinaries) {
        Clear-NativeBuildProvenance -BinaryPath $NativeBinary
    }
    if ($RequireNativeProvenance -and $NativeBuildSourceIdentity.source_dirty) {
        throw "Native provenance was required, but the source checkout is dirty"
    }
}

$PythonCommand = Get-Command python -ErrorAction SilentlyContinue
$PythonExecutable = if ($PythonCommand) { $PythonCommand.Source } else { "" }
if (-not $PythonExecutable) {
    throw "Python is required for Radia builds and for -Test"
}
$Pybind11CMakeDir = ""

$Pybind11CMakeDir = (& $PythonExecutable -c "import pybind11; print(pybind11.get_cmake_dir())").Trim()
if (-not (Test-Path "$Pybind11CMakeDir\pybind11Config.cmake")) {
    Write-Host "ERROR: pybind11 CMake package not found under $Pybind11CMakeDir" -ForegroundColor Red
    exit 1
}

# ---- Cubit install discovery ------------------------------------------------

# Picks the highest-version 'Coreform Cubit *' under Program Files or
# honors $env:CUBIT_INSTALL_DIR. Cubit absence does not block Radia core.
. "$PROJECT_DIR\tools\find_cubit.ps1"
if ($CubitInstallDir) {
    Write-Host "Cubit: $CubitVersion at $CubitInstallDir"
} else {
    Write-Host "Cubit: NOT FOUND -- Cubit-plugin .pyd build will be skipped" -ForegroundColor Yellow
}

$DetectedCubitCmakeDir = if ($CubitInstallDir) {
    Join-Path $CubitInstallDir "cmake"
} else {
    ""
}
$BuildCubitPlugin = [bool]($CubitInstallDir -and
    (Test-Path (Join-Path $CubitInstallDir "cmake\CubitConfig.cmake")))
Write-Host "Cubit native targets enabled: $BuildCubitPlugin" -ForegroundColor Gray

# Intel MKL 2026 (required for BLAS/LAPACK and HACApK). The selected
# Python environment is authoritative; MKLROOT is an explicit fallback only.

$MklCandidates = @()
$PythonLibrary = & $PythonExecutable -c "import pathlib, sys; print(pathlib.Path(sys.prefix) / 'Library')" 2>$null
if ($PythonLibrary -and (Test-Path $PythonLibrary)) {
    $MklCandidates += $PythonLibrary
}
if ($env:MKLROOT -and (Test-Path $env:MKLROOT)) {
    $MklCandidates += $env:MKLROOT
}
$MklRuntimeNames = @("mkl_rt", "mkl_core", "mkl_intel_thread", "mkl_avx2", "mkl_def")
$INTEL_MKL = ""
foreach ($candidate in ($MklCandidates | Select-Object -Unique)) {
    if (-not (Test-Path "$candidate\lib\mkl_rt.lib")) {
        continue
    }
    $missingRuntime = @($MklRuntimeNames | Where-Object {
        -not (Test-Path "$candidate\bin\$_.3.dll")
    })
    if ($missingRuntime.Count -eq 0) {
        $INTEL_MKL = $candidate
        break
    }
    Write-Host "Ignoring incomplete MKL candidate $candidate (missing $($missingRuntime -join ', '))" -ForegroundColor Yellow
}
$MklImportLibrary = "$INTEL_MKL\lib\mkl_rt.lib"
$MklRuntime = "$INTEL_MKL\bin\mkl_rt.3.dll"
if (-not (Test-Path $MklImportLibrary) -or -not (Test-Path $MklRuntime)) {
    if ($AxiFemOnly) {
    # These focused targets do not link MKL. An incremental configure on a
    # populated build dir may retain cached MKL paths, so warn, don't exit.
        Write-Host "WARNING: Intel MKL not found at $INTEL_MKL -- continuing (-AxiFemOnly does not need MKL)" -ForegroundColor Yellow
    } else {
        Write-Host "ERROR: complete Intel MKL 2026 runtime not found" -ForegroundColor Red
        Write-Host 'Install with: python -m pip install "mkl-devel>=2026,<2027"' -ForegroundColor Yellow
        exit 1
    }
}

# NGSolve (optional override via NGSOLVE_DIR / Netgen_DIR environment variables)
$NGSolveCMakeArgs = ""
$NetgenPackageDir = $null
if ($env:NGSOLVE_DIR -and (Test-Path "$env:NGSOLVE_DIR\NGSolveConfig.cmake")) {
    $NetgenCMakeDir = if ($env:Netgen_DIR -and
        (Test-Path "$env:Netgen_DIR\NetgenConfig.cmake")) {
        $env:Netgen_DIR
    } else {
        $env:NGSOLVE_DIR
    }
    $NGSolveCMakeArgs = " ^`n    -DNGSolve_DIR=`"$env:NGSOLVE_DIR`" ^`n    -DNetgen_DIR=`"$NetgenCMakeDir`""
    Write-Host "NGSolve: $env:NGSOLVE_DIR (from env)" -ForegroundColor Gray
}

$NetgenPackageDir = (& $PythonExecutable -c "import netgen,os;print(os.path.dirname(netgen.__file__))").Trim()
if (-not (Test-Path "$NetgenPackageDir\include")) {
    Write-Host "ERROR: Netgen include directory not found under $NetgenPackageDir" -ForegroundColor Red
    exit 1
}

$PythonCMakeArgs = " ^`n    -DPython3_EXECUTABLE=`"$PythonExecutable`" ^`n    -Dpybind11_DIR=`"$Pybind11CMakeDir`""

# CPU toolchain follows the MATLAB-compatible VS2022 / MSVC 14.44 policy.
$VS_PATH = $null
$CMAKE_EXE = $null
$VSWHERE = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
if (Test-Path $VSWHERE) {
    $VS_PATH = & $VSWHERE -latest -version '[17.0,18.0)' -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath 2>$null
    if ($VS_PATH) {
        $CMAKE_CANDIDATE = "$VS_PATH\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
        if (Test-Path $CMAKE_CANDIDATE) {
            $CMAKE_EXE = $CMAKE_CANDIDATE
        }
    }
}
if (-not $VS_PATH) {
    Write-Host "ERROR: Visual Studio Build Tools 2022 not found via vswhere." -ForegroundColor Red
    exit 1
}
$ToolsetVersionFile = Join-Path $VS_PATH 'VC/Auxiliary/Build/Microsoft.VCToolsVersion.default.txt'
if (-not (Test-Path -LiteralPath $ToolsetVersionFile)) {
    throw "MSVC default toolset identity is missing: $ToolsetVersionFile"
}
$ToolsetVersion = (Get-Content -LiteralPath $ToolsetVersionFile -Raw).Trim()
$CompilerPath = Join-Path $VS_PATH "VC/Tools/MSVC/$ToolsetVersion/bin/Hostx64/x64/cl.exe"
if ($ToolsetVersion -notmatch '^14\.44\.[0-9]+$' -or
    -not (Test-Path -LiteralPath $CompilerPath)) {
    throw "Build requires the default MSVC 14.44 x64 toolset; found '$ToolsetVersion' at '$VS_PATH'"
}
$CompilerVersion = (Get-Item -LiteralPath $CompilerPath).VersionInfo.FileVersion
if ($CompilerVersion -notmatch '^19\.44\.') {
    throw "Build requires cl.exe 19.44; found '$CompilerVersion' at '$CompilerPath'"
}
if (-not $CMAKE_EXE) {
    # Fallback: pip-installed cmake
    $PIP_CMAKE = & $PythonExecutable -c "import shutil; print(shutil.which('cmake') or '')" 2>$null
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
$NINJA_COMMAND = Get-Command ninja -ErrorAction SilentlyContinue
$NINJA_EXE = if ($NINJA_COMMAND) { $NINJA_COMMAND.Source } else { "" }
if (-not $NINJA_EXE -and $PythonExecutable) {
    $NINJA_EXE = & $PythonExecutable -c "import shutil; print(shutil.which('ninja') or '')" 2>$null
}
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
    if ($BuildCubitPlugin -and -not ($RadiaOnly -or $AxiFemOnly)) {
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
$netgenPkg = (& $PythonExecutable -c "import netgen,os;print(os.path.dirname(netgen.__file__))").Trim()
$pyPrefix  = (& $PythonExecutable -c "import sys;print(sys.prefix)").Trim()
$pyBasePrefix = (& $PythonExecutable -c "import sys;print(sys.base_prefix)").Trim()
$pyInc     = (& $PythonExecutable -c "import sysconfig;print(sysconfig.get_path('include'))").Trim()
$pyLib     = (& $PythonExecutable -c "import sys;print(f'python{sys.version_info.major}{sys.version_info.minor}.lib')").Trim()
$PythonImportLibrary = Join-Path $pyBasePrefix "libs\$pyLib"
if (-not (Test-Path -LiteralPath $PythonImportLibrary -PathType Leaf)) {
    throw "Python import library not found at base interpreter path: $PythonImportLibrary"
}
$axiSrc    = "$PROJECT_DIR\src\ext\axifem"
$objDir    = "$BUILD_DIR\axifem_direct"
if (-not (Test-Path $objDir)) { New-Item -ItemType Directory -Path $objDir | Out-Null }
Write-Host "  netgen pkg: $netgenPkg" -ForegroundColor Gray
Write-Host "  python    : $pyPrefix ($PythonImportLibrary)" -ForegroundColor Gray
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
link /nologo "$objDir\axifem.obj" "$objDir\axi_henrotte_fe.obj" "$objDir\axi_henrotte_fespace.obj" "$objDir\axi_henrotte_integrators.obj" "$objDir\axi_henrotte_numeric.obj" /out:"$BUILD_DIR\axifem.pyd" /implib:"$objDir\axifem.lib" /dll /machine:x64 /INCREMENTAL:NO /ignore:4273 /ignore:4217 /ignore:4049 "$PythonImportLibrary" "$netgenPkg\lib\libngsolve.lib" "$netgenPkg\lib\nglib.lib" "$netgenPkg\lib\ngcore.lib" kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib
if errorlevel 1 ( echo LINK_FAILED & exit /b 2 )
copy /Y "$BUILD_DIR\axifem.pyd" "$PROJECT_DIR\src\radia\axifem.pyd" >nul
echo Build completed.
exit /b 0
"@
} else {
$MklBatchEnvironment = @"
set LIB=$INTEL_MKL\lib;%LIB%
set INCLUDE=$INTEL_MKL\include;%INCLUDE%
set MKLROOT=$INTEL_MKL
"@

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

$MklBatchEnvironment
set RADIA_ONLY=$RadiaOnly
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
    -DCMAKE_BUILD_TYPE=Release$PythonCMakeArgs$NGSolveCMakeArgs

if errorlevel 1 (
    echo ERROR: CMake configuration failed
    exit /b 1
)

echo.
echo ========================================
echo   Building _radia_pybind
echo ========================================
"$CMAKE_EXE" --build . --config Release --target _radia_pybind -j $BuildJobs
if errorlevel 1 (
    echo ERROR: _radia_pybind build failed
    exit /b 1
)

if /I "%RUN_CPP_TESTS%"=="True" (
    echo.
    echo ========================================
    echo   Running fast C++ kernel tests
    echo ========================================
    "$CMAKE_EXE" --build . --config Release --target test_rad_beam_transfer -j $BuildJobs
    if errorlevel 1 (
        echo ERROR: beam-transfer C++ test build failed
        exit /b 1
    )
    "$CMAKE_EXE" --build . --config Release --target test_rad_beam_dynamics -j $BuildJobs
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
"$CMAKE_EXE" --build . --config Release --target radia_motor_rom -j $BuildJobs
if errorlevel 1 (
    echo ERROR: radia_motor_rom build failed
    exit /b 1
)

echo.
echo ========================================
echo   Building peec_matrices
echo ========================================
"$CMAKE_EXE" --build . --config Release --target peec_matrices -j $BuildJobs
if errorlevel 1 (
    echo ERROR: peec_matrices build failed
    exit /b 1
)

echo.
echo ========================================
echo   Building sparsesolv_ngsolve
echo ========================================
"$CMAKE_EXE" --build . --config Release --target sparsesolv_ngsolve -j $BuildJobs
if errorlevel 1 (
    echo ERROR: sparsesolv_ngsolve build failed
    exit /b 1
)

echo.
echo ========================================
echo   Building axifem
echo ========================================
"$CMAKE_EXE" --build . --config Release --target axifem -j $BuildJobs
if errorlevel 1 (
    echo ERROR: axifem build failed
    exit /b 1
)

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
set "CUBIT_DIR=$DetectedCubitCmakeDir"
set "NETGEN_DIR=$NetgenPackageDir"
rem Compact Netgen sources are in-repo (src/cubit_plugin/compact_netgen/netgen_src/).
rem No external NETGEN_SRC_DIR needed.

rem The PowerShell discovery result and paths are embedded as fixed batch
rem literals. Do not rediscover Cubit from batch state; vcvars may also change
rem delayed-expansion state, so this block intentionally contains no !VAR! use.
if /I "$BuildCubitPlugin"=="True" (
    if not exist "%CUBIT_PLUGIN_BUILD%" mkdir "%CUBIT_PLUGIN_BUILD%"
    cd /d "%CUBIT_PLUGIN_BUILD%"
    rem build-pyd: force FULL Netgen mode; disable compact_netgen detection
    rem so cubit_mesh_curver.pyd builds against pip-installed Netgen + pybind11.
    rem .ccm/.ccl in build-ccm continue to use compact_netgen without DLL deps.
    "$CMAKE_EXE" -G Ninja -DCMAKE_BUILD_TYPE=Release -DPython3_EXECUTABLE="$PythonExecutable" -Dpybind11_DIR="$Pybind11CMakeDir" -DCubit_DIR="$DetectedCubitCmakeDir" -DNETGEN_DIR="$NetgenPackageDir" -DCOMPACT_NETGEN_OVERRIDES=NONE -DCUBIT_CALLBACK_CONTRACT_TESTS=ON "$PROJECT_DIR\src\cubit_plugin"
    if errorlevel 1 (
        echo ERROR: cubit_mesh_curver configuration failed
        exit /b 1
    )
    "$CMAKE_EXE" --build . --config Release --target cubit_mesh_curver -j $BuildJobs
    if errorlevel 1 (
        echo ERROR: cubit_mesh_curver build failed
        exit /b 1
    )
    "$CMAKE_EXE" --build . --config Release --target cubit_callback_contract -j $BuildJobs
    if errorlevel 1 (
        echo ERROR: cubit_callback_contract build failed
        exit /b 1
    )
    "$CTEST_EXE" --test-dir "%CUBIT_PLUGIN_BUILD%" --output-on-failure
    if errorlevel 1 (
        echo ERROR: native Cubit callback contract failed
        exit /b 1
    )

    rem ========================================
    rem   Building cubit_mesh_export.ccm; APREPRO commands without Qt deps
    rem ========================================
    rem radia 4.80.0 removed the Qt5 RadiaComp.cpp .ccl GUI component;
    rem all GUI work is now in src/radia/panels/radia_export_menu.py using
    rem PySide6.  The .ccm has no Qt dependency and builds unconditionally
    rem when the Cubit SDK is present.
    echo.
    echo ========================================
    echo   Building cubit_mesh_export.ccm - APREPRO commands
    echo ========================================
    if not exist "$PROJECT_DIR\src\cubit_plugin\build-ccm" mkdir "$PROJECT_DIR\src\cubit_plugin\build-ccm"
    cd /d "$PROJECT_DIR\src\cubit_plugin\build-ccm"
    "$CMAKE_EXE" -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=cl -DCMAKE_CXX_COMPILER=cl -DPython3_EXECUTABLE="$PythonExecutable" -Dpybind11_DIR="$Pybind11CMakeDir" -DCubit_DIR="$DetectedCubitCmakeDir" -DNETGEN_DIR="$NetgenPackageDir" "$PROJECT_DIR\src\cubit_plugin"
    if errorlevel 1 (
        echo ERROR: cubit_mesh_export_ccm configuration failed
        exit /b 1
    )
    "$CMAKE_EXE" --build . --config Release --target cubit_mesh_export_ccm -j $BuildJobs
    if errorlevel 1 (
        echo ERROR: cubit_mesh_export_ccm build failed
        exit /b 1
    )

    cd /d "$BUILD_DIR"
) else (
    echo SKIP: Cubit SDK not found at $DetectedCubitCmakeDir
)

echo.
echo Build completed.
exit /b 0
"@
}

$BatchFile = "$BUILD_DIR\build_temp_msvc.bat"
$BuildLog = "$BUILD_DIR\build.log"
$BatchContent | Out-File -FilePath $BatchFile -Encoding ascii

try {
    Write-Host "Building..." -ForegroundColor Cyan

    # Use cmd /c directly (not Start-Process -Wait which waits for child processes)
    & cmd.exe /c "$BatchFile > `"$BuildLog`" 2>&1"
    $BuildResult = $LASTEXITCODE
    $BuildLogText = ""

    if (Test-Path $BuildLog) {
        $BuildLogText = Get-Content $BuildLog -Raw
        Get-Content $BuildLog -Tail 30 | ForEach-Object {
            if ($_ -match "error|ERROR") { Write-Host $_ -ForegroundColor Red }
            elseif ($_ -match "warning|WARNING") { Write-Host $_ -ForegroundColor Yellow }
            else { Write-Host $_ }
        }
    }
    $RequiredNativeArtifacts = @()
    if ($BuildCubitPlugin -and -not ($RadiaOnly -or $AxiFemOnly)) {
        $RequiredNativeArtifacts += "$PROJECT_DIR\src\cubit_plugin\build-pyd\cubit_mesh_curver.cp312-win_amd64.pyd"
        $RequiredNativeArtifacts += "$PROJECT_DIR\src\cubit_plugin\build-ccm\cubit_mesh_export.ccm"
    }
    $BuildResult = Get-NativeBuildResult -ExitCode $BuildResult -LogText $BuildLogText -RequiredArtifacts $RequiredNativeArtifacts

    if ($BuildResult -ne 0) { throw "Build failed with exit code $BuildResult" }

    # For -AxiFemOnly, axifem.pyd is already placed in src/radia/ by the
    # CMake POST_BUILD copy, so skip the module copy section below.
    if (-not $AxiFemOnly) {
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
            @{ src = "peec_matrices.cp312-win_amd64.pyd"; dst = "peec_matrices.pyd"; required = $true }
        )
    }

    # Raw Cubit native outputs remain in src/cubit_plugin/build-*.
    # Signed package payloads and their manifest are promoted only through
    # the reviewed release process; an ordinary build never copies or
    # touches them.

    foreach ($mod in $modules) {
        $srcPath = "$BUILD_DIR\$($mod.src)"
        $dstPath = "$PROJECT_DIR\src\radia\$($mod.dst)"
        if (Test-Path $srcPath) {
            $srcInfo = Get-Item $srcPath
            $needCopy = $true
            if (Test-Path $dstPath) {
                $dstInfo = Get-Item $dstPath
                $srcHash = (Get-FileHash -LiteralPath $srcPath -Algorithm SHA256).Hash
                $dstHash = (Get-FileHash -LiteralPath $dstPath -Algorithm SHA256).Hash
                if ($srcHash -eq $dstHash) {
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

    Write-NativeBuildProvenance `
        -BinaryPath "$PROJECT_DIR\src\radia\_radia_pybind.pyd" `
        -RepoRoot $PROJECT_DIR `
        -BuildDirectory 'build-msvc' `
        -StartIdentity $NativeBuildSourceIdentity

    # A copied .pyd is not necessarily usable: changing the pinned NGSolve
    # release can leave a loadable-looking but ABI-incompatible binary behind.
    # Probe each required extension in its own process and verify that Python
    # resolved it from this source tree rather than from stale site-packages.
    $NativeAbiModules = @("radia._radia_pybind")
    if (-not $RadiaOnly) {
        $NativeAbiModules += "radia.peec_matrices"
    }
    $NativeAbiArgs = @(
        "$PROJECT_DIR\tools\check_native_module_abi.py",
        "--source-root", "$PROJECT_DIR\src"
    )
    foreach ($NativeAbiModule in $NativeAbiModules) {
        $NativeAbiArgs += @("--module", $NativeAbiModule)
    }
    & $PythonExecutable @NativeAbiArgs
    if ($LASTEXITCODE -ne 0) {
        throw "Native module ABI check failed"
    }
    }  # end: if (-not $AxiFemOnly) -- module copy section

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

    if ($RadiaOnly -or $AxiFemOnly) {
        Write-Host "Skipping full pytest because only a focused build was requested." -ForegroundColor Yellow
    } else {
        Push-Location $PROJECT_DIR
        try { python -m pytest tests/ -v --tb=short }
        finally { Pop-Location }
    }

    $env:PYTHONPATH = $oldPythonPath
}
