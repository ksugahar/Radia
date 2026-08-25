#==============================================================================
# _build_cubit_plugin.ps1 — clean rebuild of cubit_mesh_export.ccm + cubit_mesh_export.ccl
#
# Used by tools/release_quad.py phase0. Skips the .pyd target (it requires
# pybind11 + non-compact-netgen which isn't the CI / dev path; the bundled
# .pyd in the repo is rebuilt manually when its source changes).
#
# Output:
#   src/cubit_plugin/build-ccm/cubit_mesh_export.ccm
#   src/cubit_plugin/build-ccm/cubit_mesh_export.ccl
#==============================================================================

$ErrorActionPreference = "Stop"
$pythonExecutable = (Get-Command python -ErrorAction Stop).Source
$pybind11CmakeDir = (& $pythonExecutable -c "import pybind11; print(pybind11.get_cmake_dir())").Trim()

# Visual Studio + CMake discovery via vswhere (matches Build.ps1's
# pattern).  Was hardcoded `Microsoft Visual Studio\2022\BuildTools\...`
# until 2026-05-25; that broke on LAB whose VS BuildTools installs as
# `Microsoft Visual Studio\18\BuildTools\...` (VS uses the major-version
# directory name "18", not the marketing year "2022").  Auto-discover so
# the next VS upgrade does not re-break this.
$VSWHERE = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
if (-not (Test-Path $VSWHERE)) {
    throw "vswhere.exe not found at $VSWHERE -- install Visual Studio Installer."
}
$vsPath = & $VSWHERE -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath 2>$null
if (-not $vsPath) {
    throw "Visual Studio with VC x86/x64 tools not found via vswhere."
}
$vcvarsBat = "$vsPath\VC\Auxiliary\Build\vcvars64.bat"
$cmake     = "$vsPath\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
if (-not (Test-Path $vcvarsBat)) { throw "vcvars64.bat not found at $vcvarsBat" }
if (-not (Test-Path $cmake)) {
    # Fallback: pip-installed cmake
    $pipCmake = & $pythonExecutable -c "import shutil; print(shutil.which('cmake') or '')" 2>$null
    if ($pipCmake -and (Test-Path $pipCmake)) { $cmake = $pipCmake }
    else { throw "cmake.exe not found in VS install or on PATH" }
}
$ninja = & $pythonExecutable -c "import shutil; print(shutil.which('ninja') or '')" 2>$null
if (-not $ninja -or -not (Test-Path $ninja)) {
    throw "ninja.exe not found on PATH or in the active Python environment."
}

# Auto-discover Cubit (was 2025.3, LAB now has 2025.12). Honors
# CUBIT_INSTALL_DIR env override.
. "$PSScriptRoot\find_cubit.ps1"
if (-not $CubitCmakeDir) {
    throw "Cubit not found under C:\Program Files. Set CUBIT_INSTALL_DIR to override."
}
$env:CUBIT_DIR = $CubitCmakeDir
$netgenPackageDir = (& $pythonExecutable -c "import netgen,os;print(os.path.dirname(netgen.__file__))").Trim()
if (-not (Test-Path "$netgenPackageDir\include")) {
    throw "Netgen include directory not found under $netgenPackageDir"
}
$env:NETGEN_DIR = $netgenPackageDir

# Resolve the cubit plugin source dir relative to this script.
#
# CRITICAL: keep the path on the mapped drive form (e.g. S:\...) and
# AVOID the UNC form (\\192.168.11.100\work\...).  CMake propagates
# the source path into try-compile working directories, and CMD.EXE
# refuses UNC pwds with "UNC paths are not supported. Defaulting to
# Windows directory.", killing the toolchain probe.
#
# (Resolve-Path).ProviderPath canonicalizes mapped drives back to UNC
# on this lab's setup -- DO NOT use it here.  Instead, look up which
# PSDrive backs $scriptDir and rewrite UNC -> drive letter explicitly
# if the script was invoked via UNC (e.g. by a Python caller that did
# Path(__file__).resolve()).
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
function Convert-UNCToMappedDrive($unc) {
  if ($unc -notmatch '^\\\\') { return $unc }
  foreach ($d in Get-PSDrive -PSProvider FileSystem) {
    if (-not $d.DisplayRoot) { continue }
    $root = $d.DisplayRoot.TrimEnd('\').ToLowerInvariant()
    $u = $unc.ToLowerInvariant()
    if ($u.StartsWith($root + '\') -or $u -eq $root) {
      return ($d.Name + ':' + $unc.Substring($root.Length))
    }
  }
  return $unc
}
$scriptDir = Convert-UNCToMappedDrive $scriptDir
$src = Join-Path $scriptDir "..\src\cubit_plugin"
$src = $src -replace '/', '\'

# Pull MSVC vars into pwsh.
& cmd /c "`"$vcvarsBat`" && set" | ForEach-Object {
  if ($_ -match "^([^=]+)=(.*)$") {
    Set-Item -Path "env:$($matches[1])" -Value $matches[2]
  }
}

$buildCcm = Join-Path $src "build-ccm"
if (-not (Test-Path $buildCcm)) { New-Item -ItemType Directory -Path $buildCcm | Out-Null }
Set-Location $buildCcm

& $cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=cl -DCMAKE_CXX_COMPILER=cl `
    "-DCMAKE_MAKE_PROGRAM=$ninja" `
    "-DPython3_EXECUTABLE=$pythonExecutable" `
    "-Dpybind11_DIR=$pybind11CmakeDir" `
    "-DCubit_DIR=$($env:CUBIT_DIR)" "-DNETGEN_DIR=$($env:NETGEN_DIR)" $src
if ($LASTEXITCODE -ne 0) { throw "cmake configure build-ccm failed" }

& $cmake --build . --config Release --target cubit_mesh_export_ccm -j
if ($LASTEXITCODE -ne 0) { throw "cubit_mesh_export_ccm build failed" }

# cubit_mesh_export_ccl removed in radia 4.80.0 (Qt5 GUI deleted; PySide6 toolbar at
# src/radia/panels/radia_export_menu.py replaces it).

Write-Host "BUILD OK"
