#==============================================================================
# _build_cubit_plugin.ps1 — clean rebuild of radia_cubit.ccm + radia_cubit.ccl
#
# Used by tools/release_triple.py phase0. Skips the .pyd target (it requires
# pybind11 + non-compact-netgen which isn't the CI / dev path; the bundled
# .pyd in the repo is rebuilt manually when its source changes).
#
# Output:
#   src/cubit_plugin/build-ccm/radia_cubit.ccm
#   src/cubit_plugin/build-ccm/radia_cubit.ccl
#==============================================================================

$ErrorActionPreference = "Stop"

$vcvarsBat = 'C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat'
$env:CUBIT_DIR = 'C:\Program Files\Coreform Cubit 2025.3\cmake'
$env:NETGEN_DIR = 'C:\Program Files\Python312\Lib\site-packages\netgen'
$cmake = 'C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe'

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
    "-DCubit_DIR=$($env:CUBIT_DIR)" "-DNETGEN_DIR=$($env:NETGEN_DIR)" $src
if ($LASTEXITCODE -ne 0) { throw "cmake configure build-ccm failed" }

& $cmake --build . --config Release --target radia_cubit_ccm -j
if ($LASTEXITCODE -ne 0) { throw "radia_cubit_ccm build failed" }

& $cmake --build . --config Release --target radia_cubit_ccl -j
if ($LASTEXITCODE -ne 0) { throw "radia_cubit_ccl build failed" }

Write-Host "BUILD OK"
