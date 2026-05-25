$ErrorActionPreference = "Stop"

# Visual Studio + CMake discovery via vswhere (matches Build.ps1 +
# tools/_build_cubit_plugin.ps1 pattern). Was hardcoded
# `Microsoft Visual Studio\2022\BuildTools\...` until 2026-05-25;
# that broke on LAB whose VS BuildTools installs as
# `Microsoft Visual Studio\18\BuildTools\...` (the directory name uses
# the major-version "18", not marketing year "2022"). Auto-discover so
# future VS upgrades do not re-break this.
$VSWHERE = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
if (-not (Test-Path $VSWHERE)) { throw "vswhere.exe not found at $VSWHERE" }
$vsPath = & $VSWHERE -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath 2>$null
if (-not $vsPath) { throw "Visual Studio with VC x86/x64 tools not found via vswhere." }
$vcvarsBat = "$vsPath\VC\Auxiliary\Build\vcvars64.bat"
$cmake     = "$vsPath\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
if (-not (Test-Path $vcvarsBat)) { throw "vcvars64.bat not found at $vcvarsBat" }
if (-not (Test-Path $cmake)) {
    $pipCmake = & python -c "import shutil; print(shutil.which('cmake') or '')" 2>$null
    if ($pipCmake -and (Test-Path $pipCmake)) { $cmake = $pipCmake } else { throw "cmake.exe not found" }
}

# Auto-discover the latest Cubit install (was hardcoded 2025.3 until
# 2026-05-25; LAB now has 2025.12 only). Override via CUBIT_INSTALL_DIR.
. "$PSScriptRoot\..\..\tools\find_cubit.ps1"
if (-not $CubitCmakeDir) { throw "Cubit not found; set CUBIT_INSTALL_DIR" }
$env:CUBIT_DIR = $CubitCmakeDir
$env:NETGEN_DIR = 'C:\Program Files\Python312\Lib\site-packages\netgen'
$src = 'S:\Radia\01_GitHub\src\cubit_plugin'

# Pull vcvars env into pwsh
& cmd /c "`"$vcvarsBat`" && set" | ForEach-Object {
  if ($_ -match "^([^=]+)=(.*)$") { Set-Item -Path "env:$($matches[1])" -Value $matches[2] }
}

# build-pyd
$buildPyd = Join-Path $src "build-pyd"
if (-not (Test-Path $buildPyd)) { New-Item -ItemType Directory -Path $buildPyd | Out-Null }
Set-Location $buildPyd
& $cmake -G Ninja -DCMAKE_BUILD_TYPE=Release "-DCubit_DIR=$($env:CUBIT_DIR)" "-DNETGEN_DIR=$($env:NETGEN_DIR)" $src
if ($LASTEXITCODE -ne 0) { throw "cmake configure build-pyd failed" }
& $cmake --build . --config Release --target radia_cubit_mesh -j
if ($LASTEXITCODE -ne 0) { throw "cmake build radia_cubit_mesh failed" }

# build-ccm
$buildCcm = Join-Path $src "build-ccm"
if (-not (Test-Path $buildCcm)) { New-Item -ItemType Directory -Path $buildCcm | Out-Null }
Set-Location $buildCcm
& $cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=cl -DCMAKE_CXX_COMPILER=cl "-DCubit_DIR=$($env:CUBIT_DIR)" "-DNETGEN_DIR=$($env:NETGEN_DIR)" $src
if ($LASTEXITCODE -ne 0) { throw "cmake configure build-ccm failed" }
& $cmake --build . --config Release --target radia_cubit_ccm -j
if ($LASTEXITCODE -ne 0) { throw "cmake build radia_cubit_ccm failed" }
& $cmake --build . --config Release --target radia_cubit_ccl -j
if ($LASTEXITCODE -ne 0) { throw "cmake build radia_cubit_ccl failed" }

Write-Host "BUILD OK"
