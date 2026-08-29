param(
    [switch]$Rebuild
)

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
$ninja = & python -c "import shutil; print(shutil.which('ninja') or '')" 2>$null
if (-not $ninja -or -not (Test-Path $ninja)) {
    throw "ninja.exe not found on PATH or in the active Python environment."
}
$pythonExecutable = (Get-Command python -ErrorAction Stop).Source
$pybind11CmakeDir = (& $pythonExecutable -c "import pybind11; print(pybind11.get_cmake_dir())").Trim()
if (-not (Test-Path "$pybind11CmakeDir\pybind11Config.cmake")) {
    throw "pybind11 CMake package not found under $pybind11CmakeDir"
}

# Auto-discover the latest Cubit install (was hardcoded 2025.3 until
# 2026-05-25; LAB now has 2025.12 only). Override via CUBIT_INSTALL_DIR.
. "$PSScriptRoot\..\..\tools\find_cubit.ps1"
if (-not $CubitCmakeDir) { throw "Cubit not found; set CUBIT_INSTALL_DIR" }
$env:CUBIT_DIR = $CubitCmakeDir
$env:NETGEN_DIR = (& $pythonExecutable -c `
    "import netgen, os; print(os.path.dirname(netgen.__file__))").Trim()
if (-not (Test-Path "$($env:NETGEN_DIR)\include")) {
    throw "Netgen include directory not found under $($env:NETGEN_DIR)"
}

# Keep every path relative to this script. A hard-coded repository path makes
# a clean worktree silently build another checkout.
$src = $PSScriptRoot
$repoRoot = [System.IO.Path]::GetFullPath((Join-Path $src "..\.."))
$packageDir = Join-Path $repoRoot `
    "packages\cubit-mesh-export\src\cubit_mesh_export"

function Reset-BuildDirectory([string]$Path) {
    if (-not $Rebuild -or -not (Test-Path -LiteralPath $Path)) {
        return
    }
    $resolved = [System.IO.Path]::GetFullPath($Path)
    $pluginRoot = [System.IO.Path]::GetFullPath($src).TrimEnd('\') + '\'
    if (-not $resolved.StartsWith(
            $pluginRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to remove build directory outside ${pluginRoot}: $resolved"
    }
    Write-Host "Cleaning $resolved" -ForegroundColor Yellow
    Remove-Item -LiteralPath $resolved -Recurse -Force
}

# Pull vcvars env into pwsh
& cmd /c "`"$vcvarsBat`" && set" | ForEach-Object {
  if ($_ -match "^([^=]+)=(.*)$") { Set-Item -Path "env:$($matches[1])" -Value $matches[2] }
}

# build-pyd
$buildPyd = Join-Path $src "build-pyd"
Reset-BuildDirectory $buildPyd
if (-not (Test-Path $buildPyd)) { New-Item -ItemType Directory -Path $buildPyd | Out-Null }
Push-Location $buildPyd
try {
    & $cmake -G Ninja -DCMAKE_BUILD_TYPE=Release `
        "-DCMAKE_MAKE_PROGRAM=$ninja" `
        "-DPython3_EXECUTABLE=$pythonExecutable" `
        "-Dpybind11_DIR=$pybind11CmakeDir" `
        "-DCubit_DIR=$($env:CUBIT_DIR)" `
        "-DNETGEN_DIR=$($env:NETGEN_DIR)" `
        -DCOMPACT_NETGEN_OVERRIDES=NONE $src
    if ($LASTEXITCODE -ne 0) { throw "cmake configure build-pyd failed" }
    & $cmake --build . --config Release --target cubit_mesh_curver -j
    if ($LASTEXITCODE -ne 0) { throw "cmake build cubit_mesh_curver failed" }
} finally {
    Pop-Location
}

# build-ccm
$buildCcm = Join-Path $src "build-ccm"
Reset-BuildDirectory $buildCcm
if (-not (Test-Path $buildCcm)) { New-Item -ItemType Directory -Path $buildCcm | Out-Null }
Push-Location $buildCcm
try {
    & $cmake -G Ninja -DCMAKE_BUILD_TYPE=Release `
        -DCMAKE_C_COMPILER=cl -DCMAKE_CXX_COMPILER=cl `
        "-DCMAKE_MAKE_PROGRAM=$ninja" `
        "-DPython3_EXECUTABLE=$pythonExecutable" `
        "-Dpybind11_DIR=$pybind11CmakeDir" `
        "-DCubit_DIR=$($env:CUBIT_DIR)" `
        "-DNETGEN_DIR=$($env:NETGEN_DIR)" $src
    if ($LASTEXITCODE -ne 0) { throw "cmake configure build-ccm failed" }
    & $cmake --build . --config Release --target cubit_mesh_export_ccm -j
    if ($LASTEXITCODE -ne 0) { throw "cmake build cubit_mesh_export_ccm failed" }
} finally {
    Pop-Location
}
# cubit_mesh_export_ccl removed in radia 4.80.0 (Qt5 .ccl deleted; PySide6
# toolbar at src/radia/panels/radia_export_menu.py replaces it).

$pydOutput = Get-ChildItem -LiteralPath $buildPyd `
    -Filter "cubit_mesh_curver*.pyd" -File |
    Sort-Object LastWriteTime -Descending |
    Select-Object -First 1
$ccmOutput = Get-Item -LiteralPath `
    (Join-Path $buildCcm "cubit_mesh_export.ccm") -ErrorAction SilentlyContinue
if (-not $pydOutput -or -not $ccmOutput) {
    throw "Native Cubit build did not produce both cubit_mesh_curver.pyd and cubit_mesh_export.ccm"
}
if (-not (Test-Path $packageDir)) {
    throw "cubit-mesh-export package directory not found: $packageDir"
}

$deployments = @(
    @{ Source = $pydOutput.FullName; Destination = Join-Path $packageDir "cubit_mesh_curver.pyd" },
    @{ Source = $ccmOutput.FullName; Destination = Join-Path $packageDir "cubit_mesh_export.ccm" }
)
foreach ($deployment in $deployments) {
    Copy-Item -LiteralPath $deployment.Source `
        -Destination $deployment.Destination -Force
    $sourceHash = (Get-FileHash -LiteralPath $deployment.Source -Algorithm SHA256).Hash
    $destinationHash = (Get-FileHash -LiteralPath $deployment.Destination -Algorithm SHA256).Hash
    if ($sourceHash -ne $destinationHash) {
        throw "Binary propagation hash mismatch: $($deployment.Destination)"
    }
    Write-Host "  $($deployment.Destination)" -ForegroundColor Green
    Write-Host "    sha256=$destinationHash" -ForegroundColor DarkGray
}

Write-Host "BUILD AND PROPAGATION OK" -ForegroundColor Green
