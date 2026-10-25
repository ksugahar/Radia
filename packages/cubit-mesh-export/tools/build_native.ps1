param(
    [switch]$Rebuild,
    [ValidateRange(1, 52)][int]$BuildJobs = 4
)

$ErrorActionPreference = "Stop"
$repoRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot "..\..\.."))

# Select only Visual Studio Build Tools 2022 (major 17), with no older compiler fallback.
function Find-CubitVisualStudio2022([string]$VswherePath) {
    $instances = & $VswherePath -latest -products '*' -version '[17,18)' `
        -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -format json
    if ($LASTEXITCODE -ne 0) { throw "Visual Studio Build Tools 2022 discovery failed via vswhere." }
    $selected = @($instances | ConvertFrom-Json)
    if ($selected.Count -ne 1 -or
        -not $selected[0].installationPath -or
        ([version]$selected[0].installationVersion).Major -ne 17) {
        throw "Visual Studio Build Tools 2022 (17.x) with VC x86/x64 tools is required; no older-version fallback."
    }
    $requiredCompiler = Join-Path $selected[0].installationPath 'VC\Tools\MSVC\14.44.35207\bin\Hostx64\x64\cl.exe'
    if (-not (Test-Path -LiteralPath $requiredCompiler)) {
        throw "Required MSVC toolset 14.44.35207 is missing; no alternate toolset fallback."
    }
    return $selected[0].installationPath
}
$VSWHERE = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
if (-not (Test-Path $VSWHERE)) { throw "vswhere.exe not found at $VSWHERE" }
$vsPath = Find-CubitVisualStudio2022 $VSWHERE
$vcvarsBat = "$vsPath\VC\Auxiliary\Build\vcvars64.bat"
$cmake     = "$vsPath\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
if (-not (Test-Path $vcvarsBat)) { throw "vcvars64.bat not found at $vcvarsBat" }
if (-not (Test-Path $cmake)) {
    $pipCmake = & python -c "import shutil; print(shutil.which('cmake') or '')" 2>$null
    if ($pipCmake -and (Test-Path $pipCmake)) { $cmake = $pipCmake } else { throw "cmake.exe not found" }
}
$ctest = Join-Path (Split-Path -Parent $cmake) "ctest.exe"
if (-not (Test-Path -LiteralPath $ctest)) { throw "ctest.exe not found beside selected CMake" }
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
. "$repoRoot\tools\find_cubit.ps1"
if (-not $CubitCmakeDir) { throw "Cubit not found; set CUBIT_INSTALL_DIR" }
$env:CUBIT_DIR = $CubitCmakeDir
$env:NETGEN_DIR = (& $pythonExecutable -c `
    "import netgen, os; print(os.path.dirname(netgen.__file__))").Trim()
if (-not (Test-Path "$($env:NETGEN_DIR)\include")) {
    throw "Netgen include directory not found under $($env:NETGEN_DIR)"
}

# Keep every path relative to this script. A hard-coded repository path makes
# a clean worktree silently build another checkout.
$src = Join-Path $repoRoot "src\cubit_plugin"


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
& cmd /c "`"$vcvarsBat`" -vcvars_ver=14.44.35207 && set" | ForEach-Object {
  if ($_ -match "^([^=]+)=(.*)$") { Set-Item -Path "env:$($matches[1])" -Value $matches[2] }
}

if ($LASTEXITCODE -ne 0) { throw "Visual Studio Build Tools 2022 vcvars64 initialization failed." }
$compilerPath = (Get-Command cl.exe -ErrorAction Stop).Source
$requiredCompiler = Join-Path $vsPath 'VC\Tools\MSVC\14.44.35207\bin\Hostx64\x64\cl.exe'
if ([System.IO.Path]::GetFullPath($compilerPath) -ne [System.IO.Path]::GetFullPath($requiredCompiler)) {
    throw "cl.exe does not match the pinned MSVC toolset 14.44.35207."
}

function Assert-CubitCompilerCache([string]$Directory, [string]$Compiler) {
    $cache = Join-Path $Directory 'CMakeCache.txt'
    if (-not (Test-Path -LiteralPath $cache)) { return }
    foreach ($line in Get-Content -LiteralPath $cache) {
        if ($line -match '^CMAKE_C(?:XX)?_COMPILER:[^=]+=(.+)$') {
            if ([System.IO.Path]::GetFullPath($matches[1]) -ne [System.IO.Path]::GetFullPath($Compiler)) {
                throw "CMake compiler cache differs from the pinned toolchain; rebuild explicitly with -Rebuild."
            }
        }
    }
}

# build-pyd
$buildPyd = Join-Path $src "build-pyd"
if (-not $Rebuild) { Assert-CubitCompilerCache $buildPyd $compilerPath }
Reset-BuildDirectory $buildPyd
if (-not (Test-Path $buildPyd)) { New-Item -ItemType Directory -Path $buildPyd | Out-Null }
Push-Location $buildPyd
try {
    & $cmake -G Ninja -DCMAKE_BUILD_TYPE=Release `
        "-DCMAKE_C_COMPILER=$compilerPath" "-DCMAKE_CXX_COMPILER=$compilerPath" `
        "-DCMAKE_MAKE_PROGRAM=$ninja" `
        "-DPython3_EXECUTABLE=$pythonExecutable" `
        "-Dpybind11_DIR=$pybind11CmakeDir" `
        "-DCubit_DIR=$($env:CUBIT_DIR)" `
        "-DNETGEN_DIR=$($env:NETGEN_DIR)" `
        -DCOMPACT_NETGEN_OVERRIDES=NONE -DCUBIT_CALLBACK_CONTRACT_TESTS=ON $src
    if ($LASTEXITCODE -ne 0) { throw "cmake configure build-pyd failed" }
    & $cmake --build . --config Release --target cubit_mesh_curver -j $BuildJobs
    if ($LASTEXITCODE -ne 0) { throw "cmake build cubit_mesh_curver failed" }
    & $cmake --build . --config Release --target cubit_callback_contract -j $BuildJobs
    if ($LASTEXITCODE -ne 0) { throw "cmake build cubit_callback_contract failed" }
    & $ctest --test-dir $buildPyd --output-on-failure
    if ($LASTEXITCODE -ne 0) { throw "native Cubit callback contract failed" }
} finally {
    Pop-Location
}

# build-ccm
$buildCcm = Join-Path $src "build-ccm"
if (-not $Rebuild) { Assert-CubitCompilerCache $buildCcm $compilerPath }
Reset-BuildDirectory $buildCcm
if (-not (Test-Path $buildCcm)) { New-Item -ItemType Directory -Path $buildCcm | Out-Null }
Push-Location $buildCcm
try {
    & $cmake -G Ninja -DCMAKE_BUILD_TYPE=Release `
        "-DCMAKE_C_COMPILER=$compilerPath" "-DCMAKE_CXX_COMPILER=$compilerPath" `
        "-DCMAKE_MAKE_PROGRAM=$ninja" `
        "-DPython3_EXECUTABLE=$pythonExecutable" `
        "-Dpybind11_DIR=$pybind11CmakeDir" `
        "-DCubit_DIR=$($env:CUBIT_DIR)" `
        "-DNETGEN_DIR=$($env:NETGEN_DIR)" $src
    if ($LASTEXITCODE -ne 0) { throw "cmake configure build-ccm failed" }
    & $cmake --build . --config Release --target cubit_mesh_export_ccm -j $BuildJobs
    if ($LASTEXITCODE -ne 0) { throw "cmake build cubit_mesh_export_ccm failed" }
} finally {
    Pop-Location
}

foreach ($output in @((Join-Path $buildPyd "cubit_mesh_curver.cp312-win_amd64.pyd"),
                       (Join-Path $buildCcm "cubit_mesh_export.ccm"))) {
    if (-not (Test-Path -LiteralPath $output -PathType Leaf)) { throw "Missing raw native output: $output" }
    Get-FileHash -LiteralPath $output -Algorithm SHA256 | Format-List
}
Write-Host "RAW BUILD AND CTEST OK - signed payloads require reviewed release promotion"
