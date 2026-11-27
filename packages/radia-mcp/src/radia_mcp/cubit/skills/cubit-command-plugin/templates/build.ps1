# Build the plugin in Release with a pinned MSVC toolset.
#   pwsh -File build.ps1 [-CubitInstallDir <dir>] [-Toolset 14.44.35207] [-Rebuild]
# Output: <BuildDir>\my_cubit_plugin.dll
param(
    [string]$CubitInstallDir = $(if ($env:CUBIT_INSTALL_DIR) { $env:CUBIT_INSTALL_DIR } else { 'C:\Program Files\Coreform Cubit 2026.8' }),
    [string]$ExpectedRelease = '2026.8',
    [string]$Toolset = '14.44.35207',
    [string]$BuildDir = (Join-Path $PSScriptRoot 'build'),
    [switch]$Rebuild
)
$ErrorActionPreference = 'Stop'

# The release is read from BUILD_INFO, never from the directory name.
$buildInfo = Join-Path $CubitInstallDir 'BUILD_INFO'
if (-not (Test-Path -LiteralPath $buildInfo)) { throw "BUILD_INFO not found under $CubitInstallDir" }
$short = (Select-String -LiteralPath $buildInfo -Pattern '^short=(.+)$').Matches[0].Groups[1].Value.Trim()
if ($short -ne $ExpectedRelease) { throw "Cubit release is $short, expected $ExpectedRelease" }
$cubitCmake = Join-Path $CubitInstallDir 'cmake'

# Visual Studio 2022 (17.x) with the pinned toolset; no silent fallback.
$vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
$vs = (& $vswhere -latest -products '*' -version '[17,18)' `
    -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath | Select-Object -First 1)
if (-not $vs) { throw 'Visual Studio 2022 with C++ x64 tools not found' }
$cl = Join-Path $vs "VC\Tools\MSVC\$Toolset\bin\Hostx64\x64\cl.exe"
if (-not (Test-Path -LiteralPath $cl)) { throw "MSVC toolset $Toolset not installed" }
$ide = Join-Path $vs 'Common7\IDE\CommonExtensions\Microsoft\CMake'
$cmake = Join-Path $ide 'CMake\bin\cmake.exe'
$ninja = Join-Path $ide 'Ninja\ninja.exe'
foreach ($tool in $cmake, $ninja) { if (-not (Test-Path -LiteralPath $tool)) { throw "Missing $tool" } }

# Import the vcvars environment for that toolset.
& cmd /c "`"$vs\VC\Auxiliary\Build\vcvars64.bat`" -vcvars_ver=$Toolset >nul 2>&1 && set" | ForEach-Object {
    if ($_ -match '^([^=]+)=(.*)$') { Set-Item -Path "env:$($matches[1])" -Value $matches[2] }
}
if ([IO.Path]::GetFullPath((Get-Command cl.exe).Source) -ne [IO.Path]::GetFullPath($cl)) {
    throw "vcvars did not select toolset $Toolset"
}

if ($Rebuild -and (Test-Path -LiteralPath $BuildDir)) {
    # Delete only a real directory strictly inside the source tree: the
    # separator stops a sibling like "<src>2" or the source dir itself from
    # matching, and a reparse point is never followed by the recursive delete.
    $root = [IO.Path]::GetFullPath($PSScriptRoot).TrimEnd('\', '/') + [IO.Path]::DirectorySeparatorChar
    $full = [IO.Path]::GetFullPath($BuildDir).TrimEnd('\', '/')
    if (-not ($full + [IO.Path]::DirectorySeparatorChar).StartsWith($root, [StringComparison]::OrdinalIgnoreCase) -or
        ($full + [IO.Path]::DirectorySeparatorChar) -ieq $root) {
        throw "Refusing to delete a build directory not strictly inside ${root}: $full"
    }
    # Walk every component from the trusted root down to the build dir: a
    # junction or symlink anywhere on that chain could redirect the delete.
    $relative = $full.Substring($root.Length)
    $cursor = $root.TrimEnd('\', '/')
    foreach ($part in $relative.Split([char[]]@('\', '/'), [StringSplitOptions]::RemoveEmptyEntries)) {
        $cursor = Join-Path $cursor $part
        $node = Get-Item -LiteralPath $cursor -Force
        if (-not $node.PSIsContainer -or ($node.Attributes -band [IO.FileAttributes]::ReparsePoint)) {
            throw "Refusing to delete through a non-directory or reparse point: $cursor"
        }
    }
    if (Get-ChildItem -LiteralPath $full -Recurse -Force -Attributes ReparsePoint -ErrorAction Stop) {
        throw "Refusing to delete a build directory that contains reparse points: $full"
    }
    Remove-Item -LiteralPath $full -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $BuildDir | Out-Null

& $cmake -S $PSScriptRoot -B $BuildDir -G Ninja -DCMAKE_BUILD_TYPE=Release `
    "-DCMAKE_MAKE_PROGRAM=$ninja" "-DCMAKE_CXX_COMPILER=$cl" "-DCubit_DIR=$cubitCmake"
if ($LASTEXITCODE -ne 0) { throw 'cmake configure failed' }
& $cmake --build $BuildDir --config Release
if ($LASTEXITCODE -ne 0) { throw 'cmake build failed' }

$dll = Join-Path $BuildDir 'my_cubit_plugin.dll'
if (-not (Test-Path -LiteralPath $dll)) { throw "Build produced no $dll" }
# The loader finds plugins by this marker; a Debug build carries debug=yes.
$bytes = [IO.File]::ReadAllBytes($dll)
$text = [Text.Encoding]::ASCII.GetString($bytes)
if (-not $text.Contains('cubit_plugin_verification;debug=no')) {
    throw 'Plugin marker missing or Debug build: CUBIT_PLUGIN(...) not compiled in Release'
}
Write-Host "BUILT $dll sha256=$((Get-FileHash -LiteralPath $dll -Algorithm SHA256).Hash)"
