[CmdletBinding(DefaultParameterSetName='Direct')]
param(
    [Parameter(Mandatory, ParameterSetName='Direct')][string]$Compiler,
    [Parameter(Mandatory, ParameterSetName='CMake')][string]$BuildDirectory,
    [string]$ReportPath
)
$ErrorActionPreference = 'Stop'
$cmakeVersion = $null
if ($PSCmdlet.ParameterSetName -eq 'CMake') {
    $records = @(Get-ChildItem -LiteralPath (Join-Path $BuildDirectory 'CMakeFiles') -Filter CMakeCXXCompiler.cmake -Recurse -File)
    if ($records.Count -ne 1) { throw 'Exactly one CMake C++ compiler identity is required.' }
    $record = Get-Content -LiteralPath $records[0].FullName -Raw
    $versionMatch = [regex]::Match($record, 'set\(CMAKE_CXX_COMPILER_VERSION "([^"\r\n]+)"\)')
    $compilerMatch = [regex]::Match($record, 'set\(CMAKE_CXX_COMPILER "([^"\r\n]+)"\)')
    if (-not $versionMatch.Success -or -not $compilerMatch.Success) { throw 'Missing CMake compiler/version identity.' }
    $cmakeVersion = $versionMatch.Groups[1].Value
    if ($cmakeVersion -notmatch '^19\.44\.') { throw "CMake used unsupported compiler $cmakeVersion" }
    $Compiler = $compilerMatch.Groups[1].Value
}
if (-not (Test-Path -LiteralPath $Compiler -PathType Leaf)) { throw 'Selected cl.exe is missing.' }
$oldPreference = $ErrorActionPreference
try {
    # /Bv is an identity query only; no source input and no object is produced.
    $ErrorActionPreference = 'Continue'
    $banner = (& $Compiler /Bv 2>&1 | Out-String)
} finally { $ErrorActionPreference = $oldPreference }
$actual = [regex]::Match($banner, '(?<![0-9])19\.44\.[0-9.]+')
if (-not $actual.Success) { throw 'Actual cl.exe banner is not MSVC 19.44.*.' }
$identity = [ordered]@{compiler_path=[IO.Path]::GetFullPath($Compiler);
    compiler_sha256=(Get-FileHash -LiteralPath $Compiler -Algorithm SHA256).Hash;
    compiler_banner_version=$actual.Value; cmake_compiler_version=$cmakeVersion}
$json = $identity | ConvertTo-Json
if ($ReportPath) {
    $parent = Split-Path -Parent ([IO.Path]::GetFullPath($ReportPath))
    [IO.Directory]::CreateDirectory($parent) | Out-Null
    $json | Set-Content -LiteralPath $ReportPath -Encoding utf8
}
$json
