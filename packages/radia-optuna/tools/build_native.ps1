param(
    [Parameter(Mandatory)][ValidateSet('python','matlab')][string]$Stage,
    [Parameter(Mandatory)][string]$Python,
    [Parameter(Mandatory)][string]$MatlabRoot,
    [Parameter(Mandatory)][string]$BuildDirectory,
    [string]$PythonVerification
)
$ErrorActionPreference='Stop'
$vswhere=Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio/Installer/vswhere.exe'
$vs=(& $vswhere -products Microsoft.VisualStudio.Product.BuildTools -version '[17.0,18.0)' -latest -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath).Trim()
if(!$vs -or $vs -match '["&|<>\r\n]'){throw 'Invalid VS Build Tools2022 location'}
$tools=Join-Path $vs 'VC/Tools/MSVC/14.44.35207'
$cl=Join-Path $tools 'bin/HostX64/x64/cl.exe'
if(!(Test-Path -LiteralPath $cl)){throw 'MSVC14.44.35207 missing'}
$vcvars=Join-Path $vs 'VC/Auxiliary/Build/vcvarsall.bat'
# Only this isolated build process and its children receive the compiler env.
$lines=& $env:ComSpec /d /s /c "`"$vcvars`" x64 -vcvars_ver=14.44 >nul && set"
if($LASTEXITCODE -ne 0){throw 'Cannot initialize the approved host toolchain'}
foreach($line in $lines){if($line -match '^([^=]+)=(.*)$'){[Environment]::SetEnvironmentVariable($matches[1],$matches[2],'Process')}}
$builder=Join-Path $PSScriptRoot 'build_native.py'
$arguments=@($builder,'--stage',$Stage,'--matlab-root',$MatlabRoot,'--host-compiler',$cl,'--build-dir',$BuildDirectory)
if($PythonVerification){$arguments+=@('--python-verification',$PythonVerification)}
& $Python @arguments
if($LASTEXITCODE -ne 0){exit $LASTEXITCODE}
