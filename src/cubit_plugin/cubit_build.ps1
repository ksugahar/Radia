$ErrorActionPreference = "Stop"
# Enter MSVC env via cmd
$vcvarsBat = 'C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat'
$env:CUBIT_DIR = 'C:\Program Files\Coreform Cubit 2025.3\cmake'
$env:NETGEN_DIR = 'C:\Program Files\Python312\Lib\site-packages\netgen'
$cmake = 'C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe'
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
