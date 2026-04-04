@echo off
setlocal

set "VSPATH=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools"
set "VCVARS=%VSPATH%\VC\Auxiliary\Build\vcvarsall.bat"
set "SRC=S:\Radia\01_GitHub\src\cubit_plugin"
set "BUILD=S:\Radia\01_GitHub\src\cubit_plugin\build-compact"
set "CUBIT_DIR=C:\Program Files\Coreform Cubit 2025.3\cmake"
set "NETGEN_SRC_DIR=C:\netgen_build\netgen_fork"
set "COMPACT_NETGEN_OVERRIDES=C:\compact_netgen_build\overrides"

echo === Setting up MSVC environment ===
call "%VCVARS%" amd64
if errorlevel 1 (
    echo ERROR: Failed to set up MSVC environment
    exit /b 1
)

where cl
if errorlevel 1 (
    echo ERROR: cl.exe not found
    exit /b 1
)

echo.
echo === Configuring Cubit Plugin (Compact Netgen - source integrated) ===
cmake -G Ninja ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_CXX_COMPILER=cl ^
  -DCubit_DIR="%CUBIT_DIR%" ^
  -DNETGEN_DIR= ^
  -DNETGEN_SRC_DIR="%NETGEN_SRC_DIR%" ^
  -DCOMPACT_NETGEN_OVERRIDES="%COMPACT_NETGEN_OVERRIDES%" ^
  -S "%SRC%" ^
  -B "%BUILD%"
if errorlevel 1 (
    echo ERROR: CMake configure failed
    exit /b 1
)

echo.
echo === Building Plugin ===
cmake --build "%BUILD%" --config Release --parallel
if errorlevel 1 (
    echo ERROR: Build failed
    exit /b 1
)

echo.
echo === Build successful ===
dir "%BUILD%\radia_cubit.ccm" 2>nul
dir "%BUILD%\radia_cubit.ccl" 2>nul
