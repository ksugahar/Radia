@echo off
setlocal

set "VSPATH=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools"
set "VCVARS=%VSPATH%\VC\Auxiliary\Build\vcvarsall.bat"
set "SRC=S:\Radia\01_GitHub\src\cubit_plugin"
set "BUILD=S:\Radia\01_GitHub\src\cubit_plugin\build-test"
set "CUBIT_DIR=C:\Program Files\Coreform Cubit 2025.3\cmake"
rem Disable Netgen for fallback test
set "NETGEN_DIR="

echo === Setting up MSVC environment ===
call "%VCVARS%" amd64
if errorlevel 1 (
    echo ERROR: Failed to set up MSVC environment
    exit /b 1
)

echo.
echo === Configuring Cubit Plugin ===
cmake -G Ninja ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCubit_DIR="%CUBIT_DIR%" ^
  -DNETGEN_DIR="%NETGEN_DIR%" ^
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
