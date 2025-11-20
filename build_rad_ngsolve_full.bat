@echo off
echo ========================================
echo  Building rad_ngsolve module
echo ========================================
echo.

echo [1/3] Loading Visual Studio environment...
call "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\Tools\VsDevCmd.bat"
if %ERRORLEVEL% neq 0 (
    echo [ERROR] Failed to load VS environment
    exit /b 1
)
echo       VS environment loaded
echo.

echo [2/3] Configuring with CMake...
cd /d S:\radia\01_GitHub
cmake -S . -B build -G "Visual Studio 17 2022"
if %ERRORLEVEL% neq 0 (
    echo [ERROR] CMake configuration failed
    exit /b 1
)
echo       Configuration complete
echo.

echo [3/3] Building rad_ngsolve...
cmake --build build --config Release --target rad_ngsolve
if %ERRORLEVEL% neq 0 (
    echo [ERROR] Build failed
    exit /b 1
)
echo       Build complete
echo.

echo ========================================
echo  Build SUCCESS
echo ========================================
echo.
echo Module: build\Release\rad_ngsolve.pyd
echo.
