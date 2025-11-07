@echo off
REM Build HACApK tests using CMake
REM Based on Radia's build system
REM
REM Usage:
REM   build.cmd [clean|rebuild]
REM
REM Options:
REM   clean   - Remove build directory
REM   rebuild - Clean and rebuild

setlocal enabledelayedexpansion

set BUILD_DIR=build
set BUILD_TYPE=Release

REM Parse command line arguments
if "%1"=="clean" (
	echo Cleaning build directory...
	if exist %BUILD_DIR% (
		rmdir /S /Q %BUILD_DIR%
		echo Build directory cleaned.
	)
	goto :EOF
)

if "%1"=="rebuild" (
	echo Rebuilding...
	if exist %BUILD_DIR% (
		rmdir /S /Q %BUILD_DIR%
	)
)

REM Find CMake
where cmake >nul 2>&1
if errorlevel 1 (
	echo ERROR: CMake not found in PATH
	echo.
	echo Please install CMake from: https://cmake.org/download/
	echo Or add CMake to your PATH environment variable
	exit /b 1
)

REM Create build directory
if not exist %BUILD_DIR% (
	echo Creating build directory...
	mkdir %BUILD_DIR%
)

echo.
echo ======================================================================
echo HACApK Tests Build
echo ======================================================================
echo.

REM Configure with CMake
echo Configuring with CMake...
echo.
cd %BUILD_DIR%
cmake .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_BUILD_TYPE=%BUILD_TYPE%

if errorlevel 1 (
	echo.
	echo [ERROR] CMake configuration failed
	echo.
	echo Troubleshooting:
	echo   1. Check that HACApK directory exists at:
	echo      S:\Radia\2025_10_31_HaCapK\HACApK_LH-Cimplm
	echo   2. Verify Visual Studio 2022 is installed
	echo   3. Check CMake output above for specific errors
	cd ..
	exit /b 1
)

echo.
echo ======================================================================
echo Building...
echo ======================================================================
echo.

REM Build
cmake --build . --config %BUILD_TYPE% --verbose

if errorlevel 1 (
	echo.
	echo [ERROR] Build failed
	echo.
	echo Check the error messages above for details
	cd ..
	exit /b 1
)

cd ..

echo.
echo ======================================================================
echo Build Complete
echo ======================================================================
echo.
echo Build type: %BUILD_TYPE%
echo Build directory: %BUILD_DIR%
echo.
echo Executables:
echo   %BUILD_DIR%\%BUILD_TYPE%\test_hacapk_basic.exe
echo   %BUILD_DIR%\%BUILD_TYPE%\test_hacapk_radia_concept.exe
echo.
echo To run tests:
echo   cd %BUILD_DIR%\%BUILD_TYPE%
echo   test_hacapk_basic.exe
echo   test_hacapk_radia_concept.exe
echo.
echo Or use CTest:
echo   cd %BUILD_DIR%
echo   ctest -C %BUILD_TYPE% --verbose
echo.
echo ======================================================================

endlocal
