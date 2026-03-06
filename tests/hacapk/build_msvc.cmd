@echo off
REM Manual build script using MSVC directly (no CMake required)
REM
REM Prerequisites: Visual Studio 2019 or later

setlocal

echo ======================================================================
echo Building HACApK Tests with MSVC
echo ======================================================================

REM Setup Visual Studio environment
REM Try to find vcvarsall.bat automatically
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"

if exist "%VSWHERE%" (
	for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do (
		set "VS_PATH=%%i"
	)
)

if exist "%VS_PATH%\VC\Auxiliary\Build\vcvars64.bat" (
	echo Setting up Visual Studio environment...
	call "%VS_PATH%\VC\Auxiliary\Build\vcvars64.bat"
) else (
	echo ERROR: Could not find Visual Studio installation
	echo Please run this from "x64 Native Tools Command Prompt for VS"
	exit /b 1
)

REM Create build directory
if not exist build mkdir build
cd build

REM Paths
set HACAPK_DIR=..\..\..\2025_10_31_HaCapK\HACApK_LH-Cimplm
set SRC_DIR=..

echo.
echo ======================================================================
echo Compiling HACApK C sources...
echo ======================================================================

REM Compile HACApK C files
cl /c /O2 /openmp /I"%HACAPK_DIR%" ^
	"%HACAPK_DIR%\cHACApK_base.c" ^
	"%HACAPK_DIR%\cHACApK_lib.c" ^
	"%HACAPK_DIR%\cHACApK_base_if.c"

if errorlevel 1 (
	echo [ERROR] HACApK compilation failed
	cd ..
	exit /b 1
)

echo [OK] HACApK sources compiled

REM Create static library
echo.
echo Creating HACApK library...
lib /OUT:hacapk.lib cHACApK_base.obj cHACApK_lib.obj cHACApK_base_if.obj

if errorlevel 1 (
	echo [ERROR] Library creation failed
	cd ..
	exit /b 1
)

echo [OK] hacapk.lib created

echo.
echo ======================================================================
echo Compiling test programs...
echo ======================================================================

REM Compile test_hacapk_basic
echo Compiling test_hacapk_basic.cpp...
cl /EHsc /O2 /std:c++17 /openmp /I"%HACAPK_DIR%" ^
	"%SRC_DIR%\test_hacapk_basic.cpp" ^
	hacapk.lib

if errorlevel 1 (
	echo [ERROR] test_hacapk_basic compilation failed
	cd ..
	exit /b 1
)

echo [OK] test_hacapk_basic.exe created

REM Compile test_hacapk_radia_concept
echo.
echo Compiling test_hacapk_radia_concept.cpp...
cl /EHsc /O2 /std:c++17 /openmp /I"%HACAPK_DIR%" ^
	"%SRC_DIR%\test_hacapk_radia_concept.cpp" ^
	hacapk.lib

if errorlevel 1 (
	echo [ERROR] test_hacapk_radia_concept compilation failed
	cd ..
	exit /b 1
)

echo [OK] test_hacapk_radia_concept.exe created

cd ..

echo.
echo ======================================================================
echo Build Complete
echo ======================================================================
echo.
echo Test executables created in build\ directory:
echo   - test_hacapk_basic.exe
echo   - test_hacapk_radia_concept.exe
echo.
echo To run tests:
echo   cd build
echo   test_hacapk_basic.exe
echo   test_hacapk_radia_concept.exe
echo.
echo ======================================================================

endlocal
