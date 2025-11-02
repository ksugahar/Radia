@echo off
REM Build script for creating PyPI distribution packages
REM
REM This script builds both source distribution (sdist) and wheel distribution (bdist_wheel)
REM for uploading to PyPI.
REM
REM Prerequisites:
REM 1. Run Build.ps1 first to build radia.pyd
REM 2. (Optional) Run build_ngsolve.ps1 to build rad_ngsolve.pyd
REM 3. Install build tools: pip install build twine

echo ========================================
echo Radia PyPI Package Build Script
echo ========================================

REM Check if radia.pyd exists
if not exist "dist\radia.pyd" (
    echo ERROR: dist\radia.pyd not found
    echo Please run Build.ps1 first to build the core module
    exit /b 1
)

echo [OK] Found dist\radia.pyd

REM Check for rad_ngsolve.pyd (optional)
if exist "build\Release\rad_ngsolve.pyd" (
    echo [OK] Found build\Release\rad_ngsolve.pyd
) else (
    echo [INFO] rad_ngsolve.pyd not found (optional)
)

REM Clean previous builds
echo.
echo Cleaning previous builds...
if exist "build\lib" rmdir /s /q "build\lib"
if exist "build\bdist.*" rmdir /s /q "build\bdist.*"
if exist "dist\*.whl" del /q "dist\*.whl"
if exist "dist\*.tar.gz" del /q "dist\*.tar.gz"
if exist "src\python\radia.egg-info" rmdir /s /q "src\python\radia.egg-info"
if exist "radia.egg-info" rmdir /s /q "radia.egg-info"

REM Build source distribution
echo.
echo ========================================
echo Building Source Distribution (sdist)
echo ========================================
python -m build --sdist

if errorlevel 1 (
    echo ERROR: Source distribution build failed
    exit /b 1
)

echo [OK] Source distribution created

REM Build wheel distribution
echo.
echo ========================================
echo Building Wheel Distribution (bdist_wheel)
echo ========================================
python -m build --wheel

if errorlevel 1 (
    echo ERROR: Wheel distribution build failed
    exit /b 1
)

echo [OK] Wheel distribution created

REM Show results
echo.
echo ========================================
echo Build Complete
echo ========================================
echo.
echo Distribution files created in dist/:
dir /b dist\*.whl dist\*.tar.gz

echo.
echo ========================================
echo Next Steps
echo ========================================
echo.
echo To test the package locally:
echo   pip install dist\radia-ngsolve-1.0.0-py3-none-any.whl
echo.
echo To upload to Test PyPI:
echo   python -m twine upload --repository testpypi dist/*
echo.
echo To upload to PyPI:
echo   python -m twine upload dist/*
echo.
echo ========================================
