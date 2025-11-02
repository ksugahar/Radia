@echo off
REM Upload script for PyPI
REM
REM This script uploads the built distributions to PyPI or Test PyPI
REM
REM Prerequisites:
REM 1. Run python_build.cmd first to build distributions
REM 2. Install twine: pip install twine
REM 3. Configure PyPI credentials in ~/.pypirc or use environment variables

echo ========================================
echo Radia PyPI Upload Script
echo ========================================

REM Check if distribution files exist
if not exist "dist\*.whl" (
    echo ERROR: No wheel files found in dist/
    echo Please run python_build.cmd first
    exit /b 1
)

if not exist "dist\*.tar.gz" (
    echo ERROR: No source distribution found in dist/
    echo Please run python_build.cmd first
    exit /b 1
)

echo [OK] Found distribution files

REM Show what will be uploaded
echo.
echo Files to upload:
dir /b dist\*.whl dist\*.tar.gz

echo.
echo ========================================
echo Upload Options
echo ========================================
echo.
echo 1. Upload to Test PyPI (recommended first)
echo 2. Upload to PyPI (production)
echo 3. Check package with twine
echo 4. Cancel
echo.
set /p choice="Select option (1-4): "

if "%choice%"=="1" goto testpypi
if "%choice%"=="2" goto pypi
if "%choice%"=="3" goto check
if "%choice%"=="4" goto cancel
goto cancel

:check
echo.
echo ========================================
echo Checking Package with Twine
echo ========================================
python -m twine check dist/*
if errorlevel 1 (
    echo ERROR: Package check failed
    exit /b 1
)
echo [OK] Package check passed
goto end

:testpypi
echo.
echo ========================================
echo Uploading to Test PyPI
echo ========================================
echo.
echo URL: https://test.pypi.org/
echo.
python -m twine upload --repository testpypi dist/*
if errorlevel 1 (
    echo ERROR: Upload to Test PyPI failed
    exit /b 1
)
echo.
echo [OK] Upload to Test PyPI successful
echo.
echo To install from Test PyPI:
echo   pip install --index-url https://test.pypi.org/simple/ radia
goto end

:pypi
echo.
echo ========================================
echo WARNING: Uploading to Production PyPI
echo ========================================
echo.
echo This will upload to the PRODUCTION PyPI repository.
echo Once uploaded, you CANNOT delete or modify the release.
echo.
set /p confirm="Are you sure? Type 'yes' to continue: "
if not "%confirm%"=="yes" goto cancel

echo.
echo Uploading to PyPI...
python -m twine upload dist/*
if errorlevel 1 (
    echo ERROR: Upload to PyPI failed
    exit /b 1
)
echo.
echo [OK] Upload to PyPI successful
echo.
echo Package available at: https://pypi.org/project/radia/
echo.
echo To install:
echo   pip install radia
goto end

:cancel
echo.
echo Upload cancelled
exit /b 0

:end
echo.
echo ========================================
echo Done
echo ========================================
