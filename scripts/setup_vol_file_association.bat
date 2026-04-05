@echo off
REM Setup Windows file association for .vol files
REM Run as Administrator to modify registry

echo ========================================
echo Setup .vol File Association for Netgen
echo ========================================
echo.

REM Get current directory
set "SCRIPT_DIR=%~dp0"
set "VIEWER_SCRIPT=%SCRIPT_DIR%netgen_vol_viewer.py"

echo Script location: %VIEWER_SCRIPT%
echo.

REM Check if script exists
if not exist "%VIEWER_SCRIPT%" (
    echo Error: netgen_vol_viewer.py not found!
    echo Expected location: %VIEWER_SCRIPT%
    pause
    exit /b 1
)

REM Find Python executable
where python >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo Error: Python not found in PATH!
    echo Please install Python and add to PATH.
    pause
    exit /b 1
)

echo Python found:
python --version
echo.

REM Registry setup (requires Administrator)
echo Setting up registry...
echo.
echo Creating file association: .vol = NetgenMeshFile
reg add "HKEY_CLASSES_ROOT\.vol" /ve /d "NetgenMeshFile" /f

echo.
echo Creating open command for NetgenMeshFile
reg add "HKEY_CLASSES_ROOT\NetgenMeshFile\shell\open\command" /ve /d "python.exe \"%VIEWER_SCRIPT%\" \"%%1\"" /f

echo.
echo ========================================
echo Setup Complete!
echo ========================================
echo.
echo You can now double-click .vol files to open in Netgen GUI.
echo.
echo To test:
echo   1. Navigate to a .vol file in Windows Explorer
echo   2. Double-click the file
echo   3. Netgen GUI should open automatically
echo.
echo Note: If this doesn't work, you may need to:
echo   - Run this script as Administrator
echo   - Restart Windows Explorer (or reboot)
echo.

pause
