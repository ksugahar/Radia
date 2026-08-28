@echo off
setlocal
cd /d "%~dp0.."
if not exist "build\Eqnedit64.exe" (
    echo [ERROR] build\Eqnedit64.exe is missing. Run build_eqnedt64.bat first.
    pause
    exit /b 1
)
start "" "build\Eqnedit64.exe" --debug-operations
