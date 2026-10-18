@echo off
REM Deliberately no setlocal: vcvars must configure the calling build process.
REM Never trust an inherited cl.exe or fall back to Visual Studio 2026.
set "EQN_VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
set "EQN_VSROOT="
if not exist "%EQN_VSWHERE%" (
    echo [ERROR] vswhere.exe is required for Visual Studio 2022 discovery.
    exit /b 1
)
for /f "delims=" %%I in ('powershell -NoProfile -Command "& $env:EQN_VSWHERE -latest -products '*' -version '[17.0,18.0)' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath"') do set "EQN_VSROOT=%%I"
if not defined EQN_VSROOT (
    echo [ERROR] Visual Studio 2022 C++ tools are required. No older fallback.
    exit /b 1
)
if not exist "%EQN_VSROOT%\VC\Auxiliary\Build\vcvars64.bat" (
    echo [ERROR] vcvars64.bat is missing from the selected VS2022 installation.
    exit /b 1
)
call "%EQN_VSROOT%\VC\Auxiliary\Build\vcvars64.bat" -vcvars_ver=14.44.35207 >nul 2>&1
if errorlevel 1 (
    echo [ERROR] vcvars64.bat failed to initialize MSVC 14.44.35207.
    exit /b 1
)
if not "%VisualStudioVersion%"=="17.0" (
    echo [ERROR] The selected developer environment is not Visual Studio 2022.
    exit /b 1
)
if not "%VCToolsVersion%"=="14.44.35207" (
    echo [ERROR] MSVC 14.44.35207 is required. No alternate toolset.
    exit /b 1
)
set "EQN_FIRST_CL="
for /f "delims=" %%I in ('where cl.exe 2^>nul') do if not defined EQN_FIRST_CL set "EQN_FIRST_CL=%%I"
if /I not "%EQN_FIRST_CL%"=="%VCToolsInstallDir%bin\Hostx64\x64\cl.exe" (
    echo [ERROR] PATH compiler differs from the selected VS2022 x64 tools.
    exit /b 1
)
echo [INFO] VS2022 compiler: %EQN_FIRST_CL%
pwsh -NoProfile -File "%~dp0assert_msvc1944.ps1" -Compiler "%EQN_FIRST_CL%"
if errorlevel 1 exit /b 1
exit /b 0
