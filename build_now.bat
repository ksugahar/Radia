@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
cd /d s:\Radia\01_GitHub\build-msvc
cmake --build . --config Release --target radia
if %errorlevel% neq 0 (
    echo Build failed with error %errorlevel%
    exit /b %errorlevel%
)
echo Build completed
copy /Y radia.cp312-win_amd64.pyd ..\src\radia\radia.pyd
