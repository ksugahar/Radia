@echo on
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
REM Auto-discover the latest Cubit install via the shared PowerShell helper.
REM Parses KEY=VALUE output from tools/find_cubit.ps1 and exposes CubitCmakeDir.
for /f "tokens=1,2 delims==" %%A in ('powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0..\..\tools\find_cubit.ps1"') do (
    set "%%A=%%B"
)
set "CUBIT_DIR=%CubitCmakeDir%"
if not defined CUBIT_DIR (
    echo ERROR: Cubit not found. Set CUBIT_INSTALL_DIR to override.
    exit /b 1
)
set "NETGEN_DIR=C:\Program Files\Python312\Lib\site-packages\netgen"
set "CMAKE_EXE=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
set "SRC=S:\Radia\01_GitHub\src\cubit_plugin"
mkdir "%SRC%\build-pyd" 2>/dev/null
cd /d "%SRC%\build-pyd"
"%CMAKE_EXE%" -G Ninja -DCMAKE_BUILD_TYPE=Release -DCubit_DIR="%CUBIT_DIR%" -DNETGEN_DIR="%NETGEN_DIR%" "%SRC%" && "%CMAKE_EXE%" --build . --config Release --target radia_cubit_mesh -j || exit /b 1
mkdir "%SRC%\build-ccm" 2>/dev/null
cd /d "%SRC%\build-ccm"
REM radia_cubit_ccl removed in radia 4.80.0 (Qt5 .ccl deleted; PySide6
REM toolbar at src/radia/panels/radia_export_menu.py replaces it).
"%CMAKE_EXE%" -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=cl -DCMAKE_CXX_COMPILER=cl -DCubit_DIR="%CUBIT_DIR%" -DNETGEN_DIR="%NETGEN_DIR%" "%SRC%" && "%CMAKE_EXE%" --build . --config Release --target radia_cubit_ccm -j || exit /b 1
exit /b 0
