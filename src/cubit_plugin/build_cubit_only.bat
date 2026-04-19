@echo on
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
set "CUBIT_DIR=C:\Program Files\Coreform Cubit 2025.3\cmake"
set "NETGEN_DIR=C:\Program Files\Python312\Lib\site-packages\netgen"
set "CMAKE_EXE=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
set "SRC=S:\Radia\01_GitHub\src\cubit_plugin"
mkdir "%SRC%\build-pyd" 2>/dev/null
cd /d "%SRC%\build-pyd"
"%CMAKE_EXE%" -G Ninja -DCMAKE_BUILD_TYPE=Release -DCubit_DIR="%CUBIT_DIR%" -DNETGEN_DIR="%NETGEN_DIR%" "%SRC%" && "%CMAKE_EXE%" --build . --config Release --target radia_cubit_mesh -j || exit /b 1
mkdir "%SRC%\build-ccm" 2>/dev/null
cd /d "%SRC%\build-ccm"
"%CMAKE_EXE%" -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=cl -DCMAKE_CXX_COMPILER=cl -DCubit_DIR="%CUBIT_DIR%" -DNETGEN_DIR="%NETGEN_DIR%" "%SRC%" && "%CMAKE_EXE%" --build . --config Release --target radia_cubit_ccm -j && "%CMAKE_EXE%" --build . --config Release --target radia_cubit_ccl -j || exit /b 1
exit /b 0
