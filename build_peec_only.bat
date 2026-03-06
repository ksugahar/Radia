@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
set LIB=C:\Program Files (x86)\Intel\oneAPI\mkl\latest\lib;%LIB%
set INCLUDE=C:\Program Files (x86)\Intel\oneAPI\mkl\latest\include;%INCLUDE%
set MKLROOT=C:\Program Files (x86)\Intel\oneAPI\mkl\latest
cd /d S:\Radia\01_GitHub\build-msvc
cmake --build . --config Release --target peec_matrices -j
