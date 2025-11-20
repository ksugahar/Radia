@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\Tools\VsDevCmd.bat"
cd /d S:\radia\01_GitHub
cmake --build build --config Release --target rad_ngsolve
