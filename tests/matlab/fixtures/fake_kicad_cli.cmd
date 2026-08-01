@echo off
if not "%~1"=="sch" exit /b 2
if not "%~2"=="export" exit /b 3
if not "%~3"=="netlist" exit /b 4
> "%~7" echo KiCad fake SPICE export
>> "%~7" echo V1 in 0 1
>> "%~7" echo R1 in out 1k
>> "%~7" echo C1 out 0 1u
>> "%~7" echo .tran 0 1m
>> "%~7" echo .end
exit /b 0
