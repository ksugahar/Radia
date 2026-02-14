@echo off
REM Generate Circular Coil Mesh with Port Definitions using Coreform Cubit

echo ======================================================================
echo Generating Circular Coil Mesh with Ports
echo ======================================================================
echo.

"C:\Program Files\Coreform Cubit 2025.3\bin\python3\python.exe" generate_coil_with_ports.py

echo.
echo ======================================================================
echo Mesh generation completed!
echo ======================================================================
echo.
echo Next steps:
echo   1. View mesh: gmsh circular_coil_with_ports.msh
echo   2. Run PEEC: cd .. ^& python demo_peec_dc.py
echo.

pause
