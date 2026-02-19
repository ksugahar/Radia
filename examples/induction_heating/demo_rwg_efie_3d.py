"""
RWG-EFIE 3D Induction Heating Demo

This example demonstrates the RWG-EFIE (Rao-Wilton-Glisson Electric Field
Integral Equation) solver for induction heating analysis with 3D surface
elements.

The RWG-EFIE solver has been migrated from Python to C++ with OpenMP
parallelization for improved performance (~20x speedup).

Key Features Demonstrated:
1. Loop coil mesh generation (RwgMeshLoop)
2. Spiral coil mesh generation (RwgMeshSpiral)
3. Plate workpiece mesh generation (RwgMeshRect)
4. Disk workpiece mesh generation (RwgMeshDisk)
5. Cylindrical workpiece mesh generation (RwgMeshCylinder)
6. Coupled solver setup
7. Impedance calculation
8. Magnetic field computation

Author: Radia Development Team
Date: 2026-01-09
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import radia as rad
import numpy as np

rad.FldUnits('m')


def demo_mesh_creation():
    """Demonstrate various RWG mesh creation functions."""
    print("=" * 60)
    print("Demo 1: RWG Mesh Creation")
    print("=" * 60)
    print()

    # Create different mesh types
    print("Creating various RWG meshes:")
    print()

    # Loop coil
    loop = rad.RwgMeshLoop([0, 0, 0], 0.05, [0, 0, 1], 0.002, 8, 24)
    print(f"  Loop coil mesh handle: {loop}")

    # Spiral coil
    spiral = rad.RwgMeshSpiral([0, 0, 0], 0.02, 0.05, 0.005, 3, [0, 0, 1], 0.002, 8)
    print(f"  Spiral coil mesh handle: {spiral}")

    # Rectangular plate
    rect = rad.RwgMeshRect([0, 0, -0.01], 0.1, 0.1, [0, 0, 1], 10, 10)
    print(f"  Rectangular plate mesh handle: {rect}")

    # Circular disk
    disk = rad.RwgMeshDisk([0, 0, -0.01], 0.05, [0, 0, 1], 5, 16)
    print(f"  Circular disk mesh handle: {disk}")

    # Cylindrical shell
    cylinder = rad.RwgMeshCylinder([0, 0, 0], 0.02, 0.05, [0, 0, 1], 5, 12)
    print(f"  Cylindrical shell mesh handle: {cylinder}")

    print()
    print("All mesh types created successfully!")
    print()


def demo_coupled_solver():
    """Demonstrate the coupled EFIE solver for induction heating."""
    print("=" * 60)
    print("Demo 2: Coupled EFIE Solver for Induction Heating")
    print("=" * 60)
    print()

    # Create coil mesh (spiral)
    print("Creating spiral coil mesh...")
    coil = rad.RwgMeshSpiral(
        [0, 0, 0],      # center
        0.02,           # inner radius 20mm
        0.04,           # outer radius 40mm
        0.005,          # pitch 5mm
        3,              # 3 turns
        [0, 0, 1],      # axis
        0.002,          # wire radius 2mm
        8               # 8 panels around wire
    )
    print(f"  Coil mesh handle: {coil}")

    # Create workpiece mesh (disk below coil)
    print("Creating disk workpiece mesh...")
    workpiece = rad.RwgMeshDisk(
        [0, 0, -0.01],  # center (10mm below coil)
        0.05,           # radius 50mm
        [0, 0, 1],      # normal
        5,              # 5 radial divisions
        16              # 16 angular divisions
    )
    print(f"  Workpiece mesh handle: {workpiece}")

    # Create solver
    print("Creating coupled EFIE solver...")
    solver = rad.RwgSolverCreate()
    print(f"  Solver handle: {solver}")

    # Set meshes
    rad.RwgSetCoilMesh(solver, coil)
    rad.RwgSetWorkpieceMesh(solver, workpiece)

    # Set material properties
    freq = 50000  # 50 kHz
    rad.RwgSetFrequency(solver, freq)
    rad.RwgSetCoilConductivity(solver, 5.8e7)      # Copper: 5.8e7 S/m
    rad.RwgSetWorkpieceConductivity(solver, 5e6)   # Steel: 5e6 S/m
    rad.RwgSetWorkpiecePermeability(solver, 100)   # Steel mu_r = 100

    print()
    print("Material Properties:")
    print(f"  Frequency: {freq/1000:.1f} kHz")
    print(f"  Coil conductivity: 5.8e7 S/m (copper)")
    print(f"  Workpiece conductivity: 5e6 S/m (steel)")
    print(f"  Workpiece relative permeability: 100")

    # Set voltage excitation (1V)
    rad.RwgSetVoltage(solver, 1.0, 0.0)

    # Solve
    print()
    print("Solving coupled EFIE system (OpenMP parallelized)...")
    rad.RwgSolve(solver)
    print("  Solve completed!")

    # Get results
    Z = rad.RwgGetImpedance(solver)
    P = rad.RwgGetWorkpiecePower(solver)

    print()
    print("Results:")
    print(f"  Impedance: Z = {Z[0]*1000:.4f} + j*{Z[1]*1000:.4f} mOhm")
    print(f"  |Z| = {np.sqrt(Z[0]**2 + Z[1]**2)*1000:.4f} mOhm")
    print(f"  Workpiece power: {P*1000:.6f} mW")

    # Compute B field along z-axis
    print()
    print("Magnetic field along z-axis:")
    print(f"  {'z [mm]':>10} | {'|B| [uT]':>12}")
    print("  " + "-" * 26)

    z_points = [0.01, 0.02, 0.03, 0.04, 0.05]
    for z in z_points:
        B = rad.RwgComputeB(solver, [0, 0, z])
        B_mag = np.sqrt(B[0]**2 + B[1]**2 + B[2]**2) * 1e6  # Convert to uT
        print(f"  {z*1000:>10.1f} | {B_mag:>12.3f}")

    print()


def demo_frequency_sweep():
    """Demonstrate impedance variation with frequency."""
    print("=" * 60)
    print("Demo 3: Frequency Sweep")
    print("=" * 60)
    print()

    # Create simple loop coil and disk workpiece
    coil = rad.RwgMeshLoop([0, 0, 0], 0.04, [0, 0, 1], 0.002, 6, 20)
    workpiece = rad.RwgMeshDisk([0, 0, -0.01], 0.05, [0, 0, 1], 4, 12)

    solver = rad.RwgSolverCreate()
    rad.RwgSetCoilMesh(solver, coil)
    rad.RwgSetWorkpieceMesh(solver, workpiece)
    rad.RwgSetCoilConductivity(solver, 5.8e7)
    rad.RwgSetWorkpieceConductivity(solver, 5e6)
    rad.RwgSetWorkpiecePermeability(solver, 100)
    rad.RwgSetVoltage(solver, 1.0, 0.0)

    frequencies = [10e3, 20e3, 50e3, 100e3, 200e3]

    print("Loop coil + disk workpiece:")
    print()
    print(f"  {'f [kHz]':>10} | {'R [mOhm]':>12} | {'X [mOhm]':>12} | {'P_wp [mW]':>12}")
    print("  " + "-" * 56)

    for freq in frequencies:
        rad.RwgSetFrequency(solver, freq)
        rad.RwgSolve(solver)
        Z = rad.RwgGetImpedance(solver)
        P = rad.RwgGetWorkpiecePower(solver)
        print(f"  {freq/1000:>10.0f} | {Z[0]*1000:>12.4f} | {Z[1]*1000:>12.4f} | {P*1000:>12.6f}")

    print()


if __name__ == "__main__":
    print()
    print("RWG-EFIE 3D Induction Heating Demo")
    print("C++ implementation with OpenMP parallelization")
    print("=" * 60)
    print()

    demo_mesh_creation()
    demo_coupled_solver()
    demo_frequency_sweep()

    print("All demos completed successfully!")
    print()
