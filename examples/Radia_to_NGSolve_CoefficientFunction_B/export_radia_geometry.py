"""
Export Radia Magnet Geometry to VTK Format

This script creates the Radia magnet geometry from visualize_field.py
and exports it as a VTK file using the existing radia_vtk_export module.

Usage:
	python export_radia_geometry.py

Output:
	radia_components.vtk - VTK file containing the magnet geometry


Date: 2025-11-01
"""

import sys
sys.path.insert(0, r"S:\radia\01_GitHub\build\Release")
sys.path.insert(0, r"S:\radia\01_GitHub\dist")
sys.path.insert(0, r"S:\radia\01_GitHub\src\python")

import radia as rad
from radia_vtk_export import exportGeometryToVTK

print("=" * 70)
print("Radia Geometry Export to VTK")
print("=" * 70)

# ============================================================================
# Step 1: Create Radia Magnet Geometry (from visualize_field.py)
# ============================================================================

print("\n[Step 1] Creating Radia Magnet Geometry")
print("-" * 70)

magnet_center = [0, 0, 0]
magnet_size = [20, 20, 30]

magnet = rad.ObjRecMag(magnet_center, magnet_size, [0, 0, 1.2])
rad.MatApl(magnet, rad.MatStd('NdFeB', 1.2))
rad.Solve(magnet, 0.0001, 10000)

print(f"Magnet created: object #{magnet}")
print(f"  Center: {magnet_center} mm")
print(f"  Size: {magnet_size} mm")
print(f"  Material: NdFeB, Br = 1.2 T")

# ============================================================================
# Step 2: Export to VTK using existing module
# ============================================================================

print("\n[Step 2] Exporting to VTK")
print("-" * 70)

exportGeometryToVTK(magnet, 'radia_components')

# ============================================================================
# Summary
# ============================================================================

print("\n" + "=" * 70)
print("Summary")
print("=" * 70)

print("\nOutput file: radia_components.vtk")
print("  Magnet: NdFeB rectangular magnet")
print("  Magnetization: [0, 0, 1.2] T")

print("\nVisualization:")
print("  1. Open Paraview")
print("  2. File -> Open -> radia_components.vtk")
print("  3. Click 'Apply' in the Properties panel")
print("  4. Color by 'Radia_colours' to see the magnet")

print("\n" + "=" * 70)
print("Complete")
print("=" * 70)
