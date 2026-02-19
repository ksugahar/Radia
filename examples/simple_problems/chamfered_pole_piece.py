#!/usr/bin/env python
"""
Case 2: Multiple Extrusion with Chamfer - Pole Piece Magnet
Converted from Mathematica/Wolfram Language to Python

This example demonstrates:
- Creating complex 3D geometry using ObjMltExtRtg (multiple extrusion)
- Applying chamfer to edges
- Subdividing geometry for accurate field calculation
- Calculating magnetic field from extruded magnet
"""

import sys
import os
import math

# Add parent directory to path to import radia
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'dist'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'lib', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "build", "Release"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src', 'radia'))

import radia as rad

# Clear all objects
rad.UtiDelAll()
rad.FldUnits('m')

print("=" * 70)
print("Case 2: Chamfered Pole Piece using Multiple Extrusion")
print("=" * 70)

# Geometry parameters (all dimensions in meters)
gap = 0.010        # 10 mm gap between pole pieces
thick = 0.050      # 50 mm total thickness
width = 0.040      # 40 mm width
chamfer = 0.008    # 8 mm chamfer size
lz1 = 0.020        # 20 mm height of pole piece

# Magnetization (T)
magnetization = [0, 0, 1.0]  # 1.0 T in Z direction

print(f"\nGeometry parameters:")
print(f"  Gap: {gap*1000:.0f} mm")
print(f"  Thickness: {thick*1000:.0f} mm")
print(f"  Width: {width*1000:.0f} mm")
print(f"  Chamfer: {chamfer*1000:.0f} mm")
print(f"  Height: {lz1*1000:.0f} mm")
print(f"  Magnetization: {magnetization} T")

# Define cross-sections for multiple extrusion
# This creates a pole piece with chamfered edges
# Each section is defined as [[x, y, z], [dx, dy]]

# k1: first chamfered section (narrower due to chamfer)
k1 = [
	[thick/4 - chamfer/2, 0, gap/2],
	[thick/2 - chamfer, width - 2*chamfer]
]

# k2: second section at chamfer height (full width)
k2 = [
	[thick/4, 0, gap/2 + chamfer],
	[thick/2, width]
]

# k3: final section (same width as k2)
k3 = [
	[thick/4, 0, gap/2 + lz1],
	[thick/2, width]
]

print(f"\nCross-section levels:")
print(f"  Level 1 (z={gap/2*1000:.0f} mm): size=[{(thick/2 - chamfer)*1000:.0f}, {(width - 2*chamfer)*1000:.0f}] mm (chamfered)")
print(f"  Level 2 (z={(gap/2 + chamfer)*1000:.0f} mm): size=[{thick/2*1000:.0f}, {width*1000:.0f}] mm (full)")
print(f"  Level 3 (z={(gap/2 + lz1)*1000:.0f} mm): size=[{thick/2*1000:.0f}, {width*1000:.0f}] mm (full)")

# Create multiple extrusion object with magnetization
g1 = rad.ObjMltExtRtg([k1, k2, k3], magnetization)

print(f"\nMagnet object created: ID = {g1}")
print("Note: ObjMltExtRtg creates a single element (no subdivision)")

# Calculate magnetic field at various points
print("\n" + "=" * 70)
print("Magnetic Field Calculation")
print("=" * 70)

test_points = [
	[0, 0, gap/2 + lz1/2],            # Center of pole piece
	[0, 0, gap/2 + lz1 + 0.010],      # 10mm above pole piece
	[thick/2 + 0.005, 0, gap/2 + lz1/2],  # 5mm to the side
	[0, width/2 + 0.005, gap/2 + lz1/2],  # 5mm to front
]

print(f"\n{'Point (mm)':<30} {'Bx (mT)':<12} {'By (mT)':<12} {'Bz (mT)':<12} {'|B| (mT)':<12}")
print("-" * 80)

for point in test_points:
	field = rad.Fld(g1, 'b', point)
	Bx_mT = field[0] * 1000
	By_mT = field[1] * 1000
	Bz_mT = field[2] * 1000
	B_mag = math.sqrt(Bx_mT**2 + By_mT**2 + Bz_mT**2)

	point_str = f"({point[0]*1000:6.1f}, {point[1]*1000:6.1f}, {point[2]*1000:6.1f})"
	print(f"{point_str:<30} {Bx_mT:<12.3f} {By_mT:<12.3f} {Bz_mT:<12.3f} {B_mag:<12.3f}")

print("\n" + "=" * 70)
print("Note: This example demonstrates multiple extrusion (ObjMltExtRtg)")
print("      The chamfer creates a smooth transition from narrow to full width")
print("=" * 70)
print("Calculation complete.")
print("=" * 70)

# VTS Export - Export field distribution with same filename as script
try:
	# Get script basename without extension
	script_name = os.path.splitext(os.path.basename(__file__))[0]
	vts_filename = f"{script_name}.vts"
	vts_path = os.path.join(os.path.dirname(__file__), vts_filename)

	# Geometry: pole piece with gap=10mm, thick=50mm, width=40mm, lz1=20mm
	# Extend range to cover geometry with margin
	x_range = [-0.040, 0.040]
	y_range = [-0.040, 0.040]
	z_range = [0, 0.050]

	rad.FldVTS(g1, vts_path, x_range, y_range, z_range, 21, 21, 21, 1, 0, 1.0)
	print(f"\n[VTS] Exported: {vts_filename}")
	print(f"      View with: paraview {vts_filename}")
except Exception as e:
	print(f"\n[VTS] Warning: Export failed: {e}")
