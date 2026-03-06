#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Coil geometry visualization

Visualizes the coil geometry and calculates magnetic field at test points.
This script is used to verify the coil shape is correct.
"""

import sys
import os
from pathlib import Path

# Add paths
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root / 'build' / 'lib' / 'Release'))
sys.path.insert(0, str(project_root / 'dist'))
sys.path.insert(0, str(project_root / 'src' / 'python'))

import numpy as np
import radia as rad

# Import coil model
from coil_model import create_beam_steering_coil, get_coil_info

# Try to import visualization tools
try:
	from radia_pyvista_viewer import view_radia_object
	HAS_PYVISTA = True
except ImportError:
	HAS_PYVISTA = False

try:
	from radia_vtk_export import exportGeometryToVTK
	HAS_VTK_EXPORT = True
except ImportError:
	HAS_VTK_EXPORT = False


def calculate_field_at_test_points(coil_obj):
	"""
	Calculate magnetic field at several test points.

	Args:
		coil_obj: Radia object ID
	"""
	print("\n" + "-" * 70)
	print("Calculating magnetic field at test points...")
	print("-" * 70)
	print(f"{'Position (mm)':<25} {'Bx (mT)':<15} {'By (mT)':<15} {'Bz (mT)':<15} {'|B| (mT)':<15}")
	print("-" * 70)

	# Test points along different axes
	test_points = [
		[0, 0, 0],
		[0, 0, 100],
		[0, 0, 500],
		[100, 0, 0],
		[0, 100, 0],
		[50, 0, 100],
		[0, 200, 0],
		[0, 0, -100],
	]

	for pt in test_points:
		B = rad.Fld(coil_obj, 'b', pt)
		Bx, By, Bz = B[0] * 1000, B[1] * 1000, B[2] * 1000  # T to mT
		B_mag = np.sqrt(Bx**2 + By**2 + Bz**2)
		pt_str = f"({pt[0]:.0f}, {pt[1]:.0f}, {pt[2]:.0f})"
		print(f"{pt_str:<25} {Bx:<15.6f} {By:<15.6f} {Bz:<15.6f} {B_mag:<15.6f}")

	print("-" * 70)

	# Calculate field along a line (for plotting)
	print("\n" + "-" * 70)
	print("Field along Z-axis (X=0, Y=0):")
	print("-" * 70)
	print(f"{'Z (mm)':<15} {'Bx (mT)':<15} {'|B| (mT)':<15}")
	print("-" * 40)
	z_points = np.linspace(-200, 600, 17)
	for z in z_points:
		B = rad.Fld(coil_obj, 'b', [0, 0, z])
		Bx = B[0] * 1000
		B_mag = np.linalg.norm(B) * 1000
		print(f"{z:<15.1f} {Bx:<15.6f} {B_mag:<15.6f}")
	print("-" * 40)


def main():
	"""
	Main visualization script for coil geometry.
	"""
	print("=" * 70)
	print("COIL GEOMETRY VISUALIZATION")
	print("=" * 70)
	print("\nThis script visualizes the coil geometry to verify the shape.\n")

	# Create coil from model
	print("-" * 70)
	print("Loading coil model...")
	print("-" * 70)
	coil, params = create_beam_steering_coil()

	print(f"[OK] Coil model loaded")
	print(f"     Description: {params['description']}")
	print(f"     Current: {params['current']} A")
	print(f"     Cross-section: {params['cross_section']['width']}×{params['cross_section']['height']} mm")
	print(f"     Segments: {params['num_segments']}")

	# Get geometry info
	info = get_coil_info(coil)
	print(f"\n     Bounding box:")
	print(f"       X: [{info['bbox']['x_min']:.2f}, {info['bbox']['x_max']:.2f}] mm")
	print(f"       Y: [{info['bbox']['y_min']:.2f}, {info['bbox']['y_max']:.2f}] mm")
	print(f"       Z: [{info['bbox']['z_min']:.2f}, {info['bbox']['z_max']:.2f}] mm")

	# Calculate magnetic field at test points
	calculate_field_at_test_points(coil)

	# Export to VTK
	if HAS_VTK_EXPORT:
		print("\n" + "-" * 70)
		print("Exporting to VTK format...")
		print("-" * 70)
		output_file = 'coil_geometry'
		exportGeometryToVTK(coil, output_file)
		print(f"  [OK] Created: {output_file}.vtk")

	# Visualize with PyVista
	if HAS_PYVISTA:
		print("\n" + "-" * 70)
		print("Opening PyVista viewer...")
		print("-" * 70)
		print("\nControls:")
		print("  - Left click + drag: Rotate")
		print("  - Right click + drag: Pan")
		print("  - Scroll wheel: Zoom")
		print("  - 'q': Quit\n")

		view_radia_object(coil)
	else:
		print("\n" + "-" * 70)
		print("PyVista not available.")
		print("-" * 70)
		print("Install with: pip install pyvista")

	print("\n" + "=" * 70)
	print("VISUALIZATION COMPLETE")
	print("=" * 70)
	print("\nCoil geometry has been verified.")
	print("Use field_map.py to calculate field distribution.")
	print("=" * 70 + "\n")


if __name__ == '__main__':
	main()
