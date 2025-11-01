#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Magnetic field map calculation and VTK export

This script calculates the magnetic field distribution around the complex
8-segment coil geometry and exports it to VTK format for visualization in ParaView.
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

print("=" * 70)
print("MAGNETIC FIELD MAP CALCULATION")
print("=" * 70)
print("\nComplex 8-segment coil geometry\n")


def calculate_field_grid(coil_obj, grid_params):
	"""
	Calculate magnetic field on a 3D grid.

	Args:
		coil_obj: Radia object ID
		grid_params: Dictionary with grid parameters
			- x_range: [min, max, num_points]
			- y_range: [min, max, num_points]
			- z_range: [min, max, num_points]

	Returns:
		Dictionary with grid coordinates and field values
	"""
	print("\n" + "-" * 70)
	print("Calculating field on 3D grid...")
	print("-" * 70)

	# Create grid
	x_range = grid_params['x_range']
	y_range = grid_params['y_range']
	z_range = grid_params['z_range']

	x = np.linspace(x_range[0], x_range[1], x_range[2])
	y = np.linspace(y_range[0], y_range[1], y_range[2])
	z = np.linspace(z_range[0], z_range[1], z_range[2])

	X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

	print(f"  Grid size: {x_range[2]} × {y_range[2]} × {z_range[2]} = {X.size} points")
	print(f"  X range: [{x_range[0]}, {x_range[1]}] mm")
	print(f"  Y range: [{y_range[0]}, {y_range[1]}] mm")
	print(f"  Z range: [{z_range[0]}, {z_range[1]}] mm")

	# Calculate field at each point
	Bx = np.zeros_like(X)
	By = np.zeros_like(X)
	Bz = np.zeros_like(X)

	total_points = X.size
	print(f"\n  Calculating fields at {total_points} points...")

	# Progress indicator
	progress_interval = max(1, total_points // 20)

	for idx in range(total_points):
		i = np.unravel_index(idx, X.shape)
		point = [X[i], Y[i], Z[i]]

		B = rad.Fld(coil_obj, 'b', point)
		Bx[i] = B[0] * 1000  # T to mT
		By[i] = B[1] * 1000
		Bz[i] = B[2] * 1000

		if (idx + 1) % progress_interval == 0:
			percent = (idx + 1) * 100 / total_points
			print(f"    Progress: {percent:.0f}% ({idx + 1}/{total_points})", end='\r')

	print(f"\n  [OK] Field calculation complete")

	# Calculate field magnitude
	B_mag = np.sqrt(Bx**2 + By**2 + Bz**2)

	print(f"\n  Field statistics:")
	print(f"    Bx range: [{np.min(Bx):.3f}, {np.max(Bx):.3f}] mT")
	print(f"    By range: [{np.min(By):.3f}, {np.max(By):.3f}] mT")
	print(f"    Bz range: [{np.min(Bz):.3f}, {np.max(Bz):.3f}] mT")
	print(f"    |B| range: [{np.min(B_mag):.3f}, {np.max(B_mag):.3f}] mT")

	return {
		'X': X, 'Y': Y, 'Z': Z,
		'Bx': Bx, 'By': By, 'Bz': Bz,
		'B_mag': B_mag
	}


def export_field_to_vtk(field_data, filename):
	"""
	Export field data to VTK structured grid format.

	Args:
		field_data: Dictionary with grid coordinates and field values
		filename: Output filename (without .vtk extension)
	"""
	print("\n" + "-" * 70)
	print("Exporting field data to VTK...")
	print("-" * 70)

	X = field_data['X']
	Y = field_data['Y']
	Z = field_data['Z']
	Bx = field_data['Bx']
	By = field_data['By']
	Bz = field_data['Bz']
	B_mag = field_data['B_mag']

	nx, ny, nz = X.shape

	vtk_file = f"{filename}.vtk"

	with open(vtk_file, 'w') as f:
		# Header
		f.write("# vtk DataFile Version 3.0\n")
		f.write("Magnetic field map from Radia\n")
		f.write("ASCII\n")
		f.write("DATASET STRUCTURED_GRID\n")
		f.write(f"DIMENSIONS {nx} {ny} {nz}\n")
		f.write(f"POINTS {X.size} float\n")

		# Write grid points
		for k in range(nz):
			for j in range(ny):
				for i in range(nx):
					f.write(f"{X[i,j,k]:.6f} {Y[i,j,k]:.6f} {Z[i,j,k]:.6f}\n")

		# Write field data
		f.write(f"\nPOINT_DATA {X.size}\n")

		# Vector field (B) only
		f.write("VECTORS B float\n")
		for k in range(nz):
			for j in range(ny):
				for i in range(nx):
					f.write(f"{Bx[i,j,k]:.6f} {By[i,j,k]:.6f} {Bz[i,j,k]:.6f}\n")

	file_size = os.path.getsize(vtk_file)
	file_size_mb = file_size / (1024 * 1024)

	print(f"  [OK] Created: {vtk_file}")
	print(f"       File size: {file_size_mb:.2f} MB")
	print(f"       Grid: {nx}×{ny}×{nz} = {X.size} points")
	print(f"\n  Open in ParaView to visualize:")
	print(f"    - Use 'Glyph' filter to show field vectors")
	print(f"    - Use 'Contour' filter to show field isosurfaces")
	print(f"    - Use 'Slice' filter to show field on planes")


def main():
	"""Main field map calculation script."""

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

	# Get coil bounding box
	info = get_coil_info(coil)
	bounds = info['bbox']

	print(f"\n" + "-" * 70)
	print("Coil bounding box:")
	print("-" * 70)
	print(f"  X: [{bounds['x_min']:.2f}, {bounds['x_max']:.2f}] mm (span: {info['span']['x']:.2f} mm)")
	print(f"  Y: [{bounds['y_min']:.2f}, {bounds['y_max']:.2f}] mm (span: {info['span']['y']:.2f} mm)")
	print(f"  Z: [{bounds['z_min']:.2f}, {bounds['z_max']:.2f}] mm (span: {info['span']['z']:.2f} mm)")

	# Define grid parameters with 100mm margin around coil
	margin = 100.0  # mm

	# Calculate grid range with margin
	x_min = bounds['x_min'] - margin
	x_max = bounds['x_max'] + margin
	y_min = bounds['y_min'] - margin
	y_max = bounds['y_max'] + margin
	z_min = bounds['z_min'] - margin
	z_max = bounds['z_max'] + margin

	# Number of points in each direction
	# Adjust resolution based on span to maintain reasonable aspect ratio
	nx = 31  # X direction
	ny = 51  # Y direction (typically longer)
	nz = 31  # Z direction

	grid_params = {
		'x_range': [x_min, x_max, nx],
		'y_range': [y_min, y_max, ny],
		'z_range': [z_min, z_max, nz],
	}

	print("\n" + "=" * 70)
	print("GRID CONFIGURATION")
	print("=" * 70)
	print(f"\nField evaluation region (100mm margin around coil):")
	print(f"  X: [{x_min:.2f}, {x_max:.2f}] mm")
	print(f"  Y: [{y_min:.2f}, {y_max:.2f}] mm")
	print(f"  Z: [{z_min:.2f}, {z_max:.2f}] mm")

	total_points = nx * ny * nz
	print(f"\nTotal grid points: {total_points:,}")
	print(f"Estimated calculation time: ~{total_points * 0.01:.1f} seconds")
	print("\nNote: For faster calculation, reduce grid resolution.")
	print("      For finer resolution, increase grid points (may take longer).")

	# Calculate field
	field_data = calculate_field_grid(coil, grid_params)

	# Export to VTK
	export_field_to_vtk(field_data, 'field_map')

	# Cleanup
	rad.UtiDelAll()

	print("\n" + "=" * 70)
	print("FIELD MAP CALCULATION COMPLETE")
	print("=" * 70)
	print("\nNext steps:")
	print("  1. Open 'field_map.vtk' in ParaView")
	print("  2. Apply filters to visualize the field:")
	print("     - Glyph: Show field direction with arrows")
	print("     - StreamTracer: Show field lines")
	print("     - Contour: Show constant field magnitude surfaces")
	print("     - Slice: Show field on cutting planes")
	print("=" * 70 + "\n")


if __name__ == '__main__':
	main()
