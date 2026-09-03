#!/usr/bin/env python
"""Calculate a magnetic field map and export checked GMSH post data."""

import os
from pathlib import Path

import numpy as np

# Import coil model (defines mm = 1e-3)
from coil_model import create_beam_steering_coil, mm

import radia as rad

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

	print(f"  Grid size: {x_range[2]} x {y_range[2]} x {z_range[2]} = {X.size} points")
	print(f"  X range: [{x_range[0]/mm:.1f}, {x_range[1]/mm:.1f}] mm")
	print(f"  Y range: [{y_range[0]/mm:.1f}, {y_range[1]/mm:.1f}] mm")
	print(f"  Z range: [{z_range[0]/mm:.1f}, {z_range[1]/mm:.1f}] mm")

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

	print("\n  [OK] Field calculation complete")

	# Calculate field magnitude
	B_mag = np.sqrt(Bx**2 + By**2 + Bz**2)

	print("\n  Field statistics:")
	print(f"    Bx range: [{np.min(Bx):.3f}, {np.max(Bx):.3f}] mT")
	print(f"    By range: [{np.min(By):.3f}, {np.max(By):.3f}] mT")
	print(f"    Bz range: [{np.min(Bz):.3f}, {np.max(Bz):.3f}] mT")
	print(f"    |B| range: [{np.min(B_mag):.3f}, {np.max(B_mag):.3f}] mT")

	return {
		'X': X, 'Y': Y, 'Z': Z,
		'Bx': Bx, 'By': By, 'Bz': Bz,
		'B_mag': B_mag
	}


def export_field_to_gmsh(coil_obj, grid_params, filename):
	"""Export the field on a Netgen air mesh as GMSH MSH v4.1 NodeData."""
	from netgen.csg import CSGeometry, OrthoBrick, Pnt
	from ngsolve import Mesh

	from radia.gmsh_post_export import GmshPostExport

	print("\n" + "-" * 70)
	print("Exporting checked GMSH field data...")
	print("-" * 70)

	x_min, x_max, nx = grid_params['x_range']
	y_min, y_max, ny = grid_params['y_range']
	z_min, z_max, nz = grid_params['z_range']
	spacing = max(
		(x_max - x_min) / max(nx - 1, 1),
		(y_max - y_min) / max(ny - 1, 1),
		(z_max - z_min) / max(nz - 1, 1),
	)

	geometry = CSGeometry()
	geometry.Add(
		OrthoBrick(Pnt(x_min, y_min, z_min), Pnt(x_max, y_max, z_max)).mat('air')
	)
	mesh = Mesh(geometry.GenerateMesh(maxh=spacing))
	points = [[float(value) for value in vertex.point] for vertex in mesh.vertices]
	B_mT = np.asarray(rad.Fld(coil_obj, 'b', points), dtype=float) * 1000.0

	output = f"{filename}.msh"
	post = GmshPostExport(mesh)
	post.add_vector_field('B_mT', B_mT)
	post.write(output)

	file_size_mb = os.path.getsize(output) / (1024 * 1024)
	print(f"  [OK] Created: {output}")
	print(f"       File size: {file_size_mb:.2f} MB")
	print(f"       Mesh: {mesh.nv} vertices, {mesh.ne} volume elements")
	print(f"\n  Open {Path(output).with_suffix('.geo')} in GMSH for review.")


def main():
	"""Main field map calculation script."""

	# Create coil from model
	print("-" * 70)
	print("Loading coil model...")
	print("-" * 70)
	coil, params = create_beam_steering_coil()

	print("[OK] Coil model loaded")
	print(f"     Description: {params['description']}")
	print(f"     Current: {params['current']} A")
	print(f"     Cross-section: {params['cross_section']['width']/mm:.0f}x{params['cross_section']['height']/mm:.0f} mm")
	print(f"     Segments: {params['num_segments']}")

	# Fixed evaluation grid sized for the 8-segment beam steering coil.
	# Covers the full coil extent (~1.1 m x 1.5 m x 0.5 m) with ~100 mm margin.
	x_min, x_max = -400 * mm,  800 * mm
	y_min, y_max = -800 * mm,  800 * mm
	z_min, z_max = -300 * mm,  300 * mm

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
	print("\nField evaluation region (100 mm margin around coil):")
	print(f"  X: [{x_min/mm:.2f}, {x_max/mm:.2f}] mm")
	print(f"  Y: [{y_min/mm:.2f}, {y_max/mm:.2f}] mm")
	print(f"  Z: [{z_min/mm:.2f}, {z_max/mm:.2f}] mm")

	total_points = nx * ny * nz
	print(f"\nTotal grid points: {total_points:,}")
	print(f"Estimated calculation time: ~{total_points * 0.01:.1f} seconds")
	print("\nNote: For faster calculation, reduce grid resolution.")
	print("      For finer resolution, increase grid points (may take longer).")

	# Calculate field
	calculate_field_grid(coil, grid_params)

	# Export a checked GMSH post-processing artifact.
	from ngsolve import TaskManager

	with TaskManager():
		export_field_to_gmsh(coil, grid_params, 'field_map')

	# Cleanup
	rad.UtiDelAll()

	print("\n" + "=" * 70)
	print("FIELD MAP CALCULATION COMPLETE")
	print("=" * 70)
	print("\nNext steps:")
	print("  1. Open 'field_map.geo' in GMSH")
	print("  2. Select the B_mT vector view and inspect sections or clipping planes")
	print("=" * 70 + "\n")


if __name__ == '__main__':
	main()
