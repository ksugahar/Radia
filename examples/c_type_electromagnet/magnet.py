#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Beam steering electromagnet simulation

Complete simulation with coil + magnetic yoke from Nastran mesh.
"""

import os
import sys
from pathlib import Path

# Add radia module path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root / 'build' / 'lib' / 'Release'))
sys.path.insert(0, str(project_root / 'dist'))
sys.path.insert(0, str(project_root / 'src' / 'python'))

# Add complex_coil_geometry for coil model
sys.path.insert(0, str(project_root / 'examples' / 'complex_coil_geometry'))

import numpy as np
import radia as rad

# Import racetrack coil model
from racetrack_coil_model import create_racetrack_coil, get_coil_info

# Import yoke model
from yoke_model import create_yoke_from_nastran

try:
	from radia_vtk_export import exportGeometryToVTK
	HAS_VTK_EXPORT = True
except ImportError:
	HAS_VTK_EXPORT = False
	print("Warning: VTK export not available")


def main():
	"""Main electromagnet simulation."""
	print("\n" + "=" * 70)
	print("RACETRACK ELECTROMAGNET SIMULATION")
	print("=" * 70)
	print("\nRacetrack Coil + Magnetic Yoke (Nastran mesh)\n")

	# Create racetrack coil
	print("-" * 70)
	print("Creating racetrack coil...")
	print("-" * 70)
	coil, coil_params = create_racetrack_coil(current=-2000)

	print(f"[OK] Coil created")
	print(f"     Type: {coil_params['type']}")
	print(f"     Current: {coil_params['current']} A")
	print(f"     Turns: {coil_params['turns']}")
	print(f"     Current density: {coil_params['current_density']:.6f} A/mm^2")

	# Get coil info
	coil_info = get_coil_info(coil)
	print(f"\n     Bounding box:")
	print(f"       X: [{coil_info['bbox']['x_min']:.2f}, {coil_info['bbox']['x_max']:.2f}] mm")
	print(f"       Y: [{coil_info['bbox']['y_min']:.2f}, {coil_info['bbox']['y_max']:.2f}] mm")
	print(f"       Z: [{coil_info['bbox']['z_min']:.2f}, {coil_info['bbox']['z_max']:.2f}] mm")

	# Create yoke from Nastran mesh
	nas_file = os.path.join(os.path.dirname(__file__), 'York.bdf')

	if os.path.exists(nas_file):
		yoke, yoke_info = create_yoke_from_nastran(
			nas_file,
			material_properties={'permeability': 1000, 'type': 'linear'}
		)
	else:
		print(f"\n[WARNING] Nastran file not found: {nas_file}")
		print(f"[INFO] Continuing with coil only")
		yoke = None
		yoke_info = None

	# Combine geometry
	if yoke:
		g = rad.ObjCnt([coil, yoke])
		print(f"\n[INFO] Geometry: Coils + Yoke")
		print(f"       Yoke elements: {yoke_info['num_hex_elements']} hex + {yoke_info['num_penta_elements']} penta = {yoke_info['num_total_elements']} total")
	else:
		g = coil
		print(f"\n[INFO] Geometry: Coils only (no yoke)")

	# Solve magnetostatics (only needed if magnetic materials present)
	if yoke:
		print("\n" + "=" * 70)
		print("Solving magnetostatics...")
		print("=" * 70)
		print(f"  Precision: 0.01")
		print(f"  Max iterations: 1000")

		res = rad.Solve(g, 0.01, 1000, 4)
		print(f"  Solver result: {res}")
		print(f"  Result type: {type(res)}")

		# Check for NaN in solver result
		if isinstance(res, (list, tuple)):
			has_nan = any(str(x) == 'nan' for x in res)
			if has_nan:
				print(f"  [ERROR] Solver returned NaN - geometry or material issue")
			else:
				print(f"  [OK] Solver completed (iterations: {res[-1] if len(res) > 0 else 'unknown'})")
		else:
			if res > 0:
				print(f"  [OK] Solver converged")
			else:
				print(f"  [WARNING] Solver may not have converged properly")
	else:
		print("\n[INFO] No magnetic materials - solver not needed")

	# Calculate field
	print("\n" + "=" * 70)
	print("Calculating magnetic field...")
	print("=" * 70)

	positions = [[0, 0, 0], [0, 0, 100], [0, 0, 500]]
	print(f"{'Position (mm)':<20} {'Bx (mT)':<15} {'By (mT)':<15} {'Bz (mT)':<15} {'|B| (mT)':<15}")
	print("-" * 70)

	for pos in positions:
		B = rad.Fld(g, "b", pos)
		Bx, By, Bz = B[0] * 1000, B[1] * 1000, B[2] * 1000
		B_mag = np.sqrt(Bx**2 + By**2 + Bz**2)
		pos_str = f"({pos[0]:.0f}, {pos[1]:.0f}, {pos[2]:.0f})"
		print(f"{pos_str:<20} {Bx:<15.6f} {By:<15.6f} {Bz:<15.6f} {B_mag:<15.6f}")

	# Export geometry to VTK
	if HAS_VTK_EXPORT:
		print("\n" + "=" * 70)
		print("Exporting geometry to VTK...")
		print("=" * 70)
		output_path = os.path.join(os.path.dirname(__file__), 'electromagnet')
		exportGeometryToVTK(g, output_path)
		print(f"  [OK] Created: electromagnet.vtk")

	# Export field distribution to VTK
	if yoke:
		print("\n" + "=" * 70)
		print("Calculating field distribution...")
		print("=" * 70)

		# Get bounding box of entire geometry
		bbox = rad.ObjGeoLim(g)
		print(f"  Geometry bounding box:")
		print(f"    X: [{bbox[0]:.2f}, {bbox[1]:.2f}] mm")
		print(f"    Y: [{bbox[2]:.2f}, {bbox[3]:.2f}] mm")
		print(f"    Z: [{bbox[4]:.2f}, {bbox[5]:.2f}] mm")

		# Expand bbox by 50mm in all directions
		margin = 50.0
		x_min, x_max = bbox[0] - margin, bbox[1] + margin
		y_min, y_max = bbox[2] - margin, bbox[3] + margin
		z_min, z_max = bbox[4] - margin, bbox[5] + margin

		print(f"\n  Field calculation range (bbox + 50mm):")
		print(f"    X: [{x_min:.2f}, {x_max:.2f}] mm")
		print(f"    Y: [{y_min:.2f}, {y_max:.2f}] mm")
		print(f"    Z: [{z_min:.2f}, {z_max:.2f}] mm")

		# Create grid for field calculation
		nx, ny, nz = 21, 31, 21  # Grid resolution
		x_vals = np.linspace(x_min, x_max, nx)
		y_vals = np.linspace(y_min, y_max, ny)
		z_vals = np.linspace(z_min, z_max, nz)

		print(f"\n  Grid resolution: {nx} × {ny} × {nz} = {nx*ny*nz} points")
		print(f"  Calculating magnetic field...")

		# Calculate field at grid points
		# VTK STRUCTURED_POINTS ordering: Z varies slowest, then Y, X varies fastest
		field_data = []
		total_points = nx * ny * nz
		calculated = 0

		for iz, z in enumerate(z_vals):
			for iy, y in enumerate(y_vals):
				for ix, x in enumerate(x_vals):
					B = rad.Fld(g, 'b', [x, y, z])
					field_data.append([x, y, z, B[0], B[1], B[2]])
					calculated += 1

					if calculated % 1000 == 0:
						print(f"    Progress: {calculated}/{total_points} points", end='\r')

		print(f"    Progress: {total_points}/{total_points} points")
		print(f"  [OK] Field calculation complete")

		# Export field to VTK
		field_vtk_path = os.path.join(os.path.dirname(__file__), 'field_distribution.vtk')
		with open(field_vtk_path, 'w') as f:
			f.write("# vtk DataFile Version 3.0\n")
			f.write("Magnetic field distribution\n")
			f.write("ASCII\n")
			f.write("DATASET STRUCTURED_POINTS\n")
			f.write(f"DIMENSIONS {nx} {ny} {nz}\n")
			f.write(f"ORIGIN {x_min} {y_min} {z_min}\n")
			f.write(f"SPACING {(x_max-x_min)/(nx-1)} {(y_max-y_min)/(ny-1)} {(z_max-z_min)/(nz-1)}\n")
			f.write(f"POINT_DATA {nx*ny*nz}\n")
			f.write("VECTORS B_field float\n")

			for data in field_data:
				f.write(f"{data[3]} {data[4]} {data[5]}\n")

		print(f"\n  [OK] Created: field_distribution.vtk")
		print(f"       Open in ParaView and use 'Glyph' filter to visualize vectors")

	# Try OpenGL viewer (may fail on some systems)
	print("\n" + "=" * 70)
	print("Opening OpenGL viewer...")
	print("=" * 70)
	try:
		rad.ObjDrwOpenGL(g)
	except Exception as e:
		print(f"  [INFO] OpenGL viewer not available: {e}")
		print(f"  Use ParaView to visualize electromagnet.vtk instead")

	print("\n" + "=" * 70)
	print("SIMULATION COMPLETE")
	print("=" * 70)
	print(f"\nOutput files:")
	print(f"  - electromagnet.vtk (geometry)")
	if yoke:
		print(f"  - field_distribution.vtk (magnetic field vectors)")
	print("=" * 70 + "\n")


if __name__ == '__main__':
	main()
