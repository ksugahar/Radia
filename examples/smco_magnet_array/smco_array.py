#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Samarium-Cobalt (SmCo) Magnet Array Simulation

Ported from Mathematica notebook: 2023_10_01_サマコバ/magnet.nb
Creates a hexagonal array of cylindrical SmCo magnets.
"""

import sys
import os
from pathlib import Path
import numpy as np

# Add paths
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root / 'build' / 'lib' / 'Release'))
sys.path.insert(0, str(project_root / 'dist'))

import radia as rad

try:
	sys.path.insert(0, str(project_root / 'src' / 'python'))
	from radia_vtk_export import exportGeometryToVTK
	HAS_VTK_EXPORT = True
except ImportError:
	HAS_VTK_EXPORT = False
	print("Warning: VTK export not available")


def create_meshed_disk(R, H, n_radial, n_angular, n_z=1, x0=0, y0=0, z0=0):
	"""
	Create a meshed circular disk using hexahedral elements.

	The disk is divided into:
	- n_radial rings in radial direction
	- n_angular segments in angular direction
	- n_z layers in vertical (Z) direction

	Args:
		R: Disk radius (mm)
		H: Disk height (mm)
		n_radial: Number of radial divisions
		n_angular: Number of angular divisions
		n_z: Number of vertical (Z) divisions (default: 1)
		x0, y0, z0: Center position (mm)

	Returns:
		Radia container object with hexahedral elements
	"""
	print(f"	Creating meshed disk: {n_radial} radial × {n_angular} angular × {n_z} vertical = {n_radial * n_angular * n_z} elements")

	hex_elements = []

	# Radial divisions
	r_vals = np.linspace(0, R, n_radial + 1)

	# Angular divisions
	theta_vals = np.linspace(0, 2 * np.pi, n_angular + 1)

	# Vertical (Z) divisions
	z_vals = np.linspace(z0 - H/2, z0 + H/2, n_z + 1)

	# Create hexahedral elements
	for k in range(n_z):
		z_bottom = z_vals[k]
		z_top = z_vals[k + 1]

		for i in range(n_radial):
			r_inner = r_vals[i]
			r_outer = r_vals[i + 1]

			for j in range(n_angular):
				theta1 = theta_vals[j]
				theta2 = theta_vals[j + 1]

				# Special case: innermost ring (wedge/pentahedron)
				if i == 0:
					# Center point
					x_center = x0
					y_center = y0

					# Outer edge points
					x1_outer = r_outer * np.cos(theta1) + x0
					y1_outer = r_outer * np.sin(theta1) + y0
					x2_outer = r_outer * np.cos(theta2) + x0
					y2_outer = r_outer * np.sin(theta2) + y0

					# Pentahedron (wedge): 6 vertices
					# Bottom: 3 points (center, outer1, outer2)
					# Top: 3 points (center, outer1, outer2)
					points = [
						[x_center, y_center, z_bottom],  # 1: bottom center
						[x1_outer, y1_outer, z_bottom],  # 2: bottom outer1
						[x2_outer, y2_outer, z_bottom],  # 3: bottom outer2
						[x_center, y_center, z_top],     # 4: top center
						[x1_outer, y1_outer, z_top],     # 5: top outer1
						[x2_outer, y2_outer, z_top],     # 6: top outer2
					]

					# Pentahedron faces
					faces = [
						[1, 2, 3],         # Bottom triangle
						[4, 5, 6],         # Top triangle
						[1, 2, 5, 4],      # Side face 1 (quad)
						[2, 3, 6, 5],      # Side face 2 (quad)
						[3, 1, 4, 6],      # Side face 3 (quad)
					]
				else:
					# Regular hexahedron (annular sector)
					# Inner edge points
					x1_inner = r_inner * np.cos(theta1) + x0
					y1_inner = r_inner * np.sin(theta1) + y0
					x2_inner = r_inner * np.cos(theta2) + x0
					y2_inner = r_inner * np.sin(theta2) + y0

					# Outer edge points
					x1_outer = r_outer * np.cos(theta1) + x0
					y1_outer = r_outer * np.sin(theta1) + y0
					x2_outer = r_outer * np.cos(theta2) + x0
					y2_outer = r_outer * np.sin(theta2) + y0

					# 8 vertices of hexahedron (bottom 4, top 4)
					points = [
						[x1_inner, y1_inner, z_bottom],  # 1: bottom inner1
						[x1_outer, y1_outer, z_bottom],  # 2: bottom outer1
						[x2_outer, y2_outer, z_bottom],  # 3: bottom outer2
						[x2_inner, y2_inner, z_bottom],  # 4: bottom inner2
						[x1_inner, y1_inner, z_top],     # 5: top inner1
						[x1_outer, y1_outer, z_top],     # 6: top outer1
						[x2_outer, y2_outer, z_top],     # 7: top outer2
						[x2_inner, y2_inner, z_top],     # 8: top inner2
					]

					# Hexahedron faces (standard ordering)
					faces = [
						[1, 2, 3, 4],      # Bottom
						[5, 6, 7, 8],      # Top
						[1, 2, 6, 5],      # Inner radial side
						[2, 3, 7, 6],      # Outer angular side
						[3, 4, 8, 7],      # Outer radial side
						[4, 1, 5, 8],      # Inner angular side
					]

				# Create polyhedron (no magnetization for iron base plate)
				hex_elem = rad.ObjPolyhdr(points, faces)
				hex_elements.append(hex_elem)

	# Combine into container
	disk = rad.ObjCnt(hex_elements)

	return disk


def create_smco_magnet_array(
	mag_radius=5,	  # Magnet radius (mm)
	mag_height=10,	   # Magnet height (mm)
	mag_M=[0, 0, 1],	   # Magnetization (T)
	spacing=10,		  # Magnet spacing (mm)
	array_radius=60,	 # Array radius (mm)
	base_plate_height=5  # Base plate height (mm)
):
	"""
	Create a hexagonal array of SmCo magnets on a base plate.

	Args:
		mag_radius: Individual magnet radius (mm)
		mag_height: Individual magnet height (mm)
		mag_M: Magnetization vector [Mx, My, Mz] (T)
		spacing: Distance between magnet centers (mm)
		array_radius: Radius of the entire array (mm)
		base_plate_height: Height of the base plate (mm)

	Returns:
		tuple: (geometry_object, array_info)
	"""
	print("=" * 70)
	print("Creating SmCo magnet array...")
	print("=" * 70)

	print(f"  Magnet radius: {mag_radius:.2f} mm")
	print(f"  Magnet height: {mag_height:.2f} mm")
	print(f"  Magnetization: {mag_M} T")
	print(f"  Array radius: {array_radius:.2f} mm")
	print(f"  Magnet spacing: {spacing:.2f} mm")

	# Create base plate (meshed iron disk with hexahedral elements)
	print(f"\n  Creating base plate...")
	n_radial = 6	# Number of radial divisions
	n_angular = 24  # Number of angular divisions
	n_z = 2         # Number of vertical (Z) divisions
	base_plate = create_meshed_disk(
		array_radius, base_plate_height, n_radial, n_angular, n_z, 0, 0, 0
	)
	rad.ObjDrwAtr(base_plate, [0.5, 0.5, 0.5], 0.1)  # Gray color

	# Apply iron material properties for magnetic yoke behavior
	mat = rad.MatLin(1000)  # μr = 1000 (isotropic)
	rad.MatApl(base_plate, mat)

	# Create hexagonal array of magnets
	print(f"  Creating magnet array...")
	magnets = [base_plate]
	magnet_count = 0

	# Hexagonal grid pattern
	for nx in range(-20, 21):
		for ny in range(-20, 21):
			# Hexagonal packing: offset every other row by half spacing
			x = nx * spacing + (ny % 2) * spacing / 2
			y = ny * spacing * np.sqrt(3) / 2
			# Position magnet on top of base plate (no gap)
			# Base plate top: z = base_plate_height/2
			# Magnet center: z = base_plate_height/2 + mag_height/2
			z = base_plate_height / 2 + mag_height / 2

			# Only create magnets within the array radius
			if x**2 + y**2 < array_radius**2:
				# Create cylindrical magnet directly using rad.ObjCylMag
				# rad.ObjCylMag([x,y,z], radius, height, nseg, axis, magnetization)
				magnet = rad.ObjCylMag([x, y, z], mag_radius, mag_height, 16, 'z', mag_M)
				magnets.append(magnet)
				magnet_count += 1

	print(f"  [OK] Created {magnet_count} magnets in hexagonal array")

	# Combine all objects into container
	geometry = rad.ObjCnt(magnets)

	# Set visualization color for magnets (blue)
	rad.ObjDrwAtr(geometry, [0.3, 0.3, 1.0], 0.1)

	array_info = {
		'num_magnets': magnet_count,
		'mag_radius': mag_radius,
		'mag_height': mag_height,
		'magnetization': mag_M,
		'array_radius': array_radius,
		'spacing': spacing
	}

	return geometry, array_info


def main():
	"""Main SmCo magnet array simulation."""
	print("\n" + "=" * 70)
	print("SMCO MAGNET ARRAY SIMULATION")
	print("=" * 70)
	print("\nHexagonal array of cylindrical SmCo magnets\n")

	# Create magnet array (all dimensions in mm for Radia)
	geometry, info = create_smco_magnet_array(
		mag_radius=5,	  # 5 mm radius
		mag_height=10,	   # 10 mm height
		mag_M=[0, 0, 1],	   # 1 T vertical magnetization
		spacing=10,		  # 10 mm spacing
		array_radius=60,	 # 60 mm array radius
		base_plate_height=20  # 20 mm base plate
	)

	# Solve magnetostatics (required for magnetic materials)
	print("\n" + "=" * 70)
	print("Solving magnetostatics...")
	print("=" * 70)
	print(f"  Precision: 0.01")
	print(f"  Max iterations: 1000")

	res = rad.Solve(geometry, 0.01, 1000, 4)
	print(f"  Solver result: {res}")

	# Check for convergence
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

	# Calculate field at test points
	print("\n" + "=" * 70)
	print("Calculating magnetic field...")
	print("=" * 70)

	test_points = [
		[0, 0, 0.02],	# 20 mm above center
		[0, 0, 0.05],	# 50 mm above center
		[0.03, 0, 0.02], # 30 mm off-axis, 20 mm above
	]

	print(f"{'Position (m)':<25} {'Bx (mT)':<15} {'By (mT)':<15} {'Bz (mT)':<15} {'|B| (mT)':<15}")
	print("-" * 85)

	for pos in test_points:
		B = rad.Fld(geometry, 'b', pos)
		Bx, By, Bz = B[0] * 1000, B[1] * 1000, B[2] * 1000
		B_mag = np.sqrt(Bx**2 + By**2 + Bz**2)
		pos_str = f"({pos[0]:.3f}, {pos[1]:.3f}, {pos[2]:.3f})"
		print(f"{pos_str:<25} {Bx:<15.6f} {By:<15.6f} {Bz:<15.6f} {B_mag:<15.6f}")

	# Export to VTK
	if HAS_VTK_EXPORT:
		print("\n" + "=" * 70)
		print("Exporting to VTK...")
		print("=" * 70)
		output_path = os.path.join(os.path.dirname(__file__), 'smco_array')
		exportGeometryToVTK(geometry, output_path)
		print(f"  [OK] Created: smco_array.vtk")

	# Export field distribution to VTK
	print("\n" + "=" * 70)
	print("Calculating field distribution...")
	print("=" * 70)

	# Get bounding box of entire geometry
	bbox = rad.ObjGeoLim(geometry)
	print(f"  Geometry bounding box:")
	print(f"    X: [{bbox[0]:.2f}, {bbox[1]:.2f}] mm")
	print(f"    Y: [{bbox[2]:.2f}, {bbox[3]:.2f}] mm")
	print(f"    Z: [{bbox[4]:.2f}, {bbox[5]:.2f}] mm")

	# Expand bbox by 20mm in all directions
	margin = 20.0
	x_min, x_max = bbox[0] - margin, bbox[1] + margin
	y_min, y_max = bbox[2] - margin, bbox[3] + margin
	z_min, z_max = bbox[4] - margin, bbox[5] + margin

	print(f"\n  Field calculation range (bbox + 20mm):")
	print(f"    X: [{x_min:.2f}, {x_max:.2f}] mm")
	print(f"    Y: [{y_min:.2f}, {y_max:.2f}] mm")
	print(f"    Z: [{z_min:.2f}, {z_max:.2f}] mm")

	# Create grid for field calculation
	nx, ny, nz = 21, 21, 21  # Grid resolution
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
				B = rad.Fld(geometry, 'b', [x, y, z])
				field_data.append([x, y, z, B[0], B[1], B[2]])
				calculated += 1

				if calculated % 1000 == 0:
					print(f"    Progress: {calculated}/{total_points} points", end='\r')

	print(f"    Progress: {total_points}/{total_points} points")
	print(f"  [OK] Field calculation complete")

	# Export field to VTK
	field_vtk_path = os.path.join(os.path.dirname(__file__), 'smco_field_distribution.vtk')
	with open(field_vtk_path, 'w') as f:
		f.write("# vtk DataFile Version 3.0\n")
		f.write("SmCo magnet array field distribution\n")
		f.write("ASCII\n")
		f.write("DATASET STRUCTURED_POINTS\n")
		f.write(f"DIMENSIONS {nx} {ny} {nz}\n")
		f.write(f"ORIGIN {x_min} {y_min} {z_min}\n")
		f.write(f"SPACING {(x_max-x_min)/(nx-1)} {(y_max-y_min)/(ny-1)} {(z_max-z_min)/(nz-1)}\n")
		f.write(f"POINT_DATA {nx*ny*nz}\n")
		f.write("VECTORS B_field float\n")

		for data in field_data:
			f.write(f"{data[3]} {data[4]} {data[5]}\n")

	print(f"\n  [OK] Created: smco_field_distribution.vtk")
	print(f"       Open in ParaView and use 'Glyph' filter to visualize vectors")

	# Cleanup
	rad.UtiDelAll()

	print("\n" + "=" * 70)
	print("SIMULATION COMPLETE")
	print("=" * 70)
	print(f"\nSummary:")
	print(f"  Number of magnets: {info['num_magnets']}")
	print(f"  Array radius: {info['array_radius']:.2f} mm")
	print(f"  Magnet radius: {info['mag_radius']:.2f} mm")
	print(f"  Magnet height: {info['mag_height']:.2f} mm")
	print(f"\nOutput files:")
	if HAS_VTK_EXPORT:
		print(f"  - smco_array.vtk (geometry)")
	print(f"  - smco_field_distribution.vtk (magnetic field vectors)")
	print("=" * 70 + "\n")


if __name__ == '__main__':
	main()
