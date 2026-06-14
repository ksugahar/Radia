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

# radia is an installed package (pip install -e .) -- no source-tree sys.path hack needed.
import radia as rad

mm = 1e-3  # 1 mm in meters


def create_meshed_disk(R, H, n_radial, n_angular, n_z=1, x0=0, y0=0, z0=0):
	"""
	Create a meshed circular disk using hexahedral elements.

	The disk is divided into:
	- n_radial rings in radial direction
	- n_angular segments in angular direction
	- n_z layers in vertical (Z) direction

	Args:
		R: Disk radius (m)
		H: Disk height (m)
		n_radial: Number of radial divisions
		n_angular: Number of angular divisions
		n_z: Number of vertical (Z) divisions (default: 1)
		x0, y0, z0: Center position (m)

	Returns:
		Radia container object with hexahedral elements
	"""
	print(f"	Creating meshed disk: {n_radial} radial x {n_angular} angular x {n_z} vertical = {n_radial * n_angular * n_z} elements")

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

					# Create hexahedron (no magnetization for iron base plate)
				# ObjHexahedron auto-generates standard face topology
				hex_elem = rad.ObjHexahedron(points)
				hex_elements.append(hex_elem)

	# Combine into container
	disk = rad.ObjCnt(hex_elements)

	return disk


def create_smco_magnet_array(
	mag_radius=5*mm,	      # Magnet radius (m)
	mag_height=10*mm,	      # Magnet height (m)
	mag_M=[0, 0, 1],	      # Magnetization (T)
	spacing=10*mm,		      # Magnet spacing (m)
	array_radius=60*mm,	      # Array radius (m)
	base_plate_height=5*mm    # Base plate height (m)
):
	"""
	Create a hexagonal array of SmCo magnets on a base plate.

	Args:
		mag_radius: Individual magnet radius (m)
		mag_height: Individual magnet height (m)
		mag_M: Magnetization vector [Mx, My, Mz] (T)
		spacing: Distance between magnet centers (m)
		array_radius: Radius of the entire array (m)
		base_plate_height: Height of the base plate (m)

	Returns:
		tuple: (geometry_object, array_info)
	"""
	print("=" * 70)
	print("Creating SmCo magnet array...")
	print("=" * 70)

	print(f"  Magnet radius: {mag_radius/mm:.2f} mm")
	print(f"  Magnet height: {mag_height/mm:.2f} mm")
	print(f"  Magnetization: {mag_M} T")
	print(f"  Array radius: {array_radius/mm:.2f} mm")
	print(f"  Magnet spacing: {spacing/mm:.2f} mm")

	# Create base plate (meshed iron disk with hexahedral elements)
	print(f"\n  Creating base plate...")
	n_radial = 6	# Number of radial divisions
	n_angular = 24  # Number of angular divisions
	n_z = 2         # Number of vertical (Z) divisions
	base_plate = create_meshed_disk(
		array_radius, base_plate_height, n_radial, n_angular, n_z, 0, 0, 0
	)

	# Apply iron material properties for magnetic yoke behavior
	mat = rad.MatLin(1000)  # mur = 1000 (isotropic)
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

	# Create magnet array (all dimensions in meters, Radia always uses meters)
	geometry, info = create_smco_magnet_array(
		mag_radius=5*mm,	   # 5 mm radius
		mag_height=10*mm,	  # 10 mm height
		mag_M=[0, 0, 1],	  # 1 T vertical magnetization
		spacing=10*mm,		 # 10 mm spacing
		array_radius=60*mm,	# 60 mm array radius
		base_plate_height=20*mm  # 20 mm base plate
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

	# Cleanup
	rad.UtiDelAll()

	print("\n" + "=" * 70)
	print("SIMULATION COMPLETE")
	print("=" * 70)
	print(f"\nSummary:")
	print(f"  Number of magnets: {info['num_magnets']}")
	print(f"  Array radius: {info['array_radius']/mm:.2f} mm")
	print(f"  Magnet radius: {info['mag_radius']/mm:.2f} mm")
	print(f"  Magnet height: {info['mag_height']/mm:.2f} mm")
	print("=" * 70 + "\n")


if __name__ == '__main__':
	main()
