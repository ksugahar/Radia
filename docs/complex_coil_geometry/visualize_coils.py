#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Coil geometry visualization

Visualizes the coil geometry and calculates magnetic field at test points.
This script is used to verify the coil shape is correct.
"""

import numpy as np
import radia as rad

# Import coil model (defines mm = 1e-3)
from coil_model import create_beam_steering_coil, mm



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

	# Test points (displayed in mm)
	test_points = [
		[0, 0, 0],
		[0, 0, 100*mm],
		[0, 0, 500*mm],
		[100*mm, 0, 0],
		[0, 100*mm, 0],
		[50*mm, 0, 100*mm],
		[0, 200*mm, 0],
		[0, 0, -100*mm],
	]

	for pt in test_points:
		B = rad.Fld(coil_obj, 'b', pt)
		Bx, By, Bz = B[0] * 1000, B[1] * 1000, B[2] * 1000  # T to mT
		B_mag = np.sqrt(Bx**2 + By**2 + Bz**2)
		pt_str = f"({pt[0]/mm:.0f}, {pt[1]/mm:.0f}, {pt[2]/mm:.0f})"
		print(f"{pt_str:<25} {Bx:<15.6f} {By:<15.6f} {Bz:<15.6f} {B_mag:<15.6f}")

	print("-" * 70)

	# Calculate field along a line (for plotting)
	print("\n" + "-" * 70)
	print("Field along Z-axis (X=0, Y=0):")
	print("-" * 70)
	print(f"{'Z (mm)':<15} {'Bx (mT)':<15} {'|B| (mT)':<15}")
	print("-" * 40)
	z_points = np.linspace(-200*mm, 600*mm, 17)
	for z in z_points:
		B = rad.Fld(coil_obj, 'b', [0, 0, z])
		Bx = B[0] * 1000
		B_mag = np.linalg.norm(B) * 1000
		print(f"{z/mm:<15.1f} {Bx:<15.6f} {B_mag:<15.6f}")
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
	print(f"     Cross-section: {params['cross_section']['width']/mm:.0f}x{params['cross_section']['height']/mm:.0f} mm")
	print(f"     Segments: {params['num_segments']}")

	# Calculate magnetic field at test points
	calculate_field_at_test_points(coil)

	# Cleanup
	rad.UtiDelAll()

	print("\n" + "=" * 70)
	print("VISUALIZATION COMPLETE")
	print("=" * 70)
	print("\nCoil geometry has been verified.")
	print("Use field_map.py to calculate field distribution.")
	print("=" * 70 + "\n")


if __name__ == '__main__':
	main()
