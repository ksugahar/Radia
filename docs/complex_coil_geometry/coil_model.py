#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Complex coil model definition

This module defines the 8-segment beam steering coil geometry.
The model can be imported by other scripts for visualization or field calculation.

Units: All geometry in meters (Radia always uses meters), current in Amperes.
"""

import numpy as np
import radia as rad
from radia.coil_builder import CoilBuilder

# Set unit system to meters

# Unit conversion factor: original design values are in mm
mm = 1e-3  # 1 mm in meters


def create_beam_steering_coil():
	"""
	Create the 8-segment beam steering magnet coil.

	Returns:
		tuple: (coil_object, coil_parameters)
			- coil_object: Radia object ID
			- coil_parameters: Dictionary with coil specifications
	"""
	# Coil parameters
	I = 1265.0       # Current (A)
	W = 122 * mm     # Width: 122 mm -> 0.122 m
	H = 122 * mm     # Height: 122 mm -> 0.122 m

	# Initial orientation and position
	V = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]).T
	L_start = 16.43186645 * 2 * mm
	x0 = np.array([(48 + 170) * mm, -L_start / 2, -20 * mm - W / 2])

	# Build coil using modern fluent interface
	coil_segments = (CoilBuilder(current=I)
		.set_start(x0, V)
		.set_cross_section(width=W, height=H)
		.add_straight(length=16.43186645 * 2 * mm, tilt=0)
		.add_arc(radius=121 * mm, arc_angle=64.59228189, tilt=90)
		.add_straight(length=1018.51313197 * mm, tilt=90)
		.add_arc(radius=121 * mm, arc_angle=115.40771811, tilt=-90)
		.add_straight(length=453.43186645 * 2 * mm, tilt=90)
		.add_arc(radius=121 * mm, arc_angle=115.40771811, tilt=-90)
		.add_straight(length=1018.51313197 * mm, tilt=90)
		.add_arc(radius=121 * mm, arc_angle=64.59228189, tilt=-90)
		.to_radia()
	)

	# Combine all coils
	coils_container = rad.ObjCnt(coil_segments)

	# Note: TrfZerPara was removed from API (2026-01-31).
	# Symmetry hint is optional; field calculation is still correct without it.

	# Store parameters for reference
	parameters = {
		'current': I,
		'cross_section': {'width': W, 'height': H},
		'num_segments': len(coil_segments),
		'description': '8-segment beam steering magnet coil'
	}

	return coils_container, parameters


if __name__ == '__main__':
	"""
	Test: Create coil and display basic information
	"""
	print("=" * 70)
	print("COIL MODEL TEST")
	print("=" * 70)

	# Create coil
	coil, params = create_beam_steering_coil()

	print("\nCoil Parameters:")
	print(f"  Description: {params['description']}")
	print(f"  Current: {params['current']} A")
	print(f"  Cross-section: {params['cross_section']['width']/mm:.0f}x{params['cross_section']['height']/mm:.0f} mm")
	print(f"  Segments: {params['num_segments']}")

	# Cleanup
	rad.UtiDelAll()

	print("\n" + "=" * 70)
	print("[OK] Coil model test complete")
	print("=" * 70)
