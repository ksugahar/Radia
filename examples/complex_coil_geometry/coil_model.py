#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Complex coil model definition

This module defines the 8-segment beam steering coil geometry.
The model can be imported by other scripts for visualization or field calculation.
"""

import sys
from pathlib import Path

# Add paths
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root / 'build' / 'lib' / 'Release'))
sys.path.insert(0, str(project_root / 'dist'))
sys.path.insert(0, str(project_root / 'src' / 'python'))

import numpy as np
import radia as rad
from radia_coil_builder import CoilBuilder


def create_beam_steering_coil():
	"""
	Create the 8-segment beam steering magnet coil.

	Returns:
		tuple: (coil_object, coil_parameters)
			- coil_object: Radia object ID
			- coil_parameters: Dictionary with coil specifications
	"""
	# Coil parameters
	I = 1265.0  # Current (A)
	W = 122     # Width (mm)
	H = 122     # Height (mm)

	# Initial orientation and position
	V = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]).T
	L_start = 16.43186645 * 2
	x0 = np.array([(48 + 170 / 1), -L_start / 2, -20 - W / 2])

	# Build coil using modern fluent interface
	coil_segments = (CoilBuilder(current=I)
		.set_start(x0, V)
		.set_cross_section(width=W, height=H)
		.add_straight(length=16.43186645 * 2, tilt=0)
		.add_arc(radius=121, arc_angle=64.59228189, tilt=90)
		.add_straight(length=1018.51313197, tilt=90)
		.add_arc(radius=121, arc_angle=115.40771811, tilt=-90)
		.add_straight(length=453.43186645 * 2, tilt=90)
		.add_arc(radius=121, arc_angle=115.40771811, tilt=-90)
		.add_straight(length=1018.51313197, tilt=90)
		.add_arc(radius=121, arc_angle=64.59228189, tilt=-90)
		.to_radia()
	)

	# Combine all coils
	coils_container = rad.ObjCnt(coil_segments)

	# Set visualization attributes
	rad.ObjDrwAtr(coils_container, [1, 0, 0], 0.1)  # Red, 10% opacity

	# Apply symmetry for faster computation
	rad.TrfZerPara(coils_container, [0, 0, 0], [1, 0, 0])

	# Store parameters for reference
	parameters = {
		'current': I,
		'cross_section': {'width': W, 'height': H},
		'num_segments': len(coil_segments),
		'description': '8-segment beam steering magnet coil'
	}

	return coils_container, parameters


def get_coil_info(coil_obj):
	"""
	Get information about the coil geometry.

	Args:
		coil_obj: Radia object ID

	Returns:
		dict: Coil information including bounding box
	"""
	# Get bounding box
	bbox = rad.ObjGeoLim(coil_obj)

	info = {
		'bbox': {
			'x_min': bbox[0], 'x_max': bbox[1],
			'y_min': bbox[2], 'y_max': bbox[3],
			'z_min': bbox[4], 'z_max': bbox[5]
		},
		'span': {
			'x': bbox[1] - bbox[0],
			'y': bbox[3] - bbox[2],
			'z': bbox[5] - bbox[4]
		}
	}

	return info


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
	print(f"  Cross-section: {params['cross_section']['width']}×{params['cross_section']['height']} mm")
	print(f"  Segments: {params['num_segments']}")

	# Get geometry info
	info = get_coil_info(coil)

	print("\nCoil Geometry:")
	print(f"  Bounding box:")
	print(f"    X: [{info['bbox']['x_min']:.2f}, {info['bbox']['x_max']:.2f}] mm")
	print(f"    Y: [{info['bbox']['y_min']:.2f}, {info['bbox']['y_max']:.2f}] mm")
	print(f"    Z: [{info['bbox']['z_min']:.2f}, {info['bbox']['z_max']:.2f}] mm")
	print(f"  Span:")
	print(f"    X: {info['span']['x']:.2f} mm")
	print(f"    Y: {info['span']['y']:.2f} mm")
	print(f"    Z: {info['span']['z']:.2f} mm")

	# VTK Export - Export geometry with same filename as script
	try:
		from radia_vtk_export import exportGeometryToVTK
		import os

		script_name = os.path.splitext(os.path.basename(__file__))[0]
		vtk_filename = f"{script_name}.vtk"
		vtk_path = os.path.join(os.path.dirname(__file__), vtk_filename)

		exportGeometryToVTK(coil, vtk_path)
		print(f"\n[VTK] Exported: {vtk_filename}")
		print(f"      View with: paraview {vtk_filename}")
	except ImportError:
		print("\n[VTK] Warning: radia_vtk_export not available (VTK export skipped)")
	except Exception as e:
		print(f"\n[VTK] Warning: Export failed: {e}")

	# Cleanup
	rad.UtiDelAll()

	print("\n" + "=" * 70)
	print("[OK] Coil model test complete")
	print("=" * 70)
