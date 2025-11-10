#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Racetrack coil model for electromagnet

Creates a racetrack coil matching the York.nas yoke geometry.
Based on specifications from the original Radia.py script.
"""

import sys
from pathlib import Path

# Add paths
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root / 'build' / 'lib' / 'Release'))
sys.path.insert(0, str(project_root / 'dist'))

import numpy as np
import radia as rad


def create_racetrack_coil(current=-2000):
	"""
	Create racetrack coil for electromagnet.

	Original Mathematica specification:
	coil = radObjRaceTrk[{0, 131.25, 0}, {5, 40}, {50, 62.5}, 105, 3, current/105/35]

	Args:
		current: Total current (A), default -2000 A

	Returns:
		tuple: (coil_object, coil_parameters)
	"""
	print("=" * 70)
	print("Creating racetrack coil...")
	print("=" * 70)

	# Radia racetrack coil parameters (from Mathematica specification)
	# radObjRaceTrk[center, {lx_min, lx_max}, {ly_min, ly_max}, height, nseg, J]

	center = [0, 131.25, 0]  # Center position (mm)
	lx_min = 5               # Inner semi-width X (mm)
	lx_max = 40              # Outer semi-width X (mm)
	ly_min = 50              # Inner semi-length Y (mm)
	ly_max = 62.5            # Outer semi-length Y (mm)
	height = 105             # Height Z (mm)
	nseg = 3                 # Number of segments for arc approximation

	# Current density J = total_current / (turns * cross_section_area)
	# From original: current / 105 / 35
	# This means: 105 turns and effective area of 35 mm²
	turns = 105
	area = 35  # mm²
	J = current / turns / area  # A/mm²

	print(f"  Center: {center} mm")
	print(f"  X dimensions: inner={lx_min} mm, outer={lx_max} mm")
	print(f"  Y dimensions: inner={ly_min} mm, outer={ly_max} mm")
	print(f"  Height: {height} mm")
	print(f"  Current: {current} A")
	print(f"  Turns: {turns}")
	print(f"  Current density: {J:.6f} A/mm^2")
	print(f"  Arc segments: {nseg}")

	# Create racetrack coil using Radia
	coil = rad.ObjRaceTrk(
		center,
		[lx_min, lx_max],
		[ly_min, ly_max],
		height,
		nseg,
		J
	)

	# Set visualization color (red)
	rad.ObjDrwAtr(coil, [1, 0, 0], 0.1)

	print(f"  [OK] Racetrack coil created")

	# Store parameters
	parameters = {
		'type': 'racetrack',
		'current': current,
		'turns': turns,
		'center': center,
		'dimensions': {
			'lx_min': lx_min,
			'lx_max': lx_max,
			'ly_min': ly_min,
			'ly_max': ly_max,
			'height': height
		},
		'nseg': nseg,
		'current_density': J
	}

	return coil, parameters


def get_coil_info(coil_obj):
	"""
	Get information about the coil geometry.

	Args:
		coil_obj: Radia object ID

	Returns:
		dict: Coil information including bounding box
	"""
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
	"""Test racetrack coil creation"""
	print("=" * 70)
	print("RACETRACK COIL TEST")
	print("=" * 70 + "\n")

	# Create coil
	coil, params = create_racetrack_coil()

	print("\nCoil Parameters:")
	print(f"  Type: {params['type']}")
	print(f"  Current: {params['current']} A")
	print(f"  Turns: {params['turns']}")
	print(f"  Current density: {params['current_density']:.6f} A/mm^2")

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

	# Test field calculation
	print("\nTest field calculation (coil only):")
	positions = [[0, 0, 0], [0, 131.25, 0], [0, 0, 50]]
	for pos in positions:
		B = rad.Fld(coil, 'b', pos)
		Bx, By, Bz = B[0] * 1000, B[1] * 1000, B[2] * 1000
		print(f"  {pos}: B = [{Bx:.6f}, {By:.6f}, {Bz:.6f}] mT")

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
	print("[OK] Racetrack coil test complete")
	print("=" * 70)
