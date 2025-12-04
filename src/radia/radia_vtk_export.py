#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Radia VTK Export Utilities

Functions for exporting Radia geometry to VTK format for visualization
in ParaView and other VTK-compatible tools.
"""

import csv
from itertools import accumulate

def _get_radia_length_unit():
	"""
	Get current Radia length unit by querying rad.FldUnits().

	Returns:
		tuple: (unit_name, scale_to_meters)
		       ('mm', 0.001) if Radia is using millimeters
		       ('m', 1.0) if Radia is using meters

	Raises:
		ValueError: If length unit cannot be determined
	"""
	import radia as rad

	units_str = rad.FldUnits()

	if 'Length:  mm' in units_str:
		return ('mm', 0.001)
	elif 'Length:  m' in units_str:
		return ('m', 1.0)
	else:
		raise ValueError(f"Cannot determine Radia length unit from: {units_str}")

def chunks(lst, n):
	"""
	Yield successive n-sized chunks from a list.

	Args:
		lst: List to be chunked
		n: Chunk size

	Yields:
		Chunks of size n from the input list
	"""
	for i in range(0, len(lst), n):
		yield lst[i:i + n]

def exportGeometryToVTK(obj, fileName='radia_Geometry'):
	"""
	Export Radia object geometry to VTK Legacy format file.

	Writes the geometry of a Radia object to a .vtk file for visualization
	in ParaView. The format is VTK Legacy (ASCII), consisting of polygons only.

	Args:
		obj: Radia object ID (integer)
		fileName: Output filename without extension (default: 'radia_Geometry')

	Output:
		Creates fileName.vtk in the current directory

	Example:
		>>> import radia as rad
		>>> from radia_vtk_export import exportGeometryToVTK
		>>> mag = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1])
		>>> exportGeometryToVTK(mag, 'my_magnet')
	"""
	import radia as rad

	vtkData = rad.ObjDrwVTK(obj, 'Axes->False')

	lengths = vtkData['polygons']['lengths']
	nPoly = len(lengths)
	offsets = list(accumulate(lengths))
	offsets.insert(0, 0) # prepend list with a zero
	points = vtkData['polygons']['vertices']
	nPnts = int(len(points)/3)

	# Get Radia's current length unit and convert to meters for VTK standard
	unit_name, scale_to_meters = _get_radia_length_unit()
	points = [round(num * scale_to_meters, 8) for num in points]

	# define colours array
	colors = vtkData['polygons']['colors']

	# pre-process the output lists to have chunkLength items per line
	chunkLength = 9 # this writes 9 numbers per line
	points = list(chunks(points, chunkLength))
	colors = list(chunks(colors, chunkLength))

	# write the data to file (VTK Legacy format - most compatible)
	with open(fileName + ".vtk", "w", newline="") as f:
		f.write('# vtk DataFile Version 3.0\n')
		f.write('vtk output\nASCII\nDATASET POLYDATA\n')
		f.write('POINTS ' + str(nPnts) + ' float\n')
		writer = csv.writer(f, delimiter=" ")
		writer.writerows(points)
		f.write('\n')

		# POLYGONS in classic format (most compatible with all ParaView versions)
		# Format: nPoly totalSize
		# Each line: nVertices v1 v2 v3 ...
		total_size = sum(lengths) + nPoly  # sum of (nVertices + nVertices) for each polygon
		f.write('POLYGONS ' + str(nPoly) + ' ' + str(total_size) + '\n')
		for i in range(nPoly):
			n_vertices = lengths[i]
			start = offsets[i]
			end = offsets[i+1]
			f.write(str(n_vertices))
			for j in range(start, end):
				f.write(' ' + str(j))
			f.write('\n')

		f.write('\n')
		f.write('CELL_DATA ' + str(nPoly) + '\n')
		f.write('COLOR_SCALARS Radia_colours 3\n')
		writer.writerows(colors)

	print(f"VTK file exported: {fileName}.vtk")
	print(f"  Polygons: {nPoly}")
	print(f"  Points: {nPnts}")


def exportFieldToVTK(points, field_data, fileName='field_distribution', field_name='B_field'):
	"""
	Export magnetic field distribution to VTK Legacy format file.

	Args:
		points: List of observation points [[x, y, z], ...] in Radia's current units
		field_data: List of field vectors [[Bx, By, Bz], ...] in Tesla
		fileName: Output file name (without .vtk extension)
		field_name: Name for the vector field in VTK file

	Note:
		Points are automatically converted to meters (VTK standard) based on
		Radia's current unit setting (queried via rad.FldUnits()).
	"""
	import numpy as np

	points = np.array(points)
	field_data = np.array(field_data)

	if points.shape[0] != field_data.shape[0]:
		raise ValueError(f"Points ({points.shape[0]}) and field data ({field_data.shape[0]}) must have same length")

	if points.shape[1] != 3 or field_data.shape[1] != 3:
		raise ValueError("Points and field data must be Nx3 arrays")

	# Get Radia's current length unit
	unit_name, scale_to_meters = _get_radia_length_unit()

	# Ensure .vtk extension
	if not fileName.endswith('.vtk'):
		fileName += '.vtk'

	with open(fileName, 'w') as f:
		# Header
		f.write('# vtk DataFile Version 3.0\n')
		f.write('Magnetic Field Distribution\n')
		f.write('ASCII\n')
		f.write('DATASET POLYDATA\n')

		# Points (convert to meters for VTK standard)
		f.write(f'POINTS {len(points)} float\n')
		for pt in points:
			# Convert to meters using detected unit scale
			f.write(f'{pt[0]*scale_to_meters} {pt[1]*scale_to_meters} {pt[2]*scale_to_meters}\n')

		# Point data (vector field)
		f.write(f'\nPOINT_DATA {len(points)}\n')
		f.write(f'VECTORS {field_name} float\n')
		for B in field_data:
			f.write(f'{B[0]} {B[1]} {B[2]}\n')

	print(f"VTK field file exported: {fileName}")
	print(f"  Points: {len(points)}")
	print(f"  Field: {field_name}")


if __name__ == '__main__':
	"""
	Demo: Export a simple Radia geometry to VTK format
	"""
	import sys
	import os

	# Add build directory to path
	sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'Release'))

	import radia as rad

	print("=" * 60)
	print("Radia VTK Export Demo")
	print("=" * 60)

	# Create a simple test geometry
	print("\nCreating test geometry...")

	# Rectangular magnet
	mag = rad.ObjRecMag([0, 0, 0], [30, 30, 10], [0, 0, 1])

	# Cylindrical magnet
	cyl = rad.ObjCylMag([50, 0, 0], 15, 20, 16, 'z', [0, 0, 1])

	# Container
	container = rad.ObjCnt([mag, cyl])

	# Export to VTK
	output_file = 'radia_demo_geometry'
	print(f"\nExporting geometry to {output_file}.vtk...")
	exportGeometryToVTK(container, output_file)

	print("\n" + "=" * 60)
	print("Export complete!")
	print("\nTo view in ParaView:")
	print(f"  1. Open ParaView")
	print(f"  2. File → Open → {output_file}.vtk")
	print(f"  3. Click 'Apply' in the Properties panel")
	print("=" * 60)
