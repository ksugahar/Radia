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


class RadiaVTKOutput:
	"""
	VTK output class for Radia magnetic field analysis.

	Follows NGSolve's VTKOutput pattern with:
	- Constructor sets up the output configuration
	- Do() method writes the actual file
	- Support for VTS (structured grid) format for 3D field grids

	Example:
		import radia as rad
		from radia import RadiaVTKOutput

		rad.FldUnits('m')

		# Create magnet
		vertices = [[-0.05,-0.05,-0.05], [0.05,-0.05,-0.05], ...]
		magnet = rad.ObjHexahedron(vertices, [0, 0, 954930])

		# Define observation grid
		grid_params = {
			'x_range': [-0.2, 0.2],
			'y_range': [-0.2, 0.2],
			'z_range': [0.1, 0.3],
			'nx': 41, 'ny': 41, 'nz': 21
		}

		# Export to VTS
		vtk = RadiaVTKOutput(
			obj=magnet,
			coefs=['B', 'B_magnitude'],
			filename='B_field',
			grid_params=grid_params
		)
		vtk.Do()  # Creates B_field.vts
	"""

	# Available coefficient names for field output
	AVAILABLE_COEFS = [
		'B',            # [T] Magnetic flux density (vector)
		'H',            # [A/m] Magnetic field intensity (vector)
		'B_magnitude',  # [T] |B| magnitude (scalar)
		'H_magnitude',  # [A/m] |H| magnitude (scalar)
		'Bx', 'By', 'Bz',  # [T] B components (scalar)
		'Hx', 'Hy', 'Hz',  # [A/m] H components (scalar)
	]

	def __init__(self, obj=None, coefs=None, filename='radia_output',
	             floatsize='double', legacy=False, grid_params=None):
		"""
		Initialize Radia VTK output.

		Parameters:
			obj: Radia object handle (integer)
			coefs: List of coefficient names to export
				   Default: ['B', 'B_magnitude']
				   Available: 'B', 'H', 'B_magnitude', 'H_magnitude',
				              'Bx', 'By', 'Bz', 'Hx', 'Hy', 'Hz'
			filename: Base filename without extension (default: 'radia_output')
			floatsize: 'single' or 'double' precision (default: 'double')
			legacy: Use legacy VTK format (default: False, uses VTS)
			grid_params: Dict for 3D field grid (required for field output):
				'x_range': [x_min, x_max] in current Radia units
				'y_range': [y_min, y_max] in current Radia units
				'z_range': [z_min, z_max] in current Radia units
				'nx', 'ny', 'nz': Grid points per direction (default: 21 each)
		"""
		self.obj = obj
		self.filename = filename
		self.floatsize = floatsize
		self.legacy = legacy
		self.grid_params = grid_params

		# Default coefficients
		if coefs is None:
			coefs = ['B', 'B_magnitude']
		self.coefs = coefs

		# Validate coefficients
		for coef in self.coefs:
			if coef not in self.AVAILABLE_COEFS:
				raise ValueError(f"Unknown coefficient: {coef}. Available: {self.AVAILABLE_COEFS}")

	def Do(self, time=None):
		"""
		Write VTK output file.

		Parameters:
			time: Optional time stamp (for future time series support)

		Returns:
			str: The filename that was written
		"""
		if self.grid_params is None:
			raise ValueError("grid_params is required for field output")

		if self.obj is None:
			raise ValueError("obj (Radia object handle) is required")

		return self._write_field_grid_vts(self.filename)

	def _get_float_format(self):
		"""Get format string for float output."""
		if self.floatsize == 'single':
			return '{:.6e}'
		else:
			return '{:.15e}'

	def _compute_field_grid(self, x, y, z):
		"""
		Compute B and H fields on 3D grid using FldBatch (OpenMP parallelized).

		Parameters:
			x: 1D array of x coordinates
			y: 1D array of y coordinates
			z: 1D array of z coordinates

		Returns:
			tuple: (B_data, H_data) as (nx*ny*nz, 3) arrays
		"""
		import radia as rad
		import numpy as np

		# Build point list in VTS order (i varies fastest, then j, then k)
		points = []
		for k in range(len(z)):
			for j in range(len(y)):
				for i in range(len(x)):
					points.append([x[i], y[j], z[k]])

		n_points = len(points)
		print(f"Computing field at {n_points} points using FldBatch (OpenMP)...")

		# FldBatch is OpenMP-parallelized in C++
		# Returns {'B': [...], 'H': [...]} or similar structure
		result = rad.FldBatch(self.obj, points, 0)

		# Extract B and H data
		# FldBatch returns flat list: [Bx1,By1,Bz1, Bx2,By2,Bz2, ...]
		B_flat = np.array(result['B'])
		H_flat = np.array(result['H'])

		# Reshape to (n_points, 3)
		B_data = B_flat.reshape(-1, 3)
		H_data = H_flat.reshape(-1, 3)

		return B_data, H_data

	def _write_field_grid_vts(self, filename):
		"""
		Write 3D field grid in VTS (Structured Grid) format.

		Parameters:
			filename: Base filename without extension

		Returns:
			str: The filename that was written
		"""
		import numpy as np

		params = self.grid_params
		x_min, x_max = params['x_range']
		y_min, y_max = params['y_range']
		z_min, z_max = params['z_range']
		nx = params.get('nx', 21)
		ny = params.get('ny', 21)
		nz = params.get('nz', 21)

		# Create coordinate arrays in Radia's current units
		x = np.linspace(x_min, x_max, nx)
		y = np.linspace(y_min, y_max, ny)
		z = np.linspace(z_min, z_max, nz)

		# Compute fields using FldBatch
		B_data, H_data = self._compute_field_grid(x, y, z)

		# Get unit conversion factor for coordinates
		unit_name, scale_to_meters = _get_radia_length_unit()

		fmt = self._get_float_format()
		output_file = f"{filename}.vts"

		with open(output_file, 'w') as f:
			f.write('<?xml version="1.0"?>\n')
			f.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">\n')
			f.write(f'  <StructuredGrid WholeExtent="0 {nx-1} 0 {ny-1} 0 {nz-1}">\n')
			f.write(f'    <Piece Extent="0 {nx-1} 0 {ny-1} 0 {nz-1}">\n')

			# Points (convert to meters for VTK standard)
			f.write('      <Points>\n')
			f.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
			for k in range(nz):
				for j in range(ny):
					for i in range(nx):
						# Convert to meters
						xm = x[i] * scale_to_meters
						ym = y[j] * scale_to_meters
						zm = z[k] * scale_to_meters
						f.write(f"          {fmt.format(xm)} {fmt.format(ym)} {fmt.format(zm)}\n")
			f.write('        </DataArray>\n')
			f.write('      </Points>\n')

			# Point data
			f.write('      <PointData>\n')

			# Write requested coefficients
			for coef in self.coefs:
				if coef == 'B':
					# B field vector [T]
					f.write('        <DataArray type="Float64" Name="B" NumberOfComponents="3" format="ascii">\n')
					for idx in range(len(B_data)):
						Bx, By, Bz = B_data[idx]
						f.write(f"          {fmt.format(Bx)} {fmt.format(By)} {fmt.format(Bz)}\n")
					f.write('        </DataArray>\n')

				elif coef == 'H':
					# H field vector [A/m]
					f.write('        <DataArray type="Float64" Name="H" NumberOfComponents="3" format="ascii">\n')
					for idx in range(len(H_data)):
						Hx, Hy, Hz = H_data[idx]
						f.write(f"          {fmt.format(Hx)} {fmt.format(Hy)} {fmt.format(Hz)}\n")
					f.write('        </DataArray>\n')

				elif coef == 'B_magnitude':
					# |B| [T]
					f.write('        <DataArray type="Float64" Name="B_magnitude" format="ascii">\n')
					for idx in range(len(B_data)):
						B_mag = np.linalg.norm(B_data[idx])
						f.write(f"          {fmt.format(B_mag)}\n")
					f.write('        </DataArray>\n')

				elif coef == 'H_magnitude':
					# |H| [A/m]
					f.write('        <DataArray type="Float64" Name="H_magnitude" format="ascii">\n')
					for idx in range(len(H_data)):
						H_mag = np.linalg.norm(H_data[idx])
						f.write(f"          {fmt.format(H_mag)}\n")
					f.write('        </DataArray>\n')

				elif coef == 'Bx':
					f.write('        <DataArray type="Float64" Name="Bx" format="ascii">\n')
					for idx in range(len(B_data)):
						f.write(f"          {fmt.format(B_data[idx][0])}\n")
					f.write('        </DataArray>\n')

				elif coef == 'By':
					f.write('        <DataArray type="Float64" Name="By" format="ascii">\n')
					for idx in range(len(B_data)):
						f.write(f"          {fmt.format(B_data[idx][1])}\n")
					f.write('        </DataArray>\n')

				elif coef == 'Bz':
					f.write('        <DataArray type="Float64" Name="Bz" format="ascii">\n')
					for idx in range(len(B_data)):
						f.write(f"          {fmt.format(B_data[idx][2])}\n")
					f.write('        </DataArray>\n')

				elif coef == 'Hx':
					f.write('        <DataArray type="Float64" Name="Hx" format="ascii">\n')
					for idx in range(len(H_data)):
						f.write(f"          {fmt.format(H_data[idx][0])}\n")
					f.write('        </DataArray>\n')

				elif coef == 'Hy':
					f.write('        <DataArray type="Float64" Name="Hy" format="ascii">\n')
					for idx in range(len(H_data)):
						f.write(f"          {fmt.format(H_data[idx][1])}\n")
					f.write('        </DataArray>\n')

				elif coef == 'Hz':
					f.write('        <DataArray type="Float64" Name="Hz" format="ascii">\n')
					for idx in range(len(H_data)):
						f.write(f"          {fmt.format(H_data[idx][2])}\n")
					f.write('        </DataArray>\n')

			f.write('      </PointData>\n')
			f.write('    </Piece>\n')
			f.write('  </StructuredGrid>\n')
			f.write('</VTKFile>\n')

		total_points = nx * ny * nz
		print(f"Exported: {output_file} ({nx}x{ny}x{nz} = {total_points} points)")
		print(f"  Coefficients: {', '.join(self.coefs)}")
		print(f"  Coordinates: {unit_name} -> meters (VTK standard)")

		return output_file


def export_field_grid_vts(obj, grid_params, filename='field_grid',
                          coefs=None, floatsize='double'):
	"""
	Export magnetic field on structured 3D grid to VTS format.

	This is a convenience function that creates a RadiaVTKOutput object
	and calls Do() to write the file.

	Parameters:
		obj: Radia object handle (integer)
		grid_params: Dict defining the observation grid:
			'x_range': [x_min, x_max] in current Radia units
			'y_range': [y_min, y_max] in current Radia units
			'z_range': [z_min, z_max] in current Radia units
			'nx', 'ny', 'nz': Grid points per direction (default: 21)
		filename: Output filename without extension (default: 'field_grid')
		coefs: List of coefficients to export
			   Default: ['B', 'B_magnitude']
			   Available: 'B', 'H', 'B_magnitude', 'H_magnitude',
			              'Bx', 'By', 'Bz', 'Hx', 'Hy', 'Hz'
		floatsize: 'single' or 'double' precision (default: 'double')

	Returns:
		str: The filename that was written

	Example:
		import radia as rad
		from radia import export_field_grid_vts

		rad.FldUnits('m')

		# Create magnet
		vertices = [[-0.05,-0.05,-0.05], [0.05,-0.05,-0.05], ...]
		magnet = rad.ObjHexahedron(vertices, [0, 0, 954930])

		# Define observation grid
		grid = {
			'x_range': [-0.2, 0.2],
			'y_range': [-0.2, 0.2],
			'z_range': [0.1, 0.3],
			'nx': 41, 'ny': 41, 'nz': 21
		}

		# Export to VTS
		export_field_grid_vts(magnet, grid, 'B_field', coefs=['B', 'B_magnitude'])
		# -> B_field.vts (ParaView compatible)
	"""
	vtk = RadiaVTKOutput(
		obj=obj,
		coefs=coefs,
		filename=filename,
		floatsize=floatsize,
		grid_params=grid_params
	)
	return vtk.Do()


if __name__ == '__main__':
	"""
	Demo: Export Radia geometry and field to VTK format
	"""
	import sys
	import os

	# Add build directory to path
	sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'Release'))

	import radia as rad

	print("=" * 60)
	print("Radia VTK Export Demo")
	print("=" * 60)

	# Set units to meters
	rad.FldUnits('m')

	# Create a simple test geometry
	print("\nCreating test geometry...")

	# Rectangular magnet (0.03m x 0.03m x 0.01m, 1.2T equivalent magnetization)
	mag = rad.ObjRecMag([0, 0, 0], [0.03, 0.03, 0.01], [0, 0, 954930])

	# Export geometry to Legacy VTK
	output_file = 'radia_demo_geometry'
	print(f"\nExporting geometry to {output_file}.vtk...")
	exportGeometryToVTK(mag, output_file)

	# Export field grid to VTS
	print("\nExporting field grid to VTS...")
	grid_params = {
		'x_range': [-0.05, 0.05],
		'y_range': [-0.05, 0.05],
		'z_range': [0.02, 0.08],
		'nx': 21, 'ny': 21, 'nz': 15
	}
	export_field_grid_vts(mag, grid_params, 'radia_demo_field',
	                      coefs=['B', 'B_magnitude'])

	print("\n" + "=" * 60)
	print("Export complete!")
	print("\nTo view in ParaView:")
	print("  1. Open ParaView")
	print(f"  2. File -> Open -> {output_file}.vtk (geometry)")
	print("  3. File -> Open -> radia_demo_field.vts (field)")
	print("  4. Click 'Apply' in the Properties panel")
	print("=" * 60)
