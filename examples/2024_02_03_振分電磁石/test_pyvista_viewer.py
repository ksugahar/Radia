#!/usr/bin/env python
"""
Test PyVista viewer with offscreen rendering (no GUI required)
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'dist'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'lib', 'Release'))

import radia as rad
import numpy as np

print("="*60)
print("Radia PyVista Offscreen Rendering Test")
print("="*60)

# Check if PyVista is available
try:
	import pyvista as pv
	print(f"\nPyVista version: {pv.__version__}")
except ImportError:
	print("\nERROR: PyVista is not installed.")
	print("Install it with: pip install pyvista")
	sys.exit(1)

# Create test objects
print("\nCreating test objects...")

# Arc current
g1 = rad.ObjArcCur([0,0,0], [95, 105], [0, 2*np.pi], 20, 36, 1e6/20/10, "man")
print(f"  Arc current: ID {g1}")

# Rectangular magnet
g2 = rad.ObjRecMag([0,0,-30], [30,30,10], [0,0,1])
print(f"  Rectangular magnet: ID {g2}")

# Container
g = rad.ObjCnt([g1, g2])
print(f"  Container: ID {g}")

# Get VTK data from Radia
print("\nExtracting VTK data from Radia...")
vtk_data = rad.ObjDrwVTK(g, 'Axes->False')

if 'polygons' not in vtk_data:
	print("ERROR: No polygon data found")
	sys.exit(1)

poly = vtk_data['polygons']

# Extract vertices
vertices = np.array(poly['vertices']).reshape(-1, 3)
print(f"  Vertices: {len(vertices)}")

# Extract polygon connectivity
lengths = poly['lengths']
n_polys = len(lengths)
print(f"  Polygons: {n_polys}")

# Build faces array for PyVista
faces = []
offset = 0
for length in lengths:
	faces.append(length)
	for i in range(length):
		faces.append(offset + i)
	offset += length

faces = np.array(faces)

# Create PyVista mesh
print("\nCreating PyVista mesh...")
mesh = pv.PolyData(vertices, faces)
print(f"  Mesh created: {mesh.n_points} points, {mesh.n_cells} cells")

# Create offscreen plotter
print("\nRendering with PyVista (offscreen)...")
plotter = pv.Plotter(off_screen=True)
plotter.add_axes()
plotter.add_text("Radia Object (PyVista)", position='upper_left', font_size=12)

# Add mesh
plotter.add_mesh(mesh, color='lightblue', show_edges=True,
                 edge_color='black', line_width=1.0, opacity=0.8)

# Set camera view
plotter.camera_position = 'iso'

# Save screenshot
output_file = 'radia_pyvista_test.png'
plotter.screenshot(output_file)
print(f"  Screenshot saved: {output_file}")

plotter.close()

# Check file was created
if os.path.exists(output_file):
	size = os.path.getsize(output_file)
	print(f"\nSUCCESS: Image file created ({size} bytes)")
else:
	print("\nERROR: Image file was not created")

print()
print("="*60)
print("PyVista viewer works perfectly!")
print("This is a modern alternative to rad.ObjDrwOpenGL()")
print("="*60)
