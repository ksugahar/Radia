#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Magnetic yoke model from Nastran mesh

Creates magnetic yoke geometry for Radia from Nastran .nas file.
"""

import sys
import os
from pathlib import Path

# Add paths
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root / 'build' / 'lib' / 'Release'))
sys.path.insert(0, str(project_root / 'dist'))

import numpy as np
import radia as rad

from nastran_reader import read_nastran_mesh


def create_yoke_from_nastran(nas_file, material_properties=None):
	"""
	Create magnetic yoke from Nastran mesh file.

	Args:
		nas_file: Path to .nas file
		material_properties: Dictionary with material properties
			- permeability: [mu_x, mu_y] or scalar
			- saturation (optional): For nonlinear materials

	Returns:
		tuple: (yoke_object, mesh_info)
			- yoke_object: Radia object ID
			- mesh_info: Dictionary with mesh statistics
	"""
	print("=" * 70)
	print("Creating magnetic yoke from Nastran mesh...")
	print("=" * 70)

	# Read mesh
	mesh = read_nastran_mesh(nas_file)

	nodes = mesh['nodes']
	hex_elements_data = mesh['hex_elements']
	penta_elements_data = mesh['penta_elements']
	node_id_list = mesh['node_id_list']

	print(f"\nConverting to Radia polyhedra...")

	all_polyhedra = []

	# Hexahedron face connectivity (1-indexed)
	# Nastran CHEXA node numbering:
	# Nodes: G1, G2, G3, G4, G5, G6, G7, G8
	hex_faces = [
		[1, 2, 3, 4],  # Bottom face
		[5, 6, 7, 8],  # Top face
		[1, 2, 6, 5],  # Side face 1
		[2, 3, 7, 6],  # Side face 2
		[3, 4, 8, 7],  # Side face 3
		[4, 1, 5, 8],  # Side face 4
	]

	# Pentahedron face connectivity (1-indexed)
	# Nastran CPENTA node numbering:
	# Nodes: G1, G2, G3 (bottom triangle), G4, G5, G6 (top triangle)
	penta_faces = [
		[1, 2, 3],     # Bottom triangle
		[4, 5, 6],     # Top triangle
		[1, 2, 5, 4],  # Side face 1 (quad)
		[2, 3, 6, 5],  # Side face 2 (quad)
		[3, 1, 4, 6],  # Side face 3 (quad)
	]

	# Process hexahedra
	total_elements = len(hex_elements_data) + len(penta_elements_data)
	processed = 0

	for i, elem in enumerate(hex_elements_data):
		# Get node IDs for this element
		node_ids = elem

		# Get coordinates for each node
		coords = []
		for nid in node_ids:
			try:
				idx = node_id_list.index(nid)
				coords.append(list(nodes[idx]))
			except ValueError:
				print(f"Warning: Node {nid} not found in mesh")
				coords.append([0, 0, 0])

		# Create polyhedron
		try:
			poly = rad.ObjPolyhdr(coords, hex_faces)
			all_polyhedra.append(poly)
		except Exception as e:
			print(f"Warning: Failed to create hexahedron {i+1}: {e}")

		processed += 1
		if processed % 50 == 0:
			print(f"  Progress: {processed}/{total_elements} elements", end='\r')

	# Process pentahedra
	for i, elem in enumerate(penta_elements_data):
		# Get node IDs for this element
		node_ids = elem

		# Get coordinates for each node
		coords = []
		for nid in node_ids:
			try:
				idx = node_id_list.index(nid)
				coords.append(list(nodes[idx]))
			except ValueError:
				print(f"Warning: Node {nid} not found in mesh")
				coords.append([0, 0, 0])

		# Create polyhedron
		try:
			poly = rad.ObjPolyhdr(coords, penta_faces)
			all_polyhedra.append(poly)
		except Exception as e:
			print(f"Warning: Failed to create pentahedron {i+1}: {e}")

		processed += 1
		if processed % 50 == 0:
			print(f"  Progress: {processed}/{total_elements} elements", end='\r')

	print(f"  Progress: {total_elements}/{total_elements} elements")
	print(f"  [OK] Created {len(all_polyhedra)} polyhedra ({len(hex_elements_data)} hex + {len(penta_elements_data)} penta)")

	# Combine into container
	if not all_polyhedra:
		print("  [ERROR] No elements created")
		return None, None

	yoke = rad.ObjCnt(all_polyhedra)

	# Set visualization attributes (cyan, semi-transparent)
	rad.ObjDrwAtr(yoke, [0, 1, 1], 0.1)

	# Apply material properties
	if material_properties is None:
		# Default: Linear material with high permeability
		material_properties = {
			'permeability': 1000,  # Relative permeability
			'type': 'linear'
		}

	print(f"\nApplying material properties...")
	print(f"  Type: {material_properties.get('type', 'linear')}")

	if material_properties.get('type') == 'linear':
		mu_r = material_properties['permeability']
		if isinstance(mu_r, (list, tuple)):
			mat = rad.MatLin(mu_r, [0, 0, 1])
			print(f"  Permeability: μx={mu_r[0]}, μy={mu_r[1]}")
		else:
			mat = rad.MatLin([mu_r, mu_r], [0, 0, 1])
			print(f"  Permeability: μr={mu_r}")
	else:
		# For future: nonlinear materials
		mat = rad.MatLin([1000, 1000], [0, 0, 1])
		print(f"  Using default linear material")

	rad.MatApl(yoke, mat)

	# Note: Symmetries disabled for compatibility with coil geometry
	# rad.TrfZerPerp(yoke, [0, 0, 0], [0, 1, 0])  # Y-Z plane symmetry
	# rad.TrfZerPara(yoke, [0, 0, 0], [1, 0, 0])  # X=0 plane symmetry

	print(f"  [OK] Material applied")

	mesh_info = {
		'num_nodes': len(nodes),
		'num_hex_elements': len(hex_elements_data),
		'num_penta_elements': len(penta_elements_data),
		'num_total_elements': len(hex_elements_data) + len(penta_elements_data),
		'num_polyhedra': len(all_polyhedra)
	}

	return yoke, mesh_info


if __name__ == '__main__':
	"""Test yoke model creation"""
	import os

	nas_file = os.path.join(os.path.dirname(__file__), 'York.bdf')

	if os.path.exists(nas_file):
		print("=" * 70)
		print("YOKE MODEL TEST")
		print("=" * 70 + "\n")

		# Create yoke
		yoke, info = create_yoke_from_nastran(nas_file)

		if yoke:
			print("\n" + "=" * 70)
			print("Yoke Statistics:")
			print("=" * 70)
			print(f"  Nodes: {info['num_nodes']}")
			print(f"  Hexahedra: {info['num_hex_elements']}")
			print(f"  Pentahedra: {info['num_penta_elements']}")
			print(f"  Total elements: {info['num_total_elements']}")
			print(f"  Radia polyhedra: {info['num_polyhedra']}")

			# Get bounding box
			bbox = rad.ObjGeoLim(yoke)
			print(f"\nBounding box:")
			print(f"  X: [{bbox[0]:.2f}, {bbox[1]:.2f}] mm")
			print(f"  Y: [{bbox[2]:.2f}, {bbox[3]:.2f}] mm")
			print(f"  Z: [{bbox[4]:.2f}, {bbox[5]:.2f}] mm")

			# Test field calculation
			print("\n" + "=" * 70)
			print("Testing field calculation (yoke only, no current)...")
			print("=" * 70)
			B = rad.Fld(yoke, 'b', [0, 0, 0])
			print(f"  B at origin: {B} T")
			print("  (Should be ~zero with no excitation)")

			# VTK Export - Export geometry with same filename as script
			try:
				from radia_vtk_export import exportGeometryToVTK
				import os

				script_name = os.path.splitext(os.path.basename(__file__))[0]
				vtk_filename = f"{script_name}.vtk"
				vtk_path = os.path.join(os.path.dirname(__file__), vtk_filename)

				exportGeometryToVTK(yoke, vtk_path)
				print(f"\n[VTK] Exported: {vtk_filename}")
				print(f"      View with: paraview {vtk_filename}")
			except ImportError:
				print("\n[VTK] Warning: radia_vtk_export not available (VTK export skipped)")
			except Exception as e:
				print(f"\n[VTK] Warning: Export failed: {e}")

			# Cleanup
			rad.UtiDelAll()

			print("\n" + "=" * 70)
			print("[OK] Yoke model test complete")
			print("=" * 70)
	else:
		print(f"Error: File not found: {nas_file}")
