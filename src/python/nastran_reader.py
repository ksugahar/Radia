#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Nastran Mesh Reader for Radia

Reads Nastran .nas/.bdf files and converts to Python data structures.
Supports GRID (nodes), CHEXA (hexahedron), CPENTA (pentahedron), CTETRA (tetrahedron), and CTRIA3 (triangle) elements.

For surface meshes (CTRIA3), triangles with the same material ID (PID) are grouped into a single polyhedron,
which is useful for linear magnetic analysis where only the surface representation is needed.

Usage:
	from nastran_reader import read_nastran_mesh

	mesh = read_nastran_mesh('sphere.bdf')
	nodes = mesh['nodes']  # numpy array (N, 3)
	hex_elements = mesh['hex_elements']  # list of [n1, ..., n8]
	penta_elements = mesh['penta_elements']  # list of [n1, ..., n6]
	tetra_elements = mesh['tetra_elements']  # list of [n1, ..., n4]
	tria_groups = mesh['tria_groups']  # dict {material_id: {'faces': [[n1,n2,n3], ...], 'nodes': set(...)}}

Date: 2025-11-01
"""

import re
import numpy as np


def read_nastran_mesh(filename, verbose=True):
	"""
	Read Nastran mesh file (.nas/.bdf format).

	Args:
		filename: Path to .nas or .bdf file
		verbose: Print progress messages (default: True)

	Returns:
		dict: Dictionary with mesh data
			- nodes: numpy array (N, 3) with node coordinates [mm]
			- hex_elements: list of hexahedra [n1, n2, n3, n4, n5, n6, n7, n8]
			- penta_elements: list of pentahedra [n1, n2, n3, n4, n5, n6]
			- tetra_elements: list of tetrahedra [n1, n2, n3, n4]
			- tria_groups: dict {material_id: {'faces': [[n1,n2,n3], ...], 'nodes': set(node_ids)}}
			- node_ids: dict mapping node_id to array index
			- node_id_list: list of node IDs in order
	"""
	if verbose:
		print(f"Reading Nastran mesh: {filename}")

	with open(filename, 'r') as f:
		lines = f.readlines()

	# Parse GRID entries
	nodes = {}  # node_id -> [x, y, z]
	hex_elements = []  # List of [n1, n2, n3, n4, n5, n6, n7, n8]
	penta_elements = []  # List of [n1, n2, n3, n4, n5, n6]
	tetra_elements = []  # List of [n1, n2, n3, n4]
	tria_groups = {}  # material_id -> {'faces': [[n1,n2,n3], ...], 'nodes': set(node_ids)}

	i = 0
	while i < len(lines):
		line = lines[i].rstrip('\n')

		# Parse GRID* lines (Extended format with continuation)
		if line.startswith('GRID*'):
			try:
				# GRID*   ID              CP              X1              X2
				# *       X3
				# Field width: 16 characters
				node_id = int(line[8:24].strip())
				# Skip CP field (line[24:40])
				x = float(line[40:56].strip())
				y = float(line[56:72].strip())

				# Get continuation line for Z coordinate
				i += 1
				cont_line = lines[i].rstrip('\n')
				z = float(cont_line[8:24].strip())

				nodes[node_id] = [x, y, z]
			except (ValueError, IndexError) as e:
				if verbose:
					print(f"Warning: Failed to parse GRID* line {i}: {e}")

		# Parse GRID lines (Fixed format: 8 characters per field)
		elif line.startswith('GRID') and not line.startswith('GRID*'):
			try:
				# GRID    ID      CP      X1      X2      X3
				# 01234567890123456789012345678901234567890123456789
				node_id = int(line[8:16])
				x = float(line[24:32])
				y = float(line[32:40])
				z = float(line[40:48])
				nodes[node_id] = [x, y, z]
			except ValueError as e:
				if verbose:
					print(f"Warning: Failed to parse GRID line {i+1}: {e}")

		# Parse CHEXA lines (Hexahedron with continuation)
		elif line.startswith('CHEXA'):
			try:
				# CHEXA   EID     PID     G1      G2      G3      G4      G5      G6      +
				# +              G7      G8
				# Get nodes from first line (6 nodes)
				n1 = int(line[24:32].strip())
				n2 = int(line[32:40].strip())
				n3 = int(line[40:48].strip())
				n4 = int(line[48:56].strip())
				n5 = int(line[56:64].strip())
				n6 = int(line[64:72].strip())

				# Get continuation line
				i += 1
				cont_line = lines[i].rstrip('\n')
				# +              G7      G8
				n7 = int(cont_line[14:22].strip())
				n8 = int(cont_line[22:30].strip())

				hex_elements.append([n1, n2, n3, n4, n5, n6, n7, n8])
			except (ValueError, IndexError) as e:
				if verbose:
					print(f"Warning: Failed to parse CHEXA at line {i}: {e}")

		# Parse CPENTA lines (Pentahedron, single line)
		elif line.startswith('CPENTA'):
			try:
				# CPENTA  EID     PID     G1      G2      G3      G4      G5      G6
				# All 6 nodes on one line
				n1 = int(line[24:32])
				n2 = int(line[32:40])
				n3 = int(line[40:48])
				n4 = int(line[48:56])
				n5 = int(line[56:64])
				n6 = int(line[64:72])

				penta_elements.append([n1, n2, n3, n4, n5, n6])
			except (ValueError, IndexError) as e:
				if verbose:
					print(f"Warning: Failed to parse CPENTA at line {i}: {e}")

		# Parse CTETRA lines (Tetrahedron, single line)
		elif line.startswith('CTETRA'):
			try:
				# CTETRA  EID     PID     G1      G2      G3      G4
				# All 4 nodes on one line (10-node tetra has continuation, but we only use first 4)
				n1 = int(line[24:32])
				n2 = int(line[32:40])
				n3 = int(line[40:48])
				n4 = int(line[48:56])

				tetra_elements.append([n1, n2, n3, n4])
			except (ValueError, IndexError) as e:
				if verbose:
					print(f"Warning: Failed to parse CTETRA at line {i}: {e}")

		# Parse CTRIA3 lines (Triangle, single line)
		elif line.startswith('CTRIA3'):
			try:
				# CTRIA3  EID     PID     G1      G2      G3
				# Element ID, Property ID (material), and 3 node IDs
				element_id = int(line[8:16].strip())
				material_id = int(line[16:24].strip())
				n1 = int(line[24:32].strip())
				n2 = int(line[32:40].strip())
				n3 = int(line[40:48].strip())

				# Group triangles by material ID
				if material_id not in tria_groups:
					tria_groups[material_id] = {'faces': [], 'nodes': set()}

				tria_groups[material_id]['faces'].append([n1, n2, n3])
				tria_groups[material_id]['nodes'].update([n1, n2, n3])
			except (ValueError, IndexError) as e:
				if verbose:
					print(f"Warning: Failed to parse CTRIA3 at line {i}: {e}")

		i += 1

	# Convert to numpy arrays
	node_ids = sorted(nodes.keys())
	node_id_to_idx = {nid: idx for idx, nid in enumerate(node_ids)}

	nodes_array = np.array([nodes[nid] for nid in node_ids])
	hex_array = np.array(hex_elements, dtype=int) if hex_elements else np.array([])
	penta_array = np.array(penta_elements, dtype=int) if penta_elements else np.array([])
	tetra_array = np.array(tetra_elements, dtype=int) if tetra_elements else np.array([])

	if verbose:
		print(f"  Nodes: {len(nodes_array)}")
		print(f"  Elements (CHEXA): {len(hex_elements)}")
		print(f"  Elements (CPENTA): {len(penta_elements)}")
		print(f"  Elements (CTETRA): {len(tetra_elements)}")
		if tria_groups:
			total_trias = sum(len(group['faces']) for group in tria_groups.values())
			print(f"  Elements (CTRIA3): {total_trias} triangles in {len(tria_groups)} material group(s)")
			for mat_id, group in tria_groups.items():
				print(f"    Material {mat_id}: {len(group['faces'])} triangles, {len(group['nodes'])} unique nodes")
		total_elements = len(hex_elements) + len(penta_elements) + len(tetra_elements)
		if not tria_groups:
			print(f"  Total elements: {total_elements}")

	return {
		'nodes': nodes_array,
		'hex_elements': hex_array,
		'penta_elements': penta_array,
		'tetra_elements': tetra_array,
		'tria_groups': tria_groups,
		'node_ids': node_id_to_idx,
		'node_id_list': node_ids
	}


# Element face connectivity for Radia ObjPolyhdr
# Node numbering is 1-indexed (as used in connectivity arrays)

# Hexahedron face connectivity (1-indexed)
# Nastran CHEXA node numbering:
# Nodes: G1, G2, G3, G4 (bottom), G5, G6, G7, G8 (top)
HEX_FACES = [
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
PENTA_FACES = [
	[1, 2, 3],     # Bottom triangle
	[4, 5, 6],     # Top triangle
	[1, 2, 5, 4],  # Side face 1 (quad)
	[2, 3, 6, 5],  # Side face 2 (quad)
	[3, 1, 4, 6],  # Side face 3 (quad)
]

# Tetrahedron face connectivity (1-indexed)
# Nastran CTETRA node numbering:
# Nodes: G1, G2, G3, G4
TETRA_FACES = [
	[1, 2, 3],  # Face 1
	[1, 4, 2],  # Face 2
	[2, 4, 3],  # Face 3
	[3, 4, 1],  # Face 4
]


if __name__ == '__main__':
	"""Test the Nastran reader"""
	import os
	import sys

	# Try to find a test .bdf file
	test_files = [
		'../../examples/electromagnet/York.bdf',
		'../../examples/NGSolve_CoefficientFunction_to_Radia_BackgroundField/sphere.bdf',
	]

	nas_file = None
	for f in test_files:
		if os.path.exists(f):
			nas_file = f
			break

	if nas_file:
		mesh = read_nastran_mesh(nas_file)

		print("\nMesh statistics:")
		print(f"  Total nodes: {len(mesh['nodes'])}")
		print(f"  Hexahedra (CHEXA): {len(mesh['hex_elements'])}")
		print(f"  Pentahedra (CPENTA): {len(mesh['penta_elements'])}")
		print(f"  Tetrahedra (CTETRA): {len(mesh['tetra_elements'])}")
		print(f"\nFirst 5 nodes:")
		for i in range(min(5, len(mesh['nodes']))):
			print(f"    Node {mesh['node_id_list'][i]}: {mesh['nodes'][i]}")

		if len(mesh['hex_elements']) > 0:
			print(f"\nFirst 3 hexahedra (node IDs):")
			for i in range(min(3, len(mesh['hex_elements']))):
				print(f"    Hex {i+1}: {mesh['hex_elements'][i]}")

		if len(mesh['penta_elements']) > 0:
			print(f"\nFirst 3 pentahedra (node IDs):")
			for i in range(min(3, len(mesh['penta_elements']))):
				print(f"    Penta {i+1}: {mesh['penta_elements'][i]}")

		if len(mesh['tetra_elements']) > 0:
			print(f"\nFirst 3 tetrahedra (node IDs):")
			for i in range(min(3, len(mesh['tetra_elements']))):
				print(f"    Tetra {i+1}: {mesh['tetra_elements'][i]}")
	else:
		print(f"Error: No test .bdf files found")
		print(f"Searched for: {test_files}")
