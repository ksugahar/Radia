#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Nastran mesh reader for Radia

Reads Nastran .nas/.bdf files and converts to Python data structures.
Supports GRID (nodes), CHEXA (hexahedron), and CPENTA (pentahedron) elements.
"""

import re
import numpy as np


def read_nastran_mesh(filename):
	"""
	Read Nastran mesh file (.nas/.bdf format).

	Args:
		filename: Path to .nas or .bdf file

	Returns:
		dict: Dictionary with mesh data
			- nodes: numpy array (N, 3) with node coordinates
			- hex_elements: list of hexahedra [n1, n2, n3, n4, n5, n6, n7, n8]
			- penta_elements: list of pentahedra [n1, n2, n3, n4, n5, n6]
			- node_ids: dict mapping node_id to array index
			- node_id_list: list of node IDs in order
	"""
	print(f"Reading Nastran mesh: {filename}")

	with open(filename, 'r') as f:
		lines = f.readlines()

	# Parse GRID entries
	nodes = {}  # node_id -> [x, y, z]
	hex_elements = []  # List of [n1, n2, n3, n4, n5, n6, n7, n8]
	penta_elements = []  # List of [n1, n2, n3, n4, n5, n6]

	i = 0
	while i < len(lines):
		line = lines[i].rstrip('\n')

		# Parse GRID lines (Fixed format: 8 characters per field)
		if line.startswith('GRID'):
			try:
				# GRID    ID      CP      X1      X2      X3
				# 01234567890123456789012345678901234567890123456789
				node_id = int(line[8:16])
				x = float(line[24:32])
				y = float(line[32:40])
				z = float(line[40:48])
				nodes[node_id] = [x, y, z]
			except ValueError as e:
				print(f"Warning: Failed to parse GRID line {i+1}: {e}")

		# Parse CHEXA lines (Hexahedron with continuation)
		elif line.startswith('CHEXA'):
			try:
				# CHEXA   EID     PID     G1      G2      G3      G4      G5      G6      +
				# +       G7      G8
				# Get nodes from first line (6 nodes)
				n1 = int(line[24:32])
				n2 = int(line[32:40])
				n3 = int(line[40:48])
				n4 = int(line[48:56])
				n5 = int(line[56:64])
				n6 = int(line[64:72])

				# Get continuation line
				i += 1
				cont_line = lines[i].rstrip('\n')
				# +       G7      G8
				n7 = int(cont_line[8:16])
				n8 = int(cont_line[16:24])

				hex_elements.append([n1, n2, n3, n4, n5, n6, n7, n8])
			except (ValueError, IndexError) as e:
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
				print(f"Warning: Failed to parse CPENTA at line {i}: {e}")

		i += 1

	# Convert to numpy arrays
	node_ids = sorted(nodes.keys())
	node_id_to_idx = {nid: idx for idx, nid in enumerate(node_ids)}

	nodes_array = np.array([nodes[nid] for nid in node_ids])
	hex_array = np.array(hex_elements, dtype=int) if hex_elements else np.array([])
	penta_array = np.array(penta_elements, dtype=int) if penta_elements else np.array([])

	print(f"  Nodes: {len(nodes_array)}")
	print(f"  Elements (CHEXA): {len(hex_elements)}")
	print(f"  Elements (CPENTA): {len(penta_elements)}")
	print(f"  Total elements: {len(hex_elements) + len(penta_elements)}")

	return {
		'nodes': nodes_array,
		'hex_elements': hex_array,
		'penta_elements': penta_array,
		'node_ids': node_id_to_idx,
		'node_id_list': node_ids
	}


if __name__ == '__main__':
	"""Test the Nastran reader"""
	import os

	nas_file = os.path.join(os.path.dirname(__file__), 'York.bdf')

	if os.path.exists(nas_file):
		mesh = read_nastran_mesh(nas_file)

		print("\nMesh statistics:")
		print(f"  Total nodes: {len(mesh['nodes'])}")
		print(f"  Hexahedra (CHEXA): {len(mesh['hex_elements'])}")
		print(f"  Pentahedra (CPENTA): {len(mesh['penta_elements'])}")
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
	else:
		print(f"Error: File not found: {nas_file}")
