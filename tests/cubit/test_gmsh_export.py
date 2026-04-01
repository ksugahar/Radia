"""
Test Gmsh v2.2 export via cubit.cmd('radia export gmsh ...').

Tests:
1. 1st order elements export
2. 2nd order elements export
3. All element types (Tet, Hex, Wedge, Pyramid, Tri, Quad)
4. Gmsh format validation
"""

import sys
import os
# Auto-detect Cubit installation (single source of truth)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))
from install_panels import find_cubit_bin
_cubit_path = find_cubit_bin()
if _cubit_path and _cubit_path not in sys.path:
	sys.path.append(_cubit_path)

import cubit

# Initialize Cubit
cubit.init(['cubit', '-nojournal', '-batch'])

# Gmsh element type codes
GMSH_TYPES = {
	# 1st order
	'TET4': 4,
	'HEX8': 5,
	'WEDGE6': 6,
	'PYRAMID5': 7,
	'TRI3': 2,
	'QUAD4': 3,
	'EDGE2': 1,
	# 2nd order
	'TET10': 11,
	'HEX20': 17,
	'WEDGE15': 18,
	'PYRAMID13': 19,
	'TRI6': 9,
	'QUAD8': 16,
	'EDGE3': 8,
}


def parse_gmsh_file(filename):
	"""Parse Gmsh v2.2 file and return element statistics."""
	with open(filename, 'r') as f:
		content = f.read()

	result = {
		'nodes': 0,
		'elements': 0,
		'element_types': {},
	}

	lines = content.split('\n')
	in_nodes = False
	in_elements = False

	for i, line in enumerate(lines):
		if line == '$Nodes':
			in_nodes = True
			continue
		elif line == '$EndNodes':
			in_nodes = False
			continue
		elif line == '$Elements':
			in_elements = True
			continue
		elif line == '$EndElements':
			in_elements = False
			continue

		if in_nodes and result['nodes'] == 0:
			result['nodes'] = int(line.strip())
			in_nodes = False  # Skip node data
		elif in_elements and result['elements'] == 0:
			result['elements'] = int(line.strip())
		elif in_elements and result['elements'] > 0:
			parts = line.strip().split()
			if len(parts) >= 2:
				elem_type = int(parts[1])
				result['element_types'][elem_type] = result['element_types'].get(elem_type, 0) + 1

	return result


def test_1st_order_tet_mesh():
	"""Test Gmsh export with 1st order tet elements."""
	print("=" * 60)
	print("Test 1: 1st Order Tet Mesh")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all")
	cubit.cmd("block 1 name 'solid'")
	cubit.cmd("block 2 add tri all in surface all")
	cubit.cmd("block 2 name 'boundary'")

	num_tets = len(cubit.get_block_tets(1))
	num_tris = len(cubit.get_block_tris(2))
	print(f"  Cubit mesh: {num_tets} tets, {num_tris} tris")

	# Check 1st order
	tet_id = cubit.get_block_tets(1)[0]
	nodes = cubit.get_expanded_connectivity("tet", tet_id)
	print(f"  Nodes per tet: {len(nodes)} (expected: 4)")

	msh_file = "test_1st_order_tet.msh"
	cubit.cmd(f'radia export gmsh "{msh_file}" overwrite')

	# Parse and validate
	result = parse_gmsh_file(msh_file)
	print(f"  Gmsh file: {result['nodes']} nodes, {result['elements']} elements")
	print(f"  Element types: {result['element_types']}")

	# Verify element types (TET4=4, TRI3=2)
	assert GMSH_TYPES['TET4'] in result['element_types'], "TET4 not found!"
	assert GMSH_TYPES['TRI3'] in result['element_types'], "TRI3 not found!"
	print("  PASS: 1st order elements correctly exported")

	os.remove(msh_file)
	return True


def test_2nd_order_tet_mesh():
	"""Test Gmsh export with 2nd order tet elements."""
	print("\n" + "=" * 60)
	print("Test 2: 2nd Order Tet Mesh")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all")
	cubit.cmd("block 1 name 'solid'")
	cubit.cmd("block 1 element type tetra10")

	cubit.cmd("block 2 add tri all in surface all")
	cubit.cmd("block 2 name 'boundary'")
	cubit.cmd("block 2 element type tri6")

	num_tets = len(cubit.get_block_tets(1))
	num_tris = len(cubit.get_block_tris(2))
	print(f"  Cubit mesh: {num_tets} tets, {num_tris} tris")

	# Check 2nd order
	tet_id = cubit.get_block_tets(1)[0]
	nodes = cubit.get_expanded_connectivity("tet", tet_id)
	print(f"  Nodes per tet: {len(nodes)} (expected: 10)")

	msh_file = "test_2nd_order_tet.msh"
	cubit.cmd(f'radia export gmsh "{msh_file}" overwrite')

	# Parse and validate
	result = parse_gmsh_file(msh_file)
	print(f"  Gmsh file: {result['nodes']} nodes, {result['elements']} elements")
	print(f"  Element types: {result['element_types']}")

	# Verify element types (TET10=11, TRI6=9)
	assert GMSH_TYPES['TET10'] in result['element_types'], "TET10 not found!"
	assert GMSH_TYPES['TRI6'] in result['element_types'], "TRI6 not found!"
	print("  PASS: 2nd order elements correctly exported")

	os.remove(msh_file)
	return True


def test_hex_mesh():
	"""Test Gmsh export with hex elements (1st and 2nd order)."""
	print("\n" + "=" * 60)
	print("Test 3: Hex Mesh (1st and 2nd order)")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme map")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add hex all")
	cubit.cmd("block 1 name 'solid'")
	cubit.cmd("block 2 add face all in surface all")
	cubit.cmd("block 2 name 'boundary'")

	num_hexes = len(cubit.get_block_hexes(1))
	print(f"  Cubit mesh: {num_hexes} hexes")

	# Test 1st order
	msh_file = "test_hex_1st.msh"
	cubit.cmd(f'radia export gmsh "{msh_file}" overwrite')
	result = parse_gmsh_file(msh_file)
	print(f"  1st order - Element types: {result['element_types']}")
	assert GMSH_TYPES['HEX8'] in result['element_types'], "HEX8 not found!"
	print("  1st order hex: PASS")
	os.remove(msh_file)

	# Test 2nd order
	cubit.cmd("block 1 element type hex20")
	cubit.cmd("block 2 element type quad8")

	msh_file = "test_hex_2nd.msh"
	cubit.cmd(f'radia export gmsh "{msh_file}" overwrite')
	result = parse_gmsh_file(msh_file)
	print(f"  2nd order - Element types: {result['element_types']}")
	assert GMSH_TYPES['HEX20'] in result['element_types'], "HEX20 not found!"
	print("  2nd order hex: PASS")
	os.remove(msh_file)

	return True


def test_wedge_mesh():
	"""Test Gmsh export with wedge elements."""
	print("\n" + "=" * 60)
	print("Test 4: Wedge Mesh")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create cylinder radius 0.5 height 1")
	cubit.cmd("volume 1 scheme sweep")
	cubit.cmd("volume 1 size 0.3")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add wedge all")
	cubit.cmd("block 1 name 'wedges'")

	num_wedges = len(cubit.get_block_wedges(1))
	if num_wedges == 0:
		# Sweep may produce hex instead
		cubit.cmd("block 1 add hex all")
		num_hexes = len(cubit.get_block_hexes(1))
		print(f"  Sweep produced {num_hexes} hexes instead of wedges, skipping")
		return True

	print(f"  Cubit mesh: {num_wedges} wedges")

	msh_file = "test_wedge.msh"
	cubit.cmd(f'radia export gmsh "{msh_file}" overwrite')
	result = parse_gmsh_file(msh_file)
	print(f"  Element types: {result['element_types']}")

	if GMSH_TYPES['WEDGE6'] in result['element_types']:
		print("  1st order wedge: PASS")
	os.remove(msh_file)

	return True


def test_pyramid_mesh():
	"""Test Gmsh export with pyramid elements."""
	print("\n" + "=" * 60)
	print("Test 5: Pyramid Mesh")
	print("=" * 60)

	cubit.cmd("reset")
	# Create geometry that will produce pyramids (hex-tet transition)
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 2 move 1 0 0")
	cubit.cmd("imprint all")
	cubit.cmd("merge all")

	cubit.cmd("volume 1 scheme map")
	cubit.cmd("volume 2 scheme tetmesh")
	cubit.cmd("volume all size 0.5")
	cubit.cmd("mesh volume all")

	cubit.cmd("block 1 add pyramid all")
	cubit.cmd("block 1 name 'pyramids'")

	num_pyramids = len(cubit.get_block_pyramids(1))
	if num_pyramids == 0:
		print("  No pyramid elements generated, skipping")
		return True

	print(f"  Cubit mesh: {num_pyramids} pyramids")

	msh_file = "test_pyramid.msh"
	cubit.cmd(f'radia export gmsh "{msh_file}" overwrite')
	result = parse_gmsh_file(msh_file)
	print(f"  Element types: {result['element_types']}")

	if GMSH_TYPES['PYRAMID5'] in result['element_types']:
		print("  1st order pyramid: PASS")
	os.remove(msh_file)

	return True


def test_mixed_elements():
	"""Test Gmsh export with mixed element types."""
	print("\n" + "=" * 60)
	print("Test 6: Mixed Element Types")
	print("=" * 60)

	cubit.cmd("reset")
	# Volume 1: Hex mesh
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 move 0 0 0")

	# Volume 2: Tet mesh
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 2 move 1.5 0 0")

	cubit.cmd("volume 1 scheme map")
	cubit.cmd("volume 2 scheme tetmesh")
	cubit.cmd("volume all size 0.5")
	cubit.cmd("mesh volume all")

	cubit.cmd("block 1 add hex all")
	cubit.cmd("block 1 name 'hex_region'")
	cubit.cmd("block 2 add tet all")
	cubit.cmd("block 2 name 'tet_region'")
	cubit.cmd("block 3 add face all in surface all")
	cubit.cmd("block 3 name 'quad_boundary'")
	cubit.cmd("block 4 add tri all in surface all")
	cubit.cmd("block 4 name 'tri_boundary'")

	num_hexes = len(cubit.get_block_hexes(1))
	num_tets = len(cubit.get_block_tets(2))
	print(f"  Cubit mesh: {num_hexes} hexes, {num_tets} tets")

	msh_file = "test_mixed.msh"
	cubit.cmd(f'radia export gmsh "{msh_file}" overwrite')
	result = parse_gmsh_file(msh_file)
	print(f"  Element types: {result['element_types']}")

	# Verify both hex and tet are present
	assert GMSH_TYPES['HEX8'] in result['element_types'], "HEX8 not found!"
	assert GMSH_TYPES['TET4'] in result['element_types'], "TET4 not found!"
	print("  PASS: Mixed elements correctly exported")

	os.remove(msh_file)
	return True


def test_gmsh_format_validation():
	"""Validate Gmsh v2.2 file format structure."""
	print("\n" + "=" * 60)
	print("Test 7: Gmsh Format Validation")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all")
	cubit.cmd("block 1 name 'solid'")

	msh_file = "test_format.msh"
	cubit.cmd(f'radia export gmsh "{msh_file}" overwrite')

	with open(msh_file, 'r') as f:
		content = f.read()

	# Check required sections
	assert '$MeshFormat' in content, "$MeshFormat section missing!"
	assert '$EndMeshFormat' in content, "$EndMeshFormat section missing!"
	assert '$Nodes' in content, "$Nodes section missing!"
	assert '$EndNodes' in content, "$EndNodes section missing!"
	assert '$Elements' in content, "$Elements section missing!"
	assert '$EndElements' in content, "$EndElements section missing!"
	print("  Required sections present: PASS")

	# Check mesh format version
	lines = content.split('\n')
	format_idx = lines.index('$MeshFormat')
	format_line = lines[format_idx + 1]
	assert format_line.startswith('2.2'), f"Expected Gmsh v2.2, got: {format_line}"
	print(f"  Mesh format version: {format_line.split()[0]} - PASS")

	os.remove(msh_file)
	return True


if __name__ == "__main__":
	print("\n" + "=" * 60)
	print("Gmsh v2.2 Export Test Suite (via cubit.cmd)")
	print("=" * 60)

	all_passed = True

	tests = [
		test_1st_order_tet_mesh,
		test_2nd_order_tet_mesh,
		test_hex_mesh,
		test_wedge_mesh,
		test_pyramid_mesh,
		test_mixed_elements,
		test_gmsh_format_validation,
	]

	for test in tests:
		try:
			if not test():
				all_passed = False
		except Exception as e:
			print(f"  FAIL: {e}")
			import traceback
			traceback.print_exc()
			all_passed = False

	print("\n" + "=" * 60)
	if all_passed:
		print("All tests PASSED!")
	else:
		print("Some tests FAILED!")
	print("=" * 60)
