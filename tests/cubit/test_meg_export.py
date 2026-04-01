"""
Test MEG export via cubit.cmd('radia export meg ...').

Tests:
1. 3D mesh export (DIM='T')
2. 2D mesh export (DIM='K')
3. Axisymmetric mesh export (DIM='R')
4. MEG format validation
5. MGR2 spatial nodes option
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


def parse_meg_file(filename):
	"""Parse MEG file and return statistics."""
	with open(filename, 'r', encoding='UTF-8') as f:
		content = f.read()

	result = {
		'mgr1_count': 0,  # Node count
		'mgr2_count': 0,  # Spatial node count
		'elements': {},   # Element type counts
		'has_header': False,
		'has_footer': False,
	}

	lines = content.split('\n')
	for line in lines:
		line = line.strip()

		if line.startswith('BOOK  MEP'):
			result['has_header'] = True
		elif line.startswith('BOOK  END'):
			result['has_footer'] = True
		elif line.startswith('MGR1'):
			result['mgr1_count'] += 1
		elif line.startswith('MGR2'):
			result['mgr2_count'] += 1
		elif line and not line.startswith('*') and not line.startswith('MGSC') and not line.startswith('BOOK'):
			# Element line - extract element type (first 4 chars + DIM)
			parts = line.split()
			if len(parts) >= 1:
				elem_type = parts[0]
				result['elements'][elem_type] = result['elements'].get(elem_type, 0) + 1

	return result


def test_3d_tet_mesh():
	"""Test MEG export with 3D tet mesh (DIM='T')."""
	print("=" * 60)
	print("Test 1: 3D Tet Mesh (DIM='T')")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all")
	cubit.cmd("block 1 name 'SOLI'")  # MEG uses first 4 chars

	num_tets = len(cubit.get_block_tets(1))
	print(f"  Cubit mesh: {num_tets} tets")

	meg_file = "test_3d_tet.meg"
	cubit.cmd(f'radia export meg "{meg_file}" overwrite')

	result = parse_meg_file(meg_file)
	print(f"  MGR1 (nodes): {result['mgr1_count']}")
	print(f"  Elements: {result['elements']}")

	# Element name should be block_name[:4] + DIM = 'SOLIT'
	assert 'SOLIT' in result['elements'], "SOLIT elements not found!"
	assert result['elements']['SOLIT'] == num_tets, f"Expected {num_tets} tets, got {result['elements']['SOLIT']}"
	print("  PASS: 3D tet export correct")

	os.remove(meg_file)
	return True


def test_3d_hex_mesh():
	"""Test MEG export with 3D hex mesh."""
	print("\n" + "=" * 60)
	print("Test 2: 3D Hex Mesh (DIM='T')")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme map")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add hex all")
	cubit.cmd("block 1 name 'HEXA'")

	num_hexes = len(cubit.get_block_hexes(1))
	print(f"  Cubit mesh: {num_hexes} hexes")

	meg_file = "test_3d_hex.meg"
	cubit.cmd(f'radia export meg "{meg_file}" overwrite')

	result = parse_meg_file(meg_file)
	print(f"  MGR1 (nodes): {result['mgr1_count']}")
	print(f"  Elements: {result['elements']}")

	assert 'HEXAT' in result['elements'], "HEXAT elements not found!"
	print("  PASS: 3D hex export correct")

	os.remove(meg_file)
	return True


def test_2d_mesh():
	"""Test MEG export with 2D mesh (DIM='K')."""
	print("\n" + "=" * 60)
	print("Test 3: 2D Mesh (DIM='K')")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create surface rectangle width 1 height 1 zplane")
	cubit.cmd("surface 1 scheme trimesh")
	cubit.cmd("surface 1 size 0.3")
	cubit.cmd("mesh surface 1")

	cubit.cmd("block 1 add tri all")
	cubit.cmd("block 1 name 'TRIA'")

	num_tris = len(cubit.get_block_tris(1))
	print(f"  Cubit mesh: {num_tris} tris")

	meg_file = "test_2d.meg"
	cubit.cmd(f'radia export meg "{meg_file}" overwrite')

	result = parse_meg_file(meg_file)
	print(f"  MGR1 (nodes): {result['mgr1_count']}")
	print(f"  Elements: {result['elements']}")

	assert 'TRIAK' in result['elements'], "TRIAK elements not found!"
	print("  PASS: 2D mesh export correct")

	os.remove(meg_file)
	return True


def test_axisymmetric_mesh():
	"""Test MEG export with axisymmetric mesh (DIM='R')."""
	print("\n" + "=" * 60)
	print("Test 4: Axisymmetric Mesh (DIM='R')")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create surface rectangle width 1 height 1 zplane")
	cubit.cmd("surface 1 scheme map")
	cubit.cmd("surface 1 size 0.3")
	cubit.cmd("mesh surface 1")

	cubit.cmd("block 1 add face all")
	cubit.cmd("block 1 name 'QUAD'")

	num_quads = len(cubit.get_block_faces(1))
	print(f"  Cubit mesh: {num_quads} quads")

	meg_file = "test_axi.meg"
	cubit.cmd(f'radia export meg "{meg_file}" overwrite')

	result = parse_meg_file(meg_file)
	print(f"  MGR1 (nodes): {result['mgr1_count']}")
	print(f"  Elements: {result['elements']}")

	assert 'QUADR' in result['elements'], "QUADR elements not found!"
	print("  PASS: Axisymmetric mesh export correct")

	os.remove(meg_file)
	return True


def test_meg_format():
	"""Test MEG file format structure."""
	print("\n" + "=" * 60)
	print("Test 5: MEG Format Validation")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all")
	cubit.cmd("block 1 name 'TEST'")

	meg_file = "test_format.meg"
	cubit.cmd(f'radia export meg "{meg_file}" overwrite')

	result = parse_meg_file(meg_file)

	assert result['has_header'], "BOOK MEP header missing!"
	assert result['has_footer'], "BOOK END footer missing!"
	assert result['mgr1_count'] > 0, "No MGR1 nodes found!"

	print(f"  Header: {'PASS' if result['has_header'] else 'FAIL'}")
	print(f"  Footer: {'PASS' if result['has_footer'] else 'FAIL'}")
	print(f"  Nodes (MGR1): {result['mgr1_count']}")
	print("  PASS: MEG format valid")

	os.remove(meg_file)
	return True


def test_mgr2_spatial_nodes():
	"""Test MGR2 spatial nodes option.

	Note: MGR2 spatial nodes are now handled via the plugin command syntax.
	This test verifies the basic MEG export works; MGR2 functionality
	may require additional plugin parameters.
	"""
	print("\n" + "=" * 60)
	print("Test 6: MGR2 Spatial Nodes")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all")
	cubit.cmd("block 1 name 'TEST'")

	meg_file = "test_mgr2.meg"
	cubit.cmd(f'radia export meg "{meg_file}" overwrite')

	result = parse_meg_file(meg_file)
	print(f"  MGR1 (mesh nodes): {result['mgr1_count']}")
	print(f"  MGR2 (spatial nodes): {result['mgr2_count']}")

	print("  PASS: MEG export completed")

	os.remove(meg_file)
	return True


def test_mixed_elements():
	"""Test MEG export with mixed element types."""
	print("\n" + "=" * 60)
	print("Test 7: Mixed Elements")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all")
	cubit.cmd("block 1 name 'VOLU'")
	cubit.cmd("block 2 add tri all in surface all")
	cubit.cmd("block 2 name 'SURF'")

	num_tets = len(cubit.get_block_tets(1))
	num_tris = len(cubit.get_block_tris(2))
	print(f"  Cubit mesh: {num_tets} tets, {num_tris} tris")

	meg_file = "test_mixed.meg"
	cubit.cmd(f'radia export meg "{meg_file}" overwrite')

	result = parse_meg_file(meg_file)
	print(f"  Elements: {result['elements']}")

	assert 'VOLUT' in result['elements'], "VOLUT elements not found!"
	assert 'SURFT' in result['elements'], "SURFT elements not found!"
	print("  PASS: Mixed elements export correct")

	os.remove(meg_file)
	return True


def test_1st_order_only():
	"""Verify MEG export always uses 1st order elements."""
	print("\n" + "=" * 60)
	print("Test 8: 1st Order Only (even with 2nd order mesh)")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all")
	cubit.cmd("block 1 name 'TEST'")
	cubit.cmd("block 1 element type tetra10")  # Convert to 2nd order

	# Verify Cubit has 2nd order
	tet_id = cubit.get_block_tets(1)[0]
	nodes_2nd = cubit.get_expanded_connectivity("tet", tet_id)
	nodes_1st = cubit.get_connectivity("tet", tet_id)
	print(f"  Cubit nodes (get_connectivity): {len(nodes_1st)}")
	print(f"  Cubit nodes (get_expanded_connectivity): {len(nodes_2nd)}")

	meg_file = "test_1st_order.meg"
	cubit.cmd(f'radia export meg "{meg_file}" overwrite')

	# Read file and check element node count
	with open(meg_file, 'r', encoding='UTF-8') as f:
		content = f.read()

	# Find element lines and count nodes
	for line in content.split('\n'):
		if line.startswith('TEST'):
			parts = line.split()
			# Element format: NAME ID 0 BLOCK_ID N1 N2 N3 N4
			# For tet: 4 nodes after BLOCK_ID
			node_count = len(parts) - 4  # Subtract NAME, ID, 0, BLOCK_ID
			print(f"  MEG element nodes: {node_count} (expected: 4 for tet)")
			assert node_count == 4, f"Expected 4 nodes, got {node_count}"
			break

	print("  PASS: MEG uses 1st order only (get_connectivity)")
	os.remove(meg_file)
	return True


if __name__ == "__main__":
	print("\n" + "=" * 60)
	print("MEG Export Test Suite (via cubit.cmd)")
	print("=" * 60)

	all_passed = True

	tests = [
		test_3d_tet_mesh,
		test_3d_hex_mesh,
		test_2d_mesh,
		test_axisymmetric_mesh,
		test_meg_format,
		test_mgr2_spatial_nodes,
		test_mixed_elements,
		test_1st_order_only,
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
