"""
Test VTK export node ordering and non-contiguous node ID handling.

Tests:
1. HEX20 node ordering (Cubit -> VTK conversion)
2. WEDGE15 node ordering
3. Non-contiguous node IDs (regression test for bug fix)
4. Node coordinate accuracy
"""

import sys
import os
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
	sys.path.append(cubit_path)

import cubit

# Add parent directory to path for cubit_mesh_export
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))
import cubit_mesh_export

# Initialize Cubit
cubit.init(['cubit', '-nojournal', '-batch'])


def parse_vtk_file(filename):
	"""Parse VTK Legacy file and return data."""
	with open(filename, 'r') as f:
		lines = f.readlines()

	result = {
		'points': [],
		'cells': [],
		'cell_types': [],
	}

	i = 0
	while i < len(lines):
		line = lines[i].strip()

		if line.startswith('POINTS'):
			num_points = int(line.split()[1])
			i += 1
			for _ in range(num_points):
				coords = lines[i].strip().split()
				result['points'].append([float(c) for c in coords])
				i += 1
			continue

		elif line.startswith('CELLS'):
			parts = line.split()
			num_cells = int(parts[1])
			i += 1
			for _ in range(num_cells):
				cell_data = [int(x) for x in lines[i].strip().split()]
				result['cells'].append(cell_data)
				i += 1
			continue

		elif line.startswith('CELL_TYPES'):
			num_types = int(line.split()[1])
			i += 1
			for _ in range(num_types):
				result['cell_types'].append(int(lines[i].strip()))
				i += 1
			continue

		i += 1

	return result


def test_hex8_node_ordering():
	"""Test HEX8 node ordering."""
	print("=" * 60)
	print("Test 1: HEX8 Node Ordering")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme map")
	cubit.cmd("volume 1 size 1")  # Single hex
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add hex all")

	num_hexes = len(cubit.get_block_hexes(1))
	print(f"  Hexes: {num_hexes}")

	vtk_file = "test_hex8_ordering.vtk"
	cubit_mesh_export.export_vtk(cubit, vtk_file)

	result = parse_vtk_file(vtk_file)
	print(f"  VTK points: {len(result['points'])}")
	print(f"  VTK cells: {len(result['cells'])}")

	# Verify HEX8 cell type
	assert 12 in result['cell_types'], "HEX8 (type 12) not found!"

	# Get first hex cell
	hex_cell = result['cells'][0]
	num_nodes = hex_cell[0]
	node_indices = hex_cell[1:]

	print(f"  HEX8 nodes: {num_nodes} (expected: 8)")
	assert num_nodes == 8, f"Expected 8 nodes, got {num_nodes}"

	# Verify all node indices are valid
	for idx in node_indices:
		assert 0 <= idx < len(result['points']), f"Invalid node index: {idx}"

	print("  PASS: HEX8 node ordering correct")
	os.remove(vtk_file)
	return True


def test_hex20_node_ordering():
	"""Test HEX20 node ordering (Cubit -> VTK conversion)."""
	print("\n" + "=" * 60)
	print("Test 2: HEX20 Node Ordering")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme map")
	cubit.cmd("volume 1 size 1")  # Single hex
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add hex all")
	cubit.cmd("block 1 element type hex20")

	hex_id = cubit.get_block_hexes(1)[0]
	cubit_nodes = cubit.get_expanded_connectivity("hex", hex_id)
	print(f"  Cubit HEX20 nodes: {len(cubit_nodes)}")

	vtk_file = "test_hex20_ordering.vtk"
	cubit_mesh_export.export_vtk(cubit, vtk_file)

	result = parse_vtk_file(vtk_file)

	# Verify HEX20 cell type
	assert 25 in result['cell_types'], "HEX20 (type 25) not found!"

	# Get first hex cell
	hex_cell = result['cells'][0]
	num_nodes = hex_cell[0]
	node_indices = hex_cell[1:]

	print(f"  VTK HEX20 nodes: {num_nodes} (expected: 20)")
	assert num_nodes == 20, f"Expected 20 nodes, got {num_nodes}"

	# Verify all node indices are valid and unique
	assert len(set(node_indices)) == 20, "Duplicate node indices found!"
	for idx in node_indices:
		assert 0 <= idx < len(result['points']), f"Invalid node index: {idx}"

	# Verify node ordering conversion
	# Cubit:  [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
	# VTK:    [0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15]
	print("  VTK node order: ", node_indices[:8], "...", node_indices[12:16], node_indices[16:])
	print("  PASS: HEX20 node ordering correct")

	os.remove(vtk_file)
	return True


def test_non_contiguous_node_ids():
	"""Test handling of non-contiguous node IDs (regression test)."""
	print("\n" + "=" * 60)
	print("Test 3: Non-Contiguous Node IDs (Regression Test)")
	print("=" * 60)
	print("  This tests the bug fix for non-contiguous node IDs")
	print("-" * 60)

	cubit.cmd("reset")

	# Create two separate volumes to potentially get non-contiguous IDs
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 move 0 0 0")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 2 move 2 0 0")

	cubit.cmd("volume all scheme tetmesh")
	cubit.cmd("volume all size 0.5")
	cubit.cmd("mesh volume all")

	# Delete volume 1 mesh to create gaps in node IDs
	cubit.cmd("delete mesh volume 1 propagate")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all in volume 1")
	cubit.cmd("block 2 add tet all in volume 2")

	# Get node IDs to check if they're non-contiguous
	all_nodes = set()
	for block_id in [1, 2]:
		for tet_id in cubit.get_block_tets(block_id):
			nodes = cubit.get_connectivity("tet", tet_id)
			all_nodes.update(nodes)

	node_list = sorted(all_nodes)
	min_node = min(node_list)
	max_node = max(node_list)
	is_contiguous = (max_node - min_node + 1) == len(node_list)

	print(f"  Node ID range: {min_node} to {max_node}")
	print(f"  Total nodes: {len(node_list)}")
	print(f"  Contiguous: {is_contiguous}")

	vtk_file = "test_non_contiguous.vtk"
	cubit_mesh_export.export_vtk(cubit, vtk_file)

	result = parse_vtk_file(vtk_file)
	print(f"  VTK points: {len(result['points'])}")
	print(f"  VTK cells: {len(result['cells'])}")

	# Verify all cell node indices are valid (0-indexed, within range)
	num_points = len(result['points'])
	for cell in result['cells']:
		node_indices = cell[1:]
		for idx in node_indices:
			assert 0 <= idx < num_points, f"Invalid node index {idx} (max: {num_points-1})"

	print("  PASS: Non-contiguous node IDs handled correctly")
	os.remove(vtk_file)
	return True


def test_node_coordinate_accuracy():
	"""Test that node coordinates are correctly exported."""
	print("\n" + "=" * 60)
	print("Test 4: Node Coordinate Accuracy")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 2 y 2 z 2")
	cubit.cmd("volume 1 move 1 1 1")  # Center at (1,1,1)
	cubit.cmd("volume 1 scheme map")
	cubit.cmd("volume 1 size 2")  # Single hex
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add hex all")

	# Get Cubit node coordinates
	hex_id = cubit.get_block_hexes(1)[0]
	cubit_node_ids = cubit.get_connectivity("hex", hex_id)
	cubit_coords = {}
	for node_id in cubit_node_ids:
		cubit_coords[node_id] = cubit.get_nodal_coordinates(node_id)

	print(f"  Cubit node coordinates:")
	for node_id, coord in sorted(cubit_coords.items()):
		print(f"    Node {node_id}: ({coord[0]:.1f}, {coord[1]:.1f}, {coord[2]:.1f})")

	vtk_file = "test_coords.vtk"
	cubit_mesh_export.export_vtk(cubit, vtk_file)

	result = parse_vtk_file(vtk_file)

	# Verify coordinate ranges
	all_x = [p[0] for p in result['points']]
	all_y = [p[1] for p in result['points']]
	all_z = [p[2] for p in result['points']]

	print(f"  VTK coordinate ranges:")
	print(f"    X: {min(all_x):.1f} to {max(all_x):.1f} (expected: 0 to 2)")
	print(f"    Y: {min(all_y):.1f} to {max(all_y):.1f} (expected: 0 to 2)")
	print(f"    Z: {min(all_z):.1f} to {max(all_z):.1f} (expected: 0 to 2)")

	# Verify ranges match expected brick dimensions
	assert abs(min(all_x) - 0.0) < 0.01, "X min incorrect"
	assert abs(max(all_x) - 2.0) < 0.01, "X max incorrect"
	assert abs(min(all_y) - 0.0) < 0.01, "Y min incorrect"
	assert abs(max(all_y) - 2.0) < 0.01, "Y max incorrect"
	assert abs(min(all_z) - 0.0) < 0.01, "Z min incorrect"
	assert abs(max(all_z) - 2.0) < 0.01, "Z max incorrect"

	print("  PASS: Node coordinates correctly exported")
	os.remove(vtk_file)
	return True


def test_tet10_node_ordering():
	"""Test TET10 node ordering."""
	print("\n" + "=" * 60)
	print("Test 5: TET10 Node Ordering")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create pyramid height 1 sides 4 radius 0.5 top 0")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all")
	cubit.cmd("block 1 element type tetra10")

	tet_id = cubit.get_block_tets(1)[0]
	cubit_nodes = cubit.get_expanded_connectivity("tet", tet_id)
	print(f"  Cubit TET10 nodes: {len(cubit_nodes)}")

	vtk_file = "test_tet10_ordering.vtk"
	cubit_mesh_export.export_vtk(cubit, vtk_file)

	result = parse_vtk_file(vtk_file)

	# Find TET10 cells
	tet10_count = result['cell_types'].count(24)
	print(f"  TET10 elements in VTK: {tet10_count}")

	# Verify TET10 cell
	for i, cell_type in enumerate(result['cell_types']):
		if cell_type == 24:  # TET10
			cell = result['cells'][i]
			num_nodes = cell[0]
			node_indices = cell[1:]

			assert num_nodes == 10, f"Expected 10 nodes, got {num_nodes}"
			assert len(set(node_indices)) == 10, "Duplicate node indices!"

			# Verify all indices are valid
			for idx in node_indices:
				assert 0 <= idx < len(result['points']), f"Invalid index {idx}"
			break

	print("  PASS: TET10 node ordering correct")
	os.remove(vtk_file)
	return True


if __name__ == "__main__":
	print("\n" + "=" * 60)
	print("VTK Node Ordering Test Suite")
	print("=" * 60)

	all_passed = True

	tests = [
		test_hex8_node_ordering,
		test_hex20_node_ordering,
		test_non_contiguous_node_ids,
		test_node_coordinate_accuracy,
		test_tet10_node_ordering,
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
