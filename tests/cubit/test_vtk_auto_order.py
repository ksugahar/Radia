"""
Test VTK export with automatic element order detection.

Tests:
1. 1st order mesh only
2. 2nd order mesh only
3. Mixed 1st and 2nd order mesh (if possible via Cubit API)
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

def test_1st_order_mesh():
	"""Test VTK export with 1st order elements."""
	print("=" * 60)
	print("Test 1: 1st Order Mesh")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	# Add elements to blocks
	cubit.cmd("block 1 add tet all")
	cubit.cmd("block 1 name 'tet_region'")
	cubit.cmd("block 2 add tri all in surface all")
	cubit.cmd("block 2 name 'boundary'")

	# Get element counts
	num_tets = len(cubit.get_block_tets(1))
	num_tris = len(cubit.get_block_tris(2))
	print(f"  Tets: {num_tets}, Tris: {num_tris}")

	# Check node count of first tet
	tet_id = cubit.get_block_tets(1)[0]
	nodes = cubit.get_expanded_connectivity("tet", tet_id)
	print(f"  Nodes per tet: {len(nodes)} (expected: 4)")

	# Export to VTK
	vtk_file = "test_1st_order.vtk"
	cubit_mesh_export.export_vtk(cubit, vtk_file)

	# Verify file and check cell types
	with open(vtk_file, 'r') as f:
		content = f.read()

	# Check for 1st order cell types (TET=10, TRI=5)
	lines = content.split('\n')
	cell_types_idx = None
	for i, line in enumerate(lines):
		if line.startswith('CELL_TYPES'):
			cell_types_idx = i
			break

	if cell_types_idx:
		# Read cell types
		cell_types = []
		for line in lines[cell_types_idx+1:]:
			if line.strip() and not line.startswith('CELL_DATA'):
				try:
					cell_types.append(int(line.strip()))
				except:
					break
			else:
				break

		print(f"  Cell types found: {set(cell_types)}")
		print(f"  Expected: 10 (TET4), 5 (TRI3)")

		# Verify
		assert 10 in cell_types, "TET4 (type 10) not found!"
		assert 5 in cell_types, "TRI3 (type 5) not found!"
		assert 24 not in cell_types, "TET10 (type 24) should not be present!"
		print("  PASS: 1st order cell types correctly detected")

	os.remove(vtk_file)
	return True


def test_2nd_order_mesh():
	"""Test VTK export with 2nd order elements."""
	print("\n" + "=" * 60)
	print("Test 2: 2nd Order Mesh")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	# Add elements to blocks first, then convert to 2nd order
	cubit.cmd("block 1 add tet all")
	cubit.cmd("block 1 name 'tet_region'")
	cubit.cmd("block 1 element type tetra10")  # Convert to 2nd order

	cubit.cmd("block 2 add tri all in surface all")
	cubit.cmd("block 2 name 'boundary'")
	cubit.cmd("block 2 element type tri6")  # Convert to 2nd order

	# Get element counts
	num_tets = len(cubit.get_block_tets(1))
	num_tris = len(cubit.get_block_tris(2))
	print(f"  Tets: {num_tets}, Tris: {num_tris}")

	# Check node count of first tet
	tet_id = cubit.get_block_tets(1)[0]
	nodes = cubit.get_expanded_connectivity("tet", tet_id)
	print(f"  Nodes per tet: {len(nodes)} (expected: 10)")

	# Export to VTK
	vtk_file = "test_2nd_order.vtk"
	cubit_mesh_export.export_vtk(cubit, vtk_file)

	# Verify file and check cell types
	with open(vtk_file, 'r') as f:
		content = f.read()

	# Check for 2nd order cell types (TET=24, TRI=22)
	lines = content.split('\n')
	cell_types_idx = None
	for i, line in enumerate(lines):
		if line.startswith('CELL_TYPES'):
			cell_types_idx = i
			break

	if cell_types_idx:
		# Read cell types
		cell_types = []
		for line in lines[cell_types_idx+1:]:
			if line.strip() and not line.startswith('CELL_DATA'):
				try:
					cell_types.append(int(line.strip()))
				except:
					break
			else:
				break

		print(f"  Cell types found: {set(cell_types)}")
		print(f"  Expected: 24 (TET10), 22 (TRI6)")

		# Verify
		assert 24 in cell_types, "TET10 (type 24) not found!"
		assert 22 in cell_types, "TRI6 (type 22) not found!"
		assert 10 not in cell_types, "TET4 (type 10) should not be present!"
		print("  PASS: 2nd order cell types correctly detected")

	os.remove(vtk_file)
	return True


def test_hex_elements():
	"""Test VTK export with hex elements (1st and 2nd order)."""
	print("\n" + "=" * 60)
	print("Test 3: Hex Elements (1st and 2nd order)")
	print("=" * 60)

	# 1st order hex
	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme map")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add hex all")
	cubit.cmd("block 1 name 'hex_region'")

	hex_id = cubit.get_block_hexes(1)[0]
	nodes_1st = cubit.get_expanded_connectivity("hex", hex_id)
	print(f"  1st order hex nodes: {len(nodes_1st)} (expected: 8)")

	vtk_file = "test_hex_1st.vtk"
	cubit_mesh_export.export_vtk(cubit, vtk_file)

	with open(vtk_file, 'r') as f:
		content = f.read()
	assert '12\n' in content, "HEX8 (type 12) not found!"
	print("  1st order hex: PASS")
	os.remove(vtk_file)

	# 2nd order hex - use block element type
	cubit.cmd("block 1 element type hex20")

	hex_id = cubit.get_block_hexes(1)[0]
	nodes_2nd = cubit.get_expanded_connectivity("hex", hex_id)
	print(f"  2nd order hex nodes: {len(nodes_2nd)} (expected: 20)")

	vtk_file = "test_hex_2nd.vtk"
	cubit_mesh_export.export_vtk(cubit, vtk_file)

	with open(vtk_file, 'r') as f:
		content = f.read()
	assert '25\n' in content, "HEX20 (type 25) not found!"
	print("  2nd order hex: PASS")
	os.remove(vtk_file)

	return True


def test_wedge_elements():
	"""Test VTK export with wedge elements."""
	print("\n" + "=" * 60)
	print("Test 4: Wedge Elements")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create cylinder radius 0.5 height 1")
	cubit.cmd("volume 1 scheme sweep")
	cubit.cmd("volume 1 size 0.3")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add wedge all")
	cubit.cmd("block 1 name 'wedge_region'")

	num_wedges = len(cubit.get_block_wedges(1))
	if num_wedges > 0:
		wedge_id = cubit.get_block_wedges(1)[0]
		nodes = cubit.get_expanded_connectivity("wedge", wedge_id)
		print(f"  Wedges: {num_wedges}, Nodes per wedge: {len(nodes)}")

		vtk_file = "test_wedge.vtk"
		cubit_mesh_export.export_vtk(cubit, vtk_file)

		with open(vtk_file, 'r') as f:
			content = f.read()

		if len(nodes) == 6:
			assert '13\n' in content, "WEDGE6 (type 13) not found!"
			print("  1st order wedge: PASS")
		elif len(nodes) == 15:
			assert '26\n' in content, "WEDGE15 (type 26) not found!"
			print("  2nd order wedge: PASS")

		os.remove(vtk_file)
	else:
		print("  No wedge elements generated, skipping")

	return True


def test_cells_size_calculation():
	"""Verify CELLS size is correctly calculated for mixed node counts."""
	print("\n" + "=" * 60)
	print("Test 5: CELLS Size Calculation")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all")

	num_tets = len(cubit.get_block_tets(1))

	# 1st order: each tet has 4 nodes, CELLS line format: "4 n1 n2 n3 n4"
	# Size = num_tets * (1 + 4) = num_tets * 5
	expected_size_1st = num_tets * 5

	vtk_file = "test_cells_size.vtk"
	cubit_mesh_export.export_vtk(cubit, vtk_file)

	with open(vtk_file, 'r') as f:
		for line in f:
			if line.startswith('CELLS'):
				parts = line.split()
				actual_num_cells = int(parts[1])
				actual_size = int(parts[2])
				print(f"  CELLS line: {line.strip()}")
				print(f"  Expected: CELLS {num_tets} {expected_size_1st}")
				assert actual_num_cells == num_tets, f"Cell count mismatch: {actual_num_cells} vs {num_tets}"
				assert actual_size == expected_size_1st, f"Size mismatch: {actual_size} vs {expected_size_1st}"
				print("  PASS: CELLS size correctly calculated")
				break

	os.remove(vtk_file)

	# Now test 2nd order - use block element type
	cubit.cmd("block 1 element type tetra10")

	# 2nd order: each tet has 10 nodes, CELLS line format: "10 n1 n2 ... n10"
	# Size = num_tets * (1 + 10) = num_tets * 11
	expected_size_2nd = num_tets * 11

	cubit_mesh_export.export_vtk(cubit, vtk_file)

	with open(vtk_file, 'r') as f:
		for line in f:
			if line.startswith('CELLS'):
				parts = line.split()
				actual_size = int(parts[2])
				print(f"  CELLS line (2nd order): {line.strip()}")
				print(f"  Expected: CELLS {num_tets} {expected_size_2nd}")
				assert actual_size == expected_size_2nd, f"Size mismatch: {actual_size} vs {expected_size_2nd}"
				print("  PASS: 2nd order CELLS size correctly calculated")
				break

	os.remove(vtk_file)
	return True


if __name__ == "__main__":
	print("\n" + "=" * 60)
	print("VTK Auto-Order Detection Test Suite")
	print("=" * 60)

	all_passed = True

	try:
		test_1st_order_mesh()
	except Exception as e:
		print(f"  FAIL: {e}")
		all_passed = False

	try:
		test_2nd_order_mesh()
	except Exception as e:
		print(f"  FAIL: {e}")
		all_passed = False

	try:
		test_hex_elements()
	except Exception as e:
		print(f"  FAIL: {e}")
		all_passed = False

	try:
		test_wedge_elements()
	except Exception as e:
		print(f"  FAIL: {e}")
		all_passed = False

	try:
		test_cells_size_calculation()
	except Exception as e:
		print(f"  FAIL: {e}")
		all_passed = False

	print("\n" + "=" * 60)
	if all_passed:
		print("All tests PASSED!")
	else:
		print("Some tests FAILED!")
	print("=" * 60)
