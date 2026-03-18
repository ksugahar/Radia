"""
Test export_vtu() function (VTK XML format).

Tests:
1. 1st order tet mesh export
2. 2nd order tet mesh export
3. Hex mesh export (1st and 2nd order)
4. Mixed element types
5. VTU format validation
6. Non-contiguous node IDs
7. Block ID cell data
"""

import sys
import os
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
	sys.path.append(cubit_path)

import cubit
import xml.etree.ElementTree as ET

# Add parent directory to path for cubit_mesh_export
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))
import cubit_mesh_export

# Initialize Cubit
cubit.init(['cubit', '-nojournal', '-batch'])

# VTK cell type codes for reference
VTK_TYPES = {
	'VERTEX': 1,
	'LINE': 3,
	'QUADRATIC_EDGE': 21,
	'TRIANGLE': 5,
	'QUADRATIC_TRIANGLE': 22,
	'QUAD': 9,
	'QUADRATIC_QUAD': 23,
	'TETRA': 10,
	'QUADRATIC_TETRA': 24,
	'HEXAHEDRON': 12,
	'QUADRATIC_HEXAHEDRON': 25,
	'WEDGE': 13,
	'QUADRATIC_WEDGE': 26,
	'PYRAMID': 14,
	'QUADRATIC_PYRAMID': 27,
}


def parse_vtu_file(filename):
	"""Parse VTU file and return data."""
	tree = ET.parse(filename)
	root = tree.getroot()

	result = {
		'num_points': 0,
		'num_cells': 0,
		'cell_types': [],
		'has_block_ids': False,
		'has_node_ids': False,
		'valid_xml': True,
	}

	# Find UnstructuredGrid/Piece
	piece = root.find('.//Piece')
	if piece is not None:
		result['num_points'] = int(piece.get('NumberOfPoints', 0))
		result['num_cells'] = int(piece.get('NumberOfCells', 0))

	# Parse cell types
	types_array = root.find('.//Cells/DataArray[@Name="types"]')
	if types_array is not None and types_array.text:
		result['cell_types'] = [int(t) for t in types_array.text.split()]

	# Check for BlockID cell data
	block_id_array = root.find('.//CellData/DataArray[@Name="BlockID"]')
	if block_id_array is not None:
		result['has_block_ids'] = True

	# Check for NodeID point data
	node_id_array = root.find('.//PointData/DataArray[@Name="NodeID"]')
	if node_id_array is not None:
		result['has_node_ids'] = True

	return result


def test_1st_order_tet_mesh():
	"""Test VTU export with 1st order tet mesh."""
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

	num_tets = len(cubit.get_block_tets(1))
	print(f"  Cubit mesh: {num_tets} tets")

	vtu_file = "test_1st_order_tet.vtu"
	cubit_mesh_export.export_vtu(cubit, vtu_file)

	result = parse_vtu_file(vtu_file)
	print(f"  VTU file: {result['num_points']} points, {result['num_cells']} cells")
	print(f"  Cell types: {set(result['cell_types'])}")

	# Verify tet type (10 = VTK_TETRA)
	assert VTK_TYPES['TETRA'] in result['cell_types'], "VTK_TETRA not found!"
	assert result['num_cells'] == num_tets, f"Expected {num_tets} cells, got {result['num_cells']}"

	print("  PASS: 1st order tet export correct")
	os.remove(vtu_file)
	return True


def test_2nd_order_tet_mesh():
	"""Test VTU export with 2nd order tet mesh."""
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

	num_tets = len(cubit.get_block_tets(1))
	print(f"  Cubit mesh: {num_tets} tets (tetra10)")

	# Verify 2nd order in Cubit
	tet_id = cubit.get_block_tets(1)[0]
	nodes = cubit.get_expanded_connectivity("tet", tet_id)
	print(f"  Nodes per tet: {len(nodes)} (expected: 10)")

	vtu_file = "test_2nd_order_tet.vtu"
	cubit_mesh_export.export_vtu(cubit, vtu_file)

	result = parse_vtu_file(vtu_file)
	print(f"  VTU file: {result['num_points']} points, {result['num_cells']} cells")
	print(f"  Cell types: {set(result['cell_types'])}")

	# Verify quadratic tet type (24 = VTK_QUADRATIC_TETRA)
	assert VTK_TYPES['QUADRATIC_TETRA'] in result['cell_types'], "VTK_QUADRATIC_TETRA not found!"

	print("  PASS: 2nd order tet export correct")
	os.remove(vtu_file)
	return True


def test_hex_mesh():
	"""Test VTU export with hex mesh (1st and 2nd order)."""
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

	num_hexes = len(cubit.get_block_hexes(1))
	print(f"  Cubit mesh: {num_hexes} hexes")

	# Test 1st order
	vtu_file = "test_hex_1st.vtu"
	cubit_mesh_export.export_vtu(cubit, vtu_file)
	result = parse_vtu_file(vtu_file)
	print(f"  1st order - Cell types: {set(result['cell_types'])}")
	assert VTK_TYPES['HEXAHEDRON'] in result['cell_types'], "VTK_HEXAHEDRON not found!"
	print("  1st order hex: PASS")
	os.remove(vtu_file)

	# Test 2nd order
	cubit.cmd("block 1 element type hex20")
	vtu_file = "test_hex_2nd.vtu"
	cubit_mesh_export.export_vtu(cubit, vtu_file)
	result = parse_vtu_file(vtu_file)
	print(f"  2nd order - Cell types: {set(result['cell_types'])}")
	assert VTK_TYPES['QUADRATIC_HEXAHEDRON'] in result['cell_types'], "VTK_QUADRATIC_HEXAHEDRON not found!"
	print("  2nd order hex: PASS")
	os.remove(vtu_file)

	return True


def test_mixed_elements():
	"""Test VTU export with mixed element types."""
	print("\n" + "=" * 60)
	print("Test 4: Mixed Element Types")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 move 0 0 0")
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

	num_hexes = len(cubit.get_block_hexes(1))
	num_tets = len(cubit.get_block_tets(2))
	print(f"  Cubit mesh: {num_hexes} hexes, {num_tets} tets")

	vtu_file = "test_mixed.vtu"
	cubit_mesh_export.export_vtu(cubit, vtu_file)

	result = parse_vtu_file(vtu_file)
	print(f"  VTU file: {result['num_points']} points, {result['num_cells']} cells")
	print(f"  Cell types: {set(result['cell_types'])}")

	assert VTK_TYPES['HEXAHEDRON'] in result['cell_types'], "VTK_HEXAHEDRON not found!"
	assert VTK_TYPES['TETRA'] in result['cell_types'], "VTK_TETRA not found!"
	assert result['num_cells'] == num_hexes + num_tets

	print("  PASS: Mixed elements export correct")
	os.remove(vtu_file)
	return True


def test_vtu_format_validation():
	"""Test VTU file format validation."""
	print("\n" + "=" * 60)
	print("Test 5: VTU Format Validation")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all")

	vtu_file = "test_format.vtu"
	cubit_mesh_export.export_vtu(cubit, vtu_file)

	# Read file content
	with open(vtu_file, 'r') as f:
		content = f.read()

	# Check XML structure
	assert '<?xml version="1.0"?>' in content, "XML declaration missing!"
	assert '<VTKFile type="UnstructuredGrid"' in content, "VTKFile header missing!"
	assert '<UnstructuredGrid>' in content, "UnstructuredGrid missing!"
	assert '<Piece' in content, "Piece missing!"
	assert '<Points>' in content, "Points missing!"
	assert '<Cells>' in content, "Cells missing!"
	assert 'connectivity' in content, "connectivity array missing!"
	assert 'offsets' in content, "offsets array missing!"
	assert 'types' in content, "types array missing!"

	print("  Required elements present: PASS")

	# Validate XML parsing
	result = parse_vtu_file(vtu_file)
	assert result['valid_xml'], "Invalid XML!"
	print("  Valid XML: PASS")

	# Check data arrays
	assert result['has_block_ids'], "BlockID cell data missing!"
	assert result['has_node_ids'], "NodeID point data missing!"
	print("  Cell/Point data present: PASS")

	os.remove(vtu_file)
	return True


def test_non_contiguous_node_ids():
	"""Test VTU export with non-contiguous node IDs."""
	print("\n" + "=" * 60)
	print("Test 6: Non-Contiguous Node IDs")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 move 0 0 0")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 2 move 2 0 0")

	cubit.cmd("volume all scheme tetmesh")
	cubit.cmd("volume all size 0.5")
	cubit.cmd("mesh volume all")

	# Delete and remesh to create gaps in node IDs
	cubit.cmd("delete mesh volume 1 propagate")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all in volume 1")
	cubit.cmd("block 2 add tet all in volume 2")

	# Check if node IDs are non-contiguous
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

	vtu_file = "test_non_contiguous.vtu"
	cubit_mesh_export.export_vtu(cubit, vtu_file)

	result = parse_vtu_file(vtu_file)
	print(f"  VTU file: {result['num_points']} points, {result['num_cells']} cells")

	# VTU should have correct number of points
	assert result['num_points'] == len(node_list), f"Expected {len(node_list)} points, got {result['num_points']}"

	print("  PASS: Non-contiguous node IDs handled correctly")
	os.remove(vtu_file)
	return True


def test_block_id_cell_data():
	"""Test BlockID cell data in VTU export."""
	print("\n" + "=" * 60)
	print("Test 7: BlockID Cell Data")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 move 0 0 0")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 2 move 1.5 0 0")

	cubit.cmd("volume all scheme tetmesh")
	cubit.cmd("volume all size 0.5")
	cubit.cmd("mesh volume all")

	cubit.cmd("block 1 add tet all in volume 1")
	cubit.cmd("block 1 name 'material1'")
	cubit.cmd("block 2 add tet all in volume 2")
	cubit.cmd("block 2 name 'material2'")

	num_tets_1 = len(cubit.get_block_tets(1))
	num_tets_2 = len(cubit.get_block_tets(2))
	print(f"  Block 1: {num_tets_1} tets")
	print(f"  Block 2: {num_tets_2} tets")

	vtu_file = "test_block_ids.vtu"
	cubit_mesh_export.export_vtu(cubit, vtu_file)

	# Parse and check BlockID array
	tree = ET.parse(vtu_file)
	root = tree.getroot()
	block_id_array = root.find('.//CellData/DataArray[@Name="BlockID"]')

	assert block_id_array is not None, "BlockID array not found!"

	block_ids = [int(b) for b in block_id_array.text.split()]
	print(f"  BlockID values: {set(block_ids)}")

	# Verify block IDs are present
	assert 1 in block_ids, "Block ID 1 not found!"
	assert 2 in block_ids, "Block ID 2 not found!"
	assert block_ids.count(1) + block_ids.count(2) == num_tets_1 + num_tets_2

	print("  PASS: BlockID cell data correct")
	os.remove(vtu_file)
	return True


def test_2d_surface_mesh():
	"""Test VTU export with 2D surface mesh."""
	print("\n" + "=" * 60)
	print("Test 8: 2D Surface Mesh")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create surface rectangle width 1 height 1 zplane")
	cubit.cmd("surface 1 scheme trimesh")
	cubit.cmd("surface 1 size 0.3")
	cubit.cmd("mesh surface 1")

	cubit.cmd("block 1 add tri all")
	cubit.cmd("block 1 name 'surface'")

	num_tris = len(cubit.get_block_tris(1))
	print(f"  Cubit mesh: {num_tris} tris")

	vtu_file = "test_2d_surface.vtu"
	cubit_mesh_export.export_vtu(cubit, vtu_file)

	result = parse_vtu_file(vtu_file)
	print(f"  VTU file: {result['num_points']} points, {result['num_cells']} cells")
	print(f"  Cell types: {set(result['cell_types'])}")

	assert VTK_TYPES['TRIANGLE'] in result['cell_types'], "VTK_TRIANGLE not found!"
	assert result['num_cells'] == num_tris

	print("  PASS: 2D surface mesh export correct")
	os.remove(vtu_file)
	return True


if __name__ == "__main__":
	print("\n" + "=" * 60)
	print("export_vtu() Test Suite")
	print("=" * 60)

	all_passed = True

	tests = [
		test_1st_order_tet_mesh,
		test_2nd_order_tet_mesh,
		test_hex_mesh,
		test_mixed_elements,
		test_vtu_format_validation,
		test_non_contiguous_node_ids,
		test_block_id_cell_data,
		test_2d_surface_mesh,
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
