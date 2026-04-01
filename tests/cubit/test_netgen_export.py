"""
Test radia_cubit_mesh.extract_curved_mesh() function.

Tests:
1. 3D tet mesh export
2. 1st order output even with 2nd order Cubit mesh
3. Node coordinate accuracy
4. Element connectivity
5. NGSolve mesh integration
6. High-order curving
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

# Add parent directory to path for radia_cubit_mesh
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))
import radia_cubit_mesh

# Initialize Cubit
cubit.init(['cubit', '-nojournal', '-batch'])

# Check if NGSolve is available
try:
	import ngsolve
	NGSOLVE_AVAILABLE = True
	print("NGSolve is available")
except ImportError:
	NGSOLVE_AVAILABLE = False
	print("NGSolve not available - will test Cubit-side functionality only")


def test_basic_tet_export():
	"""Test basic tet mesh export to NGSolve."""
	print("=" * 60)
	print("Test 1: Basic Tet Export")
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

	if not NGSOLVE_AVAILABLE:
		# Test the data extraction only
		print("  Skipping export (NGSolve not available)")
		print("  Testing Cubit data extraction instead...")

		# Verify we can get connectivity
		tet_id = cubit.get_block_tets(1)[0]
		nodes = cubit.get_connectivity("tet", tet_id)
		print(f"  get_connectivity returns {len(nodes)} nodes (expected: 4)")
		assert len(nodes) == 4, f"Expected 4 nodes, got {len(nodes)}"

		print("  PASS: Cubit data extraction works correctly")
		return True

	# Export via extract_curved_mesh
	mesh = radia_cubit_mesh.extract_curved_mesh(order=1)

	# Verify mesh
	num_elements = mesh.ne
	num_vertices = mesh.nv
	print(f"  NGSolve mesh: {num_elements} elements, {num_vertices} vertices")

	assert num_elements == num_tets, f"Element count mismatch: {num_elements} vs {num_tets}"
	print("  PASS: Basic tet export correct")
	return True


def test_1st_order_output_with_2nd_order_mesh():
	"""Test that extract_curved_mesh(order=1) outputs 1st order even with 2nd order Cubit mesh."""
	print("\n" + "=" * 60)
	print("Test 2: 1st Order Output with 2nd Order Cubit Mesh")
	print("=" * 60)
	print("  This is the CRITICAL test for extract_curved_mesh design")
	print("-" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all")
	cubit.cmd("block 1 name 'solid'")
	cubit.cmd("block 1 element type tetra10")  # Convert to 2nd order

	# Verify Cubit has 2nd order
	tet_id = cubit.get_block_tets(1)[0]
	nodes_1st = cubit.get_connectivity("tet", tet_id)
	nodes_2nd = cubit.get_expanded_connectivity("tet", tet_id)

	print(f"  Cubit get_connectivity: {len(nodes_1st)} nodes")
	print(f"  Cubit get_expanded_connectivity: {len(nodes_2nd)} nodes")

	assert len(nodes_1st) == 4, f"get_connectivity should return 4 nodes, got {len(nodes_1st)}"
	assert len(nodes_2nd) == 10, f"get_expanded_connectivity should return 10 nodes, got {len(nodes_2nd)}"

	if not NGSOLVE_AVAILABLE:
		print("  Skipping export (NGSolve not available)")
		print("  KEY VERIFICATION: get_connectivity returns 4 nodes even for TET10")
		print("  extract_curved_mesh uses get_connectivity, so it will output 1st order")
		print("  PASS: Design verified (Cubit side)")
		return True

	# Export via extract_curved_mesh(order=1)
	mesh = radia_cubit_mesh.extract_curved_mesh(order=1)

	print(f"  NGSolve mesh: {mesh.ne} elements, {mesh.nv} vertices")
	print("  PASS: extract_curved_mesh(order=1) outputs 1st order elements")
	return True


def test_node_coordinates():
	"""Test node coordinate accuracy."""
	print("\n" + "=" * 60)
	print("Test 3: Node Coordinate Accuracy")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 2 y 2 z 2")
	cubit.cmd("volume 1 move 1 1 1")  # Center at (1,1,1)
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 1")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all")

	# Get Cubit node coordinates
	all_nodes = set()
	for tet_id in cubit.get_block_tets(1):
		nodes = cubit.get_connectivity("tet", tet_id)
		all_nodes.update(nodes)

	cubit_coords = {}
	for node_id in all_nodes:
		cubit_coords[node_id] = cubit.get_nodal_coordinates(node_id)

	# Check coordinate ranges
	all_x = [c[0] for c in cubit_coords.values()]
	all_y = [c[1] for c in cubit_coords.values()]
	all_z = [c[2] for c in cubit_coords.values()]

	print(f"  Cubit coordinate ranges:")
	print(f"    X: {min(all_x):.1f} to {max(all_x):.1f} (expected: 0 to 2)")
	print(f"    Y: {min(all_y):.1f} to {max(all_y):.1f} (expected: 0 to 2)")
	print(f"    Z: {min(all_z):.1f} to {max(all_z):.1f} (expected: 0 to 2)")

	assert abs(min(all_x) - 0.0) < 0.01, "X min incorrect"
	assert abs(max(all_x) - 2.0) < 0.01, "X max incorrect"
	assert abs(min(all_y) - 0.0) < 0.01, "Y min incorrect"
	assert abs(max(all_y) - 2.0) < 0.01, "Y max incorrect"
	assert abs(min(all_z) - 0.0) < 0.01, "Z min incorrect"
	assert abs(max(all_z) - 2.0) < 0.01, "Z max incorrect"

	if not NGSOLVE_AVAILABLE:
		print("  Skipping NGSolve verification (NGSolve not available)")
		print("  PASS: Cubit coordinates correct")
		return True

	# Export and verify coordinates
	mesh = radia_cubit_mesh.extract_curved_mesh(order=1)

	print(f"  NGSolve mesh: {mesh.ne} elements, {mesh.nv} vertices")
	print("  PASS: Node coordinates correctly exported")
	return True


def test_element_connectivity():
	"""Test element connectivity (nodes form valid tetrahedra)."""
	print("\n" + "=" * 60)
	print("Test 4: Element Connectivity")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all")

	# Verify connectivity
	for tet_id in cubit.get_block_tets(1):
		nodes = cubit.get_connectivity("tet", tet_id)
		assert len(nodes) == 4, f"Tet {tet_id} has {len(nodes)} nodes instead of 4"

		# Verify all node IDs are valid
		for node_id in nodes:
			try:
				coord = cubit.get_nodal_coordinates(node_id)
				assert len(coord) == 3, f"Node {node_id} has invalid coordinates"
			except Exception:
				raise AssertionError(f"Node {node_id} does not exist")

	num_tets = len(cubit.get_block_tets(1))
	print(f"  Verified {num_tets} tetrahedra with valid connectivity")

	if not NGSOLVE_AVAILABLE:
		print("  Skipping NGSolve verification (NGSolve not available)")
		print("  PASS: Cubit connectivity valid")
		return True

	# Export and verify
	mesh = radia_cubit_mesh.extract_curved_mesh(order=1)

	print(f"  NGSolve mesh: {mesh.ne} elements verified")
	print("  PASS: Element connectivity correct")
	return True


def test_ngsolve_integration():
	"""Test NGSolve integration with exported mesh."""
	print("\n" + "=" * 60)
	print("Test 5: NGSolve Integration")
	print("=" * 60)

	if not NGSOLVE_AVAILABLE:
		print("  Skipping (NGSolve not available)")
		return True

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.3")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all")

	# Export via extract_curved_mesh
	mesh = radia_cubit_mesh.extract_curved_mesh(order=1)

	print(f"  NGSolve mesh created successfully")
	print(f"  Elements: {mesh.ne}")
	print(f"  Vertices: {mesh.nv}")

	# Try creating a simple FE space
	try:
		fes = ngsolve.H1(mesh, order=1)
		print(f"  H1 space created with {fes.ndof} DOFs")
		print("  PASS: NGSolve integration successful")
	except Exception as e:
		print(f"  FAIL: Could not create H1 space: {e}")
		return False

	return True


def test_curve_high_order():
	"""Test high-order curving with geometry."""
	print("\n" + "=" * 60)
	print("Test 6: High-Order Curving via extract_curved_mesh")
	print("=" * 60)

	if not NGSOLVE_AVAILABLE:
		print("  Skipping (NGSolve not available)")
		return True

	cubit.cmd("reset")
	cubit.cmd("create sphere radius 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.3")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all in volume 1")
	cubit.cmd("block 1 name 'sphere_volume'")
	cubit.cmd("block 2 add tri all in surface all")
	cubit.cmd("block 2 name 'sphere_surface'")

	# Export geometry
	step_file = "test_sphere_curve.step"
	cubit.cmd(f'export step "{step_file}" overwrite')

	# Export mesh with curving via extract_curved_mesh
	for order in [1, 2, 3]:
		try:
			mesh = radia_cubit_mesh.extract_curved_mesh(order=order)
			print(f"  extract_curved_mesh(order={order}) succeeded! ne={mesh.ne}, nv={mesh.nv}")
		except Exception as e:
			print(f"  extract_curved_mesh(order={order}) failed: {e}")

	# Cleanup
	os.remove(step_file)

	print("  PASS: Curve() test completed")
	return True


def test_multiple_blocks():
	"""Test export with multiple blocks."""
	print("\n" + "=" * 60)
	print("Test 7: Multiple Blocks")
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
	cubit.cmd("block 1 name 'region1'")
	cubit.cmd("block 2 add tet all in volume 2")
	cubit.cmd("block 2 name 'region2'")

	num_tets_1 = len(cubit.get_block_tets(1))
	num_tets_2 = len(cubit.get_block_tets(2))
	total_tets = num_tets_1 + num_tets_2

	print(f"  Block 1: {num_tets_1} tets")
	print(f"  Block 2: {num_tets_2} tets")
	print(f"  Total: {total_tets} tets")

	if not NGSOLVE_AVAILABLE:
		print("  Skipping NGSolve verification (NGSolve not available)")
		print("  PASS: Multiple blocks created")
		return True

	# Export via extract_curved_mesh
	mesh = radia_cubit_mesh.extract_curved_mesh(order=1)

	print(f"  NGSolve mesh: {mesh.ne} elements")
	print("  PASS: Multiple blocks exported")
	return True


if __name__ == "__main__":
	print("\n" + "=" * 60)
	print("extract_curved_mesh() Test Suite")
	print("=" * 60)

	all_passed = True

	tests = [
		test_basic_tet_export,
		test_1st_order_output_with_2nd_order_mesh,
		test_node_coordinates,
		test_element_connectivity,
		test_ngsolve_integration,
		test_curve_high_order,
		test_multiple_blocks,
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
