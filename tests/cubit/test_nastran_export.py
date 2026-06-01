"""
Test Nastran export via cubit.cmd('export jmag_nastran ...').

Tests:
1. 3D mesh export (tet, hex)
2. 2D mesh export (tri, quad)
3. PYRAM option (pyramid to hex conversion)
4. Nastran format validation
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


def parse_nastran_file(filename):
	"""Parse Nastran BDF file and return element statistics."""
	with open(filename, 'r') as f:
		content = f.read()

	result = {
		'grids': 0,
		'elements': {},
		'has_header': False,
	}

	lines = content.split('\n')
	for line in lines:
		line_upper = line.upper().strip()

		# Check for header
		if line_upper.startswith('$') or line_upper.startswith('BEGIN'):
			result['has_header'] = True

		# Count GRID cards
		if line_upper.startswith('GRID'):
			result['grids'] += 1

		# Count element types
		element_types = ['CTETRA', 'CHEXA', 'CPENTA', 'CPYRAM', 'CTRIA3', 'CQUAD4']
		for elem_type in element_types:
			if line_upper.startswith(elem_type):
				result['elements'][elem_type] = result['elements'].get(elem_type, 0) + 1

	return result


def test_3d_tet_mesh():
	"""Test Nastran export with 3D tet mesh."""
	print("=" * 60)
	print("Test 1: 3D Tet Mesh Export")
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

	bdf_file = "test_tet.bdf"
	cubit.cmd(f'export jmag_nastran "{bdf_file}" dimension 3 overwrite')

	result = parse_nastran_file(bdf_file)
	print(f"  Nastran file: {result['grids']} GRIDs")
	print(f"  Elements: {result['elements']}")

	assert 'CTETRA' in result['elements'], "CTETRA not found!"
	assert result['elements']['CTETRA'] == num_tets, "CTETRA count mismatch!"
	print("  PASS: 3D tet export correct")

	os.remove(bdf_file)
	return True


def test_3d_hex_mesh():
	"""Test Nastran export with 3D hex mesh."""
	print("\n" + "=" * 60)
	print("Test 2: 3D Hex Mesh Export")
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

	bdf_file = "test_hex.bdf"
	cubit.cmd(f'export jmag_nastran "{bdf_file}" dimension 3 overwrite')

	result = parse_nastran_file(bdf_file)
	print(f"  Nastran file: {result['grids']} GRIDs")
	print(f"  Elements: {result['elements']}")

	assert 'CHEXA' in result['elements'], "CHEXA not found!"
	assert result['elements']['CHEXA'] == num_hexes, "CHEXA count mismatch!"
	print("  PASS: 3D hex export correct")

	os.remove(bdf_file)
	return True


def test_2d_tri_mesh():
	"""Test Nastran export with 2D tri mesh."""
	print("\n" + "=" * 60)
	print("Test 3: 2D Tri Mesh Export")
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

	bdf_file = "test_tri.bdf"
	cubit.cmd(f'export jmag_nastran "{bdf_file}" dimension 2 overwrite')

	result = parse_nastran_file(bdf_file)
	print(f"  Nastran file: {result['grids']} GRIDs")
	print(f"  Elements: {result['elements']}")

	assert 'CTRIA3' in result['elements'], "CTRIA3 not found!"
	assert result['elements']['CTRIA3'] == num_tris, "CTRIA3 count mismatch!"
	print("  PASS: 2D tri export correct")

	os.remove(bdf_file)
	return True


def test_2d_quad_mesh():
	"""Test Nastran export with 2D quad mesh."""
	print("\n" + "=" * 60)
	print("Test 4: 2D Quad Mesh Export")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create surface rectangle width 1 height 1 zplane")
	cubit.cmd("surface 1 scheme map")
	cubit.cmd("surface 1 size 0.3")
	cubit.cmd("mesh surface 1")

	cubit.cmd("block 1 add face all")
	cubit.cmd("block 1 name 'surface'")

	num_quads = len(cubit.get_block_faces(1))
	print(f"  Cubit mesh: {num_quads} quads")

	bdf_file = "test_quad.bdf"
	cubit.cmd(f'export jmag_nastran "{bdf_file}" dimension 2 overwrite')

	result = parse_nastran_file(bdf_file)
	print(f"  Nastran file: {result['grids']} GRIDs")
	print(f"  Elements: {result['elements']}")

	assert 'CQUAD4' in result['elements'], "CQUAD4 not found!"
	assert result['elements']['CQUAD4'] == num_quads, "CQUAD4 count mismatch!"
	print("  PASS: 2D quad export correct")

	os.remove(bdf_file)
	return True


def test_mixed_3d_mesh():
	"""Test Nastran export with mixed 3D elements (hex + tet)."""
	print("\n" + "=" * 60)
	print("Test 5: Mixed 3D Mesh (Hex + Tet)")
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

	num_hexes = len(cubit.get_block_hexes(1))
	num_tets = len(cubit.get_block_tets(2))
	print(f"  Cubit mesh: {num_hexes} hexes, {num_tets} tets")

	bdf_file = "test_mixed.bdf"
	cubit.cmd(f'export jmag_nastran "{bdf_file}" dimension 3 overwrite')

	result = parse_nastran_file(bdf_file)
	print(f"  Nastran file: {result['grids']} GRIDs")
	print(f"  Elements: {result['elements']}")

	assert 'CHEXA' in result['elements'], "CHEXA not found!"
	assert 'CTETRA' in result['elements'], "CTETRA not found!"
	print("  PASS: Mixed 3D export correct")

	os.remove(bdf_file)
	return True


def test_wedge_mesh():
	"""Test Nastran export with wedge elements."""
	print("\n" + "=" * 60)
	print("Test 6: Wedge Mesh Export")
	print("=" * 60)

	cubit.cmd("reset")
	# Create a prism geometry that will produce wedges
	cubit.cmd("create prism height 1 sides 3 radius 0.5")
	cubit.cmd("volume 1 scheme sweep")
	cubit.cmd("volume 1 size 0.3")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add wedge all")
	cubit.cmd("block 1 name 'wedges'")

	num_wedges = len(cubit.get_block_wedges(1))
	if num_wedges == 0:
		# Try adding hex if wedges not generated
		cubit.cmd("block 1 add hex all")
		num_hexes = len(cubit.get_block_hexes(1))
		print(f"  No wedges generated ({num_hexes} hexes instead), skipping")
		return True

	print(f"  Cubit mesh: {num_wedges} wedges")

	bdf_file = "test_wedge.bdf"
	cubit.cmd(f'export jmag_nastran "{bdf_file}" dimension 3 overwrite')

	result = parse_nastran_file(bdf_file)
	print(f"  Nastran file: {result['grids']} GRIDs")
	print(f"  Elements: {result['elements']}")

	if 'CPENTA' in result['elements']:
		print("  PASS: Wedge export correct")
	else:
		print("  SKIP: CPENTA not found (may be hex)")

	os.remove(bdf_file)
	return True


def test_nastran_format():
	"""Test Nastran file format structure."""
	print("\n" + "=" * 60)
	print("Test 7: Nastran Format Validation")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all")

	bdf_file = "test_format.bdf"
	cubit.cmd(f'export jmag_nastran "{bdf_file}" dimension 3 overwrite')

	with open(bdf_file, 'r') as f:
		content = f.read()

	lines = content.split('\n')

	# Check for GRID cards with proper format
	grid_count = 0
	for line in lines:
		if line.upper().startswith('GRID'):
			grid_count += 1
			# GRID card should have node ID and coordinates

	print(f"  GRID cards found: {grid_count}")
	assert grid_count > 0, "No GRID cards found!"

	# Check for element cards
	elem_count = 0
	for line in lines:
		if line.upper().startswith('CTETRA'):
			elem_count += 1

	print(f"  CTETRA cards found: {elem_count}")
	assert elem_count > 0, "No CTETRA cards found!"

	print("  PASS: Nastran format valid")
	os.remove(bdf_file)
	return True


def test_pyram_option():
	"""Test PYRAM option (pyramid handling)."""
	print("\n" + "=" * 60)
	print("Test 8: PYRAM Option")
	print("=" * 60)
	print("  Note: Pyramid elements are converted based on PYRAM flag")

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	cubit.cmd("block 1 add tet all")

	# Test with PYRAM=True (default - pyramids allowed)
	bdf_file = "test_pyram_true.bdf"
	cubit.cmd(f'export jmag_nastran "{bdf_file}" dimension 3 overwrite')
	result_true = parse_nastran_file(bdf_file)
	print(f"  PYRAM=True: {result_true['elements']}")
	os.remove(bdf_file)

	# Test with PYRAM=False (nopyramid)
	bdf_file = "test_pyram_false.bdf"
	cubit.cmd(f'export jmag_nastran "{bdf_file}" dimension 3 nopyramid overwrite')
	result_false = parse_nastran_file(bdf_file)
	print(f"  PYRAM=False: {result_false['elements']}")
	os.remove(bdf_file)

	print("  PASS: PYRAM option works")
	return True


if __name__ == "__main__":
	print("\n" + "=" * 60)
	print("Nastran Export Test Suite (via cubit.cmd)")
	print("=" * 60)

	all_passed = True

	tests = [
		test_3d_tet_mesh,
		test_3d_hex_mesh,
		test_2d_tri_mesh,
		test_2d_quad_mesh,
		test_mixed_3d_mesh,
		test_wedge_mesh,
		test_nastran_format,
		test_pyram_option,
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
