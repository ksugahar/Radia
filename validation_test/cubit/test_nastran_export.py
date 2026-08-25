"""
Test Nastran export via Cubit's tool-neutral ``export nastran_bdf`` command.

Tests:
1. 3D mesh export (tet, hex)
2. 2D mesh export (tri, quad)
3. PYRAM option (pyramid to hex conversion)
4. Nastran format, property, sideset, nodeset, and alias validation
"""

import os
import sys
# Auto-detect Cubit installation (single source of truth)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))
from install_panels import find_cubit_bin
_cubit_path = find_cubit_bin()
if _cubit_path and _cubit_path not in sys.path:
	sys.path.append(_cubit_path)

_plugin_dir = os.path.join(os.path.dirname(_cubit_path), 'bin', 'plugins') \
	if _cubit_path else None
if _plugin_dir and os.path.isdir(_plugin_dir):
	os.environ['CUBIT_PLUGIN_DIR'] = _plugin_dir

import cubit

# Initialize Cubit
cubit.init(['cubit', '-nojournal', '-batch',
	'-commandplugindir', _plugin_dir or ''])


def parse_nastran_file(filename):
	"""Parse Nastran BDF file and return element statistics."""
	with open(filename, 'r') as f:
		content = f.read()

	result = {
		'grids': 0,
		'grid_coordinates': [],
		'elements': {},
		'element_pids': {},
		'properties': {},
		'property_ids': {},
		'sets': {},
		'has_header': False,
	}

	lines = content.split('\n')
	for line_index, line in enumerate(lines):
		line_upper = line.upper().strip()

		# Check for header
		if line_upper.startswith(('$', 'BEGIN')):
			result['has_header'] = True

		# Count GRID cards
		if line_upper.startswith('GRID'):
			result['grids'] += 1
			if line_upper.startswith('GRID*') and line_index + 1 < len(lines):
				result['grid_coordinates'].append((
					float(line[40:56]),
					float(line[56:72]),
					float(lines[line_index + 1][8:24]),
				))

		# Count element types
		element_types = ['CTETRA', 'CHEXA', 'CPENTA', 'CPYRAM',
		                 'CTRIA3', 'CTRIA6', 'CQUAD4', 'CQUAD8']
		for elem_type in element_types:
			if line_upper.startswith(elem_type):
				result['elements'][elem_type] = result['elements'].get(elem_type, 0) + 1
				pid = int(line[16:24])
				result['element_pids'].setdefault(elem_type, set()).add(pid)

		for property_type in ['PSOLID', 'PSHELL']:
			if line_upper.startswith(property_type):
				result['properties'][property_type] = \
					result['properties'].get(property_type, 0) + 1
				pid = int(line[8:16])
				mid = int(line[16:24])
				result['property_ids'][pid] = (property_type, mid)

		if line_upper.startswith('SET1'):
			set_id = int(line[8:16])
			result['sets'][set_id] = result['sets'].get(set_id, 0) + 1

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
	cubit.cmd(f'export nastran_bdf "{bdf_file}" dimension 3 overwrite')

	result = parse_nastran_file(bdf_file)
	print(f"  Nastran file: {result['grids']} GRIDs")
	print(f"  Elements: {result['elements']}")

	assert 'CTETRA' in result['elements'], "CTETRA not found!"
	assert result['elements']['CTETRA'] == num_tets, "CTETRA count mismatch!"
	assert result['properties'].get('PSOLID') == 1
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
	cubit.cmd(f'export nastran_bdf "{bdf_file}" dimension 3 overwrite')

	result = parse_nastran_file(bdf_file)
	print(f"  Nastran file: {result['grids']} GRIDs")
	print(f"  Elements: {result['elements']}")

	assert 'CHEXA' in result['elements'], "CHEXA not found!"
	assert result['elements']['CHEXA'] == num_hexes, "CHEXA count mismatch!"
	assert result['properties'].get('PSOLID') == 1
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
	cubit.cmd("move surface 1 z 2")
	cubit.cmd("surface 1 scheme trimesh")
	cubit.cmd("surface 1 size 0.3")
	cubit.cmd("mesh surface 1")

	cubit.cmd("block 1 add tri all")
	cubit.cmd("block 1 name 'surface'")

	num_tris = len(cubit.get_block_tris(1))
	print(f"  Cubit mesh: {num_tris} tris")

	bdf_file = "test_tri.bdf"
	cubit.cmd(f'export nastran_bdf "{bdf_file}" dimension 2 overwrite')

	result = parse_nastran_file(bdf_file)
	print(f"  Nastran file: {result['grids']} GRIDs")
	print(f"  Elements: {result['elements']}")

	assert 'CTRIA3' in result['elements'], "CTRIA3 not found!"
	assert result['elements']['CTRIA3'] == num_tris, "CTRIA3 count mismatch!"
	assert result['properties'].get('PSHELL') == 1
	assert 'PSOLID' not in result['properties']
	assert result['grid_coordinates']
	assert all(abs(coords[2] - 2.0) < 1e-12
		for coords in result['grid_coordinates'])
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
	cubit.cmd(f'export nastran_bdf "{bdf_file}" dimension 2 overwrite')

	result = parse_nastran_file(bdf_file)
	print(f"  Nastran file: {result['grids']} GRIDs")
	print(f"  Elements: {result['elements']}")

	assert 'CQUAD4' in result['elements'], "CQUAD4 not found!"
	assert result['elements']['CQUAD4'] == num_quads, "CQUAD4 count mismatch!"
	assert result['properties'].get('PSHELL') == 1
	assert 'PSOLID' not in result['properties']
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
	cubit.cmd(f'export nastran_bdf "{bdf_file}" dimension 3 overwrite')

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

	for cmd in [
		"reset",
		"create Cylinder height 1 radius 0.5",
		"create sphere radius 0.5",
		"move Volume 2 z 0.5 include_merged",
		"unite volume all",
		"compress",
		"create boundary_layer 1",
		"modify boundary_layer 1 uniform height 0.02 growth 1.2 layers 3",
		"modify boundary_layer 1 add surface 1 volume 1 surface 3 volume 1",
		"modify boundary_layer 1 continuity on",
		"volume 1 scheme tetmesh",
		"volume 1 size 0.2",
		"mesh volume 1",
		"block 1 add volume 1",
		"block 1 name 'wedge_region'",
	]:
		cubit.cmd(cmd)

	num_wedges = len(cubit.parse_cubit_list("wedge", "all"))
	assert num_wedges > 0, "Wedge validation model generated no wedges"

	print(f"  Cubit mesh: {num_wedges} wedges")

	bdf_file = "test_wedge.bdf"
	cubit.cmd(f'export nastran_bdf "{bdf_file}" dimension 3 overwrite')

	result = parse_nastran_file(bdf_file)
	print(f"  Nastran file: {result['grids']} GRIDs")
	print(f"  Elements: {result['elements']}")

	assert result['elements'].get('CPENTA') == num_wedges, \
		"CPENTA count mismatch"
	print("  PASS: Wedge export correct")

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
	cubit.cmd(f'export nastran_bdf "{bdf_file}" dimension 3 overwrite')

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

	for cmd in [
		"reset",
		"brick x 2 y 1 z 1",
		"webcut volume 1 with plane xplane imprint merge",
		"volume 1 scheme map",
		"volume 1 size 1",
		"mesh volume 1",
		"volume 2 scheme tetmesh",
		"volume 2 size 5",
		"mesh volume 2",
		"block 1 add volume 1",
		"block 1 name 'map'",
		"block 2 add volume 2",
		"block 2 name 'tet_pyramid'",
	]:
		cubit.cmd(cmd)
	num_pyramids = len(cubit.parse_cubit_list("pyramid", "all"))
	assert num_pyramids > 0, "Pyramid validation model generated no pyramids"

	# Test with PYRAM=True (default - pyramids allowed)
	bdf_file = "test_pyram_true.bdf"
	cubit.cmd(f'export nastran_bdf "{bdf_file}" dimension 3 overwrite')
	result_true = parse_nastran_file(bdf_file)
	print(f"  PYRAM=True: {result_true['elements']}")
	assert result_true['elements'].get('CPYRAM') == num_pyramids
	os.remove(bdf_file)

	# Test with PYRAM=False (nopyramid)
	bdf_file = "test_pyram_false.bdf"
	cubit.cmd(f'export nastran_bdf "{bdf_file}" dimension 3 nopyramid overwrite')
	result_false = parse_nastran_file(bdf_file)
	print(f"  PYRAM=False: {result_false['elements']}")
	assert 'CPYRAM' not in result_false['elements']
	assert result_false['elements'].get('CHEXA', 0) == \
		result_true['elements'].get('CHEXA', 0) + num_pyramids
	os.remove(bdf_file)

	print("  PASS: PYRAM option works")
	return True


def test_groups_properties_and_alias():
	"""Check machine-readable groups, collision-free PIDs, and legacy alias."""
	print("\n" + "=" * 60)
	print("Test 9: Groups, Properties, and Compatibility Alias")
	print("=" * 60)

	for cmd in [
		"reset",
		"create brick x 1 y 1 z 1",
		"volume 1 scheme tetmesh",
		"volume 1 size 0.5",
		"mesh volume 1",
		"block 1 add tet all",
		"block 1 name 'solid'",
		"sideset 1 add surface 1",
		"sideset 1 name 'boundary'",
		"nodeset 1 add node all",
		"nodeset 1 name 'all_nodes'",
	]:
		cubit.cmd(cmd)

	primary_file = "test_groups_primary.bdf"
	alias_file = "test_groups_alias.bdf"
	cubit.cmd(f'export nastran_bdf "{primary_file}" dimension 3 overwrite')
	cubit.cmd(f'export jmag_nastran "{alias_file}" dimension 3 overwrite')

	primary = parse_nastran_file(primary_file)
	alias = parse_nastran_file(alias_file)
	assert primary['elements'] == alias['elements']
	assert primary['property_ids'] == alias['property_ids']
	assert primary['sets'] == alias['sets']

	assert primary['element_pids']['CTETRA'] == {1}
	assert primary['property_ids'][1][0] == 'PSOLID'
	assert primary['element_pids']['CTRIA3'] == {2}
	assert primary['property_ids'][2] == ('PSHELL', 1)
	assert primary['sets'] == {1: 1}
	print("  PASS: block/sideset property IDs do not collide, SET1 is present")
	print("  PASS: deprecated jmag_nastran alias matches nastran_bdf")

	os.remove(primary_file)
	os.remove(alias_file)
	return True


def test_dimension_filter_rejects_volume_only_2d():
	"""A 3D block must never become a zero-thickness 2D BDF by accident."""
	print("\n" + "=" * 60)
	print("Test 10: Dimension Filter Rejects Volume-only 2D Export")
	print("=" * 60)

	for cmd in [
		"reset",
		"create brick x 1 y 1 z 1",
		"volume 1 scheme tetmesh",
		"volume 1 size 0.5",
		"mesh volume 1",
		"block 1 add tet all",
	]:
		cubit.cmd(cmd)

	bdf_file = "test_volume_as_2d.bdf"
	if os.path.exists(bdf_file):
		os.remove(bdf_file)
	cubit.cmd(f'export nastran_bdf "{bdf_file}" dimension 2 overwrite')
	assert not os.path.exists(bdf_file), \
		"Volume-only dimension=2 export silently wrote a BDF"
	print("  PASS: incompatible dimension is rejected without writing a file")
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
		test_groups_properties_and_alias,
		test_dimension_filter_rejects_volume_only_2d,
	]

	for test in tests:
		try:
			if not test():
				all_passed = False
		except Exception as e:  # noqa: BLE001 - continue the validation suite
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
	cubit.destroy()
	sys.exit(0 if all_passed else 1)
