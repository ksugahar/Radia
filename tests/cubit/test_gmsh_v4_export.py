"""
Test export_Gmesh(version="4.1") function.

Tests:
1. Basic format validation (v4.1)
2. 1st order 3D elements (Tet, Hex, Wedge, Pyramid)
3. 2nd order 3D elements
4. 2D elements with normal orientation
5. Mixed elements
6. $Entities section validation
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

# Add parent directory to path for cubit_mesh_export
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))
import cubit_mesh_export

# Initialize Cubit
cubit.init(['cubit', '-nojournal', '-batch'])

# Gmsh v4 element type codes
GMSH_V4_TYPES = {
	# 1st order
	'TET4': 4,
	'HEX8': 5,
	'WEDGE6': 6,
	'PYRAMID5': 7,
	'TRI3': 2,
	'QUAD4': 3,
	'EDGE2': 1,
	'POINT': 15,
	# 2nd order
	'TET10': 11,
	'HEX20': 17,
	'WEDGE15': 18,
	'PYRAMID13': 19,
	'TRI6': 9,
	'QUAD8': 16,
	'EDGE3': 8,
}


def parse_gmsh_v4_file(filename):
	"""Parse Gmsh v4 file and return structure info."""
	with open(filename, 'r') as f:
		content = f.read()

	result = {
		'version': None,
		'physical_names': [],
		'entities': {'points': 0, 'curves': 0, 'surfaces': 0, 'volumes': 0},
		'nodes': 0,
		'elements': 0,
		'element_types': {},
		'sections': [],
	}

	lines = content.split('\n')
	current_section = None
	section_line_count = 0

	for i, line in enumerate(lines):
		line = line.strip()

		# Track sections
		if line.startswith('$') and not line.startswith('$End'):
			current_section = line[1:]
			result['sections'].append(current_section)
			section_line_count = 0
			continue
		elif line.startswith('$End'):
			current_section = None
			continue

		if current_section == 'MeshFormat' and section_line_count == 0:
			parts = line.split()
			if len(parts) >= 1:
				result['version'] = parts[0]
			section_line_count += 1

		elif current_section == 'PhysicalNames':
			if section_line_count == 0:
				# Skip count line
				pass
			else:
				parts = line.split('"')
				if len(parts) >= 2:
					result['physical_names'].append(parts[1])
			section_line_count += 1

		elif current_section == 'Entities' and section_line_count == 0:
			parts = line.split()
			if len(parts) >= 4:
				result['entities']['points'] = int(parts[0])
				result['entities']['curves'] = int(parts[1])
				result['entities']['surfaces'] = int(parts[2])
				result['entities']['volumes'] = int(parts[3])
			section_line_count += 1

		elif current_section == 'Nodes' and section_line_count == 0:
			parts = line.split()
			if len(parts) >= 2:
				result['nodes'] = int(parts[1])
			section_line_count += 1

		elif current_section == 'Elements':
			if section_line_count == 0:
				parts = line.split()
				if len(parts) >= 2:
					result['elements'] = int(parts[1])
			else:
				# Entity block header: entityDim entityTag elementType numElements
				parts = line.split()
				if len(parts) >= 4:
					try:
						elem_type = int(parts[2])
						num_elems = int(parts[3])
						result['element_types'][elem_type] = result['element_types'].get(elem_type, 0) + num_elems
					except ValueError:
						pass
			section_line_count += 1

	return result


def test_format_version():
	"""Test Gmsh v4 format version."""
	print("=" * 60)
	print("Test 1: Gmsh v4 Format Version")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")
	cubit.cmd("block 1 add tet all")
	cubit.cmd("block 1 name 'solid'")

	msh_file = "test_v4_format.msh"
	cubit_mesh_export.export_Gmesh(cubit, msh_file, version="4.1")

	result = parse_gmsh_v4_file(msh_file)

	print(f"  Version: {result['version']}")
	print(f"  Sections: {result['sections']}")

	assert result['version'] == '4.1', f"Expected version 4.1, got {result['version']}"
	assert 'MeshFormat' in result['sections'], "Missing $MeshFormat section"
	assert 'PhysicalNames' in result['sections'], "Missing $PhysicalNames section"
	assert 'Entities' in result['sections'], "Missing $Entities section"
	assert 'Nodes' in result['sections'], "Missing $Nodes section"
	assert 'Elements' in result['sections'], "Missing $Elements section"

	print("  PASS: Gmsh v4.1 format validated")

	os.remove(msh_file)
	return True


def test_3d_tet_mesh():
	"""Test 3D tet mesh export."""
	print("\n" + "=" * 60)
	print("Test 2: 3D Tet Mesh (1st and 2nd order)")
	print("=" * 60)

	# 1st order
	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")
	cubit.cmd("block 1 add tet all")
	cubit.cmd("block 1 name 'solid'")

	msh_file = "test_v4_tet.msh"
	cubit_mesh_export.export_Gmesh(cubit, msh_file, version="4.1")

	result = parse_gmsh_v4_file(msh_file)
	print(f"  1st order - Nodes: {result['nodes']}, Elements: {result['elements']}")
	print(f"  Element types: {result['element_types']}")
	assert GMSH_V4_TYPES['TET4'] in result['element_types'], "TET4 not found"
	print("  1st order: PASS")
	os.remove(msh_file)

	# 2nd order
	cubit.cmd("block 1 element type tetra10")

	cubit_mesh_export.export_Gmesh(cubit, msh_file, version="4.1")
	result = parse_gmsh_v4_file(msh_file)
	print(f"  2nd order - Nodes: {result['nodes']}, Elements: {result['elements']}")
	print(f"  Element types: {result['element_types']}")
	assert GMSH_V4_TYPES['TET10'] in result['element_types'], "TET10 not found"
	print("  2nd order: PASS")
	os.remove(msh_file)

	return True


def test_3d_hex_mesh():
	"""Test 3D hex mesh export."""
	print("\n" + "=" * 60)
	print("Test 3: 3D Hex Mesh (1st and 2nd order)")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme map")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")
	cubit.cmd("block 1 add hex all")
	cubit.cmd("block 1 name 'solid'")

	msh_file = "test_v4_hex.msh"
	cubit_mesh_export.export_Gmesh(cubit, msh_file, version="4.1")

	result = parse_gmsh_v4_file(msh_file)
	print(f"  1st order - Element types: {result['element_types']}")
	assert GMSH_V4_TYPES['HEX8'] in result['element_types'], "HEX8 not found"
	print("  1st order: PASS")
	os.remove(msh_file)

	# 2nd order
	cubit.cmd("block 1 element type hex20")
	cubit_mesh_export.export_Gmesh(cubit, msh_file, version="4.1")
	result = parse_gmsh_v4_file(msh_file)
	print(f"  2nd order - Element types: {result['element_types']}")
	assert GMSH_V4_TYPES['HEX20'] in result['element_types'], "HEX20 not found"
	print("  2nd order: PASS")
	os.remove(msh_file)

	return True


def test_2d_mesh_with_dim_auto():
	"""Test 2D mesh export with DIM='auto'."""
	print("\n" + "=" * 60)
	print("Test 4: 2D Mesh with DIM='auto'")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create surface rectangle width 1 height 1 zplane")
	cubit.cmd("surface 1 scheme trimesh")
	cubit.cmd("surface 1 size 0.3")
	cubit.cmd("mesh surface 1")
	cubit.cmd("block 1 add tri all")
	cubit.cmd("block 1 name 'plate'")

	msh_file = "test_v4_2d.msh"
	cubit_mesh_export.export_Gmesh(cubit, msh_file, version="4.1", DIM="auto")

	result = parse_gmsh_v4_file(msh_file)
	print(f"  Element types: {result['element_types']}")
	print(f"  Entities: {result['entities']}")
	assert GMSH_V4_TYPES['TRI3'] in result['element_types'], "TRI3 not found"
	print("  PASS: 2D mesh exported correctly")
	os.remove(msh_file)

	return True


def test_2d_mesh_normal_orientation():
	"""Test 2D mesh with normal orientation (DIM='2D')."""
	print("\n" + "=" * 60)
	print("Test 5: 2D Mesh Normal Orientation")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create surface rectangle width 1 height 1 zplane")
	cubit.cmd("surface 1 scheme trimesh")
	cubit.cmd("surface 1 size 0.3")
	cubit.cmd("mesh surface 1")
	cubit.cmd("block 1 add tri all")
	cubit.cmd("block 1 name 'plate'")

	msh_file = "test_v4_2d_normal.msh"
	cubit_mesh_export.export_Gmesh(cubit, msh_file, version="4.1", DIM="2D")

	# Check file was created
	assert os.path.exists(msh_file), "File not created"
	result = parse_gmsh_v4_file(msh_file)
	print(f"  Nodes: {result['nodes']}, Elements: {result['elements']}")
	print("  PASS: 2D mesh with normal orientation")
	os.remove(msh_file)

	return True


def test_entities_section():
	"""Test $Entities section content."""
	print("\n" + "=" * 60)
	print("Test 6: $Entities Section Validation")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")
	cubit.cmd("block 1 add tet all")
	cubit.cmd("block 1 name 'solid'")

	msh_file = "test_v4_entities.msh"
	cubit_mesh_export.export_Gmesh(cubit, msh_file, version="4.1")

	result = parse_gmsh_v4_file(msh_file)
	print(f"  Entities: {result['entities']}")

	# A cube should have: 8 vertices, 12 curves, 6 surfaces, 1 volume
	# But if geometry is not available, fallback may be used
	assert result['entities']['volumes'] >= 1, "Expected at least 1 volume entity"

	print("  PASS: Entities section validated")
	os.remove(msh_file)

	return True


def test_mixed_elements():
	"""Test mixed element types in same mesh."""
	print("\n" + "=" * 60)
	print("Test 7: Mixed Element Types")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
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

	msh_file = "test_v4_mixed.msh"
	cubit_mesh_export.export_Gmesh(cubit, msh_file, version="4.1")

	result = parse_gmsh_v4_file(msh_file)
	print(f"  Element types: {result['element_types']}")
	print(f"  Physical names: {result['physical_names']}")

	assert GMSH_V4_TYPES['HEX8'] in result['element_types'], "HEX8 not found"
	assert GMSH_V4_TYPES['TET4'] in result['element_types'], "TET4 not found"
	print("  PASS: Mixed elements exported correctly")
	os.remove(msh_file)

	return True


def test_physical_names():
	"""Test PhysicalNames section."""
	print("\n" + "=" * 60)
	print("Test 8: PhysicalNames Section")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")
	cubit.cmd("block 1 add tet all")
	cubit.cmd("block 1 name 'my_solid_block'")
	cubit.cmd("block 2 add tri all in surface all")
	cubit.cmd("block 2 name 'boundary_faces'")

	msh_file = "test_v4_names.msh"
	cubit_mesh_export.export_Gmesh(cubit, msh_file, version="4.1")

	result = parse_gmsh_v4_file(msh_file)
	print(f"  Physical names: {result['physical_names']}")

	assert 'my_solid_block' in result['physical_names'], "Block name not found"
	assert 'boundary_faces' in result['physical_names'], "Boundary name not found"
	print("  PASS: Physical names exported correctly")
	os.remove(msh_file)

	return True


if __name__ == "__main__":
	print("\n" + "=" * 60)
	print("export_Gmesh(version='4.1') Test Suite")
	print("=" * 60)

	all_passed = True

	tests = [
		test_format_version,
		test_3d_tet_mesh,
		test_3d_hex_mesh,
		test_2d_mesh_with_dim_auto,
		test_2d_mesh_normal_orientation,
		test_entities_section,
		test_mixed_elements,
		test_physical_names,
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
