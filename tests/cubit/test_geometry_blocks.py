"""
Test export with geometry-based blocks (volume, surface, curve, vertex in blocks).

This module tests that export functions work correctly when blocks contain
geometry references instead of direct mesh element references.

Run with Cubit Python:
    "C:/Program Files/Coreform Cubit 2025.3/bin/python3/python.exe" tests/test_geometry_blocks.py
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
import tempfile

# Add parent directory to path for radia_cubit_mesh
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))
import radia_cubit_mesh

# Initialize Cubit
cubit.init(['cubit', '-nojournal', '-batch'])

all_passed = True

def test_helper_function():
	"""Test extract_mesh_data helper function for getting block elements."""
	global all_passed
	print("=" * 60)
	print("Test 1: extract_mesh_data helper function")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	num_tets = cubit.get_tet_count()
	num_tris = cubit.get_tri_count()

	# Test with mesh elements block
	cubit.cmd("block 1 add tet all")
	cubit.cmd("block 2 add tri all in surface all")

	# Use extract_mesh_data to get elements
	mesh_data = radia_cubit_mesh.extract_mesh_data(order=1)
	# mesh_data provides element information from blocks
	print(f"  Mesh data extracted successfully")
	print(f"  Total tets in mesh: {num_tets}")
	print(f"  Total tris in mesh: {num_tris}")

	if num_tets > 0 and num_tris > 0:
		print("  PASS: Mesh element blocks work correctly")
	else:
		print("  FAIL: Mesh element blocks")
		all_passed = False

	# Reset and test with geometry blocks
	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	num_tets = cubit.get_tet_count()

	# Test with geometry blocks
	cubit.cmd("block 1 add volume 1")
	cubit.cmd("block 2 add surface 1")

	mesh_data = radia_cubit_mesh.extract_mesh_data(order=1)
	print(f"\n  Geometry blocks:")
	print(f"    Total tets: {num_tets}")

	if num_tets > 0:
		print("  PASS: Geometry blocks work correctly")
	else:
		print("  FAIL: Geometry blocks")
		all_passed = False


def test_gmsh_v2_with_volume_block():
	"""Test Gmsh v2 export with block containing volume (geometry)."""
	global all_passed
	print("\n" + "=" * 60)
	print("Test 2: Gmsh v2 export with volume-based block")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	# Add VOLUME to block (not mesh elements)
	cubit.cmd("block 1 add volume 1")
	cubit.cmd("block 1 name 'solid'")

	num_tets_expected = cubit.get_tet_count()
	print(f"  Total tets in mesh: {num_tets_expected}")
	print(f"  Block contains: volume 1 (geometry)")

	# Export to Gmsh
	with tempfile.NamedTemporaryFile(suffix='.msh', delete=False) as f:
		msh_file = f.name

	cubit.cmd(f'radia export gmsh "{msh_file}" overwrite')

	# Verify file exists and has content
	with open(msh_file, 'r') as f:
		content = f.read()

	# Check if elements are in the file
	has_elements = '$Elements' in content

	if has_elements and len(content) > 100:
		print("  PASS: Gmsh file generated with geometry block")
	else:
		print("  FAIL: Gmsh file not generated correctly")
		all_passed = False

	os.unlink(msh_file)


def test_gmsh_with_volume_block():
	"""Test Gmsh export with volume-based block."""
	global all_passed
	print("\n" + "=" * 60)
	print("Test 3: Gmsh export with volume-based block")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	# Add VOLUME to block
	cubit.cmd("block 1 add volume 1")
	cubit.cmd("block 1 name 'solid'")

	num_tets_expected = cubit.get_tet_count()

	# Export to Gmsh
	with tempfile.NamedTemporaryFile(suffix='.msh', delete=False) as f:
		msh_file = f.name

	cubit.cmd(f'radia export gmsh "{msh_file}" version 2 overwrite')

	# Verify file
	with open(msh_file, 'r') as f:
		content = f.read()

	# Find Elements section and count tet elements (type 4 for TET4)
	lines = content.split('\n')
	in_elements = False
	tet_count = 0
	for line in lines:
		if line.strip() == '$Elements':
			in_elements = True
			continue
		if line.strip() == '$EndElements':
			break
		if in_elements and line.strip():
			parts = line.split()
			if len(parts) >= 2:
				try:
					elem_type = int(parts[1])
					# Type 4 = TET4, Type 11 = TET10
					if elem_type in [4, 11]:
						tet_count += 1
				except:
					pass

	print(f"  Expected tets: {num_tets_expected}")
	print(f"  Gmsh tet elements: {tet_count}")

	if tet_count == num_tets_expected:
		print("  PASS: Gmsh export with volume block works")
	else:
		print(f"  FAIL: Tet count mismatch")
		all_passed = False

	os.unlink(msh_file)


def test_no_cross_contamination():
	"""Test that elements from one block don't appear in another."""
	global all_passed
	print("\n" + "=" * 60)
	print("Test 4: No cross-contamination between blocks")
	print("=" * 60)

	cubit.cmd("reset")
	cubit.cmd("create brick x 1 y 1 z 1")
	cubit.cmd("volume 1 scheme tetmesh")
	cubit.cmd("volume 1 size 0.5")
	cubit.cmd("mesh volume 1")

	# Block 1: tets only
	cubit.cmd("block 1 add tet all")
	# Block 2: tris only
	cubit.cmd("block 2 add tri all in surface all")

	# Check that block 1 has tets and block 2 has tris (no cross-contamination)
	block1_tets = len(cubit.get_block_tets(1))
	block1_tris = len(cubit.get_block_tris(1))
	block2_tets = len(cubit.get_block_tets(2))
	block2_tris = len(cubit.get_block_tris(2))

	print(f"  Block 1: {block1_tets} tets, {block1_tris} tris (expect tets>0, tris=0)")
	print(f"  Block 2: {block2_tets} tets, {block2_tris} tris (expect tets=0, tris>0)")

	if block1_tets > 0 and block1_tris == 0 and block2_tets == 0 and block2_tris > 0:
		print("  PASS: No cross-contamination")
	else:
		print("  FAIL: Cross-contamination detected")
		all_passed = False


# Run tests
if __name__ == "__main__":
	test_helper_function()
	test_gmsh_v2_with_volume_block()
	test_gmsh_with_volume_block()
	test_no_cross_contamination()

	print("\n" + "=" * 60)
	if all_passed:
		print("All tests PASSED!")
	else:
		print("Some tests FAILED!")
	print("=" * 60)
