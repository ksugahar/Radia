"""
Test script to investigate Cubit API behavior when blocks contain geometry.

This script tests various API functions to find the best way to get
mesh elements from blocks that contain geometry (volume, surface, curve, vertex)
instead of direct mesh element references.

Run with Cubit Python:
    "C:/Program Files/Coreform Cubit 2025.12/bin/python3/python.exe" tests/test_block_geometry_api.py
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
cubit.init(['cubit', '-nojournal', '-batch'])

print("=" * 70)
print("Testing Cubit API for blocks with geometry")
print("=" * 70)

# Create a simple geometry and mesh
cubit.cmd("reset")
cubit.cmd("brick x 1 y 1 z 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.5")
cubit.cmd("mesh volume 1")

print("\n--- Mesh statistics ---")
print(f"Total nodes: {cubit.get_node_count()}")
print(f"Total tets: {cubit.get_tet_count()}")
print(f"Total tris: {cubit.get_tri_count()}")

# Test 1: Block with volume (geometry)
print("\n" + "=" * 70)
print("Test 1: Block with VOLUME (geometry)")
print("=" * 70)
cubit.cmd("block 1 add volume 1")
cubit.cmd("block 1 name 'volume_block'")

block_id = 1
print(f"\nBlock {block_id} contents:")
print(f"  get_block_volumes({block_id}): {cubit.get_block_volumes(block_id)}")
print(f"  get_block_tets({block_id}): {cubit.get_block_tets(block_id)}")
print(f"  get_block_hexes({block_id}): {cubit.get_block_hexes(block_id)}")
print(f"  get_block_tris({block_id}): {cubit.get_block_tris(block_id)}")

# Test parse_cubit_list
print(f"\nUsing parse_cubit_list:")
try:
	tets_in_block = cubit.parse_cubit_list("tet", f"in block {block_id}")
	print(f"  parse_cubit_list('tet', 'in block {block_id}'): {len(tets_in_block)} tets")
	print(f"    First 5: {tets_in_block[:5]}")
except Exception as e:
	print(f"  parse_cubit_list failed: {e}")

try:
	hexes_in_block = cubit.parse_cubit_list("hex", f"in block {block_id}")
	print(f"  parse_cubit_list('hex', 'in block {block_id}'): {len(hexes_in_block)} hexes")
except Exception as e:
	print(f"  parse_cubit_list failed: {e}")

try:
	tris_in_block = cubit.parse_cubit_list("tri", f"in block {block_id}")
	print(f"  parse_cubit_list('tri', 'in block {block_id}'): {len(tris_in_block)} tris")
except Exception as e:
	print(f"  parse_cubit_list failed: {e}")

try:
	nodes_in_block = cubit.parse_cubit_list("node", f"in block {block_id}")
	print(f"  parse_cubit_list('node', 'in block {block_id}'): {len(nodes_in_block)} nodes")
except Exception as e:
	print(f"  parse_cubit_list failed: {e}")

# Test get_block_elements_and_nodes (if available)
print(f"\nUsing get_block_elements_and_nodes:")
try:
	result = cubit.get_block_elements_and_nodes(block_id)
	print(f"  Result type: {type(result)}")
	print(f"  Result: {result}")
except AttributeError:
	print("  Function not available in this version")
except Exception as e:
	print(f"  Error: {e}")


# Test 2: Block with surface (geometry)
print("\n" + "=" * 70)
print("Test 2: Block with SURFACE (geometry)")
print("=" * 70)
cubit.cmd("block 2 add surface 1")
cubit.cmd("block 2 name 'surface_block'")

block_id = 2
print(f"\nBlock {block_id} contents:")
print(f"  get_block_surfaces({block_id}): {cubit.get_block_surfaces(block_id)}")
print(f"  get_block_tris({block_id}): {cubit.get_block_tris(block_id)}")
print(f"  get_block_faces({block_id}): {cubit.get_block_faces(block_id)}")

print(f"\nUsing parse_cubit_list:")
try:
	tris_in_block = cubit.parse_cubit_list("tri", f"in block {block_id}")
	print(f"  parse_cubit_list('tri', 'in block {block_id}'): {len(tris_in_block)} tris")
except Exception as e:
	print(f"  parse_cubit_list failed: {e}")

try:
	nodes_in_block = cubit.parse_cubit_list("node", f"in block {block_id}")
	print(f"  parse_cubit_list('node', 'in block {block_id}'): {len(nodes_in_block)} nodes")
except Exception as e:
	print(f"  parse_cubit_list failed: {e}")


# Test 3: Block with curve (geometry)
print("\n" + "=" * 70)
print("Test 3: Block with CURVE (geometry)")
print("=" * 70)
cubit.cmd("block 3 add curve 1")
cubit.cmd("block 3 name 'curve_block'")

block_id = 3
print(f"\nBlock {block_id} contents:")
try:
	print(f"  get_block_curves({block_id}): {cubit.get_block_curves(block_id)}")
except:
	print(f"  get_block_curves not available")
print(f"  get_block_edges({block_id}): {cubit.get_block_edges(block_id)}")

print(f"\nUsing parse_cubit_list:")
try:
	edges_in_block = cubit.parse_cubit_list("edge", f"in block {block_id}")
	print(f"  parse_cubit_list('edge', 'in block {block_id}'): {len(edges_in_block)} edges")
except Exception as e:
	print(f"  parse_cubit_list failed: {e}")

try:
	nodes_in_block = cubit.parse_cubit_list("node", f"in block {block_id}")
	print(f"  parse_cubit_list('node', 'in block {block_id}'): {len(nodes_in_block)} nodes")
except Exception as e:
	print(f"  parse_cubit_list failed: {e}")


# Test 4: Block with vertex (geometry)
print("\n" + "=" * 70)
print("Test 4: Block with VERTEX (geometry)")
print("=" * 70)
cubit.cmd("block 4 add vertex 1")
cubit.cmd("block 4 name 'vertex_block'")

block_id = 4
print(f"\nBlock {block_id} contents:")
try:
	print(f"  get_block_vertices({block_id}): {cubit.get_block_vertices(block_id)}")
except:
	print(f"  get_block_vertices not available")
print(f"  get_block_nodes({block_id}): {cubit.get_block_nodes(block_id)}")

print(f"\nUsing parse_cubit_list:")
try:
	nodes_in_block = cubit.parse_cubit_list("node", f"in block {block_id}")
	print(f"  parse_cubit_list('node', 'in block {block_id}'): {len(nodes_in_block)} nodes")
except Exception as e:
	print(f"  parse_cubit_list failed: {e}")


# Test 5: Block with mesh elements (direct)
print("\n" + "=" * 70)
print("Test 5: Block with MESH ELEMENTS (direct)")
print("=" * 70)
cubit.cmd("block 5 add tet all")
cubit.cmd("block 5 name 'mesh_block'")

block_id = 5
print(f"\nBlock {block_id} contents:")
print(f"  get_block_volumes({block_id}): {cubit.get_block_volumes(block_id)}")
print(f"  get_block_tets({block_id}): {len(cubit.get_block_tets(block_id))} tets")

print(f"\nUsing parse_cubit_list:")
try:
	tets_in_block = cubit.parse_cubit_list("tet", f"in block {block_id}")
	print(f"  parse_cubit_list('tet', 'in block {block_id}'): {len(tets_in_block)} tets")
except Exception as e:
	print(f"  parse_cubit_list failed: {e}")


# Summary
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
Key findings:
1. get_block_tets/hexes/tris returns empty for geometry-based blocks
2. parse_cubit_list('tet', 'in block X') should work for both cases
3. This mimics Cubit's 'draw tet in block X' behavior
""")

print("\nDone!")
