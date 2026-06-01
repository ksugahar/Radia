"""
Test the warning for mixed element types in blocks.

Run with Cubit Python:
    "C:/Program Files/Coreform Cubit 2025.12/bin/python3/python.exe" tests/test_mixed_element_warning.py
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

# Initialize Cubit
cubit.init(['cubit', '-nojournal', '-batch'])

print("=" * 70)
print("Test: Mixed Element Type Warning")
print("=" * 70)

# Test 1: Single element type in block (no warning expected)
print("\n" + "-" * 70)
print("Test 1: Single element type in block (no warning)")
print("-" * 70)

cubit.cmd("reset")
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.5")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 name 'solid'")

with tempfile.NamedTemporaryFile(suffix='.msh', delete=False) as f:
	vtk_file = f.name

print("  Calling export gmsh (should see no warning):")
cubit.cmd(f'export gmsh "{vtk_file}" overwrite')
os.unlink(vtk_file)

# Test 2: Mixed element types in block (warning expected)
print("\n" + "-" * 70)
print("Test 2: Mixed element types in block (warning expected)")
print("-" * 70)

cubit.cmd("reset")
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 move 0 0 0")
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 2 move 1.5 0 0")

cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 2 scheme map")
cubit.cmd("volume all size 0.5")
cubit.cmd("mesh volume all")

# Add both tets and hexes to the same block
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 add hex all")
cubit.cmd("block 1 name 'mixed'")

with tempfile.NamedTemporaryFile(suffix='.msh', delete=False) as f:
	vtk_file = f.name

print("  Calling export gmsh (should see warning about mixed types):")
cubit.cmd(f'export gmsh "{vtk_file}" overwrite')
os.unlink(vtk_file)

# Test 3: Multiple separate blocks (no warning expected)
print("\n" + "-" * 70)
print("Test 3: Separate blocks for each element type (no warning)")
print("-" * 70)

cubit.cmd("reset")
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 move 0 0 0")
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 2 move 1.5 0 0")

cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 2 scheme map")
cubit.cmd("volume all size 0.5")
cubit.cmd("mesh volume all")

cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 name 'tets'")
cubit.cmd("block 2 add hex all")
cubit.cmd("block 2 name 'hexes'")

with tempfile.NamedTemporaryFile(suffix='.msh', delete=False) as f:
	vtk_file = f.name

print("  Calling export gmsh (should see no warning):")
cubit.cmd(f'export gmsh "{vtk_file}" overwrite')
os.unlink(vtk_file)

print("\n" + "=" * 70)
print("Test complete!")
print("=" * 70)
