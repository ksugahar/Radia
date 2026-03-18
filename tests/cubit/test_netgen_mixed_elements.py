"""
Test script for export_NetgenMesh() with mixed element types (Hex, Wedge, Pyramid, Tet)

Run this script:
  python test_ngsolve_mixed_elements.py
"""

import sys
import os

# IMPORTANT: Import NGSolve/Netgen BEFORE Cubit to avoid DLL conflicts
import ngsolve

cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
	sys.path.append(cubit_path)

print("=" * 60)
print("Testing export_NetgenMesh() with mixed elements")
print("=" * 60)

# Step 1: Import and initialize Cubit
print("\nStep 1: Importing Cubit...")
import cubit
cubit.init(['cubit', '-nojournal', '-batch'])
print("  [OK] Cubit initialized")

# Step 2: Create geometry with multiple volumes
print("\nStep 2: Creating geometry with hex, wedge, pyramid, tet...")
cubit.cmd("reset")

# Create brick for hex mesh
cubit.cmd("brick x 1 y 1 z 1")
cubit.cmd("volume 1 move 0 0 0.5")

# Create brick for tet/pyramid mesh
cubit.cmd("brick x 1 y 1 z 1")
cubit.cmd("volume 2 move 0 0 1.5")

# Create wedge by sweeping (cylinder approximation)
cubit.cmd("create cylinder radius 0.4 height 1")
cubit.cmd("volume 3 move 0 0 -0.5")

# Imprint and merge
cubit.cmd("imprint all")
cubit.cmd("merge all")
cubit.cmd("compress")

# Mesh volume 1 with hex
cubit.cmd("volume 1 scheme map")
cubit.cmd("volume 1 size 0.5")
cubit.cmd("mesh volume 1")

# Mesh volume 2 with tet (will have some pyramids at interface)
cubit.cmd("volume 2 scheme tetmesh")
cubit.cmd("volume 2 size 0.5")
cubit.cmd("mesh volume 2")

# Mesh volume 3 with sweep (creates wedges)
cubit.cmd("volume 3 scheme sweep")
cubit.cmd("volume 3 size 0.3")
cubit.cmd("mesh volume 3")

# Add elements to blocks (1st order only for this test)
cubit.cmd("block 1 add hex all")
cubit.cmd("block 1 name 'hex_region'")

cubit.cmd("block 2 add tet all")
cubit.cmd("block 2 name 'tet_region'")

cubit.cmd("block 3 add wedge all")
cubit.cmd("block 3 name 'wedge_region'")

cubit.cmd("block 4 add pyramid all")
cubit.cmd("block 4 name 'pyramid_region'")

# Add boundary elements
cubit.cmd("block 5 add tri all in surface all")
cubit.cmd("block 5 name 'tri_boundary'")

cubit.cmd("block 6 add face all in surface all")
cubit.cmd("block 6 name 'quad_boundary'")

# Get mesh statistics
num_hex = len(cubit.get_block_hexes(1))
num_tet = len(cubit.get_block_tets(2))
num_wedge = len(cubit.get_block_wedges(3))
num_pyramid = len(cubit.get_block_pyramids(4))
num_tri = len(cubit.get_block_tris(5))
num_quad = len(cubit.get_block_faces(6))

print(f"  Created mesh:")
print(f"    Hex:     {num_hex}")
print(f"    Tet:     {num_tet}")
print(f"    Wedge:   {num_wedge}")
print(f"    Pyramid: {num_pyramid}")
print(f"    Tri:     {num_tri}")
print(f"    Quad:    {num_quad}")

# Step 3: Export geometry to STEP
print("\nStep 3: Exporting geometry to STEP...")
step_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_mixed.step")
cubit.cmd(f'export step "{step_file}" overwrite')
print(f"  [OK] Exported to {step_file}")

# Step 4: Use export_Netgen function
print("\nStep 4: Using export_NetgenMesh() function...")
import cubit_mesh_export

ngmesh = cubit_mesh_export.export_NetgenMesh(cubit, geometry_file=step_file)
print(f"  [OK] Created netgen.Mesh: ne={ngmesh.ne}")

# Step 5: Convert to ngsolve.Mesh
print("\nStep 5: Converting to ngsolve.Mesh...")
ngs_mesh = ngsolve.Mesh(ngmesh)
print(f"  [OK] Converted: nv={ngs_mesh.nv}, ne={ngs_mesh.ne}")

# Step 6: Test Curve() with different orders
print("\nStep 6: Testing mesh.Curve() with different orders...")
for order in [2, 3, 4]:
	try:
		ngs_mesh.Curve(order)
		print(f"  [OK] mesh.Curve({order}) succeeded!")
	except Exception as e:
		print(f"  [FAIL] mesh.Curve({order}) error: {e}")

# Step 7: Check materials and boundaries
print("\nStep 7: Checking materials and boundaries...")
print(f"  Materials: {ngs_mesh.GetMaterials()}")
print(f"  Boundaries: {ngs_mesh.GetBoundaries()}")

# Step 8: Verify element counts
print("\nStep 8: Verifying element counts...")
total_vol_elements = num_hex + num_tet + num_wedge + num_pyramid
print(f"  Cubit total 3D elements: {total_vol_elements}")
print(f"  NGSolve ne: {ngs_mesh.ne}")
if ngs_mesh.ne == total_vol_elements:
	print("  [OK] Element counts match!")
else:
	print(f"  [WARN] Element count mismatch: expected {total_vol_elements}, got {ngs_mesh.ne}")

print("\n" + "=" * 60)
print("Test complete!")
print("=" * 60)
