"""
Netgen Export Example: Hex Cylinder with export_curved()

Key use case: Cubit's hex meshing + NGSolve high-order elements.
Demonstrates that export_curved() works with hexahedral meshes too.

Run: python netgen_hex_cylinder_setgeominfo.py
"""

import sys
import os
import math

cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
	sys.path.append(cubit_path)

work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(work_dir))
sys.path.insert(0, repo_root)

from ngsolve import Integrate, CF, BND
import cubit
import cubit_mesh_export

# ============================================================
# Parameters
# ============================================================
R = 0.5  # Radius
H = 2.0  # Height

print("=== Netgen Export: Hex Cylinder with export_curved() ===")
print(f"Parameters: R={R}, H={H}")
print("(Cubit hex meshing + NGSolve high-order)")
print()

# ============================================================
# Step 1: Create and mesh geometry in Cubit (hex sweep)
# ============================================================
print("Step 1: Create and mesh geometry in Cubit")

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create cylinder height {H} radius {R}")

# Hex mesh using sweep
cubit.cmd("volume 1 scheme sweep")
cubit.cmd("surface all size 0.15")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add hex all")
cubit.cmd('block 1 name "domain"')
cubit.cmd("block 2 add face all")
cubit.cmd('block 2 name "boundary"')

print(f"  Hexes: {cubit.get_hex_count()}")

# ============================================================
# Step 2: Export with automatic curving
# ============================================================
expected_area = 2 * math.pi * R * H + 2 * math.pi * R * R
expected_vol = math.pi * R * R * H

for order in [2, 3]:
	print(f"\nexport_curved(order={order})")
	mesh = cubit_mesh_export.export_curved(cubit, order=order)

	area = Integrate(CF(1), mesh, VOL_or_BND=BND)
	vol = Integrate(CF(1), mesh)

	print(f"  Area: {area:.6f} (err {abs(area-expected_area)/expected_area*100:.4f}%)")
	print(f"  Vol:  {vol:.6f} (err {abs(vol-expected_vol)/expected_vol*100:.4f}%)")

print("\n=== DONE ===")
