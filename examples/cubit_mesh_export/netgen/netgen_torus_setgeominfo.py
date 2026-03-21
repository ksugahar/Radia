"""
Netgen Export Example: Torus with export_curved()

Demonstrates automatic high-order curving for a torus using export_curved().

Run: python netgen_torus_setgeominfo.py
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
R_major = 1.0  # Major radius
R_minor = 0.3  # Minor radius

print("=== Netgen Export: Torus with export_curved() ===")
print(f"Parameters: R_major={R_major}, R_minor={R_minor}")
print()

# ============================================================
# Step 1: Create and mesh geometry in Cubit
# ============================================================
print("Step 1: Create and mesh geometry in Cubit")

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create torus major radius {R_major} minor radius {R_minor}")

cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.08")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add tet all")
cubit.cmd('block 1 name "domain"')
cubit.cmd("block 2 add tri all")
cubit.cmd('block 2 name "boundary"')

print(f"  Tets: {cubit.get_tet_count()}")

# ============================================================
# Step 2: Export with automatic curving
# ============================================================
expected_area = 4 * math.pi**2 * R_major * R_minor
expected_vol = 2 * math.pi**2 * R_major * R_minor**2

for order in [2, 3]:
	print(f"\nexport_curved(order={order})")
	mesh = cubit_mesh_export.export_curved(cubit, order=order)

	area = Integrate(CF(1), mesh, VOL_or_BND=BND)
	vol = Integrate(CF(1), mesh)

	print(f"  Area: {area:.6f} (err {abs(area-expected_area)/expected_area*100:.4f}%)")
	print(f"  Vol:  {vol:.6f} (err {abs(vol-expected_vol)/expected_vol*100:.4f}%)")

print("\n=== DONE ===")
