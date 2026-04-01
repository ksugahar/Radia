"""
Netgen Export Example: Sphere with export_NGSolveCurvedMesh()

Demonstrates automatic high-order curving for a sphere using export_NGSolveCurvedMesh().

Run: python netgen_sphere_setgeominfo.py
"""

import sys
import os
import math

# Auto-detect Cubit installation
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))
from install_panels import find_cubit_bin as _fcb
_cubit_path = _fcb()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)

work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(work_dir))
sys.path.insert(0, repo_root)

from ngsolve import Mesh, Integrate, CF, BND
import cubit
from cubit_netgen_bridge import extract_curved_mesh

# ============================================================
# Parameters
# ============================================================
R = 0.5  # Radius

print("=== Netgen Export: Sphere with export_NGSolveCurvedMesh() ===")
print(f"Parameters: R={R}")
print()

# ============================================================
# Step 1: Create and mesh geometry in Cubit
# ============================================================
print("Step 1: Create and mesh geometry in Cubit")

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create sphere radius {R}")

cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.1")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add tet all")
cubit.cmd('block 1 name "domain"')
cubit.cmd("block 2 add tri all")
cubit.cmd('block 2 name "boundary"')

print(f"  Tets: {cubit.get_tet_count()}")

# ============================================================
# Step 2: Export with automatic curving
# ============================================================
expected_area = 4 * math.pi * R * R
expected_vol = 4 / 3 * math.pi * R * R * R

for order in [2, 3]:
	print(f"\nexport_NGSolveCurvedMesh(order={order})")
	mesh = Mesh(extract_curved_mesh(cubit, order=order))

	area = Integrate(CF(1), mesh, VOL_or_BND=BND)
	vol = Integrate(CF(1), mesh)

	print(f"  Area: {area:.6f} (err {abs(area-expected_area)/expected_area*100:.4f}%)")
	print(f"  Vol:  {vol:.6f} (err {abs(vol-expected_vol)/expected_vol*100:.4f}%)")

print("\n=== DONE ===")
