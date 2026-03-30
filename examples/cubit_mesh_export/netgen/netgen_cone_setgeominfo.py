"""
Netgen Export Example: Cone with export_NGSolveCurvedMesh()

Demonstrates automatic high-order curving for a cone using export_NGSolveCurvedMesh().

Run: python netgen_cone_setgeominfo.py
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

from ngsolve import Integrate, CF, BND
import cubit
import cubit_mesh_export

# ============================================================
# Parameters
# ============================================================
R = 0.5  # Base radius
H = 2.0  # Height

print("=== Netgen Export: Cone with export_NGSolveCurvedMesh() ===")
print(f"Parameters: base_radius={R}, height={H}")
print()

# ============================================================
# Step 1: Create and mesh geometry in Cubit
# ============================================================
print("Step 1: Create and mesh geometry in Cubit")

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create frustum height {H} radius {R} top 0")

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
slant_height = math.sqrt(R*R + H*H)
expected_area = math.pi * R * R + math.pi * R * slant_height
expected_vol = (1/3) * math.pi * R * R * H

for order in [2, 3]:
	print(f"\nexport_NGSolveCurvedMesh(order={order})")
	mesh = cubit_mesh_export.export_NGSolveCurvedMesh(cubit, order=order)

	area = Integrate(CF(1), mesh, VOL_or_BND=BND)
	vol = Integrate(CF(1), mesh)

	print(f"  Area: {area:.6f} (err {abs(area-expected_area)/expected_area*100:.4f}%)")
	print(f"  Vol:  {vol:.6f} (err {abs(vol-expected_vol)/expected_vol*100:.4f}%)")

print("\n=== DONE ===")
