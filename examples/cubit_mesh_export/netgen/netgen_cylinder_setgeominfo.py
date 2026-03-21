"""
Netgen Export Example: Cylinder with export_NGSolveCurvedMesh()

Demonstrates the unified high-order curving workflow using export_NGSolveCurvedMesh(),
which works directly with Cubit's ACIS kernel for UV parameter computation
and mesh.Curve(order). No STEP files needed.

Run: python netgen_cylinder_setgeominfo.py
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

print("=== Netgen Export: Cylinder with export_NGSolveCurvedMesh() ===")
print(f"Parameters: R={R}, H={H}")
print()

# ============================================================
# Step 1: Create and mesh geometry in Cubit
# ============================================================
print("Step 1: Create and mesh geometry in Cubit")

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create cylinder height {H} radius {R}")

cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.15")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add tet all")
cubit.cmd('block 1 name "domain"')
cubit.cmd("block 2 add tri all")
cubit.cmd('block 2 name "boundary"')

print(f"  Tets: {cubit.get_tet_count()}")

# ============================================================
# Step 2: Export with automatic curving (order=2)
# ============================================================
print("\nStep 2: export_NGSolveCurvedMesh(order=2)")

mesh2 = cubit_mesh_export.export_NGSolveCurvedMesh(cubit, order=2)
print(f"  Boundaries: {mesh2.GetBoundaries()}")

# ============================================================
# Step 2b: Accuracy comparison
# ============================================================
print("\nStep 2b: Accuracy")

expected_area = 2*math.pi*R*H + 2*math.pi*R*R
expected_vol = math.pi*R*R*H

area2 = Integrate(CF(1), mesh2, VOL_or_BND=BND)
vol2 = Integrate(CF(1), mesh2)

print(f"  Expected:  Area={expected_area:.6f}, Vol={expected_vol:.6f}")
print(f"  order=2:   Area={area2:.6f} ({abs(area2-expected_area)/expected_area*100:.4f}% err), "
      f"Vol={vol2:.6f} ({abs(vol2-expected_vol)/expected_vol*100:.4f}% err)")

# ============================================================
# Step 3: Compare with order=3
# ============================================================
print("\nStep 3: export_NGSolveCurvedMesh(order=3)")

mesh3 = cubit_mesh_export.export_NGSolveCurvedMesh(cubit, order=3)

area3 = Integrate(CF(1), mesh3, VOL_or_BND=BND)
vol3 = Integrate(CF(1), mesh3)

print(f"  order=3:   Area={area3:.6f} ({abs(area3-expected_area)/expected_area*100:.4f}% err), "
      f"Vol={vol3:.6f} ({abs(vol3-expected_vol)/expected_vol*100:.4f}% err)")

print("\n=== DONE ===")
