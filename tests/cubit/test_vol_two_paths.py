"""
Debug: compare two .vol generation paths from the SAME Cubit model.

Path 1 (jou -> export netgen):
  Cubit creates mesh -> "export netgen" C++ command -> .vol (compact_netgen)

Path 2 (cub5 -> import cubit -> extract_curved_mesh):
  Cubit creates mesh -> save cub5 -> fresh Python loads cub5 ->
  extract_curved_mesh() -> ng_mesh.Save() -> .vol (full Netgen)

Both paths use the same sphere model. Volume error should match.

Usage: python tests/cubit/test_vol_two_paths.py
"""
import sys
import os
import math
import subprocess
import json
import tempfile

_test_dir = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.dirname(os.path.dirname(_test_dir))

RADIUS = 0.05
V_EXACT = (4.0 / 3.0) * math.pi * RADIUS ** 3
OUT_DIR = os.path.join(_test_dir, 'export_netgen_test')
os.makedirs(OUT_DIR, exist_ok=True)

# ================================================================
# Step 1: Create mesh in Cubit, export netgen + save cub5
# ================================================================
print("Step 1: Create mesh in Cubit")

sys.path.insert(0, os.path.join(_repo_root, 'src', 'radia'))
from install_panels import find_cubit_bin
_cubit_path = find_cubit_bin()
if _cubit_path: sys.path.append(_cubit_path)
os.environ['CUBIT_PLUGIN_DIR'] = r'C:\Program Files\Coreform Cubit 2025.3\bin\plugins'

# Import NGSolve first
import netgen.meshing
from ngsolve import Mesh, Integrate, CF

import cubit
cubit.init(['cubit', '-nojournal', '-batch',
            '-commandplugindir', os.environ['CUBIT_PLUGIN_DIR']])

cubit.cmd("reset")
cubit.cmd(f"create sphere radius {RADIUS}")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size auto factor 5")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add volume 1")
cubit.cmd('block 1 name "sphere"')
print(f"  Tets: {cubit.get_tet_count()}")

# Save cub5 for Path 2
cub5_path = os.path.join(OUT_DIR, "sphere_debug.cub5")
cubit.cmd(f'save cub5 "{cub5_path}" overwrite')
print(f"  Saved: {cub5_path}")

for order in [2, 3]:
    print(f"\n{'='*70}")
    print(f"Order {order}")
    print(f"{'='*70}")

    # --- Path 1: export netgen (compact_netgen) ---
    path1 = os.path.join(OUT_DIR, f"debug_path1_o{order}.vol")
    cubit.cmd(f'export netgen "{path1}" order {order} overwrite')
    mesh1 = Mesh(path1)
    vol1 = Integrate(CF(1), mesh1)
    err1 = (vol1 - V_EXACT) / V_EXACT * 100
    print(f"  Path 1 (export netgen):      V={vol1:.10e}  err={err1:+.6e}%")

    # --- Path 2: extract_curved_mesh (import cubit) ---
    # Run in same process since cubit is already imported
    from cubit_mesh_export import extract_curved_mesh
    ng_mesh2 = extract_curved_mesh(cubit, order=order)
    mesh2_mem = Mesh(ng_mesh2)
    vol2_mem = Integrate(CF(1), mesh2_mem)
    err2_mem = (vol2_mem - V_EXACT) / V_EXACT * 100

    path2 = os.path.join(OUT_DIR, f"debug_path2_o{order}.vol")
    ng_mesh2.Save(path2)
    mesh2_reload = Mesh(path2)
    vol2_reload = Integrate(CF(1), mesh2_reload)
    err2_reload = (vol2_reload - V_EXACT) / V_EXACT * 100

    print(f"  Path 2 (extract, in-mem):    V={vol2_mem:.10e}  err={err2_mem:+.6e}%")
    print(f"  Path 2 (extract, reload):    V={vol2_reload:.10e}  err={err2_reload:+.6e}%")

    # --- Structural comparison ---
    def count_curved(path):
        with open(path) as f:
            lines = f.readlines()
        info = {}
        for i, line in enumerate(lines):
            if line.strip() == 'curvedelements':
                nedges = int(lines[i+1].strip())
                nho = sum(1 for j in range(i+2, i+2+nedges) if int(lines[j].strip()) > 1)
                info['nedges'] = nedges
                info['n_curved_edges'] = nho
                break
        # Count mesh stats from header
        for i, line in enumerate(lines):
            if line.strip() == 'volumeelements':
                info['ne'] = int(lines[i+1].strip())
                break
        for i, line in enumerate(lines):
            if line.strip() == 'surfaceelements':
                info['nse'] = int(lines[i+1].strip())
                break
        for i, line in enumerate(lines):
            if line.strip() == 'edgesegmentsgi2':
                info['nseg'] = int(lines[i+1].strip())
                break
        for i, line in enumerate(lines):
            if line.strip() == 'points':
                info['np'] = int(lines[i+1].strip())
                break
        return info

    info1 = count_curved(path1)
    info2 = count_curved(path2)

    print(f"\n  .vol structure:")
    print(f"    {'':25} {'Path1':>10} {'Path2':>10} {'Match':>6}")
    for key in ['np', 'ne', 'nse', 'nseg', 'nedges', 'n_curved_edges']:
        v1 = info1.get(key, 'N/A')
        v2 = info2.get(key, 'N/A')
        m = "OK" if v1 == v2 else "DIFF"
        print(f"    {key:25} {str(v1):>10} {str(v2):>10} {m:>6}")

print(f"\n{'='*70}")
print("Root cause: 'nseg' difference = segment elements on curves")
print("  Path 1 (NetgenCurver) does not add AddSegment()")
print("  Path 2 (NetgenCurverPure) adds segments from edge_elements")
print("  -> topology.Update() produces different edge count")
print("  -> BuildCurvedElements curves fewer edges in Path 1")
