"""
Compare two export paths for volume accuracy:
  Path A: extract_curved_mesh() (full Netgen, import cubit) -> .vol Save
  Path B: export netgen (compact_netgen, C++ plugin) -> .vol

Both should produce identical curving. Diagnose any difference.
"""
import sys
import os
import math

# Setup
_test_dir = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.dirname(os.path.dirname(_test_dir))
sys.path.insert(0, os.path.join(_repo_root, 'src', 'radia'))
from install_panels import find_cubit_bin
_cubit_path = find_cubit_bin()
if _cubit_path: sys.path.append(_cubit_path)
os.environ['CUBIT_PLUGIN_DIR'] = r'C:\Program Files\Coreform Cubit 2025.3\bin\plugins'

# NGSolve FIRST
import netgen.meshing
from ngsolve import Mesh, Integrate, CF

# Then Cubit
import cubit
cubit.init(['cubit', '-nojournal', '-batch',
            '-commandplugindir', os.environ['CUBIT_PLUGIN_DIR']])
from cubit_mesh_export import extract_curved_mesh

RADIUS = 0.05
V_EXACT = (4.0 / 3.0) * math.pi * RADIUS ** 3
OUT_DIR = os.path.join(_test_dir, 'export_netgen_test')
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 70)
print("Compare: extract_curved_mesh vs export netgen")
print("=" * 70)
print(f"Sphere R={RADIUS}, V_exact={V_EXACT:.10e}")

# Create mesh once
cubit.cmd("reset")
cubit.cmd(f"create sphere radius {RADIUS}")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size auto factor 5")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add volume 1")
cubit.cmd('block 1 name "sphere"')
print(f"Tets: {cubit.get_tet_count()}")

for order in [2, 3]:
    print(f"\n{'='*70}")
    print(f"Order {order}")
    print(f"{'='*70}")

    # --- Path A: extract_curved_mesh -> Save .vol ---
    ng_mesh_a = extract_curved_mesh(cubit, order=order)
    mesh_a_mem = Mesh(ng_mesh_a)
    vol_a_mem = Integrate(CF(1), mesh_a_mem)

    path_a = os.path.join(OUT_DIR, f"compare_pathA_o{order}.vol")
    ng_mesh_a.Save(path_a)

    mesh_a_reload = Mesh(path_a)
    vol_a_reload = Integrate(CF(1), mesh_a_reload)

    # --- Path B: export netgen (compact_netgen) ---
    path_b = os.path.join(OUT_DIR, f"compare_pathB_o{order}.vol")
    cubit.cmd(f'export netgen "{path_b}" order {order} overwrite')
    mesh_b = Mesh(path_b)
    vol_b = Integrate(CF(1), mesh_b)

    # --- Compare ---
    err_a_mem = (vol_a_mem - V_EXACT) / V_EXACT * 100
    err_a_reload = (vol_a_reload - V_EXACT) / V_EXACT * 100
    err_b = (vol_b - V_EXACT) / V_EXACT * 100

    print(f"  Path A (in-memory):  V={vol_a_mem:.10e}  err={err_a_mem:+.6e}%")
    print(f"  Path A (reload):     V={vol_a_reload:.10e}  err={err_a_reload:+.6e}%")
    print(f"  Path B (compact):    V={vol_b:.10e}  err={err_b:+.6e}%")

    # --- Structural comparison of .vol files ---
    def parse_curvedelements(path):
        with open(path) as f:
            lines = f.readlines()
        info = {}
        for i, line in enumerate(lines):
            if line.strip() == 'curvedelements':
                # edgeorder
                nedges = int(lines[i+1].strip())
                eo = [int(lines[i+2+j].strip()) for j in range(nedges)]
                info['nedges'] = nedges
                info['n_curved_edges'] = sum(1 for x in eo if x > 1)
                # faceorder
                off = i + 2 + nedges
                nfaces = int(lines[off].strip())
                fo = [int(lines[off+1+j].strip()) for j in range(nfaces)]
                info['nfaces'] = nfaces
                info['n_curved_faces'] = sum(1 for x in fo if x > 1)
                # edgecoeffsindex
                off2 = off + 1 + nfaces
                neci = int(lines[off2].strip())
                eci = [int(lines[off2+1+j].strip()) for j in range(neci)]
                info['edgecoeffsindex_size'] = neci
                info['total_edgecoeffs'] = eci[-1] if eci else 0
                # facecoeffsindex
                off3 = off2 + 1 + neci
                nfci = int(lines[off3].strip())
                fci = [int(lines[off3+1+j].strip()) for j in range(nfci)]
                info['facecoeffsindex_size'] = nfci
                info['total_facecoeffs'] = fci[-1] if fci else 0
                break
        else:
            info['nedges'] = 0
        return info

    info_a = parse_curvedelements(path_a)
    info_b = parse_curvedelements(path_b)

    print(f"\n  Structural comparison:")
    for key in sorted(set(list(info_a.keys()) + list(info_b.keys()))):
        va = info_a.get(key, 'N/A')
        vb = info_b.get(key, 'N/A')
        match = "OK" if va == vb else "DIFF"
        print(f"    {key:<25} A={va:<10} B={vb:<10} {match}")

    # --- Compare mesh element counts ---
    print(f"\n  Mesh comparison:")
    print(f"    Elements:    A={mesh_a_reload.ne:<8} B={mesh_b.ne:<8}")
    print(f"    Vertices:    A={mesh_a_reload.nv:<8} B={mesh_b.nv:<8}")
    print(f"    Materials A: {set(mesh_a_reload.GetMaterials())}")
    print(f"    Materials B: {set(mesh_b.GetMaterials())}")
