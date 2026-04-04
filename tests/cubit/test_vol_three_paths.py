"""
Compare three .vol generation paths for a sphere (R=0.05m).

Path A: export netgen (compact_netgen C++ plugin)
  Cubit -> "export netgen" APREPRO command -> .vol

Path B: extract_curved_mesh (full Netgen, import cubit)
  Cubit -> extract_curved_mesh(cubit, order) -> ng_mesh.Save() -> .vol

Path C: CSG Sphere -> Netgen GenerateMesh -> Curve(order) -> Save()
  Independent Netgen mesh with exact sphere geometry (CSG kernel).
  Different mesh, but volume/area are valid reference values.

Usage: python tests/cubit/test_vol_three_paths.py
"""
import sys
import os
import math

_test_dir = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.dirname(os.path.dirname(_test_dir))
sys.path.insert(0, os.path.join(_repo_root, 'src', 'radia'))
from install_panels import find_cubit_bin
_cubit_path = find_cubit_bin()
if _cubit_path: sys.path.append(_cubit_path)
os.environ['CUBIT_PLUGIN_DIR'] = r'C:\Program Files\Coreform Cubit 2025.3\bin\plugins'

# NGSolve FIRST
import netgen.meshing
from ngsolve import Mesh, Integrate, CF, BND

from netgen.csg import CSGeometry, Sphere, Pnt

import cubit
cubit.init(['cubit', '-nojournal', '-batch',
            '-commandplugindir', os.environ['CUBIT_PLUGIN_DIR']])
from cubit_mesh_export import extract_curved_mesh

RADIUS = 0.05
V_EXACT = (4.0 / 3.0) * math.pi * RADIUS ** 3
A_EXACT = 4.0 * math.pi * RADIUS ** 2
OUT_DIR = os.path.join(_test_dir, 'export_netgen_test')
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 75)
print("Three-path comparison: sphere R=0.05m")
print("=" * 75)
print(f"V_exact = {V_EXACT:.10e} m^3")
print(f"A_exact = {A_EXACT:.10e} m^2")

# ================================================================
# Create mesh in Cubit
# ================================================================
cubit.cmd("reset")
cubit.cmd(f"create sphere radius {RADIUS}")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size auto factor 5")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add volume 1")
cubit.cmd('block 1 name "sphere"')
print(f"\nCubit mesh: {cubit.get_tet_count()} tets, {cubit.get_tri_count()} surface tris")

# Save cub5 for potential reuse
cub5_path = os.path.join(OUT_DIR, "sphere_debug.cub5")
cubit.cmd(f'save cub5 "{cub5_path}" overwrite')

# ================================================================
# Results table
# ================================================================
header = f"{'Path':<40} {'Volume':>14} {'V_err%':>12} {'Area':>14} {'A_err%':>12} {'ne':>7} {'nse':>7} {'nseg':>5} {'nedges':>7} {'n_ho_e':>7}"

def measure(mesh_or_path, label):
    """Measure volume and area from mesh or .vol path."""
    if isinstance(mesh_or_path, str):
        m = Mesh(mesh_or_path)
    else:
        m = mesh_or_path
    vol = Integrate(CF(1), m)
    area = Integrate(CF(1), m, BND)
    v_err = (vol - V_EXACT) / V_EXACT * 100
    a_err = (area - A_EXACT) / A_EXACT * 100
    return vol, v_err, area, a_err, m.ne


def count_vol_structure(path):
    """Count structural elements in .vol file."""
    info = {'ne': 0, 'nse': 0, 'nseg': 0, 'nedges': 0, 'n_ho_e': 0}
    with open(path) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        s = line.strip()
        if s == 'volumeelements':
            info['ne'] = int(lines[i + 1].strip())
        elif s == 'surfaceelements':
            info['nse'] = int(lines[i + 1].strip())
        elif s == 'edgesegmentsgi2':
            info['nseg'] = int(lines[i + 1].strip())
        elif s == 'curvedelements':
            nedges = int(lines[i + 1].strip())
            info['nedges'] = nedges
            nho = sum(1 for j in range(i + 2, i + 2 + nedges)
                      if int(lines[j].strip()) > 1)
            info['n_ho_e'] = nho
    return info


def print_row(label, vol, v_err, area, a_err, ne, si):
    print(f"{label:<40} {vol:>14.8e} {v_err:>+12.6e} {area:>14.8e} {a_err:>+12.6e} {ne:>7} {si.get('nse','-'):>7} {si.get('nseg','-'):>5} {si.get('nedges','-'):>7} {si.get('n_ho_e','-'):>7}")


for order in [2, 3]:
    print(f"\n{'='*75}")
    print(f"Order {order}")
    print(f"{'='*75}")
    print(header)
    print("-" * len(header))

    # --- Path A: export netgen (compact_netgen) ---
    path_a = os.path.join(OUT_DIR, f"3path_A_o{order}.vol")
    cubit.cmd(f'export netgen "{path_a}" order {order} overwrite')
    vol_a, verr_a, area_a, aerr_a, ne_a = measure(path_a, "A")
    si_a = count_vol_structure(path_a)
    print_row(f"A: export netgen (compact)", vol_a, verr_a, area_a, aerr_a, ne_a, si_a)

    # --- Path B: extract_curved_mesh -> Save ---
    ng_mesh_b = extract_curved_mesh(cubit, order=order)
    path_b = os.path.join(OUT_DIR, f"3path_B_o{order}.vol")
    ng_mesh_b.Save(path_b)

    # In-memory
    mesh_b_mem = Mesh(ng_mesh_b)
    vol_b, verr_b, area_b, aerr_b, ne_b = (
        Integrate(CF(1), mesh_b_mem),
        (Integrate(CF(1), mesh_b_mem) - V_EXACT) / V_EXACT * 100,
        Integrate(CF(1), mesh_b_mem, BND),
        (Integrate(CF(1), mesh_b_mem, BND) - A_EXACT) / A_EXACT * 100,
        mesh_b_mem.ne)
    si_b_mem = {'nse': '-', 'nseg': '-', 'nedges': '-', 'n_ho_e': '-'}
    print_row(f"B: extract_curved (in-mem)", vol_b, verr_b, area_b, aerr_b, ne_b, si_b_mem)

    # Reloaded
    vol_br, verr_br, area_br, aerr_br, ne_br = measure(path_b, "B reload")
    si_b = count_vol_structure(path_b)
    print_row(f"B: extract_curved (reload)", vol_br, verr_br, area_br, aerr_br, ne_br, si_b)

    # --- Path C: STEP -> netgen.occ -> Curve ---
    try:
        geo_c = CSGeometry()
        geo_c.Add(Sphere(Pnt(0, 0, 0), RADIUS).bc('surface'))
        ng_mesh_c = geo_c.GenerateMesh(maxh=0.0058)  # ~10100 tets to match Cubit ~10359
        ng_mesh_c.Curve(order)
        path_c = os.path.join(OUT_DIR, f"3path_C_o{order}.vol")
        ng_mesh_c.Save(path_c)

        mesh_c = Mesh(ng_mesh_c)
        vol_c = Integrate(CF(1), mesh_c)
        area_c = Integrate(CF(1), mesh_c, BND)
        verr_c = (vol_c - V_EXACT) / V_EXACT * 100
        aerr_c = (area_c - A_EXACT) / A_EXACT * 100
        ne_c = mesh_c.ne
        si_c = count_vol_structure(path_c)
        print_row(f"C: CSG->Netgen->Curve({order})", vol_c, verr_c, area_c, aerr_c, ne_c, si_c)

        # Also reload
        vol_cr, verr_cr, area_cr, aerr_cr, ne_cr = measure(path_c, "C reload")
        print_row(f"C: CSG->Netgen (reload)", vol_cr, verr_cr, area_cr, aerr_cr, ne_cr, si_c)
    except Exception as e:
        import traceback
        print(f"C: CSG->Netgen FAILED: {e}")
        traceback.print_exc()

    # --- Exact ---
    print("-" * len(header))
    print(f"{'Exact':40} {V_EXACT:>14.8e} {'0':>12} {A_EXACT:>14.8e} {'0':>12}")

print(f"\n{'='*75}")
print("Root cause summary:")
print("  Path A: get_surface_tris() returns UNSHARED node IDs (rendering mesh)")
print("         -> surface elements don't share nodes with volume elements")
print("         -> topology differs, fewer edges curved -> volume error")
print("  Path B: parse_cubit_list('tri') + get_connectivity returns SHARED nodes")
print("         -> correct topology, proper curving")
print("  Path C: CSG sphere + Netgen mesher (exact geometry, different mesh)")
print("         -> volume/area reference (independent validation)")
