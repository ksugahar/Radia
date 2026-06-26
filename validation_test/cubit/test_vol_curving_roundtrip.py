"""
Test: .vol curving round-trip (Save + Load curvedelements).

Workflow:
  1. Cubit: sphere -> tet mesh
  2. export netgen .vol with order=2,3 curving
  3. Load .vol -> NGSolve Mesh
  4. Check volume accuracy and curvedelements section presence

This tests the Netgen fork's curvedelements Load parser (meshclass.cpp patch).

Run inside Cubit Python or with system Python (needs Cubit + NGSolve):
  python validation_test/cubit/test_vol_curving_roundtrip.py
"""
import sys
import os
import math

# --- Setup paths ---
_test_dir = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.dirname(os.path.dirname(_test_dir))
sys.path.insert(0, os.path.join(_repo_root, 'src', 'radia'))
from install_panels import find_cubit_bin

_cubit_path = find_cubit_bin()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)

# Import NGSolve FIRST (DLL conflict avoidance)
import netgen.meshing
from ngsolve import Mesh, Integrate, CF
from ngsolve import TaskManager

# Then Cubit
_plugin_dir = os.path.join(os.path.dirname(_cubit_path), 'bin', 'plugins') \
    if _cubit_path else None
if _plugin_dir and os.path.isdir(_plugin_dir):
    os.environ['CUBIT_PLUGIN_DIR'] = _plugin_dir

import cubit
cubit.init(['cubit', '-nojournal', '-batch',
            '-commandplugindir', _plugin_dir or ''])

# ============================================================
# Parameters
# ============================================================
RADIUS = 0.05  # meters
V_EXACT = (4.0 / 3.0) * math.pi * RADIUS ** 3
OUT_DIR = os.path.join(_test_dir, 'export_netgen_test')
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print("Test: .vol curving round-trip")
print("=" * 60)
print(f"Sphere R={RADIUS} m, V_exact={V_EXACT:.10e} m^3")
print()

# ============================================================
# Step 1: Create sphere mesh in Cubit
# ============================================================
cubit.cmd("reset")
cubit.cmd(f"create sphere radius {RADIUS}")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size auto factor 5")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add volume 1")
cubit.cmd('block 1 name "sphere"')
print(f"Cubit mesh: {cubit.get_tet_count()} tets")

# ============================================================
# Step 2: export netgen .vol for each order, load and check
# ============================================================
results = {}

for order in [2, 3]:
    print(f"\n--- Order {order} ---")

    # Export .vol via C++ plugin
    vol_path = os.path.join(OUT_DIR, f"sphere_roundtrip_o{order}.vol")
    cubit.cmd(f'export netgen "{vol_path}" order {order} overwrite')
    sz = os.path.getsize(vol_path)
    print(f"  Exported: {vol_path} ({sz:,} bytes)")

    # Check if curvedelements section exists
    has_curved = False
    with open(vol_path, 'r') as f:
        for line in f:
            if line.strip() == 'curvedelements':
                has_curved = True
                break
    print(f"  curvedelements section: {'YES' if has_curved else 'NO'}")

    # Load from .vol
    mesh_reload = Mesh(vol_path)
    with TaskManager():
        vol_reload = Integrate(CF(1), mesh_reload)
        err_reload = (vol_reload - V_EXACT) / V_EXACT * 100
        curve_order = mesh_reload.GetCurveOrder()
        print(f"  Loaded:     V={vol_reload:.10e}  err={err_reload:+.6e}%")
        print(f"  CurveOrder: {curve_order} (expected {order})")

        results[order] = {
            'vol_reload': vol_reload, 'err_reload': err_reload,
            'has_curved': has_curved,
            'curve_order': curve_order,
        }

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 60)
print("Summary")
print("=" * 60)
print(f"{'Order':<6} {'Loaded err%':<18} {'curvedelems':<12}")
for order, r in results.items():
    print(f"{order:<6} {r['err_reload']:>+.6e}%    {'YES' if r['has_curved'] else 'NO':<12}")

# Verify order 3 is better than order 2
if 2 in results and 3 in results:
    if abs(results[3]['err_reload']) < abs(results[2]['err_reload']):
        print("\nOrder 3 is more accurate than order 2: PASS")
    else:
        print(f"\nWARNING: Order 3 ({results[3]['err_reload']:+.6e}%) "
              f"is NOT more accurate than order 2 ({results[2]['err_reload']:+.6e}%)")

# Verify curvedelements section is present
all_curved = all(r['has_curved'] for r in results.values())
if all_curved:
    print("curvedelements section present in all .vol files: PASS")
else:
    print("curvedelements section MISSING in some .vol files: FAIL")

# Verify CurveOrder is actually loaded (requires Netgen fork with curvedelements parser)
curve_ok = all(r['curve_order'] == order for order, r in results.items())
if curve_ok:
    print("CurveOrder matches export order for all .vol files: PASS")
else:
    for order, r in results.items():
        if r['curve_order'] != order:
            print(f"WARNING: Order {order} .vol loaded as CurveOrder={r['curve_order']} "
                  f"(expected {order}). Netgen fork with curvedelements Load parser required.")
