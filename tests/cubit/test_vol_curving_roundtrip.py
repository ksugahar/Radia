"""
Test: .vol curving round-trip (Save + Load curvedelements).

Workflow:
  1. Cubit: sphere -> tet mesh
  2. extract_curved_mesh(order=2) -> NGSolve Mesh (in-memory, curved)
  3. Save .vol (includes curvedelements section)
  4. Reload .vol -> NGSolve Mesh
  5. Compare volume accuracy: in-memory vs reloaded

This tests the Netgen fork's curvedelements Load parser (meshclass.cpp patch).

Run inside Cubit Python or with system Python (needs Cubit + NGSolve):
  python tests/cubit/test_vol_curving_roundtrip.py
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

# Then Cubit
import cubit
cubit.init(['cubit', '-nojournal', '-batch'])

# cubit_mesh_export (extract_curved_mesh)
from cubit_mesh_export import extract_curved_mesh

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
# Step 2: extract_curved_mesh for each order, save .vol, reload
# ============================================================
results = {}

for order in [2, 3]:
    print(f"\n--- Order {order} ---")

    # Get curved mesh in memory
    ng_mesh = extract_curved_mesh(cubit, order=order)
    mesh_mem = Mesh(ng_mesh)
    vol_mem = Integrate(CF(1), mesh_mem)
    err_mem = (vol_mem - V_EXACT) / V_EXACT * 100
    print(f"  In-memory:  V={vol_mem:.10e}  err={err_mem:+.6e}%")

    # Save to .vol
    vol_path = os.path.join(OUT_DIR, f"sphere_roundtrip_o{order}.vol")
    ng_mesh.Save(vol_path)
    sz = os.path.getsize(vol_path)
    print(f"  Saved: {vol_path} ({sz:,} bytes)")

    # Check if curvedelements section exists
    has_curved = False
    with open(vol_path, 'r') as f:
        for line in f:
            if line.strip() == 'curvedelements':
                has_curved = True
                break
    print(f"  curvedelements section: {'YES' if has_curved else 'NO'}")

    # Reload from .vol
    mesh_reload = Mesh(vol_path)
    vol_reload = Integrate(CF(1), mesh_reload)
    err_reload = (vol_reload - V_EXACT) / V_EXACT * 100
    print(f"  Reloaded:   V={vol_reload:.10e}  err={err_reload:+.6e}%")

    # Compare
    diff = abs(vol_mem - vol_reload)
    print(f"  Difference: {diff:.6e} ({diff/V_EXACT*100:.6e}%)")

    results[order] = {
        'vol_mem': vol_mem, 'err_mem': err_mem,
        'vol_reload': vol_reload, 'err_reload': err_reload,
        'has_curved': has_curved, 'diff': diff,
    }

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 60)
print("Summary")
print("=" * 60)
print(f"{'Order':<6} {'In-memory err%':<18} {'Reloaded err%':<18} {'Match?':<8} {'curvedelems':<12}")
for order, r in results.items():
    match = abs(r['diff']) < 1e-6 * V_EXACT  # float64 rounding
    print(f"{order:<6} {r['err_mem']:>+.6e}%    {r['err_reload']:>+.6e}%    {'OK' if match else 'FAIL':<8} {'YES' if r['has_curved'] else 'NO':<12}")

# Verify order 3 is better than order 2
if 2 in results and 3 in results:
    if abs(results[3]['err_reload']) < abs(results[2]['err_reload']):
        print("\nOrder 3 reloaded is more accurate than order 2: PASS")
    else:
        print(f"\nWARNING: Order 3 reloaded ({results[3]['err_reload']:+.6e}%) "
              f"is NOT more accurate than order 2 ({results[2]['err_reload']:+.6e}%)")

# Verify round-trip preserves curving
all_match = all(abs(r['diff']) < 1e-6 * V_EXACT for r in results.values())  # float64 rounding
if all_match:
    print("Round-trip curving preservation: PASS")
else:
    print("Round-trip curving preservation: FAIL")
    for order, r in results.items():
        if abs(r['diff']) > 1e-6 * V_EXACT:
            print(f"  Order {order}: diff={r['diff']:.6e}")
