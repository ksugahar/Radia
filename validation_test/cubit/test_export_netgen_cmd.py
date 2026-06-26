"""
Test 'export netgen' APREPRO command (compact Netgen, static link).

Tests:
1. Sphere tet mesh (order 1-4)
2. Brick hex mesh (order 2)
3. Multi-block with labels
4. Sideset boundary labels
5. Volume accuracy (sphere, order 1 vs 2)
6. NGSolve readback verification

Usage:
  python validation_test/cubit/test_export_netgen_cmd.py
"""

import sys
import os
import math
import tempfile

# --- Setup paths ---
_test_dir = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.dirname(os.path.dirname(_test_dir))
sys.path.insert(0, os.path.join(_repo_root, 'src', 'radia'))
from install_panels import find_cubit_bin

_cubit_path = find_cubit_bin()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)

_plugin_dir = os.path.join(os.path.dirname(_cubit_path), 'bin', 'plugins') \
    if _cubit_path else None
if _plugin_dir and os.path.isdir(_plugin_dir):
    os.environ['CUBIT_PLUGIN_DIR'] = _plugin_dir

import cubit
cubit.init(['cubit', '-nojournal', '-batch',
            '-commandplugindir', _plugin_dir or ''])

# Check if NGSolve is available
try:
    import netgen.meshing
    from ngsolve import Mesh, Integrate, CF
    NGSOLVE_AVAILABLE = True
    print("NGSolve available")
except ImportError:
    NGSOLVE_AVAILABLE = False
    print("NGSolve not available - file-only tests")

OUT_DIR = os.path.join(_test_dir, 'export_netgen_test')
os.makedirs(OUT_DIR, exist_ok=True)

RADIUS = 0.05
V_EXACT = (4.0 / 3.0) * math.pi * RADIUS ** 3
passed = 0
failed = 0


def check_vol_file(path):
    """Check .vol file exists and has reasonable content."""
    assert os.path.exists(path), f"File not found: {path}"
    sz = os.path.getsize(path)
    assert sz > 100, f"File too small: {sz} bytes"
    # Check header — Netgen .vol starts with optional comment then 'mesh3d'
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            assert line == 'mesh3d', f"Bad header: {line!r}"
            break
    return sz


def run_test(name, func):
    global passed, failed
    print(f"\n{'=' * 60}")
    print(f"Test: {name}")
    print('=' * 60)
    try:
        func()
        print(f"  PASS")
        passed += 1
    except Exception as e:
        print(f"  FAIL: {e}")
        import traceback
        traceback.print_exc()
        failed += 1


# ================================================================
# Test 1: Sphere tet, order 1-4
# ================================================================
def test_sphere_tet_orders():
    cubit.cmd("reset")
    cubit.cmd(f"create sphere radius {RADIUS}")
    cubit.cmd("volume 1 scheme tetmesh")
    cubit.cmd(f"volume 1 size auto factor 5")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add volume 1")
    cubit.cmd('block 1 name "sphere"')

    for order in [1, 2, 3, 4]:
        vol_path = os.path.join(OUT_DIR, f"sphere_o{order}.vol")
        cubit.cmd(f'export netgen "{vol_path}" order {order} overwrite')
        sz = check_vol_file(vol_path)
        print(f"  order {order}: {sz:,} bytes")


# ================================================================
# Test 2: Brick hex mesh, order 2
# ================================================================
def test_brick_hex():
    cubit.cmd("reset")
    cubit.cmd("create brick x 0.1 y 0.1 z 0.1")
    cubit.cmd("volume 1 scheme map")
    cubit.cmd("volume 1 size 0.025")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add volume 1")
    cubit.cmd('block 1 name "brick"')

    vol_path = os.path.join(OUT_DIR, "brick_hex_o2.vol")
    cubit.cmd(f'export netgen "{vol_path}" order 2 overwrite')
    sz = check_vol_file(vol_path)
    print(f"  hex order 2: {sz:,} bytes")


# ================================================================
# Test 3: Multi-block with material labels
# ================================================================
def test_multi_block_labels():
    cubit.cmd("reset")
    cubit.cmd("create brick x 0.1 y 0.1 z 0.1")
    cubit.cmd("volume 1 move 0 0 0")
    cubit.cmd("create brick x 0.1 y 0.1 z 0.1")
    cubit.cmd("volume 2 move 0.15 0 0")
    cubit.cmd("volume all scheme tetmesh")
    cubit.cmd("volume all size auto factor 5")
    cubit.cmd("mesh volume all")
    cubit.cmd("block 1 add volume 1")
    cubit.cmd('block 1 name "iron"')
    cubit.cmd("block 2 add volume 2")
    cubit.cmd('block 2 name "air"')

    vol_path = os.path.join(OUT_DIR, "multi_block.vol")
    cubit.cmd(f'export netgen "{vol_path}" order 2 overwrite')
    sz = check_vol_file(vol_path)
    print(f"  multi-block: {sz:,} bytes")

    if NGSOLVE_AVAILABLE:
        mesh = Mesh(vol_path)
        mats = set(mesh.GetMaterials())
        print(f"  Materials: {mats}")
        assert 'iron' in mats or 'volume_1' in mats, f"Expected 'iron' in materials, got {mats}"


# ================================================================
# Test 4: Sideset boundary labels
# ================================================================
def test_sideset_labels():
    cubit.cmd("reset")
    cubit.cmd("create brick x 0.1 y 0.1 z 0.1")
    cubit.cmd("volume 1 scheme tetmesh")
    cubit.cmd("volume 1 size auto factor 5")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add volume 1")
    cubit.cmd('block 1 name "solid"')
    cubit.cmd("sideset 1 add surface 1")
    cubit.cmd('sideset 1 name "bottom"')
    cubit.cmd("sideset 2 add surface 2")
    cubit.cmd('sideset 2 name "top"')

    vol_path = os.path.join(OUT_DIR, "sideset_labels.vol")
    cubit.cmd(f'export netgen "{vol_path}" order 2 overwrite')
    sz = check_vol_file(vol_path)
    print(f"  sideset: {sz:,} bytes")

    if NGSOLVE_AVAILABLE:
        mesh = Mesh(vol_path)
        bnds = set(mesh.GetBoundaries())
        print(f"  Boundaries: {bnds}")
        # Check that at least some named boundaries exist
        assert len(bnds) > 0, "No boundaries found"


# ================================================================
# Test 5: Volume accuracy (sphere order 1 vs 2)
# ================================================================
def test_volume_accuracy():
    if not NGSOLVE_AVAILABLE:
        print("  Skipping (NGSolve required)")
        return

    cubit.cmd("reset")
    cubit.cmd(f"create sphere radius {RADIUS}")
    cubit.cmd("volume 1 scheme tetmesh")
    cubit.cmd(f"volume 1 size auto factor 5")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add volume 1")
    cubit.cmd('block 1 name "sphere"')

    results = {}
    for order in [1, 2]:
        vol_path = os.path.join(OUT_DIR, f"sphere_vol_o{order}.vol")
        cubit.cmd(f'export netgen "{vol_path}" order {order} overwrite')

        mesh = Mesh(vol_path)
        vol = Integrate(CF(1), mesh)
        err = (vol - V_EXACT) / V_EXACT * 100
        results[order] = (vol, err)
        print(f"  order {order}: V={vol:.6e}, err={err:+.2e}%")

    # Order 2 should be more accurate than order 1
    assert abs(results[2][1]) < abs(results[1][1]), \
        f"Order 2 ({results[2][1]:.2e}%) should be more accurate than order 1 ({results[1][1]:.2e}%)"
    print(f"  Volume improvement: order 1 -> order 2")


# ================================================================
# Test 6: NGSolve readback + FEM solve
# ================================================================
def test_ngsolve_readback():
    if not NGSOLVE_AVAILABLE:
        print("  Skipping (NGSolve required)")
        return

    from ngsolve import H1, GridFunction, BilinearForm, LinearForm, grad, dx
    from ngsolve import TaskManager

    cubit.cmd("reset")
    cubit.cmd("create brick x 0.1 y 0.1 z 0.1")
    cubit.cmd("volume 1 scheme tetmesh")
    cubit.cmd("volume 1 size auto factor 5")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add volume 1")
    cubit.cmd('block 1 name "solid"')

    vol_path = os.path.join(OUT_DIR, "fem_test.vol")
    cubit.cmd(f'export netgen "{vol_path}" order 2 overwrite')

    mesh = Mesh(vol_path)
    print(f"  Mesh: {mesh.ne} elements, {mesh.nv} vertices")

    fes = H1(mesh, order=2, dirichlet=".*")
    print(f"  H1 space: {fes.ndof} DOFs")

    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += grad(u) * grad(v) * dx
    with TaskManager():
        a.Assemble()

        f = LinearForm(fes)
        f += 1 * v * dx
        f.Assemble()

        gfu = GridFunction(fes)
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec
        print(f"  Poisson solve OK, max(u) = {max(gfu.vec):.6e}")


# ================================================================
# Main
# ================================================================
if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("export netgen Command Test Suite (compact Netgen)")
    print("=" * 60)

    run_test("Sphere tet orders 1-4", test_sphere_tet_orders)
    run_test("Brick hex order 2", test_brick_hex)
    run_test("Multi-block labels", test_multi_block_labels)
    run_test("Sideset boundary labels", test_sideset_labels)
    run_test("Volume accuracy (sphere)", test_volume_accuracy)
    run_test("NGSolve readback + FEM", test_ngsolve_readback)

    print(f"\n{'=' * 60}")
    print(f"Results: {passed} passed, {failed} failed")
    print('=' * 60)
    sys.exit(1 if failed > 0 else 0)
