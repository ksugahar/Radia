"""
High-order volume convergence test for mixed element types.

Tests:
  1. p-convergence: export netgen .vol at order 1-5, volume must converge to CAD
  2. gmsh cross-check: export .msh order 1 and 2, compare volume vs netgen .vol
  3. Path A==B: export netgen (C++) vs extract_curved_mesh (Python) volume match

Cases:
  1. Flat brick (hex + tet + pyramid) -- exact volume expected
  2. Cylinder hex (scheme map) -- curved surface
  3. Cylinder webcut hex (multi-volume) -- curved, 4 volumes
  4. Cylinder boundary-layer hex -- curved, BL mesh
  5. Loft circle-to-rectangle hex -- complex shape

Usage:
  python tests/cubit/test_ho_volume_convergence.py
"""
import sys
import os
import math
import numpy as np

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

# Import NGSolve FIRST (DLL conflict avoidance)
import netgen.meshing
from ngsolve import Mesh, Integrate, CF

import cubit
cubit.init(['cubit', '-nojournal', '-batch',
            '-commandplugindir', _plugin_dir or ''])
from cubit_mesh_export import extract_curved_mesh

OUT_DIR = os.path.join(_test_dir, 'ho_convergence_test')
os.makedirs(OUT_DIR, exist_ok=True)


# ================================================================
# Test case definitions
# ================================================================
test_cases = [
    ("case01_flat_brick", [
        "reset",
        "brick x 2 y 1 z 1",
        "webcut volume 1 with plane xplane imprint merge",
        "volume 1 scheme map", "volume 1 size 1", "mesh volume 1",
        "volume 2 scheme tetmesh", "volume 2 size 5", "mesh volume 2",
        "block 1 add hex all", 'block 1 name "map"',
        "block 2 add tet all", 'block 2 name "tet"',
        "block 3 add pyramid all", 'block 3 name "pyram"',
        "Volume all scale 0.001",
    ], True),   # is_flat=True (no curved surfaces)
    ("case02_cylinder_map", [
        "reset",
        "cylinder height 1 radius 1",
        "volume 1 size 0.5", "mesh volume 1",
        "block 1 add hex all", 'block 1 name "cylinder_hex"',
        "Volume all scale 0.001",
    ], False),
    ("case03_cylinder_webcut", [
        "reset",
        "cylinder height 1 radius 1",
        "webcut volume all with plane xplane offset 0 imprint merge",
        "webcut volume all with plane yplane offset 0 imprint merge",
        "imprint all", "merge all",
        "volume all size 1", "mesh volume all",
        "block 1 add hex all", 'block 1 name "cylinder_hex"',
        "Volume all scale 0.001",
    ], False),
    ("case04_cylinder_BL", [
        "reset",
        "cylinder height 1 radius 1",
        "create boundary_layer 1",
        "modify boundary_layer 1 uniform height 0.03 growth 1.2 layers 3",
        "modify boundary_layer 1 add surface 1 volume 1",
        "modify boundary_layer 1 continuity on",
        "volume 1 size 0.6", "mesh volume 1",
        "block 1 add hex all", 'block 1 name "cylinder_hex"',
        "Volume all scale 0.001",
    ], False),
    ("case05_loft", [
        "reset",
        "create surface rectangle width 1 height 1 zplane",
        "create surface circle radius 1 zplane",
        "move Surface 1 z 1 include_merged",
        "create volume loft surface 2 1",
        "delete body 1 2", "compress",
        "volume 1 size auto factor 9", "mesh volume 1",
        "block 1 add hex all", 'block 1 name "map"',
        "Volume all scale 0.001",
    ], False),
]


# ================================================================
# HEX20 serendipity volume quadrature
# ================================================================
def hex20_volume(coords):
    """Volume of HEX20 using 3x3x3 Gauss quadrature."""
    gp = math.sqrt(3.0 / 5.0)
    gw = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]
    gx = [-gp, 0.0, gp]

    # Corner reference coordinates [-1,1]^3
    xc = np.array([-1, 1, 1, -1, -1, 1, 1, -1], dtype=float)
    yc = np.array([-1, -1, 1, 1, -1, -1, 1, 1], dtype=float)
    zc = np.array([-1, -1, -1, -1, 1, 1, 1, 1], dtype=float)

    vol = 0.0
    for i in range(3):
        for j in range(3):
            for k in range(3):
                xi, eta, zeta = gx[i], gx[j], gx[k]
                w = gw[i] * gw[j] * gw[k]

                dNdxi = np.zeros(20)
                dNdeta = np.zeros(20)
                dNdzeta = np.zeros(20)

                # Corner nodes (0-7): serendipity formula
                for n in range(8):
                    x0, y0, z0 = xc[n], yc[n], zc[n]
                    a = 1 + x0 * xi
                    b = 1 + y0 * eta
                    c = 1 + z0 * zeta
                    s = x0 * xi + y0 * eta + z0 * zeta - 2
                    dNdxi[n] = (1 / 8.0) * x0 * b * c * s + (1 / 8.0) * a * b * c * x0
                    dNdeta[n] = (1 / 8.0) * a * y0 * c * s + (1 / 8.0) * a * b * c * y0
                    dNdzeta[n] = (1 / 8.0) * a * b * z0 * s + (1 / 8.0) * a * b * c * z0

                # Mid-edge nodes along xi (8, 10, 12, 14)
                for idx, e0, e1 in [(8, -1, -1), (10, 1, -1),
                                     (12, -1, 1), (14, 1, 1)]:
                    aa = 1 - xi * xi
                    bb = 1 + e0 * eta
                    cc = 1 + e1 * zeta
                    dNdxi[idx] = (1 / 4.0) * (-2 * xi) * bb * cc
                    dNdeta[idx] = (1 / 4.0) * aa * e0 * cc
                    dNdzeta[idx] = (1 / 4.0) * aa * bb * e1

                # Mid-edge nodes along eta (9, 11, 13, 15)
                for idx, x0, z0 in [(9, 1, -1), (11, -1, -1),
                                     (13, 1, 1), (15, -1, 1)]:
                    aa = 1 + x0 * xi
                    bb = 1 - eta * eta
                    cc = 1 + z0 * zeta
                    dNdxi[idx] = (1 / 4.0) * x0 * bb * cc
                    dNdeta[idx] = (1 / 4.0) * aa * (-2 * eta) * cc
                    dNdzeta[idx] = (1 / 4.0) * aa * bb * z0

                # Mid-edge nodes along zeta (16, 17, 18, 19)
                for idx, x0, y0 in [(16, -1, -1), (17, 1, -1),
                                     (18, 1, 1), (19, -1, 1)]:
                    aa = 1 + x0 * xi
                    bb = 1 + y0 * eta
                    cc = 1 - zeta * zeta
                    dNdxi[idx] = (1 / 4.0) * x0 * bb * cc
                    dNdeta[idx] = (1 / 4.0) * aa * y0 * cc
                    dNdzeta[idx] = (1 / 4.0) * aa * bb * (-2 * zeta)

                J = np.zeros((3, 3))
                for n in range(20):
                    J[0] += dNdxi[n] * coords[n]
                    J[1] += dNdeta[n] * coords[n]
                    J[2] += dNdzeta[n] * coords[n]

                vol += w * abs(np.linalg.det(J))
    return vol


# ================================================================
# TET10 volume quadrature
# ================================================================
def tet10_volume(coords):
    """Volume of TET10 using 4-point Gauss quadrature."""
    a = 0.1381966011250105
    b = 0.5854101966249685
    gauss_pts = [(a, a, a, b), (a, a, b, a), (a, b, a, a), (b, a, a, a)]

    vol = 0.0
    for L1, L2, L3, L4 in gauss_pts:
        dNdL1 = np.array([
            4 * L1 - 1, 0, 0, -(4 * L4 - 1),
            4 * L2, 0, 4 * L3, 4 * (L4 - L1), -4 * L2, -4 * L3
        ])
        dNdL2 = np.array([
            0, 4 * L2 - 1, 0, -(4 * L4 - 1),
            4 * L1, 4 * L3, 0, -4 * L1, 4 * (L4 - L2), -4 * L3
        ])
        dNdL3 = np.array([
            0, 0, 4 * L3 - 1, -(4 * L4 - 1),
            0, 4 * L2, 4 * L1, -4 * L1, -4 * L2, 4 * (L4 - L3)
        ])
        J = np.array([
            sum(dNdL1[k] * coords[k] for k in range(10)),
            sum(dNdL2[k] * coords[k] for k in range(10)),
            sum(dNdL3[k] * coords[k] for k in range(10)),
        ])
        vol += abs(np.linalg.det(J)) / 24.0
    return vol


# ================================================================
# Simple volume routines (linear elements + PYRAMID13 fallback)
# ================================================================
def tet4_volume(pts):
    v1, v2, v3 = pts[1] - pts[0], pts[2] - pts[0], pts[3] - pts[0]
    return abs(np.dot(v1, np.cross(v2, v3))) / 6.0


def hex8_volume(pts):
    gp = 1.0 / math.sqrt(3.0)
    gauss = [(-gp, -gp, -gp), (gp, -gp, -gp), (-gp, gp, -gp), (gp, gp, -gp),
             (-gp, -gp, gp), (gp, -gp, gp), (-gp, gp, gp), (gp, gp, gp)]
    vol = 0.0
    for xi, eta, zeta in gauss:
        dNdxi = np.array([
            -(1 - eta) * (1 - zeta), (1 - eta) * (1 - zeta),
            (1 + eta) * (1 - zeta), -(1 + eta) * (1 - zeta),
            -(1 - eta) * (1 + zeta), (1 - eta) * (1 + zeta),
            (1 + eta) * (1 + zeta), -(1 + eta) * (1 + zeta)]) / 8.0
        dNdeta = np.array([
            -(1 - xi) * (1 - zeta), -(1 + xi) * (1 - zeta),
            (1 + xi) * (1 - zeta), (1 - xi) * (1 - zeta),
            -(1 - xi) * (1 + zeta), -(1 + xi) * (1 + zeta),
            (1 + xi) * (1 + zeta), (1 - xi) * (1 + zeta)]) / 8.0
        dNdzeta = np.array([
            -(1 - xi) * (1 - eta), -(1 + xi) * (1 - eta),
            -(1 + xi) * (1 + eta), -(1 - xi) * (1 + eta),
            (1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
            (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)]) / 8.0
        J = np.zeros((3, 3))
        for n in range(8):
            J[0] += dNdxi[n] * pts[n]
            J[1] += dNdeta[n] * pts[n]
            J[2] += dNdzeta[n] * pts[n]
        vol += abs(np.linalg.det(J))
    return vol


def pyramid5_volume(pts):
    v1 = tet4_volume(np.array([pts[0], pts[1], pts[2], pts[4]]))
    v2 = tet4_volume(np.array([pts[0], pts[2], pts[3], pts[4]]))
    return v1 + v2


# ================================================================
# GMSH v2.2 parser
# ================================================================
def parse_gmsh_v22(filename):
    """Parse GMSH v2.2 and compute total volume from all 3D elements."""
    with open(filename, 'r') as f:
        content = f.read()

    nodes = {}
    elements = []
    lines = content.split('\n')
    section = None
    n_nodes = n_elems = 0

    for line in lines:
        line = line.strip()
        if line.startswith('$'):
            if line == '$Nodes':
                section = 'nodes'
            elif line == '$Elements':
                section = 'elements'
            elif line.startswith('$End'):
                section = None
            continue
        if section == 'nodes':
            if n_nodes == 0:
                n_nodes = int(line)
            else:
                parts = line.split()
                if len(parts) >= 4:
                    nodes[int(parts[0])] = np.array(
                        [float(parts[1]), float(parts[2]), float(parts[3])])
        elif section == 'elements':
            if n_elems == 0:
                n_elems = int(line)
            else:
                parts = line.split()
                if len(parts) >= 4:
                    etype = int(parts[1])
                    ntags = int(parts[2])
                    conn = [int(p) for p in parts[3 + ntags:]]
                    elements.append((etype, conn))

    type_names = {
        4: 'TET4', 11: 'TET10', 5: 'HEX8', 17: 'HEX20',
        6: 'WEDGE6', 18: 'WEDGE15', 7: 'PYRAMID5', 19: 'PYRAMID13',
    }
    vol = 0.0
    counts = {}
    for etype, conn in elements:
        try:
            if etype == 4:
                vol += tet4_volume(np.array([nodes[c] for c in conn[:4]]))
            elif etype == 11:
                vol += tet10_volume([nodes[c] for c in conn[:10]])
            elif etype == 5:
                vol += hex8_volume(np.array([nodes[c] for c in conn[:8]]))
            elif etype == 17:
                vol += hex20_volume([nodes[c] for c in conn[:20]])
            elif etype == 7:
                vol += pyramid5_volume(np.array([nodes[c] for c in conn[:5]]))
            elif etype == 19:
                # PYRAMID13: use linear fallback (proper quadrature needs
                # diagonal mid-node not present in GMSH PYRAMID13)
                vol += pyramid5_volume(np.array([nodes[c] for c in conn[:5]]))
            else:
                continue
        except Exception:
            continue
        name = type_names.get(etype, f'type{etype}')
        counts[name] = counts.get(name, 0) + 1

    return {'volume': vol, 'n_nodes': len(nodes), 'counts': counts}


# ================================================================
# Test 1: p-convergence (netgen .vol, order 1-5)
# ================================================================
def test_p_convergence(case_name, cad_volume, is_flat):
    """Export .vol at orders 1-5, load in NGSolve, check volume convergence."""
    print(f"\n  p-convergence test (CAD={cad_volume:.6e}):")
    print(f"  {'Order':>5} {'Volume':>14} {'Error':>14} {'Verdict':>8}")
    print(f"  {'-----':>5} {'-' * 14:>14} {'-' * 14:>14} {'--------':>8}")

    errors = []
    all_pass = True

    for order in range(1, 6):
        vol_path = os.path.join(OUT_DIR, f"{case_name}_o{order}.vol")
        cmd = f'export netgen "{vol_path}" order {order} overwrite'
        try:
            cubit.cmd(cmd)
            mesh = Mesh(vol_path)
            vol = Integrate(CF(1), mesh)
        except Exception as e:
            print(f"  {order:>5} FAILED: {e}")
            errors.append(None)
            all_pass = False
            continue

        err_pct = (vol - cad_volume) / cad_volume * 100.0
        verdict = "OK"

        if is_flat:
            # Flat geometry: volume must be exact at all orders
            if abs(err_pct) > 0.01:
                verdict = "FAIL"
                all_pass = False
        elif len(errors) > 0 and errors[-1] is not None:
            prev_err = errors[-1]
            # Allow 10x regression (numerical noise at high order)
            if abs(err_pct) > abs(prev_err) * 10 and abs(err_pct) > 1e-3:
                verdict = "REGRESS"
                all_pass = False

        print(f"  {order:>5} {vol:>14.6e} {err_pct:>+13.4e}% {verdict:>8}")
        errors.append(err_pct)

    # Convergence check: order 2 must be better than order 1 (curved cases)
    if not is_flat and len(errors) >= 2:
        if errors[0] is not None and errors[1] is not None:
            if abs(errors[1]) >= abs(errors[0]):
                print(f"  order 2 not better than order 1: FAIL")
                all_pass = False

    return all_pass


# ================================================================
# Test 2: gmsh vs netgen volume at order 2
# ================================================================
def test_gmsh_vs_netgen(case_name, cad_volume, is_flat):
    """Compare .vol and .msh volume at order 1 and 2."""
    print(f"\n  gmsh vs netgen (CAD={cad_volume:.6e}):")

    results = {}
    for order in [1, 2]:
        msh_path = os.path.join(OUT_DIR, f"{case_name}_o{order}.msh")
        cmd = f'export gmsh "{msh_path}" order {order} version 2 dimension 3 overwrite'
        cubit.cmd(cmd)
        r = parse_gmsh_v22(msh_path)
        r['error_pct'] = (r['volume'] - cad_volume) / cad_volume * 100.0
        results[f'gmsh_o{order}'] = r

    # Load netgen vol order 2 (already exported in p-convergence test)
    vol_path = os.path.join(OUT_DIR, f"{case_name}_o2.vol")
    try:
        mesh = Mesh(vol_path)
        ng_vol = Integrate(CF(1), mesh)
        ng_err = (ng_vol - cad_volume) / cad_volume * 100.0
    except Exception:
        ng_vol = ng_err = None

    print(f"  {'Source':<20} {'Volume':>14} {'Error':>14} {'Elements':>30}")
    print(f"  {'-' * 20:<20} {'-' * 14:>14} {'-' * 14:>14} {'-' * 30:>30}")
    for key in ['gmsh_o1', 'gmsh_o2']:
        r = results[key]
        elem_str = ', '.join(f"{k}={v}" for k, v in r['counts'].items())
        print(f"  {key:<20} {r['volume']:>14.6e} {r['error_pct']:>+13.4e}% {elem_str:>30}")
    if ng_vol is not None:
        print(f"  {'netgen_o2':<20} {ng_vol:>14.6e} {ng_err:>+13.4e}%")

    all_pass = True
    g1 = results['gmsh_o1']
    g2 = results['gmsh_o2']

    # Check 1: gmsh o2 should improve over o1 (curved surfaces only)
    if not is_flat:
        if abs(g2['error_pct']) >= abs(g1['error_pct']):
            print(f"  gmsh o2 not better than o1: [FAIL]")
            all_pass = False
        else:
            print(f"  gmsh o2 better than o1: [PASS]")

    # Check 2: gmsh o2 and netgen o2 must agree within 1%
    if ng_vol is not None:
        diff_pct = abs(g2['volume'] - ng_vol) / cad_volume * 100.0
        cross_pass = diff_pct < 1.0
        print(f"  gmsh_o2 vs netgen_o2 diff: {diff_pct:.4e}% "
              f"[{'PASS' if cross_pass else 'FAIL'}]")
        if not cross_pass:
            all_pass = False

    # Check 3: flat geometry must be exact
    if is_flat and abs(g1['error_pct']) > 0.01:
        print(f"  flat geometry o1 not exact: [FAIL]")
        all_pass = False

    return all_pass


# ================================================================
# Test 3: Path A (export netgen C++) vs Path B (extract_curved_mesh Python)
# ================================================================
def test_path_ab(case_name, cad_volume):
    """Compare .vol from Path A and Path B at order=2."""
    print(f"\n  Path A vs B (order=2, CAD={cad_volume:.6e}):")

    ORDER = 2

    # Path A: export netgen (already tested, re-use .vol)
    path_a = os.path.join(OUT_DIR, f"{case_name}_o{ORDER}.vol")
    if not os.path.exists(path_a):
        cubit.cmd(f'export netgen "{path_a}" order {ORDER} overwrite')
    mesh_a = Mesh(path_a)
    vol_a = Integrate(CF(1), mesh_a)

    # Path B: extract_curved_mesh (Python)
    path_b = os.path.join(OUT_DIR, f"{case_name}_pathB.vol")
    try:
        ng_b = extract_curved_mesh(cubit, order=ORDER)
        ng_b.Save(path_b)
        mesh_b = Mesh(path_b)
        vol_b = Integrate(CF(1), mesh_b)
    except Exception as e:
        print(f"  Path B failed: {e}")
        return False

    err_a = (vol_a - cad_volume) / cad_volume * 100.0
    err_b = (vol_b - cad_volume) / cad_volume * 100.0
    diff_pct = abs(vol_a - vol_b) / cad_volume * 100.0

    print(f"  {'Path A':<10} vol={vol_a:.6e}  err={err_a:+.4e}%")
    print(f"  {'Path B':<10} vol={vol_b:.6e}  err={err_b:+.4e}%")
    print(f"  A vs B diff: {diff_pct:.4e}%", end="")

    # Volume must agree within 0.1%
    ok = diff_pct < 0.1
    print(f" [{'PASS' if ok else 'FAIL'}]")
    return ok


# ================================================================
# Main
# ================================================================
def main():
    print("=" * 70)
    print("High-Order Volume Convergence Test")
    print("  1. p-convergence (netgen .vol, order 1-5)")
    print("  2. gmsh cross-check (.msh order 1 vs 2 vs netgen)")
    print("  3. Path A (C++) vs Path B (Python) volume match")
    print("=" * 70)

    overall_pass = True

    for case_name, commands, is_flat in test_cases:
        print(f"\n{'=' * 70}")
        print(f"  {case_name}{' (flat)' if is_flat else ''}")
        print(f"{'=' * 70}")

        for cmd in commands:
            cubit.cmd(cmd)

        vol_ids = cubit.parse_cubit_list("volume", "all")
        cad_volume = sum(cubit.volume(vid).volume() for vid in vol_ids)

        n_hex = len(cubit.parse_cubit_list("hex", "all"))
        n_tet = len(cubit.parse_cubit_list("tet", "all"))
        n_pyr = len(cubit.parse_cubit_list("pyramid", "all"))
        print(f"  Elements: hex={n_hex} tet={n_tet} pyramid={n_pyr}")
        print(f"  CAD volume: {cad_volume:.6e}")

        p_pass = test_p_convergence(case_name, cad_volume, is_flat)
        g_pass = test_gmsh_vs_netgen(case_name, cad_volume, is_flat)
        ab_pass = test_path_ab(case_name, cad_volume)

        case_pass = p_pass and g_pass and ab_pass
        print(f"\n  Result: [{'PASS' if case_pass else 'FAIL'}]")
        if not case_pass:
            overall_pass = False

    print(f"\n{'=' * 70}")
    if overall_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED — see details above")
    print(f"{'=' * 70}")

    return 0 if overall_pass else 1


if __name__ == "__main__":
    sys.exit(main())
