"""
High-order volume convergence test for mixed element types.

Tests:
  1. p-convergence: export netgen .vol at order 1-5, volume must converge to CAD

Cases:
  1. Flat brick (hex + tet + pyramid) -- exact volume expected
  2. Cylinder hex (scheme map) -- curved surface
  3. Cylinder webcut hex (multi-volume) -- curved, 4 volumes
  4. Cylinder boundary-layer hex -- curved, BL mesh
  5. Loft circle-to-rectangle hex -- complex shape

Usage:
  python validation_test/cubit/test_ho_volume_convergence.py
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
from ngsolve import Mesh, Integrate, CF, BND, BBND
from ngsolve import TaskManager

import cubit
cubit.init(['cubit', '-nojournal', '-batch',
            '-commandplugindir', _plugin_dir or ''])
OUT_DIR = os.path.join(_test_dir, 'ho_convergence_test')
os.makedirs(OUT_DIR, exist_ok=True)


# ================================================================
# Test case definitions
# ================================================================
test_cases = [
# (name, commands, is_flat, area_warn)
# area_warn=True: area p-convergence is WARN-only (ACIS loft trimming limitation)
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
    ], True, False),
    ("case02_cylinder_map", [
        "reset",
        "cylinder height 1 radius 1",
        "volume 1 size 0.5", "mesh volume 1",
        "block 1 add hex all", 'block 1 name "cylinder_hex"',
        "Volume all scale 0.001",
    ], False, False),
    ("case03_cylinder_webcut", [
        "reset",
        "cylinder height 1 radius 1",
        "webcut volume all with plane xplane offset 0 imprint merge",
        "webcut volume all with plane yplane offset 0 imprint merge",
        "imprint all", "merge all",
        "volume all size 1", "mesh volume all",
        "block 1 add hex all", 'block 1 name "cylinder_hex"',
        "Volume all scale 0.001",
    ], False, False),
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
    ], False, False),
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
    ], False, True),   # area_warn: ACIS loft surface 5-curve trimming
    ("case06_cyl_sphere_BL", [
        "reset",
        "create Cylinder height 1 radius 0.5",
        "create sphere radius 0.5",
        "move Volume 2 z 0.5 include_merged",
        "unite volume all", "compress",
        "create boundary_layer 1",
        "modify boundary_layer 1 uniform height 0.02 growth 1.2 layers 3",
        "modify boundary_layer 1 add surface 1 volume 1 surface 3 volume 1",
        "modify boundary_layer 1 continuity on",
        "volume 1 scheme tetmesh", "volume 1 size 0.2",
        "mesh volume 1",
        "block 1 add tet all", 'block 1 name "tet"',
        "block 2 add wedge all", 'block 2 name "wedge"',
        "Volume all scale 0.001",
    ], False, False),
    ("case07_sweep_torus", [
        "reset",
        "create curve arc radius 0.5 center location 2 0 0 "
        "normal 0 1 0 start angle 0 stop angle 360",
        "create surface curve 1",
        "sweep surface 1 zaxis angle 360",
        "volume 1 scheme tetmesh", "volume 1 size auto factor 5",
        "mesh volume 1",
        "block 1 add volume 1", 'block 1 name "torus"',
        "Volume all scale 0.001",
    ], False, False),
    ("case08_half_torus_webcut", [
        "reset",
        "create torus major radius 2 minor radius 0.5",
        "webcut volume 1 with plane xplane noimprint nomerge",
        "delete volume 2",
        "volume 1 scheme tetmesh", "volume 1 size auto factor 5",
        "mesh volume 1",
        "block 1 add volume 1", 'block 1 name "coil"',
        'sideset 1 add surface 2', 'sideset 1 name "source"',
        'sideset 2 add surface 3', 'sideset 2 name "sink"',
        "Volume all scale 0.001",
    ], False, False),
]
# NOTE: Racetrack (boolean) coil tested separately in test_racetrack_coil.py
# Boolean trimming surfaces saturate at ~0.025% at high order (ACIS limitation).


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


def wedge6_volume(pts):
    """Volume of linear wedge using 2-point Gauss in zeta, 1-point in triangle."""
    gp = 1.0 / math.sqrt(3.0)
    gauss = [(1/3.0, 1/3.0, -gp), (1/3.0, 1/3.0, gp)]
    vol = 0.0
    for L1, L2, zeta in gauss:
        L3 = 1 - L1 - L2
        dNdL1 = np.array([(1-zeta)/2, 0, -(1-zeta)/2,
                           (1+zeta)/2, 0, -(1+zeta)/2])
        dNdL2 = np.array([0, (1-zeta)/2, -(1-zeta)/2,
                           0, (1+zeta)/2, -(1+zeta)/2])
        dNdzeta = np.array([-L1/2, -L2/2, -L3/2, L1/2, L2/2, L3/2])
        J = np.zeros((3, 3))
        for i in range(6):
            J[0] += dNdL1[i] * pts[i]
            J[1] += dNdL2[i] * pts[i]
            J[2] += dNdzeta[i] * pts[i]
        vol += 0.5 * abs(np.linalg.det(J))  # 0.5 = ref triangle area
    return vol


def wedge15_volume(coords):
    """Volume of WEDGE15 using 3-point triangle x 3-point line Gauss quadrature.

    GMSH WEDGE15 serendipity shape functions in (u,v,w) reference space:
      u,v: triangle (u+v<=1), w: [-1,1]
      L0=1-u-v, L1=u, L2=v

    Corner (serendipity, NOT simple product):
      Bottom: N_i = L_i*(1-w)*(2*L_i - 2 - w)/2
      Top:    N_i = L_i*(1+w)*(2*L_i - 2 + w)/2
    Mid-edge triangle: N = 4*La*Lb*(1 +/- w)/2
    Mid-edge vertical: N = L_i*(1-w^2)
    """
    tri_pts = [(1 / 6.0, 1 / 6.0), (2 / 3.0, 1 / 6.0), (1 / 6.0, 2 / 3.0)]
    tri_w = 1 / 6.0

    gp = math.sqrt(3.0 / 5.0)
    line_pts = [-gp, 0, gp]
    line_w = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]

    # dL/du, dL/dv for L0=1-u-v, L1=u, L2=v
    dLdu = [-1, 1, 0]
    dLdv = [-1, 0, 1]

    vol = 0.0
    for u, v in tri_pts:
        Lv = [1 - u - v, u, v]
        for li in range(3):
            w = line_pts[li]
            wt = tri_w * line_w[li]
            h0 = 1 - w * w

            dNdu = np.zeros(15)
            dNdv = np.zeros(15)
            dNdw = np.zeros(15)

            # Bottom corners: N = Li*(1-w)*(2*Li-2-w)/2
            for i in range(3):
                Li = Lv[i]
                a = 1 - w
                dNdu[i] = dLdu[i] * a * (4 * Li - 2 - w) / 2
                dNdv[i] = dLdv[i] * a * (4 * Li - 2 - w) / 2
                dNdw[i] = Li * (-2 * Li + 1 + 2 * w) / 2

            # Top corners: N = Li*(1+w)*(2*Li-2+w)/2
            for i in range(3):
                Li = Lv[i]
                a = 1 + w
                dNdu[3 + i] = dLdu[i] * a * (4 * Li - 2 + w) / 2
                dNdv[3 + i] = dLdv[i] * a * (4 * Li - 2 + w) / 2
                dNdw[3 + i] = Li * (2 * Li - 1 + 2 * w) / 2

            # Mid-edge bottom (6,7,8): N = 4*La*Lb*(1-w)/2
            hm = (1 - w) / 2
            for idx, ia, ib in [(6, 0, 1), (7, 1, 2), (8, 2, 0)]:
                La, Lb = Lv[ia], Lv[ib]
                dNdu[idx] = 4 * (dLdu[ia] * Lb + La * dLdu[ib]) * hm
                dNdv[idx] = 4 * (dLdv[ia] * Lb + La * dLdv[ib]) * hm
                dNdw[idx] = 4 * La * Lb * (-0.5)

            # Mid-edge top (9,10,11): N = 4*La*Lb*(1+w)/2
            hp = (1 + w) / 2
            for idx, ia, ib in [(9, 0, 1), (10, 1, 2), (11, 2, 0)]:
                La, Lb = Lv[ia], Lv[ib]
                dNdu[idx] = 4 * (dLdu[ia] * Lb + La * dLdu[ib]) * hp
                dNdv[idx] = 4 * (dLdv[ia] * Lb + La * dLdv[ib]) * hp
                dNdw[idx] = 4 * La * Lb * (0.5)

            # Mid-edge vertical (12,13,14): N = Li*(1-w^2)
            for idx, i in [(12, 0), (13, 1), (14, 2)]:
                dNdu[idx] = dLdu[i] * h0
                dNdv[idx] = dLdv[i] * h0
                dNdw[idx] = Lv[i] * (-2 * w)

            J = np.zeros((3, 3))
            for n in range(15):
                J[0] += dNdu[n] * coords[n]
                J[1] += dNdv[n] * coords[n]
                J[2] += dNdw[n] * coords[n]
            vol += wt * abs(np.linalg.det(J))
    return vol


def pyramid5_volume(pts):
    v1 = tet4_volume(np.array([pts[0], pts[1], pts[2], pts[4]]))
    v2 = tet4_volume(np.array([pts[0], pts[2], pts[3], pts[4]]))
    return v1 + v2


# ================================================================
# Test 1: p-convergence (netgen .vol, order 1-5, volume + area)
# ================================================================
def test_p_convergence(case_name, cad_volume, cad_area, is_flat, area_warn=False):
    """Export .vol at orders 1-5, check volume and area convergence."""
    print(f"\n  p-convergence (CAD vol={cad_volume:.6e}, area={cad_area:.6e}):")
    print(f"  {'Order':>5} {'Volume':>14} {'V_err':>14} {'Area':>14} {'A_err':>14} {'Verdict':>8}")
    print(f"  {'-----':>5} {'-'*14:>14} {'-'*14:>14} {'-'*14:>14} {'-'*14:>14} {'--------':>8}")

    v_errors = []
    a_errors = []
    all_pass = True

    for order in range(1, 6):
        vol_path = os.path.join(OUT_DIR, f"{case_name}_o{order}.vol")
        cmd = f'export netgen "{vol_path}" order {order} overwrite'
        try:
            cubit.cmd(cmd)
            mesh = Mesh(vol_path)
            with TaskManager():
                vol = Integrate(CF(1), mesh)
                area = Integrate(CF(1), mesh, BND)
        except Exception as e:
            print(f"  {order:>5} FAILED: {e}")
            v_errors.append(None)
            a_errors.append(None)
            all_pass = False
            continue

        v_err = (vol - cad_volume) / cad_volume * 100.0
        a_err = (area - cad_area) / cad_area * 100.0
        verdict = "OK"

        if is_flat:
            if abs(v_err) > 0.01 or abs(a_err) > 0.01:
                verdict = "FAIL"
                all_pass = False
        else:
            # Check volume regression
            if v_errors and v_errors[-1] is not None:
                if abs(v_err) > abs(v_errors[-1]) * 100 and abs(v_err) > 1e-2:
                    verdict = "REGRESS"
                    all_pass = False
            # Check area regression
            if a_errors and a_errors[-1] is not None:
                if abs(a_err) > abs(a_errors[-1]) * 100 and abs(a_err) > 1e-2:
                    if area_warn:
                        verdict = "WARN(area)"
                    else:
                        verdict = "REGRESS"
                        all_pass = False

        print(f"  {order:>5} {vol:>14.6e} {v_err:>+13.4e}% {area:>14.6e} {a_err:>+13.4e}% {verdict:>8}")
        v_errors.append(v_err)
        a_errors.append(a_err)

    # Convergence check: order 2 must improve over order 1
    if not is_flat:
        if len(v_errors) >= 2 and v_errors[0] is not None and v_errors[1] is not None:
            if abs(v_errors[1]) >= abs(v_errors[0]):
                print(f"  volume: order 2 not better than order 1: FAIL")
                all_pass = False
        if len(a_errors) >= 2 and a_errors[0] is not None and a_errors[1] is not None:
            if abs(a_errors[1]) >= abs(a_errors[0]):
                if area_warn:
                    print(f"  area: order 2 not better than order 1: WARN (ACIS loft)")
                else:
                    print(f"  area: order 2 not better than order 1: FAIL")
                    all_pass = False

    return all_pass


# ================================================================
# Main
# ================================================================
def main():
    print("=" * 70)
    print("High-Order Volume Convergence Test")
    print("  p-convergence (netgen .vol, order 1-5)")
    print("=" * 70)

    overall_pass = True

    for case_name, commands, is_flat, area_warn in test_cases:
        print(f"\n{'=' * 70}")
        print(f"  {case_name}{' (flat)' if is_flat else ''}")
        print(f"{'=' * 70}")

        for cmd in commands:
            cubit.cmd(cmd)

        vol_ids = cubit.parse_cubit_list("volume", "all")
        cad_volume = sum(cubit.volume(vid).volume() for vid in vol_ids)
        surf_ids = cubit.parse_cubit_list("surface", "all")
        cad_area = sum(cubit.surface(sid).area() for sid in surf_ids)

        n_hex = len(cubit.parse_cubit_list("hex", "all"))
        n_tet = len(cubit.parse_cubit_list("tet", "all"))
        n_pyr = len(cubit.parse_cubit_list("pyramid", "all"))
        print(f"  Elements: hex={n_hex} tet={n_tet} pyramid={n_pyr}")
        print(f"  CAD volume: {cad_volume:.6e}, area: {cad_area:.6e}")

        case_pass = test_p_convergence(case_name, cad_volume, cad_area, is_flat, area_warn)
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
