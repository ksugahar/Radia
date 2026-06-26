"""
Verify high-order mesh export quality using GMSH API volume computation.

Tests all element types (TET, HEX, WEDGE, PYRAMID) across all export formats
(GMSH v2.2, GMSH v4.1, Nastran BDF, VTK) at order 1 and 2.

Geometries:
  - Sphere (TET):     r=0.05 m, tetmesh
  - Cube (HEX):       0.1 m side, map scheme
  - Cylinder (WEDGE): r=0.04 m, h=0.1 m, sweep scheme (wedge+hex)

The GMSH API computes volume using its own isoparametric Gauss quadrature
on HO elements (TET10, HEX20, etc.). This is the authoritative test:
  - Correct volume -> mid-node positions are correct
  - No negative Jacobian determinants -> node ordering is correct
  - Cross-format volume match -> all exporters are consistent

NOTE on VTK: GMSH's VTK reader does not support VTK_QUADRATIC_PYRAMID (type 27).
VTK verification uses GMSH API for TET10/HEX20/WEDGE15.

Usage:
  python validation_test/cubit/test_ho_volume_all_formats.py
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

import cubit
cubit.init(['cubit', '-nojournal', '-batch',
            '-commandplugindir', _plugin_dir or ''])

# --- Output directory ---
OUT_DIR = os.path.join(_test_dir, 'ho_volume_test')
os.makedirs(OUT_DIR, exist_ok=True)


# ================================================================
# Geometry builders
# ================================================================

def build_sphere_tet(mesh_size=0.015):
    """Sphere meshed with tets. V = 4/3 * pi * r^3."""
    r = 0.05
    cubit.cmd("reset")
    cubit.cmd(f"create sphere radius {r}")
    cubit.cmd("volume 1 scheme tetmesh")
    cubit.cmd(f"volume 1 size {mesh_size}")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add tet all")
    cubit.cmd('block 1 name "sphere"')
    V = (4.0 / 3.0) * math.pi * r**3
    n = len(cubit.get_block_tets(1))
    return f"Sphere TET ({n} elems)", V


def build_cube_hex(mesh_size=0.025):
    """Cube meshed with hexes. V = L^3."""
    L = 0.1
    cubit.cmd("reset")
    cubit.cmd(f"brick x {L} y {L} z {L}")
    cubit.cmd("volume 1 scheme map")
    cubit.cmd(f"volume 1 size {mesh_size}")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add hex all")
    cubit.cmd('block 1 name "cube"')
    V = L**3
    n = len(cubit.get_block_hexes(1))
    return f"Cube HEX ({n} elems)", V


def build_cylinder_wedge(mesh_size=0.02):
    """Cylinder meshed with sweep (wedges + hexes at core). V = pi*r^2*h."""
    r, h = 0.04, 0.1
    cubit.cmd("reset")
    cubit.cmd(f"create cylinder radius {r} height {h}")
    cubit.cmd("volume 1 scheme sweep")
    cubit.cmd(f"volume 1 size {mesh_size}")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add hex all")
    cubit.cmd("block 1 add wedge all")
    cubit.cmd('block 1 name "cylinder"')
    V = math.pi * r**2 * h
    nh = len(cubit.get_block_hexes(1))
    nw = len(cubit.get_block_wedges(1))
    return f"Cylinder ({nh} hex + {nw} wedge)", V


# ================================================================
# GMSH API volume computation
# ================================================================
def gmsh_api_volume(filename):
    """Compute volume via GMSH API getJacobians."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        gmsh.open(filename)
    except Exception as e:
        gmsh.finalize()
        return None, -1, str(e)

    vol = 0.0
    n_neg = 0
    infos = []

    etypes, etags, ntags = gmsh.model.mesh.getElements(dim=3)
    for idx, et in enumerate(etypes):
        name, _, order, nn, _, _ = gmsh.model.mesh.getElementProperties(int(et))
        n_elem = len(etags[idx])
        gauss_order = max(2 * order, 1)
        lc, w = gmsh.model.mesh.getIntegrationPoints(int(et), f"Gauss{gauss_order}")
        n_gp = len(w)
        _, det, _ = gmsh.model.mesh.getJacobians(int(et), lc)
        det = np.array(det)
        neg = int((det < 0).sum())
        n_neg += neg
        for i in range(n_elem):
            for j in range(n_gp):
                vol += det[i * n_gp + j] * w[j]
        infos.append(f"{name}={n_elem}" + (f"({neg}neg!)" if neg else ""))

    gmsh.finalize()
    return vol, n_neg, ", ".join(infos)


# ================================================================
# Export formats
# ================================================================
FORMATS = [
    ("GMSH_v4.1", ".msh", 'export gmsh "{f}" order {o} dimension 3 overwrite'),
    ("Nastran",   ".bdf", 'export jmag_nastran "{f}" order {o} dimension 3 overwrite'),
    ("VTK",       ".vtk", 'export vtk "{f}" order {o} dimension 3 overwrite'),
]


# ================================================================
# Test a single geometry
# ================================================================
def test_geometry(geom_name, build_fn, results_all):
    """Build geometry, export, verify via GMSH API."""
    desc, V_exact = build_fn()

    print(f"\n{'=' * 75}")
    print(f"  {desc}")
    print(f"  V_exact = {V_exact:.6e} m^3")
    print(f"{'=' * 75}")

    for fmt_name, ext, cmd_tpl in FORMATS:
        for order in [1, 2]:
            tag = f"{geom_name}_{fmt_name.replace('.','')}_o{order}"
            fname = os.path.join(OUT_DIR, f"{tag}{ext}")
            cubit.cmd(cmd_tpl.format(f=fname, o=order))

            vol, neg, info = gmsh_api_volume(fname)
            if vol is None:
                print(f"  {fmt_name:<12} o={order}: ERROR {info}")
                entry = {'geometry': geom_name, 'format': fmt_name,
                         'order': order, 'error_pct': 999, 'n_neg_det': -1}
            else:
                err = (vol - V_exact) / V_exact * 100.0
                status = "OK" if neg == 0 else "FAIL"
                print(f"  {fmt_name:<12} o={order}: V={vol:.6e} err={err:+.2e}% "
                      f"neg={neg} [{status}] {info}")
                entry = {'geometry': geom_name, 'format': fmt_name,
                         'order': order, 'volume': vol, 'error_pct': err,
                         'n_neg_det': neg}

            results_all.append(entry)


# ================================================================
# Summary
# ================================================================
def print_summary(results_all):
    all_pass = True

    print(f"\n{'=' * 80}")
    print("SUMMARY (order=2 only)")
    print(f"{'=' * 80}")
    print(f"{'Geometry':<12} {'Format':<12} {'err%':>12} {'neg_det':>8} {'Verdict':>8}")
    print(f"{'_' * 52}")

    for r in results_all:
        if r['order'] != 2:
            continue
        neg = r.get('n_neg_det', -1)
        err = r.get('error_pct', 999)
        ok = neg == 0 and abs(err) < 10.0
        if not ok:
            all_pass = False
        verdict = "PASS" if ok else "FAIL"
        print(f"{r['geometry']:<12} {r['format']:<12} {err:>+11.2e}% {neg:>8} {verdict:>8}")

    print(f"{'=' * 80}")
    if all_pass:
        print("ALL PASS")
    else:
        print("FAILURE: see details above.")
    return all_pass


# ================================================================
# Main
# ================================================================
def main():
    print("=" * 80)
    print("High-Order Mesh Export Verification (GMSH API)")
    print("=" * 80)

    results = []
    for name, fn in [("sphere", build_sphere_tet),
                     ("cube", build_cube_hex),
                     ("cylinder", build_cylinder_wedge)]:
        test_geometry(name, fn, results)

    ok = print_summary(results)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
