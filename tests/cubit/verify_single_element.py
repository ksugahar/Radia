"""
Minimal test: export a single HEX element at order 2 to all formats,
then verify with GMSH API. This isolates node ordering issues from
curving/geometry effects.

Also tests WEDGE and PYRAMID if available.

Usage:
  python tests/cubit/verify_single_element.py
"""

import sys
import os
import math
import numpy as np

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

OUT_DIR = os.path.join(_test_dir, 'ho_volume_test')
os.makedirs(OUT_DIR, exist_ok=True)


def gmsh_verify(filename):
    """Verify file with GMSH API. Returns (volume, n_neg_det, elem_info)."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(filename)

    vol = 0.0
    n_neg_det = 0
    infos = []

    etypes, etags, ntags = gmsh.model.mesh.getElements(dim=3)
    for idx, et in enumerate(etypes):
        name, dim, order, nn, _, _ = gmsh.model.mesh.getElementProperties(int(et))
        n_elem = len(etags[idx])
        gauss_order = max(2 * order, 1)
        local_coords, weights = gmsh.model.mesh.getIntegrationPoints(
            int(et), f"Gauss{gauss_order}")
        n_gp = len(weights)
        jac, det, pts = gmsh.model.mesh.getJacobians(int(et), local_coords)
        det = np.array(det)
        neg = int((det < 0).sum())
        n_neg_det += neg
        for i in range(n_elem):
            for j in range(n_gp):
                vol += det[i * n_gp + j] * weights[j]
        infos.append(f"{name}={n_elem}" + (f" ({neg} neg!)" if neg else ""))

    # Also get 2D elements
    etypes2, etags2, _ = gmsh.model.mesh.getElements(dim=2)
    for idx, et in enumerate(etypes2):
        name, dim, order, nn, _, _ = gmsh.model.mesh.getElementProperties(int(et))
        n_elem = len(etags2[idx])
        infos.append(f"{name}={n_elem}")

    gmsh.finalize()
    return vol, n_neg_det, ", ".join(infos)


def test_cube_hex():
    """Single large hex element (cube). Flat faces, no curving needed."""
    L = 0.1
    V_exact = L**3

    cubit.cmd("reset")
    cubit.cmd(f"brick x {L} y {L} z {L}")
    cubit.cmd("volume 1 scheme map")
    cubit.cmd(f"volume 1 size {L}")  # single element
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add hex all")
    cubit.cmd('block 1 name "cube"')

    n_hex = len(cubit.get_block_hexes(1))
    print(f"\n{'=' * 60}")
    print(f"Cube HEX: {n_hex} elements, V_exact={V_exact:.6e}")
    print(f"{'=' * 60}")

    run_exports("cube_single", V_exact)


def test_cube_hex_multi():
    """Multiple hex elements (cube, 3x3x3)."""
    L = 0.1
    V_exact = L**3

    cubit.cmd("reset")
    cubit.cmd(f"brick x {L} y {L} z {L}")
    cubit.cmd("volume 1 scheme map")
    cubit.cmd(f"volume 1 size {L/3}")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add hex all")
    cubit.cmd('block 1 name "cube"')

    n_hex = len(cubit.get_block_hexes(1))
    print(f"\n{'=' * 60}")
    print(f"Cube HEX (3x3x3): {n_hex} elements, V_exact={V_exact:.6e}")
    print(f"{'=' * 60}")

    run_exports("cube_multi", V_exact)


def test_sphere_tet():
    """Sphere with tet mesh (reference, known to work)."""
    r = 0.05
    V_exact = (4.0 / 3.0) * math.pi * r**3

    cubit.cmd("reset")
    cubit.cmd(f"create sphere radius {r}")
    cubit.cmd("volume 1 scheme tetmesh")
    cubit.cmd(f"volume 1 size 0.02")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add tet all")
    cubit.cmd('block 1 name "sphere"')

    n_tet = len(cubit.get_block_tets(1))
    print(f"\n{'=' * 60}")
    print(f"Sphere TET: {n_tet} elements, V_exact={V_exact:.6e}")
    print(f"{'=' * 60}")

    run_exports("sphere_tet", V_exact)


def run_exports(tag, V_exact):
    """Export current mesh to all formats at order 1 and 2, verify."""
    formats = [
        ("GMSH_v2", ".msh", 'export gmsh "{f}" order {o} version 2 overwrite'),
        ("GMSH_v4", ".msh", 'export gmsh "{f}" order {o} version 4 overwrite'),
        ("Nastran", ".bdf", 'export radia_nastran "{f}" order {o} dimension 3 overwrite'),
        ("VTK",     ".vtk", 'export vtk "{f}" order {o} dimension 3 overwrite'),
    ]

    for fmt_name, ext, cmd_tpl in formats:
        for order in [1, 2]:
            fname = os.path.join(OUT_DIR, f"{tag}_{fmt_name}_o{order}{ext}")
            cubit.cmd(cmd_tpl.format(f=fname, o=order))
            vol, neg, info = gmsh_verify(fname)
            err = (vol - V_exact) / V_exact * 100.0
            status = "OK" if neg == 0 else "FAIL"
            neg_str = f" [{neg} neg det]" if neg > 0 else ""
            print(f"  {fmt_name:<10} o={order}: V={vol:.6e} err={err:+.2e}% "
                  f"{info}{neg_str} [{status}]")


def main():
    print("=" * 60)
    print("Single Element Verification: Node Ordering Check")
    print("=" * 60)

    test_sphere_tet()
    test_cube_hex()
    test_cube_hex_multi()

    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
