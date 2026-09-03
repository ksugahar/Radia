"""
Test Radia-NGSolve volume panel pipeline with hex sphere mesh.

Creates a hex-meshed sphere in Cubit, exports via cubit_mesh_export
at order=1,2,3, integrates volume in NGSolve, compares with analytical.

This tests the FULL panel pipeline:
  Cubit hex mesh -> export netgen order=N -> NGSolve Integrate

Usage:
  python validation_test/cubit/test_ngsolve_volume_hex_sphere.py
"""

import math
import os
import sys

from ngsolve import CF, Integrate, TaskManager

# NGSolve MUST be imported before cubit (DLL conflict avoidance)
from ngsolve import Mesh as NGMesh

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

cubit.init(['cubit', '-nojournal', '-batch', '-nographics',
            '-commandplugindir', _plugin_dir or ''])

OUT_DIR = os.path.join(_test_dir, "ngsolve_hex_sphere_test")
os.makedirs(OUT_DIR, exist_ok=True)

RADIUS = 0.05
V_EXACT = (4.0 / 3.0) * math.pi * RADIUS**3


def build_hex_sphere(mesh_size=0.015):
    """Create a hex-meshed sphere using Cubit's sphere scheme.

    Cubit's built-in sphere hex scheme creates a structured hex mesh
    from a cube core + 6 projected caps.
    """
    cubit.cmd("reset")
    cubit.cmd(f"create sphere radius {RADIUS}")
    # Use Cubit's built-in sphere hex scheme
    cubit.cmd(f"volume 1 size {mesh_size}")
    cubit.cmd("volume 1 scheme sphere")
    cubit.cmd("mesh volume 1")
    n_hex = len(cubit.parse_cubit_list("hex", "all"))
    n_tet = len(cubit.parse_cubit_list("tet", "all"))
    n_nodes = cubit.get_node_count()
    print(f"Sphere mesh: {n_hex} hexes, {n_tet} tets, {n_nodes} nodes")
    return n_hex, n_tet, n_nodes


def build_tet_sphere(mesh_size=0.015):
    """Create a tet-meshed sphere for comparison."""
    cubit.cmd("reset")
    cubit.cmd(f"create sphere radius {RADIUS}")
    cubit.cmd("volume 1 scheme tetmesh")
    cubit.cmd(f"volume 1 size {mesh_size}")
    cubit.cmd("mesh volume 1")
    n_tet = len(cubit.parse_cubit_list("tet", "all"))
    n_nodes = cubit.get_node_count()
    print(f"Sphere mesh: {n_tet} tets, {n_nodes} nodes")
    return 0, n_tet, n_nodes


def _evaluate_volume(mesh_type, build_func, mesh_size, orders):
    """Test volume accuracy for a given mesh type and orders."""
    print(f"\n{'=' * 60}")
    print(f"  {mesh_type} mesh (size={mesh_size})")
    print(f"{'=' * 60}")

    n_hex, _n_tet, _n_nodes = build_func(mesh_size)

    # CAD volume from Cubit
    v = cubit.volume(1)
    cad_vol = v.volume()
    cad_err = (cad_vol - V_EXACT) / V_EXACT * 100
    print(f"  CAD volume: {cad_vol:.6e} (err={cad_err:+.4f}%)")

    # Assign a block (elements must be in a block to be exported), then export a
    # high-order Netgen .vol per order and load it AS-IS.  A high-order .vol already
    # carries its curved mid-side nodes -- do NOT call mesh.Curve() on a loaded .vol
    # (that re-curves from absent CAD and resets every element to straight-sided).
    elem = "hex" if n_hex else "tet"
    cubit.cmd("delete block all")
    cubit.cmd(f"block 1 add {elem} all")

    results = []
    for order in orders:
        out = os.path.join(OUT_DIR, f"{mesh_type.lower()}_sphere_o{order}.vol")
        try:
            cubit.cmd(f'export netgen "{out.replace(os.sep, "/")}" order {order} overwrite')
            mesh = NGMesh(out)
            vol = Integrate(CF(1), mesh)
            err = (vol - V_EXACT) / V_EXACT * 100
            status = "OK"
        except Exception as exc:  # noqa: BLE001 - report any Cubit/export failure
            vol = 0
            err = -100
            status = f"FAIL: {exc}"

        print(f"  order={order}: V={vol:.6e}, err={err:+.4f}%  [{status}]")
        results.append({
            'order': order, 'volume': vol, 'error_pct': err, 'status': status
        })

    return results


def main():
    print("=" * 60)
    print("Radia-NGSolve Volume Panel: Hex Sphere Test")
    print(f"Sphere r={RADIUS} m, V_exact={V_EXACT:.6e} m^3")
    print("=" * 60)

    all_results = {}

    # Tet mesh (baseline)
    with TaskManager():
        all_results['tet'] = _evaluate_volume(
            "Tet", build_tet_sphere, 0.015, [1, 2, 3])

        # Hex mesh (finer to get enough elements for convergence)
        all_results['hex'] = _evaluate_volume(
            "Hex", build_hex_sphere, 0.008, [1, 2, 3])

    # Summary
    print(f"\n{'=' * 60}")
    print(f"{'Mesh':<6} {'Order':>5} {'Volume':>14} {'Error':>10} {'Status':>8}")
    print(f"{'-' * 60}")
    for mtype, results in all_results.items():
        for r in results:
            s = "OK" if "OK" in r['status'] else "FAIL"
            print(f"{mtype:<6} {r['order']:>5} {r['volume']:>14.6e} "
                  f"{r['error_pct']:>+9.4f}% {s:>8}")
    print(f"{'=' * 60}")

    # Pass/fail
    # Pass criteria: order>=2 must have < 1% error (proves curving works)
    # order=1 is allowed to have larger error (linear approximation)
    all_ok = True
    for results in all_results.values():
        for r in results:
            if "OK" not in r['status']:
                all_ok = False
            if r['order'] >= 2 and abs(r['error_pct']) > 1.0:
                all_ok = False
    if all_ok:
        print("\nALL PASS")
    else:
        print("\nFAILURES detected")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
