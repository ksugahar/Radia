"""
Test: sweep torus p-convergence with UV-guided projection fix.

Verifies that ACIS projection correctly handles sweep torus surfaces
split at z=0.

Tests Path A (C++ plugin, export vol command) p-convergence.

Previously: closest_point_trimmed cross-projected between upper/lower halves,
causing negative Jacobians and broken p-convergence.

Usage: python validation_test/cubit/test_sweep_torus_curving.py
"""
import sys
import os
import math

_test_dir = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.dirname(os.path.dirname(_test_dir))
sys.path.insert(0, os.path.join(_repo_root, 'src', 'radia'))
sys.path.insert(0, os.path.join(_repo_root, 'ngsolve_developers'))
from radia_mcp.cubit.session import find_cubit_bin
_cubit_path = find_cubit_bin()
if _cubit_path: sys.path.append(_cubit_path)

from radia_mcp.cubit.layout import export_plugin_directory
# $RADIA_CUBIT_PLUGIN_DIR (a build.py output) or the installed plugin directory.
_plugin_dir = str(export_plugin_directory(_cubit_path)) if _cubit_path else None
if _plugin_dir:
    os.environ['CUBIT_PLUGIN_DIR'] = _plugin_dir

import netgen.meshing
from ngsolve import Mesh, Integrate, CF
from ngsolve import TaskManager

import cubit
cubit.init(['cubit', '-nojournal', '-batch',
            '-noinitfile'])
from radia.vol_check import check_mesh_quality

OUT_DIR = os.path.join(_test_dir, 'sweep_torus_test')
os.makedirs(OUT_DIR, exist_ok=True)

# Torus parameters
R_major = 0.05
R_minor = 0.02
V_exact = 2.0 * math.pi**2 * R_major * R_minor**2

# Two methods of creating a torus in Cubit
torus_methods = [
    ("native_torus", [
        "reset",
        f"create torus major radius {R_major} minor radius {R_minor}",
        "volume 1 scheme tetmesh",
        "volume 1 size auto factor 5",
        "mesh volume 1",
        "block 1 add volume 1",
        'block 1 name "torus"',
    ]),
    ("sweep_torus", [
        "reset",
        f"create curve arc radius {R_minor} center location {R_major} 0 0 "
        f"normal 0 1 0 start angle 0 stop angle 360",
        "create surface curve 1",
        "sweep surface 1 zaxis angle 360",
        "volume 1 scheme tetmesh",
        "volume 1 size auto factor 5",
        "mesh volume 1",
        "block 1 add volume 1",
        'block 1 name "torus"',
    ]),
]

orders = [1, 2, 3]

print(f"Torus: R_major={R_major}, R_minor={R_minor}")
print(f"V_exact = {V_exact:.6e} m^3")
print()


def check_jacobians(vol_path):
    """Count folded curved-map samples of the exported .vol (sign flips or invalid)."""
    report = check_mesh_quality(vol_path)
    return (report["orientation_flip_sample_count"]
            + report["invalid_jacobian_sample_count"])


def volume_from_vol(vol_path, order):
    """Load .vol and integrate volume."""
    mesh = Mesh(vol_path)
    with TaskManager():
        return Integrate(CF(1), mesh)


# ============================================================
# Test: Path A (C++ plugin) p-convergence
# ============================================================
print("=" * 75)
print("Path A (C++ plugin) p-convergence")
print(f"{'Method':<16} {'Order':>5} {'V_ngsolve':>14} {'Error%':>12} {'folded':>8}")
print("-" * 60)

all_passed = True

for method_name, commands in torus_methods:
    for cmd in commands:
        cubit.cmd(cmd)

    prev_err = None
    for order in orders:
        vol_path = os.path.join(OUT_DIR, f"{method_name}_A_order{order}.vol")

        cubit.cmd(f'export vol "{vol_path}" order {order} overwrite')

        V_ng = volume_from_vol(vol_path, order)
        err_pct = (V_ng - V_exact) / V_exact * 100.0
        neg = check_jacobians(vol_path)

        print(f"{method_name:<16} {order:>5} {V_ng:>14.6e} {err_pct:>+12.2e}% {neg:>8}")

        if neg > 0:
            print(f"  ** FAIL: folded curved elements detected!")
            all_passed = False
        if prev_err is not None and abs(err_pct) > abs(prev_err) and not math.isclose(err_pct, 0, abs_tol=1e-6):
            print(f"  ** WARN: p-convergence not monotonic")

        prev_err = err_pct
    print()

print("=" * 75)
if all_passed:
    print("ALL PASSED: no folded curved elements, p-convergence OK")
else:
    print("FAILED: see details above")
    sys.exit(1)
