"""
Test: racetrack coil p-convergence.

Racetrack = straight segments + arc turns, swept along path.
This is the typical accelerator coil geometry.

Usage: python validation_test/cubit/test_racetrack_coil.py
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

_plugin_dir = os.path.join(os.path.dirname(_cubit_path), 'bin', 'plugins') \
    if _cubit_path else None
if _plugin_dir and os.path.isdir(_plugin_dir):
    os.environ['CUBIT_PLUGIN_DIR'] = _plugin_dir

import netgen.meshing
from ngsolve import Mesh, Integrate, CF
from ngsolve import TaskManager

import cubit
cubit.init(['cubit', '-nojournal', '-batch',
            '-commandplugindir', _plugin_dir or ''])

OUT_DIR = os.path.join(_test_dir, 'racetrack_test')
os.makedirs(OUT_DIR, exist_ok=True)

# Racetrack parameters
R_arc = 0.05       # arc radius (half-width)
L_straight = 0.10  # straight section length
R_wire = 0.01      # wire cross-section radius

orders = [1, 2, 3, 4, 5]

# ============================================================
# Method 1: Sweep circle along racetrack path
# ============================================================
racetrack_methods = [
    ("sweep_racetrack", [
        "reset",
        # Create racetrack by subtracting inner torus from outer torus (pill shape)
        # Outer: cylinder + two hemispheres
        f"create cylinder height {L_straight} radius {R_arc + R_wire}",
        # Inner: cylinder - wire radius
        f"create cylinder height {L_straight} radius {R_arc - R_wire}",
        # Hollow cylinder
        "subtract volume 2 from volume 1",
        # Cap with half-torus at each end
        f"create torus major radius {R_arc} minor radius {R_wire}",
        f"move volume 3 x {L_straight/2} include_merged",
        f"create torus major radius {R_arc} minor radius {R_wire}",
        f"move volume 4 x {-L_straight/2} include_merged",
        # Unite all
        "unite volume all",
        "compress",
        # Webcut for source/sink
        f"webcut volume 1 with plane yplane offset 0 noimprint nomerge",
        "delete volume 2",
        # Mesh
        "volume 1 scheme tetmesh",
        "volume 1 size auto factor 5",
        "mesh volume 1",
        "block 1 add volume 1",
        'block 1 name "coil"',
        'sideset 1 add surface 2',
        'sideset 1 name "source"',
        'sideset 2 add surface 3',
        'sideset 2 name "sink"',
    ]),
    ("native_torus_racetrack", [
        # Simpler approximation: just use native torus + webcut to make C-shape
        # This tests the webcut pattern with a gap
        "reset",
        f"create torus major radius {R_arc + L_straight/2} minor radius {R_wire}",
        # Cut 350 degrees to leave a gap
        f"webcut volume 1 with plane yplane offset 0 noimprint nomerge",
        "delete volume 2",
        # Mesh
        "volume 1 scheme tetmesh",
        "volume 1 size auto factor 5",
        "mesh volume 1",
        "block 1 add volume 1",
        'block 1 name "coil"',
        'sideset 1 add surface 2',
        'sideset 1 name "source"',
        'sideset 2 add surface 3',
        'sideset 2 name "sink"',
    ]),
]

print("Racetrack Coil p-Convergence Test")
print(f"R_arc={R_arc}, L_straight={L_straight}, R_wire={R_wire}")
print()

all_passed = True

for method_name, commands in racetrack_methods:
    print(f"=== {method_name} ===")

    for cmd in commands:
        cubit.cmd(cmd)

    # Get CAD volume for reference
    vol_ids = list(cubit.get_entities("volume"))
    cad_vol = sum(cubit.volume(vid).volume() for vid in vol_ids)

    print(f"CAD volume: {cad_vol:.6e} m^3")
    print(f"{'Order':>5} {'V_ngsolve':>14} {'Error%':>12}")
    print("-" * 35)

    prev_err = None
    for order in orders:
        vol_path = os.path.join(OUT_DIR, f"{method_name}_order{order}.vol")

        try:
            cubit.cmd(f'export netgen "{vol_path}" order {order} overwrite')
            mesh = Mesh(vol_path)
            with TaskManager():
                V_ng = Integrate(CF(1), mesh)
                err_pct = (V_ng - cad_vol) / cad_vol * 100.0
        except Exception as e:
            print(f"  {order:>3}  FAILED: {e}")
            all_passed = False
            continue

        print(f"  {order:>3} {V_ng:>14.6e} {err_pct:>+12.4e}%")

        if prev_err is not None and abs(err_pct) > abs(prev_err) and abs(err_pct) > 1e-4:
            print(f"       ** WARN: not monotonic (prev={prev_err:+.4e}%)")

        prev_err = err_pct

    print()

if all_passed:
    print("ALL PASSED")
else:
    print("SOME FAILURES")
    sys.exit(1)
