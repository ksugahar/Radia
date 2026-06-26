"""
Test: BEM inductance p-convergence for half-torus coil.

Exports half-torus at mesh order 1-4, computes BEM inductance via
saddle-point EFIE, checks convergence.

Analytical reference: L = mu_0 * R * (ln(8R/a) - 2) for full torus,
half-torus = L_full / 2 (approximately).

Usage: python validation_test/cubit/test_inductance_p_convergence.py
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

# NGSolve FIRST, then scipy BEFORE cubit.init() (Cubit pollutes sys.path)
import netgen.meshing
from ngsolve import Mesh, Integrate, CF, BND
from ngsolve import TaskManager
import scipy.sparse  # noqa: F401 — must import before cubit.init()

# Cubit
import cubit
cubit.init(['cubit', '-nojournal', '-batch',
            '-commandplugindir', _plugin_dir or ''])

# Inductance extractor
sys.path.insert(0, os.path.join(_repo_root, 'src', 'radia', 'panels'))
from calc_inductance import extract_inductance_vol

OUT_DIR = os.path.join(_test_dir, 'inductance_p_conv')
os.makedirs(OUT_DIR, exist_ok=True)

# Parameters
R_major = 0.11    # major radius [m]
R_minor = 0.01    # minor radius (wire) [m]
MU_0 = 4e-7 * math.pi

L_full = MU_0 * R_major * (math.log(8 * R_major / R_minor) - 2.0)
L_half_approx = L_full / 2.0

# Two coil creation methods
coil_methods = [
    ("native_torus_webcut", [
        "reset",
        f"create torus major radius {R_major} minor radius {R_minor}",
        "webcut volume 1 with plane xplane noimprint nomerge",
        "delete volume 2",
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
    ("boolean_half_torus", [
        # Same half-torus shape, but built via boolean (sphere subtract)
        # to test ACIS boolean trimming surfaces
        "reset",
        # Full torus via boolean: cylinder rotate - inner cylinder rotate
        f"create cylinder height {2*R_minor} radius {R_major + R_minor}",
        f"create cylinder height {2*R_minor} radius {R_major - R_minor}",
        "subtract volume 2 from volume 1",
        # Revolve doesn't work easily, so use sphere intersection instead:
        # Create torus by sweep, then webcut
        f"create torus major radius {R_major} minor radius {R_minor}",
        # Boolean: subtract a box to keep only x>0 half (same as webcut)
        f"create brick x {4*R_major} y {4*R_major} z {4*R_major}",
        f"move volume 4 x {-2*R_major} include_merged",
        "subtract volume 4 from volume 3",
        # Clean up: delete the hollow cylinder, keep only boolean torus half
        "delete volume 1",
        "compress",
        # Mesh
        "volume 1 scheme tetmesh",
        "volume 1 size auto factor 5",
        "mesh volume 1",
        "block 1 add volume 1",
        'block 1 name "coil"',
    ]),
]

orders = [1, 2, 3, 4]

print(f"R_major={R_major}, R_minor={R_minor}")
print(f"L_full (analytical) = {L_full*1e9:.2f} nH")
print(f"L_half (approx)     = {L_half_approx*1e9:.2f} nH")
print()

for method_name, cmds in coil_methods:
    print(f"=== {method_name} ===")
    for cmd in cmds:
        cubit.cmd(cmd)

    # Find source/sink sidesets — may not exist for boolean method
    has_sidesets = len(list(cubit.get_sideset_id_list())) >= 2

    # If no sidesets, find webcut faces and add them
    if not has_sidesets:
        # After webcut, terminal faces are the planar surfaces
        all_surfs = list(cubit.parse_cubit_list("surface", "all"))
        # Find surfaces with smallest area (terminal faces)
        areas = [(sid, cubit.surface(sid).area()) for sid in all_surfs]
        areas.sort(key=lambda x: x[1])
        if len(areas) >= 2:
            cubit.cmd(f'sideset 1 add surface {areas[0][0]}')
            cubit.cmd('sideset 1 name "source"')
            cubit.cmd(f'sideset 2 add surface {areas[1][0]}')
            cubit.cmd('sideset 2 name "sink"')

    cad_vol = sum(cubit.volume(vid).volume()
                  for vid in cubit.get_entities("volume"))

    print(f"CAD volume: {cad_vol:.6e}")
    print(f"{'Order':>5} {'V_err%':>10} {'L [nH]':>10} {'n_DOF':>8}")
    print("-" * 40)

    for order in orders:
        vol_path = os.path.join(OUT_DIR, f"{method_name}_o{order}.vol")
        cubit.cmd(f'export netgen "{vol_path}" order {order} overwrite')

        mesh = Mesh(vol_path)
        with TaskManager():
            vol = Integrate(CF(1), mesh)
            v_err = (vol - cad_vol) / cad_vol * 100

            try:
                result = extract_inductance_vol(
                    vol_path, source_label="source", sink_label="sink", fes_order=0)
                if 'error' in result:
                    L_nH = float('nan')
                    ndof = 0
                    sys.stderr.write(f"  {method_name} order {order}: {result['error']}\n")
                else:
                    L_nH = result['inductance_H'] * 1e9
                    ndof = result['n_dofs']
            except Exception as e:
                L_nH = float('nan')
                ndof = 0
                sys.stderr.write(f"  {method_name} order {order}: {e}\n")

            print(f"  {order:>3} {v_err:>+10.4e}% {L_nH:>10.2f} {ndof:>8}")

    print()
