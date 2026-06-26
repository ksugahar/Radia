"""
Multi-geometry volume accuracy test (export netgen C++ command).

Shapes:
  1. Sphere (R=0.05) — single curved surface
  2. Cylinder (R=0.03, H=0.1) — curved + flat surfaces
  3. Torus (R1=0.05, R2=0.02) — doubly curved
  4. Brick (0.1 x 0.05 x 0.02) — all flat (curving should be trivial)
  5. Two-volume: sphere + brick (multi-material)

Usage: python validation_test/cubit/test_vol_multi_geometry.py
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
os.environ['CUBIT_PLUGIN_DIR'] = os.path.join(_cubit_path, 'plugins') if _cubit_path else ''

import netgen.meshing
from ngsolve import Mesh, Integrate, CF, BND, TaskManager
import cubit
cubit.init(['cubit', '-nojournal', '-batch',
            '-commandplugindir', os.environ['CUBIT_PLUGIN_DIR']])
OUT_DIR = os.path.join(_test_dir, 'export_netgen_test')
os.makedirs(OUT_DIR, exist_ok=True)

ORDER = 2

# ================================================================
# Test cases: (name, setup_commands, V_exact, A_exact)
# ================================================================
R_sph = 0.05
R_cyl = 0.03
H_cyl = 0.1
R1_tor = 0.05
R2_tor = 0.02
Lx, Ly, Lz = 0.1, 0.05, 0.02

test_cases = [
    ("sphere",
     [f"create sphere radius {R_sph}",
      "volume 1 scheme tetmesh",
      "volume 1 size auto factor 5",
      "mesh volume 1",
      "block 1 add volume 1",
      'block 1 name "sphere"'],
     (4.0/3.0)*math.pi*R_sph**3,
     4.0*math.pi*R_sph**2),

    ("cylinder",
     [f"create cylinder height {H_cyl} radius {R_cyl}",
      "volume 1 scheme tetmesh",
      "volume 1 size auto factor 5",
      "mesh volume 1",
      "block 1 add volume 1",
      'block 1 name "cylinder"'],
     math.pi*R_cyl**2*H_cyl,
     2*math.pi*R_cyl*H_cyl + 2*math.pi*R_cyl**2),

    ("torus",
     [f"create torus major radius {R1_tor} minor radius {R2_tor}",
      "volume 1 scheme tetmesh",
      "volume 1 size auto factor 5",
      "mesh volume 1",
      "block 1 add volume 1",
      'block 1 name "torus"'],
     2*math.pi**2*R1_tor*R2_tor**2,
     4*math.pi**2*R1_tor*R2_tor),

    ("brick",
     [f"create brick x {Lx} y {Ly} z {Lz}",
      "volume 1 scheme tetmesh",
      "volume 1 size auto factor 5",
      "mesh volume 1",
      "block 1 add volume 1",
      'block 1 name "brick"'],
     Lx*Ly*Lz,
     2*(Lx*Ly + Ly*Lz + Lz*Lx)),

    ("two_volume",
     [f"create sphere radius {R_sph}",
      f"create brick x {3*R_sph} y {3*R_sph} z {3*R_sph}",
      f"volume 2 move {4*R_sph} 0 0",
      "volume all scheme tetmesh",
      "volume all size auto factor 5",
      "mesh volume all",
      "block 1 add volume 1",
      'block 1 name "iron"',
      "block 2 add volume 2",
      'block 2 name "air"'],
     (4.0/3.0)*math.pi*R_sph**3 + (3*R_sph)**3,
     None),  # complex boundary, skip area check

    ("loft_rect",
     [# Bottom: rectangle in XY plane
      "create vertex 0 0 0",
      "create vertex 0.1 0 0",
      "create vertex 0.1 0.04 0",
      "create vertex 0 0.04 0",
      "create surface vertex 1 2 3 4",
      # Top: smaller rotated rectangle at z=0.08
      "create vertex 0.02 0.01 0.08",
      "create vertex 0.08 0.005 0.08",
      "create vertex 0.08 0.035 0.08",
      "create vertex 0.02 0.03 0.08",
      "create surface vertex 5 6 7 8",
      "create volume loft surface 1 2",
      "volume all scheme tetmesh",
      "volume all size auto factor 5",
      "mesh volume all",
      "block 1 add volume all",
      'block 1 name "loft_rect"'],
     None,  # no analytical volume
     None),

    ("loft_circle",
     [# Bottom: circle R=0.04 at z=0
      "create surface circle radius 0.04 zplane",
      # Top: circle R=0.02 at z=0.1, offset in x
      "create surface circle radius 0.02 zplane",
      "surface 2 move 0.01 0 0.1",
      "create volume loft surface 1 2",
      "volume all scheme tetmesh",
      "volume all size auto factor 5",
      "mesh volume all",
      "block 1 add volume all",
      'block 1 name "loft_circle"'],
     None,
     None),

    ("cone",
     ["create frustum height 0.08 radius 0.04 top 0.01",
      "volume 1 scheme tetmesh",
      "volume 1 size auto factor 5",
      "mesh volume 1",
      "block 1 add volume 1",
      'block 1 name "cone"'],
     (1.0/3.0)*math.pi*0.08*(0.04**2 + 0.04*0.01 + 0.01**2),
     None),  # complex area

    # --- Hex meshes with curved surfaces ---
    ("hex_cylinder",
     [f"create cylinder height {H_cyl} radius {R_cyl}",
      "volume 1 size auto factor 5",
      "mesh volume 1",
      "block 1 add volume 1",
      'block 1 name "hex_cyl"'],
     math.pi*R_cyl**2*H_cyl,
     2*math.pi*R_cyl*H_cyl + 2*math.pi*R_cyl**2),

    ("hex_sphere",
     [f"create sphere radius {R_sph}",
      "volume 1 size auto factor 5",
      "mesh volume 1",
      "block 1 add volume 1",
      'block 1 name "hex_sph"'],
     (4.0/3.0)*math.pi*R_sph**3,
     4.0*math.pi*R_sph**2),

]

print("=" * 90)
print(f"Multi-geometry volume accuracy test (order {ORDER})")
print("=" * 90)

header = f"{'Shape':<14} {'V_err%':>12} {'A_err%':>12} {'ne':>7} {'nse':>6} {'nedges':>7} {'n_ho_e':>7} {'Result':>8}"
print(header)
print("-" * len(header))

all_pass = True

for name, cmds, v_exact, a_exact in test_cases:
    cubit.cmd("reset")
    for cmd in cmds:
        cubit.cmd(cmd)

    # --- Export netgen ---
    vol_path = os.path.join(OUT_DIR, f"multi_{name}.vol")
    cubit.cmd(f'export netgen "{vol_path}" order {ORDER} overwrite')
    with TaskManager():
        mesh = Mesh(vol_path)
        vol = Integrate(CF(1), mesh)
        area = Integrate(CF(1), mesh, BND)
    verr = (vol - v_exact) / v_exact * 100 if v_exact else float('nan')
    aerr = (area - a_exact) / a_exact * 100 if a_exact else float('nan')

    # --- Structural info ---
    def vol_info(path):
        info = {}
        with open(path) as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            s = line.strip()
            if s == 'surfaceelements':
                info['nse'] = int(lines[i+1].strip())
            elif s == 'curvedelements':
                ned = int(lines[i+1].strip())
                info['nedges'] = ned
                info['n_ho_e'] = sum(1 for j in range(i+2, i+2+ned)
                                     if int(lines[j].strip()) > 1)
        return info

    si = vol_info(vol_path)

    # Volume error check (1% threshold for curved, exact for flat)
    ok = True
    if v_exact and abs(verr) > 1.0:
        ok = False
    tag = "OK" if ok else "FAIL"
    if not ok:
        all_pass = False

    print(f"{name:<14} {verr:>+12.6e} {aerr:>+12.6e} {mesh.ne:>7} {si.get('nse','?'):>6} {si.get('nedges','?'):>7} {si.get('n_ho_e','?'):>7} {tag:>8}")

print("-" * len(header))
if all_pass:
    print("ALL SHAPES: volume accuracy OK  (PASS)")
else:
    print("SOME SHAPES FAILED: investigate errors above")
