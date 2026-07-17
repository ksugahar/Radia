"""Generate the 50 mm-workpiece FEM-Kelvin mesh for the ESIM
cross-formulation validation lane (Cubit-embedded Python).

This is the SA-26-070 chapter-5 benchmark workpiece size
(dia. 50 mm x height 25 mm), meshed as a HCurl A + Kelvin domain so
`calc_fem_kelvin.py --impedance esim` can drive it.  It is the
FEM-side counterpart of the BIE-side `samples/ih_bem_sample_p1.vol`
(regenerated from `samples/ih_bem_sample.jou`); the two are
INDEPENDENT surface meshes of the same 50 mm cylinder, so a P_wp
comparison between them isolates the outer-field formulation.

Same coil (gapped torus R=30mm / r=3mm / 355deg sweep), same air
sphere R=60mm, and the same `add_kelvin_cubit` construction as
`samples/ih_fem_kelvin_demo.py`; ONLY the workpiece cylinder differs
(demo is 20mm x 20mm, this is 50mm dia x 25mm).  The wp cavity
surface is sized at 3 mm to match the BIE mesh density.

Run via the Cubit Python API (NOT `play`, which cannot import the
Kelvin helper module).  `make_meshes.py` in this directory drives it;
to run standalone:

    python make_meshes.py --only fem50

Output: <lane>/ih_fem_kelvin_50mm.vol  (gitignored; regenerated on
demand).  Coil STEP is NOT re-exported -- the shared
`samples/ih_fem_kelvin_demo_coil.step` is the identical coil.

The Kelvin helper is imported from the DEPLOYED Cubit plugin
(`<Cubit>/bin/plugins/cubit_helpers/add_kelvin.py`), which is the
self-contained, Cubit-embedded-Python version.  Set
`RADIA_50MM_VOL_OUT` to override the output path.
"""
import os
import sys

import cubit

_HELPERS = os.path.join(
    os.environ.get("CUBIT_PLUGIN_DIR",
                   r"C:\Program Files\Coreform Cubit 2025.12\bin\plugins"),
    "cubit_helpers")
if os.path.isdir(_HELPERS) and _HELPERS not in sys.path:
    sys.path.insert(0, _HELPERS)
try:
    from add_kelvin import add_kelvin_cubit
except ImportError:
    from cubit_mesh_export.cubit_helpers.add_kelvin import add_kelvin_cubit

# `__file__` is undefined when this script is exec()'d by make_meshes.py
# inside a single cubit.init session; RADIA_50MM_VOL_OUT then supplies
# the output path.  Standalone (python make_ih_fem_kelvin_50mm.py under
# a cubit-enabled interpreter) it falls back to the script directory.
try:
    _HERE = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _HERE = os.getcwd()
OUT_VOL = os.environ.get(
    "RADIA_50MM_VOL_OUT",
    os.path.join(_HERE, "ih_fem_kelvin_50mm.vol")).replace("\\", "/")

cubit.cmd("reset")

R_air = 0.060

# === 1. Coil: gapped torus by sweep (identical to demo) ===
cubit.cmd("create vertex 0.030000 0 0")
cubit.cmd("create vertex 0.033000 0 0")
cubit.cmd("create vertex 0.030000 0 0.003000")
cubit.cmd("create curve arc center vertex 1 vertex 2 vertex 3 "
          "normal 0 1 0 full")
cubit.cmd("create surface curve 1")
cubit.cmd("sweep surface 1 axis 0 0 0 0 0 1 angle 355")
coil_vid = cubit.get_last_id("volume")

# === 2. Workpiece cylinder: SA-26-070 chapter-5 size (dia 50 x 25) ===
cubit.cmd("create cylinder radius 0.025 height 0.025")
wp_vid = cubit.get_last_id("volume")

# === 3. Air sphere ===
cubit.cmd("create sphere radius %g" % R_air)
air_vid = cubit.get_last_id("volume")

# === 4. Webcut air sphere at z=0 (copy-mesh seam for Kelvin) ===
id_before = cubit.get_last_id("volume")
cubit.cmd("webcut volume %d with plane zplane" % air_vid)
air_top = air_vid
air_bot = id_before + 1

# === 5. Subtract coil + wp from BOTH air hemispheres ===
cubit.cmd("subtract volume %d %d from volume %d %d keep_tool"
          % (coil_vid, wp_vid, air_top, air_bot))

# === 6. Imprint + merge ===
cubit.cmd("imprint volume %d %d %d %d"
          % (coil_vid, wp_vid, air_top, air_bot))
cubit.cmd("merge volume %d %d %d %d"
          % (coil_vid, wp_vid, air_top, air_bot))

# === 7. Mesh (workpiece NOT meshed -- HOLE policy) ===
cubit.cmd("curve in volume %d with length < 0.05 interval 12" % coil_vid)
cubit.cmd("curve in volume %d with length > 0.05 interval 48" % coil_vid)
cubit.cmd("volume %d scheme tetmesh" % coil_vid)
cubit.cmd("volume %d size 0.003" % coil_vid)
cubit.cmd("mesh volume %d" % coil_vid)

# wp cavity surface at 3 mm to match the BIE mesh density
cubit.cmd("surface in volume %d size 0.003" % wp_vid)
cubit.cmd("volume %d %d scheme tetmesh" % (air_top, air_bot))
cubit.cmd("volume %d %d size 0.008" % (air_top, air_bot))
cubit.cmd("mesh volume %d %d" % (air_top, air_bot))

# === 8. Blocks (NO 'work' block per WP HOLE policy) ===
cubit.cmd("set duplicate block elements off")
cubit.cmd("block 1 add volume %d" % coil_vid)
cubit.cmd('block 1 name "coil"')
cubit.cmd("block 2 add volume %d %d" % (air_top, air_bot))
cubit.cmd('block 2 name "air"')

# === 9. Sidesets ===
cubit.cmd('group "coil_gaps" add surface in volume %d with area < 0.0001'
          % coil_vid)
cubit.cmd('sideset 1 add surface in coil_gaps with y_coord > -0.001')
cubit.cmd('sideset 1 name "source"')
cubit.cmd('sideset 2 add surface in coil_gaps with y_coord < -0.001')
cubit.cmd('sideset 2 name "sink"')
cubit.cmd("sideset 3 add surface in volume %d" % wp_vid)
cubit.cmd('sideset 3 name "sibc"')

# === 10. Kelvin open boundary ===
info = add_kelvin_cubit(R=R_air, symmetry=["z"])

# === 11. Export ===
cubit.cmd('export netgen "%s" order 1 overwrite' % OUT_VOL)

print()
print("make_ih_fem_kelvin_50mm.py: done")
print("  Kelvin centre: (%g, %g, %g)" % info["center"])
print("  Mesh .vol:     %s" % OUT_VOL)
