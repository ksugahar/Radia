"""Closed (360 deg) torus variant of ih_fem_kelvin_skin.py.

(2) A-V sanity check. On a closed torus T0 is topologically ill-defined
(no valid Dirichlet phi for a single-valued potential). Historical
note: the PEEC Biot-Savart scattered path was retired 2026-04-24; this
demoted sample is kept only for the analytical closed-torus L reference
below. Expected closed-torus L:

  DC (uniform current, GMD = a * exp(-1/4)):
    L = mu_0 * R * (ln(8R/a) - 7/4) = 99.23 nH   (R=30mm, a=3mm)
  HF (skin limit, GMD = a):
    L = mu_0 * R * (ln(8R/a) - 2)   = 89.81 nH
"""
import os
import sys

import cubit

# See README.md in this directory for why this file lives here (demoted
# from panels/samples/ on 2026-04-23). Locate add_kelvin: repo-relative
# first, then pip-installed radia package.
_here = os.path.dirname(os.path.abspath(__file__))
_panels_repo = os.path.normpath(
    os.path.join(_here, "..", "..", "..", "src", "radia", "panels"))
if os.path.isfile(os.path.join(_panels_repo, "add_kelvin.py")):
    _panels_dir = _panels_repo
else:
    try:
        import radia as _radia
        _panels_dir = os.path.join(os.path.dirname(_radia.__file__), "panels")
    except ImportError:
        raise RuntimeError(
            "Cannot find add_kelvin.py.  Install 'radia' (pip install radia) "
            "or run this script from the Radia source tree.")
if _panels_dir not in sys.path:
    sys.path.insert(0, _panels_dir)
from add_kelvin import add_kelvin_cubit

cubit.cmd("reset")

R = 0.060
SAMPLES_DIR = os.getcwd()

cubit.cmd("create vertex 0.030000 0 0")
cubit.cmd("create vertex 0.033000 0 0")
cubit.cmd("create vertex 0.030000 0 0.003000")
cubit.cmd("create curve arc center vertex 1 vertex 2 vertex 3 normal 0 1 0 full")
cubit.cmd("create surface curve 1")
cubit.cmd("sweep surface 1 axis 0 0 0 0 0 1 angle 360")
coil_vid = cubit.get_last_id("volume")

cubit.cmd("create cylinder radius 0.025000 height 0.025000")
wp_vid = cubit.get_last_id("volume")

cubit.cmd("create sphere radius %g" % R)
air_vid = cubit.get_last_id("volume")

id_before = cubit.get_last_id("volume")
cubit.cmd("webcut volume %d with plane zplane" % air_vid)
air_top = air_vid
air_bot = id_before + 1

cubit.cmd("subtract volume %d %d from volume %d %d keep_tool" % (coil_vid, wp_vid, air_top, air_bot))

cubit.cmd("imprint volume %d %d %d %d" % (coil_vid, wp_vid, air_top, air_bot))
cubit.cmd("merge volume %d %d %d %d" % (coil_vid, wp_vid, air_top, air_bot))

step_path = os.path.join(SAMPLES_DIR, "ih_closed_torus_coil.step")
cubit.cmd('export step "%s" volume %d overwrite' % (step_path, coil_vid))
print("STEP exported: %s" % step_path)

# Coil volume is NOT meshed: calc_fem_coilmesh treats it as a hole in
# the air domain, with SIBC Robin BC on the coil_surface sideset.
# Only air + (later) kelvin are meshed with HCurl DOFs.
cubit.cmd("curve in volume %d with length < 0.05 interval 60" % coil_vid)
cubit.cmd("curve in volume %d with length > 0.05 interval 180" % coil_vid)

cubit.cmd("surface in volume %d size 0.002" % wp_vid)
cubit.cmd("volume %d %d scheme tetmesh" % (air_top, air_bot))
cubit.cmd("volume %d %d size 0.020" % (air_top, air_bot))
cubit.cmd("mesh volume %d %d" % (air_top, air_bot))

cubit.cmd("block 2 add volume %d %d" % (air_top, air_bot))
cubit.cmd('block 2 name "air"')

cubit.cmd("sideset 1 add surface in volume %d" % wp_vid)
cubit.cmd('sideset 1 name "sibc"')

# Coil outer surface for calc_fem_coilmesh (coil SIBC + K_phi source).
# Closed torus has no gap faces, so all coil surfaces carry K.
cubit.cmd("sideset 2 add surface in volume %d" % coil_vid)
cubit.cmd('sideset 2 name "coil_surface"')

info = add_kelvin_cubit(R=R, symmetry=["z"])

vol_path = os.path.join(SAMPLES_DIR, "ih_closed_torus.vol")
cubit.cmd('export netgen "%s" order 1 overwrite' % vol_path)

print("")
print("ih_closed_torus.py: done")
print("  Kelvin center: (%g, %g, %g)" % info["center"])
print("  .vol exported:  %s" % vol_path)
print("  .step exported: %s" % step_path)
