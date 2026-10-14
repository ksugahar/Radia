"""FEM-Kelvin sample: gapped torus coil + cylindrical workpiece + air + Kelvin.

Cubit Python script (playback from Cubit GUI or batch mode).

This script creates the physical domain (coil + workpiece + air sphere),
then calls add_kelvin_cubit() to add the Kelvin open-boundary sphere
automatically.  The user does NOT need to create the Kelvin geometry.

Required labels (consumed by calc_fem_kelvin.py):
  block 1 = "coil"      (gapped torus, source J)
  block 2 = "air"       (interior sphere with coil+workpiece holes)
  + block "kelvin"       (added by add_kelvin_cubit)
  sideset 1 = "source"  (one gap face on the coil)
  sideset 2 = "sink"    (other gap face on the coil)
  sideset 3 = "sibc"    (hole boundary where workpiece was subtracted)
  + sideset "kelvin_int"/"kelvin_ext" (added by add_kelvin_cubit)
  + nodeset 100 = "GND"  (added by add_kelvin_cubit)

HOLE APPROACH: workpiece is subtracted from air and NOT meshed.
The SIBC Robin BC is applied on the air-side hole boundary.

Geometry:
  Coil:      gapped torus, R_major = 0.030 m, R_minor = 0.003 m, gap = 5 deg
  Workpiece: cylinder R = 0.025 m, H = 0.025 m
  Air sphere R = 0.060 m (~2x coil major radius)
  Kelvin:    exterior sphere (same R=0.060), offset 0.180 m in +x

Run: Solve -> Radia-NGSolve -> Induction Heating -> Method: FEM
"""
import math
import os
import sys

import cubit

# --- Locate add_kelvin.py from the Radia panels folder. ------------------
# This file was demoted from src/radia/panels/samples/ on 2026-04-23 and
# preserved here on 2026-06-29 (see README.md). The loader below
# tries two locations in order:
#   1. Repo-relative (when running from a checked-out Radia source tree)
#   2. pip-installed radia package (fallback for non-repo environments)
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

R = 0.060  # Air / Kelvin sphere radius

# === 1. Coil: gapped torus by sweep =========================================
cubit.cmd("create vertex 0.030000 0 0")
cubit.cmd("create vertex 0.033000 0 0")
cubit.cmd("create vertex 0.030000 0 0.003000")
cubit.cmd("create curve arc center vertex 1 vertex 2 vertex 3 normal 0 1 0 full")
cubit.cmd("create surface curve 1")
cubit.cmd("sweep surface 1 axis 0 0 0 0 0 1 angle 355")
coil_vid = cubit.get_last_id("volume")

# === 2. Workpiece cylinder ===================================================
cubit.cmd("create cylinder radius 0.025000 height 0.025000")
wp_vid = cubit.get_last_id("volume")

# === 3. Air sphere ===========================================================
cubit.cmd("create sphere radius %g" % R)
air_vid = cubit.get_last_id("volume")

# === 4. Webcut air sphere (equator curves for copy mesh) =====================
id_before = cubit.get_last_id("volume")
cubit.cmd("webcut volume %d with plane zplane" % air_vid)
air_top = air_vid
air_bot = id_before + 1

# === 5. Subtract coil + workpiece from air ===================================
cubit.cmd("subtract volume %d %d from volume %d %d keep_tool"
          % (coil_vid, wp_vid, air_top, air_bot))

# === 6. Imprint + merge =====================================================
cubit.cmd("imprint volume %d %d %d %d" % (coil_vid, wp_vid, air_top, air_bot))
cubit.cmd("merge volume %d %d %d %d" % (coil_vid, wp_vid, air_top, air_bot))

# === 7. Mesh (workpiece NOT meshed — hole approach) =========================
cubit.cmd("curve in volume %d with length < 0.05 interval 12" % coil_vid)
cubit.cmd("curve in volume %d with length > 0.05 interval 48" % coil_vid)
cubit.cmd("volume %d scheme tetmesh" % coil_vid)
cubit.cmd("volume %d size 0.003" % coil_vid)
cubit.cmd("mesh volume %d" % coil_vid)

# Workpiece is NOT meshed: SIBC treats it as a hole in the air domain.

cubit.cmd("volume %d %d scheme tetmesh" % (air_top, air_bot))
cubit.cmd("volume %d %d size 0.020" % (air_top, air_bot))
cubit.cmd("mesh volume %d %d" % (air_top, air_bot))

# === 8. Blocks (NO workpiece block) =========================================
cubit.cmd("block 1 add volume %d" % coil_vid)
cubit.cmd('block 1 name "coil"')
cubit.cmd("block 2 add volume %d %d" % (air_top, air_bot))
cubit.cmd('block 2 name "air"')

# === 9. Sidesets =============================================================
cubit.cmd('group "coil_gaps" add surface in volume %d with area < 0.0001'
          % coil_vid)
cubit.cmd('sideset 1 add surface in coil_gaps with y_coord > -0.001')
cubit.cmd('sideset 1 name "source"')
cubit.cmd('sideset 2 add surface in coil_gaps with y_coord < -0.001')
cubit.cmd('sideset 2 name "sink"')

# Workpiece hole boundary (SIBC Robin BC on air-side surface).
# After subtract + imprint/merge, surfaces of wp_vid are shared with air.
cubit.cmd("sideset 3 add surface in volume %d" % wp_vid)
cubit.cmd('sideset 3 name "sibc"')

# === 10. Add Kelvin open boundary (ONE LINE) =================================
info = add_kelvin_cubit(R=R, symmetry=["z"])

print("")
print("ih_fem_kelvin_sample.py: done")
print("  Kelvin center: (%g, %g, %g)" % info["center"])
print("  Next: export netgen \"ih_fem_kelvin_sample.vol\" order 1 overwrite")
