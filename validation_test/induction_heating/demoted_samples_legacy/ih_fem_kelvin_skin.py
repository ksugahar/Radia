"""FEM-Kelvin sample with SKIN-DEPTH-RESOLVED coil mesh.

Cubit Python script (playback from Cubit GUI or batch mode).

Same geometry as ih_fem_kelvin_sample.py but with fine coil mesh to
support --coil-sigma 5.8e7 (eddy current in coil volume).

Cu at 7 kHz: delta = 0.79 mm.  Cross-section circles: 60 intervals
(edge ~ 0.31 mm, ~2.5 edges per skin depth).  Interior grades to ~1 mm.

Required labels (consumed by calc_fem_kelvin.py):
  block 1 = "coil"      (gapped torus, source J + optional eddy)
  block 2 = "air"       (interior sphere with coil+workpiece holes)
  + block "kelvin"       (added by add_kelvin_cubit)
  sideset 1 = "source"  (one gap face on the coil)
  sideset 2 = "sink"    (other gap face on the coil)
  sideset 3 = "sibc"    (hole boundary where workpiece was subtracted)
  + sideset "kelvin_int"/"kelvin_ext" (added by add_kelvin_cubit)
  + nodeset 100 = "GND"  (added by add_kelvin_cubit)

HOLE APPROACH: workpiece is subtracted from air and NOT meshed.

Geometry:
  Coil:      gapped torus, R_major = 0.030 m, R_minor = 0.003 m, gap = 5 deg
             L_air = mu_0 * R * (ln(8R/a) - 2) * (1 - gap/360) = 88.55 nH
  Workpiece: cylinder R = 0.025 m, H = 0.025 m
  Air sphere R = 0.060 m (~2x coil major radius)
  Kelvin:    exterior sphere (same R=0.060), offset 0.180 m in +x

Run: Solve -> Radia-NGSolve -> Induction Heating -> Method: FEM
  with --coil-sigma 5.8e7 for FEM truth,
  or  --peec-step .step --peec-nwinc 3 for PEEC comparison.
"""
import math
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

R = 0.060  # Air / Kelvin sphere radius
# __file__ may resolve to the *installed* package directory when Cubit
# plays back via startup.py.  Use CWD (= the directory the user launched
# the script from) as the output location for .vol / .step.
SAMPLES_DIR = os.getcwd()

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
cubit.cmd("subtract volume %d %d from volume %d %d keep_tool" % (coil_vid, wp_vid, air_top, air_bot))

# === 6. Imprint + merge =====================================================
cubit.cmd("imprint volume %d %d %d %d" % (coil_vid, wp_vid, air_top, air_bot))
cubit.cmd("merge volume %d %d %d %d" % (coil_vid, wp_vid, air_top, air_bot))

# === 7. STEP export (coil only, for PEEC comparison) ========================
step_path = os.path.join(SAMPLES_DIR, "ih_fem_kelvin_skin_coil.step")
cubit.cmd('export step "%s" volume %d overwrite' % (step_path, coil_vid))
print("STEP exported: %s" % step_path)

# === 8. Mesh: SKIN-DEPTH-RESOLVED coil ======================================
# Cross-section circles (length ~ 0.019 m < 0.05):
#   60 intervals -> edge ~ 0.31 mm  (delta = 0.79 mm -> ~2.5 edges/delta)
# Toroidal arcs (length ~ 0.186 m > 0.05):
#   180 intervals -> edge ~ 1.03 mm
# Volume interior size 0.001 (1 mm) -> natural grading from 0.3 mm surface

cubit.cmd("curve in volume %d with length < 0.05 interval 60" % coil_vid)
cubit.cmd("curve in volume %d with length > 0.05 interval 180" % coil_vid)
cubit.cmd("volume %d scheme tetmesh" % coil_vid)
cubit.cmd("volume %d size 0.001" % coil_vid)
cubit.cmd("mesh volume %d" % coil_vid)

print("Coil meshed (volume %d)" % coil_vid)

# Workpiece is NOT meshed: SIBC treats it as a hole in the air domain.

# Air mesh: fine near workpiece hole (2mm gap to coil) + coarse far field.
# The SIBC surface (workpiece hole boundary) needs fine elements so that
# the PEEC reduced-A formulation can resolve the Biot-Savart A_s field.
cubit.cmd("surface in volume %d size 0.002" % wp_vid)
cubit.cmd("volume %d %d scheme tetmesh" % (air_top, air_bot))
cubit.cmd("volume %d %d size 0.020" % (air_top, air_bot))
cubit.cmd("mesh volume %d %d" % (air_top, air_bot))

# === 9. Blocks (NO workpiece block) =========================================
cubit.cmd("block 1 add volume %d" % coil_vid)
cubit.cmd('block 1 name "coil"')
cubit.cmd("block 2 add volume %d %d" % (air_top, air_bot))
cubit.cmd('block 2 name "air"')

# === 10. Sidesets ============================================================
cubit.cmd('group "coil_gaps" add surface in volume %d with area < 0.0001' % coil_vid)
cubit.cmd('sideset 1 add surface in coil_gaps with y_coord > -0.001')
cubit.cmd('sideset 1 name "source"')
cubit.cmd('sideset 2 add surface in coil_gaps with y_coord < -0.001')
cubit.cmd('sideset 2 name "sink"')

# Workpiece hole boundary (SIBC Robin BC on air-side surface).
cubit.cmd("sideset 3 add surface in volume %d" % wp_vid)
cubit.cmd('sideset 3 name "sibc"')

# === 11. Add Kelvin open boundary (ONE LINE) =================================
info = add_kelvin_cubit(R=R, symmetry=["z"])

# === 12. Export .vol =========================================================
vol_path = os.path.join(SAMPLES_DIR, "ih_fem_kelvin_skin.vol")
cubit.cmd('export netgen "%s" order 1 overwrite' % vol_path)

print("")
print("ih_fem_kelvin_skin.py: done")
print("  Kelvin center: (%g, %g, %g)" % info["center"])
print("  .vol exported:  %s" % vol_path)
print("  .step exported: %s" % step_path)
print("")
print("Next: compare_fem_vs_peec_skin.py --vol '%s' --step '%s'" % (vol_path, step_path))
