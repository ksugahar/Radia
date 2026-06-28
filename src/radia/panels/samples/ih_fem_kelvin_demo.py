"""PEEC + FEM-Kelvin demo sample (self-contained).

This is the *canonical end-to-end working* PEEC-FEM-Kelvin sample for
the IH panel.  Unlike the `*.jou` siblings, this `.py` script EXPLICITLY
calls `add_kelvin_cubit()` so the exported `.vol` is guaranteed to
contain the `kelvin` material block and the `kelvin_int` / `kelvin_ext`
sidesets that `calc_fem_kelvin.py --require-kelvin` insists on.

Why a `.py` and not a `.jou`: Cubit journal files cannot import
Python modules.  The `add_kelvin_cubit()` helper that builds the
Kelvin shell + 1:1 copy-mesh + GND nodeset is a Python function in
`cubit-mesh-export`, so any sample that needs Kelvin must be a `.py`
script that runs in Cubit's embedded Python.  See
`validation_test/induction_heating/demoted_samples_legacy/README.md` for the
.jou-vs-.py rationale (2026-04-23 demotion of confusing .jou+.py
pairs).

Geometry (DELIBERATELY SMALL for fast end-to-end demo):
  Coil:        gapped torus, R_major = 0.030 m, R_minor = 0.003 m,
               gap angle = 5 deg (sweep 355 deg)
  Workpiece:   Cu cylinder, R = 0.010 m, H = 0.020 m (axis along z)
  Air sphere:  R = 0.060 m (~ 2x coil major radius), webcut at z=0
  Kelvin:      added by add_kelvin_cubit(R=0.060, symmetry=["z"]) --
               appears at offset (3 R, 0, 0) by default

Mesh sizing (coarse — fast):
  Coil cross-section circle: 12 intervals (~0.5 mm edge)
  Coil toroidal arc:         48 intervals (~3 mm edge)
  Coil volume:               3 mm tetmesh
  Air sphere hemispheres:    8 mm tetmesh
  Kelvin shell hemispheres:  inherited from air via 1:1 copy-mesh

Total tet count ~ 30k-50k.  Mesh + add_kelvin + export ~ 1-2 min on
typical desktop.  Trade off mesh fineness for runtime: this is the
"is the pipeline wired correctly" demo, NOT a converged production
result -- for production use `ih_fem_kelvin_skin_fine.jou` and add
the launcher's add_kelvin step manually.

Workflow once this script is played in Cubit (GUI or batch):

    # Step 1: play this script in Cubit (builds + add_kelvin + ready for export)
    play "ih_fem_kelvin_demo.py"

    # Step 2: export .vol (Cubit APREPRO command from the Radia plugin)
    export netgen "ih_fem_kelvin_demo.vol" order 1 overwrite

    # Step 3: outside Cubit -- PEEC coil L+R from the .step.
    # Use the canonical `_peec_ind.json` output suffix so the .gitignore
    # rule `*_peec_ind.json` keeps the result file untracked.
    python -m radia.panels.calc_inductance \
        --coil-step ih_fem_kelvin_demo_coil.step \
        --coil-solver peec \
        --frequency 50000 \
        --coil-only \
        --output ih_fem_kelvin_demo_peec_ind.json

    # Step 4: FEM-Kelvin workpiece P_wp + R_wp from the .vol.
    # `_fem_kelvin.json` suffix matches the .gitignore convention.
    # NOTE: --solver pardiso is REQUIRED for Periodic Kelvin BC.
    # The default `auto` selects `ams` which raises
    # "ams is not supported with Periodic Kelvin BC", and `bddc`
    # produces a NaN HCurl solution on this geometry (see
    # memory project_calc_fem_kelvin_bddc_periodic_nan_2026_05_13).
    python -m radia.panels.calc_fem_kelvin \
        --vol ih_fem_kelvin_demo.vol \
        --peec-step ih_fem_kelvin_demo_coil.step \
        --frequency 50000 \
        --material copper \
        --solver pardiso \
        --output ih_fem_kelvin_demo_fem_kelvin.json

Steps 3 and 4 produce the two halves of the PEEC + FEM-Kelvin result:
PEEC gives L_coil + R_coil_AC, FEM-Kelvin gives the workpiece eddy
current dissipation + induced L back-reaction.  Combine externally
for full circuit parameters.

Required labels (consumed by calc_fem_kelvin.py):
  block 1 = "coil"      (gapped torus, source J)
  block 2 = "air"       (interior sphere with coil+workpiece holes)
  block (added) = "kelvin"      -- by add_kelvin_cubit
  sideset 1 = "source"  (one gap face on the coil)
  sideset 2 = "sink"    (other gap face on the coil)
  sideset 3 = "sibc"    (hole boundary where workpiece was subtracted)
  sideset (added) = "kelvin_int", "kelvin_ext"  -- by add_kelvin_cubit
  nodeset (added) = 100 (GND vertex at Kelvin centre)  -- by add_kelvin_cubit
"""
import math
import os
import sys

import cubit

# --- Locate add_kelvin_cubit ----------------------------------------------
# Two resolution paths (matching the demoted-samples loader pattern):
#   1. Repo-relative `src/radia/panels/` (development tree)
#   2. Installed `radia.panels.add_kelvin` shim (pip-installed PyPI radia)
# In Cubit-embedded Python, `radia` is usually NOT importable, so we
# fall back to the deployed `<Cubit>/bin/plugins/cubit_helpers/` path
# which the .ccm command prepends to sys.path at panel startup.
#
# `__file__` is set when run as `python script.py`, but NOT when
# Cubit's `play` (or `coreform_cubit -batch script.py`) executes the
# file via `exec()` -- those run with `__name__ == '__main__'` but no
# `__file__`.  Fall back to CWD which Cubit sets to the script's dir.
try:
    _here = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _here = os.getcwd()
# This file lives at <repo>/src/radia/panels/samples/, so the panels
# directory is exactly its parent (one level up).
_panels_repo = os.path.dirname(_here)
if os.path.isfile(os.path.join(_panels_repo, "add_kelvin.py")):
    if _panels_repo not in sys.path:
        sys.path.insert(0, _panels_repo)
# Also try the deployed Cubit plugin location -- the .ccm prepends
# this to sys.path at panel startup, but in headless batch we may
# have started without that path being injected.
_cubit_plugin_helpers = os.path.join(
    os.environ.get("CUBIT_PLUGIN_DIR",
                   r"C:\Program Files\Coreform Cubit 2025.12\bin\plugins"),
    "cubit_helpers")
if os.path.isdir(_cubit_plugin_helpers) and _cubit_plugin_helpers not in sys.path:
    sys.path.insert(0, _cubit_plugin_helpers)
try:
    from add_kelvin import add_kelvin_cubit
except ImportError:
    # Last-resort fallback for pip-installed cubit-mesh-export
    from cubit_mesh_export.cubit_helpers.add_kelvin import add_kelvin_cubit


cubit.cmd("reset")

R_air = 0.060   # Air / Kelvin sphere radius [m]

# === 1. Coil: gapped torus by sweep =========================================
# Major R = 30 mm, minor r = 3 mm, gap = 5 deg (sweep 355 deg).
cubit.cmd("create vertex 0.030000 0 0")
cubit.cmd("create vertex 0.033000 0 0")
cubit.cmd("create vertex 0.030000 0 0.003000")
cubit.cmd("create curve arc center vertex 1 vertex 2 vertex 3 "
          "normal 0 1 0 full")
cubit.cmd("create surface curve 1")
cubit.cmd("sweep surface 1 axis 0 0 0 0 0 1 angle 355")
coil_vid = cubit.get_last_id("volume")

# === 2. Workpiece cylinder =================================================
# Cu billet centered at origin, axis along z, R=10 mm, H=20 mm.
cubit.cmd("create cylinder radius 0.010 height 0.020")
wp_vid = cubit.get_last_id("volume")

# === 3. Air sphere ==========================================================
# R = 60 mm (~ 2x coil major radius -- enough clearance for the SIBC
# hole boundary at the wp surface to feel like it's in free space).
cubit.cmd("create sphere radius %g" % R_air)
air_vid = cubit.get_last_id("volume")

# === 4. Webcut air sphere at z=0 ===========================================
# Mesh-seam circle on the outer surface -- required by add_kelvin_cubit's
# copy-mesh step (it needs source-curve / target-curve anchors).  Same
# pattern as ih_fem_kelvin_skin_fine.jou and em_1-1_full.jou.
id_before = cubit.get_last_id("volume")
cubit.cmd("webcut volume %d with plane zplane" % air_vid)
air_top = air_vid
air_bot = id_before + 1

# === 5. Subtract coil + wp from BOTH air hemispheres =======================
# `keep_tool` retains the coil + wp at their original IDs and modifies
# the air bodies (carves the coil + wp shapes out of the air).
cubit.cmd("subtract volume %d %d from volume %d %d keep_tool"
          % (coil_vid, wp_vid, air_top, air_bot))

# === 6. Imprint + merge ====================================================
# Conformal mesh across the z=0 seam AND around the coil + wp cavities.
cubit.cmd("imprint volume %d %d %d %d"
          % (coil_vid, wp_vid, air_top, air_bot))
cubit.cmd("merge volume %d %d %d %d"
          % (coil_vid, wp_vid, air_top, air_bot))

# === 7. Coil-only STEP for the PEEC step ===================================
# `calc_inductance.py --step` reads this to extract filaments + L_coil.
out_step = os.path.join(_here, "ih_fem_kelvin_demo_coil.step")
cubit.cmd('export step "%s" volume %d overwrite' % (out_step, coil_vid))

# === 8. Mesh (workpiece NOT meshed -- HOLE policy) ==========================
# Coil: 12 intervals on cross-section circles + 48 on toroidal arcs
# gives a ~3 mm tet target.  Air: 8 mm bulk.  Total ~30k-50k tets.
cubit.cmd("curve in volume %d with length < 0.05 interval 12" % coil_vid)
cubit.cmd("curve in volume %d with length > 0.05 interval 48" % coil_vid)
cubit.cmd("volume %d scheme tetmesh" % coil_vid)
cubit.cmd("volume %d size 0.003" % coil_vid)
cubit.cmd("mesh volume %d" % coil_vid)

# Workpiece is intentionally NOT meshed -- SIBC treats it as a HOLE in
# the air domain (see CLAUDE.md "WP HOLE policy").
cubit.cmd("volume %d %d scheme tetmesh" % (air_top, air_bot))
cubit.cmd("volume %d %d size 0.008" % (air_top, air_bot))
cubit.cmd("mesh volume %d %d" % (air_top, air_bot))

# === 9. Blocks (NO 'work' block per WP HOLE policy) =========================
cubit.cmd("set duplicate block elements off")
cubit.cmd("block 1 add volume %d" % coil_vid)
cubit.cmd('block 1 name "coil"')
cubit.cmd("block 2 add volume %d %d" % (air_top, air_bot))
cubit.cmd('block 2 name "air"')

# === 10. Sidesets ===========================================================
# source / sink: the two coil gap faces (5-degree wedge at y ~ 0).  We
# pick them by area (small cap face = pi * r_minor^2 ~ 2.83e-5 m^2)
# and y sign (one positive, one negative).  Cubit `with` predicates are
# evaluated at parse time and immune to surface ID renumbering (per
# CLAUDE.md "Journal File Portability Policy").
cubit.cmd('group "coil_gaps" add surface in volume %d with area < 0.0001'
          % coil_vid)
cubit.cmd('sideset 1 add surface in coil_gaps with y_coord > -0.001')
cubit.cmd('sideset 1 name "source"')
cubit.cmd('sideset 2 add surface in coil_gaps with y_coord < -0.001')
cubit.cmd('sideset 2 name "sink"')

# sibc: the air-side surface of the wp cavity (= every face of the
# subtracted wp volume).  After subtract + imprint + merge, these
# faces are shared between wp_vid and air_top/air_bot; tagging on
# wp_vid propagates correctly to the BND triangles exported with
# the air block (verified by the legacy validation fixture
# ih_fem_kelvin_sample.py, which uses this exact predicate).
cubit.cmd("sideset 3 add surface in volume %d" % wp_vid)
cubit.cmd('sideset 3 name "sibc"')

# === 11. Add Kelvin open boundary (one call) ===============================
# symmetry=["z"] = full-sphere geometry with the z=0 mesh seam retained
# (NOT a domain-reduction; both hemispheres are kept in the mesh).
# The Kelvin shell is created at offset (3*R, 0, 0) by default (offset
# direction auto-detected from the symmetry list).  Inside the helper:
#   - fresh sphere R = 0.060 m at (0.180, 0, 0)
#   - webcut at z=0 (matching the air mesh seam)
#   - 1:1 copy-mesh from air outer surface to Kelvin shell inner surface
#   - blocks "kelvin", sidesets "kelvin_int" / "kelvin_ext"
#   - nodeset 100 = "GND" at Kelvin centre (anchors the Omega gauge)
info = add_kelvin_cubit(R=R_air, symmetry=["z"])

# === 12. Export .vol with the Radia Cubit plugin =====================
# This bakes blocks (coil + air + kelvin), sidesets (source / sink /
# sibc / kelvin_int / kelvin_ext), and the Kelvin identifications into
# a self-contained .vol that `calc_fem_kelvin.py --require-kelvin`
# can consume directly.  Order 1 is enough for the demo (the panel
# tests are at order 1; higher orders are for production p-convergence).
out_vol = os.path.join(_here, "ih_fem_kelvin_demo.vol").replace("\\", "/")
cubit.cmd('export netgen "%s" order 1 overwrite' % out_vol)

print()
print("ih_fem_kelvin_demo.py: done")
print("  Kelvin centre: (%g, %g, %g)" % info["center"])
print("  Coil STEP:     %s" % out_step)
print("  Mesh .vol:     %s" % out_vol)
print()
print("Next steps (run outside Cubit, in radia-enabled Python 3.12):")
print("  python -m radia.panels.calc_inductance "
      "--coil-step ih_fem_kelvin_demo_coil.step "
      "--coil-solver peec --frequency 50000 --coil-only "
      "--output ih_fem_kelvin_demo_peec_ind.json")
print("  python -m radia.panels.calc_fem_kelvin "
      "--vol ih_fem_kelvin_demo.vol "
      "--peec-step ih_fem_kelvin_demo_coil.step "
      "--frequency 50000 --material copper --solver pardiso "
      "--output ih_fem_kelvin_demo_fem_kelvin.json")
print()
print("(--solver pardiso is REQUIRED for Periodic Kelvin BC: auto -> ams")
print(" raises, and bddc produces a NaN HCurl solution on this geometry.)")
