"""FEM-Kelvin sample: gapped torus coil + cylindrical workpiece + air + Kelvin sphere.

Cubit Python script (playback from Cubit GUI or batch mode).

Kelvin = TWO IDENTICAL SPHERES (same radius R), OFFSET in space.
  - Interior sphere: coil + workpiece + air (physical domain)
  - Exterior sphere: Kelvin domain (air only, material modulated)
  - Periodic BC via TRANSLATION (offset vector)
  - copy mesh surface from interior to exterior (MANDATORY for 1:1 DOFs)

Required labels (consumed by calc_fem_kelvin.py):
  block 1 = "coil"      (gapped torus, source J)
  block 2 = "workpiece" (cylinder)
  block 3 = "air"       (interior sphere with coil+workpiece holes)
  block 4 = "kelvin"    (exterior sphere, offset)
  sideset 1 = "source"  (one gap face on the coil)
  sideset 2 = "sink"    (other gap face on the coil)
  sideset 3 = "sibc"    (workpiece outer surface)
  nodeset 100 = "GND"   (Kelvin sphere center = physical infinity)

Geometry:
  Coil:      gapped torus, R_major = 0.030 m, R_minor = 0.003 m, gap = 5 deg
  Workpiece: cylinder R = 0.025 m, H = 0.025 m
  Sphere R:  0.060 m (both interior and exterior, SAME radius)
  Offset:    0.200 m in x direction (> 3R, no overlap)

Run: Solve -> Radia-NGSolve -> Induction Heating -> Method: FEM
"""
import math
import cubit

cubit.cmd("reset")

R = 0.060       # Kelvin sphere radius (both interior and exterior)
offset_x = 0.200  # Exterior sphere offset in x

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

# === 3. Interior air sphere ==================================================
cubit.cmd(f"create sphere radius {R}")
air_vid = cubit.get_last_id("volume")

# === 4. Exterior (Kelvin) sphere — SAME radius, offset in x =================
cubit.cmd(f"create sphere radius {R}")
kelvin_vid = cubit.get_last_id("volume")
cubit.cmd(f"move volume {kelvin_vid} x {offset_x}")

# === 5. Webcut BOTH spheres with zplane ======================================
# ACIS spheres have 0 curves — webcut creates equator curves needed for copy mesh.
id_before = cubit.get_last_id("volume")
cubit.cmd(f"webcut volume {air_vid} with plane zplane")
air_top = air_vid
air_bot = id_before + 1

id_before2 = cubit.get_last_id("volume")
cubit.cmd(f"webcut volume {kelvin_vid} with plane zplane")
# zplane at z=0 passes through the kelvin sphere center (0.200, 0, 0)
# -> same hemisphere topology as the interior sphere cut.
kelvin_top = kelvin_vid
kelvin_bot = cubit.get_last_id("volume")

# === 6. Subtract coil + workpiece from air ===================================
cubit.cmd(f"subtract volume {coil_vid} {wp_vid} from volume {air_top} {air_bot} keep_tool")

# === 7. Imprint + merge =====================================================
# Imprint air halves with each other (equator), kelvin halves with each other.
# Do NOT imprint interior with exterior (they don't touch).
cubit.cmd(f"imprint volume {air_top} {air_bot}")
cubit.cmd(f"merge volume {air_top} {air_bot}")
cubit.cmd(f"imprint volume {kelvin_top} {kelvin_bot}")
cubit.cmd(f"merge volume {kelvin_top} {kelvin_bot}")
# Also imprint coil/wp with air
cubit.cmd(f"imprint volume {coil_vid} {wp_vid} {air_top} {air_bot}")
cubit.cmd(f"merge volume {coil_vid} {wp_vid} {air_top} {air_bot}")

# === 8. Mesh coil + workpiece + air (NOT kelvin yet) ========================
cubit.cmd(f"curve in volume {coil_vid} with length < 0.05 interval 12")
cubit.cmd(f"curve in volume {coil_vid} with length > 0.05 interval 48")
cubit.cmd(f"volume {coil_vid} scheme tetmesh")
cubit.cmd(f"volume {coil_vid} size 0.003")
cubit.cmd(f"mesh volume {coil_vid}")

cubit.cmd(f"volume {wp_vid} scheme tetmesh")
cubit.cmd(f"volume {wp_vid} size 0.005")
cubit.cmd(f"mesh volume {wp_vid}")

cubit.cmd(f"volume {air_top} {air_bot} scheme tetmesh")
cubit.cmd(f"volume {air_top} {air_bot} size 0.020")
cubit.cmd(f"mesh volume {air_top} {air_bot}")

# === 9. Copy mesh: interior sphere -> exterior sphere ========================
# The interior sphere's outer surface is now meshed (from air meshing).
# Copy that mesh onto the exterior sphere's surface for 1:1 DOF matching.
A_hemi = 2.0 * math.pi * R**2  # hemisphere area for identification

for inner_vid, outer_vid in [(air_top, kelvin_top), (air_bot, kelvin_bot)]:
    # Find hemisphere surface (largest area) in each volume
    inner_surfs = cubit.get_relatives("volume", inner_vid, "surface")
    outer_surfs = cubit.get_relatives("volume", outer_vid, "surface")

    # Interior: the hemisphere (outer boundary of air) is the largest surface
    # that is NOT shared with coil/workpiece
    src = max(inner_surfs, key=lambda s: cubit.surface(s).area())
    dst = max(outer_surfs, key=lambda s: cubit.surface(s).area())

    # Pick longest curve (equator) for orientation reference
    src_curves = cubit.get_relatives("surface", src, "curve")
    dst_curves = cubit.get_relatives("surface", dst, "curve")

    if not src_curves or not dst_curves:
        print(f"WARNING: no curves on surfaces {src}/{dst}, skipping copy mesh")
        continue

    src_c = max(src_curves, key=lambda c: cubit.curve(c).length())
    dst_c = max(dst_curves, key=lambda c: cubit.curve(c).length())
    src_v = cubit.get_relatives("curve", src_c, "vertex")[0]
    dst_v = cubit.get_relatives("curve", dst_c, "vertex")[0]

    print(f"  copy mesh: surface {src} (A={cubit.surface(src).area():.6f}) "
          f"-> surface {dst} (A={cubit.surface(dst).area():.6f})")

    cubit.cmd(
        f"copy mesh surface {src} onto surface {dst} "
        f"source curve {src_c} source vertex {src_v} "
        f"target curve {dst_c} target vertex {dst_v}")

# === 10. Mesh kelvin volumes (boundary already constrained from copy) ========
# Do NOT set volume size — it would override the copied surface mesh.
cubit.cmd(f"volume {kelvin_top} {kelvin_bot} scheme tetmesh")
cubit.cmd(f"mesh volume {kelvin_top} {kelvin_bot}")

# === 11. Blocks ==============================================================
cubit.cmd(f"block 1 add volume {coil_vid}")
cubit.cmd('block 1 name "coil"')
cubit.cmd(f"block 2 add volume {wp_vid}")
cubit.cmd('block 2 name "workpiece"')
cubit.cmd(f"block 3 add volume {air_top} {air_bot}")
cubit.cmd('block 3 name "air"')
cubit.cmd(f"block 4 add volume {kelvin_top} {kelvin_bot}")
cubit.cmd('block 4 name "kelvin"')

# === 12. Sidesets ============================================================
cubit.cmd(f'group "coil_gaps" add surface in volume {coil_vid} with area < 0.0001')
cubit.cmd('sideset 1 add surface in coil_gaps with y_coord > -0.001')
cubit.cmd('sideset 1 name "source"')
cubit.cmd('sideset 2 add surface in coil_gaps with y_coord < -0.001')
cubit.cmd('sideset 2 name "sink"')

cubit.cmd(f"sideset 3 add surface in volume {wp_vid}")
cubit.cmd('sideset 3 name "sibc"')

# Kelvin periodic boundaries: interior sphere outer surface + exterior sphere surface.
# After subtract, the interior sphere outer surface = largest surface in each air half.
# The exterior sphere surface = all surfaces in each kelvin half.
for a_vid in [air_top, air_bot]:
    surfs = cubit.get_relatives("volume", a_vid, "surface")
    # The outer hemisphere surface has area ~ 2*pi*R^2
    outer_s = max(surfs, key=lambda s: cubit.surface(s).area())
    cubit.cmd(f"sideset 4 add surface {outer_s}")
cubit.cmd('sideset 4 name "kelvin_int"')

for k_vid in [kelvin_top, kelvin_bot]:
    surfs = cubit.get_relatives("volume", k_vid, "surface")
    # All surfaces of the kelvin sphere
    for s in surfs:
        cubit.cmd(f"sideset 5 add surface {s}")
cubit.cmd('sideset 5 name "kelvin_ext"')

# === 13. GND nodeset =========================================================
# GND = Kelvin sphere center = physical infinity
cubit.cmd(f"create vertex {offset_x} 0 0")
gnd_vid = cubit.get_last_id("vertex")
cubit.cmd(f"nodeset 100 add vertex {gnd_vid}")
cubit.cmd('nodeset 100 name "GND"')

# Hide air + kelvin
cubit.cmd(f"volume {air_top} {air_bot} visibility off")
cubit.cmd(f"volume {kelvin_top} {kelvin_bot} visibility off")

print("ih_fem_kelvin_sample.py: done (offset Kelvin)")
