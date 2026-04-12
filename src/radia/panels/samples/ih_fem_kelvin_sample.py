"""FEM-Kelvin sample: gapped torus coil + cylindrical workpiece + air + Kelvin shell.

Cubit Python script (playback from Cubit GUI or batch mode).
Uses 'copy mesh surface' to ensure 1:1 triangle correspondence between
kelvin_int and kelvin_ext for NGSolve periodic BC.

Required labels (consumed by calc_fem_kelvin.py):
  block 1 = "coil"      (gapped torus, source J)
  block 2 = "workpiece" (cylinder)
  block 3 = "air"       (carved inner air half-spheres)
  block 4 = "kelvin"    (outer hemispherical shells, mesh-copied)
  sideset 1 = "source"  (one gap face on the coil)
  sideset 2 = "sink"    (other gap face on the coil)
  sideset 3 = "sibc"    (workpiece outer surface)
  nodeset 100 = "GND"   (origin vertex)

Geometry:
  Coil:      gapped torus, R_major = 0.030 m, R_minor = 0.003 m, gap = 5 deg
  Workpiece: cylinder R = 0.025 m, H = 0.025 m
  Air:       sphere R = 0.060 m
  Kelvin:    shell R = 0.060 m to R = 0.120 m

Run: Solve -> Radia-NGSolve -> Induction Heating -> Method: FEM
"""
import math
import cubit

cubit.cmd("reset")

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

# === 3. Air sphere (inner) ===================================================
R_inner = 0.060
cubit.cmd(f"create sphere radius {R_inner}")
air_vid = cubit.get_last_id("volume")

# === 4. Kelvin outer sphere ==================================================
R_outer = 0.120
cubit.cmd(f"create sphere radius {R_outer}")
kelvin_vid = cubit.get_last_id("volume")

# === 5. Webcut air + kelvin with zplane ======================================
id_before = cubit.get_last_id("volume")
cubit.cmd(f"webcut volume {air_vid} {kelvin_vid} with plane zplane")
air_top = air_vid
air_bot = id_before + 1
kelvin_top = kelvin_vid
kelvin_bot = id_before + 2

# === 6. Subtract air from kelvin -> shells ===================================
cubit.cmd(f"subtract volume {air_top} from volume {kelvin_top} keep_tool")
cubit.cmd(f"subtract volume {air_bot} from volume {kelvin_bot} keep_tool")

# === 7. Subtract coil + workpiece from air ===================================
cubit.cmd(f"subtract volume {coil_vid} {wp_vid} from volume {air_top} {air_bot} keep_tool")

# === 8. Imprint + merge =====================================================
cubit.cmd("imprint all")
cubit.cmd("merge all")

# === 9. Mesh coil + workpiece + air (NOT kelvin yet) ========================
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

# === 10. Copy mesh: kelvin inner -> kelvin outer =============================
# After meshing air, the shared air|kelvin interface is triangulated.
# Copy that triangulation onto the kelvin outer surface so periodic BC
# has exact 1:1 node correspondence.
#
# For each hemisphere pair (top, bottom):
#   - src_surf = kelvin inner hemi (shared with air, area ~ 2*pi*R_inner^2)
#   - dst_surf = kelvin outer hemi (area ~ 2*pi*R_outer^2)
# Identify by area: inner hemi is smaller, outer hemi is larger.
# Equator flat face has area ~ pi*(R_outer^2 - R_inner^2).

A_inner_hemi = 2.0 * math.pi * R_inner**2  # ~ 0.0226
A_outer_hemi = 2.0 * math.pi * R_outer**2  # ~ 0.0905
A_equator = math.pi * (R_outer**2 - R_inner**2)  # ~ 0.0339

for k_vid in [kelvin_top, kelvin_bot]:
    surfs = cubit.get_relatives("volume", k_vid, "surface")
    # Classify surfaces by area
    surf_areas = [(s, cubit.surface(s).area()) for s in surfs]
    surf_areas.sort(key=lambda x: x[1])

    # Find inner (closest to A_inner_hemi) and outer (closest to A_outer_hemi)
    src_surf = min(surf_areas, key=lambda x: abs(x[1] - A_inner_hemi))[0]
    dst_surf = min(surf_areas, key=lambda x: abs(x[1] - A_outer_hemi))[0]

    # Mesh the source surface first (it should already be meshed from air)
    # Copy mesh requires source curve + vertex for orientation mapping
    src_curves = cubit.get_relatives("surface", src_surf, "curve")
    dst_curves = cubit.get_relatives("surface", dst_surf, "curve")

    if not src_curves or not dst_curves:
        print(f"WARNING: no curves on kelvin surfaces {src_surf}/{dst_surf}")
        continue

    # Pick the longest curve (equator circle) for orientation reference
    src_c = max(src_curves, key=lambda c: cubit.curve(c).length())
    dst_c = max(dst_curves, key=lambda c: cubit.curve(c).length())

    src_verts = cubit.get_relatives("curve", src_c, "vertex")
    dst_verts = cubit.get_relatives("curve", dst_c, "vertex")

    if not src_verts or not dst_verts:
        print(f"WARNING: no vertices on curves {src_c}/{dst_c}")
        continue

    src_v = src_verts[0]
    dst_v = dst_verts[0]

    print(f"  copy mesh: surface {src_surf} (A={cubit.surface(src_surf).area():.6f}) "
          f"-> surface {dst_surf} (A={cubit.surface(dst_surf).area():.6f})")
    print(f"    src curve {src_c} vertex {src_v} -> dst curve {dst_c} vertex {dst_v}")

    cubit.cmd(
        f"copy mesh surface {src_surf} onto surface {dst_surf} "
        f"source curve {src_c} source vertex {src_v} "
        f"target curve {dst_c} target vertex {dst_v}")

# === 11. Mesh kelvin volumes (boundary already constrained) ==================
# Do NOT set volume size — it would override the copied surface mesh.
cubit.cmd(f"volume {kelvin_top} {kelvin_bot} scheme tetmesh")
cubit.cmd(f"mesh volume {kelvin_top} {kelvin_bot}")

# === 12. Blocks ==============================================================
cubit.cmd(f"block 1 add volume {coil_vid}")
cubit.cmd('block 1 name "coil"')
cubit.cmd(f"block 2 add volume {wp_vid}")
cubit.cmd('block 2 name "workpiece"')
cubit.cmd(f"block 3 add volume {air_top} {air_bot}")
cubit.cmd('block 3 name "air"')
cubit.cmd(f"block 4 add volume {kelvin_top} {kelvin_bot}")
cubit.cmd('block 4 name "kelvin"')

# === 13. Sidesets ============================================================
cubit.cmd(f'group "coil_gaps" add surface in volume {coil_vid} with area < 0.0001')
cubit.cmd('sideset 1 add surface in coil_gaps with y_coord > -0.001')
cubit.cmd('sideset 1 name "source"')
cubit.cmd('sideset 2 add surface in coil_gaps with y_coord < -0.001')
cubit.cmd('sideset 2 name "sink"')

cubit.cmd(f"sideset 3 add surface in volume {wp_vid}")
cubit.cmd('sideset 3 name "sibc"')

# === 14. GND nodeset =========================================================
cubit.cmd("create vertex 0 0 0")
gnd_vid = cubit.get_last_id("vertex")
cubit.cmd(f"nodeset 100 add vertex {gnd_vid}")
cubit.cmd('nodeset 100 name "GND"')

# Hide air + kelvin
cubit.cmd(f"volume {air_top} {air_bot} visibility off")
cubit.cmd(f"volume {kelvin_top} {kelvin_bot} visibility off")

print("ih_fem_kelvin_sample.py: done")
