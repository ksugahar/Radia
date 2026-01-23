#!/usr/bin/env python
"""
Verify Coil Geometry: Racetrack Coil

Creates racetrack coil geometry using Netgen OCC for visualization.
This is Step 0 before mesh generation.

Usage:
  python verify_coil.py           # Export to STEP
  py -i verify_coil.py            # Interactive mode

Output:
  coil.step   : STEP file for ParaView/CAD viewing
"""

import sys
import os
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))

# Use custom NGSolve with GUI support and patches (PR#231, PR#232)
NGSOLVE_PATH = r'S:\NGSolve\01_GitHub\install_ksugahar\lib\site-packages'
sys.path.insert(0, NGSOLVE_PATH)

from netgen.occ import *

# ============================================================
# Coil parameters from ELF_magic.mei
# ============================================================
# AA O 0 131.25 0
# AA SQRING 60 72.5 105 35 5 3
#
# ELF SQRING: lx=60, ly=72.5, lz=105, lt=35, r=5, Mr=3
# - lx=60: inner opening X (mm) - full width, not half
# - ly=72.5: inner opening Y (mm) - full width, not half
# - lz=105: coil height Z (mm)
# - lt=35: conductor thickness (mm)
# - r=5: inner corner radius (mm)
# - Mr=3: corner divisions (not used in Radia 'auto' mode)
#
# ELF coordinate system (from ELF_magic.mei):
# - IMA X: mirror about X axis (Y -> -Y)
# - IMA -Z: mirror about -Z axis (Z -> -Z)
# - GSC 0.001: scale factor (mm -> m)

# Coil center from ELF (AA O 0 131.25 0)
# ELF reads the Cubit-transformed yoke, so ELF coordinates = post-transform coordinates
# The coil wraps around the pole piece at Y = 131.25 mm
cx_mm, cy_mm, cz_mm = 0, 131.25, 0  # mm (ELF specification)

# SQRING parameters (mm)
lx_mm = 60      # inner opening X (full width)
ly_mm = 72.5    # inner opening Y (full width)
lz_mm = 105     # coil height Z
lt_mm = 35      # conductor thickness
r_mm = 5        # inner corner radius

# Convert to meters
cx, cy, cz = cx_mm / 1000, cy_mm / 1000, cz_mm / 1000  # m
lx, ly = lx_mm / 1000, ly_mm / 1000  # m (full widths)
h = lz_mm / 1000  # m
conductor_width = lt_mm / 1000  # m
r_inner = r_mm / 1000  # m
r_outer = r_inner + conductor_width  # m

print("=" * 60)
print("VERIFY COIL GEOMETRY: RACETRACK COIL")
print("=" * 60)
print()
print("ELF_magic.mei parameters:")
print(f"  AA O {cx_mm} {cy_mm} {cz_mm}")
print(f"  AA SQRING {lx_mm} {ly_mm} {lz_mm} {lt_mm} {r_mm} 3")
print()
print("Coil parameters (m):")
print(f"  Center: ({cx:.4f}, {cy:.4f}, {cz:.4f})")
print(f"  Inner opening: {lx:.4f} x {ly:.4f}")
print(f"  Coil height: {h:.4f}")
print(f"  Conductor thickness: {conductor_width:.4f}")
print(f"  Inner corner radius: {r_inner:.4f}")
print(f"  Outer corner radius: {r_outer:.4f}")
print()

# ============================================================
# Step 1: Create racetrack profile (2D cross-section in XY plane)
# ============================================================
print("Step 1: Create racetrack profile")

# Racetrack = rounded rectangle (straight sections + arc corners)
# Inner boundary: lx x ly with r_inner corner radius
# Outer boundary: (lx + 2*lt) x (ly + 2*lt) with r_outer corner radius

# Method: Create 2D rounded rectangles and extrude
# Use WorkPlane to create 2D profiles

# Outer profile: half-widths are (lx/2 + lt) and (ly/2 + lt)
outer_half_x = lx / 2 + conductor_width
outer_half_y = ly / 2 + conductor_width

# Inner profile: half-widths are lx/2 and ly/2
inner_half_x = lx / 2
inner_half_y = ly / 2

print(f"  Outer half-size: {outer_half_x*1000:.1f} x {outer_half_y*1000:.1f} mm")
print(f"  Inner half-size: {inner_half_x*1000:.1f} x {inner_half_y*1000:.1f} mm")
print(f"  Outer full-size: {outer_half_x*2*1000:.1f} x {outer_half_y*2*1000:.1f} mm")
print(f"  Inner full-size: {inner_half_x*2*1000:.1f} x {inner_half_y*2*1000:.1f} mm")

# Create outer rounded rectangle using WorkPlane
# Start at bottom-left corner (after radius) and go counter-clockwise
wp_outer = WorkPlane()
wp_outer.MoveTo(-outer_half_x + r_outer, -outer_half_y)
wp_outer.LineTo(outer_half_x - r_outer, -outer_half_y)
wp_outer.Arc(r_outer, 90)  # Bottom-right corner
wp_outer.LineTo(outer_half_x, outer_half_y - r_outer)
wp_outer.Arc(r_outer, 90)  # Top-right corner
wp_outer.LineTo(-outer_half_x + r_outer, outer_half_y)
wp_outer.Arc(r_outer, 90)  # Top-left corner
wp_outer.LineTo(-outer_half_x, -outer_half_y + r_outer)
wp_outer.Arc(r_outer, 90)  # Bottom-left corner
wp_outer.Close()
outer_face = wp_outer.Face()

# Create inner rounded rectangle
wp_inner = WorkPlane()
wp_inner.MoveTo(-inner_half_x + r_inner, -inner_half_y)
wp_inner.LineTo(inner_half_x - r_inner, -inner_half_y)
wp_inner.Arc(r_inner, 90)  # Bottom-right corner
wp_inner.LineTo(inner_half_x, inner_half_y - r_inner)
wp_inner.Arc(r_inner, 90)  # Top-right corner
wp_inner.LineTo(-inner_half_x + r_inner, inner_half_y)
wp_inner.Arc(r_inner, 90)  # Top-left corner
wp_inner.LineTo(-inner_half_x, -inner_half_y + r_inner)
wp_inner.Arc(r_inner, 90)  # Bottom-left corner
wp_inner.Close()
inner_face = wp_inner.Face()

# Racetrack cross-section = outer - inner
racetrack_face = outer_face - inner_face

print(f"  [OK] Racetrack 2D profile created (rounded corners)")

# ============================================================
# Step 2: Extrude to 3D
# ============================================================
print("\nStep 2: Extrude to 3D")

# Extrude in Z direction
coil_local = racetrack_face.Extrude(h)

print(f"  Extruded height: {h*1000:.1f} mm")

# ============================================================
# Step 3: Move to coil center position
# ============================================================
print("\nStep 3: Position coil")

# The coil is currently at Z = [0, h] after extrusion
# We need to center it at Z = 0, so move by -h/2 in Z
# Final position: center at (cx, cy, cz) with Z centered at cz
coil_3d = coil_local.Move((cx, cy, cz - h/2))

print(f"  Coil center: ({cx*1000:.2f}, {cy*1000:.2f}, {cz*1000:.2f}) mm")
print(f"  Z offset: {-h/2*1000:.2f} mm (to center at Z=0)")

# ============================================================
# Step 4: Compute bounding box
# ============================================================
print("\nStep 4: Coil bounding box (mm)")

# Get bounding box from OCC shape
bbox = coil_3d.bounding_box
x_min, y_min, z_min = bbox[0]
x_max, y_max, z_max = bbox[1]

print(f"  X: [{x_min*1000:.1f}, {x_max*1000:.1f}] (size: {(x_max-x_min)*1000:.1f})")
print(f"  Y: [{y_min*1000:.1f}, {y_max*1000:.1f}] (size: {(y_max-y_min)*1000:.1f})")
print(f"  Z: [{z_min*1000:.1f}, {z_max*1000:.1f}] (size: {(z_max-z_min)*1000:.1f})")

# ============================================================
# Step 5: Export to STEP
# ============================================================
print("\nStep 5: Export to STEP")

step_file = os.path.join(work_dir, 'coil.step')

# Create OCCGeometry
geo = OCCGeometry(coil_3d)

# Export to STEP using OCC shape's method
try:
    # Try WriteStep if available
    from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_AsIs
    from OCC.Core.IFSelect import IFSelect_RetDone

    writer = STEPControl_Writer()
    writer.Transfer(coil_3d.IsShape(), STEPControl_AsIs)
    status = writer.Write(step_file)
    if status == IFSelect_RetDone:
        print(f"  Saved: {step_file}")
    else:
        print(f"  [WARNING] STEP export failed")
except ImportError:
    # Use netgen's built-in method if available
    try:
        coil_3d.WriteStep(step_file)
        print(f"  Saved: {step_file}")
    except AttributeError:
        print("  [INFO] STEP export not available")
        print("  Use Draw(coil_3d) in interactive mode for visualization")

# ============================================================
# Step 6: Visualize in Netgen GUI
# ============================================================
print("\nStep 6: Visualize in Netgen GUI")

# Import Draw from ngsolve (NGSolve's Draw, not webgui)
from ngsolve import Draw
import netgen.gui  # Activates GUI window

# Draw the OCC geometry
Draw(geo)
print("  GUI opened - close window to exit")

# ============================================================
# Summary
# ============================================================
print()
print("=" * 60)
print("COIL GEOMETRY COMPLETE")
print("=" * 60)
print()
print("Racetrack coil specifications:")
print(f"  Inner opening: {lx*1000:.1f} x {ly*1000:.1f} mm")
print(f"  Outer size: {(lx + 2*conductor_width)*1000:.1f} x {(ly + 2*conductor_width)*1000:.1f} mm")
print(f"  Conductor width: {conductor_width*1000:.1f} mm")
print(f"  Corner radii: inner={r_inner*1000:.1f}, outer={r_outer*1000:.1f} mm")
print(f"  Height: {h*1000:.1f} mm")
print(f"  Center: ({cx*1000:.2f}, {cy*1000:.2f}, {cz*1000:.2f}) mm")
print()
print("Output files:")
print(f"  - coil.step : STEP file (view in ParaView/CAD)")
print()
print("Next step:")
print("  python generate_hex_mesh.py")
print("=" * 60)

# For interactive visualization
coil = coil_3d  # Alias for convenience
geo = OCCGeometry(coil_3d)
