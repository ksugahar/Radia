#!/usr/bin/env python
"""
Verify Combined Geometry: C-Type Yoke + Racetrack Coil

Loads both yoke.step and coil.step to visualize the complete electromagnet model.
- Yoke: Cyan (water blue)
- Coil: Red

All units are in meters (m).

Usage:
  python verify_geometry.py           # Visualize combined model
  py -i verify_geometry.py            # Interactive mode

Requirements:
  - yoke.step from generate_hex_mesh.py
  - coil.step from verify_coil.py
"""

import sys
import os

work_dir = os.path.dirname(os.path.abspath(__file__))

# Use custom NGSolve with GUI support and patches (PR#231, PR#232)
NGSOLVE_PATH = r'S:\NGSolve\01_GitHub\install_ksugahar\lib\site-packages'
sys.path.insert(0, NGSOLVE_PATH)

from netgen.occ import *

print("=" * 60)
print("VERIFY COMBINED GEOMETRY: C-TYPE YOKE + RACETRACK COIL")
print("=" * 60)
print("All units in meters (m)")
print()

# ============================================================
# Step 1: Load yoke geometry
# ============================================================
print("Step 1: Load yoke geometry")

yoke_step = os.path.join(work_dir, 'yoke_united.step')
if not os.path.exists(yoke_step):
    print(f"  [ERROR] yoke_united.step not found!")
    print(f"  Run: python generate_hex_mesh.py first")
    sys.exit(1)

yoke_raw = OCCGeometry(yoke_step).shape
print(f"  Loaded: {yoke_step}")

# yoke.step is in mm (from Cubit), scale to meters
# Scale(origin, scale_factor)
yoke = yoke_raw.Scale(Pnt(0, 0, 0), 0.001)
print(f"  Scaled from mm to m")

# Yoke bounding box (in meters)
bbox_yoke = yoke.bounding_box
print(f"  Yoke bounding box (m):")
print(f"    X: [{bbox_yoke[0][0]:.4f}, {bbox_yoke[1][0]:.4f}]")
print(f"    Y: [{bbox_yoke[0][1]:.4f}, {bbox_yoke[1][1]:.4f}]")
print(f"    Z: [{bbox_yoke[0][2]:.4f}, {bbox_yoke[1][2]:.4f}]")

# Set yoke color to cyan (water blue)
# RGB: (0, 1, 1) = cyan
for face in yoke.faces:
    face.col = (0.3, 0.7, 1.0)  # Light blue / cyan

print("  Color: Cyan (light blue)")

# ============================================================
# Step 2: Load coil geometry
# ============================================================
print("\nStep 2: Load coil geometry")

coil_step = os.path.join(work_dir, 'coil.step')
if not os.path.exists(coil_step):
    print(f"  [ERROR] coil.step not found!")
    print(f"  Run: python verify_coil.py first")
    sys.exit(1)

coil = OCCGeometry(coil_step).shape
print(f"  Loaded: {coil_step}")

# Coil bounding box (already in meters)
bbox_coil = coil.bounding_box
print(f"  Coil bounding box (m):")
print(f"    X: [{bbox_coil[0][0]:.4f}, {bbox_coil[1][0]:.4f}]")
print(f"    Y: [{bbox_coil[0][1]:.4f}, {bbox_coil[1][1]:.4f}]")
print(f"    Z: [{bbox_coil[0][2]:.4f}, {bbox_coil[1][2]:.4f}]")

# Set coil color to red
# RGB: (1, 0, 0) = red
for face in coil.faces:
    face.col = (1.0, 0.2, 0.2)  # Red

print("  Color: Red")

# ============================================================
# Step 3: Combine geometries
# ============================================================
print("\nStep 3: Combine geometries")

# Create a Compound of both shapes (not boolean union)
# This keeps them as separate bodies for visualization
combined = Glue([yoke, coil])

print(f"  Combined yoke and coil into single geometry")

# Combined bounding box (in meters)
bbox_combined = combined.bounding_box
print(f"  Combined bounding box (m):")
print(f"    X: [{bbox_combined[0][0]:.4f}, {bbox_combined[1][0]:.4f}]")
print(f"    Y: [{bbox_combined[0][1]:.4f}, {bbox_combined[1][1]:.4f}]")
print(f"    Z: [{bbox_combined[0][2]:.4f}, {bbox_combined[1][2]:.4f}]")

# ============================================================
# Step 4: Export combined STEP
# ============================================================
print("\nStep 4: Export combined STEP")

combined_step = os.path.join(work_dir, 'electromagnet.step')
try:
    combined.WriteStep(combined_step)
    print(f"  Saved: {combined_step}")
except AttributeError:
    print("  [INFO] STEP export not available")

# ============================================================
# Step 5: Visualize in Netgen GUI
# ============================================================
print("\nStep 5: Visualize in Netgen GUI")

# Create OCCGeometry for visualization
geo = OCCGeometry(combined)

# Import Draw from ngsolve
from ngsolve import Draw
import netgen.gui  # Activates GUI window

# Draw the combined geometry
Draw(geo)
print("  GUI opened - close window to exit")
print("  Yoke: Cyan (light blue)")
print("  Coil: Red")

# ============================================================
# Summary
# ============================================================
print()
print("=" * 60)
print("COMBINED GEOMETRY COMPLETE")
print("=" * 60)
print()
print("Components:")
print("  - Yoke: C-type electromagnet yoke (Cyan)")
print("  - Coil: Racetrack coil (Red)")
print()
print("Output files:")
print(f"  - electromagnet.step : Combined STEP file (m)")
print()
print("Next step:")
print("  python generate_hex_mesh.py  (for Radia hex mesh)")
print("  python generate_tet_mesh.py  (for NGSolve tet mesh)")
print("=" * 60)

# For interactive visualization
electromagnet = combined
