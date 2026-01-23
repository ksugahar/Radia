#!/usr/bin/env python
"""
Verify Yoke Geometry Only

Loads yoke.step to visualize the C-type electromagnet yoke.
All units are in meters (m).

Usage:
  python verify_yoke.py           # Visualize yoke
  py -i verify_yoke.py            # Interactive mode
"""

import sys
import os

work_dir = os.path.dirname(os.path.abspath(__file__))

# Use custom NGSolve with GUI support and patches (PR#231, PR#232)
NGSOLVE_PATH = r'S:\NGSolve\01_GitHub\install_ksugahar\lib\site-packages'
sys.path.insert(0, NGSOLVE_PATH)

from netgen.occ import *

print("=" * 60)
print("VERIFY YOKE GEOMETRY: C-TYPE ELECTROMAGNET")
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
yoke = yoke_raw.Scale(Pnt(0, 0, 0), 0.001)
print(f"  Scaled from mm to m")

# ============================================================
# Step 2: Bounding box
# ============================================================
print("\nStep 2: Yoke bounding box (m)")

bbox = yoke.bounding_box
x_min, y_min, z_min = bbox[0]
x_max, y_max, z_max = bbox[1]

print(f"  X: [{x_min:.4f}, {x_max:.4f}] (size: {(x_max-x_min):.4f})")
print(f"  Y: [{y_min:.4f}, {y_max:.4f}] (size: {(y_max-y_min):.4f})")
print(f"  Z: [{z_min:.4f}, {z_max:.4f}] (size: {(z_max-z_min):.4f})")

# ============================================================
# Step 3: Set color to cyan
# ============================================================
print("\nStep 3: Set color")

# Set color on original (mm) shape for GUI display
for face in yoke_raw.faces:
    face.col = (0.3, 0.7, 1.0)  # Cyan (light blue)

# Also set on scaled shape
for face in yoke.faces:
    face.col = (0.3, 0.7, 1.0)  # Cyan (light blue)

print("  Color: Cyan (light blue)")

# ============================================================
# Step 4: Visualize in Netgen GUI
# ============================================================
print("\nStep 4: Visualize in Netgen GUI")

# Use original mm scale for better visibility in GUI
# (meters scale is too small for GUI visualization)
geo_mm = OCCGeometry(yoke_raw)  # Original mm scale

from ngsolve import Draw
import netgen.gui

Draw(geo_mm)
print("  GUI opened - close window to exit")
print("  NOTE: Displayed in mm scale for visibility")

# ============================================================
# Summary
# ============================================================
print()
print("=" * 60)
print("YOKE GEOMETRY LOADED")
print("=" * 60)
print()
print("Yoke size (m):")
print(f"  X: {(x_max-x_min):.4f}")
print(f"  Y: {(y_max-y_min):.4f}")
print(f"  Z: {(z_max-z_min):.4f}")
print()
print("Center (m):")
print(f"  X: {(x_min+x_max)/2:.4f}")
print(f"  Y: {(y_min+y_max)/2:.4f}")
print(f"  Z: {(z_min+z_max)/2:.4f}")
print("=" * 60)

# For interactive mode
yoke_shape = yoke
