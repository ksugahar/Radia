#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ParaView Visualization Script for Electromagnet Simulation

This script loads the Radia geometry and magnetic field distribution
and creates a visualization with glyphs showing the magnetic field vectors.

Usage:
    pvpython visualize_paraview.py

Requirements:
    - ParaView installation with pvpython
    - Radia_model.vtk (geometry from Radia)
    - field_distribution.vtk (magnetic field data)
"""

import os
import sys

# Import ParaView modules
try:
    from paraview.simple import *
except ImportError:
    print("[ERROR] Could not import paraview.simple")
    print("This script must be run with pvpython (ParaView Python)")
    print("Usage: pvpython visualize_paraview.py")
    sys.exit(1)

print("=" * 70)
print("PARAVIEW VISUALIZATION - ELECTROMAGNET SIMULATION")
print("=" * 70)

# Get script directory
script_dir = os.path.dirname(os.path.abspath(__file__))

# File paths
radia_model_file = os.path.join(script_dir, 'Radia_model.vtk')
field_file = os.path.join(script_dir, 'field_distribution.vtk')

# Check files exist
if not os.path.exists(radia_model_file):
    print(f"\n[ERROR] File not found: {radia_model_file}")
    print("Please run main_simulation_workflow.py first to generate VTK files")
    sys.exit(1)

if not os.path.exists(field_file):
    print(f"\n[ERROR] File not found: {field_file}")
    print("Please run main_simulation_workflow.py first to generate VTK files")
    sys.exit(1)

print(f"\n[Step 1/5] Loading Radia geometry model...")
print(f"  File: {os.path.basename(radia_model_file)}")

# Load Radia geometry (coil + yoke)
radia_model = LegacyVTKReader(FileNames=[radia_model_file])
radia_display = Show(radia_model)
radia_display.Representation = 'Surface'
radia_display.ColorArrayName = ['CELL_DATA', 'colors']

print(f"  [OK] Loaded Radia geometry")

# ========================================================================
# Load Field Distribution
# ========================================================================
print(f"\n[Step 2/5] Loading magnetic field distribution...")
print(f"  File: {os.path.basename(field_file)}")

field_data = LegacyVTKReader(FileNames=[field_file])
field_display = Show(field_data)
field_display.Representation = 'Points'
field_display.PointSize = 2.0

print(f"  [OK] Loaded field data")

# ========================================================================
# Create Glyph Filter for Field Vectors
# ========================================================================
print(f"\n[Step 3/5] Creating glyph visualization...")

# Create glyph filter (arrows showing magnetic field vectors)
glyph = Glyph(Input=field_data, GlyphType='Arrow')
glyph.OrientationArray = ['POINTS', 'B_field']
glyph.ScaleArray = ['POINTS', 'B_magnitude']
glyph.ScaleFactor = 50.0  # Adjust arrow size
glyph.GlyphMode = 'All Points'

# Show glyphs
glyph_display = Show(glyph)
glyph_display.Representation = 'Surface'

# Color by magnetic field magnitude
ColorBy(glyph_display, ('POINTS', 'B_magnitude'))

# Get color transfer function
lut = GetColorTransferFunction('B_magnitude')
lut.RescaleTransferFunction(0.0, 1.0)  # Tesla range

# Show color bar
glyph_display.SetScalarBarVisibility(GetActiveView(), True)

print(f"  [OK] Glyph filter created")
print(f"  Arrow scale factor: {glyph.ScaleFactor}")

# ========================================================================
# Setup View and Camera
# ========================================================================
print(f"\n[Step 4/5] Setting up view...")

# Get render view
renderView = GetActiveViewOrCreate('RenderView')

# Set background color (white)
renderView.Background = [1.0, 1.0, 1.0]

# Reset camera to show all objects
renderView.ResetCamera()

# Adjust camera for better viewing angle
camera = renderView.GetActiveCamera()
camera.Elevation(30)
camera.Azimuth(45)
renderView.ResetCamera()

print(f"  [OK] View configured")

# ========================================================================
# Render and Save Screenshot
# ========================================================================
print(f"\n[Step 5/5] Rendering visualization...")

# Render
Render()

# Save screenshot
screenshot_file = os.path.join(script_dir, 'electromagnet.png')
SaveScreenshot(screenshot_file, renderView, ImageResolution=[1920, 1080])

print(f"  [OK] Rendered")
print(f"  Screenshot saved: {os.path.basename(screenshot_file)}")

# ========================================================================
# Summary
# ========================================================================
print("\n" + "=" * 70)
print("VISUALIZATION COMPLETE")
print("=" * 70)
print(f"\nVisualization includes:")
print(f"  1. Radia geometry (coil + yoke) with colors")
print(f"  2. Magnetic field points")
print(f"  3. Glyph arrows showing field vectors")
print(f"  4. Color mapping by field magnitude (Tesla)")
print(f"\nScreenshot: {os.path.basename(screenshot_file)}")
print(f"\nTo view interactively:")
print(f"  1. Open ParaView")
print(f"  2. Load Radia_model.vtk")
print(f"  3. Load field_distribution.vtk")
print(f"  4. Apply Glyph filter to field_distribution.vtk")
print(f"  5. Set Orientation Array = B_field")
print(f"  6. Set Scale Array = B_magnitude")
print(f"  7. Color by B_magnitude")
print("=" * 70)
