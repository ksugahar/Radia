#!/usr/bin/env python
"""
ParaView Publication-Quality Figure Generation

Demonstrates:
- Automated high-resolution figure generation with pvpython
- Vector field visualization (glyphs)
- Contour/isosurface rendering
- Camera positioning
- Publication-ready PNG/SVG export

Use case: Journal paper figures, presentation slides

Requirements:
- ParaView installation with pvpython
- Run as: pvpython generate_paraview_figure.py
  OR: python generate_paraview_figure.py (generates VTS, manual ParaView)
"""

import sys
import os

# Check if running under pvpython
try:
    from paraview.simple import *
    HAS_PARAVIEW = True
except ImportError:
    HAS_PARAVIEW = False
    print("Note: ParaView Python not available")
    print("This script will generate VTS file only")
    print("Open the VTS file in ParaView GUI manually")

# Add Radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
import radia as rad


def generate_field_vts():
    """Generate field VTS file with Radia."""
    print("="*60)
    print("Generating Radia Field VTS")
    print("="*60)


    # Create magnet
    print("\n[1] Creating magnet...")
    magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])
    print("    Magnet: 40x40x20 mm, Br=1.2T (NdFeB)")

    # Export field
    print("\n[2] Exporting field to VTS...")
    vts_file = 'field_for_paraview.vts'
    rad.FldVTS(magnet, vts_file,
               [-0.1, 0.1], [-0.1, 0.1], [0.02, 0.15],
               41, 41, 27, 1, 0, 1.0)
    print(f"    VTS file created: {vts_file}")

    return vts_file


def generate_paraview_figure(vts_file):
    """Generate publication figure with ParaView Python."""
    if not HAS_PARAVIEW:
        print("\n[SKIP] ParaView Python not available")
        print(f"\nTo generate figure manually:")
        print(f"1. Open ParaView")
        print(f"2. File > Open > {vts_file}")
        print(f"3. Apply filters (Glyph, Contour, etc.)")
        print(f"4. File > Save Screenshot (3000x2000, PNG or SVG)")
        return

    print("\n" + "="*60)
    print("Generating ParaView Figure")
    print("="*60)

    # Load VTS
    print("\n[1] Loading VTS file...")
    reader = XMLStructuredGridReader(FileName=[vts_file])
    reader.UpdatePipeline()

    # Figure 1: Scalar field with contours
    print("\n[2] Figure 1: B magnitude with contours")
    view = CreateView('RenderView')
    view.ViewSize = [1200, 900]
    view.Background = [1, 1, 1]  # White background

    # Display mesh
    display = Show(reader, view)
    ColorBy(display, ('POINTS', 'B_magnitude'))

    # Add contour lines
    contour = Contour(Input=reader)
    contour.ContourBy = ['POINTS', 'B_magnitude']
    contour.Isosurfaces = [0.1, 0.2, 0.3, 0.4, 0.5]

    contour_display = Show(contour, view)
    contour_display.LineWidth = 2.0

    # Camera
    view.CameraPosition = [0.3, 0.3, 0.4]
    view.CameraFocalPoint = [0, 0, 0.08]
    view.CameraViewUp = [0, 0, 1]

    # Save
    print("    Saving figure1_contours.png...")
    SaveScreenshot('figure1_contours.png', view, ImageResolution=[3000, 2000])

    Delete(view)

    # Figure 2: Vector field with arrows
    print("\n[3] Figure 2: B field vectors")
    view = CreateView('RenderView')
    view.ViewSize = [1200, 900]
    view.Background = [1, 1, 1]

    # Glyph (arrows)
    glyph = Glyph(Input=reader, GlyphType='Arrow')
    glyph.OrientationArray = ['POINTS', 'B_field']
    glyph.ScaleArray = ['POINTS', 'B_magnitude']
    glyph.ScaleFactor = 0.02
    glyph.GlyphMode = 'Every Nth Point'
    glyph.Stride = [2, 2, 2]

    glyph_display = Show(glyph, view)
    ColorBy(glyph_display, ('POINTS', 'B_magnitude'))

    # Camera
    view.CameraPosition = [0.3, 0.2, 0.5]
    view.CameraFocalPoint = [0, 0, 0.08]
    view.CameraViewUp = [0, 0, 1]

    # Save
    print("    Saving figure2_vectors.png...")
    SaveScreenshot('figure2_vectors.png', view, ImageResolution=[3000, 2000])

    Delete(view)

    # Figure 3: Slice with streamlines
    print("\n[4] Figure 3: Z-plane slice with streamlines")
    view = CreateView('RenderView')
    view.ViewSize = [1200, 900]
    view.Background = [1, 1, 1]

    # Slice
    slice_filter = Slice(Input=reader)
    slice_filter.SliceType = 'Plane'
    slice_filter.SliceType.Origin = [0, 0, 0.05]
    slice_filter.SliceType.Normal = [0, 0, 1]

    slice_display = Show(slice_filter, view)
    ColorBy(slice_display, ('POINTS', 'B_magnitude'))

    # Stream tracer
    stream = StreamTracer(Input=reader, SeedType='Line')
    stream.Vectors = ['POINTS', 'B_field']
    stream.MaximumStreamlineLength = 0.2
    stream.SeedType.Point1 = [-0.08, -0.08, 0.05]
    stream.SeedType.Point2 = [0.08, 0.08, 0.05]
    stream.SeedType.Resolution = 20

    stream_display = Show(stream, view)
    stream_display.LineWidth = 2.0

    # Camera (top view)
    view.CameraPosition = [0, 0, 0.5]
    view.CameraFocalPoint = [0, 0, 0.05]
    view.CameraViewUp = [0, 1, 0]

    # Save
    print("    Saving figure3_streamlines.png...")
    SaveScreenshot('figure3_streamlines.png', view, ImageResolution=[3000, 2000])

    Delete(view)

    print("\n[OK] ParaView figures generated")
    print("\nGenerated files:")
    print("  - figure1_contours.png (3000x2000)")
    print("  - figure2_vectors.png (3000x2000)")
    print("  - figure3_streamlines.png (3000x2000)")


def main():
    # Generate VTS file
    vts_file = generate_field_vts()

    # Generate ParaView figures (if pvpython available)
    generate_paraview_figure(vts_file)

    print("\n" + "="*60)
    print("Figure Generation Complete")
    print("="*60)


if __name__ == '__main__':
    main()
