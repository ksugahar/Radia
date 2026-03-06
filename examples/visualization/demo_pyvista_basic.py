#!/usr/bin/env python
"""
PyVista Basic Visualization Demo for Radia-NGSolve

Demonstrates:
- Field visualization with PyVista
- Vector field arrows
- Scalar field colormaps
- Interactive exploration (インタラクティブウィンドウ)

Use case: Daily development, quick field checks, parameter studies

Note: PyVistaは通常のPythonスクリプト（.py）から使用可能です。
      Jupyterノートブック不要で、インタラクティブウィンドウが開きます。
      マウスでカメラ操作、ズーム、回転が可能です。
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import radia as rad
import pyvista as pv
import numpy as np

def main():
    print("="*60)
    print("PyVista Basic Visualization Demo")
    print("="*60)


    # Create rectangular permanent magnet
    print("\n[1] Creating magnet...")
    magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])
    print("    Magnet: 40x40x20 mm, Br=1.2T (NdFeB)")

    # Export field to VTS
    print("\n[2] Exporting field to VTS...")
    vts_file = 'field_output.vts'
    rad.FldVTS(magnet, vts_file,
               [-0.1, 0.1], [-0.1, 0.1], [0.02, 0.15],
               41, 41, 27, 1, 0, 1.0)
    print(f"    VTS file created: {vts_file}")

    # Load with PyVista
    print("\n[3] Loading VTS with PyVista...")
    grid = pv.read(vts_file)
    print(f"    Grid: {grid.n_points} points, {grid.n_cells} cells")

    # Visualization 1: Scalar field (magnitude)
    print("\n[4] Visualization 1: B magnitude (scalar)")
    plotter = pv.Plotter(window_size=[800, 600])
    plotter.add_mesh(grid, scalars='B_magnitude',
                     cmap='coolwarm',
                     opacity=0.8,
                     show_edges=False)
    plotter.add_scalar_bar('B magnitude [T]')
    plotter.show_grid()
    plotter.add_axes()
    plotter.camera_position = 'iso'
    plotter.show()

    # Visualization 2: Vector field (arrows)
    print("\n[5] Visualization 2: B field vectors")
    plotter = pv.Plotter(window_size=[800, 600])

    # Add surface mesh with transparency
    plotter.add_mesh(grid, scalars='B_magnitude',
                     cmap='viridis',
                     opacity=0.3)

    # Add vector arrows
    arrows = grid.glyph(orient='B_field', scale='B_magnitude', factor=0.02)
    plotter.add_mesh(arrows, color='black')

    plotter.add_scalar_bar('B magnitude [T]')
    plotter.add_axes()
    plotter.camera_position = 'iso'
    plotter.show()

    # Visualization 3: Slice
    print("\n[6] Visualization 3: Z-plane slice")
    plotter = pv.Plotter(window_size=[800, 600])

    # Create slice at z=0.05m
    slice_z = grid.slice(normal='z', origin=[0, 0, 0.05])

    plotter.add_mesh(slice_z, scalars='B_magnitude',
                     cmap='rainbow',
                     show_edges=True)

    # Add arrows on slice
    arrows = slice_z.glyph(orient='B_field', scale='B_magnitude', factor=0.02)
    plotter.add_mesh(arrows, color='white')

    plotter.add_scalar_bar('B magnitude [T]')
    plotter.add_axes()
    plotter.camera_position = 'xy'
    plotter.show()

    # Cleanup
    os.remove(vts_file)
    print("\n[OK] Demo complete")


if __name__ == '__main__':
    main()
