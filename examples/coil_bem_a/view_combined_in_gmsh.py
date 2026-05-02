"""Open the combined BEM-A mesh + PEEC filament .msh in an INTERACTIVE
GMSH window.  Stays open until the user closes the window (gmsh.fltk.run
blocks).  Rotate / zoom / pan with the mouse; toggle physical groups via
Tools -> Visibility (Ctrl+Shift+V).

Run from a console with X / RDP available::

    python view_combined_in_gmsh.py
"""
import os
import sys

import gmsh

HERE = os.path.dirname(os.path.abspath(__file__))
MSH = os.path.join(HERE, "rect_torus_loft_combined.msh")


def main():
    gmsh.initialize()
    gmsh.open(MSH)

    # Display options to make the filaments visible against the surface
    gmsh.option.setNumber("Mesh.SurfaceFaces", 1)
    gmsh.option.setNumber("Mesh.SurfaceEdges", 1)
    gmsh.option.setNumber("Mesh.Lines", 1)         # show 1D elements
    gmsh.option.setNumber("Mesh.LineWidth", 1.5)   # mesh edge thickness
    gmsh.option.setNumber("Geometry.LineWidth", 2) # filament thickness
    gmsh.option.setNumber("Mesh.ColorCarousel", 2)  # color by physical group
    gmsh.option.setNumber("General.RotationX", 35)
    gmsh.option.setNumber("General.RotationY", 0)
    gmsh.option.setNumber("General.RotationZ", 30)

    print(f"Opened: {MSH}")
    print()
    print("Mouse: rotate / zoom / pan")
    print("Tools -> Visibility (Ctrl+Shift+V): toggle physical groups")
    print("  - body  : 298 surface tris (BEM-A lateral mesh)")
    print("  - source: 2 cap tris")
    print("  - sink  : 2 cap tris")
    print("  - peec_fil : 304 line segments (PEEC perimeter filaments,")
    print("               16 filaments x 19 segments each)")
    print()
    print("Close the window to exit.")
    gmsh.fltk.run()
    gmsh.finalize()


if __name__ == "__main__":
    main()
