#!/usr/bin/env python
"""
GMSH Workflow for Radia-NGSolve

Demonstrates:
- GMSH mesh generation (Python API)
- Direct import to NGSolve
- Radia field computation
- Visualization options

Advantages of GMSH:
- [OK] Open source (free)
- [OK] Cross-platform (Windows/Linux/Mac)
- [OK] Python API + GUI
- [OK] Native NGSolve support
- [OK] Surface elements automatically included
- [OK] Can import CAD (STEP, IGES, etc.)

Requirements:
- gmsh (pip install gmsh)
- ngsolve
- radia
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import gmsh
from ngsolve import Mesh
import radia as rad


def create_box_mesh_with_gmsh():
    """Create a simple box mesh with GMSH."""
    print("="*60)
    print("GMSH Mesh Generation")
    print("="*60)

    # Initialize GMSH
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)

    # Create new model
    gmsh.model.add("box")

    # Create box geometry
    print("\n[1] Creating geometry...")
    lc = 0.01  # Characteristic length (mesh size)

    # Define points
    p1 = gmsh.model.geo.addPoint(-0.05, -0.05, -0.05, lc)
    p2 = gmsh.model.geo.addPoint(0.05, -0.05, -0.05, lc)
    p3 = gmsh.model.geo.addPoint(0.05, 0.05, -0.05, lc)
    p4 = gmsh.model.geo.addPoint(-0.05, 0.05, -0.05, lc)
    p5 = gmsh.model.geo.addPoint(-0.05, -0.05, 0.05, lc)
    p6 = gmsh.model.geo.addPoint(0.05, -0.05, 0.05, lc)
    p7 = gmsh.model.geo.addPoint(0.05, 0.05, 0.05, lc)
    p8 = gmsh.model.geo.addPoint(-0.05, 0.05, 0.05, lc)

    # Define lines (bottom face)
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    # Define lines (top face)
    l5 = gmsh.model.geo.addLine(p5, p6)
    l6 = gmsh.model.geo.addLine(p6, p7)
    l7 = gmsh.model.geo.addLine(p7, p8)
    l8 = gmsh.model.geo.addLine(p8, p5)

    # Define lines (vertical)
    l9 = gmsh.model.geo.addLine(p1, p5)
    l10 = gmsh.model.geo.addLine(p2, p6)
    l11 = gmsh.model.geo.addLine(p3, p7)
    l12 = gmsh.model.geo.addLine(p4, p8)

    # Define faces
    # Bottom face
    cl1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    s1 = gmsh.model.geo.addPlaneSurface([cl1])

    # Top face
    cl2 = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])
    s2 = gmsh.model.geo.addPlaneSurface([cl2])

    # Side faces
    cl3 = gmsh.model.geo.addCurveLoop([l1, l10, -l5, -l9])
    s3 = gmsh.model.geo.addPlaneSurface([cl3])

    cl4 = gmsh.model.geo.addCurveLoop([l2, l11, -l6, -l10])
    s4 = gmsh.model.geo.addPlaneSurface([cl4])

    cl5 = gmsh.model.geo.addCurveLoop([l3, l12, -l7, -l11])
    s5 = gmsh.model.geo.addPlaneSurface([cl5])

    cl6 = gmsh.model.geo.addCurveLoop([l4, l9, -l8, -l12])
    s6 = gmsh.model.geo.addPlaneSurface([cl6])

    # Define volume
    sl = gmsh.model.geo.addSurfaceLoop([s1, s2, s3, s4, s5, s6])
    v = gmsh.model.geo.addVolume([sl])

    # Synchronize
    gmsh.model.geo.synchronize()

    # Physical groups (for NGSolve materials)
    gmsh.model.addPhysicalGroup(3, [v], 1)  # Volume
    gmsh.model.setPhysicalName(3, 1, "magnet")

    gmsh.model.addPhysicalGroup(2, [s1, s2, s3, s4, s5, s6], 1)  # Surfaces
    gmsh.model.setPhysicalName(2, 1, "boundary")

    print("    Geometry created")

    # Generate mesh
    print("\n[2] Generating mesh...")
    gmsh.model.mesh.generate(3)  # 3D mesh

    # Get mesh statistics
    nodes = gmsh.model.mesh.getNodes()
    elements = gmsh.model.mesh.getElements()

    print(f"    Nodes: {len(nodes[0])}")
    print(f"    Elements: {len(elements[1][0])} (volume)")

    # Save to file
    msh_file = 'gmsh_box.msh'
    gmsh.write(msh_file)
    print(f"    Saved: {msh_file}")

    # Finalize
    gmsh.finalize()

    return msh_file


def use_with_ngsolve_and_radia(msh_file):
    """Use GMSH mesh with NGSolve and Radia."""
    print("\n" + "="*60)
    print("NGSolve + Radia Integration")
    print("="*60)

    # Import to NGSolve
    print("\n[3] Importing to NGSolve...")
    mesh = Mesh(msh_file)
    print(f"    NGSolve mesh loaded")
    print(f"    Vertices: {mesh.nv}")
    print(f"    Elements: {mesh.ne}")

    # Check surface elements
    print(f"    Surface elements: {mesh.ngmesh.nse}")
    if mesh.ngmesh.nse > 0:
        print("    [OK] Surface elements present - Netgen GUI compatible")

    # Convert to Radia (example)
    print("\n[4] Radia field computation...")

    # Create Radia magnet
    magnet = rad.ObjRecMag([0, 0, 0], [0.1, 0.1, 0.1], [0, 0, 954930])
    print("    Radia magnet created")

    # Compute field at mesh vertices (example)
    import numpy as np

    # Get a few vertex positions
    sample_points = []
    for i in range(min(5, mesh.nv)):
        v = mesh.vertices[i]
        sample_points.append([v.point[0], v.point[1], v.point[2]])

    print("\n    Field at sample points:")
    for i, pt in enumerate(sample_points):
        B = rad.Fld(magnet, 'b', pt)
        B_mag = np.linalg.norm(B)
        print(f"      Point {i+1}: |B| = {B_mag:.6f} T")

    print("\n[OK] GMSH workflow complete")


def main():
    print("="*60)
    print("GMSH Workflow Demonstration")
    print("="*60)
    print("\nAdvantages of GMSH:")
    print("  - Open source (free)")
    print("  - Cross-platform")
    print("  - Python API + GUI")
    print("  - Native NGSolve support")
    print("  - Surface elements included")

    # Generate mesh with GMSH
    msh_file = create_box_mesh_with_gmsh()

    # Use with NGSolve and Radia
    use_with_ngsolve_and_radia(msh_file)

    print("\n" + "="*60)
    print("Next steps:")
    print("="*60)
    print("\n1. View in GMSH GUI:")
    print(f"   gmsh {msh_file}")
    print("\n2. Convert to .vol (optional):")
    print("   from ngsolve import Mesh")
    print(f"   mesh = Mesh('{msh_file}')")
    print("   mesh.ngmesh.Save('gmsh_box.vol')")
    print("\n3. View .vol in Netgen GUI:")
    print("   python scripts/netgen_vol_viewer.py gmsh_box.vol")


if __name__ == '__main__':
    main()
