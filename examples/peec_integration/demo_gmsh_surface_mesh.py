#!/usr/bin/env python
"""
GMSH Surface Mesh for PEEC Conductors

CRITICAL: PEEC requires SURFACE MESH ONLY, not volume mesh.

Workflow:
    GMSH surface mesh → .msh → NGSolve → Radia PEEC

Key Points:
    - Conductors: Surface elements (Tri3, Quad4) only
    - Skin effect: Handled by SIBC (Surface Impedance Boundary Condition)
    - No volume discretization needed for conductors
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import gmsh
import numpy as np


def create_coil_surface_mesh(radius, width, height, n_radial=8, n_axial=4):
    """
    Create surface mesh for rectangular cross-section coil.

    Args:
        radius: Coil mean radius (m)
        width: Cross-section width (radial direction) (m)
        height: Cross-section height (axial direction) (m)
        n_radial: Number of divisions in radial direction
        n_axial: Number of divisions in axial direction

    Returns:
        Path to generated .msh file
    """
    print("=" * 60)
    print("GMSH Surface Mesh for PEEC Conductor")
    print("=" * 60)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("coil_surface")

    # Coil geometry parameters
    r_inner = radius - width / 2
    r_outer = radius + width / 2
    z_bottom = -height / 2
    z_top = height / 2

    print(f"\n[1] Creating coil geometry:")
    print(f"    Mean radius: {radius} m")
    print(f"    Width (radial): {width} m")
    print(f"    Height (axial): {height} m")
    print(f"    Inner radius: {r_inner} m")
    print(f"    Outer radius: {r_outer} m")

    # Create rectangular cross-section surface (in XZ plane)
    # Will be revolved around Z axis to create coil surface
    p1 = gmsh.model.geo.addPoint(r_inner, 0, z_bottom)
    p2 = gmsh.model.geo.addPoint(r_outer, 0, z_bottom)
    p3 = gmsh.model.geo.addPoint(r_outer, 0, z_top)
    p4 = gmsh.model.geo.addPoint(r_inner, 0, z_top)

    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surf = gmsh.model.geo.addPlaneSurface([loop])

    # Revolve to create coil surface
    print("\n[2] Revolving surface to create coil...")
    axis_point = [0, 0, 0]
    axis_direction = [0, 0, 1]  # Z axis
    angle = 2 * np.pi  # Full revolution

    gmsh.model.geo.revolve(
        [(2, surf)],  # Surface to revolve
        *axis_point,
        *axis_direction,
        angle
    )

    gmsh.model.geo.synchronize()

    # Set mesh size
    print(f"\n[3] Setting mesh parameters:")
    print(f"    Radial divisions: {n_radial}")
    print(f"    Axial divisions: {n_axial}")

    # Characteristic length based on divisions
    lc_radial = width / n_radial
    lc_axial = height / n_axial
    lc = min(lc_radial, lc_axial)

    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc)

    # Define physical group for PEEC (surface only!)
    print("\n[4] Defining physical surface group...")
    surfaces = gmsh.model.getEntities(2)
    if surfaces:
        surface_tags = [s[1] for s in surfaces]
        gmsh.model.addPhysicalGroup(2, surface_tags, 1)
        gmsh.model.setPhysicalName(2, 1, "conductor")
        print(f"    Surfaces: {len(surface_tags)}")

    # Generate SURFACE mesh only (dim=2)
    print("\n[5] Generating surface mesh...")
    gmsh.model.mesh.generate(2)  # 2D surface mesh, NOT 3D volume

    # Get statistics
    nodes = gmsh.model.mesh.getNodes()
    surf_elements = gmsh.model.mesh.getElements(2)  # Surface elements

    print(f"    Nodes: {len(nodes[0])}")
    if surf_elements[1]:
        total_surf = sum(len(e) for e in surf_elements[1])
        print(f"    Surface elements: {total_surf}")
        print("    ✅ Surface mesh only - correct for PEEC")

    # Check for volume elements (should be ZERO for PEEC)
    vol_elements = gmsh.model.mesh.getElements(3)
    if vol_elements[1] and any(len(e) > 0 for e in vol_elements[1]):
        print("    ⚠️  WARNING: Volume elements found - PEEC only needs surface!")
    else:
        print("    ✅ No volume elements - correct for PEEC")

    # Save mesh
    msh_file = 'coil_surface.msh'
    gmsh.write(msh_file)
    print(f"\n[6] Saved surface mesh: {msh_file}")
    print("    This .msh file contains SURFACE ELEMENTS ONLY")
    print("    Ready for PEEC conductor import")

    gmsh.finalize()

    return msh_file


def main():
    print("\nPEEC Conductor Surface Mesh Generation\n")
    print("Key Concept:")
    print("  - PEEC models conductors using SURFACE currents only")
    print("  - Skin effect handled by SIBC (Surface Impedance)")
    print("  - No volume discretization needed")
    print()

    # Create surface mesh for coil
    msh_file = create_coil_surface_mesh(
        radius=0.05,     # 50mm mean radius
        width=0.002,     # 2mm width
        height=0.002,    # 2mm height
        n_radial=4,      # 4 divisions radially
        n_axial=4        # 4 divisions axially
    )

    print("\n" + "=" * 60)
    print("Next Steps")
    print("=" * 60)
    print("\n1. View in GMSH:")
    print(f"   gmsh {msh_file}")

    print("\n2. Import to NGSolve:")
    print(f"   from ngsolve import Mesh")
    print(f"   mesh = Mesh('{msh_file}')")

    print("\n3. Convert to Radia PEEC (planned API):")
    print(f"   from peec_mesh_import import surface_mesh_to_peec")
    print(f"   conductor = surface_mesh_to_peec(mesh, sigma=5.8e7)")

    print("\n4. Alternative - Use CndLoop for simple geometry:")
    print("   coil = rad.CndLoop([0,0,0], 0.05, [0,0,1], 'r',")
    print("                       0.002, 0.002, 5.8e7, 8, 36)")


if __name__ == '__main__':
    main()
