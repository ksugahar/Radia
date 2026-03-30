#!/usr/bin/env python
"""
Circular Coil Model in GMSH

Creates a toroidal coil (circular coil with rectangular cross-section) using GMSH.

Workflow:
    1. Define rectangular cross-section in XZ plane
    2. Revolve around Z axis to create toroidal coil
    3. Generate surface mesh (for PEEC conductor)
    4. Visualize in GMSH GUI

Parameters:
    - Mean radius: Coil center radius (m)
    - Wire width: Cross-section width (radial direction) (m)
    - Wire height: Cross-section height (axial direction) (m)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import gmsh
import numpy as np


def create_circular_coil(mean_radius=0.05, wire_width=0.004, wire_height=0.004,
                         mesh_size=0.001, show_gui=True):
    """
    Create circular coil with rectangular cross-section.

    Args:
        mean_radius: Coil mean radius (m)
        wire_width: Wire cross-section width (radial) (m)
        wire_height: Wire cross-section height (axial) (m)
        mesh_size: Mesh element size (m)
        show_gui: Show GMSH GUI (True/False)

    Returns:
        Path to generated .msh file
    """
    print("=" * 60)
    print("Circular Coil Model in GMSH")
    print("=" * 60)

    # Initialize GMSH
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("circular_coil")

    # Coil geometry parameters
    r_inner = mean_radius - wire_width / 2
    r_outer = mean_radius + wire_width / 2
    z_bottom = -wire_height / 2
    z_top = wire_height / 2

    print(f"\n[1] Coil Parameters:")
    print(f"    Mean radius: {mean_radius * 1000:.1f} mm")
    print(f"    Wire width (radial): {wire_width * 1000:.2f} mm")
    print(f"    Wire height (axial): {wire_height * 1000:.2f} mm")
    print(f"    Inner radius: {r_inner * 1000:.1f} mm")
    print(f"    Outer radius: {r_outer * 1000:.1f} mm")

    # Create rectangular cross-section in XZ plane
    print("\n[2] Creating rectangular cross-section (XZ plane)...")
    p1 = gmsh.model.geo.addPoint(r_inner, 0, z_bottom)
    p2 = gmsh.model.geo.addPoint(r_outer, 0, z_bottom)
    p3 = gmsh.model.geo.addPoint(r_outer, 0, z_top)
    p4 = gmsh.model.geo.addPoint(r_inner, 0, z_top)

    l1 = gmsh.model.geo.addLine(p1, p2)  # Bottom edge
    l2 = gmsh.model.geo.addLine(p2, p3)  # Outer edge
    l3 = gmsh.model.geo.addLine(p3, p4)  # Top edge
    l4 = gmsh.model.geo.addLine(p4, p1)  # Inner edge

    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surf = gmsh.model.geo.addPlaneSurface([loop])

    # Revolve around Z axis to create toroidal coil
    print("[3] Revolving around Z axis to create toroidal coil...")
    axis_point = [0, 0, 0]      # Origin
    axis_direction = [0, 0, 1]   # Z axis
    angle = 2 * np.pi            # Full revolution (360 degrees)

    gmsh.model.geo.revolve(
        [(2, surf)],      # Surface entity to revolve
        *axis_point,
        *axis_direction,
        angle
    )

    gmsh.model.geo.synchronize()

    # Set mesh size
    print(f"\n[4] Setting mesh parameters:")
    print(f"    Mesh size: {mesh_size * 1000:.3f} mm")

    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size)

    # For better quality on curved surfaces
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 1)
    gmsh.option.setNumber("Mesh.MinimumCirclePoints", 20)

    # Define physical group for PEEC (surface only!)
    print("\n[5] Defining physical surface group (PEEC conductor)...")
    surfaces = gmsh.model.getEntities(2)
    if surfaces:
        surface_tags = [s[1] for s in surfaces]
        gmsh.model.addPhysicalGroup(2, surface_tags, 1)
        gmsh.model.setPhysicalName(2, 1, "conductor")
        print(f"    Surfaces: {len(surface_tags)}")

    # Generate SURFACE mesh (dim=2 for PEEC)
    print("\n[6] Generating surface mesh...")
    gmsh.model.mesh.generate(2)  # 2D surface mesh ONLY

    # Get statistics
    nodes = gmsh.model.mesh.getNodes()
    surf_elements = gmsh.model.mesh.getElements(2)

    print(f"\n[7] Mesh Statistics:")
    print(f"    Nodes: {len(nodes[0])}")
    if surf_elements[1]:
        total_surf = sum(len(e) for e in surf_elements[1])
        print(f"    Surface elements: {total_surf}")

        # Element type breakdown
        for i, elem_type in enumerate(surf_elements[0]):
            elem_name = gmsh.model.mesh.getElementProperties(elem_type)[0]
            num_elems = len(surf_elements[1][i])
            print(f"      {elem_name}: {num_elems}")

    # Check for volume elements (should be ZERO for PEEC)
    vol_elements = gmsh.model.mesh.getElements(3)
    if vol_elements[1] and any(len(e) > 0 for e in vol_elements[1]):
        print("    ⚠️  WARNING: Volume elements found - PEEC only needs surface!")
    else:
        print("    [OK] No volume elements - correct for PEEC")

    # Save mesh
    msh_file = 'circular_coil.msh'
    gmsh.write(msh_file)
    print(f"\n[8] Saved mesh: {msh_file}")
    print("    [OK] Surface mesh only - ready for PEEC")

    # Show in GMSH GUI
    if show_gui:
        print("\n[9] Opening GMSH GUI...")
        print("    Visualization:")
        print("      - Rotate: Left drag")
        print("      - Zoom: Mouse wheel")
        print("      - Pan: Middle drag")
        print("      - Toggle mesh: Press '0'")
        print("      - Show edges: Press 'e'")

        gmsh.fltk.run()

    gmsh.finalize()

    return msh_file


def create_coil_with_field_data(mean_radius=0.05, wire_width=0.004,
                                 wire_height=0.004, current_density=1.5e6):
    """
    Create circular coil with current density field data.

    Args:
        mean_radius: Coil mean radius (m)
        wire_width: Wire width (m)
        wire_height: Wire height (m)
        current_density: Surface current density (A/m^2)
    """
    print("=" * 60)
    print("Circular Coil with Current Density Field")
    print("=" * 60)

    # Create mesh
    msh_file = create_circular_coil(mean_radius, wire_width, wire_height,
                                     mesh_size=0.001, show_gui=False)

    # Re-open to add field data
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.open(msh_file)

    print("\n[10] Adding current density field data...")

    # Create view for current density
    view = gmsh.view.add("Current Density [A/m^2]")

    # Get surface elements
    surf_elements = gmsh.model.mesh.getElements(2)

    data_list = []
    for elem_type_idx, elem_type in enumerate(surf_elements[0]):
        elem_name = gmsh.model.mesh.getElementProperties(elem_type)[0]

        for elem_tag in surf_elements[1][elem_type_idx]:
            # Get element nodes
            elem_nodes = gmsh.model.mesh.getElement(elem_type, elem_tag)[1]

            # Get node coordinates
            coords = []
            for node in elem_nodes:
                coord = gmsh.model.mesh.getNode(node)[0]
                coords.extend(coord)

            # Add scalar value (constant current density)
            coords.append(current_density)
            data_list.extend(coords)

    # Add data to view
    if "Triangle" in elem_name or "Tri" in elem_name:
        gmsh.view.addListData(view, "ST", len(surf_elements[1][0]), data_list)
    elif "Quad" in elem_name:
        gmsh.view.addListData(view, "SQ", len(surf_elements[1][0]), data_list)

    # Configure colormap
    gmsh.view.option.setNumber(view, "ColormapNumber", 2)  # Jet colormap
    gmsh.view.option.setNumber(view, "RangeType", 2)       # Custom range
    gmsh.view.option.setNumber(view, "CustomMin", 0)
    gmsh.view.option.setNumber(view, "CustomMax", current_density * 1.1)

    # Save with field data
    msh_field_file = 'circular_coil_with_field.msh'
    gmsh.write(msh_field_file)
    print(f"    Saved: {msh_field_file}")
    print(f"    Current density: {current_density:.2e} A/m^2")

    # Show in GMSH GUI
    print("\n[11] Opening GMSH GUI with field data...")
    print("    View controls:")
    print("      - Press '1' to toggle current density view")
    print("      - Tools -> Visibility: Control view display")
    print("      - Tools -> Options -> View: Colormap settings")

    gmsh.fltk.run()
    gmsh.finalize()


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Create circular coil in GMSH')
    parser.add_argument('--radius', type=float, default=50.0,
                        help='Mean radius (mm, default: 50)')
    parser.add_argument('--width', type=float, default=4.0,
                        help='Wire width (mm, default: 4)')
    parser.add_argument('--height', type=float, default=4.0,
                        help='Wire height (mm, default: 4)')
    parser.add_argument('--mesh-size', type=float, default=1.0,
                        help='Mesh size (mm, default: 1)')
    parser.add_argument('--with-field', action='store_true',
                        help='Add current density field data')
    parser.add_argument('--current', type=float, default=1.5e6,
                        help='Current density (A/m^2, default: 1.5e6)')
    parser.add_argument('--no-gui', action='store_true',
                        help='Do not show GMSH GUI')

    args = parser.parse_args()

    # Convert mm to m
    mean_radius = args.radius / 1000.0
    wire_width = args.width / 1000.0
    wire_height = args.height / 1000.0
    mesh_size = args.mesh_size / 1000.0

    if args.with_field:
        create_coil_with_field_data(mean_radius, wire_width, wire_height,
                                     args.current)
    else:
        create_circular_coil(mean_radius, wire_width, wire_height,
                            mesh_size, show_gui=not args.no_gui)

    print("\n" + "=" * 60)
    print("Next Steps")
    print("=" * 60)
    print("\n1. Import to NGSolve:")
    print("   from ngsolve import Mesh")
    print("   mesh = Mesh('circular_coil.msh')")

    print("\n2. Convert to Radia PEEC (future API):")
    print("   from peec_mesh_import import surface_mesh_to_peec")
    print("   conductor = surface_mesh_to_peec(mesh, sigma=5.8e7)")

    print("\n3. Current workaround - Use CndLoop:")
    print("   coil = rad.CndLoop([0,0,0], 0.05, [0,0,1], 'r',")
    print("                       0.004, 0.004, 5.8e7, 8, 36)")


if __name__ == '__main__':
    main()
