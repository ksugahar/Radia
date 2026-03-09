"""
Demo: C-type electromagnet hex mesh using GmshBuilder.

Creates a C-shaped iron core with structured hexahedral mesh,
then solves with Radia magnetostatics.

Uses GMSH OCC kernel as backend with Gmsh-based high-level API.

Requires: pip install gmsh
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

from radia.gmsh_builder import GmshBuilder


def demo_c_magnet_mesh():
    """Create C-type magnet with hex mesh."""

    # C-type magnet dimensions (meters)
    # Bottom yoke: 100mm x 20mm x 30mm
    # Left pole:   20mm x 20mm x 20mm
    # Right pole:  20mm x 20mm x 20mm
    yoke_length = 0.10
    yoke_height = 0.02
    yoke_depth = 0.03
    pole_width = 0.02
    pole_height = 0.02
    gap = 0.02  # air gap between poles

    with GmshBuilder(verbose=True) as hb:
        # 1. Create yoke (bottom bar)
        yoke = hb.add_box(
            center=[0, 0, 0],
            size=[yoke_length, yoke_depth, yoke_height])

        # 2. Create left pole
        left_pole = hb.add_box(
            center=[-(yoke_length / 2 - pole_width / 2), 0,
                    yoke_height / 2 + pole_height / 2],
            size=[pole_width, yoke_depth, pole_height])

        # 3. Create right pole
        right_pole = hb.add_box(
            center=[(yoke_length / 2 - pole_width / 2), 0,
                    yoke_height / 2 + pole_height / 2],
            size=[pole_width, yoke_depth, pole_height])

        # 4. Set mesh divisions
        hb.set_divisions(yoke, nx=10, ny=3, nz=2)
        hb.set_divisions(left_pole, nx=2, ny=3, nz=2)
        hb.set_divisions(right_pole, nx=2, ny=3, nz=2)

        # 5. Generate hex mesh
        hb.generate()

        # 6. Show mesh info
        info = hb.get_mesh_info()
        print(f"\nMesh statistics:")
        print(f"  Nodes: {info['n_nodes']}")
        print(f"  Hex elements: {info['n_hex']}")
        print(f"  Tet elements: {info['n_tet']}")

        # 7. Export
        output_file = os.path.join(os.path.dirname(__file__),
                                    'c_magnet.msh')
        hb.export(output_file)
        print(f"\nExported: {output_file}")

        # 8. Convert to Radia (optional - needs radia)
        try:
            import radia as rad
            rad.UtiDelAll()
            container = hb.to_radia(mu_r=1000)
            print(f"\nRadia container: {container}")
            print("C-type magnet ready for Radia.Solve()")
        except ImportError:
            print("\nRadia not installed - mesh export only")


def demo_webcut():
    """Demo webcut: split a box into 4 pieces."""
    print("\n--- Webcut Demo ---")

    with GmshBuilder(verbose=True) as hb:
        # Create a box
        box = hb.add_box([0, 0, 0], [0.06, 0.06, 0.06])

        # Webcut along x (split into left/right)
        parts_x = hb.webcut_plane(box, 'x', 0.0)
        print(f"After x-cut: {len(parts_x)} pieces")

        # Webcut left piece along z (split into bottom/top)
        parts_z = hb.webcut_plane(parts_x[0], 'z', 0.0)
        print(f"After z-cut: {len(hb.get_volumes())} total pieces")

        # Mesh all pieces
        hb.set_divisions_all(3, 3, 3)
        hb.generate()

        info = hb.get_mesh_info()
        print(f"Total hex: {info['n_hex']} "
              f"(3 pieces x 27 each = {3 * 27})")


if __name__ == '__main__':
    demo_c_magnet_mesh()
    demo_webcut()
