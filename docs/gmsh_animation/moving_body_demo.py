"""
Demo: Stator (STEP, fixed) + Rotor (.vol, translating) -> GMSH animation.

Creates:
  1. A stator geometry (STEP file, static)
  2. A rotor mesh (Netgen .vol with curving)
  3. A GMSH .msh v4.1 with multiple time steps (rotor displaced at each step)
  4. A .geo companion file for GMSH settings

The resulting .msh can be opened in GMSH and animated using
  Tools -> Options -> General -> Time step slider
or GMSH's built-in animation (View -> Animation).

Usage:
  python moving_body_demo.py

Requirements:
  pip install ngsolve
  gmsh (for visualization)
"""

import math
import os
import sys
import numpy as np

# Add radia to path
_this_dir = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.dirname(os.path.dirname(_this_dir))
sys.path.insert(0, os.path.join(_repo_root, "src"))
sys.path.insert(0, _this_dir)

from ngsolve import Mesh, CF, Integrate, BND, TaskManager
from netgen.csg import CSGeometry, Cylinder, Plane, Pnt, Vec, OrthoBrick

from tet10_orientation_fix import fix_msh_tet_orientation


def create_stator_step(filename):
    """Create a simple stator-like geometry and export as STEP.

    Stator = hollow cylinder (outer R=0.08, inner R=0.052, height=0.1).
    """
    try:
        from netgen.occ import OCCGeometry, Box, Cylinder as OCCCylinder, Axes
        ax = Axes(Pnt(0, 0, 0), Vec(0, 0, 1))
        outer = OCCCylinder(ax, 0.08, 0.1)
        inner = OCCCylinder(ax, 0.052, 0.1)
        stator = outer - inner
        geo = OCCGeometry(stator)
        geo.Glue()
        # Export STEP
        from netgen.occ import step as occ_step
        occ_step.OCCGeometryWrite(geo, filename)
        return True
    except Exception:
        pass

    # Fallback: create STEP via simple CSG description file
    # (STEP export requires OCC which may not be available in fork)
    return False


def create_rotor_vol(filename, order=2):
    """Create a rotor mesh (solid cylinder R=0.05, H=0.1) as .vol.

    Uses Netgen CSG (no Cubit needed).
    """
    geo = CSGeometry()
    cyl = Cylinder(Pnt(0, 0, 0), Pnt(0, 0, 0.1), 0.05)
    bot = Plane(Pnt(0, 0, 0), Vec(0, 0, -1))
    top = Plane(Pnt(0, 0, 0.1), Vec(0, 0, 1))
    rotor = cyl * bot * top
    geo.Add(rotor.bc("rotor_surface").mat("rotor"))

    with TaskManager():
        ng_mesh = geo.GenerateMesh(maxh=0.015)
        ng_mesh.Curve(order)
        ng_mesh.Save(filename)
        return filename


def create_animation_msh(rotor_vol, output_msh, n_steps=20,
                          displacement=0.15, direction=(1, 0, 0)):
    """Create GMSH .msh v4.1 with time-stepped displaced rotor mesh.

    Each time step shifts the rotor by (step/n_steps * displacement * direction).
    Displacement field is stored as $NodeData for visualization.

    Args:
        rotor_vol: Path to rotor .vol file.
        output_msh: Output .msh path.
        n_steps: Number of animation frames.
        displacement: Total displacement magnitude [m].
        direction: Displacement direction (normalized).
    """
    mesh = Mesh(rotor_vol)
    ne = mesh.ne

    # Normalize direction
    d = np.array(direction, dtype=float)
    d = d / np.linalg.norm(d)

    from radia.gmsh_post_export import GmshPostExport

    print(f"Rotor: {ne} elements")
    print(f"Displacement: {displacement} m in direction {d}")
    print(f"Frames: {n_steps}")

    # Write the mesh structure (TET10 for a Curve(2) .vol)
    with TaskManager():
        post = GmshPostExport(mesh)
        post.write_mesh(output_msh)

    # GmshPostExport writes GMSH tet corners in NGSolve el.vertices order,
    # which is an odd (inverted) permutation of the GMSH-positive ordering
    # (all-negative getJacobians). Reorient in place; see
    # tet10_orientation_fix.py. The root fix belongs in
    # src/radia/gmsh_post_export.py (_NGSOLVE_TO_GMSH_NODE_ORDER['TET']).
    stats = fix_msh_tet_orientation(output_msh)
    print(f"Orientation fix: {stats['fixed']} of {stats['elements']} "
          f"tets reoriented")

    # Displacement NodeData must cover ALL nodes of the written .msh
    # (vertices + TET10 mid-edge nodes), not just mesh.nv vertices --
    # otherwise GMSH leaves mid-edge nodes behind in displacement views.
    with open(output_msh, "r") as f:
        lines = f.read().splitlines()
    n_msh_nodes = int(lines[lines.index("$Nodes") + 1].split()[1])

    # Append $NodeData sections for each time step
    with open(output_msh, "a") as f:
        for step in range(n_steps + 1):
            t = step / n_steps
            dx = displacement * t * d

            # Displacement field (3-component vector per node)
            f.write("$NodeData\n")
            f.write("1\n")  # num string tags
            f.write(f'"Displacement [m]"\n')
            f.write("1\n")  # num real tags
            f.write(f"{t:.6f}\n")  # time value
            f.write("3\n")  # num int tags
            f.write(f"{step}\n")  # timestep
            f.write("3\n")  # num components (vector)
            f.write(f"{n_msh_nodes}\n")  # num nodes

            for vid in range(n_msh_nodes):
                f.write(f"{vid + 1} {dx[0]:.10e} {dx[1]:.10e} {dx[2]:.10e}\n")

            f.write("$EndNodeData\n")

    print(f"Written: {output_msh} ({n_steps + 1} time steps)")
    return output_msh


def create_geo_file(output_geo, stator_step, rotor_msh):
    """Create .geo companion file for GMSH visualization settings."""
    with open(output_geo, "w") as f:
        f.write(f'// Animation: stator (STEP) + moving rotor (.msh)\n')
        f.write(f'// Open this file in GMSH, then use View -> Animation\n\n')
        f.write(f'Merge "{os.path.basename(stator_step)}";\n')
        f.write(f'Merge "{os.path.basename(rotor_msh)}";\n\n')
        f.write(f'// Display settings\n')
        f.write(f'Mesh.SurfaceEdges = 1;\n')
        f.write(f'Mesh.VolumeEdges = 0;\n')
        f.write(f'Mesh.NumSubEdges = 4;  // curved element rendering\n\n')
        f.write(f'// Displacement view settings\n')
        f.write(f'View[0].VectorType = 5;  // displacement\n')
        f.write(f'View[0].DisplacementFactor = 1.0;\n')
        f.write(f'View[0].ShowElement = 1;\n')
        f.write(f'View[0].ShowScale = 0;\n\n')
        f.write(f'// Animation\n')
        f.write(f'General.TranslationX = 0;\n')
        f.write(f'General.TranslationY = 0;\n')
    print(f"Written: {output_geo}")


def main():
    out_dir = os.path.join(_this_dir)
    os.makedirs(out_dir, exist_ok=True)

    rotor_vol = os.path.join(out_dir, "rotor.vol")
    rotor_msh = os.path.join(out_dir, "rotor_animation.msh")
    stator_step = os.path.join(out_dir, "stator.step")
    geo_file = os.path.join(out_dir, "animation.geo")

    print("=" * 60)
    print("Stator + Moving Rotor Animation Demo")
    print("=" * 60)

    # Create rotor mesh
    print("\n1. Creating rotor mesh (.vol with curving)...")
    create_rotor_vol(rotor_vol, order=2)
    m = Mesh(rotor_vol)
    vol = Integrate(CF(1), m)
    area = Integrate(CF(1), m, BND)
    print(f"   Rotor: {m.ne} elements, V={vol:.6e} m^3, A={area:.6e} m^2")

    # Create stator STEP (if OCC available)
    print("\n2. Creating stator geometry (STEP)...")
    has_step = create_stator_step(stator_step)
    if has_step:
        print(f"   Stator STEP: {stator_step}")
    else:
        print("   OCC not available -- stator STEP skipped")
        print("   (demo works without stator, rotor animation still visible)")

    # Create animation .msh
    print("\n3. Creating animation .msh (GMSH v4.1)...")
    create_animation_msh(rotor_vol, rotor_msh,
                          n_steps=20, displacement=0.15,
                          direction=(1, 0, 0))

    # Create .geo companion
    print("\n4. Creating .geo companion...")
    if has_step:
        create_geo_file(geo_file, stator_step, rotor_msh)
    else:
        # Just rotor
        with open(geo_file, "w") as f:
            f.write(f'Merge "{os.path.basename(rotor_msh)}";\n')
            f.write(f'Mesh.NumSubEdges = 4;\n')
            f.write(f'View[0].VectorType = 5;\n')
            f.write(f'View[0].DisplacementFactor = 1.0;\n')
            f.write(f'View[0].ShowElement = 1;\n')
            f.write(f'View[0].ShowScale = 0;\n')
        print(f"   Written: {geo_file}")

    print("\n" + "=" * 60)
    print("To view animation:")
    print(f"  gmsh {geo_file}")
    print("  -> Tools -> Options -> General -> Time step slider")
    print("  -> or View -> Animation -> Play")
    print("=" * 60)


if __name__ == "__main__":
    main()
