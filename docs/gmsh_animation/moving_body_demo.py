"""
Demo: Stator (STEP, fixed) + Rotor (.vol, translating) -> GMSH animation.

Creates:
  1. A stator geometry (STEP file, static)
  2. A rotor mesh (Netgen .vol with curving)
  3. A GMSH .msh v4.2 with multiple time steps (rotor displaced at each step)
The tracked animation.geo (+ .geo.opt) merges the STEP and the .msh for display.

The resulting .msh can be opened in GMSH and animated using
  Tools -> Options -> General -> Time step slider
or GMSH's built-in animation (View -> Animation).

Usage:
  python moving_body_demo.py

Requirements:
  ngsolve (with netgen.occ for the STEP), radia
  Gmsh 5 + MSH 4.2 (for visualization)
"""

import math
import os
import numpy as np

_this_dir = os.path.dirname(os.path.abspath(__file__))

from ngsolve import Mesh, CF, Integrate, BND, TaskManager
from netgen.csg import CSGeometry, Cylinder, Plane, Pnt, Vec, OrthoBrick

from tet10_orientation_fix import fix_msh_tet_orientation


def create_stator_step(filename):
    """Create a simple stator-like geometry and export as STEP.

    Stator = hollow cylinder (outer R=0.08, inner R=0.052, height=0.1).
    """
    from netgen.occ import Cylinder as OCCCylinder, Pnt as OPnt, Z
    stator = OCCCylinder(OPnt(0, 0, 0), Z, 0.08, 0.1) - OCCCylinder(OPnt(0, 0, 0), Z, 0.052, 0.1)
    stator.WriteStep(filename)
    if not os.path.isfile(filename) or os.path.getsize(filename) == 0:
        raise RuntimeError(f"stator STEP was not written: {filename}")


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
    """Create GMSH .msh v4.2 with time-stepped displaced rotor mesh.

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

    from radia.gmsh_post_export import GmshPostExport, DataField, _emit_data_section

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

            _emit_data_section(f, DataField(
                kind='$NodeData', name='Displacement [m]', ncomp=3,
                items=[(vid + 1, tuple(dx)) for vid in range(n_msh_nodes)],
                time=t, timestep=step))

    print(f"Written: {output_msh} ({n_steps + 1} time steps)")
    return output_msh


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

    print("\n2. Creating stator geometry (STEP)...")
    create_stator_step(stator_step)
    print(f"   Stator STEP: {stator_step}")

    print("\n3. Creating animation .msh (GMSH v4.2)...")
    create_animation_msh(rotor_vol, rotor_msh,
                          n_steps=20, displacement=0.15,
                          direction=(1, 0, 0))

    # The tracked animation.geo (+ .geo.opt) is the curated launch file; it is not rewritten here.
    if not os.path.isfile(geo_file):
        raise FileNotFoundError(f"tracked launch file missing: {geo_file}")

    print("\n" + "=" * 60)
    print("To view animation:")
    print(f"  gmsh {geo_file}")
    print("  -> Tools -> Options -> General -> Time step slider")
    print("  -> or View -> Animation -> Play")
    print("=" * 60)


if __name__ == "__main__":
    main()
