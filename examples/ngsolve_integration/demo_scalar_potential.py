"""Demo: Scalar potential solver combining Radia (MMM/MSC) + NGSolve FEM.

This example demonstrates the two-scalar-potential (phi-reduced/phi)
method for computing B field from permanent magnets in the presence
of soft iron.

Problem setup:
    - NdFeB permanent magnet (rectangular block)
    - Soft iron yoke nearby
    - Air surrounding both

Method:
    1. Radia computes H_s from the permanent magnet (MMM/MSC integral
       method, analytical for ObjRecMag)
    2. NGSolve solves for the correction potential phi that accounts
       for the iron magnetization (FEM)
    3. Total field: H = H_s - grad(phi), B = mu * H

Reference:
    Simkin & Trowbridge, IJNME, 1979.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import radia as rad


def build_radia_magnet():
    """Create a permanent magnet in Radia."""
    rad.UtiDelAll()

    # NdFeB magnet: 10mm x 10mm x 10mm, Br = 1.2 T -> M = 954930 A/m
    center = [0, 0, 0]
    size = [0.01, 0.01, 0.01]
    magnetization = [0, 0, 954930]  # A/m, z-directed
    mag = rad.ObjRecMag(center, size, magnetization)

    return mag


def build_mesh_with_iron():
    """Create NGSolve mesh with an iron block near the magnet."""
    from netgen.occ import Box, Pnt, OCCGeometry, Glue
    from ngsolve import Mesh, TaskManager

    # Air box: 100mm x 100mm x 100mm centered at origin
    air_box = Box(Pnt(-0.05, -0.05, -0.05), Pnt(0.05, 0.05, 0.05))

    # Iron yoke: 10mm x 10mm x 10mm, positioned 5mm above magnet
    iron_box = Box(Pnt(-0.005, -0.005, 0.01), Pnt(0.005, 0.005, 0.02))

    # Subtract iron from air to get clean interface
    air_region = air_box - iron_box

    # Assign material names
    air_region.name = 'air'
    iron_box.name = 'iron'

    # Glue together
    shape = Glue([air_region, iron_box])

    # Set outer boundary name
    for face in air_box.faces:
        face.name = 'outer'

    geo = OCCGeometry(shape)
    with TaskManager():
        mesh = geo.GenerateMesh(maxh=0.01)
        return Mesh(mesh)


def main():
    print("=== Scalar Potential Solver: Radia + NGSolve ===")
    print()

    # Step 1: Create Radia permanent magnet model
    print("1. Building Radia permanent magnet model...")
    mag = build_radia_magnet()

    # Verify Radia field at a test point
    B_radia = rad.Fld(mag, 'b', [0, 0, 0.02])
    print(f"   Radia B at (0,0,20mm): [{B_radia[0]:.6f}, "
          f"{B_radia[1]:.6f}, {B_radia[2]:.6f}] T")

    # Step 2: Create mesh
    print("2. Creating NGSolve mesh with iron region...")
    mesh = build_mesh_with_iron()
    print(f"   Materials: {list(mesh.GetMaterials())}")
    print(f"   Boundaries: {list(mesh.GetBoundaries())[:5]}...")
    print(f"   Elements: {mesh.ne}")

    # Step 3: Solve with single-potential method
    print("3. Solving with single-potential method (mu_r=1000)...")
    from radia.scalar_potential_solver import ScalarPotentialSolver

    solver = ScalarPotentialSolver(mesh, iron_domains='iron',
                                  mu_r=1000, order=2)
    solver.set_source_from_radia(mag, resolution=31)
    phi_gf = solver.solve(method='single')

    # Step 4: Get results
    print("4. Getting results...")
    B_cf = solver.get_B()
    H_cf = solver.get_H()

    # Evaluate at a test point in air (above the iron)
    from ngsolve import Integrate, CF, sqrt

    # Compare Radia B vs FEM-corrected B at several points
    test_points = [
        [0, 0, 0.03],    # above iron
        [0.02, 0, 0],    # side of magnet
        [0, 0, -0.02],   # below magnet
    ]

    print()
    print("   Comparison: Radia (no iron) vs FEM-corrected (with iron)")
    print("   " + "-" * 65)
    print(f"   {'Point':25s} {'Radia |B|':>12s} {'FEM |B|':>12s} {'Ratio':>10s}")
    print("   " + "-" * 65)

    for pt in test_points:
        B_r = rad.Fld(mag, 'b', pt)
        B_r_mag = np.sqrt(sum(b**2 for b in B_r[:3]))

        # Evaluate FEM field at point using mesh
        try:
            mip = mesh(pt[0], pt[1], pt[2])
            B_fem = B_cf(mip)
            B_fem_mag = np.sqrt(sum(b**2 for b in B_fem))
            ratio = B_fem_mag / B_r_mag if B_r_mag > 0 else float('inf')
            print(f"   ({pt[0]*1000:5.1f},{pt[1]*1000:5.1f},{pt[2]*1000:5.1f}) mm"
                  f"    {B_r_mag:10.6f} T  {B_fem_mag:10.6f} T  {ratio:8.4f}")
        except Exception as e:
            print(f"   ({pt[0]*1000:5.1f},{pt[1]*1000:5.1f},{pt[2]*1000:5.1f}) mm"
                  f"    {B_r_mag:10.6f} T  (point outside mesh)")

    # Step 5: Project to HDiv (optional)
    print()
    print("5. Projecting B to HDiv space...")
    B_hdiv = solver.project_to_hdiv()
    print(f"   HDiv DOFs: {B_hdiv.space.ndof}")

    # Step 6: VTK export (optional)
    print("6. Exporting to VTK...")
    from ngsolve import VTKOutput
    from ngsolve import TaskManager
    vtk_file = os.path.join(os.path.dirname(__file__),
                            'demo_scalar_potential')
    try:
        vtk = VTKOutput(mesh, coefs=[B_cf, H_cf],
                        names=['B', 'H'],
                        filename=vtk_file)
        vtk.Do()
        print(f"   Written: {vtk_file}.vtu")
    except Exception as e:
        print(f"   VTK export skipped: {e}")

    print()
    print("=== Done ===")
    print()
    print("The iron yoke concentrates the magnetic flux.")
    print("B field near the iron should be enhanced compared to")
    print("the bare magnet (Radia-only) result.")

    rad.UtiDelAll()


if __name__ == '__main__':
    main()
