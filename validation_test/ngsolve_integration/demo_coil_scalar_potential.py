"""Demo: Coil source field with phi-reduced scalar potential (Radia + NGSolve).

This example demonstrates the phi-reduced scalar potential method where
the source field H_s comes from a **coil** (not a permanent magnet).

Physical model:
    - Rectangular solenoid coil (100 turns, 1A)
    - Soft iron core inside the coil bore
    - Air surrounding both

Method:
    1. Model the coil in Radia using equivalent magnetization:
       A solenoid with N turns, current I, height h has M = NI/h (A/m)
       along the axis.  ObjRecMag blocks with this M produce exactly
       the same external field as the real coil winding.
    2. Radia computes H_s (coil field in free space, no iron).
    3. NGSolve solves for the correction potential phi that accounts
       for the iron core magnetization.
    4. Total field: H = H_s - grad(phi), B = mu * H.

The key physics:
    The jump [H_s . n] across the iron-air interface drives the FEM
    correction potential.  Without iron, H_s is the complete solution.
    With iron, the correction potential redistributes the flux through
    the high-permeability core.

Reference:
    Simkin & Trowbridge, IJNME, vol. 14, pp. 423-440, 1979.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import radia as rad


def build_radia_coil(N=100, I=1.0, h=0.020,
                     bore_xy=0.020, outer_xy=0.030):
    """Create a rectangular solenoid coil in Radia using equivalent magnetization.

    A solenoid winding with N turns carrying current I over height h
    is equivalent to a uniformly magnetized hollow prism with
    M_z = N * I / h (A/m).

    Parameters
    ----------
    N : int
        Number of turns.
    I : float
        Current in amperes.
    h : float
        Coil height in meters.
    bore_xy : float
        Inner bore dimension (square cross-section) in meters.
    outer_xy : float
        Outer winding dimension (square cross-section) in meters.

    Returns
    -------
    int
        Radia container object handle.
    """
    rad.UtiDelAll()

    M_z = N * I / h  # Equivalent magnetization (A/m)
    mag = [0, 0, M_z]

    a_in = bore_xy / 2.0
    a_out = outer_xy / 2.0
    h2 = h / 2.0

    # Decompose hollow rectangular prism into 4 wall blocks
    # Left wall: x in [-a_out, -a_in], y in [-a_out, a_out]
    left = rad.magnet_box([-(a_in + a_out) / 2, 0, 0],
                         [a_out - a_in, outer_xy, h], mag)
    # Right wall: x in [a_in, a_out], y in [-a_out, a_out]
    right = rad.magnet_box([(a_in + a_out) / 2, 0, 0],
                          [a_out - a_in, outer_xy, h], mag)
    # Front wall: x in [-a_in, a_in], y in [a_in, a_out]
    front = rad.magnet_box([0, (a_in + a_out) / 2, 0],
                          [bore_xy, a_out - a_in, h], mag)
    # Back wall: x in [-a_in, a_in], y in [-a_out, -a_in]
    back = rad.magnet_box([0, -(a_in + a_out) / 2, 0],
                         [bore_xy, a_out - a_in, h], mag)

    coil = rad.ObjCnt([left, right, front, back])
    return coil


def build_mesh_with_iron_core(bore_xy=0.020, outer_xy=0.030,
                              coil_h=0.020, core_xy=0.018,
                              core_h=0.030, domain=0.08, maxh=0.008):
    """Create NGSolve mesh with an iron core inside the coil bore.

    Parameters
    ----------
    bore_xy : float
        Coil bore dimension (m).
    outer_xy : float
        Coil outer dimension (m).
    coil_h : float
        Coil height (m).
    core_xy : float
        Iron core cross-section (m), must be < bore_xy.
    core_h : float
        Iron core length (m), can extend beyond coil.
    domain : float
        Air domain half-size (m).
    maxh : float
        Maximum mesh element size (m).

    Returns
    -------
    ngsolve.Mesh
    """
    from netgen.occ import Box, Pnt, OCCGeometry, Glue

    # Air box
    air_box = Box(Pnt(-domain, -domain, -domain),
                  Pnt(domain, domain, domain))

    # Iron core (centered, extends in z beyond the coil)
    c2 = core_xy / 2.0
    ch2 = core_h / 2.0
    iron_box = Box(Pnt(-c2, -c2, -ch2), Pnt(c2, c2, ch2))

    # Subtract iron from air
    air_region = air_box - iron_box
    air_region.name = 'air'
    iron_box.name = 'iron'

    shape = Glue([air_region, iron_box])

    geo = OCCGeometry(shape)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh)

        from ngsolve import Mesh
        return Mesh(ngmesh)


def main():
    print("=== Coil Source Field: Phi-Reduced Scalar Potential ===")
    print()

    # Coil parameters
    N = 100       # turns
    I_coil = 1.0  # amperes
    h_coil = 0.020  # 20 mm height
    bore = 0.020    # 20 mm bore (square)
    outer = 0.030   # 30 mm outer (square)
    M_z = N * I_coil / h_coil

    print("Coil parameters:")
    print(f"  Turns: {N}, Current: {I_coil} A, NI = {N * I_coil} A-turns")
    print(f"  Height: {h_coil * 1000:.0f} mm")
    print(f"  Bore: {bore * 1000:.0f} x {bore * 1000:.0f} mm")
    print(f"  Outer: {outer * 1000:.0f} x {outer * 1000:.0f} mm")
    print(f"  Equivalent magnetization M_z = {M_z:.0f} A/m")
    print()

    # Iron core parameters
    core_xy = 0.018   # 18 mm (slightly smaller than bore)
    core_h = 0.030    # 30 mm (extends 5mm beyond coil on each side)
    mu_r = 1000.0

    print("Iron core parameters:")
    print(f"  Cross-section: {core_xy * 1000:.0f} x {core_xy * 1000:.0f} mm")
    print(f"  Length: {core_h * 1000:.0f} mm")
    print(f"  mu_r = {mu_r:.0f}")
    print()

    # Step 1: Build Radia coil model
    print("1. Building Radia coil model (equivalent magnetization)...")
    coil = build_radia_coil(N, I_coil, h_coil, bore, outer)

    # Verify coil field at test points (no iron)
    test_pts = [
        [0, 0, 0],        # center of coil
        [0, 0, 0.015],    # near top of core
        [0, 0, 0.025],    # above coil
        [0.02, 0, 0],     # side (outside bore)
    ]
    print("   Coil field H_s (no iron):")
    for pt in test_pts:
        B = rad.Fld(coil, 'b', pt)
        H = rad.Fld(coil, 'h', pt)
        B_mag = np.sqrt(sum(b ** 2 for b in B[:3]))
        H_mag = np.sqrt(sum(v ** 2 for v in H[:3]))
        label = f"({pt[0]*1000:5.1f},{pt[1]*1000:5.1f},{pt[2]*1000:5.1f}) mm"
        print(f"   {label}  |B|={B_mag:.6f} T  |H|={H_mag:.1f} A/m")

    # Step 2: Create mesh
    print()
    print("2. Creating NGSolve mesh with iron core...")
    mesh = build_mesh_with_iron_core(
        bore_xy=bore, outer_xy=outer, coil_h=h_coil,
        core_xy=core_xy, core_h=core_h,
        domain=0.06, maxh=0.006
    )
    print(f"   Materials: {list(mesh.GetMaterials())}")
    print(f"   Elements: {mesh.ne}")

    # Step 3: Solve with scalar potential
    print()
    print("3. Solving phi-reduced scalar potential (Radia coil H_s as source)...")
    from radia.scalar_potential_solver import ScalarPotentialSolver

    solver = ScalarPotentialSolver(mesh, iron_domains='iron',
                                  mu_r=mu_r, order=2)
    solver.set_source_from_radia(coil, resolution=31)
    phi_gf = solver.solve(method='single')

    # Step 4: Compare fields
    print()
    print("4. Comparing coil-only vs FEM-corrected (with iron core):")
    B_cf = solver.get_B()
    H_cf = solver.get_H()

    from ngsolve import CF, sqrt as ngsqrt

    print()
    print("   " + "-" * 70)
    print(f"   {'Point':30s} {'Coil |B|':>12s} {'FEM |B|':>12s} {'Ratio':>10s}")
    print("   " + "-" * 70)

    compare_pts = [
        ([0, 0, 0],       "center of coil"),
        ([0, 0, 0.010],   "mid-core, z=10mm"),
        ([0, 0, 0.020],   "core end, z=20mm"),
        ([0, 0, 0.030],   "above core, z=30mm"),
        ([0.02, 0, 0],    "side, x=20mm"),
        ([0, 0, -0.020],  "below core, z=-20mm"),
    ]

    for pt, desc in compare_pts:
        B_radia = rad.Fld(coil, 'b', pt)
        B_r_mag = np.sqrt(sum(b ** 2 for b in B_radia[:3]))

        try:
            mip = mesh(pt[0], pt[1], pt[2])
            B_fem = B_cf(mip)
            B_fem_mag = np.sqrt(sum(b ** 2 for b in B_fem))
            ratio = B_fem_mag / B_r_mag if B_r_mag > 1e-12 else float('inf')
            print(f"   {desc:30s} {B_r_mag:10.6f} T  {B_fem_mag:10.6f} T  {ratio:8.3f}x")
        except Exception:
            print(f"   {desc:30s} {B_r_mag:10.6f} T  (outside mesh)")

    # Step 5: Project to HDiv
    print()
    print("5. Projecting B to HDiv space (div(B)=0)...")
    B_hdiv = solver.project_to_hdiv()
    print(f"   HDiv DOFs: {B_hdiv.space.ndof}")

    # Step 6: GMSH export
    print()
    print("6. Exporting to GMSH...")
    from radia.gmsh_post_export import GmshPostExport
    gmsh_file = os.path.join(os.path.dirname(__file__),
                             'demo_coil_scalar_potential.msh')
    try:
        post = GmshPostExport(mesh)
        post.add_vector_field('B', B_cf)
        post.add_vector_field('H', H_cf)
        post.write(gmsh_file)
        print(f"   Written: {gmsh_file}")
    except Exception as e:
        print(f"   GMSH export skipped: {e}")

    print()
    print("=== Done ===")
    print()
    print("The iron core concentrates flux inside the coil bore.")
    print("At the coil center, B should be significantly enhanced")
    print("(roughly mu_r / (mu_r + demagnetization) times the coil-only field).")
    print("Outside the core, the field is also modified by the iron.")

    rad.UtiDelAll()


if __name__ == '__main__':
    main()
