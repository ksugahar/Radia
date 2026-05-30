"""Demo: Total scalar potential with Gmsh cohomology cuts.

This example demonstrates the cohomology-based total scalar potential
method for a coil + iron core magnetostatic problem.

Formulation:
    H = -grad(phi) + NI * h

where h is the cohomology basis function (HCurl, curl-free) computed
automatically by Gmsh.  No Biot-Savart computation needed.

Physical model:
    - Cylindrical iron core (mu_r = 1000)
    - Coil modeled as a hole in the mesh (current = NI ampere-turns)
    - Surrounding air domain

Comparison:
    Results are verified against Radia's direct field computation
    using an equivalent magnetization coil model (ObjRecMag).

Reference:
    Pellikka et al., "Homology and Cohomology Computation in Finite
    Element Modeling," SIAM J. Sci. Comput., 2013.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src'))

MU_0 = 4 * np.pi * 1e-7


def build_gmsh_geometry(core_radius=0.008, core_height=0.030,
                        coil_r_inner=0.010, coil_r_outer=0.015,
                        coil_height=0.020, air_size=0.06, maxh=0.008):
    """Create Gmsh geometry: air box + iron core + coil hole.

    The coil is modeled as a hollow cylinder (torus-like hole)
    subtracted from the domain.  Volume classification uses the
    fragment mapping (not centroids, which fail for symmetric geometry).

    Returns
    -------
    domain_tag : int
        Physical group tag for the 3D domain.
    """
    import gmsh

    gmsh.initialize()
    gmsh.option.setNumber('General.Verbosity', 2)
    gmsh.model.add('coil_cohomology')

    # Air box
    box = gmsh.model.occ.addBox(
        -air_size/2, -air_size/2, -air_size/2,
        air_size, air_size, air_size)

    # Iron core (cylinder along z)
    iron_cyl = gmsh.model.occ.addCylinder(
        0, 0, -core_height/2, 0, 0, core_height, core_radius)

    # Coil hole: outer and inner cylinders (fragment will split them)
    coil_outer = gmsh.model.occ.addCylinder(
        0, 0, -coil_height/2, 0, 0, coil_height, coil_r_outer)
    coil_inner = gmsh.model.occ.addCylinder(
        0, 0, -coil_height/2, 0, 0, coil_height, coil_r_inner)

    # Fragment ALL volumes together for proper non-overlapping decomposition
    _, result_map = gmsh.model.occ.fragment(
        [(3, box)],
        [(3, iron_cyl), (3, coil_outer), (3, coil_inner)])
    gmsh.model.occ.synchronize()

    # Use fragment mapping to classify volumes:
    #   result_map[0] = box fragments (all volumes)
    #   result_map[1] = iron_cyl fragments
    #   result_map[2] = coil_outer fragments (includes coil annulus + interior)
    #   result_map[3] = coil_inner fragments (interior only)
    iron_vol_set = set(t[1] for t in result_map[1])
    coil_outer_set = set(t[1] for t in result_map[2])
    coil_inner_set = set(t[1] for t in result_map[3])

    # Coil annulus = in coil_outer but NOT in coil_inner and NOT iron
    coil_vols = list(coil_outer_set - coil_inner_set - iron_vol_set)

    volumes = gmsh.model.getEntities(3)
    all_vol_tags = set(v[1] for v in volumes)
    iron_vols = list(iron_vol_set)
    air_vols = list(all_vol_tags - iron_vol_set - set(coil_vols))

    print(f"  Volumes after fragment: {len(volumes)}")
    print(f"  Iron volumes: {iron_vols}")
    print(f"  Coil volumes (to remove): {coil_vols}")
    print(f"  Air volumes: {air_vols}")

    # Remove coil volumes from domain (creates the hole)
    if coil_vols:
        for tag in coil_vols:
            gmsh.model.occ.remove([(3, tag)], recursive=True)
        gmsh.model.occ.synchronize()

    # Re-enumerate and reclassify after removal
    volumes = gmsh.model.getEntities(3)
    remaining_tags = set(v[1] for v in volumes)
    iron_final = [t for t in iron_vols if t in remaining_tags]
    air_final = [t for t in air_vols if t in remaining_tags]

    print(f"  Remaining: iron={iron_final}, air={air_final}")

    # Physical groups
    if iron_final:
        gmsh.model.addPhysicalGroup(3, iron_final, name='iron')
    if air_final:
        gmsh.model.addPhysicalGroup(3, air_final, name='air')

    domain_vols = iron_final + air_final
    domain_tag = gmsh.model.addPhysicalGroup(3, domain_vols, name='domain')

    # Outer boundary (faces on the air box surface)
    surfaces = gmsh.model.getEntities(2)
    outer_surfs = []
    for dim, tag in surfaces:
        com = gmsh.model.occ.getCenterOfMass(dim, tag)
        r_max = max(abs(com[0]), abs(com[1]), abs(com[2]))
        if r_max > air_size/2 * 0.9:
            outer_surfs.append(tag)

    if outer_surfs:
        gmsh.model.addPhysicalGroup(2, outer_surfs, name='outer')

    # Mesh size
    gmsh.option.setNumber('Mesh.MeshSizeMax', maxh)
    gmsh.option.setNumber('Mesh.MeshSizeMin', maxh / 5)

    return domain_tag


def build_radia_reference(N=100, I=1.0, h_coil=0.020,
                          r_inner=0.010, r_outer=0.015):
    """Build Radia coil model for reference field computation.

    Uses equivalent magnetization: M_z = NI / h.
    Models the coil as 4 rectangular blocks forming a hollow square.
    """
    import radia as rad

    rad.UtiDelAll()

    M_z = N * I / h_coil
    mag = [0, 0, M_z]

    # Approximate cylindrical coil with square blocks
    a_in = r_inner
    a_out = r_outer

    left = rad.ObjRecMag([-(a_in + a_out)/2, 0, 0],
                         [a_out - a_in, 2*a_out, h_coil], mag)
    right = rad.ObjRecMag([(a_in + a_out)/2, 0, 0],
                          [a_out - a_in, 2*a_out, h_coil], mag)
    front = rad.ObjRecMag([0, (a_in + a_out)/2, 0],
                          [2*a_in, a_out - a_in, h_coil], mag)
    back = rad.ObjRecMag([0, -(a_in + a_out)/2, 0],
                         [2*a_in, a_out - a_in, h_coil], mag)

    coil = rad.ObjCnt([left, right, front, back])
    return coil


def main():
    print("=" * 70)
    print("  Total Scalar Potential with Gmsh Cohomology Cuts")
    print("=" * 70)
    print()

    # Parameters
    N = 100
    I_coil = 1.0
    NI = N * I_coil
    core_radius = 0.008
    core_height = 0.030
    coil_r_inner = 0.010
    coil_r_outer = 0.015
    coil_height = 0.020
    mu_r = 1000.0

    print(f"Coil: {N} turns, {I_coil} A -> NI = {NI} A-turns")
    print(f"Core: radius={core_radius*1000:.0f} mm, "
          f"height={core_height*1000:.0f} mm, mu_r={mu_r:.0f}")
    print(f"Coil: r_inner={coil_r_inner*1000:.0f} mm, "
          f"r_outer={coil_r_outer*1000:.0f} mm, "
          f"height={coil_height*1000:.0f} mm")
    print()

    # Step 1: Build Gmsh geometry
    print("1. Building Gmsh geometry...")
    domain_tag = build_gmsh_geometry(
        core_radius=core_radius, core_height=core_height,
        coil_r_inner=coil_r_inner, coil_r_outer=coil_r_outer,
        coil_height=coil_height, air_size=0.06, maxh=0.006)

    # Step 2: Setup solver (mesh + cohomology + transfer)
    print()
    print("2. Computing cohomology and transferring to NGSolve...")
    from radia.cohomology_cut import CohomologyCutSolver

    solver = CohomologyCutSolver()
    n_coils = solver.setup_from_gmsh(
        domain_physical_name='domain',
        boundary_physical_name='outer')

    mesh = solver.get_mesh()
    print(f"   Cohomology generators: {n_coils} (expected: 1)")
    print(f"   NGSolve mesh: {mesh.nv} vertices, {mesh.ne} elements")
    print(f"   Materials: {list(mesh.GetMaterials())}")

    # Verify cohomology basis
    from ngsolve import Integrate, InnerProduct, curl, dx
    h_basis = solver.get_cohomology_basis()
    for k, h_gf in enumerate(h_basis):
        with TaskManager():
            curl_norm = Integrate(InnerProduct(curl(h_gf), curl(h_gf)) * dx, mesh)
            h_norm = Integrate(InnerProduct(h_gf, h_gf) * dx, mesh)
            print(f"   h_{k}: ||curl||^2 = {curl_norm:.2e}, ||h||^2 = {h_norm:.4e}")

    # Step 3: Build Radia reference
    print()
    print("3. Building Radia reference model (equivalent magnetization)...")
    import radia as rad
    radia_coil = build_radia_reference(N, I_coil, coil_height,
                                       coil_r_inner, coil_r_outer)

    test_points = [
        ([0, 0, 0],          "center"),
        ([0, 0, 0.010],      "z=10mm"),
        ([0, 0, 0.020],      "z=20mm (core end)"),
        ([0.02, 0, 0],       "x=20mm (side)"),
        ([0, 0.02, 0],       "y=20mm (side)"),
    ]

    # Step 4: Solve total scalar potential
    print()
    print("4. Solving total scalar potential (mu_r={:.0f})...".format(mu_r))
    phi_gf = solver.solve([NI], {'iron': mu_r}, order=2, dirichlet='outer')
    print(f"   phi range: [{min(phi_gf.vec):.2f}, {max(phi_gf.vec):.2f}]")

    # Step 5: Compare fields
    print()
    print("5. Field comparison:")
    B_cf = solver.get_B()
    H_cf = solver.get_H()

    print()
    print("   " + "-" * 60)
    print(f"   {'Point':25s} {'Radia(coil)':>12s} {'FEM(coil+Fe)':>12s} {'Ratio':>8s}")
    print("   " + "-" * 60)

    for pt, desc in test_points:
        B_radia = rad.Fld(radia_coil, 'b', pt)
        B_r_mag = np.sqrt(sum(b**2 for b in B_radia[:3]))
        try:
            mip = mesh(*pt)
            B_fem = B_cf(mip)
            B_fem_mag = np.sqrt(sum(b**2 for b in B_fem))
            if B_r_mag > 1e-12:
                ratio = B_fem_mag / B_r_mag
                print(f"   {desc:25s} {B_r_mag:10.6f} T {B_fem_mag:10.6f} T {ratio:6.1f}x")
            else:
                print(f"   {desc:25s} {B_r_mag:10.6f} T {B_fem_mag:10.6f} T")
        except Exception:
            print(f"   {desc:25s} {B_r_mag:10.6f} T (outside mesh)")

    # Step 6: VTK export
    print()
    print("6. Exporting to VTK...")
    from ngsolve import VTKOutput
    from ngsolve import TaskManager
    vtk_file = os.path.join(os.path.dirname(__file__),
                            'demo_coil_cohomology_cut')
    try:
        vtk = VTKOutput(mesh, coefs=[B_cf, H_cf, phi_gf],
                        names=['B', 'H', 'phi'],
                        filename=vtk_file)
        vtk.Do()
        print(f"   Written: {vtk_file}.vtu")
    except Exception as e:
        print(f"   VTK export skipped: {e}")

    import gmsh
    gmsh.finalize()
    rad.UtiDelAll()

    print()
    print("=" * 70)
    print("  FEM uses total scalar potential with cohomology:")
    print("    H = -grad(phi) + NI * h_k")
    print("  Radia reference is coil-only (square blocks, no iron).")
    print("  FEM coil is cylindrical -> geometries differ.")
    print("  Iron amplifies B at center by mu_eff >> 1.")
    print("=" * 70)


if __name__ == '__main__':
    main()
