#!/usr/bin/env python
"""
C-Type Electromagnet (mu_r=1000, linear) - Cohomology Total Scalar Potential

Full model using:
  - GmshBuilder for geometry construction and mesh control
  - Gmsh cohomology for coil current (no Biot-Savart needed)
  - Kelvin transform for open boundary
  - NGSolve H1 FEM solve

Geometry (all in meters):
  Main Leg:     x[-25,25] y[100,162.5] z[-105,105] mm
  Yoke Back Lo: x[-25,25] y[20,100]    z[-105,-55] mm
  Yoke Back Hi: x[-25,25] y[20,100]    z[55,105]   mm
  Pole Bottom:  x[-25,25] y[-20,20]    z[-105,-5]  mm
  Pole Top:     x[-25,25] y[-20,20]    z[5,105]    mm
  Gap:          z in [-5, 5] mm, field at origin

Coil: racetrack (2000 AT), rectangular tube hole for cohomology
Expected B at gap center: ~228 mT
"""

import sys
import os
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(script_dir))
sys.path.insert(0, os.path.join(repo_root, 'src', 'radia'))
sys.path.insert(0, os.path.join(repo_root, 'src'))

MU_0 = 4 * np.pi * 1e-7

# =========================================================================
# Geometry constants (meters)
# =========================================================================
# (name, center_x, center_y, center_z, size_x, size_y, size_z)
IRON_PARTS = [
    ('main_leg',     0.0, 0.13125, 0.0,    0.050, 0.0625, 0.210),
    ('yoke_back_lo', 0.0, 0.060,  -0.080,  0.050, 0.080,  0.050),
    ('yoke_back_hi', 0.0, 0.060,   0.080,  0.050, 0.080,  0.050),
    ('pole_bottom',  0.0, -0.000, -0.055,   0.050, 0.040,  0.100),
    ('pole_top',     0.0, -0.000,  0.055,   0.050, 0.040,  0.100),
]

# Rectangular tube coil hole (loops around main leg in x-y plane)
TUBE_OUTER_CENTER = (0.0, 0.1325, 0.0)
TUBE_OUTER_SIZE = (0.110, 0.115, 0.008)
TUBE_INNER_CENTER = (0.0, 0.1325, 0.0)
TUBE_INNER_SIZE = (0.100, 0.105, 0.008)

# Air sphere + Kelvin shell
SPHERE_CENTER = [0.0, 0.070, 0.0]
AIR_R = 0.300
KELVIN_R = 0.600

NI = 2000.0
MU_R_IRON = 1000.0


def build_ctype_geometry(gb, maxh_gap=0.005, maxh_iron=0.012, maxh_air=0.150):
    """Build C-type geometry and configure mesh sizing using GmshBuilder.

    Parameters
    ----------
    gb : GmshBuilder
        Active GmshBuilder context.
    maxh_gap : float
        Element size in gap region (m).
    maxh_iron : float
        Element size in iron region (m).
    maxh_air : float
        Element size in air/Kelvin region (m).

    Returns
    -------
    dict
        Volume classification: {name: [vol_ids]}.
    """
    # --- Iron pieces: 5 boxes fused into one ---
    iron_ids = []
    for name, cx, cy, cz, sx, sy, sz in IRON_PARTS:
        iron_ids.append(gb.add_box([cx, cy, cz], [sx, sy, sz]))
    iron = gb.fuse(iron_ids)

    # --- Coil hole: rectangular tube (outer - inner) ---
    outer = gb.add_box(list(TUBE_OUTER_CENTER), list(TUBE_OUTER_SIZE))
    inner = gb.add_box(list(TUBE_INNER_CENTER), list(TUBE_INNER_SIZE))
    tube = gb.cut(outer, inner)

    # --- Air sphere + Kelvin shell ---
    air = gb.add_sphere(SPHERE_CENTER, AIR_R)
    kelvin = gb.add_sphere(SPHERE_CENTER, KELVIN_R)

    # --- Fragment with tracking ---
    frag_map = gb.fragment_tracked([iron, tube, air, kelvin])

    iron_ids = set(frag_map[iron])
    tube_ids = set(frag_map[tube])
    air_ids = set(frag_map[air])
    kelvin_ids = set(frag_map[kelvin])

    # Classify: remove overlaps
    iron_final = sorted(iron_ids)
    tube_final = sorted(tube_ids - iron_ids)
    air_final = sorted(air_ids - iron_ids - tube_ids)
    kelvin_final = sorted(kelvin_ids - air_ids)

    print(f"  Iron: {len(iron_final)}, Air: {len(air_final)}, "
          f"Kelvin: {len(kelvin_final)}, Tube(hole): {len(tube_final)}")

    # --- Physical groups ---
    gb.add_block_by_name(iron_final, 'iron')
    gb.add_block_by_name(air_final, 'air')
    gb.add_block_by_name(kelvin_final, 'kelvin')
    gb.add_block_by_name(tube_final, 'coil_hole')

    domain_ids = iron_final + air_final + kelvin_final
    gb.add_block_by_name(domain_ids, 'domain')

    # Outer boundary: kelvin surfaces NOT shared with air
    kelvin_surfs = set()
    for vid in kelvin_final:
        kelvin_surfs.update(gb.get_surfaces(vid))
    air_surfs = set()
    for vid in air_final:
        air_surfs.update(gb.get_surfaces(vid))
    outer_surfs = sorted(kelvin_surfs - air_surfs)

    if outer_surfs:
        gb.add_sideset_by_name(outer_surfs, 'outer')
    print(f"  Outer boundary: {len(outer_surfs)} surfaces")

    # --- Mesh control: Box fields for gap and iron regions ---
    f_gap = gb.add_field_box(
        -0.030, 0.030, -0.025, 0.025, -0.010, 0.010,
        size_in=maxh_gap, size_out=maxh_air)

    f_iron = gb.add_field_box(
        -0.030, 0.030, -0.025, 0.170, -0.110, 0.110,
        size_in=maxh_iron, size_out=maxh_air)

    f_min = gb.add_field_min([f_gap, f_iron])
    gb.set_background_field(f_min)

    gb.set_min_max_size(maxh_gap / 3, maxh_air)
    gb.set_algorithm_3d(1)  # Delaunay

    return {
        'iron': iron_final, 'air': air_final,
        'kelvin': kelvin_final, 'tube': tube_final,
    }


def build_radia_reference():
    """Build Radia model (iron blocks + racetrack coil) for comparison."""
    import radia as rad

    coil_dir = os.path.join(repo_root, 'examples',
                            'c_type_electromagnet', 'mu=1000')
    sys.path.insert(0, coil_dir)
    from coil_model import create_racetrack_coil

    rad.UtiDelAll()

    blocks = []
    for name, cx, cy, cz, dx, dy, dz in IRON_PARTS:
        x0, y0, z0 = cx - dx/2, cy - dy/2, cz - dz/2
        nx, ny, nz = 2, 2, max(2, int(round(dz / dx * 2)))
        sx, sy, sz = dx / nx, dy / ny, dz / nz
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    bx = x0 + (ix + 0.5) * sx
                    by = y0 + (iy + 0.5) * sy
                    bz = z0 + (iz + 0.5) * sz
                    block = rad.ObjRecMag([bx, by, bz],
                                         [sx, sy, sz], [0, 0, 0])
                    mat = rad.MatLin(MU_R_IRON)
                    rad.MatApl(block, mat)
                    blocks.append(block)

    iron = rad.ObjCnt(blocks)
    coil = create_racetrack_coil(NI)
    model = rad.ObjCnt([iron, coil])

    print(f"   Radia model: {len(blocks)} iron blocks")
    result = rad.Solve(model, 0.0001, 1000, 0)
    print(f"   Solve result: max|M|={result[0]:.2f}, max|dM|={result[1]:.4e}")

    return model


def main():
    print("=" * 70)
    print("  C-Type Electromagnet: Cohomology Total Scalar Potential + Kelvin")
    print("=" * 70)
    print()
    print(f"NI = {NI:.0f} A-turns, mu_r = {MU_R_IRON:.0f}")
    print(f"Expected B at gap center: ~228 mT")
    print()

    from radia.gmsh_builder import GmshBuilder
    from radia.cohomology_cut import CohomologyCutSolver

    # =====================================================================
    # Step 1-2: Build geometry + cohomology (requires active Gmsh session)
    # =====================================================================
    print("1. Building geometry with GmshBuilder...", flush=True)

    with GmshBuilder(model_name='c_type_cohomology', verbose=False) as gb:
        vols = build_ctype_geometry(gb,
                                    maxh_gap=0.003,
                                    maxh_iron=0.008,
                                    maxh_air=0.060)
        print()

        print("2. Computing cohomology + transferring to NGSolve...")
        solver = CohomologyCutSolver()
        n_coils = solver.setup_from_gmsh(
            domain_physical_name='domain',
            boundary_physical_name='outer')
    # Gmsh finalized here -- NGSolve mesh is independent

    mesh = solver.get_mesh()
    mesh.Curve(3)
    print(f"   Cohomology generators: {n_coils} (expected: 1)")
    print(f"   Mesh: {mesh.nv} vertices, {mesh.ne} elements")
    print(f"   Materials: {sorted(set(mesh.GetMaterials()))}")
    print(f"   Mesh curved to order 3")

    from ngsolve import Integrate, InnerProduct, curl, dx
    for k, h_gf in enumerate(solver.get_cohomology_basis()):
        curl_norm = Integrate(InnerProduct(curl(h_gf), curl(h_gf)) * dx, mesh)
        h_norm = Integrate(InnerProduct(h_gf, h_gf) * dx, mesh)
        print(f"   h_{k}: ||curl||^2 = {curl_norm:.2e}, ||h||^2 = {h_norm:.4e}")

    if n_coils != 1:
        print(f"   WARNING: Expected 1 generator, got {n_coils}")
        return

    # =====================================================================
    # Step 3: Solve FEM with Kelvin transform
    # =====================================================================
    print()
    print("3. Solving total scalar potential with Kelvin transform...")
    phi_gf = solver.solve(
        [NI],
        mu_r_dict={'iron': MU_R_IRON},
        order=3,
        dirichlet='outer',
        kelvin_region='kelvin',
        kelvin_radius=AIR_R,
        kelvin_center=SPHERE_CENTER)

    print(f"   phi range: [{min(phi_gf.vec):.2f}, {max(phi_gf.vec):.2f}]")
    B_cf = solver.get_B()
    H_cf = solver.get_H()

    try:
        mip = mesh(0, 0, 0)
        B_fem = B_cf(mip)
        B_z = B_fem[2]
        B_mag = np.sqrt(sum(b**2 for b in B_fem))
        print(f"   B at origin: Bx={B_fem[0]*1e3:.2f}, "
              f"By={B_fem[1]*1e3:.2f}, Bz={B_fem[2]*1e3:.2f} mT")
        print(f"   |B| = {B_mag*1e3:.2f} mT")
    except Exception as e:
        print(f"   Cannot evaluate at origin: {e}")
        B_z = 0.0

    # =====================================================================
    # Step 4: Radia reference
    # =====================================================================
    print()
    print("4. Building Radia reference model...")
    try:
        import radia as rad
        radia_model = build_radia_reference()

        test_points = [
            ([0, 0, 0],          "gap center"),
            ([0, 0, 0.003],      "z=3mm"),
            ([0, 0, -0.003],     "z=-3mm"),
            ([0, 0.010, 0],      "y=10mm"),
            ([0.010, 0, 0],      "x=10mm"),
        ]

        print()
        print("   " + "-" * 65)
        print(f"   {'Point':18s} {'FEM B (mT)':>12s} {'Radia B (mT)':>14s} "
              f"{'Diff':>8s}")
        print("   " + "-" * 65)

        for pt, desc in test_points:
            B_radia = np.array(rad.Fld(radia_model, 'b', pt))
            B_r_mag = np.linalg.norm(B_radia)
            try:
                mip = mesh(*pt)
                B_fem_pt = np.array(B_cf(mip))
                B_f_mag = np.linalg.norm(B_fem_pt)
                if B_r_mag > 1e-6:
                    pct = (B_f_mag - B_r_mag) / B_r_mag * 100
                    print(f"   {desc:18s} {B_f_mag*1e3:10.2f}   "
                          f"{B_r_mag*1e3:12.2f}   {pct:+6.1f}%")
                else:
                    print(f"   {desc:18s} {B_f_mag*1e3:10.2f}   "
                          f"{B_r_mag*1e3:12.2f}")
            except Exception:
                print(f"   {desc:18s} {'(outside)':>10s}   "
                      f"{B_r_mag*1e3:12.2f}")

        rad.UtiDelAll()

    except ImportError:
        print("   Radia not available, skipping comparison")

    # =====================================================================
    # Step 5: VTK export
    # =====================================================================
    print()
    print("5. Exporting to VTK...")
    from ngsolve import VTKOutput
    vtk_file = os.path.join(script_dir, 'demo_ctype_cohomology')
    try:
        vtk = VTKOutput(mesh, coefs=[B_cf, H_cf, phi_gf],
                        names=['B', 'H', 'phi'],
                        filename=vtk_file)
        vtk.Do()
        print(f"   Written: {vtk_file}.vtu")
    except Exception as e:
        print(f"   VTK export skipped: {e}")

    # =====================================================================
    # Summary
    # =====================================================================
    print()
    print("=" * 70)
    print("  Summary")
    print("=" * 70)
    print(f"  FEM mesh: {mesh.nv} vertices, {mesh.ne} elements")
    print(f"  Cohomology generators: {n_coils}")
    print(f"  FEM B at gap center: {B_z*1e3:.2f} mT (expected ~228 mT)")
    print()
    print("  Method: H = -grad(phi) + NI * h_k")
    print("  Kelvin transform on outer shell (R = {:.0f} mm)".format(AIR_R*1e3))
    print("  No Biot-Savart computation needed")
    print("=" * 70)


if __name__ == '__main__':
    main()
