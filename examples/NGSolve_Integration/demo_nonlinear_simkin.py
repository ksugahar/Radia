#!/usr/bin/env python
"""
C-Type Electromagnet - Nonlinear Simkin-Trowbridge Solver

Same geometry as demo_ctype_simkin.py but with nonlinear B-H curve
instead of constant mu_r=1000.

Uses ScalarPotentialSolver.solve_nonlinear() which performs Picard
iteration: solve linear -> evaluate H -> update mu from B-H -> repeat.
"""

import sys
import os
import time
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(script_dir))
sys.path.insert(0, os.path.join(repo_root, 'src', 'radia'))
sys.path.insert(0, os.path.join(repo_root, 'src'))

MU_0 = 4 * np.pi * 1e-7

# =========================================================================
# Geometry constants (meters) -- same as other C-type demos
# =========================================================================
IRON_PARTS = [
    ('main_leg',     0.0, 0.13125, 0.0,    0.050, 0.0625, 0.210),
    ('yoke_back_lo', 0.0, 0.060,  -0.080,  0.050, 0.080,  0.050),
    ('yoke_back_hi', 0.0, 0.060,   0.080,  0.050, 0.080,  0.050),
    ('pole_bottom',  0.0, -0.000, -0.055,   0.050, 0.040,  0.100),
    ('pole_top',     0.0, -0.000,  0.055,   0.050, 0.040,  0.100),
]

SPHERE_CENTER = [0.0, 0.070, 0.0]
AIR_R = 0.300
KELVIN_R = 0.600

NI = 20000.0


def load_bh_curve():
    """Load B-H curve from BH.txt (from ELF nonlinear example)."""
    bh_file = os.path.join(repo_root, 'examples',
                           'c_type_electromagnet', 'nonlinear', 'BH.txt')
    data = np.loadtxt(bh_file)
    return data.tolist()  # [[H, B], ...]


def build_simkin_geometry(gb, maxh_gap=0.003, maxh_iron=0.008, maxh_air=0.060):
    """Build C-type geometry WITHOUT coil hole (Simkin method)."""
    iron_ids = []
    for name, cx, cy, cz, sx, sy, sz in IRON_PARTS:
        iron_ids.append(gb.add_box([cx, cy, cz], [sx, sy, sz]))
    iron = gb.fuse(iron_ids)

    air = gb.add_sphere(SPHERE_CENTER, AIR_R)
    kelvin = gb.add_sphere(SPHERE_CENTER, KELVIN_R)

    frag_map = gb.fragment_tracked([iron, air, kelvin])

    iron_ids = set(frag_map[iron])
    air_ids = set(frag_map[air])
    kelvin_ids = set(frag_map[kelvin])

    iron_final = sorted(iron_ids)
    air_final = sorted(air_ids - iron_ids)
    kelvin_final = sorted(kelvin_ids - air_ids)

    print(f"  Iron: {len(iron_final)}, Air: {len(air_final)}, "
          f"Kelvin: {len(kelvin_final)}")

    gb.add_block_by_name(iron_final, 'iron')
    gb.add_block_by_name(air_final, 'air')
    gb.add_block_by_name(kelvin_final, 'kelvin')

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

    f_gap = gb.add_field_box(
        -0.030, 0.030, -0.025, 0.025, -0.010, 0.010,
        size_in=maxh_gap, size_out=maxh_air)
    f_iron = gb.add_field_box(
        -0.030, 0.030, -0.025, 0.170, -0.110, 0.110,
        size_in=maxh_iron, size_out=maxh_air)
    f_min = gb.add_field_min([f_gap, f_iron])
    gb.set_background_field(f_min)
    gb.set_min_max_size(maxh_gap / 3, maxh_air)
    gb.set_algorithm_3d(1)


def build_radia_coil():
    """Build Radia racetrack coil for H_s source field."""
    coil_dir = os.path.join(repo_root, 'examples',
                            'c_type_electromagnet', 'mu=1000')
    if coil_dir not in sys.path:
        sys.path.insert(0, coil_dir)
    from coil_model import create_racetrack_coil
    return create_racetrack_coil(NI)


def build_radia_nonlinear_reference(bh_data):
    """Build Radia nonlinear reference model."""
    import radia as rad

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
                    mat = rad.MatSatIsoTab(bh_data)
                    rad.MatApl(block, mat)
                    blocks.append(block)

    iron = rad.ObjCnt(blocks)
    coil = build_radia_coil()
    model = rad.ObjCnt([iron, coil])

    print(f"   Radia: {len(blocks)} blocks (nonlinear)")
    result = rad.Solve(model, 0.0001, 1000, 0)
    print(f"   Solve: max|M|={result[0]:.2f}, max|dM|={result[1]:.4e}")
    return model


def main():
    print("=" * 70)
    print("  C-Type Electromagnet: NONLINEAR Simkin + Kelvin")
    print("=" * 70)
    print()
    print(f"NI = {NI:.0f} A-turns")
    print(f"Material: Nonlinear B-H curve (from ELF)")
    print()

    import radia as rad
    from radia.gmsh_builder import GmshBuilder
    from radia.scalar_potential_solver import ScalarPotentialSolver

    # Load B-H curve
    bh_data = load_bh_curve()
    print(f"B-H curve: {len(bh_data)} points")
    print(f"  H range: [{bh_data[0][0]:.0f}, {bh_data[-1][0]:.0f}] A/m")
    print(f"  B range: [{bh_data[0][1]:.3f}, {bh_data[-1][1]:.3f}] T")
    print()

    # =====================================================================
    # Step 1: Build geometry + mesh
    # =====================================================================
    print("1. Building geometry...", flush=True)

    with GmshBuilder(model_name='c_type_nonlinear', verbose=False) as gb:
        build_simkin_geometry(gb, maxh_gap=0.004, maxh_iron=0.012,
                              maxh_air=0.070)
        print()
        print("2. Generating mesh...")
        gb.generate(element_type='tet')
        mesh = gb.to_ngsolve_volume()

    fem_order = 2
    mesh.Curve(fem_order)
    print(f"   Mesh: {mesh.nv} vertices, {mesh.ne} elements")
    print()

    # =====================================================================
    # Step 2: Solve nonlinear
    # =====================================================================
    print("3. Computing H_s from Radia coil...", flush=True)
    rad.UtiDelAll()
    coil = build_radia_coil()

    solver = ScalarPotentialSolver(
        mesh, iron_domains='iron', order=fem_order,
        kelvin_region='kelvin', kelvin_radius=AIR_R,
        kelvin_center=SPHERE_CENTER)
    solver.set_source_from_radia(coil, resolution=51)
    print()

    print("4. Newton + SymbolicEnergy (B-H curve)...")
    t0 = time.time()
    phi_gf = solver.solve_nonlinear_newton(
        bh_data=bh_data,
        tol=1e-4,
        maxiter=50,
        dirichlet='outer')
    t_newton = time.time() - t0
    print(f"   Newton solve time: {t_newton:.1f}s")

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
    print()

    # =====================================================================
    # Step 3: Radia nonlinear reference
    # =====================================================================
    print("5. Radia nonlinear reference...")
    radia_model = build_radia_nonlinear_reference(bh_data)

    test_points = [
        ([0, 0, 0],          "gap center"),
        ([0, 0, 0.003],      "z=3mm"),
        ([0, 0, -0.003],     "z=-3mm"),
        ([0, 0.010, 0],      "y=10mm"),
    ]

    print()
    print("   " + "-" * 60)
    print(f"   {'Point':16s} {'FEM (mT)':>10s} {'Radia (mT)':>12s} {'Diff':>8s}")
    print("   " + "-" * 60)

    for pt, desc in test_points:
        B_radia = np.array(rad.Fld(radia_model, 'b', pt))
        B_r_mag = np.linalg.norm(B_radia)
        try:
            mip = mesh(*pt)
            B_fem_pt = np.array(B_cf(mip))
            B_f_mag = np.linalg.norm(B_fem_pt)
            pct = (B_f_mag - B_r_mag) / B_r_mag * 100 if B_r_mag > 1e-6 else 0
            print(f"   {desc:16s} {B_f_mag*1e3:8.2f}   {B_r_mag*1e3:10.2f}   {pct:+6.1f}%")
        except Exception:
            print(f"   {desc:16s} {'(outside)':>8s}   {B_r_mag*1e3:10.2f}")

    rad.UtiDelAll()

    # =====================================================================
    # Step 4: VTK export
    # =====================================================================
    print()
    print("6. Exporting to VTK...")
    from ngsolve import VTKOutput
    vtk_file = os.path.join(script_dir, 'demo_nonlinear_simkin')
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
    print(f"  FEM B at gap center: {B_z*1e3:.2f} mT")
    print()
    print("  Method: H = H_s - grad(phi)  [Simkin 1979]")
    print("  Material: Nonlinear B-H curve (Newton + SymbolicEnergy)")
    print("  Open boundary: Kelvin transform (R = {:.0f} mm)".format(AIR_R*1e3))
    print("=" * 70)


if __name__ == '__main__':
    main()
