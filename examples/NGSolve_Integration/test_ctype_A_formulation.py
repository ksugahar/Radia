#!/usr/bin/env python
"""
Test: A_r reduced vector potential for C-type electromagnet.

Compares VectorPotentialSolver (A_r) with ScalarPotentialSolver (Simkin).
Both use the same mesh, coil, and material -- only the formulation differs.

    Simkin: H = H_s - grad(phi),  B = mu * H
    A_r:    B = B_s + curl(A_r),  H = nu * B

Uses frustum pole tip geometry matching ELF model.

References:
    ELF/Tosca:  Bz = -995 mT at origin (NI=20000, nonlinear)
    Radia BEM:  Bz ~ -945 mT (finer mesh)
    Linear:     Bz ~ -232 mT (NI=2000, mu_r=1000)
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
# Geometry constants (same as test_elf_yoke_simkin.py)
# =========================================================================
YOKE_X_MIN = -0.025
YOKE_X_MAX = 0.025
YOKE_X_EXTENT = YOKE_X_MAX - YOKE_X_MIN

SPHERE_CENTER = [0.0, 0.070, 0.0]
AIR_R = 0.300

YOKE_POLYGON_MAIN = [
    (-0.020, -0.105),
    (0.1625, -0.105),
    (0.1625, 0.105),
    (-0.020, 0.105),
    (-0.020, 0.013),
    (0.020, 0.013),
    (0.020, 0.055),
    (0.100, 0.055),
    (0.100, -0.055),
    (0.020, -0.055),
    (0.020, -0.013),
    (-0.020, -0.013),
]


def load_bh_curve():
    """Load B-H curve from BH.txt. Returns list of [H, B] pairs."""
    bh_file = os.path.join(repo_root, 'examples',
                           'c_type_electromagnet', 'nonlinear', 'BH.txt')
    data = np.loadtxt(bh_file)
    return data.tolist()


def build_radia_coil(current_at):
    """Build Radia racetrack coil."""
    coil_dir = os.path.join(repo_root, 'examples',
                            'c_type_electromagnet', 'nonlinear')
    if coil_dir not in sys.path:
        sys.path.insert(0, coil_dir)
    from coil_model import create_racetrack_coil
    return create_racetrack_coil(current_at)


def _extrude_polygon(gb, polygon_yz):
    """Extrude a Y-Z polygon in X direction."""
    n = len(polygon_yz)
    point_tags = []
    for y, z in polygon_yz:
        pt = gb.add_point(YOKE_X_MIN, y, z)
        point_tags.append(pt)
    line_tags = []
    for i in range(n):
        j = (i + 1) % n
        line_tags.append(gb.add_line(point_tags[i], point_tags[j]))
    wire = gb.add_wire(line_tags)
    surf = gb.add_plane_surface([wire])
    extruded = gb.extrude_surface(surf, YOKE_X_EXTENT, 0, 0)
    return extruded[0]


def _make_rect_wire(gb, x_min, x_max, y_min, y_max, z):
    """Create a rectangular wire loop at given Z position."""
    p1 = gb.add_point(x_min, y_min, z)
    p2 = gb.add_point(x_max, y_min, z)
    p3 = gb.add_point(x_max, y_max, z)
    p4 = gb.add_point(x_min, y_max, z)
    l1 = gb.add_line(p1, p2)
    l2 = gb.add_line(p2, p3)
    l3 = gb.add_line(p3, p4)
    l4 = gb.add_line(p4, p1)
    return gb.add_wire([l1, l2, l3, l4])


def build_geometry_and_mesh(maxh_gap=0.004, maxh_iron=0.010, maxh_air=0.060):
    """Build C-type electromagnet mesh with frustum poles (no Kelvin)."""
    from radia.gmsh_builder import GmshBuilder

    with GmshBuilder(model_name='a_form_test', verbose=False) as gb:
        # Iron yoke with frustum pole tips
        main_vol = _extrude_polygon(gb, YOKE_POLYGON_MAIN)

        w_base_bot = _make_rect_wire(gb, -0.025, 0.025, -0.020, 0.020, -0.013)
        w_tip_bot = _make_rect_wire(gb, -0.017, 0.017, -0.012, 0.012, -0.005)
        bot_frustum = gb.add_thru_sections([w_base_bot, w_tip_bot],
                                           make_solid=True)

        w_tip_top = _make_rect_wire(gb, -0.017, 0.017, -0.012, 0.012, 0.005)
        w_base_top = _make_rect_wire(gb, -0.025, 0.025, -0.020, 0.020, 0.013)
        top_frustum = gb.add_thru_sections([w_tip_top, w_base_top],
                                           make_solid=True)

        iron_vol = gb.fuse([main_vol, bot_frustum, top_frustum])

        # Air sphere (large, no Kelvin for A-formulation initial test)
        air = gb.add_sphere(SPHERE_CENTER, AIR_R)
        frag_map = gb.fragment_tracked([iron_vol, air])

        iron_ids = set(frag_map[iron_vol])
        air_ids = set(frag_map[air])

        iron_final = sorted(iron_ids)
        air_final = sorted(air_ids - iron_ids)

        gb.add_block_by_name(iron_final, 'iron')
        gb.add_block_by_name(air_final, 'air')

        air_surfs = set()
        for vid in air_final:
            air_surfs.update(gb.get_surfaces(vid))
        iron_surfs = set()
        for vid in iron_final:
            iron_surfs.update(gb.get_surfaces(vid))
        outer_surfs = sorted(air_surfs - iron_surfs)
        if outer_surfs:
            gb.add_sideset_by_name(outer_surfs, 'outer')

        print(f"  Iron: {len(iron_final)}, Air: {len(air_final)}")
        print(f"  Outer boundary: {len(outer_surfs)} surfaces")

        # Mesh size
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

        gb.generate(element_type='tet')
        mesh = gb.to_ngsolve_volume()

    return mesh


def run_comparison(ni=2000.0, nonlinear=False, order=2, resolution=41):
    """Run both A_r and Simkin on the same mesh and compare."""
    import radia as rad
    from radia.vector_potential_solver import VectorPotentialSolver
    from radia.scalar_potential_solver import ScalarPotentialSolver

    mat_label = "Nonlinear" if nonlinear else f"Linear mu_r=1000"
    print("=" * 70)
    print(f"  C-type electromagnet: {mat_label}, NI={ni:.0f}")
    print(f"  Comparing A_r (VectorPotential) vs phi_r (Simkin)")
    print("=" * 70)

    bh_data = load_bh_curve()

    # Build mesh
    print("\n1. Building geometry (frustum poles)...")
    mesh = build_geometry_and_mesh(maxh_gap=0.006, maxh_iron=0.015,
                                   maxh_air=0.060)
    mesh.Curve(order)
    print(f"   Mesh: {mesh.nv} vertices, {mesh.ne} elements")
    print(f"   Materials: {sorted(set(mesh.GetMaterials()))}")

    # Create Radia coil
    print("\n2. Creating Radia coil...")
    rad.UtiDelAll()
    coil = build_radia_coil(ni)

    mip_origin = mesh(0, 0, 0)
    results = {}

    # --- A_r formulation ---
    print("\n3. VectorPotentialSolver (A_r)...")
    solver_A = VectorPotentialSolver(
        mesh, iron_domains='iron', mu_r=1000.0, order=order)
    solver_A.set_source_from_radia(coil, resolution=resolution)

    if nonlinear:
        print("   Newton solve (coenergy)...")
        t0 = time.time()
        solver_A.solve_nonlinear_newton(
            bh_data=bh_data, tol=1e-4, maxiter=80,
            dirichlet='outer')
        t_A = time.time() - t0
    else:
        print("   Linear solve...")
        t0 = time.time()
        solver_A.solve_linear(dirichlet='outer')
        t_A = time.time() - t0

    B_A = solver_A.get_B()
    B_A_origin = B_A(mip_origin)
    Bz_A = B_A_origin[2]
    print(f"   Bz = {Bz_A*1e3:.2f} mT, time = {t_A:.1f}s")
    results['A_r'] = Bz_A

    # --- Simkin formulation ---
    print("\n4. ScalarPotentialSolver (Simkin phi_r)...")
    solver_S = ScalarPotentialSolver(
        mesh, iron_domains='iron', mu_r=1000.0, order=order)
    solver_S.set_source_from_radia(coil, resolution=resolution)

    if nonlinear:
        print("   Newton solve (energy)...")
        t0 = time.time()
        solver_S.solve_nonlinear_newton(
            bh_data=bh_data, tol=1e-4, maxiter=80,
            dirichlet='outer')
        t_S = time.time() - t0
    else:
        print("   Linear solve...")
        t0 = time.time()
        solver_S.solve_single_potential(dirichlet='outer')
        t_S = time.time() - t0

    B_S = solver_S.get_B()
    B_S_origin = B_S(mip_origin)
    Bz_S = B_S_origin[2]
    print(f"   Bz = {Bz_S*1e3:.2f} mT, time = {t_S:.1f}s")
    results['Simkin'] = Bz_S

    rad.UtiDelAll()

    # Summary
    print("\n" + "=" * 70)
    print("  RESULTS")
    print("=" * 70)
    print(f"  {'Method':<20} {'Bz (mT)':>12} {'Time (s)':>10}")
    print(f"  {'-'*20} {'-'*12} {'-'*10}")
    print(f"  {'A_r (HCurl)':<20} {Bz_A*1e3:>12.2f} {t_A:>10.1f}")
    print(f"  {'Simkin (H1)':<20} {Bz_S*1e3:>12.2f} {t_S:>10.1f}")

    diff_pct = abs(Bz_A - Bz_S) / abs(Bz_S) * 100 if abs(Bz_S) > 1e-10 else 0
    print(f"\n  A_r vs Simkin difference: {diff_pct:.1f}%")

    if nonlinear and ni == 20000:
        Bz_ref = -0.945  # BEM reference
        err_A = abs(Bz_A * 1e3 - Bz_ref * 1e3) / abs(Bz_ref * 1e3) * 100
        err_S = abs(Bz_S * 1e3 - Bz_ref * 1e3) / abs(Bz_ref * 1e3) * 100
        print(f"  vs BEM ~{Bz_ref*1e3:.0f} mT: A_r={err_A:.1f}%, "
              f"Simkin={err_S:.1f}%")
    elif not nonlinear and ni == 2000:
        Bz_ref = -0.232
        err_A = abs(Bz_A * 1e3 - Bz_ref * 1e3) / abs(Bz_ref * 1e3) * 100
        err_S = abs(Bz_S * 1e3 - Bz_ref * 1e3) / abs(Bz_ref * 1e3) * 100
        print(f"  vs BEM ~{Bz_ref*1e3:.0f} mT: A_r={err_A:.1f}%, "
              f"Simkin={err_S:.1f}%")

    print("=" * 70)
    return results


def main():
    print("C-type Electromagnet: A_r vs Simkin comparison")
    print("Frustum pole tips, no Kelvin (large air sphere)\n")

    # Linear test first (fast)
    print("\n" + "#" * 70)
    print("# LINEAR (mu_r=1000, NI=2000)")
    print("#" * 70)
    results_lin = run_comparison(ni=2000.0, nonlinear=False,
                                 order=2, resolution=41)

    # Nonlinear test
    print("\n" + "#" * 70)
    print("# NONLINEAR (B-H curve, NI=20000)")
    print("#" * 70)
    results_nl = run_comparison(ni=20000.0, nonlinear=True,
                                order=2, resolution=41)


if __name__ == '__main__':
    main()
