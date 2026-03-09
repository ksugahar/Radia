#!/usr/bin/env python
"""
Test: 2-way comparison of hysteresis solvers.

Compares energy-based hysteresis (B-input Play model) across:
  1. ScalarPotentialSolver.solve_hysteresis()  -- Simkin reduced potential
  2. VectorPotentialSolver.solve_hysteresis()   -- A_r formulation

Both use the same mesh, coil source, Hantila polarization method, and
under-relaxation for convergence stability.

DC magnetization: NI=20000 AT, verifies Bz agreement within 10%.
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
# Geometry constants (simple C-type, same as demo_hysteresis_simkin.py)
# =========================================================================
IRON_PARTS = [
    ('main_leg',      0.0, 0.13125, 0.0,    0.050, 0.0625, 0.210),
    ('yoke_back_lo',  0.0, 0.060,  -0.080,  0.050, 0.080,  0.050),
    ('yoke_back_hi',  0.0, 0.060,   0.080,  0.050, 0.080,  0.050),
    ('pole_bottom',   0.0, -0.000, -0.055,   0.050, 0.040,  0.100),
    ('pole_top',      0.0, -0.000,  0.055,   0.050, 0.040,  0.100),
]

SPHERE_CENTER = [0.0, 0.070, 0.0]
AIR_R = 0.300
KELVIN_R = 0.600

NI = 20000.0


def make_hysteresis_factory():
    """Return a factory for energy-based hysteresis materials."""
    import radia as rad

    K = 10
    As = np.full(K, 8000.0)
    Js = np.full(K, 0.18)
    chi = np.linspace(0, 200, K)
    eps = 1e-8

    tables = []
    for k in range(K):
        r = np.linspace(0, Js[k], 200)
        f = As[k] * np.tan(np.pi / 2 * r / Js[k])
        f[-1] = f[-2] * 2
        tables.append((r, f))

    def factory():
        return rad.MatEnergyHysteresis(K, chi, tables, eps)

    return factory


def build_geometry(gb, maxh_gap=0.008, maxh_iron=0.020, maxh_air=0.090):
    """Build C-type geometry with Kelvin shell."""
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


def build_radia_coil(ni_val):
    """Build Radia racetrack coil for H_s source field."""
    coil_dir = os.path.join(repo_root, 'examples',
                            'c_type_electromagnet', 'mu=1000')
    if coil_dir not in sys.path:
        sys.path.insert(0, coil_dir)
    from coil_model import create_racetrack_coil
    return create_racetrack_coil(ni_val)


def main():
    print("=" * 70)
    print("  Hysteresis Solver Comparison")
    print("  Energy-based B-input Play Model, NI={:.0f} AT".format(NI))
    print("=" * 70)
    print()

    import radia as rad
    from radia.gmsh_builder import GmshBuilder
    from radia.scalar_potential_solver import ScalarPotentialSolver
    from radia.vector_potential_solver import VectorPotentialSolver

    # =====================================================================
    # Step 1: Build geometry + mesh (shared by both methods)
    # =====================================================================
    print("1. Building geometry and mesh...")
    with GmshBuilder(model_name='hys_comparison', verbose=False) as gb:
        build_geometry(gb)
        gb.generate(element_type='tet')
        mesh = gb.to_ngsolve_volume()

    fem_order = 2
    mesh.Curve(fem_order)
    print(f"   Mesh: {mesh.nv} vertices, {mesh.ne} elements")
    print(f"   Materials: {sorted(set(mesh.GetMaterials()))}")
    print()

    # =====================================================================
    # Step 2: Create Radia coil (shared source)
    # =====================================================================
    print("2. Creating Radia coil (NI={:.0f})...".format(NI))
    rad.UtiDelAll()
    coil = build_radia_coil(NI)
    print()

    mat_factory = make_hysteresis_factory()
    mip_origin = mesh(0, 0, 0)
    results = {}

    hantila_kwargs = dict(
        tol=1e-3, maxiter=60, alpha=500.0, dirichlet='outer',
        verbose=True, relax=0.5
    )

    # =====================================================================
    # Method 1: ScalarPotentialSolver (Simkin reduced potential)
    # =====================================================================
    print("3. ScalarPotentialSolver (Simkin)...")
    solver_S = ScalarPotentialSolver(
        mesh, iron_domains='iron', order=fem_order,
        kelvin_region='kelvin', kelvin_radius=AIR_R,
        kelvin_center=SPHERE_CENTER)
    solver_S.set_source_from_radia(coil, resolution=41)

    t0 = time.time()
    solver_S.solve_hysteresis(mat_factory, **hantila_kwargs)
    t_S = time.time() - t0

    B_S = solver_S.get_B()
    Bz_S = float(B_S(mip_origin)[2])
    results['Simkin'] = Bz_S
    print(f"   Bz = {Bz_S*1e3:.2f} mT, time = {t_S:.1f}s\n")

    # =====================================================================
    # Method 2: VectorPotentialSolver (A_r)
    # =====================================================================
    print("4. VectorPotentialSolver (A_r)...")
    solver_A = VectorPotentialSolver(
        mesh, iron_domains='iron', mu_r=1000.0, order=fem_order,
        kelvin_region='kelvin', kelvin_radius=AIR_R,
        kelvin_center=SPHERE_CENTER)
    solver_A.set_source_from_radia(coil, resolution=41)

    t0 = time.time()
    solver_A.solve_hysteresis(mat_factory, solver='auto', **hantila_kwargs)
    t_A = time.time() - t0

    B_A = solver_A.get_B()
    Bz_A = float(B_A(mip_origin)[2])
    results['A_r'] = Bz_A
    print(f"   Bz = {Bz_A*1e3:.2f} mT, time = {t_A:.1f}s\n")

    rad.UtiDelAll()

    # =====================================================================
    # Summary
    # =====================================================================
    print("=" * 70)
    print("  RESULTS")
    print("=" * 70)
    print(f"  {'Method':<30} {'Bz (mT)':>12} {'Time (s)':>10}")
    print(f"  {'-'*30} {'-'*12} {'-'*10}")
    print(f"  {'Simkin (H1, reduced)':<30} {Bz_S*1e3:>12.2f} {t_S:>10.1f}")
    print(f"  {'A_r (HCurl)':<30} {Bz_A*1e3:>12.2f} {t_A:>10.1f}")

    ref = abs(Bz_S)
    if ref > 1e-10:
        diff_pct = abs(Bz_S - Bz_A) / ref * 100
        print(f"\n  Simkin vs A_r difference: {diff_pct:.1f}%")

        if diff_pct < 10.0:
            print("  PASS: Methods agree within 10%")
        else:
            print("  FAIL: Methods differ by more than 10%")

    print("=" * 70)
    return results


if __name__ == '__main__':
    main()
