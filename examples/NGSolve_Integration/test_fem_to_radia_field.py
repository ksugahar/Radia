#!/usr/bin/env python
"""
Test: FEM hysteresis solve -> Radia analytical field evaluation.

Pipeline:
  1. VectorPotentialSolver.solve_hysteresis() -> M per element
  2. solver.to_radia(coil) -> Radia objects with solved M
  3. rad.Fld(combined, 'b', points) -> exact analytical B in gap
  4. Compare vs FEM direct B = B_s + curl(A_r)

The Radia analytical evaluation uses surface charge formulas that are
exact for constant M per element, with no mesh needed in the gap.
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
# Geometry constants (C-type electromagnet)
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
    print("  FEM -> Radia Analytical Field Evaluation")
    print("  Energy-based Hysteresis, NI={:.0f} AT".format(NI))
    print("=" * 70)
    print()

    import radia as rad
    from radia.gmsh_builder import GmshBuilder
    from radia.vector_potential_solver import VectorPotentialSolver

    # =================================================================
    # Step 1: Build geometry + mesh
    # =================================================================
    print("1. Building geometry and mesh...")
    with GmshBuilder(model_name='fem_to_radia', verbose=False) as gb:
        build_geometry(gb)
        gb.generate(element_type='tet')
        mesh = gb.to_ngsolve_volume()

    fem_order = 2
    mesh.Curve(fem_order)
    print(f"   Mesh: {mesh.nv} vertices, {mesh.ne} elements")
    print()

    # =================================================================
    # Step 2: Create Radia coil
    # =================================================================
    print("2. Creating Radia coil (NI={:.0f})...".format(NI))
    rad.UtiDelAll()
    coil = build_radia_coil(NI)
    print()

    # =================================================================
    # Step 3: Solve hysteresis with A_r formulation
    # =================================================================
    print("3. Solving hysteresis (A_r formulation)...")
    solver = VectorPotentialSolver(
        mesh, iron_domains='iron', mu_r=1000.0, order=fem_order,
        kelvin_region='kelvin', kelvin_radius=AIR_R,
        kelvin_center=SPHERE_CENTER)
    solver.set_source_from_radia(coil, resolution=41)

    mat_factory = make_hysteresis_factory()
    t0 = time.time()
    solver.solve_hysteresis(mat_factory, tol=1e-3, maxiter=60,
                            alpha=500.0, dirichlet='outer',
                            verbose=True, relax=0.5)
    t_solve = time.time() - t0
    print(f"   Solve time: {t_solve:.1f}s")

    # =================================================================
    # Step 4: Extract M and convert to Radia
    # =================================================================
    print()
    print("4. Converting solved M to Radia objects...")
    M_dict = solver.get_M_per_element()
    print(f"   Iron elements with M: {len(M_dict)}")

    M_mags = [np.linalg.norm(m) for m in M_dict.values()]
    print(f"   |M| range: {min(M_mags):.0f} - {max(M_mags):.0f} A/m")
    print(f"   |M| mean:  {np.mean(M_mags):.0f} A/m")

    t0 = time.time()
    combined = solver.to_radia(coil=coil)
    t_convert = time.time() - t0
    print(f"   Conversion time: {t_convert:.1f}s")
    print()

    # =================================================================
    # Step 5: Compare B field at gap center
    # =================================================================
    print("5. Comparing B at gap center (0, 0, 0)...")
    B_cf = solver.get_B()
    mip_origin = mesh(0, 0, 0)
    Bz_fem = float(B_cf(mip_origin)[2])

    B_radia = rad.Fld(combined, 'b', [0, 0, 0])
    Bz_radia = float(B_radia[2])

    print(f"   FEM direct (B_s + curl A_r): Bz = {Bz_fem*1e3:.2f} mT")
    print(f"   Radia analytical:            Bz = {Bz_radia*1e3:.2f} mT")
    diff_pct = abs(Bz_fem - Bz_radia) / max(abs(Bz_fem), 1e-10) * 100
    print(f"   Difference: {diff_pct:.1f}%")
    print()

    # =================================================================
    # Step 6: Scan along gap centerline (y-axis)
    # =================================================================
    print("6. Scanning B along gap centerline (y-axis)...")
    y_points = np.linspace(-0.020, 0.020, 41)
    scan_points = np.column_stack([
        np.zeros_like(y_points),
        y_points,
        np.zeros_like(y_points)
    ])

    # Radia batch evaluation
    t0 = time.time()
    B_radia_scan = np.array(rad.Fld(combined, 'b', scan_points))
    t_radia = time.time() - t0

    # FEM evaluation
    t0 = time.time()
    B_fem_scan = np.zeros((len(y_points), 3))
    for i, y in enumerate(y_points):
        try:
            mip = mesh(0, y, 0)
            B_fem_scan[i] = [float(B_cf(mip)[d]) for d in range(3)]
        except Exception:
            B_fem_scan[i] = [np.nan, np.nan, np.nan]
    t_fem = time.time() - t0

    print(f"   Radia batch ({len(y_points)} pts): {t_radia*1e3:.1f} ms")
    print(f"   FEM per-point ({len(y_points)} pts): {t_fem*1e3:.1f} ms")
    print()

    # Compare Bz along scan
    Bz_radia_line = B_radia_scan[:, 2]
    Bz_fem_line = B_fem_scan[:, 2]
    valid = ~np.isnan(Bz_fem_line)
    if np.any(valid):
        max_diff = np.max(np.abs(Bz_radia_line[valid] - Bz_fem_line[valid]))
        rms_diff = np.sqrt(np.mean(
            (Bz_radia_line[valid] - Bz_fem_line[valid])**2))
        print(f"   Gap Bz scan: max diff = {max_diff*1e3:.2f} mT, "
              f"RMS diff = {rms_diff*1e3:.2f} mT")

    # =================================================================
    # Step 7: Save scan data
    # =================================================================
    output_file = os.path.join(script_dir, 'fem_vs_radia_field_scan.csv')
    header = 'y_m,Bz_fem_T,Bz_radia_T'
    data = np.column_stack([y_points, Bz_fem_line, Bz_radia_line])
    np.savetxt(output_file, data, delimiter=',', header=header, comments='')
    print(f"   Saved scan data to {os.path.basename(output_file)}")
    print()

    # =================================================================
    # Summary
    # =================================================================
    print("=" * 70)
    print("  RESULTS")
    print("=" * 70)
    print(f"  {'Metric':<35} {'Value':>15}")
    print(f"  {'-'*35} {'-'*15}")
    print(f"  {'FEM Bz at origin (mT)':<35} {Bz_fem*1e3:>15.2f}")
    print(f"  {'Radia Bz at origin (mT)':<35} {Bz_radia*1e3:>15.2f}")
    print(f"  {'Difference (%)':<35} {diff_pct:>15.1f}")
    print(f"  {'Iron elements':<35} {len(M_dict):>15d}")
    print(f"  {'FEM solve time (s)':<35} {t_solve:>15.1f}")
    print(f"  {'M->Radia conversion time (s)':<35} {t_convert:>15.1f}")
    print()

    if diff_pct < 10.0:
        print("  PASS: FEM and Radia agree within 10%")
    else:
        print("  FAIL: FEM and Radia differ by more than 10%")
    print("=" * 70)

    rad.UtiDelAll()
    return diff_pct


if __name__ == '__main__':
    main()
