#!/usr/bin/env python
"""
Test: ELF yoke geometry with tapered pole tips + Simkin nonlinear solver.

Purpose: Match the exact ELF C-type yoke geometry including 3D tapered
pole tips (frustum shape) near the gap.

ELF pole tip geometry (from .meg file analysis):
  At Z=+-13mm (taper base): X=[-25,25], Y=[-20,20] (full rectangular)
  At Z=+-5mm (gap face):    X=[-17,17], Y=[-12,12] (tapered inward)
  Taper: 45-degree chamfer in both X and Y (8mm reduction per side)
  Gap: 10mm (Z = -5 to +5)

Geometry variants:
  'rect'  - Original rectangular poles (Y=[-20,20] at gap), extrusion
  'taper' - Y-Z taper only (Y=[-12,12] at gap), same X extrusion
  'frustum' - Full 3D frustum (X and Y taper), matches ELF exactly

References:
  ELF/Tosca 52-elem (MIMA):  Bz = -995 mT at origin (NI=20000)
  Radia BEM (finer mesh):    Bz ~ -945 mT
  Radia BEM (linear, 2000AT): Bz ~ -232 mT
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
# ELF yoke outer boundary (C-shape cross-section in Y-Z plane)
# All coordinates from ELF .meg file analysis (in meters)
# =========================================================================
# Geometry (mm -> m):
#   Back yoke:   Y = [100, 162.5], Z = [-105, 105]
#   Upper arm:   Y = [20, 100],    Z = [55, 105]
#   Lower arm:   Y = [20, 100],    Z = [-105, -55]
#   Pole (main): Y = [-20, 20],    Z = [13, 105] and [-105, -13]
#   Pole tips:   Tapered frustum from Z=+-13 to Z=+-5
#                X: [-25,25] -> [-17,17], Y: [-20,20] -> [-12,12]
#   Gap:         10mm (Z = -5 to +5)
#   Coil window: Y = [20, 100],    Z = [-55, 55] (air)

# Original rectangular polygon (Y, Z) -- poles extend to Z=+-5 at full width
YOKE_POLYGON_RECT = [
    (-0.020, -0.105),   # 1: bottom-left (pole bottom, outer)
    ( 0.1625, -0.105),  # 2: bottom-right (back yoke, bottom)
    ( 0.1625,  0.105),  # 3: top-right (back yoke, top)
    (-0.020,  0.105),   # 4: top-left (pole top, outer)
    (-0.020,  0.005),   # 5: gap top (pole top inner face)
    ( 0.020,  0.005),   # 6: gap top right
    ( 0.020,  0.055),   # 7: coil window top-left
    ( 0.100,  0.055),   # 8: coil window top-right
    ( 0.100, -0.055),   # 9: coil window bottom-right
    ( 0.020, -0.055),   # 10: coil window bottom-left
    ( 0.020, -0.005),   # 11: gap bottom right
    (-0.020, -0.005),   # 12: gap bottom (pole bottom inner face)
]

# Tapered polygon (Y, Z) -- pole tips taper from Y=+-20 to Y=+-12 at gap
YOKE_POLYGON_TAPER = [
    (-0.020, -0.105),   # 1: bottom-left (pole bottom, outer)
    ( 0.1625, -0.105),  # 2: bottom-right (back yoke, bottom)
    ( 0.1625,  0.105),  # 3: top-right (back yoke, top)
    (-0.020,  0.105),   # 4: top-left (pole top, outer)
    (-0.020,  0.013),   # 5: taper start (top, Y=-20)
    (-0.012,  0.005),   # 6: gap top (tapered to Y=-12)
    ( 0.012,  0.005),   # 7: gap top right (tapered to Y=+12)
    ( 0.020,  0.013),   # 8: taper end (top, Y=+20)
    ( 0.020,  0.055),   # 9: coil window top-left
    ( 0.100,  0.055),   # 10: coil window top-right
    ( 0.100, -0.055),   # 11: coil window bottom-right
    ( 0.020, -0.055),   # 12: coil window bottom-left
    ( 0.020, -0.013),   # 13: taper start (bottom, Y=+20)
    ( 0.012, -0.005),   # 14: gap bottom right (tapered to Y=+12)
    (-0.012, -0.005),   # 15: gap bottom (tapered to Y=-12)
    (-0.020, -0.013),   # 16: taper end (bottom, Y=-20)
]

# Main body polygon (Y, Z) -- poles end at Z=+-13, no pole tip
# Used as base for 3D frustum approach
YOKE_POLYGON_MAIN = [
    (-0.020, -0.105),   # bottom-left (pole bottom, outer)
    ( 0.1625, -0.105),  # back yoke, bottom
    ( 0.1625,  0.105),  # back yoke, top
    (-0.020,  0.105),   # pole top, outer
    (-0.020,  0.013),   # pole top, taper start (Z=13mm)
    ( 0.020,  0.013),   # inner edge at Z=13mm
    ( 0.020,  0.055),   # coil window top-left
    ( 0.100,  0.055),   # coil window top-right
    ( 0.100, -0.055),   # coil window bottom-right
    ( 0.020, -0.055),   # coil window bottom-left
    ( 0.020, -0.013),   # inner edge at Z=-13mm
    (-0.020, -0.013),   # pole bottom, taper start (Z=-13mm)
]

YOKE_X_MIN = -0.025  # -25 mm
YOKE_X_MAX =  0.025  # +25 mm
YOKE_X_EXTENT = YOKE_X_MAX - YOKE_X_MIN  # 50 mm

SPHERE_CENTER = [0.0, 0.070, 0.0]
AIR_R = 0.300
KELVIN_R = 0.600


def load_bh_curve():
    """Load B-H curve from BH.txt."""
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
    """Extrude a Y-Z polygon in X direction. Returns volume ID."""
    n = len(polygon_yz)
    point_tags = []
    for y, z in polygon_yz:
        pt = gb.add_point(YOKE_X_MIN, y, z)
        point_tags.append(pt)

    line_tags = []
    for i in range(n):
        j = (i + 1) % n
        line = gb.add_line(point_tags[i], point_tags[j])
        line_tags.append(line)

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


def build_elf_yoke_polygon(gb, geometry='frustum'):
    """Build C-type yoke matching ELF geometry.

    Parameters
    ----------
    geometry : str
        'rect'    - Rectangular poles (original, wrong)
        'taper'   - Y-Z taper only (2D extrusion)
        'frustum' - Full 3D frustum pole tips (matches ELF)
    """
    if geometry == 'rect':
        return _extrude_polygon(gb, YOKE_POLYGON_RECT)

    elif geometry == 'taper':
        return _extrude_polygon(gb, YOKE_POLYGON_TAPER)

    elif geometry == 'frustum':
        # Main body (poles end at Z=+-13mm)
        main_vol = _extrude_polygon(gb, YOKE_POLYGON_MAIN)

        # Bottom frustum: Z=-13 to Z=-5, XY tapering
        # Base at Z=-13mm: [-25,25] x [-20,20]
        # Tip at Z=-5mm:   [-17,17] x [-12,12]
        w_base_bot = _make_rect_wire(gb,
                                     -0.025, 0.025, -0.020, 0.020, -0.013)
        w_tip_bot = _make_rect_wire(gb,
                                    -0.017, 0.017, -0.012, 0.012, -0.005)
        bot_frustum = gb.add_thru_sections(
            [w_base_bot, w_tip_bot], make_solid=True)

        # Top frustum: Z=+5 to Z=+13, XY tapering
        w_tip_top = _make_rect_wire(gb,
                                    -0.017, 0.017, -0.012, 0.012, 0.005)
        w_base_top = _make_rect_wire(gb,
                                     -0.025, 0.025, -0.020, 0.020, 0.013)
        top_frustum = gb.add_thru_sections(
            [w_tip_top, w_base_top], make_solid=True)

        # Fuse all volumes
        yoke_vol = gb.fuse([main_vol, bot_frustum, top_frustum])
        return yoke_vol

    else:
        raise ValueError(f"Unknown geometry type: {geometry}")


def build_elf_yoke_geometry(gb, use_kelvin=True,
                            maxh_gap=0.004, maxh_iron=0.010,
                            maxh_air=0.060, geometry='frustum'):
    """Build C-type yoke geometry with air/kelvin domains."""
    iron_vol = build_elf_yoke_polygon(gb, geometry=geometry)

    if use_kelvin:
        air = gb.add_sphere(SPHERE_CENTER, AIR_R)
        kelvin = gb.add_sphere(SPHERE_CENTER, KELVIN_R)
        frag_map = gb.fragment_tracked([iron_vol, air, kelvin])

        iron_ids = set(frag_map[iron_vol])
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

        print(f"  Iron: {len(iron_final)}, Air: {len(air_final)}, "
              f"Kelvin: {len(kelvin_final)}")
        print(f"  Outer boundary: {len(outer_surfs)} surfaces")
    else:
        air_r = 0.600
        air = gb.add_sphere(SPHERE_CENTER, air_r)
        frag_map = gb.fragment_tracked([iron_vol, air])

        iron_ids = set(frag_map[iron_vol])
        air_ids = set(frag_map[air])

        iron_final = sorted(iron_ids)
        air_final = sorted(air_ids - iron_ids)

        gb.add_block_by_name(iron_final, 'iron')
        gb.add_block_by_name(air_final, 'air')

        air_surfs_all = set()
        for vid in air_final:
            air_surfs_all.update(gb.get_surfaces(vid))
        iron_surfs = set()
        for vid in iron_final:
            iron_surfs.update(gb.get_surfaces(vid))
        outer_surfs = sorted(air_surfs_all - iron_surfs)
        if outer_surfs:
            gb.add_sideset_by_name(outer_surfs, 'outer')

        print(f"  Iron: {len(iron_final)}, Air: {len(air_final)}")
        print(f"  Outer boundary: {len(outer_surfs)} surfaces")

    # Mesh size fields
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


def run_test(use_kelvin=True, ni=20000.0, nonlinear=True,
             maxh_gap=0.004, maxh_iron=0.010, resolution=51,
             geometry='frustum'):
    """Run Simkin solve with ELF polygon yoke geometry."""
    import radia as rad
    from radia.gmsh_builder import GmshBuilder
    from radia.scalar_potential_solver import ScalarPotentialSolver

    label = "Kelvin" if use_kelvin else "No-Kelvin"
    mat_label = "Nonlinear" if nonlinear else "Linear mu_r=1000"
    print("=" * 70)
    print(f"  ELF Yoke ({geometry}): {mat_label}, {label}, NI={ni:.0f}")
    print(f"  VoxelCoefficient resolution={resolution}")
    print("=" * 70)
    print()

    bh_data = load_bh_curve()
    print(f"B-H curve: {len(bh_data)} points")

    # Build geometry + mesh
    print(f"\n1. Building ELF yoke geometry ({geometry})...")
    with GmshBuilder(model_name='elf_yoke_simkin', verbose=False) as gb:
        build_elf_yoke_geometry(gb, use_kelvin=use_kelvin,
                                maxh_gap=maxh_gap, maxh_iron=maxh_iron,
                                geometry=geometry)
        print("\n2. Generating mesh...")
        gb.generate(element_type='tet')
        mesh = gb.to_ngsolve_volume()

    fem_order = 2
    mesh.Curve(fem_order)
    print(f"   Mesh: {mesh.nv} vertices, {mesh.ne} elements")
    mats = sorted(set(mesh.GetMaterials()))
    print(f"   Materials: {mats}")
    print()

    # Solve
    print("3. Computing H_s from Radia coil...")
    rad.UtiDelAll()
    coil = build_radia_coil(ni)

    if use_kelvin:
        solver = ScalarPotentialSolver(
            mesh, iron_domains='iron', order=fem_order,
            kelvin_region='kelvin', kelvin_radius=AIR_R,
            kelvin_center=SPHERE_CENTER)
    else:
        solver = ScalarPotentialSolver(
            mesh, iron_domains='iron', order=fem_order)

    solver.set_source_from_radia(coil, resolution=resolution)
    print()

    if nonlinear:
        print("4. Newton + SymbolicEnergy solve...")
        t0 = time.time()
        phi_gf = solver.solve_nonlinear_newton(
            bh_data=bh_data,
            tol=1e-4,
            maxiter=80,
            dirichlet='outer')
        t_solve = time.time() - t0
    else:
        print("4. Linear solve (mu_r=1000)...")
        solver._mu_r = 1000.0
        t0 = time.time()
        phi_gf = solver.solve(dirichlet='outer')
        t_solve = time.time() - t0

    print(f"   Solve time: {t_solve:.1f}s")

    B_cf = solver.get_B()
    try:
        mip = mesh(0, 0, 0)
        B_fem = B_cf(mip)
        Bz = B_fem[2]
        B_mag = np.sqrt(sum(b**2 for b in B_fem))
        print(f"   B at origin: Bz = {Bz*1e3:.2f} mT")
        print(f"   |B| = {B_mag*1e3:.2f} mT")
    except Exception as e:
        print(f"   Cannot evaluate at origin: {e}")
        Bz = 0.0

    if nonlinear and ni == 20000:
        print(f"\n   Expected (ELF/Tosca): Bz = -995 mT")
        if abs(Bz) > 1e-6:
            print(f"   Error: {abs(Bz*1e3 - (-995)) / 995 * 100:.1f}%")
    elif not nonlinear and ni == 2000:
        print(f"\n   Expected (BEM): Bz ~ -232 mT")
        if abs(Bz) > 1e-6:
            print(f"   Error: {abs(Bz*1e3 - (-232)) / 232 * 100:.1f}%")

    rad.UtiDelAll()
    return Bz


def main():
    print("Testing ELF yoke geometry variants with Simkin solver")
    print("Comparing: rect vs taper (2D) vs frustum (3D)\n")

    results = {}
    mesh_params = dict(maxh_gap=0.004, maxh_iron=0.010)

    for geom in ['rect', 'taper', 'frustum']:
        print("\n" + "#" * 70)
        print(f"# Nonlinear NI=20000, Kelvin, geometry={geom}")
        print("#" * 70)
        try:
            Bz = run_test(use_kelvin=True, ni=20000.0, nonlinear=True,
                          geometry=geom, **mesh_params)
            results[geom] = Bz
        except Exception as e:
            print(f"   FAILED: {e}")
            import traceback
            traceback.print_exc()
            results[geom] = None

    # Summary
    print("\n\n" + "=" * 70)
    print("  SUMMARY: Geometry Comparison (Nonlinear NI=20000)")
    print("=" * 70)
    ref_mT = -945.0
    for geom in ['rect', 'taper', 'frustum']:
        Bz = results.get(geom)
        if Bz is not None:
            err = abs(Bz * 1e3 - ref_mT) / abs(ref_mT) * 100
            print(f"  {geom:8s}: Bz = {Bz*1e3:8.2f} mT  "
                  f"(error vs ~{ref_mT:.0f} mT: {err:.1f}%)")
        else:
            print(f"  {geom:8s}: FAILED")
    print("=" * 70)


if __name__ == '__main__':
    main()
