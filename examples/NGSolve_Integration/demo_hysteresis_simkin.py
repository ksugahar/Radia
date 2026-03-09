#!/usr/bin/env python
"""
C-Type Electromagnet - Hysteresis Simkin-Trowbridge Solver

Same geometry as demo_ctype_simkin.py but with energy-based
hysteresis material (B-input Play model, Egger formulation).

Uses ScalarPotentialSolver.solve_hysteresis() which performs Picard
iteration using rad.MatMvsH() (Inverse: H -> M with internal state).

Demonstrates:
  1. DC magnetization (single NI level)
  2. AC minor loop (quasi-static NI cycling)
  3. B-H loop at gap center
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

NI_MAX = 20000.0


def make_hysteresis_factory():
    """Return a factory function that creates energy-based hysteresis materials.

    K=10 play operators, designed to mimic typical soft iron with
    moderate hysteresis (coercivity ~50 A/m).
    """
    import radia as rad

    K = 10
    As = np.full(K, 8000.0)     # Saturation slope [A/m]
    Js = np.full(K, 0.18)       # Saturation polarization per operator [T]
    chi = np.linspace(0, 200, K)  # Pinning strength [A/m]
    eps = 1e-8

    # Generate table-based shape functions from analytical parameters
    tables = []
    for k in range(K):
        r = np.linspace(0, Js[k], 200)
        f = As[k] * np.tan(np.pi / 2 * r / Js[k])
        f[-1] = f[-2] * 2
        tables.append((r, f))

    print(f"   Energy hysteresis: K={K}, As={As[0]:.0f}, Js={Js[0]:.2f}")
    print(f"   Pinning: chi = [0, ..., {chi[-1]:.0f}] A/m")

    def factory():
        return rad.MatEnergyHysteresis(K, chi, tables, eps)

    return factory


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


def build_radia_coil(NI):
    """Build Radia racetrack coil for H_s source field."""
    coil_dir = os.path.join(repo_root, 'examples',
                            'c_type_electromagnet', 'mu=1000')
    if coil_dir not in sys.path:
        sys.path.insert(0, coil_dir)
    from coil_model import create_racetrack_coil
    return create_racetrack_coil(NI)


def main():
    print("=" * 70)
    print("  C-Type Electromagnet: HYSTERESIS Simkin + Kelvin")
    print("  Energy-based B-input Play Model (Egger formulation)")
    print("=" * 70)
    print()

    import radia as rad
    from radia.gmsh_builder import GmshBuilder
    from radia.scalar_potential_solver import ScalarPotentialSolver

    # =====================================================================
    # Step 1: Build geometry + mesh
    # =====================================================================
    print("1. Building geometry...", flush=True)

    with GmshBuilder(model_name='c_type_hysteresis', verbose=False) as gb:
        build_simkin_geometry(gb, maxh_gap=0.008, maxh_iron=0.020,
                              maxh_air=0.090)
        print()
        print("2. Generating mesh...")
        gb.generate(element_type='tet')
        mesh = gb.to_ngsolve_volume()

    fem_order = 2
    mesh.Curve(fem_order)
    print(f"   Mesh: {mesh.nv} vertices, {mesh.ne} elements")
    print()

    # =====================================================================
    # Step 2: Create hysteresis material factory
    # =====================================================================
    print("3. Creating hysteresis material factory...")
    mat_factory = make_hysteresis_factory()
    print()

    # =====================================================================
    # Step 3: DC magnetization at NI_MAX
    # =====================================================================
    print("4. DC magnetization (NI = {:.0f} AT)...".format(NI_MAX))
    rad.UtiDelAll()
    coil = build_radia_coil(NI_MAX)

    solver = ScalarPotentialSolver(
        mesh, iron_domains='iron', order=fem_order,
        kelvin_region='kelvin', kelvin_radius=AIR_R,
        kelvin_center=SPHERE_CENTER)
    solver.set_source_from_radia(coil, resolution=51)

    t0 = time.time()
    phi_gf = solver.solve_hysteresis(
        mat_factory,
        tol=1e-3,
        maxiter=50,
        alpha=500.0,
        dirichlet='outer')
    t_dc = time.time() - t0
    print(f"   DC solve: {t_dc:.1f}s")

    B_cf = solver.get_B()
    try:
        mip = mesh(0, 0, 0)
        B_dc = np.array([float(B_cf(mip)[i]) for i in range(3)])
        print(f"   B at gap center: Bz={B_dc[2]*1e3:.2f} mT, "
              f"|B|={np.linalg.norm(B_dc)*1e3:.2f} mT")
    except Exception as e:
        print(f"   Cannot evaluate at origin: {e}")
        B_dc = np.zeros(3)
    print()

    # =====================================================================
    # Step 4: AC loop (quasi-static NI cycling)
    # Per-element material handles persist across steps via solver._hys_handles
    # =====================================================================
    print("5. AC quasi-static loop...")
    n_steps_up = 4
    n_steps_down = 8
    n_steps_back = 4

    # Up: 0 -> NI_MAX
    NI_up = np.linspace(0, NI_MAX, n_steps_up + 1)[1:]
    # Down: NI_MAX -> -NI_MAX
    NI_down = np.linspace(NI_MAX, -NI_MAX, n_steps_down + 1)[1:]
    # Back up: -NI_MAX -> NI_MAX
    NI_back = np.linspace(-NI_MAX, NI_MAX, n_steps_back + 1)[1:]

    NI_sequence = np.concatenate([NI_up, NI_down, NI_back])

    # Track B at gap center
    Bz_history = []
    NI_history = []

    # Use a SINGLE solver with persistent per-element handles
    ac_solver = ScalarPotentialSolver(
        mesh, iron_domains='iron', order=fem_order,
        kelvin_region='kelvin', kelvin_radius=AIR_R,
        kelvin_center=SPHERE_CENTER)

    for i, ni_val in enumerate(NI_sequence):
        # Update coil source (do NOT call UtiDelAll — preserves material handles)
        coil = build_radia_coil(ni_val)
        ac_solver.set_source_from_radia(coil, resolution=31)

        ac_solver.solve_hysteresis(
            mat_factory,
            tol=5e-3,
            maxiter=30,
            alpha=500.0,
            dirichlet='outer',
            verbose=False)

        B_cf_step = ac_solver.get_B()
        try:
            mip = mesh(0, 0, 0)
            Bz = float(B_cf_step(mip)[2])
        except Exception:
            Bz = 0.0

        Bz_history.append(Bz)
        NI_history.append(ni_val)

        if (i + 1) % 4 == 0 or i == 0:
            print(f"   Step {i+1:3d}/{len(NI_sequence)}: "
                  f"NI={ni_val:+8.0f} AT, Bz={Bz*1e3:+8.2f} mT")

    print()

    # =====================================================================
    # Step 5: Plot B-H loop
    # =====================================================================
    print("6. Saving B-NI loop plot...")
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        ax.plot(NI_history, [bz * 1e3 for bz in Bz_history], 'b.-',
                linewidth=1.5, markersize=4)
        ax.set_xlabel('NI [A-turns]')
        ax.set_ylabel('Bz at gap center [mT]')
        ax.set_title('C-Type Electromagnet: Hysteresis Loop (Simkin + Energy Play)')
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='k', linewidth=0.5)
        ax.axvline(x=0, color='k', linewidth=0.5)

        plot_file = os.path.join(script_dir, 'demo_hysteresis_simkin.png')
        fig.savefig(plot_file, dpi=150, bbox_inches='tight')
        print(f"   Written: {plot_file}")
        plt.close(fig)
    except Exception as e:
        print(f"   Plot skipped: {e}")

    # =====================================================================
    # Summary
    # =====================================================================
    print()
    print("=" * 70)
    print("  Summary")
    print("=" * 70)
    print(f"  FEM mesh: {mesh.nv} vertices, {mesh.ne} elements")
    print(f"  DC magnetization: Bz = {B_dc[2]*1e3:.2f} mT at NI = {NI_MAX:.0f} AT")
    if len(Bz_history) > 0:
        Bz_arr = np.array(Bz_history) * 1e3
        print(f"  AC loop: Bz range = [{Bz_arr.min():.2f}, {Bz_arr.max():.2f}] mT")
    print()
    print("  Method: H = H_s - grad(phi)  [Simkin 1979]")
    print("  Material: Energy-based hysteresis (B-input Play, Egger)")
    print("  Open boundary: Kelvin transform (R = {:.0f} mm)".format(AIR_R*1e3))
    print("=" * 70)


if __name__ == '__main__':
    main()
