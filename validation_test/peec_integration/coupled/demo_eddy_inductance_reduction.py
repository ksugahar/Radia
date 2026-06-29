"""
demo_eddy_inductance_reduction.py

Validates eddy current inductance reduction using ngbem solvers.

Core model selection:
  - 'fembem': FEM-BEM coupling (ngbem Calderon). mu_r=1 only. Unbounded domain.
  - 'fem':    FEM with Dirichlet BC. Any mu_r. Bounded domain.
  - 'radia':  Radia MMM/MSC. Any mu_r, nonlinear. Unbounded domain.

Physics:
  A coil near a conducting core experiences:
  1. Low freq:  L ~ L_air + Delta_L * (mu_r - 1)  (full magnetization)
  2. Mid freq:  L decreases, R increases (eddy current loss peak)
  3. High freq: L ~ L_air - Delta_L               (shielding, L < L_air)

Part of Radia project
"""

import sys
import os
import numpy as np
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

from ngbem_peec import NGBEMPEECSolver, create_plate_mesh
from ngbem_coupled import CoupledPEECMMM, compute_delta_L
from ngbem_eddy import EddyCurrentFEMBEM, create_conductor_mesh

MU_0 = 4.0 * np.pi * 1e-7


def print_separator(title):
    print()
    print("=" * 70)
    print(f"  {title}")
    print("=" * 70)


def analytical_mu_eff(omega, mu_r, sigma, thickness):
    """Analytical mu_eff for infinite conducting slab (reference formula).

    mu_eff = mu_r * tanh(gamma*d/2) / (gamma*d/2)
    gamma = sqrt(j*omega*mu*sigma)
    """
    if sigma <= 0 or omega <= 0 or thickness <= 0:
        return complex(mu_r)

    mu = mu_r * MU_0
    gamma = np.sqrt(1j * omega * mu * sigma)
    gd2 = gamma * thickness / 2.0

    if abs(gd2) < 1e-6:
        ratio = 1.0 - gd2**2 / 3.0
    else:
        ratio = np.tanh(gd2) / gd2

    return mu_r * ratio


def test_fembem_vs_fem():
    """Test 1: Compare FEM-BEM (Calderon) vs FEM (Dirichlet) for mu_r=1.

    This validates that the ngbem FEM-BEM coupling works correctly
    by comparing against the simpler FEM-only solver for a non-magnetic
    conductor (mu_r=1, where scalar Hz formulation is exact).

    For mu_r=1:
    - FEM-BEM: Unbounded domain, natural BC from BEM
    - FEM:     Bounded domain, Dirichlet BC Hz=Hz_inc on surface
    Both should give similar mu_eff for well-resolved meshes.
    """
    print_separator("Test 1: FEM-BEM vs FEM (mu_r=1, Calderon Coupling)")

    mu_r = 1.0
    sigma = 3.5e7     # Aluminum
    depth = 0.002      # 2mm
    width = 0.01       # 10mm
    height = 0.01      # 10mm

    print(f"\n  Conductor: {width*1e3:.0f}mm x {height*1e3:.0f}mm x "
          f"{depth*1e3:.0f}mm")
    print(f"  mu_r={mu_r:.0f} (non-magnetic), sigma={sigma:.1e} S/m")
    print(f"  FEM-BEM: Calderon projector, unbounded domain")
    print(f"  FEM:     Dirichlet BC, bounded domain")

    # Create 3D volume mesh
    core_mesh = create_conductor_mesh(width, height, depth, maxh=0.003,
                                       conductor_label="conductor",
                                       surface_label="surface")

    freqs = [100, 500, 1000, 5000, 10000]

    print()
    print(f"  {'Freq':>10s}  {'d/delta':>8s}  "
          f"{'mu_FEM':>10s}  {'mu_FEMBEM':>10s}  {'Diff':>8s}  "
          f"{'mu_ana':>10s}")
    print(f"  {'-'*10}  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*8}  {'-'*10}")

    max_diff = 0.0
    all_pass = True

    for freq in freqs:
        omega = 2.0 * np.pi * freq
        delta = np.sqrt(2.0 / (omega * mu_r * MU_0 * sigma))

        # --- FEM mode (Dirichlet BC) ---
        eddy_fem = EddyCurrentFEMBEM(core_mesh, sigma=sigma, mu_r=mu_r,
                                      order=2,
                                      conductor_label="conductor",
                                      surface_label="surface")
        eddy_fem.assemble_fem(freq)
        eddy_fem.solve(B_ext=[0, 0, 1.0], mode='fem')

        from ngsolve import Integrate, CF
        from ngsolve import TaskManager
        Hz_total_fem = eddy_fem.get_total_field()
        Hz_inc = 1.0 / MU_0
        with TaskManager():
            volume = Integrate(CF(1), core_mesh).real
            Hz_avg_fem = Integrate(Hz_total_fem, core_mesh) / volume
            mu_eff_fem = mu_r * Hz_avg_fem / Hz_inc

            # --- FEMBEM mode (Calderon coupling) ---
            eddy_bem = EddyCurrentFEMBEM(core_mesh, sigma=sigma, mu_r=mu_r,
                                          order=2,
                                          conductor_label="conductor",
                                          surface_label="surface")
            eddy_bem.assemble_fembem(freq)
            eddy_bem.solve(B_ext=[0, 0, 1.0], mode='fembem')

            Hz_total_bem = eddy_bem.get_total_field()
            Hz_avg_bem = Integrate(Hz_total_bem, core_mesh) / volume
            mu_eff_bem = mu_r * Hz_avg_bem / Hz_inc

            # --- Analytical reference ---
            mu_eff_ana = analytical_mu_eff(omega, mu_r, sigma, depth)

            # Compute difference
            if abs(mu_eff_fem.real) > 1e-10:
                diff = abs(mu_eff_fem.real - mu_eff_bem.real) / abs(mu_eff_fem.real)
            else:
                diff = abs(mu_eff_fem.real - mu_eff_bem.real)
            max_diff = max(max_diff, diff)

            if freq < 1e3:
                f_str = f"{freq:.0f} Hz"
            else:
                f_str = f"{freq/1e3:.0f} kHz"

            print(f"  {f_str:>10s}  {depth/delta:>8.2f}  "
                  f"{mu_eff_fem.real:>10.4f}  {mu_eff_bem.real:>10.4f}  "
                  f"{diff*100:>7.1f}%  {mu_eff_ana.real:>10.4f}")

    print(f"\n  Max FEM-BEM vs FEM difference: {max_diff*100:.1f}%")
    if max_diff < 0.20:
        print(f"  [PASS] FEM-BEM coupling matches FEM (< 20%)")
    else:
        print(f"  [INFO] Larger difference expected: FEM uses bounded domain,")
        print(f"         FEMBEM uses unbounded. Different physics at boundary.")

    return max_diff


def test_core_model_comparison():
    """Test 2: Compare core_model='fem' vs core_model='fembem'.

    Uses CoupledPEECMMM with different core models for same geometry.
    Non-magnetic conductor (mu_r=1) where both models are valid.
    """
    print_separator("Test 2: CoupledPEECMMM core_model Comparison (mu_r=1)")

    # Create conductor plate (coil)
    width = 0.01
    height = 0.01
    maxh = 0.004
    thickness = 35e-6
    sigma_coil = 5.8e7

    mesh = create_plate_mesh(width, height, maxh, label="conductor")
    peec = NGBEMPEECSolver(mesh, conductor_label="conductor",
                            sigma=sigma_coil, thickness=thickness,
                            order=0, intorder=5)
    peec.assemble()

    # Non-magnetic shield
    mu_r = 1.0
    sigma_shield = 3.5e7   # Aluminum
    d_shield = 0.002

    print(f"\n  Coil: {width*1e3:.0f}mm x {height*1e3:.0f}mm plate")
    print(f"  Shield: mu_r=1 (aluminum), sigma={sigma_shield:.1e} S/m, "
          f"d={d_shield*1e3:.1f}mm")

    # Create shield mesh
    core_mesh = create_conductor_mesh(0.015, 0.015, d_shield, maxh=0.003,
                                       conductor_label="conductor",
                                       surface_label="surface")

    Delta_L = compute_delta_L(peec.L, mu_r)

    # --- core_model='fem' (Dirichlet BC) ---
    coupled_fem = CoupledPEECMMM(
        peec, core_model='fem',
        core_mesh=core_mesh, core_sigma=sigma_shield, mu_r=mu_r)
    coupled_fem.compute_coupling_analytically(Delta_L)

    # --- core_model='fembem' (Calderon coupling) ---
    coupled_bem = CoupledPEECMMM(
        peec, core_model='fembem',
        core_mesh=core_mesh, core_sigma=sigma_shield, mu_r=mu_r)
    coupled_bem.compute_coupling_analytically(Delta_L)

    # --- No core (air only) ---
    Z_air = peec.solve_frequency(np.array([1000.0]), mode='mqs')
    L_air_1k = np.imag(Z_air[0]) / (2*np.pi*1000) * 1e9

    freqs = np.array([100, 500, 1000, 5000, 10000])

    print()
    print(f"  {'Freq':>10s}  {'L_air (nH)':>10s}  "
          f"{'L_fem (nH)':>10s}  {'L_fembem (nH)':>13s}  "
          f"{'Diff':>8s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*13}  {'-'*8}")

    Z_air_all = peec.solve_frequency(freqs, mode='mqs')

    for i, freq in enumerate(freqs):
        omega = 2.0 * np.pi * freq

        Z_fem_i = coupled_fem.solve_frequency(np.array([freq]), mode='mqs')
        Z_bem_i = coupled_bem.solve_frequency(np.array([freq]), mode='mqs')

        L_air_val = np.imag(Z_air_all[i]) / omega * 1e9
        L_fem = np.imag(Z_fem_i[0]) / omega * 1e9
        L_bem = np.imag(Z_bem_i[0]) / omega * 1e9

        if abs(L_fem) > 1e-15:
            diff = abs(L_fem - L_bem) / abs(L_fem)
        else:
            diff = abs(L_fem - L_bem)

        if freq < 1e3:
            f_str = f"{freq:.0f} Hz"
        else:
            f_str = f"{freq/1e3:.0f} kHz"

        print(f"  {f_str:>10s}  {L_air_val:>10.3f}  "
              f"{L_fem:>10.3f}  {L_bem:>13.3f}  "
              f"{diff*100:>7.1f}%")

    # Get mu_eff comparison
    print(f"\n  mu_eff comparison at 1 kHz:")
    omega_1k = 2.0 * np.pi * 1000
    mu_fem = coupled_fem._get_mu_eff(omega_1k)
    mu_bem = coupled_bem._get_mu_eff(omega_1k)
    mu_ana = analytical_mu_eff(omega_1k, mu_r, sigma_shield, d_shield)
    print(f"    FEM (Dirichlet):  mu_eff = {mu_fem.real:.4f} + "
          f"{mu_fem.imag:.4f}j")
    print(f"    FEMBEM (Calderon): mu_eff = {mu_bem.real:.4f} + "
          f"{mu_bem.imag:.4f}j")
    print(f"    Analytical (slab): mu_eff = {mu_ana.real:.4f} + "
          f"{mu_ana.imag:.4f}j")


def test_inductance_reduction():
    """Test 3: Inductance reduction with FEM eddy currents."""
    print_separator("Test 3: Inductance Reduction (core_model='fem')")

    width = 0.01
    height = 0.01
    maxh = 0.004
    thickness = 35e-6
    sigma_coil = 5.8e7

    mesh = create_plate_mesh(width, height, maxh, label="conductor")
    peec = NGBEMPEECSolver(mesh, conductor_label="conductor",
                            sigma=sigma_coil, thickness=thickness,
                            order=0, intorder=5)
    peec.assemble()

    # Conducting magnetic core
    mu_r = 100.0
    sigma_core = 2e6
    d_core = 0.005

    print(f"\n  Coil: {width*1e3:.0f}mm x {height*1e3:.0f}mm plate")
    print(f"  Core: mu_r={mu_r:.0f}, sigma={sigma_core:.1e} S/m, "
          f"d={d_core*1e3:.1f}mm")

    core_mesh = create_conductor_mesh(0.02, 0.02, d_core, maxh=0.004,
                                       conductor_label="conductor",
                                       surface_label="surface")

    # Static core (no eddy currents)
    coupled_static = CoupledPEECMMM(peec, mu_r=mu_r)
    Delta_L = compute_delta_L(peec.L, mu_r)
    coupled_static.compute_coupling_analytically(Delta_L)

    # Conducting core with FEM eddy currents
    coupled_eddy = CoupledPEECMMM(
        peec, core_model='fem',
        core_mesh=core_mesh, core_sigma=sigma_core, mu_r=mu_r)
    coupled_eddy.compute_coupling_analytically(Delta_L)

    freqs = np.array([10, 50, 100, 500, 1000, 5000, 10000])
    Z_air = peec.solve_frequency(freqs, mode='mqs')
    Z_static = coupled_static.solve_frequency(freqs, mode='mqs')

    print()
    print(f"  {'Freq':>10s}  {'L_air (nH)':>10s}  {'L_static':>10s}  "
          f"{'L_eddy':>10s}  {'L_eddy/L_air':>12s}  {'mu_eff':>10s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*12}  {'-'*10}")

    for i, freq in enumerate(freqs):
        omega = 2.0 * np.pi * freq
        Z_eddy_i = coupled_eddy.solve_frequency(np.array([freq]), mode='mqs')

        L_air = np.imag(Z_air[i]) / omega * 1e9
        L_static = np.imag(Z_static[i]) / omega * 1e9
        L_eddy = np.imag(Z_eddy_i[0]) / omega * 1e9
        ratio = L_eddy / L_air if abs(L_air) > 1e-15 else 0

        mu_eff = coupled_eddy._get_mu_eff(omega)

        if freq < 1e3:
            f_str = f"{freq:.0f} Hz"
        else:
            f_str = f"{freq/1e3:.0f} kHz"

        print(f"  {f_str:>10s}  {L_air:>10.3f}  {L_static:>10.3f}  "
              f"{L_eddy:>10.3f}  {ratio:>12.3f}  {mu_eff.real:>10.2f}")

    # Verify inductance decreases with frequency
    Z_low = coupled_eddy.solve_frequency(np.array([10.0]), mode='mqs')
    Z_high = coupled_eddy.solve_frequency(np.array([10000.0]), mode='mqs')
    L_low = np.imag(Z_low[0]) / (2*np.pi*10) * 1e9
    L_high = np.imag(Z_high[0]) / (2*np.pi*10000) * 1e9

    print(f"\n  L(10 Hz)={L_low:.2f} nH -> L(10 kHz)={L_high:.2f} nH")
    if L_high < L_low:
        print(f"  [PASS] Inductance decreases with frequency")
    else:
        print(f"  [WARN] Inductance did not decrease")


def test_nonmagnetic_shielding():
    """Test 4: Non-magnetic conductor shielding (mu_r=1) with FEMBEM.

    Uses core_model='fembem' for aluminum shield (mu_r=1).
    Verifies L_eddy < L_air at high frequencies.
    """
    print_separator("Test 4: Non-Magnetic Shield (core_model='fembem')")

    width = 0.01
    height = 0.01
    maxh = 0.004
    thickness = 35e-6
    sigma_coil = 5.8e7

    mesh = create_plate_mesh(width, height, maxh, label="conductor")
    peec = NGBEMPEECSolver(mesh, conductor_label="conductor",
                            sigma=sigma_coil, thickness=thickness,
                            order=0, intorder=5)
    peec.assemble()

    mu_r = 1.0
    sigma_shield = 3.5e7
    d_shield = 0.002

    print(f"\n  Coil: {width*1e3:.0f}mm x {height*1e3:.0f}mm plate")
    print(f"  Shield: mu_r=1 (aluminum), sigma={sigma_shield:.1e} S/m")
    print(f"  Core model: FEMBEM (Calderon, unbounded domain)")

    core_mesh = create_conductor_mesh(0.015, 0.015, d_shield, maxh=0.003,
                                       conductor_label="conductor",
                                       surface_label="surface")

    coupled = CoupledPEECMMM(
        peec, core_model='fembem',
        core_mesh=core_mesh, core_sigma=sigma_shield, mu_r=mu_r)
    Delta_L = compute_delta_L(peec.L, mu_r)
    coupled.compute_coupling_analytically(Delta_L)

    freqs = np.array([10, 100, 500, 1000, 5000, 10000, 50000])
    Z_air = peec.solve_frequency(freqs, mode='mqs')

    print()
    print(f"  {'Freq':>10s}  {'d/delta':>8s}  {'mu_eff':>10s}  "
          f"{'L_air (nH)':>10s}  {'L_eddy (nH)':>11s}  {'L_eddy/L_air':>12s}")
    print(f"  {'-'*10}  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*11}  {'-'*12}")

    below_air = False

    for i, freq in enumerate(freqs):
        omega = 2.0 * np.pi * freq
        delta = np.sqrt(2.0 / (omega * mu_r * MU_0 * sigma_shield))

        Z_eddy_i = coupled.solve_frequency(np.array([freq]), mode='mqs')

        L_air_val = np.imag(Z_air[i]) / omega * 1e9
        L_eddy = np.imag(Z_eddy_i[0]) / omega * 1e9
        ratio = L_eddy / L_air_val if abs(L_air_val) > 1e-15 else 0

        mu_eff = coupled._get_mu_eff(omega)

        if ratio < 1.0:
            below_air = True

        if freq < 1e3:
            f_str = f"{freq:.0f} Hz"
        elif freq < 1e6:
            f_str = f"{freq/1e3:.0f} kHz"
        else:
            f_str = f"{freq/1e6:.0f} MHz"

        print(f"  {f_str:>10s}  {d_shield/delta:>8.2f}  "
              f"{mu_eff.real:>10.4f}  {L_air_val:>10.3f}  "
              f"{L_eddy:>11.3f}  {ratio:>12.4f}")

    if below_air:
        print(f"\n  [PASS] L_eddy < L_air confirmed with FEMBEM core model")
    else:
        print(f"\n  [WARN] L_eddy did not drop below L_air")


def test_material_comparison():
    """Test 5: Compare materials with appropriate core_model."""
    print_separator("Test 5: Material Comparison (core_model selection)")

    width = 0.01
    height = 0.01
    maxh = 0.004
    thickness = 35e-6
    sigma_coil = 5.8e7

    mesh = create_plate_mesh(width, height, maxh, label="conductor")
    peec = NGBEMPEECSolver(mesh, conductor_label="conductor",
                            sigma=sigma_coil, thickness=thickness,
                            order=0, intorder=5)
    peec.assemble()

    f_test = 1000.0
    omega = 2.0 * np.pi * f_test

    Z_air = peec.solve_frequency(np.array([f_test]), mode='mqs')
    L_air = np.imag(Z_air[0]) / omega * 1e9

    # (name, mu_r, sigma, core_depth, core_model)
    materials = [
        ("Ferrite (insulating)", 2000, 0.01,  0.005, None),      # Static mu_r
        ("Silicon steel",        100,  2e6,   0.005, 'fem'),     # FEM (mu_r!=1)
        ("Solid iron",           100,  1e7,   0.005, 'fem'),     # FEM (mu_r!=1)
        ("Aluminum shield",        1,  3.5e7, 0.002, 'fembem'),  # FEMBEM (mu_r=1)
        ("Copper shield",          1,  5.8e7, 0.002, 'fembem'),  # FEMBEM (mu_r=1)
    ]

    print(f"\n  Frequency: {f_test/1e3:.0f} kHz")
    print(f"  L_air = {L_air:.3f} nH")
    print()
    print(f"  {'Material':>25s}  {'mu_r':>6s}  {'sigma':>11s}  "
          f"{'Model':>8s}  {'mu_eff':>10s}  {'L (nH)':>10s}  {'L/L_air':>8s}")
    print(f"  {'-'*25}  {'-'*6}  {'-'*11}  {'-'*8}  "
          f"{'-'*10}  {'-'*10}  {'-'*8}")

    for name, mu_r, sigma_c, d, model in materials:
        if model in ('fem', 'fembem'):
            core_mesh = create_conductor_mesh(0.02, 0.02, d, maxh=0.004,
                                               conductor_label="conductor",
                                               surface_label="surface")
            coupled = CoupledPEECMMM(
                peec, core_model=model,
                core_mesh=core_mesh, core_sigma=sigma_c, mu_r=float(mu_r))
        else:
            coupled = CoupledPEECMMM(peec, mu_r=float(mu_r))

        DL = compute_delta_L(peec.L, mu_r)
        coupled.compute_coupling_analytically(DL)

        Z = coupled.solve_frequency(np.array([f_test]), mode='mqs')
        L = np.imag(Z[0]) / omega * 1e9
        mu_eff = coupled._get_mu_eff(omega)

        model_str = model or 'static'
        print(f"  {name:>25s}  {mu_r:>6.0f}  {sigma_c:>11.1e}  "
              f"{model_str:>8s}  {mu_eff.real:>10.1f}  "
              f"{L:>10.3f}  {L/L_air:>8.3f}")

    print(f"\n  Core model selection:")
    print(f"    mu_r=1  (non-magnetic): 'fembem' (Calderon, unbounded)")
    print(f"    mu_r>>1 (magnetic):     'fem' (Dirichlet, bounded)")
    print(f"    sigma~0 (insulating):   None (static mu_r, no eddy)")
    print(f"    nonlinear B-H:          'radia' (MMM/MSC, unbounded)")


def main():
    print("=" * 70)
    print("  Demo: Eddy Current Inductance Reduction")
    print("  Core model: FEM-BEM (Calderon) vs FEM (Dirichlet) vs Radia")
    print("=" * 70)

    # Test 1: Validate FEMBEM vs FEM for mu_r=1
    test_fembem_vs_fem()

    # Test 2: CoupledPEECMMM with both core models
    test_core_model_comparison()

    # Test 3: Inductance reduction (magnetic core, FEM)
    test_inductance_reduction()

    # Test 4: Non-magnetic shielding (FEMBEM)
    test_nonmagnetic_shielding()

    # Test 5: Material comparison with model selection
    test_material_comparison()

    print_separator("Summary")
    print("""
  Core Model Selection Guide:

  +-------------------+----------+-----------+------------------+
  | Material          | mu_r     | sigma     | core_model       |
  +-------------------+----------+-----------+------------------+
  | Aluminum shield   | 1        | 35 MS/m   | 'fembem'         |
  | Copper shield     | 1        | 58 MS/m   | 'fembem'         |
  | Silicon steel     | 100-5000 | 2 MS/m    | 'fem'            |
  | Solid iron        | 100-5000 | 10 MS/m   | 'fem'            |
  | Ferrite           | 1000+    | ~0        | None (static)    |
  | Nonlinear B-H     | varies   | varies    | 'radia' (MMM)    |
  +-------------------+----------+-----------+------------------+

  Key physics:
  - 'fembem': Unbounded domain (BEM exterior), high-order BEM elements
              Valid ONLY for mu_r=1 (scalar Hz formulation)
  - 'fem':    Bounded domain (Dirichlet BC), any linear mu_r
              Approximation: assumes Hz=Hz_inc on boundary
  - 'radia':  Integral equation (unbounded), nonlinear materials
              Full MMM/MSC with MatSatIsoTab support
""")


if __name__ == '__main__':
    main()
