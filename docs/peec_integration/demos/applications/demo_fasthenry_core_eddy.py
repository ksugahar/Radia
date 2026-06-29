"""
demo_fasthenry_core_eddy.py

Demonstration of eddy current effects in magnetic core via ESIM + PEEC.

Problem: Wire above conductive iron core at various frequencies.
The core has sigma > 0, so eddy currents attenuate the field penetration.

Three models compared:
  1. Static MMM (CoupledPEECSolver) - Delta_L = constant (no eddy currents)
  2. ESIM mu_eff(f) - Delta_L(f) = Delta_L_static * mu_eff(f)/mu_r
  3. Analytical Dowell limit for comparison

Physics:
  - At low freq: skin depth >> core size -> full penetration -> L ~ L_static
  - At high freq: skin depth << core size -> shielding -> L decreases
  - ESIM captures this transition correctly via 1D cell problem

Part of Radia project
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad
from fasthenry_parser import FastHenryParser
from peec_coupled import CoupledPEECSolver
from peec_topology import PEECCircuitSolver

MU_0 = 4e-7 * np.pi

# Try importing ESIM (requires scipy)
try:
    from esim_cell_problem import ESIMFiniteSlabSolver
    ESIM_AVAILABLE = True
except ImportError:
    ESIM_AVAILABLE = False
    print("WARNING: scipy not available, ESIM comparison disabled")


def compute_static_coupling():
    """Compute static Delta_L using Radia MMM."""
    rad.UtiDelAll()

    # Wire: 100 mm along x, 1x1 mm cross-section
    inp = """\
.Units mm
.default sigma=1e6

N1 x=0 y=0 z=0
N2 x=100 y=0 z=0

E1 N1 N2 w=1 h=1

.external N1 N2

.magnetic
  type=box
  center=50,10,0
  size=60,10,10
  divisions=3,1,1
  mu_r=1000
.endmagnetic

.end
"""
    parser = FastHenryParser()
    parser.parse_string(inp)
    result = parser.solve(freqs=np.array([1.0]))

    L_air = result['topology']['L'][0, 0]
    Delta_L_static = result['Delta_L'][0, 0]

    return L_air, Delta_L_static


def compute_esim_mu_eff(freqs, mu_r, sigma, half_thickness):
    """
    Compute frequency-dependent effective permeability using ESIM.

    At each frequency, solve the 1D cell problem inside the core slab
    to get how much of the core is actually magnetized.

    mu_eff(f) = mu_r * (average B inside) / (mu_r * mu_0 * H_surface)
              = effectively reduced permeability due to skin effect
    """
    if not ESIM_AVAILABLE:
        return np.full(len(freqs), mu_r)

    mu_eff = np.zeros(len(freqs))

    for i, freq in enumerate(freqs):
        if freq < 0.1:
            # DC: full penetration
            mu_eff[i] = mu_r
            continue

        solver = ESIMFiniteSlabSolver(
            half_thickness=half_thickness,
            sigma=sigma,
            frequency=freq,
            mu_r=mu_r,
            n_nodes=200,
        )

        # Solve for H(z) distribution inside slab
        H0 = 1.0  # unit surface field
        result = solver.solve(H0)
        H_dist = result['H_solution']
        z = result['mesh_points']

        # mu_eff = average|H(z)| / |H(0)|
        # This gives how much of the core is effectively magnetized
        # At DC: H uniform -> mu_eff = mu_r
        # At high freq: H decays inside -> mu_eff < mu_r
        H_avg = np.trapezoid(np.abs(H_dist), z) / half_thickness
        mu_eff[i] = mu_r * (H_avg / abs(H0))

    return mu_eff


def compute_analytical_mu_eff(freqs, mu_r, sigma, half_thickness):
    """
    Analytical mu_eff using Dowell/tanh formula for linear slab.

    For a slab of thickness 2a in uniform tangential field:
      H(z) = H0 * cosh(gamma*(a-z)) / cosh(gamma*a)
      <H> = H0 * tanh(gamma*a) / (gamma*a)

    mu_eff = mu_r * |tanh(gamma*a) / (gamma*a)|
    """
    mu_eff = np.zeros(len(freqs))

    for i, freq in enumerate(freqs):
        if freq < 0.1:
            mu_eff[i] = mu_r
            continue

        omega = 2 * np.pi * freq
        mu = mu_r * MU_0
        gamma = np.sqrt(1j * omega * mu * sigma)
        ga = gamma * half_thickness

        # tanh(ga)/ga gives average field ratio
        if abs(ga) < 1e-6:
            ratio = 1.0
        else:
            ratio = abs(np.tanh(ga) / ga)

        mu_eff[i] = mu_r * ratio

    return mu_eff


def main():
    print()
    print("*" * 65)
    print("  Eddy Current in Magnetic Core: ESIM + PEEC Demo")
    print("*" * 65)
    print()

    # Core parameters
    mu_r = 1000
    sigma_core = 1e6  # S/m (iron-like conductivity)
    core_thickness = 0.010  # 10 mm total
    half_thickness = core_thickness / 2.0

    print(f"  Core: 60x10x10 mm iron block")
    print(f"  mu_r = {mu_r}, sigma = {sigma_core:.0e} S/m")
    print(f"  half_thickness = {half_thickness*1e3:.1f} mm")
    print()

    # Skin depth reference
    print("  Skin depth reference:")
    ref_freqs = [50, 100, 1000, 10000, 100000]
    for f in ref_freqs:
        delta = np.sqrt(2.0 / (2 * np.pi * f * mu_r * MU_0 * sigma_core))
        ratio = delta / half_thickness
        print(f"    f={f:>7.0f} Hz: delta={delta*1e3:.3f} mm, "
              f"delta/a={ratio:.3f} {'(full penetration)' if ratio > 2 else '(partial)' if ratio > 0.1 else '(skin only)'}")
    print()

    # Step 1: Static coupling (frequency-independent)
    print("=" * 65)
    print("Step 1: Static MMM Coupling (Delta_L = const)")
    print("=" * 65)
    L_air, Delta_L_static = compute_static_coupling()
    print(f"  L_air = {L_air*1e9:.2f} nH")
    print(f"  Delta_L_static = {Delta_L_static*1e9:.2f} nH "
          f"({Delta_L_static/L_air*100:.1f}% of L_air)")
    print()

    # Step 2: Frequency-dependent mu_eff
    freqs = np.logspace(1, 6, 51)  # 10 Hz to 1 MHz

    print("=" * 65)
    print("Step 2: ESIM mu_eff(f) vs Analytical")
    print("=" * 65)

    mu_eff_esim = compute_esim_mu_eff(freqs, mu_r, sigma_core, half_thickness)
    mu_eff_analytical = compute_analytical_mu_eff(freqs, mu_r, sigma_core, half_thickness)

    # Step 3: Compute frequency-dependent Delta_L and total L
    # Delta_L(f) = Delta_L_static * mu_eff(f) / mu_r
    Delta_L_f_esim = Delta_L_static * mu_eff_esim / mu_r
    Delta_L_f_analytical = Delta_L_static * mu_eff_analytical / mu_r
    L_total_static = L_air + Delta_L_static  # no eddy current
    L_total_esim = L_air + Delta_L_f_esim    # with eddy current
    L_total_analytical = L_air + Delta_L_f_analytical

    # Print comparison at selected frequencies
    print(f"\n  {'f(Hz)':>10s}  {'delta/a':>8s}  {'mu_eff':>8s}  "
          f"{'mu_ana':>8s}  {'L_static':>10s}  {'L_esim':>10s}  {'L_ana':>10s}")
    print(f"  {'-'*10}  {'-'*8}  {'-'*8}  {'-'*8}  "
          f"{'-'*10}  {'-'*10}  {'-'*10}")

    display_freqs = [10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 1000000]
    for f_target in display_freqs:
        idx = np.argmin(np.abs(freqs - f_target))
        f = freqs[idx]
        delta = np.sqrt(2.0 / (2 * np.pi * f * mu_r * MU_0 * sigma_core))
        ratio = delta / half_thickness
        print(f"  {f:10.0f}  {ratio:8.3f}  {mu_eff_esim[idx]:8.1f}  "
              f"{mu_eff_analytical[idx]:8.1f}  "
              f"{L_total_static*1e9:10.2f}  "
              f"{L_total_esim[idx]*1e9:10.2f}  "
              f"{L_total_analytical[idx]*1e9:10.2f}")

    print(f"\n  (All inductances in nH)")

    # Step 4: Show the key physics
    print()
    print("=" * 65)
    print("KEY PHYSICS")
    print("=" * 65)
    print(f"""
  Static MMM (no eddy current in core):
    L = L_air + Delta_L = {L_total_static*1e9:.2f} nH (constant for all f)

  With core eddy current (ESIM):
    L(f) = L_air + Delta_L * mu_eff(f)/mu_r
    L(10 Hz)    = {L_total_esim[0]*1e9:.2f} nH  (full penetration, ~ static)
    L(1 kHz)    = {L_total_esim[np.argmin(np.abs(freqs-1000))]*1e9:.2f} nH  (partial penetration)
    L(100 kHz)  = {L_total_esim[np.argmin(np.abs(freqs-1e5))]*1e9:.2f} nH  (skin only)
    L(1 MHz)    = {L_total_esim[-1]*1e9:.2f} nH  (nearly air)

  Transition frequency (delta = half_thickness):
    f_transition = 2/(mu*sigma*a^2) = {2.0/(mu_r*MU_0*sigma_core*half_thickness**2):.0f} Hz

  At high frequency: L -> L_air ({L_air*1e9:.2f} nH)
  Core is completely shielded, as if no magnetic material present.
""")

    print("=" * 65)
    print("IMPLEMENTATION APPROACH")
    print("=" * 65)
    print("""
  Current (demo_fasthenry_magnetic_core.py):
    CoupledPEECSolver  ->  Static Delta_L  ->  L = const
    (+ mu"_r loss tangent for ferrite)

  With core eddy currents (this demo):
    Step 1: CoupledPEECSolver -> Delta_L_static (one-time Radia Solve)
    Step 2: ESIMFiniteSlabSolver -> mu_eff(f) at each frequency
    Step 3: Delta_L(f) = Delta_L_static * mu_eff(f) / mu_r
    Step 4: L(f) = L_air + Delta_L(f)

  This gives frequency-dependent inductance from core eddy currents
  using only Radia MMM (1 static solve) + ESIM (fast 1D solve per freq).
  No ngbem required.

  Applicable to:
    - Iron cores (sigma ~ 1e6 S/m, mu_r ~ 1000)
    - Laminated cores (effective sigma from lamination)
    - Ferrite at high frequency (sigma ~ 0.01 S/m, but mu"_r more important)
""")


if __name__ == '__main__':
    main()
