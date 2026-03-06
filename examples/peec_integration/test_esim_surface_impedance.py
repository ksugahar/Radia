"""
Test ESIM for Surface Impedance calculation (not R_ac/R_dc).

ESIM (Effective Surface Impedance Method) computes Z_s, which is used in SIBC.
This is the correct approach for skin effect in PEEC.

The key insight:
- ESIM with dH/dz(a) = 0 BC computes SURFACE IMPEDANCE Z_s
- Z_s relates surface E field to surface H field: Z_s = E(0) / H(0)
- This Z_s is used directly in SIBC, not converted to R_ac/R_dc

For PEEC conductor resistance:
1. Compute Z_s from ESIM cell problem
2. Use Z_s in the conductor impedance: Z_conductor = Z_s * (perimeter / length)
3. Or use: R_conductor = Re(Z_s) * (perimeter / length)

This approach is correct for:
- Thick conductors (skin depth << dimension)
- Transition region (skin depth ~ dimension)
- Nonlinear materials (mu(H) varies)
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

from esim_cell_problem import ESIMFiniteSlabSolver, ESIMCellProblemSolver
from scipy.constants import mu_0


def analytical_surface_impedance(sigma, omega, mu):
    """
    Analytical surface impedance for semi-infinite conductor.

    Z_s = sqrt(j*omega*mu / sigma) = (1+j) / (sigma * delta)

    where delta = sqrt(2 / (omega * mu * sigma))
    """
    if omega < 1e-10:
        return 0.0

    delta = np.sqrt(2.0 / (omega * mu * sigma))
    return (1 + 1j) / (sigma * delta)


def analytical_finite_slab_Zs(xi, rho):
    """
    Analytical surface impedance for finite slab with symmetric BC.

    For slab with half-thickness a and symmetric excitation:
    Z_s = rho * gamma * tanh(gamma * a)

    where gamma = (1+j) / delta, xi = a / delta

    Returns Z_s (Ohm per unit area).
    """
    if xi < 0.001:
        # DC limit: Z_s -> 0 (no skin effect)
        return 0.0

    gamma_a = complex(1, 1) * xi  # gamma * a
    Z_s_normalized = gamma_a * np.tanh(gamma_a)  # Z_s / (rho/a)

    # Z_s = rho * gamma * tanh(gamma*a) = (rho/a) * (gamma*a) * tanh(gamma*a)
    # If we normalize by rho/a, we get gamma*a * tanh(gamma*a)
    # But we want Z_s in Ohm, so:
    # Z_s = (rho/a) * Z_s_normalized

    # Wait, let me redo this:
    # gamma = (1+j) / delta
    # gamma * a = (1+j) * (a/delta) = (1+j) * xi
    # Z_s = rho * gamma * tanh(gamma*a)
    #     = rho * (1+j) / delta * tanh((1+j)*xi)
    #     = rho * (1+j) * xi / a * tanh((1+j)*xi)
    #     = (rho/a) * (1+j) * xi * tanh((1+j)*xi)

    # For unit comparison, assume a = 1:
    # Z_s / rho = (1+j) * xi * tanh((1+j)*xi)

    return Z_s_normalized  # This is Z_s / (rho/a)


def test_esim_surface_impedance():
    """Test ESIM surface impedance calculation."""
    print("=" * 70)
    print("ESIM Surface Impedance Test")
    print("=" * 70)

    print("""
    ESIM computes the surface impedance Z_s, NOT R_ac/R_dc.

    Z_s = E(0) / H(0) = rho * gamma * tanh(gamma * a)

    For PEEC, the conductor resistance is:
        R = Re(Z_s) * (perimeter * length) / area
          = Re(Z_s) * perimeter / thickness

    This is different from R_ac/R_dc which compares AC to DC resistance
    for the same total current.
    """)

    # Copper parameters
    sigma = 5.8e7  # S/m
    rho = 1.0 / sigma
    mu_r = 1.0
    half_thickness = 0.001  # 1 mm

    frequencies = [100, 1000, 10000, 50000, 100000]

    print(f"\nCopper conductor: sigma = {sigma:.2e} S/m")
    print(f"Half-thickness a = {half_thickness*1000:.1f} mm")
    print()

    print(f"{'Freq [Hz]':>12} {'xi':>8} {'|Z_s| ESIM':>14} {'|Z_s| Anal':>14} {'Ratio':>10}")
    print("-" * 62)

    for freq in frequencies:
        omega = 2 * np.pi * freq
        delta = np.sqrt(2 * rho / (omega * mu_0 * mu_r))
        xi = half_thickness / delta

        # ESIM finite slab solver
        solver = ESIMFiniteSlabSolver(
            half_thickness=half_thickness,
            sigma=sigma,
            frequency=freq,
            mu_r=mu_r,
            n_nodes=200
        )
        result = solver.solve(1.0)  # H0 = 1 A/m
        Z_esim = result['Z']

        # Analytical formula for finite slab
        Z_anal_normalized = analytical_finite_slab_Zs(xi, rho)
        Z_anal = Z_anal_normalized * (rho / half_thickness)

        ratio = np.abs(Z_esim) / np.abs(Z_anal) if np.abs(Z_anal) > 1e-20 else 0

        print(f"{freq:>12.0f} {xi:>8.3f} {np.abs(Z_esim)*1e6:>14.4f} {np.abs(Z_anal)*1e6:>14.4f} {ratio:>10.4f}")

    print("\n(Z_s in micro-Ohm)")


def test_esim_vs_semi_infinite():
    """Compare finite slab ESIM with semi-infinite for thick conductors."""
    print("\n" + "=" * 70)
    print("Finite Slab vs Semi-Infinite (Thick Conductor)")
    print("=" * 70)

    sigma = 5.8e7
    rho = 1.0 / sigma
    half_thickness = 0.01  # 10 mm (thick)
    freq = 100000  # 100 kHz

    omega = 2 * np.pi * freq
    delta = np.sqrt(2 * rho / (omega * mu_0))
    xi = half_thickness / delta

    print(f"\nThick conductor: a = {half_thickness*1000:.1f} mm")
    print(f"Skin depth: delta = {delta*1000:.3f} mm")
    print(f"xi = a/delta = {xi:.2f} (thick regime)")
    print()

    # Semi-infinite analytical
    Z_semi = analytical_surface_impedance(sigma, omega, mu_0)

    # Finite slab ESIM
    solver = ESIMFiniteSlabSolver(
        half_thickness=half_thickness,
        sigma=sigma,
        frequency=freq,
        mu_r=1.0,
        n_nodes=200
    )
    result = solver.solve(1.0)
    Z_finite = result['Z']

    print(f"Semi-infinite Z_s: {Z_semi.real*1e6:.4f} + j{Z_semi.imag*1e6:.4f} micro-Ohm")
    print(f"Finite slab Z_s:   {Z_finite.real*1e6:.4f} + j{Z_finite.imag*1e6:.4f} micro-Ohm")
    print(f"Ratio |Z_finite/Z_semi|: {np.abs(Z_finite/Z_semi):.4f}")

    print("""
    For thick conductors (xi >> 1), finite slab should approach semi-infinite.
    The ratio should be close to 1.0.
    """)


def test_peec_conductor_resistance():
    """
    Demonstrate how to use ESIM Z_s for PEEC conductor resistance.
    """
    print("\n" + "=" * 70)
    print("PEEC Conductor Resistance from ESIM Z_s")
    print("=" * 70)

    print("""
    For a rectangular conductor with width w and height h:

    Cross-section area: A = w * h
    Perimeter: P = 2 * (w + h)
    Length: L

    DC resistance:
        R_dc = rho * L / A

    AC resistance using surface impedance:
        The current flows in a skin layer around the perimeter.
        Effective conducting area = P * delta (for high frequency)

        R_ac = Re(Z_s) * L * P / P = Re(Z_s) * L

        Wait, that's not right. Let me reconsider.

    Actually, for a conductor with current flowing through it:
        Power loss = Re(Z_s) * |H_surface|^2 * (surface area)
                   = Re(Z_s) * (I / P)^2 * (P * L)
                   = Re(Z_s) * I^2 * L / P

        R = Power / I^2 = Re(Z_s) * L / P

    Hmm, this gives R decreasing with P, which is wrong for DC.

    The issue: Z_s is defined per unit area of the conductor surface.
    For a conductor with current I, H_surface = I / P (Ampere's law).

    Power loss:
        P_loss = Re(Z_s) * |H_s|^2 * Area_surface
               = Re(Z_s) * (I/P)^2 * (P * L)
               = Re(Z_s) * L * I^2 / P

        R = P_loss / I^2 = Re(Z_s) * L / P

    At DC: Z_s = 0 (from tanh formula), so R -> 0? That's wrong!

    The problem: ESIM with dH/dz=0 BC gives Z_s -> 0 at DC.
    This is correct for surface impedance, but DC current flows through
    the BULK, not just the surface.

    SOLUTION:
    For conductors, we need:
        R_total = R_bulk + R_skin

    where:
        R_bulk = rho * L / A  (always present)
        R_skin = contribution from skin effect

    The correct formula:
        R_ac = R_dc * F(xi)

    where F(xi) is the AC/DC ratio function (Dowell or Kelvin).

    For ESIM approach with nonlinear materials:
        Solve the cell problem numerically to get Z_s(H0)
        Then use: R = (R_bulk + Z_s_contribution) depending on geometry

    Actually, for SIBC in PEEC:
        The surface impedance Z_s is used to modify the boundary condition,
        not directly as the resistance.

        The SIBC condition: E_t = Z_s * H_t (at conductor surface)

        In PEEC, this translates to additional voltage drop at the surface.
    """)

    # Example calculation
    sigma = 5.8e7
    rho = 1.0 / sigma
    width = 0.002  # 2 mm
    height = 0.001  # 1 mm
    length = 0.1  # 100 mm
    freq = 50000  # 50 kHz

    area = width * height
    perimeter = 2 * (width + height)
    half_thickness = min(width, height) / 2

    omega = 2 * np.pi * freq
    delta = np.sqrt(2 * rho / (omega * mu_0))
    xi = half_thickness / delta

    print(f"\nRectangular conductor:")
    print(f"  Width: {width*1000:.1f} mm")
    print(f"  Height: {height*1000:.1f} mm")
    print(f"  Length: {length*1000:.1f} mm")
    print(f"  Frequency: {freq/1000:.1f} kHz")
    print(f"  Skin depth: {delta*1000:.4f} mm")
    print(f"  xi = min_dim/(2*delta) = {xi:.3f}")
    print()

    # DC resistance
    R_dc = rho * length / area
    print(f"DC resistance: R_dc = {R_dc*1000:.6f} mOhm")

    # AC resistance using Dowell
    def dowell_ratio(xi):
        if xi < 0.01:
            return 1.0
        sh2, sn2 = np.sinh(2*xi), np.sin(2*xi)
        ch2, cs2 = np.cosh(2*xi), np.cos(2*xi)
        return xi * (sh2 + sn2) / (ch2 - cs2)

    R_ac_dowell = R_dc * dowell_ratio(xi)
    print(f"AC resistance (Dowell): R_ac = {R_ac_dowell*1000:.6f} mOhm")
    print(f"R_ac/R_dc = {dowell_ratio(xi):.4f}")

    # ESIM surface impedance
    solver = ESIMFiniteSlabSolver(
        half_thickness=half_thickness,
        sigma=sigma,
        frequency=freq,
        mu_r=1.0,
        n_nodes=200
    )
    result = solver.solve(1.0)
    Z_s = result['Z']

    print(f"\nESIM surface impedance: Z_s = {Z_s.real*1e6:.4f} + j{Z_s.imag*1e6:.4f} micro-Ohm")

    # The correct relationship between Z_s and R_ac/R_dc for planar geometry:
    # R_ac/R_dc = Re(Z_s) * a / rho = Re((gamma*a) * tanh(gamma*a))
    R_ratio_from_Zs = Z_s.real * half_thickness / rho
    print(f"R_ac/R_dc from Z_s: {R_ratio_from_Zs:.4f}")

    print("""
    Note: The R_ac/R_dc from Z_s doesn't match Dowell at small xi.
    This is because Z_s -> 0 at DC (no surface effect), while R_dc is finite.

    CONCLUSION:
    For PEEC, use Dowell's formula (or Kelvin for round) for R_ac/R_dc.
    Use ESIM Z_s for SIBC boundary conditions in high-frequency regime.
    For nonlinear materials, ESIM provides Z_s(H0) table for iteration.
    """)


def main():
    test_esim_surface_impedance()
    test_esim_vs_semi_infinite()
    test_peec_conductor_resistance()

    print("\n" + "=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    print("""
    ESIM Purpose:
    1. Compute surface impedance Z_s for SIBC boundary conditions
    2. Handle nonlinear materials where mu(H) varies
    3. Generate Z_s(H0) lookup tables for 3D solver iteration

    For PEEC R_ac/R_dc:
    1. Linear materials: Use Dowell (rectangular) or Kelvin (round) formulas
    2. Transition region: Same formulas are accurate
    3. Nonlinear materials: Use ESIM to get Z_s, then apply to geometry

    Key Insight:
    ESIM's dH/dz=0 BC models the surface impedance phenomenon.
    At DC, Z_s = 0 (no surface effect), but R_dc is still finite (bulk resistance).
    The total resistance = bulk + surface contributions.
    """)


if __name__ == '__main__':
    main()
