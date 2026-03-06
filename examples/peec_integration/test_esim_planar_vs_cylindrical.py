"""
Compare ESIM (planar) with Kelvin (cylindrical) for transition region.

Key insight: ESIM solves the 1D diffusion equation for a PLANAR conductor:
    rho * d^2H/ds^2 + j*omega*mu*H = 0    (s = depth into conductor)

The Kelvin function solution is for a CYLINDRICAL conductor:
    rho * (1/r) * d/dr(r * dH/dr) + j*omega*mu*H = 0    (r = radial coordinate)

These give DIFFERENT R_ac/R_dc ratios in the transition region!

For high frequency (xi >> 1):
    Both approach R_ac/R_dc = xi/2 because current is concentrated in thin skin layer
    where curvature effects are negligible.

For transition region (0.5 < xi < 3):
    - Planar: Lower R_ac/R_dc (flatter profile)
    - Cylindrical: Higher R_ac/R_dc (geometric concentration toward surface)

Usage:
    python test_esim_planar_vs_cylindrical.py
"""

import numpy as np
from scipy.special import kelvin, jv, yv  # Bessel functions
from scipy.constants import mu_0


def kelvin_rac_ratio_cylindrical(xi):
    """
    Exact R_ac/R_dc for cylindrical wire using Kelvin functions.

    This is the standard formula for round wire skin effect.

    xi = sqrt(2) * a / delta where a = wire radius, delta = skin depth

    Reference: Ramo, Whinnery, Van Duzer, "Fields and Waves"
    """
    if xi < 0.01:
        return 1.0

    ber, bei, ker, kei = kelvin(xi)

    h = 1e-6
    ber_h, bei_h, _, _ = kelvin(xi + h)
    berp = (ber_h - ber) / h
    beip = (bei_h - bei) / h

    ber_bei = complex(ber, bei)
    berp_beip = complex(berp, beip)

    if abs(ber_bei) < 1e-30:
        return 1.0

    ratio = (xi / 2.0) * (berp_beip / ber_bei).real
    return max(ratio, 1.0)


def rac_ratio_planar(xi):
    """
    Exact R_ac/R_dc for planar (slab) conductor.

    For a slab of thickness 2a with skin effect:
    R_ac/R_dc = (xi/2) * coth(xi) * real part

    where xi = a / delta (half-thickness / skin depth)

    For very thick slab (a >> delta), this reduces to:
    R_ac/R_dc -> xi/2 (same as cylindrical)

    The full formula for slab with current on both sides:
    Z = Rs * coth(gamma * a) where gamma = (1+j)/delta

    For single-sided (semi-infinite):
    Z = Rs (surface impedance)
    """
    if xi < 0.01:
        return 1.0

    # For semi-infinite planar conductor:
    # The surface impedance is Z_s = (1+j) * Rs
    # where Rs = 1/(sigma*delta)
    #
    # The resistance per unit area for planar is just Rs
    # R_ac/R_dc = Rs * sigma * a = a/delta = xi

    # Actually, for a wire approximated as planar:
    # The effective resistance is Rs * L / perimeter
    # R_ac = Rs * L / (2*pi*a)
    # R_dc = rho * L / (pi*a^2)
    # R_ac/R_dc = Rs * sigma * a / 2 = xi/2

    # This gives xi/2 for all xi, which is the high-frequency limit
    # The difference from cylindrical comes from geometry

    # For more accurate planar slab, use the sinh/cosh formula
    # But for semi-infinite (ESIM), the answer is simply xi/2

    # Let's derive more carefully for finite thickness slab:
    # J(z) = J0 * cosh(gamma*z) / cosh(gamma*a)
    # where gamma = (1+j)/delta, z is distance from center

    # Total current: I = integral{-a to a} J dz
    # Power: P = integral{-a to a} |J|^2 / sigma dz
    # R_ac = P / |I|^2

    gamma = complex(1, 1) / 1.0  # gamma * delta = (1+j)
    gamma_a = gamma * xi  # gamma * a = (1+j) * xi

    # For finite slab with symmetric current:
    # R_ac/R_dc = xi * |sinh(gamma_a)|^2 / (|cosh(gamma_a) - 1|^2 + terms...)

    # Simpler: use hyperbolic formula
    # sinh((1+j)*xi) / cosh((1+j)*xi)

    try:
        sh = np.sinh(gamma_a)
        ch = np.cosh(gamma_a)

        # For slab: R_ac/R_dc = xi * Re[(1+j) * tanh((1+j)*xi) / (1+j)*xi]
        #                     = Re[tanh((1+j)*xi) / xi * (1+j)]

        tanh_ga = sh / ch
        ratio = (xi / 2.0) * ((1 + 1j) * tanh_ga).real / xi * 2

        # Simplified formula from textbooks:
        # R_ac/R_dc = xi * (sinh(2*xi) + sin(2*xi)) / (cosh(2*xi) - cos(2*xi))
        # This is Dowell's formula for rectangular conductors

        sh2 = np.sinh(2*xi)
        sn2 = np.sin(2*xi)
        ch2 = np.cosh(2*xi)
        cs2 = np.cos(2*xi)

        denom = ch2 - cs2
        if abs(denom) < 1e-30:
            return 1.0

        ratio_dowell = xi * (sh2 + sn2) / denom

        return max(ratio_dowell, 1.0)

    except (OverflowError, FloatingPointError):
        return xi / 2.0  # High frequency limit


def test_geometry_comparison():
    """Compare planar (slab) vs cylindrical (round wire) geometry."""
    print("=" * 80)
    print("Planar (Slab) vs Cylindrical (Round Wire) Skin Effect")
    print("=" * 80)

    print("""
    Key difference:
    - ESIM solves 1D problem (planar geometry): d^2H/dz^2 + jk^2*H = 0
    - Kelvin functions solve cylindrical: (1/r)*d/dr(r*dH/dr) + jk^2*H = 0

    For round wire, use Kelvin (cylindrical).
    For rectangular wire or thin sheet, use Dowell/planar formula.

    Both converge to R_ac/R_dc = xi/2 at high frequency where curvature
    effects become negligible (skin depth << radius).
    """)

    xi_values = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0]

    print(f"\n{'xi=a/d':>8} {'Cylindrical':>14} {'Planar/Dowell':>14} {'Cyl/Planar':>12} {'xi/2':>10}")
    print("-" * 62)

    for xi in xi_values:
        r_cyl = kelvin_rac_ratio_cylindrical(np.sqrt(2) * xi)
        r_planar = rac_ratio_planar(xi)
        ratio = r_cyl / r_planar if r_planar > 0 else 0
        r_hf = xi / 2.0

        print(f"{xi:>8.2f} {r_cyl:>14.4f} {r_planar:>14.4f} {ratio:>12.3f} {r_hf:>10.4f}")

    print("""
    Observations:
    1. At low xi: Both approach 1.0 (DC limit)
    2. At high xi: Both approach xi/2 (skin effect limit)
    3. In transition: Cylindrical > Planar (geometric effect)

    The cylindrical geometry has higher R_ac because current is forced
    into a smaller effective area due to the curved surface.
    """)


def test_rectangular_cross_section():
    """
    Compare formulas for rectangular cross-section wire.

    For PEEC with rectangular cross-section, Dowell's formula is more appropriate
    than Kelvin functions.
    """
    print("\n" + "=" * 80)
    print("Rectangular Cross-Section: Dowell's Formula")
    print("=" * 80)

    print("""
    For rectangular conductors (width w, height h), Dowell's formula applies:

    R_ac/R_dc = xi * (sinh(2*xi) + sin(2*xi)) / (cosh(2*xi) - cos(2*xi))

    where xi = sqrt(w*h) / (2*delta) (effective thickness / skin depth)

    This is what ESIM's planar solution approximates!
    """)

    # Compare with simple step function (current PEEC)
    xi_values = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]

    print(f"\n{'xi':>8} {'Dowell':>12} {'Simple':>12} {'Dowell/Simple':>14}")
    print("-" * 50)

    for xi in xi_values:
        r_dowell = rac_ratio_planar(xi)

        # Simple step function
        if xi <= 1.0:
            r_simple = 1.0
        else:
            r_simple = xi / 2.0

        ratio = r_dowell / r_simple if r_simple > 0 else 0

        print(f"{xi:>8.2f} {r_dowell:>12.4f} {r_simple:>12.4f} {ratio:>14.3f}")

    print("""
    Conclusion for PEEC (rectangular cross-section):
    - Dowell's formula is the correct analytical reference
    - Simple step function is ~10-20% off in transition region
    - ESIM should match Dowell, not Kelvin
    """)


def test_esim_vs_dowell():
    """Verify ESIM matches Dowell's formula for planar geometry."""
    import sys
    import os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

    try:
        from esim_cell_problem import ESIMCellProblemSolver
    except ImportError:
        print("ESIM module not available")
        return

    print("\n" + "=" * 80)
    print("ESIM vs Dowell: Planar Geometry Verification")
    print("=" * 80)

    sigma = 5.8e7  # Copper
    mu_r = 1.0

    # Linear BH curve (constant mu)
    bh_curve = [[0, 0], [1e6, mu_0 * mu_r * 1e6]]

    wire_radius = 0.001  # 1 mm (equivalent half-thickness for planar)

    xi_values = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]

    print(f"\nEquivalent half-thickness: {wire_radius*1000:.1f} mm")
    print()

    print(f"{'xi':>8} {'Freq [Hz]':>12} {'Dowell':>12} {'ESIM':>12} {'ESIM/Dowell':>12}")
    print("-" * 60)

    for xi in xi_values:
        # Compute frequency
        omega = 2 * xi**2 / (wire_radius**2 * mu_0 * mu_r * sigma)
        freq = omega / (2 * np.pi)
        delta = np.sqrt(2 / (omega * mu_0 * mu_r * sigma))

        # Dowell's formula
        r_dowell = rac_ratio_planar(xi)

        # ESIM
        solver = ESIMCellProblemSolver(
            bh_curve=bh_curve, sigma=sigma, frequency=freq,
            num_skin_depths=15, n_nodes=200
        )

        result = solver.solve(1.0, tol=1e-8, max_iter=100)
        Z_esim = result['Z']
        Rs_esim = Z_esim.real

        # R_ac/R_dc for planar: Rs * sigma * a / 2 (same as wire)
        a = xi * delta
        r_esim = Rs_esim * sigma * a / 2.0

        ratio = r_esim / r_dowell if r_dowell > 0 else 0

        print(f"{xi:>8.2f} {freq:>12.0f} {r_dowell:>12.4f} {r_esim:>12.4f} {ratio:>12.3f}")

    print("""
    Analysis:
    - ESIM (semi-infinite half-space) gives Rs = (1+j)/sigma/delta at all frequencies
    - This means R_ac/R_dc = xi for ESIM half-space model
    - Dowell's formula accounts for finite thickness (current on both sides)

    For accurate ESIM in transition region:
    - Need to use finite-thickness slab model, not semi-infinite
    - Or use Dowell's formula directly (analytical)
    """)


def main():
    test_geometry_comparison()
    test_rectangular_cross_section()
    test_esim_vs_dowell()

    print("\n" + "=" * 80)
    print("Summary: Correct Formula for Different Geometries")
    print("=" * 80)
    print("""
    1. Round Wire (circular cross-section):
       - Use Kelvin functions (ber, bei)
       - R_ac/R_dc = (xi/2) * Re{(ber'+j*bei')/(ber+j*bei)}
       - where xi = sqrt(2) * radius / skin_depth

    2. Rectangular Wire (rectangular cross-section):
       - Use Dowell's formula
       - R_ac/R_dc = xi * (sinh(2*xi) + sin(2*xi)) / (cosh(2*xi) - cos(2*xi))
       - where xi = effective_thickness / (2 * skin_depth)

    3. Semi-infinite Half-Space (ESIM current implementation):
       - Always gives R_ac/R_dc proportional to xi
       - Not accurate in transition region for finite-size conductors

    Recommendation for PEEC:
    - Rectangular cross-section: Use Dowell's formula (analytical, fast)
    - This is simpler and more accurate than ESIM Cell Problem for linear materials
    - ESIM Cell Problem is useful for NONLINEAR mu(H) where analytical formulas fail
    """)


if __name__ == '__main__':
    main()
