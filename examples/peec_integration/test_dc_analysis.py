"""
Analyze the DC limit issue for symmetric geometry.

At DC (xi -> 0), we get R_ac/R_dc = 4/3 instead of 1.0.
This suggests a normalization problem.

Let's trace through the DC calculation step by step.
"""

import numpy as np


def analyze_dc_limit():
    """Trace through DC calculation."""
    print("=" * 70)
    print("DC Limit Analysis for Symmetric Geometry")
    print("=" * 70)

    a = 1.0  # Half-thickness
    n_points = 100
    z = np.linspace(0, a, n_points)

    print("""
    Geometry:
    - Half-thickness a = 1.0
    - Surface at z = 0 with H(0) = H0 = 1
    - Center at z = a with dH/dz = 0 (symmetry)

    At DC, the field distribution is:
    H(z) = H0 * (1 - z/a)  [linear from H0 to 0]

    But wait! This is WRONG for symmetric geometry.
    With dH/dz(a) = 0, H(a) is NOT necessarily 0.
    """)

    # DC solution with dH/dz(a) = 0:
    # d^2H/dz^2 = 0 (no omega term)
    # Solution: H = A*z + B
    # BC1: H(0) = H0 => B = H0
    # BC2: dH/dz(a) = A = 0
    # So: H(z) = H0 (constant!)

    print("""
    DC SOLUTION with symmetric BC:
    d^2H/dz^2 = 0 => H = A*z + B
    BC: H(0) = H0 => B = H0
    BC: dH/dz(a) = 0 => A = 0
    => H(z) = H0 (constant everywhere!)

    This means:
    - J = -dH/dz = 0 (NO current flows!)
    - I = H(0) - H(a) = H0 - H0 = 0

    This is the issue! At DC, the symmetric geometry with dH/dz(a)=0
    has ZERO current flow!

    This BC models a conductor that is INSULATED at the center,
    not a conductor with current flowing through.
    """)

    print("=" * 70)
    print("CORRECT DC Model for R_ac/R_dc")
    print("=" * 70)

    print("""
    For R_ac/R_dc, we need a geometry where current FLOWS at DC.

    Option 1: Single-sided (Dowell's geometry)
    - H(0) = 0, H(h) = H0
    - Current I = H0 (flows through full thickness)

    Option 2: Symmetric with fixed total current
    - H(0) = H0, H(a) = 0
    - Current I = H0 (flows through half-thickness)
    - BC: H(a) = 0 instead of dH/dz(a) = 0

    Let's compute with Option 2 BC.
    """)


def symmetric_fixed_endpoint(xi, n_points=2000):
    """
    Symmetric geometry with H(a) = 0 (not dH/dz(a) = 0).

    This models current flowing from surface to center,
    where the center is connected to a ground/return path.

    PDE: rho * d^2H/dz^2 + j*omega*mu*H = 0
    BC: H(0) = H0, H(a) = 0

    Solution: H(z) = H0 * sinh(gamma*(a-z)) / sinh(gamma*a)
    """
    a = 1.0

    if xi < 0.001:
        # DC: H(z) = H0 * (a-z)/a
        # J = H0/a
        # I = H0
        # P = |H0/a|^2 * a = H0^2/a
        # R = P/I^2 = 1/a
        # R_dc = 1/a
        # R_ac/R_dc = 1.0
        return 1.0

    gamma = complex(1, 1) * xi / a
    z = np.linspace(0, a, n_points)

    # H(z) = sinh(gamma*(a-z)) / sinh(gamma*a)
    sinh_ga = np.sinh(gamma * a)
    H = np.sinh(gamma * (a - z)) / sinh_ga

    # J(z) = -dH/dz = gamma * cosh(gamma*(a-z)) / sinh(gamma*a)
    J = gamma * np.cosh(gamma * (a - z)) / sinh_ga

    # Current: I = H(0) - H(a) = 1 - 0 = 1
    I = 1.0

    # Power
    power_ac = np.trapezoid(np.abs(J)**2, z)

    # DC power: J_dc = 1/a, power_dc = 1/a
    power_dc = 1.0 / a

    R_ratio = (power_ac / np.abs(I)**2) / power_dc

    return R_ratio.real


def dowell_formula(xi):
    """Dowell's formula."""
    if xi < 0.01:
        return 1.0 + xi**4 / 45

    sh2 = np.sinh(2*xi)
    sn2 = np.sin(2*xi)
    ch2 = np.cosh(2*xi)
    cs2 = np.cos(2*xi)

    return xi * (sh2 + sn2) / (ch2 - cs2)


def surface_impedance_coth(xi):
    """Surface impedance with coth (single-sided)."""
    if xi < 0.001:
        return 1.0

    gamma_a = complex(1, 1) * xi
    return (gamma_a / np.tanh(gamma_a)).real  # gamma*a * coth(gamma*a)


def main():
    analyze_dc_limit()

    print("\n" + "=" * 70)
    print("Comparison: Symmetric H(a)=0 vs Dowell")
    print("=" * 70)

    print("""
    Both geometries have H(0) = H0, H(endpoint) = 0.
    They should give the SAME R_ac/R_dc formula!

    - Dowell: H(0) = 0, H(h) = H0 (h = full thickness)
    - Symmetric H(a)=0: H(0) = H0, H(a) = 0 (a = half thickness)

    These are the SAME geometry with reversed z-direction!
    """)

    xi_values = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]

    print(f"\n{'xi':>8} {'Dowell':>12} {'Symm H(a)=0':>12} {'Surf coth':>12}")
    print("-" * 48)

    for xi in xi_values:
        r_dowell = dowell_formula(xi)
        r_symm = symmetric_fixed_endpoint(xi)
        r_coth = surface_impedance_coth(xi)

        print(f"{xi:>8.2f} {r_dowell:>12.4f} {r_symm:>12.4f} {r_coth:>12.4f}")

    print("\n" + "=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    print("""
    The symmetric geometry with dH/dz(a) = 0 models a conductor where
    current enters from the surface but CANNOT exit at the center.
    This is NOT the correct model for R_ac/R_dc!

    CORRECT GEOMETRIES for R_ac/R_dc:

    1. Single-sided (Dowell):
       - H(0) = 0, H(h) = H0
       - R_ac/R_dc = xi * (sinh(2xi)+sin(2xi)) / (cosh(2xi)-cos(2xi))
       - This is the INDUSTRY STANDARD formula

    2. Symmetric with H(a) = 0:
       - H(0) = H0, H(a) = 0
       - SAME formula as Dowell (reversed z-axis)

    The ESIMFiniteSlabSolver with dH/dz(a) = 0 BC is for a DIFFERENT problem:
    - Models surface impedance of a thick conductor
    - NOT for R_ac/R_dc calculation

    For R_ac/R_dc in PEEC, use Dowell's formula directly!
    """)


if __name__ == '__main__':
    main()
