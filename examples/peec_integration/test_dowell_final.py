"""
FINAL correct derivation of Dowell's formula.

The key insight: Dowell's formula uses a DIFFERENT geometry than what we assumed.

Dowell considers a conductor with:
- Width w (into page)
- Height h = 2a (total thickness)
- Current I flowing in x-direction
- H-field in y-direction

The boundary conditions are:
- H(bottom surface) = 0
- H(top surface) = I/w (total current divided by width)

This is a SINGLE-SIDED excitation, not symmetric!

For SYMMETRIC excitation (current from both surfaces):
- H(bottom) = H0
- H(top) = -H0 (or use symmetry: dH/dz = 0 at center)

Let's derive both cases.
"""

import numpy as np


def dowell_formula(xi):
    """Dowell's formula for SINGLE-SIDED excitation (foil winding)."""
    if xi < 0.01:
        return 1.0 + xi**4 / 45

    sh2 = np.sinh(2*xi)
    sn2 = np.sin(2*xi)
    ch2 = np.cosh(2*xi)
    cs2 = np.cos(2*xi)

    return xi * (sh2 + sn2) / (ch2 - cs2)


def single_sided_analytical(xi, n_points=1000):
    """
    Single-sided excitation (Dowell's geometry).

    H(0) = 0 (bottom)
    H(h) = H0 = I/w (top, where current enters)

    Solution: H(z) = H0 * sinh(gamma*z) / sinh(gamma*h)

    where h = total thickness (not half-thickness!)
    and xi = h / delta (full thickness / skin depth)

    Wait - Dowell uses xi = h/delta where h is the FULL thickness.
    Let me reread Dowell's paper definition...

    Actually in standard Dowell formulation:
    - d is the conductor thickness
    - delta is skin depth
    - The formula uses d/delta

    For our half-thickness formulation:
    - a is half-thickness
    - The formula should use 2*a/delta = 2*xi_half

    Let me use the FULL thickness formulation.
    """
    h = 1.0  # Full thickness (normalized)

    if xi < 0.001:
        return 1.0

    gamma = complex(1, 1) * xi / h  # gamma such that gamma*h = (1+j)*xi

    z = np.linspace(0, h, n_points)

    # H(z) = H0 * sinh(gamma*z) / sinh(gamma*h)
    sinh_gh = np.sinh(gamma * h)
    H = np.sinh(gamma * z) / sinh_gh

    # J(z) = -dH/dz = -H0 * gamma * cosh(gamma*z) / sinh(gamma*h)
    J = gamma * np.cosh(gamma * z) / sinh_gh

    # Power = integral{rho * |J|^2} dz
    power_ac = np.trapezoid(np.abs(J)**2, z)

    # DC case: J_dc = H0 / h (uniform current)
    # power_dc = |J_dc|^2 * h = 1/h
    power_dc = 1.0 / h

    # Current: I = H(h) - H(0) = 1 - 0 = 1 = H0
    I_ac = 1.0

    # R_ac / R_dc = (power_ac / |I|^2) / (power_dc / 1)
    R_ratio = power_ac / power_dc

    return R_ratio.real


def symmetric_excitation_analytical(xi, n_points=1000):
    """
    Symmetric excitation (current from both surfaces).

    H(0) = H0 (surface)
    dH/dz(a) = 0 (center, symmetry)

    Solution: H(z) = H0 * cosh(gamma*(a-z)) / cosh(gamma*a)

    This is our ESIMFiniteSlabSolver geometry.
    """
    a = 1.0  # Half-thickness

    if xi < 0.001:
        return 1.0

    gamma = complex(1, 1) * xi / a

    z = np.linspace(0, a, n_points)

    # H(z) = H0 * cosh(gamma*(a-z)) / cosh(gamma*a)
    cosh_ga = np.cosh(gamma * a)
    H = np.cosh(gamma * (a - z)) / cosh_ga

    # J(z) = -dH/dz = H0 * gamma * sinh(gamma*(a-z)) / cosh(gamma*a)
    J = gamma * np.sinh(gamma * (a - z)) / cosh_ga

    # Power = integral{rho * |J|^2} dz
    power_ac = np.trapezoid(np.abs(J)**2, z)

    # Current into this half: I_half = H(0) - H(a) = H0 * (1 - 1/cosh(gamma*a))
    I_half = H[0] - H[-1]

    # DC case: For uniform current in half-slab
    # J_dc = I_half_dc / a, but I_half_dc = H0 (since H goes from H0 to 0)
    # Wait, this is getting confusing. Let me think differently.

    # For symmetric excitation with total conductor thickness 2a:
    # Total current = 2 * I_half
    # Total power = 2 * power_half
    # Total DC power = 2 * power_dc_half

    # Actually the simpler approach:
    # R_ac/R_dc should compare AC power dissipation to DC power dissipation
    # for the SAME total current.

    # At DC: H(z) = H0 * (a-z)/a (linear from H0 at z=0 to 0 at z=a)
    #        J_dc = H0/a (uniform)
    #        I_dc = H0 (current = H at surface)
    #        power_dc = |J_dc|^2 * a = H0^2 / a

    power_dc = 1.0 / a  # For H0 = 1

    # For AC:
    # I_ac = H(0) = H0 = 1 (current at surface)
    I_ac = 1.0  # H(0) = H0 = 1

    R_ac = power_ac / np.abs(I_ac)**2
    R_dc = power_dc / np.abs(I_ac)**2

    R_ratio = R_ac / R_dc

    return R_ratio.real


def correct_surface_impedance(xi):
    """
    Correct surface impedance formulation.

    For SYMMETRIC excitation:
    Z_s = E(0) / H(0)

    E(0) = rho * J(0) = rho * H0 * gamma * sinh(gamma*a) / cosh(gamma*a)
         = rho * H0 * gamma * tanh(gamma*a)

    Z_s = rho * gamma * tanh(gamma*a)

    Z_dc = rho * J_dc / H0 = rho * (H0/a) / H0 = rho / a

    R_ratio = Re(Z_s) / Z_dc = Re(gamma * a * tanh(gamma*a))
    """
    if xi < 0.001:
        return 1.0

    gamma_a = complex(1, 1) * xi
    result = (gamma_a * np.tanh(gamma_a)).real

    return result


def main():
    print("=" * 80)
    print("FINAL Dowell Formula Derivation")
    print("=" * 80)

    print("""
    Key Insight: Dowell's formula is for SINGLE-SIDED excitation (foil windings)
    where current enters from one side only.

    Geometry Comparison:

    1. SINGLE-SIDED (Dowell):
       - H(0) = 0 (bottom)
       - H(h) = H0 (top, current side)
       - h = full conductor thickness
       - xi = h / delta

    2. SYMMETRIC (ESIMFiniteSlabSolver):
       - H(0) = H0 (surface)
       - dH/dz(a) = 0 (center)
       - a = half-thickness
       - xi = a / delta

    For symmetric case with half-thickness a:
    - Effective xi_Dowell = 2 * xi_symmetric (full thickness / delta)
    """)

    xi_values = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]

    print(f"\n{'xi':>6} {'Dowell':>10} {'Single':>10} {'Symm':>10} {'SurfImp':>10}")
    print("-" * 50)

    for xi in xi_values:
        r_dowell = dowell_formula(xi)
        r_single = single_sided_analytical(xi)
        r_symm = symmetric_excitation_analytical(xi)
        r_surf = correct_surface_impedance(xi)

        print(f"{xi:>6.2f} {r_dowell:>10.4f} {r_single:>10.4f} {r_symm:>10.4f} {r_surf:>10.4f}")

    print("\n" + "=" * 80)
    print("Matching Dowell with Symmetric Geometry")
    print("=" * 80)

    print("""
    For symmetric geometry (half-thickness a, xi = a/delta):
    To match Dowell's formula, we need to use xi_full = 2*xi

    Let's verify: Dowell(2*xi) should match symmetric geometry result
    """)

    print(f"\n{'xi_half':>8} {'Dowell(2xi)':>12} {'Symmetric':>12} {'Ratio':>10}")
    print("-" * 46)

    for xi_half in xi_values:
        xi_full = 2 * xi_half
        r_dowell_2xi = dowell_formula(xi_full)
        r_symm = symmetric_excitation_analytical(xi_half)

        ratio = r_symm / r_dowell_2xi if r_dowell_2xi > 0 else 0

        print(f"{xi_half:>8.2f} {r_dowell_2xi:>12.4f} {r_symm:>12.4f} {ratio:>10.4f}")


if __name__ == '__main__':
    main()
