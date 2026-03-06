"""
Correct derivation of R_ac/R_dc for finite slab.

The issue: We have two different definitions of "R_ac/R_dc"

1. Surface impedance approach (relates E to H at surface):
   Z_s = E(0) / H(0) = rho * gamma * tanh(gamma*a)
   Z_dc_surf = rho / a
   R_ratio_surf = Re(Z_s) / Z_dc_surf

2. Power loss approach (actual resistance in circuit):
   P = Re(Z) * |I|^2
   R = P / |I|^2
   R_dc = rho * L / A = rho / a (for unit length, unit width)
   R_ac = actual AC resistance

For circuits, we need approach 2 (power-based).

The correct formula for power-based R_ac/R_dc:

For a slab with current I entering from top and bottom surfaces:
- H(0) = I/2 (half the current on each surface)
- Total power = 2 * integral{rho * |J|^2} dz from 0 to a

Let me derive this properly.
"""

import numpy as np


def dowell_formula(xi):
    """Dowell's formula - the CORRECT R_ac/R_dc for planar conductor."""
    if xi < 0.01:
        return 1.0 + xi**4 / 45

    sh2 = np.sinh(2*xi)
    sn2 = np.sin(2*xi)
    ch2 = np.cosh(2*xi)
    cs2 = np.cos(2*xi)

    return xi * (sh2 + sn2) / (ch2 - cs2)


def analytical_H_field(z, a, gamma):
    """Analytical H(z) for finite slab with symmetry BC."""
    return np.cosh(gamma * (a - z)) / np.cosh(gamma * a)


def analytical_J_field(z, a, gamma):
    """Analytical J(z) = -dH/dz."""
    return gamma * np.sinh(gamma * (a - z)) / np.cosh(gamma * a)


def power_loss_integral(xi, n_points=1000):
    """
    Compute R_ac/R_dc using power loss integral.

    P_ac = integral{rho * |J|^2} dz from 0 to a
    P_dc = rho * |I|^2 / a  where I = integral{J_dc} dz = J_dc * a

    For H0 = 1: I_dc = H0 = 1 (current per unit width)

    At DC: J_dc = I/a = H0/a, P_dc = rho * (H0/a)^2 * a = rho * H0^2 / a

    For AC: J(z) = H0 * gamma * sinh(gamma*(a-z)) / cosh(gamma*a)
            I_ac = H0 - H(a) = H0 * (1 - 1/cosh(gamma*a))
            P_ac = rho * integral{|J|^2} dz

    R_ac/R_dc = P_ac / |I_ac|^2  /  (P_dc / |I_dc|^2)
    """
    a = 1.0  # Normalized half-thickness

    if xi < 0.001:
        return 1.0

    gamma = complex(1, 1) * xi / a  # gamma = (1+j) * xi / a

    # Numerical integration
    z = np.linspace(0, a, n_points)
    dz = z[1] - z[0]

    # H(z) = cosh(gamma*(a-z)) / cosh(gamma*a)
    H = analytical_H_field(z, a, gamma)

    # J(z) = gamma * sinh(gamma*(a-z)) / cosh(gamma*a)
    J = analytical_J_field(z, a, gamma)

    # Power integral: integral{|J|^2} dz
    power_ac = np.trapz(np.abs(J)**2, z)

    # Current: I = H(0) - H(a)
    I_ac = H[0] - H[-1]

    # DC: J_dc = 1/a, I_dc = 1
    # power_dc = 1/a (for rho=1, a=1)
    power_dc = 1.0 / a

    # R_ac = power_ac / |I_ac|^2 (for rho=1)
    # R_dc = power_dc / |I_dc|^2 = 1/a / 1 = 1/a

    R_ac_over_Iac2 = power_ac / (np.abs(I_ac)**2)
    R_dc_over_Idc2 = power_dc / 1.0

    R_ratio = R_ac_over_Iac2 / R_dc_over_Idc2

    return R_ratio.real


def surface_impedance_ratio(xi):
    """
    R_ac/R_dc from surface impedance.

    Z_s = rho * gamma * tanh(gamma*a)
    Z_dc = rho / a

    R_ratio = Re(Z_s) / Z_dc = Re(gamma*a * tanh(gamma*a))
    """
    if xi < 0.001:
        return 1.0

    gamma_a = complex(1, 1) * xi
    Z_s_normalized = gamma_a * np.tanh(gamma_a)

    return Z_s_normalized.real


def derived_formula(xi):
    """
    Correct formula derived from analytical integration.

    From Dowell's paper, the analytical integration gives:
    R_ac/R_dc = xi * (sinh(2xi) + sin(2xi)) / (cosh(2xi) - cos(2xi))

    Let's verify this equals Re(gamma*a * tanh(gamma*a))
    """
    if xi < 0.01:
        return 1.0 + xi**4 / 45

    # Using hyperbolic/trig identities...
    # Let z = (1+j)*xi
    # tanh(z) = sinh(z)/cosh(z)
    # |cosh(z)|^2 = cosh^2(Re(z))*cos^2(Im(z)) + sinh^2(Re(z))*sin^2(Im(z))
    #             = (cosh(2*Re(z)) + cos(2*Im(z))) / 2  for z = x + jy
    # For z = (1+j)*xi, Re(z) = Im(z) = xi
    # |cosh(z)|^2 = (cosh(2xi) + cos(2xi)) / 2

    # Similarly:
    # Re(z * tanh(z)) = ?

    # Let me just compute numerically
    z = (1+1j) * xi
    result = (z * np.tanh(z)).real

    return result


def main():
    print("=" * 70)
    print("Correct R_ac/R_dc Derivation for Finite Slab")
    print("=" * 70)

    print("""
    Comparing three approaches:

    1. Dowell's formula: xi * (sinh(2xi) + sin(2xi)) / (cosh(2xi) - cos(2xi))
       This is the industry standard formula.

    2. Surface impedance: Re(gamma*a * tanh(gamma*a))
       This is Z_s / Z_dc where Z_s = rho*gamma*tanh(gamma*a), Z_dc = rho/a

    3. Power loss integral: P_ac/|I_ac|^2 / (P_dc/|I_dc|^2)
       This uses the actual power dissipation.

    Key insight: Dowell's formula is derived from SURFACE IMPEDANCE approach,
    not from power loss approach!
    """)

    xi_values = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]

    print(f"\n{'xi':>6} {'Dowell':>10} {'SurfImp':>10} {'Power':>10} {'Surf/Dow':>10} {'Pow/Dow':>10}")
    print("-" * 62)

    for xi in xi_values:
        r_dowell = dowell_formula(xi)
        r_surf = surface_impedance_ratio(xi)
        r_power = power_loss_integral(xi)

        ratio_surf = r_surf / r_dowell if r_dowell > 0 else 0
        ratio_power = r_power / r_dowell if r_dowell > 0 else 0

        print(f"{xi:>6.2f} {r_dowell:>10.4f} {r_surf:>10.4f} {r_power:>10.4f} "
              f"{ratio_surf:>10.4f} {ratio_power:>10.4f}")

    print("\n" + "=" * 70)
    print("Analysis")
    print("=" * 70)

    print("""
    Observation:
    - Surface impedance Re(gamma*a * tanh(gamma*a)) does NOT equal Dowell's formula!
    - There is a factor difference, especially at small xi.

    The issue is in the DEFINITION:
    - Dowell's formula is for R_ac/R_dc where the "resistance" is defined
      in terms of the TOTAL current flowing through the conductor.

    Let me reconsider the derivation...
    """)


if __name__ == '__main__':
    main()
