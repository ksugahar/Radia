"""
Correct symmetric R_ac/R_dc formula.

The issue: Re(gamma*a * tanh(gamma*a)) gives values < 1.0 for small xi.
This is physically impossible since R_ac >= R_dc.

The problem is in the DEFINITION of what "R_ac/R_dc" means for symmetric geometry.

Let me reconsider the physics:

For SYMMETRIC geometry (conductor with current from both surfaces):
- Total thickness 2a
- Current I enters from both top and bottom
- By symmetry, current in top half = current in bottom half = I/2
- H field at both surfaces = I/(2w) where w is width

The resistance of this structure:
- At DC: R_dc = rho * L / A = rho * L / (2a * w)
  For unit length L=1 and unit width w=1: R_dc = rho / (2a)

- At AC: Resistance from surface impedance
  Z_s = E(0) / H(0) = rho * gamma * tanh(gamma*a)

  But this Z_s is for ONE surface. The total impedance of the symmetric
  conductor has two surfaces in PARALLEL (top and bottom).

  Z_total = Z_s / 2 (parallel combination of two surfaces)

  Z_dc = rho / (2a) (for the full conductor)

  R_ac / R_dc = Re(Z_total) / Z_dc
              = Re(Z_s/2) / (rho/(2a))
              = Re(rho * gamma * tanh(gamma*a) / 2) * 2a / rho
              = Re(gamma * a * tanh(gamma*a))

Wait, that's the same formula! Let me think again...

Actually, the issue is that at LOW frequency:
tanh(gamma*a) ~ gamma*a (for small gamma*a)
gamma*a * tanh(gamma*a) ~ (gamma*a)^2 = 2j * xi^2

This is purely imaginary for small xi! The real part is ~0.

This means the surface impedance formulation gives very small resistance
at low frequency for symmetric geometry. But this contradicts the DC limit.

THE KEY INSIGHT:
At DC (omega=0), gamma = 0, and tanh(0) = 0.
So Z_s = 0 at DC! This is correct because there's no potential drop
at the surface when current is uniformly distributed.

The DC resistance comes from the BULK of the conductor, not the surface.

So R_ac/R_dc using surface impedance is:
- At low frequency: ~ (gamma*a)^2 ~ small
- At high frequency: ~ gamma*a ~ xi

The surface impedance ratio goes to ZERO at DC, not 1.0!

This is because at DC, the "surface impedance" concept doesn't apply.
Current flows uniformly through the conductor.

For practical R_ac/R_dc, we need to add the DC resistance contribution:

R_ac = R_skin + R_dc  (simplified model)

Or more correctly, use the TOTAL impedance including bulk effects.
"""

import numpy as np


def symmetric_surface_impedance(xi):
    """
    Pure surface impedance: Re(gamma*a * tanh(gamma*a))
    This goes to 0 at DC, which is correct for surface impedance!
    """
    if xi < 0.001:
        return 0.0  # Correct DC limit for surface impedance

    gamma_a = complex(1, 1) * xi
    return (gamma_a * np.tanh(gamma_a)).real


def symmetric_total_resistance(xi):
    """
    CORRECT R_ac/R_dc for symmetric geometry.

    The total resistance has two components:
    1. Skin effect resistance (from non-uniform current)
    2. Bulk DC resistance

    At DC: Only bulk resistance contributes
    At high freq: Skin resistance dominates

    For a conductor with symmetric excitation:
    R_ac/R_dc = max(1.0, surface_impedance_ratio)

    Actually, the correct formula uses the POWER dissipation method,
    not the surface impedance method!
    """
    # Correct approach: power-based with proper normalization
    return power_based_Rac_Rdc(xi)


def power_based_Rac_Rdc(xi, n_points=2000):
    """
    Power-based R_ac/R_dc with CORRECT normalization.

    P_ac = integral{rho * |J|^2} dz
    P_dc = rho * |J_dc|^2 * a = rho * |I|^2 / a

    R_ac = P_ac / |I|^2
    R_dc = P_dc / |I|^2 = rho / a

    R_ac/R_dc = (P_ac / |I|^2) / (rho / a) = P_ac * a / (rho * |I|^2)

    For rho = 1, a = 1:
    R_ac/R_dc = P_ac / |I|^2
    """
    a = 1.0

    if xi < 0.001:
        return 1.0

    gamma = complex(1, 1) * xi / a

    z = np.linspace(0, a, n_points)

    # H(z) = cosh(gamma*(a-z)) / cosh(gamma*a)  (normalized to H0=1)
    cosh_ga = np.cosh(gamma * a)
    H = np.cosh(gamma * (a - z)) / cosh_ga

    # J(z) = gamma * sinh(gamma*(a-z)) / cosh(gamma*a)
    J = gamma * np.sinh(gamma * (a - z)) / cosh_ga

    # Current: I = H(0) - H(a) = 1 - 1/cosh(gamma*a)
    # But wait! For R_ac/R_dc, we should normalize to same current at DC and AC.

    # At DC: H(z) = H0 * (1 - z/a) (linear), J = H0/a = const
    #        I_dc = H(0) - H(a) = H0 - 0 = H0 = 1

    # At AC: I_ac = H(0) - H(a) = 1 - 1/cosh(gamma*a)

    I_ac = 1.0 - 1.0/cosh_ga

    # Power (per unit rho)
    power_ac = np.trapezoid(np.abs(J)**2, z)

    # DC power: J_dc = I_dc/a = 1/a, power_dc = |J_dc|^2 * a = 1/a
    power_dc = 1.0 / a

    # R_ac/R_dc at same current would be:
    # R = P / I^2
    # R_ac/R_dc = (P_ac/|I_ac|^2) / (P_dc/|I_dc|^2)
    #           = P_ac / |I_ac|^2 * |I_dc|^2 / P_dc
    #           = P_ac / |I_ac|^2 * 1 / (1/a)
    #           = P_ac * a / |I_ac|^2

    R_ratio = power_ac * a / (np.abs(I_ac)**2)

    return R_ratio.real


def compare_formulas():
    """Compare different R_ac/R_dc formulations."""
    print("=" * 70)
    print("Symmetric Geometry R_ac/R_dc Comparison")
    print("=" * 70)

    print("""
    Geometry: Slab with half-thickness a, symmetric excitation
    - H(0) = H0 (surface)
    - dH/dz(a) = 0 (center, symmetry)

    Formulas compared:
    1. Surface impedance: Re(gamma*a * tanh(gamma*a))
       - Goes to 0 at DC (correct for surface Z)
       - Used for high-frequency regime

    2. Power-based: P_ac * a / |I_ac|^2
       - Approaches 1.0 at DC
       - Includes bulk resistance contribution
    """)

    xi_values = [0.01, 0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]

    print(f"\n{'xi':>8} {'SurfImp':>12} {'Power':>12}")
    print("-" * 36)

    for xi in xi_values:
        r_surf = symmetric_surface_impedance(xi)
        r_power = power_based_Rac_Rdc(xi)

        print(f"{xi:>8.2f} {r_surf:>12.6f} {r_power:>12.4f}")


def correct_symmetric_formula():
    """
    Derive the CORRECT analytical formula for symmetric R_ac/R_dc.

    From power integral:
    P_ac = |gamma|^2 * integral{|sinh(gamma*(a-z))|^2} dz / |cosh(gamma*a)|^2

    Let u = gamma*(a-z), du = -gamma*dz
    = |gamma|^2 / |cosh(gamma*a)|^2 * integral{|sinh(u)|^2} du / |gamma|

    |sinh(u)|^2 = (cosh(2*Re(u)) - cos(2*Im(u))) / 2

    For u = gamma*(a-z) where gamma = (1+j)*xi/a:
    Re(u) = (xi/a)*(a-z) = xi*(1 - z/a)
    Im(u) = (xi/a)*(a-z) = xi*(1 - z/a)

    This gets complicated. Let's use numerical verification.
    """
    print("\n" + "=" * 70)
    print("Numerical Verification of Correct Formula")
    print("=" * 70)

    # At high xi, both should approach xi
    # At low xi, power-based should approach 1.0

    print("\nHigh xi limit (should approach xi):")
    for xi in [5.0, 10.0, 20.0]:
        r = power_based_Rac_Rdc(xi)
        print(f"  xi = {xi}: R_ac/R_dc = {r:.4f}, ratio to xi = {r/xi:.4f}")

    print("\nLow xi limit (should approach 1.0):")
    for xi in [0.1, 0.05, 0.01]:
        r = power_based_Rac_Rdc(xi)
        print(f"  xi = {xi}: R_ac/R_dc = {r:.4f}")


if __name__ == '__main__':
    compare_formulas()
    correct_symmetric_formula()
