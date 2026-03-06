"""
Derive the correct R_ac/R_dc formula for SYMMETRIC excitation.

For symmetric excitation geometry:
- Conductor thickness: 2a
- Surface at z = 0 with H = H0
- Center at z = a with dH/dz = 0 (symmetry)

Current enters from BOTH surfaces, so total current I_total = 2 * I_half

The solution is:
    H(z) = H0 * cosh(gamma*(a-z)) / cosh(gamma*a)

where gamma = (1+j) / delta

Let's derive R_ac/R_dc for this geometry CORRECTLY.
"""

import numpy as np


def symmetric_formula_correct(xi):
    """
    CORRECT R_ac/R_dc for symmetric excitation.

    The key is proper normalization.

    For H(z) = H0 * cosh(gamma*(a-z)) / cosh(gamma*a):

    J(z) = -dH/dz = H0 * gamma * sinh(gamma*(a-z)) / cosh(gamma*a)

    Current per half = I_half = integral{J} dz from 0 to a
                     = H0 * gamma * integral{sinh(gamma*(a-z))} dz / cosh(gamma*a)
                     = H0 * gamma * [cosh(gamma*a) - 1] / (gamma * cosh(gamma*a))
                     = H0 * (1 - 1/cosh(gamma*a))
                     = H0 * (1 - sech(gamma*a))

    At low frequency (gamma*a -> 0):
    sech(gamma*a) -> 1 - (gamma*a)^2/2
    I_half -> H0 * (gamma*a)^2 / 2

    Power per half = integral{rho * |J|^2} dz from 0 to a
    At low frequency: J ~ H0 * gamma * (a-z)
    power_half ~ |H0|^2 * |gamma|^2 * integral{(a-z)^2} dz = |H0|^2 * |gamma|^2 * a^3/3

    DC case:
    J_dc = I_dc / a
    I_dc = H0 (since H goes from H0 to 0)
    power_dc = |J_dc|^2 * a = |H0|^2 / a

    But wait - at DC, dH/dz = const = H0/a, so J_dc = H0/a
    I_dc = integral{J_dc} dz = H0/a * a = H0

    Actually, let me think about this differently.

    The SURFACE IMPEDANCE approach is CORRECT for AC resistance!

    Z_s = E(0) / H(0) where E(0) = rho * J(0)

    For symmetric geometry:
    J(0) = H0 * gamma * sinh(gamma*a) / cosh(gamma*a) = H0 * gamma * tanh(gamma*a)
    E(0) = rho * H0 * gamma * tanh(gamma*a)
    Z_s = rho * gamma * tanh(gamma*a)

    Z_dc = R_dc = rho * L / A for a conductor
    For unit length and width: A = 2a (full thickness)
    Z_dc = rho / (2a)

    Wait, that's the issue! For SYMMETRIC excitation, current flows through HALF
    the conductor (from surface to center), not the full thickness!

    Actually no. In symmetric excitation:
    - Current enters from top surface, flows to center, reflects due to symmetry
    - Same from bottom surface

    The effective "thickness" for resistance is a (half-thickness), not 2a.

    So: Z_dc = rho / a

    R_ac/R_dc = Re(Z_s) / Z_dc = Re(gamma * a * tanh(gamma*a))
              = Re((1+j)*xi * tanh((1+j)*xi))

    This is exactly what we computed before!
    """
    if xi < 0.001:
        return 1.0

    gamma_a = complex(1, 1) * xi
    return (gamma_a * np.tanh(gamma_a)).real


def dowell_formula(xi):
    """Dowell's formula."""
    if xi < 0.01:
        return 1.0 + xi**4 / 45

    sh2 = np.sinh(2*xi)
    sn2 = np.sin(2*xi)
    ch2 = np.cosh(2*xi)
    cs2 = np.cos(2*xi)

    return xi * (sh2 + sn2) / (ch2 - cs2)


def verify_formula_identity():
    """
    Check if Re((1+j)*xi * tanh((1+j)*xi)) = xi * (sinh(2xi) + sin(2xi)) / (cosh(2xi) - cos(2xi))

    Let z = (1+j)*xi

    tanh(z) = sinh(z) / cosh(z)

    sinh((1+j)*xi) = sinh(xi)*cos(xi) + j*cosh(xi)*sin(xi)
    cosh((1+j)*xi) = cosh(xi)*cos(xi) + j*sinh(xi)*sin(xi)

    |cosh(z)|^2 = cosh^2(xi)*cos^2(xi) + sinh^2(xi)*sin^2(xi)
                = (cosh(2xi) + cos(2xi)) / 2    [using identities]

    sinh(z)*conj(cosh(z)) = (sinh*cos + j*cosh*sin)(cosh*cos - j*sinh*sin)
                          = sinh*cos*cosh*cos + cosh*sin*sinh*sin
                            + j(cosh*sin*cosh*cos - sinh*cos*sinh*sin)
                          = sinh*cosh*(cos^2 + sin^2) + j*cosh*sinh*(cos^2 - sin^2) ... hmm

    Actually let's use:
    tanh(z) * conj(cosh(z)^2) = sinh(z) * conj(cosh(z))

    z * tanh(z) = (1+j)*xi * sinh(z) / cosh(z)

    Real part of z*tanh(z):
    = Re((1+j)*xi * (sinh*cos + j*cosh*sin) / (cosh*cos + j*sinh*sin))

    Let's compute the denominator:
    cosh*cos + j*sinh*sin
    |denom|^2 = cosh^2*cos^2 + sinh^2*sin^2 = (cosh(2xi) + cos(2xi))/2

    denom^(-1) = (cosh*cos - j*sinh*sin) / |denom|^2

    numerator * denom^(-1):
    = (sinh*cos + j*cosh*sin)(cosh*cos - j*sinh*sin) / |denom|^2
    = [sinh*cos*cosh*cos + cosh*sin*sinh*sin + j(...)] / |denom|^2
    = [sinh*cosh*(cos^2+sin^2) + j(...)] / |denom|^2
    = [sinh*cosh + j*...] / |denom|^2
    = [(sinh(2xi)/2) + j*...] / ((cosh(2xi)+cos(2xi))/2)
    = sinh(2xi) / (cosh(2xi)+cos(2xi)) + j*...

    So:
    tanh(z) = sinh(z)/cosh(z) = sinh(2xi)/(cosh(2xi)+cos(2xi)) + j*...

    Hmm, this doesn't match Dowell's denominator (cosh(2xi) - cos(2xi)).

    Let me recalculate more carefully.
    """
    xi = 1.0

    # Using numpy for verification
    z = (1+1j) * xi
    tanh_z = np.tanh(z)
    product = z * tanh_z

    print(f"xi = {xi}")
    print(f"z = {z}")
    print(f"tanh(z) = {tanh_z}")
    print(f"z * tanh(z) = {product}")
    print(f"Re(z * tanh(z)) = {product.real}")

    # Dowell
    dowell = dowell_formula(xi)
    print(f"Dowell formula = {dowell}")

    # The formulas are NOT equal!
    print(f"\nRatio = {product.real / dowell}")


def find_correct_formula():
    """
    Let's work backwards from Dowell's formula to find what it represents.

    Dowell: R_ac/R_dc = xi * (sinh(2xi) + sin(2xi)) / (cosh(2xi) - cos(2xi))

    This can be rewritten using complex exponentials.

    Let q = (1+j)*xi
    cosh(2xi) - cos(2xi) = |sinh(q)|^2 * 2   (need to verify)
    sinh(2xi) + sin(2xi) = ... ?

    Actually, Dowell's formula comes from SINGLE-SIDED geometry where:
    H(0) = 0, H(h) = H0

    The solution is H(z) = H0 * sinh(gamma*z) / sinh(gamma*h)

    Let's verify this gives Dowell.
    """
    print("\n" + "=" * 60)
    print("Single-Sided Geometry Verification")
    print("=" * 60)

    xi_values = [0.5, 1.0, 2.0, 5.0]

    for xi in xi_values:
        gamma_h = complex(1, 1) * xi  # gamma * h where h is full thickness

        # Surface impedance for single-sided:
        # J(h) = H0 * gamma * cosh(gamma*h) / sinh(gamma*h) = H0 * gamma * coth(gamma*h)
        # E(h) = rho * J(h) = rho * H0 * gamma * coth(gamma*h)
        # Z_s = E(h) / H(h) = rho * gamma * coth(gamma*h)
        # Z_dc = rho / h
        # R_ratio = Re(gamma*h * coth(gamma*h))

        coth_gh = 1.0 / np.tanh(gamma_h)
        R_surface_imp = (gamma_h * coth_gh).real

        R_dowell = dowell_formula(xi)

        print(f"xi = {xi}: Surface Imp = {R_surface_imp:.4f}, Dowell = {R_dowell:.4f}, "
              f"Ratio = {R_surface_imp/R_dowell:.4f}")


def main():
    print("=" * 70)
    print("Understanding Dowell's Formula vs Surface Impedance")
    print("=" * 70)

    verify_formula_identity()
    find_correct_formula()

    print("\n" + "=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    print("""
    The formulas are:

    1. SINGLE-SIDED (Dowell's original):
       R_ac/R_dc = Re(gamma*h * coth(gamma*h))
                 = xi * (sinh(2xi) + sin(2xi)) / (cosh(2xi) - cos(2xi))
       where h = full thickness, xi = h/delta

    2. SYMMETRIC (our ESIMFiniteSlabSolver):
       R_ac/R_dc = Re(gamma*a * tanh(gamma*a))
       where a = half-thickness, xi = a/delta

    These are DIFFERENT formulas for DIFFERENT geometries!

    For symmetric geometry, we should use tanh, not Dowell's formula.
    """)

    print("\n" + "=" * 70)
    print("Symmetric Formula Table")
    print("=" * 70)

    print(f"\n{'xi=a/delta':>12} {'R_ac/R_dc':>12}")
    print("-" * 28)

    xi_values = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]
    for xi in xi_values:
        r = symmetric_formula_correct(xi)
        print(f"{xi:>12.2f} {r:>12.4f}")


if __name__ == '__main__':
    main()
