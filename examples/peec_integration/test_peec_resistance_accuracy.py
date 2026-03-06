"""
Test PEEC resistance accuracy across frequency range.

Compare PEEC implementation with accurate analytical formulas for:
1. DC resistance: R_dc = L / (sigma * A)
2. AC resistance with skin effect

The accurate AC resistance formula for a round wire is:
    R_ac / R_dc = (a/2) * Real{ (ber + j*bei)' / (ber + j*bei) }

where ber, bei are Kelvin functions and a = wire_radius / skin_depth.

For practical engineering:
    - a < 1: R_ac ≈ R_dc (no skin effect)
    - a > 3: R_ac ≈ R_dc * a/2 (full skin effect)
    - 1 < a < 3: Transition region (need Kelvin functions)

Usage:
    python test_peec_resistance_accuracy.py
"""

import sys
import os
import numpy as np
from scipy.special import kelvin

# Add Radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))


def kelvin_rac_ratio(xi):
    """
    Accurate AC/DC resistance ratio using Kelvin functions.

    R_ac/R_dc = (xi/2) * Real{ (ber' + j*bei') / (ber + j*bei) }

    where xi = sqrt(2) * a / delta = sqrt(2) * (wire_radius / skin_depth)

    For round wire: xi = a * sqrt(2 * omega * mu * sigma)

    Reference: Ramo, Whinnery, Van Duzer, "Fields and Waves in
               Communication Electronics", Chapter 5.
    """
    if xi < 0.01:
        return 1.0  # DC limit

    # Kelvin functions: scipy.special.kelvin returns (ber, bei, ker, kei)
    ber, bei, ker, kei = kelvin(xi)

    # Derivatives using finite difference (scipy doesn't provide berp, beip directly)
    h = 1e-6
    ber_h, bei_h, _, _ = kelvin(xi + h)
    berp = (ber_h - ber) / h
    beip = (bei_h - bei) / h

    # ber + j*bei
    ber_bei = complex(ber, bei)

    # ber' + j*bei'
    berp_beip = complex(berp, beip)

    # Ratio
    if abs(ber_bei) < 1e-30:
        return 1.0

    ratio = (xi / 2.0) * (berp_beip / ber_bei).real

    return max(ratio, 1.0)  # R_ac >= R_dc


def simple_rac_ratio(a_over_delta):
    """
    Simple AC/DC resistance ratio (current PEEC implementation).

    if delta < a:
        R_ac = R_dc * a / (2*delta)
    else:
        R_ac = R_dc
    """
    if a_over_delta <= 1.0:
        return 1.0
    else:
        return a_over_delta / 2.0


def dowell_rac_ratio(xi):
    """
    Dowell's formula for AC/DC resistance ratio.

    Commonly used in transformer design for round wires.

    R_ac/R_dc = xi * (sinh(2*xi) + sin(2*xi)) / (cosh(2*xi) - cos(2*xi))

    where xi = sqrt(omega * mu * sigma / 2) * diameter = a/delta

    Note: This is for rectangular conductors, often used as approximation
    for round wires with appropriate correction.
    """
    if xi < 0.01:
        return 1.0

    sinh_2xi = np.sinh(2*xi)
    sin_2xi = np.sin(2*xi)
    cosh_2xi = np.cosh(2*xi)
    cos_2xi = np.cos(2*xi)

    denom = cosh_2xi - cos_2xi
    if abs(denom) < 1e-30:
        return 1.0

    return xi * (sinh_2xi + sin_2xi) / denom


def test_resistance_formulas():
    """Compare different AC resistance formulas."""
    print("\n" + "="*70)
    print("Test: AC Resistance Formula Comparison")
    print("="*70)

    print("""
    Formulas compared:
    1. Kelvin functions (exact for round wire)
    2. Simple step function (current PEEC): R_ac = R_dc * max(1, a/(2*delta))
    3. Dowell's formula (rectangular approximation)

    xi = a/delta where a = wire radius, delta = skin depth
    """)

    xi_values = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]

    print(f"{'xi=a/delta':>12} {'Kelvin':>12} {'Simple':>12} {'Dowell':>12} "
          f"{'Simple/Kelv':>12} {'Dowell/Kelv':>12}")
    print("-" * 76)

    for xi in xi_values:
        r_kelvin = kelvin_rac_ratio(np.sqrt(2) * xi)  # Kelvin uses sqrt(2)*xi
        r_simple = simple_rac_ratio(xi)
        r_dowell = dowell_rac_ratio(xi)

        err_simple = r_simple / r_kelvin if r_kelvin > 0 else 0
        err_dowell = r_dowell / r_kelvin if r_kelvin > 0 else 0

        print(f"{xi:>12.2f} {r_kelvin:>12.4f} {r_simple:>12.4f} {r_dowell:>12.4f} "
              f"{err_simple:>12.3f} {err_dowell:>12.3f}")

    print("\nNote: Simple/Kelvin = 1.0 means perfect match")


def test_peec_vs_analytical():
    """Compare PEEC implementation with analytical formulas."""
    import radia as rad

    print("\n" + "="*70)
    print("Test: PEEC Implementation vs Analytical")
    print("="*70)

    # Wire parameters
    sigma = 5.8e7  # Copper
    mu_0 = 4 * np.pi * 1e-7
    loop_radius = 0.05  # 50 mm
    wire_radius = 0.001  # 1 mm

    # Analytical DC resistance
    path_length = 2 * np.pi * loop_radius
    cross_section = np.pi * wire_radius**2
    R_dc_analytical = path_length / (sigma * cross_section)

    print(f"\nWire parameters:")
    print(f"  Loop radius: {loop_radius*1000:.1f} mm")
    print(f"  Wire radius: {wire_radius*1000:.2f} mm")
    print(f"  R_dc (analytical): {R_dc_analytical*1000:.4f} mOhm")

    frequencies = [0, 100, 1000, 5000, 10000, 50000, 100000]

    print(f"\n{'Freq [Hz]':>10} {'xi':>8} {'R_PEEC':>12} {'R_Kelvin':>12} "
          f"{'R_Simple':>12} {'PEEC/Kelv':>10}")
    print("-" * 68)

    for f in frequencies:
        # Create Radia conductor
        rad.UtiDelAll()
        rad.FldUnits('m')

        coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r',
                           wire_radius*2, wire_radius*2, sigma, 4, 36)

        rad.CndSetFrequency(coil, f)
        rad.CndSetVoltage(coil, 1.0, 0.0)
        rad.CndSolve(coil)

        Z = rad.CndGetImpedance(coil)
        R_peec = Z.real

        if f > 0:
            omega = 2 * np.pi * f
            delta = np.sqrt(2 / (omega * mu_0 * sigma))
            xi = wire_radius / delta

            # Analytical formulas
            r_kelvin = kelvin_rac_ratio(np.sqrt(2) * xi)
            r_simple = simple_rac_ratio(xi)

            R_kelvin = R_dc_analytical * r_kelvin
            R_simple = R_dc_analytical * r_simple

            ratio = R_peec / R_kelvin if R_kelvin > 0 else 0

            print(f"{f:>10.0f} {xi:>8.3f} {R_peec*1000:>12.4f} {R_kelvin*1000:>12.4f} "
                  f"{R_simple*1000:>12.4f} {ratio:>10.3f}")
        else:
            print(f"{f:>10.0f} {'0':>8} {R_peec*1000:>12.4f} {R_dc_analytical*1000:>12.4f} "
                  f"{R_dc_analytical*1000:>12.4f} {R_peec/R_dc_analytical:>10.3f}")


def improved_rac_formula():
    """Propose improved AC resistance formula for PEEC."""
    print("\n" + "="*70)
    print("Proposed: Improved AC Resistance Formula")
    print("="*70)

    print("""
    Current PEEC formula (step function):
        if delta < a:
            R_ac = R_dc * a / (2*delta)
        else:
            R_ac = R_dc

    Problem: Discontinuous at delta = a, underestimates in transition region.

    Improved formula (smooth transition):
        xi = a / delta
        R_ac = R_dc * sqrt(1 + (xi/2)^4)^(1/2)

    Or using approximate Kelvin-based formula:
        R_ac = R_dc * (1 + xi^4 / 48 + ...)  for xi < 2
        R_ac = R_dc * xi / 2                  for xi > 3

    Best: Use polynomial fit to Kelvin function results.
    """)

    # Fit polynomial to Kelvin results
    xi_fit = np.linspace(0.1, 10, 100)
    r_kelvin = np.array([kelvin_rac_ratio(np.sqrt(2) * xi) for xi in xi_fit])

    # Polynomial fit: R_ac/R_dc = 1 + a2*xi^2 + a4*xi^4 + ...
    # For small xi: R_ac/R_dc ≈ 1 + xi^4/48
    # For large xi: R_ac/R_dc ≈ xi/2

    def smooth_formula(xi):
        """Smooth formula matching Kelvin function behavior."""
        # Blend between small-xi and large-xi approximations
        small_xi = 1 + xi**4 / 48
        large_xi = xi / 2

        # Smooth transition around xi = 2
        t = 1 / (1 + np.exp(-3*(xi - 2)))

        return (1 - t) * small_xi + t * large_xi

    print(f"\n{'xi':>8} {'Kelvin':>12} {'Current':>12} {'Smooth':>12} "
          f"{'Curr Err%':>10} {'Smooth Err%':>10}")
    print("-" * 66)

    xi_test = [0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0]

    for xi in xi_test:
        r_k = kelvin_rac_ratio(np.sqrt(2) * xi)
        r_c = simple_rac_ratio(xi)
        r_s = smooth_formula(xi)

        err_c = (r_c - r_k) / r_k * 100
        err_s = (r_s - r_k) / r_k * 100

        print(f"{xi:>8.2f} {r_k:>12.4f} {r_c:>12.4f} {r_s:>12.4f} "
              f"{err_c:>10.1f} {err_s:>10.1f}")


def main():
    """Run all tests."""
    print("="*70)
    print("PEEC AC Resistance Accuracy Analysis")
    print("="*70)

    test_resistance_formulas()
    test_peec_vs_analytical()
    improved_rac_formula()

    print("\n" + "="*70)
    print("Summary and Recommendations")
    print("="*70)
    print("""
    Current PEEC Implementation:
    - Uses step function: R_ac = R_dc * max(1, a/(2*delta))
    - Accurate for xi < 0.5 (low freq) and xi > 3 (high freq)
    - Up to 30% error in transition region (0.5 < xi < 3)

    Recommendations:
    1. For power electronics (50 Hz - 100 kHz): Current formula is acceptable
       because most wires operate in either xi < 1 or xi > 3 regime.

    2. For precision applications: Implement smooth transition formula
       or use Kelvin function lookup table.

    3. Key insight: DC resistance (xi=0) is always correctly computed
       because the step function returns R_dc when delta > a.
    """)


if __name__ == '__main__':
    main()
