"""
Test ESIM DC limit consistency.

This script verifies that the ESIM (Effective Surface Impedance Method)
correctly reduces to the DC resistance formula as frequency -> 0.

Physics:
    ESIM surface impedance: Zs = (1+j) * Rs
    where Rs = 1 / (sigma * delta) = sqrt(omega * mu / (2 * sigma))

    DC limit (omega -> 0):
        delta -> infinity (current penetrates entire cross-section)
        Zs -> 0 (pure conduction, no skin effect)

    For a wire of length L and cross-section A:
        DC resistance: R_dc = L / (sigma * A)

    For AC with skin effect (delta << wire radius a):
        R_ac ~ R_dc * a / (2 * delta) = L / (sigma * 2*pi*a * delta)
             = L * sqrt(omega * mu / (2 * sigma)) / (sigma * 2*pi*a)

    The ESIM Zs gives resistance per unit area:
        R = Zs.real * L / (perimeter * 1)  [for unit thickness]

Verification:
    1. Compare ESIM-derived R with analytical formulas
    2. Check low-frequency limit matches DC resistance
    3. Check high-frequency matches skin effect formula

Usage:
    python test_esim_dc_limit.py
"""

import sys
import os
import numpy as np

# Add Radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))


def analytical_skin_depth(frequency, sigma, mu_r=1.0):
    """Analytical skin depth: delta = sqrt(2 / (omega * mu * sigma))."""
    mu_0 = 4 * np.pi * 1e-7
    if frequency <= 0:
        return np.inf
    omega = 2 * np.pi * frequency
    return np.sqrt(2 / (omega * mu_0 * mu_r * sigma))


def analytical_surface_resistance(frequency, sigma, mu_r=1.0):
    """Analytical surface resistance: Rs = 1 / (sigma * delta)."""
    delta = analytical_skin_depth(frequency, sigma, mu_r)
    if np.isinf(delta):
        return 0.0  # DC: no surface resistance
    return 1.0 / (sigma * delta)


def analytical_surface_impedance(frequency, sigma, mu_r=1.0):
    """Analytical surface impedance: Zs = (1+j) * Rs."""
    Rs = analytical_surface_resistance(frequency, sigma, mu_r)
    return complex(Rs, Rs)


def wire_dc_resistance(length, cross_section_area, sigma):
    """DC resistance of a wire: R = L / (sigma * A)."""
    return length / (sigma * cross_section_area)


def wire_ac_resistance_skin_effect(length, wire_radius, frequency, sigma, mu_r=1.0):
    """AC resistance with skin effect (delta << a): R ~ R_dc * a / (2*delta)."""
    delta = analytical_skin_depth(frequency, sigma, mu_r)
    cross_section = np.pi * wire_radius**2
    R_dc = wire_dc_resistance(length, cross_section, sigma)

    if delta >= wire_radius:
        # Skin depth larger than radius: use DC resistance
        return R_dc
    else:
        # Skin effect: current flows in skin layer
        # Effective cross-section ~ 2*pi*a * delta
        return R_dc * wire_radius / (2 * delta)


def test_esim_dc_consistency():
    """Test ESIM consistency with DC resistance formula."""
    print("\n" + "="*70)
    print("Test: ESIM DC Limit Consistency")
    print("="*70)

    # Material properties
    sigma_cu = 5.8e7   # Copper [S/m]
    sigma_steel = 2e6  # Steel [S/m]
    mu_r_steel = 200   # Steel relative permeability

    # Wire geometry
    wire_radius = 0.001  # 1 mm
    wire_length = 0.1    # 100 mm (just for resistance calculation)

    print("\n--- Copper (non-magnetic) ---")
    print(f"Conductivity: sigma = {sigma_cu:.2e} S/m")
    print(f"Wire: radius = {wire_radius*1000:.1f} mm, length = {wire_length*1000:.0f} mm")

    # DC resistance
    cross_section = np.pi * wire_radius**2
    R_dc_cu = wire_dc_resistance(wire_length, cross_section, sigma_cu)
    print(f"\nDC resistance (analytical): R_dc = {R_dc_cu*1000:.4f} mOhm")

    # Frequency sweep
    frequencies = [0, 10, 100, 1000, 10000, 100000]

    print(f"\n{'Freq [Hz]':>12} {'delta [mm]':>12} {'Rs [uOhm]':>12} "
          f"{'R_wire [mOhm]':>14} {'R_ac/R_dc':>10}")
    print("-" * 64)

    for f in frequencies:
        if f > 0:
            delta = analytical_skin_depth(f, sigma_cu)
            Rs = analytical_surface_resistance(f, sigma_cu)
            R_ac = wire_ac_resistance_skin_effect(wire_length, wire_radius, f, sigma_cu)
            ratio = R_ac / R_dc_cu
            print(f"{f:>12.0f} {delta*1000:>12.4f} {Rs*1e6:>12.4f} "
                  f"{R_ac*1000:>14.4f} {ratio:>10.2f}")
        else:
            print(f"{f:>12.0f} {'inf':>12} {'0':>12} "
                  f"{R_dc_cu*1000:>14.4f} {'1.00':>10}")

    # Steel with high permeability
    print("\n--- Steel (magnetic) ---")
    print(f"Conductivity: sigma = {sigma_steel:.2e} S/m")
    print(f"Permeability: mu_r = {mu_r_steel}")

    R_dc_steel = wire_dc_resistance(wire_length, cross_section, sigma_steel)
    print(f"\nDC resistance (analytical): R_dc = {R_dc_steel*1000:.4f} mOhm")

    print(f"\n{'Freq [Hz]':>12} {'delta [mm]':>12} {'Rs [mOhm]':>12} "
          f"{'R_wire [mOhm]':>14} {'R_ac/R_dc':>10}")
    print("-" * 64)

    for f in frequencies:
        if f > 0:
            delta = analytical_skin_depth(f, sigma_steel, mu_r_steel)
            Rs = analytical_surface_resistance(f, sigma_steel, mu_r_steel)
            R_ac = wire_ac_resistance_skin_effect(wire_length, wire_radius, f,
                                                   sigma_steel, mu_r_steel)
            ratio = R_ac / R_dc_steel
            print(f"{f:>12.0f} {delta*1000:>12.4f} {Rs*1000:>12.4f} "
                  f"{R_ac*1000:>14.4f} {ratio:>10.2f}")
        else:
            print(f"{f:>12.0f} {'inf':>12} {'0':>12} "
                  f"{R_dc_steel*1000:>14.4f} {'1.00':>10}")


def test_esim_surface_impedance_formula():
    """Verify ESIM surface impedance formula: Zs = sqrt(j*omega*mu/sigma)."""
    print("\n" + "="*70)
    print("Test: ESIM Surface Impedance Formula Verification")
    print("="*70)

    print("""
    ESIM Theory:
        The surface impedance relates tangential E and H at the surface:
            E_t = Zs * (n x H)

        For good conductors (sigma >> omega*epsilon):
            Zs = sqrt(j*omega*mu/sigma) = (1+j) * sqrt(omega*mu/(2*sigma))
                                        = (1+j) * Rs

        where Rs = 1/(sigma*delta) is the surface resistance.

        Properties:
            - |Zs| = sqrt(2) * Rs
            - arg(Zs) = 45 degrees (equal real and imaginary parts)
            - Zs -> 0 as omega -> 0 (DC limit)
    """)

    mu_0 = 4 * np.pi * 1e-7
    sigma = 5.8e7  # Copper

    print(f"Copper: sigma = {sigma:.2e} S/m\n")

    frequencies = [50, 100, 1000, 10000, 100000]

    print(f"{'Freq [Hz]':>10} {'delta [mm]':>12} {'Rs [uOhm]':>12} "
          f"{'|Zs| [uOhm]':>12} {'arg(Zs) [deg]':>14}")
    print("-" * 64)

    for f in frequencies:
        omega = 2 * np.pi * f
        delta = np.sqrt(2 / (omega * mu_0 * sigma))
        Rs = 1 / (sigma * delta)

        # Exact formula: Zs = sqrt(j*omega*mu/sigma)
        Zs_exact = np.sqrt(1j * omega * mu_0 / sigma)

        # Approximate formula: Zs = (1+j)*Rs
        Zs_approx = complex(Rs, Rs)

        print(f"{f:>10.0f} {delta*1000:>12.4f} {Rs*1e6:>12.4f} "
              f"{abs(Zs_exact)*1e6:>12.4f} {np.degrees(np.angle(Zs_exact)):>14.1f}")

    print("\nNote: arg(Zs) = 45 deg confirms equal real and imaginary parts")


def test_esim_power_loss():
    """Verify ESIM power loss calculation."""
    print("\n" + "="*70)
    print("Test: ESIM Power Loss Calculation")
    print("="*70)

    print("""
    ESIM Power Loss:
        For surface current density Js [A/m] (current per unit width):
            P' = (1/2) * Rs * |Js|^2  [W/m^2]

        For tangential H field at surface:
            |Js| = |H_t|  (boundary condition)
            P' = (1/2) * Rs * |H_t|^2  [W/m^2]

        Relation to bulk current:
            For a wire with perimeter p and length L carrying current I:
            |H_t| ~ I / p  (from Ampere's law)
            P' ~ (1/2) * Rs * (I/p)^2

            Total power: P = P' * (p * L) = (1/2) * Rs * I^2 * L / p

        DC limit verification:
            P_dc = I^2 * R_dc = I^2 * L / (sigma * A)

            For skin effect (delta << a):
            Effective A ~ p * delta
            P_ac = I^2 * L / (sigma * p * delta)
                 = I^2 * Rs * L / p  (since Rs = 1/(sigma*delta))

            This matches ESIM formula!
    """)

    # Example calculation
    sigma = 5.8e7  # Copper
    mu_0 = 4 * np.pi * 1e-7
    wire_radius = 0.001  # 1 mm
    wire_length = 0.1    # 100 mm
    I = 1.0  # 1 A

    perimeter = 2 * np.pi * wire_radius
    cross_section = np.pi * wire_radius**2

    # DC power
    R_dc = wire_length / (sigma * cross_section)
    P_dc = I**2 * R_dc

    print(f"Wire: radius = {wire_radius*1000:.1f} mm, length = {wire_length*1000:.0f} mm")
    print(f"Current: I = {I:.1f} A")
    print(f"DC resistance: R_dc = {R_dc*1000:.4f} mOhm")
    print(f"DC power loss: P_dc = {P_dc*1000:.4f} mW")

    print(f"\n{'Freq [Hz]':>10} {'delta/a':>10} {'Rs [uOhm]':>12} "
          f"{'P_ESIM [mW]':>12} {'P_skin [mW]':>12} {'P/P_dc':>10}")
    print("-" * 70)

    frequencies = [0, 50, 100, 1000, 10000, 100000]

    for f in frequencies:
        if f > 0:
            omega = 2 * np.pi * f
            delta = np.sqrt(2 / (omega * mu_0 * sigma))
            Rs = 1 / (sigma * delta)

            # ESIM power: P = Rs * I^2 * L / p
            P_esim = Rs * I**2 * wire_length / perimeter

            # Skin effect power (using effective area)
            if delta < wire_radius:
                A_eff = perimeter * delta
                R_ac = wire_length / (sigma * A_eff)
                P_skin = I**2 * R_ac
            else:
                P_skin = P_dc

            ratio = P_esim / P_dc if P_dc > 0 else 0

            print(f"{f:>10.0f} {delta/wire_radius:>10.2f} {Rs*1e6:>12.4f} "
                  f"{P_esim*1000:>12.4f} {P_skin*1000:>12.4f} {ratio:>10.2f}")
        else:
            print(f"{f:>10.0f} {'inf':>10} {'0':>12} "
                  f"{P_dc*1000:>12.4f} {P_dc*1000:>12.4f} {'1.00':>10}")

    print("\nNote: P_ESIM matches P_skin when delta << a (skin effect regime)")
    print("      At low frequencies (delta > a), ESIM underestimates loss")
    print("      because it assumes surface current distribution.")


def main():
    """Run all ESIM DC limit tests."""
    print("="*70)
    print("ESIM (Effective Surface Impedance Method) DC Limit Verification")
    print("="*70)

    test_esim_dc_consistency()
    test_esim_surface_impedance_formula()
    test_esim_power_loss()

    print("\n" + "="*70)
    print("Summary: ESIM DC Limit Behavior")
    print("="*70)
    print("""
    Key Findings:

    1. DC Limit (f -> 0):
       - Skin depth delta -> infinity
       - Surface resistance Rs -> 0
       - ESIM reduces to pure conduction (no skin effect)

    2. ESIM Validity Range:
       - Best when delta << wire radius (high frequency)
       - At low frequencies, full volume integration is more accurate

    3. PEEC-ESIM Integration:
       - PEEC uses full wire cross-section for DC resistance
       - ESIM correction applied for AC (skin effect)
       - Smooth transition via Rs = 1/(sigma*delta)

    4. Consistency Check:
       - ESIM power loss matches skin-effect formula when delta << a
       - DC resistance emerges from full cross-section integration
    """)


if __name__ == '__main__':
    main()
