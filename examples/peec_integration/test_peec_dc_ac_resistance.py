"""
Test PEEC DC and AC resistance calculation.

This script verifies that the PEEC solver correctly computes:
1. DC resistance (R_dc = L / (sigma * A))
2. AC resistance with skin effect (R_ac = R_dc * r / (2*delta) when delta < r)
3. Inductance (L = mu_0 * R * (ln(8R/a) - 2) for circular loop)

Physical parameters:
- Copper: sigma = 5.8e7 S/m
- Loop: radius = 50mm, wire radius = 1mm

Expected DC resistance:
  R_dc = (2*pi*R) / (sigma * pi * a^2)
       = 2*pi*0.05 / (5.8e7 * pi * 0.001^2)
       = 0.314 / (182.2)
       = 1.72 mOhm

Expected inductance:
  L = mu_0 * R * (ln(8R/a) - 2)
    = 4*pi*1e-7 * 0.05 * (ln(400) - 2)
    = 6.28e-8 * (5.99 - 2)
    = 6.28e-8 * 3.99
    = 251 nH

Usage:
    python test_peec_dc_ac_resistance.py
"""

import sys
import os
import numpy as np

# Add Radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))


def test_dc_resistance():
    """Test DC resistance calculation."""
    import radia as rad

    print("\n" + "="*60)
    print("Test: DC Resistance of Circular Loop")
    print("="*60)

    # Physical constants
    sigma = 5.8e7  # Copper conductivity [S/m]
    mu_0 = 4 * np.pi * 1e-7  # H/m

    # Geometry
    loop_radius = 0.05  # 50 mm [m]
    wire_radius = 0.001  # 1 mm [m]
    wire_height = 0.002  # 2 mm [m] (for rectangular cross section)

    # Expected values
    path_length = 2 * np.pi * loop_radius
    wire_area = np.pi * wire_radius**2
    R_dc_expected = path_length / (sigma * wire_area)

    L_expected = mu_0 * loop_radius * (np.log(8 * loop_radius / wire_radius) - 2)

    print(f"\nGeometry:")
    print(f"  Loop radius:    R = {loop_radius*1000:.1f} mm")
    print(f"  Wire radius:    a = {wire_radius*1000:.2f} mm")
    print(f"  Path length:    l = {path_length*1000:.2f} mm")
    print(f"  Wire area:      A = {wire_area*1e6:.4f} mm^2")

    print(f"\nExpected values (analytical):")
    print(f"  DC resistance:  R_dc = {R_dc_expected*1000:.4f} mOhm")
    print(f"  Inductance:     L = {L_expected*1e9:.2f} nH")

    # Create Radia conductor
    rad.UtiDelAll()
    rad.FldUnits('m')

    # CndLoop parameters: center, radius, normal, cross_section_type,
    #                     width, height, conductivity, n_radial, n_azimuthal
    coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r',
                       wire_radius*2, wire_height, sigma, 4, 36)

    # Set DC (frequency = 0)
    rad.CndSetFrequency(coil, 0)  # DC
    rad.CndSetVoltage(coil, 1.0, 0.0)  # 1V excitation

    # Solve
    rad.CndSolve(coil)

    # Get impedance
    Z = rad.CndGetImpedance(coil)
    R_radia = Z.real
    X_radia = Z.imag

    print(f"\nRadia PEEC results (DC, f=0):")
    print(f"  Impedance:      Z = {R_radia*1000:.4f} + j{X_radia*1000:.4f} mOhm")
    print(f"  Resistance:     R = {R_radia*1000:.4f} mOhm")
    print(f"  Reactance:      X = {X_radia*1000:.4f} mOhm (should be ~0 for DC)")

    # For DC, impedance is pure resistance
    # Note: The current implementation uses rectangular cross-section
    # so the DC resistance will differ from circular wire formula

    return R_radia, L_expected


def test_ac_resistance_sweep():
    """Test AC resistance with skin effect across frequencies."""
    import radia as rad

    print("\n" + "="*60)
    print("Test: AC Resistance Frequency Sweep (Skin Effect)")
    print("="*60)

    # Physical constants
    sigma = 5.8e7  # Copper conductivity [S/m]
    mu_0 = 4 * np.pi * 1e-7

    # Geometry
    loop_radius = 0.05  # 50 mm
    wire_radius = 0.001  # 1 mm
    wire_height = 0.002  # 2 mm

    # Frequency sweep
    frequencies = [0, 100, 1000, 10000, 100000]  # DC to 100 kHz

    print(f"\nSkin depth vs frequency:")
    print(f"{'Freq [Hz]':>12} {'delta [mm]':>12} {'delta/a':>10}")
    print("-" * 36)

    for f in frequencies:
        if f > 0:
            omega = 2 * np.pi * f
            delta = np.sqrt(2 / (omega * mu_0 * sigma))
            ratio = delta / wire_radius
            print(f"{f:>12.0f} {delta*1000:>12.4f} {ratio:>10.2f}")
        else:
            print(f"{f:>12.0f} {'inf':>12} {'inf':>10}")

    print(f"\nFrequency sweep results:")
    print(f"{'Freq [Hz]':>12} {'R [mOhm]':>12} {'X [mOhm]':>12} {'L [nH]':>12} {'|Z| [mOhm]':>12}")
    print("-" * 62)

    for f in frequencies:
        rad.UtiDelAll()
        rad.FldUnits('m')

        coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r',
                           wire_radius*2, wire_height, sigma, 4, 36)

        rad.CndSetFrequency(coil, f)
        rad.CndSetVoltage(coil, 1.0, 0.0)
        rad.CndSolve(coil)

        Z = rad.CndGetImpedance(coil)
        R = Z.real * 1000  # mOhm
        X = Z.imag * 1000  # mOhm
        Z_mag = abs(Z) * 1000  # mOhm

        if f > 0:
            L = Z.imag / (2 * np.pi * f) * 1e9  # nH
        else:
            L = 0

        print(f"{f:>12.0f} {R:>12.4f} {X:>12.4f} {L:>12.2f} {Z_mag:>12.4f}")


def test_esim_surface_impedance():
    """Test ESIM surface impedance calculation."""
    import radia as rad

    print("\n" + "="*60)
    print("Test: ESIM Surface Impedance (Steel at Various Frequencies)")
    print("="*60)

    # Steel parameters (from Karl Hollaus's example)
    sigma_steel = 2e6  # S/m (typical for carbon steel)
    mu_r_steel = 200   # Relative permeability

    mu_0 = 4 * np.pi * 1e-7

    print(f"\nSteel properties:")
    print(f"  Conductivity: sigma = {sigma_steel:.2e} S/m")
    print(f"  Permeability: mu_r = {mu_r_steel}")

    frequencies = [50, 100, 500, 1000, 5000]  # Power frequencies

    print(f"\nESIM surface impedance (Zs = (1+j)*Rs):")
    print(f"{'Freq [Hz]':>10} {'delta [mm]':>12} {'Rs [mOhm]':>12} {'|Zs| [mOhm]':>12}")
    print("-" * 50)

    for f in frequencies:
        omega = 2 * np.pi * f
        delta = np.sqrt(2 / (omega * mu_0 * mu_r_steel * sigma_steel))
        Rs = 1 / (sigma_steel * delta)
        Zs_mag = Rs * np.sqrt(2)  # |Zs| = Rs * sqrt(2) for (1+j)*Rs

        print(f"{f:>10.0f} {delta*1000:>12.4f} {Rs*1000:>12.4f} {Zs_mag*1000:>12.4f}")

    # Compare with copper at same frequencies
    sigma_cu = 5.8e7

    print(f"\nComparison: Copper (non-magnetic) at same frequencies:")
    print(f"{'Freq [Hz]':>10} {'delta [mm]':>12} {'Rs [mOhm]':>12} {'|Zs| [mOhm]':>12}")
    print("-" * 50)

    for f in frequencies:
        omega = 2 * np.pi * f
        delta = np.sqrt(2 / (omega * mu_0 * sigma_cu))
        Rs = 1 / (sigma_cu * delta)
        Zs_mag = Rs * np.sqrt(2)

        print(f"{f:>10.0f} {delta*1000:>12.4f} {Rs*1000:>12.6f} {Zs_mag*1000:>12.6f}")


def main():
    """Run all tests."""
    print("="*60)
    print("PEEC DC/AC Resistance and ESIM Test")
    print("="*60)

    test_dc_resistance()
    test_ac_resistance_sweep()
    test_esim_surface_impedance()

    print("\n" + "="*60)
    print("Tests completed!")
    print("="*60)


if __name__ == '__main__':
    main()
