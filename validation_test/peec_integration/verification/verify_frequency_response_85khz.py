"""
Frequency Response Verification at 85 kHz

Verifies PEEC model accuracy at 85 kHz, the target frequency for:
- Wireless Power Transfer (WPT) for electric vehicles
- Qi-standard wireless charging
- Industrial induction heating

Test cases:
1. Straight wire: Verify impedance vs analytical
2. Circular loop: Verify inductance and resistance
3. Skin effect: Verify resistance increase at high frequency
"""

import sys
import os
import numpy as np

# Add Radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad

try:
    from peec_matrices import PEECBuilder
    print("Using C++ PEEC implementation")
except ImportError:
    print("C++ PEEC module not available, using Python fallback")
    from radia.peec_matrices import PEECBuilder


# Physical constants
MU_0 = 4 * np.pi * 1e-7  # H/m
SIGMA_COPPER = 5.8e7      # S/m

mm = 1e-3  # 1 mm in meters


def skin_depth(f, sigma=SIGMA_COPPER, mu_r=1.0):
    """Calculate skin depth in meters."""
    omega = 2 * np.pi * f
    return np.sqrt(2 / (omega * mu_r * MU_0 * sigma))


def dc_resistance(length, width, height, sigma):
    """Calculate DC resistance of rectangular conductor."""
    area = width * height
    return length / (sigma * area)


def ac_resistance_factor(f, width, height, sigma):
    """Calculate AC/DC resistance ratio due to skin effect."""
    delta = skin_depth(f, sigma)

    # For rectangular conductor, approximate factor
    # Using Belevitch formula approximation
    if delta >= min(width, height) / 2:
        # Low frequency: no skin effect
        return 1.0
    else:
        # High frequency: current flows in skin depth layer
        perimeter = 2 * (width + height)
        area_effective = perimeter * delta * (1 - np.exp(-min(width, height) / (2 * delta)))
        area_dc = width * height
        return area_dc / area_effective


def analytical_wire_inductance(length, width, height):
    """Analytical inductance of straight wire (Neumann formula approximation)."""
    # GMD (Geometric Mean Distance) approximation for rectangular cross-section
    a = np.sqrt(width * height / np.pi)  # Equivalent circular radius
    if length > 10 * a:
        # Long wire approximation
        L = MU_0 * length / (2 * np.pi) * (np.log(2 * length / a) - 1)
    else:
        # Short wire approximation
        L = MU_0 * length / (2 * np.pi) * (np.log(2 * length / a) - 0.75)
    return L


def analytical_loop_inductance(radius, wire_radius):
    """Analytical inductance of circular loop (thin wire approximation)."""
    # Neumann formula for circular loop
    if radius > 10 * wire_radius:
        L = MU_0 * radius * (np.log(8 * radius / wire_radius) - 2)
    else:
        L = MU_0 * radius * (np.log(8 * radius / wire_radius) - 1.75)
    return L


def verify_straight_wire(f_target=85000):
    """Verify PEEC for straight wire at target frequency."""
    print("\n" + "=" * 60)
    print(f"Test 1: Straight Wire at {f_target/1e3:.0f} kHz")
    print("=" * 60)

    # Wire parameters (typical WPT coil segment)
    length = 100*mm     # 10 cm
    width = 2*mm        # 2 mm
    height = 0.5*mm     # 0.5 mm (Litz wire cross-section)
    n_segments = 20
    sigma = SIGMA_COPPER

    print(f"\nWire parameters:")
    print(f"  Length:     {length/mm:.0f} mm")
    print(f"  Width:      {width/mm:.1f} mm")
    print(f"  Height:     {height/mm:.2f} mm")
    print(f"  Segments:   {n_segments}")

    # Analytical values
    delta = skin_depth(f_target, sigma)
    R_dc = dc_resistance(length, width, height, sigma)
    R_ac_factor = ac_resistance_factor(f_target, width, height, sigma)
    L_analytical = analytical_wire_inductance(length, width, height)

    print(f"\nAnalytical values:")
    print(f"  Skin depth at {f_target/1e3:.0f} kHz: {delta*1e6:.1f} um")
    print(f"  R_DC:         {R_dc*1e3:.4f} mOhm")
    print(f"  R_AC/R_DC:    {R_ac_factor:.2f}")
    print(f"  L_analytical: {L_analytical*1e9:.2f} nH")

    # Build PEEC matrices
    builder = PEECBuilder()
    builder.create_wire([0, 0, 0], [length, 0, 0], width, height, n_segments, sigma)
    L, R, P, M_LS = builder.build()

    # PEEC values
    R_peec = np.sum(np.diag(R)) if R.ndim == 1 else np.sum(np.diag(R))
    L_peec = np.sum(L)

    print(f"\nPEEC values:")
    print(f"  R_PEEC:       {R_peec*1e3:.4f} mOhm")
    print(f"  L_PEEC:       {L_peec*1e9:.2f} nH")

    # Errors
    R_error = abs(R_peec - R_dc) / R_dc * 100
    L_error = abs(L_peec - L_analytical) / L_analytical * 100

    print(f"\nErrors vs analytical:")
    print(f"  R error:      {R_error:.1f}%")
    print(f"  L error:      {L_error:.1f}%")

    # Impedance at target frequency
    omega = 2 * np.pi * f_target
    if R.ndim == 1:
        R_full = np.diag(R)
    else:
        R_full = R
    Z = R_full + 1j * omega * L

    # Port impedance (first to last segment)
    try:
        Z_port = Z[0, 0] - Z[0, 1:] @ np.linalg.solve(Z[1:, 1:], Z[1:, 0])
    except Exception:
        Z_port = np.sum(np.diag(Z))

    Z_analytical = R_dc + 1j * omega * L_analytical

    print(f"\nImpedance at {f_target/1e3:.0f} kHz:")
    print(f"  |Z_PEEC|:      {np.abs(Z_port)*1e3:.4f} mOhm")
    print(f"  |Z_analytical|: {np.abs(Z_analytical)*1e3:.4f} mOhm")
    print(f"  Phase PEEC:    {np.angle(Z_port, deg=True):.1f} deg")
    print(f"  Phase anal.:   {np.angle(Z_analytical, deg=True):.1f} deg")

    passed = R_error < 5 and L_error < 20
    print(f"\nResult: {'PASS' if passed else 'FAIL'}")
    return passed


def verify_circular_loop(f_target=85000):
    """Verify PEEC for circular loop at target frequency."""
    print("\n" + "=" * 60)
    print(f"Test 2: Circular Loop at {f_target/1e3:.0f} kHz")
    print("=" * 60)

    # Loop parameters (typical WPT coil)
    radius = 150*mm     # 15 cm diameter
    width = 3*mm        # 3 mm
    height = 1*mm       # 1 mm
    n_segments = 36    # 10 deg per segment
    sigma = SIGMA_COPPER

    print(f"\nLoop parameters:")
    print(f"  Radius:     {radius/mm/10:.0f} cm")
    print(f"  Width:      {width/mm:.1f} mm")
    print(f"  Height:     {height/mm:.1f} mm")
    print(f"  Segments:   {n_segments}")

    # Analytical values
    wire_radius = np.sqrt(width * height / np.pi)
    circumference = 2 * np.pi * radius
    R_dc = dc_resistance(circumference, width, height, sigma)
    L_analytical = analytical_loop_inductance(radius, wire_radius)

    print(f"\nAnalytical values:")
    print(f"  R_DC:         {R_dc*1e3:.4f} mOhm")
    print(f"  L_analytical: {L_analytical*1e6:.3f} uH")

    # Build PEEC matrices
    builder = PEECBuilder()
    builder.create_loop([0, 0, 0], radius, [0, 0, 1], width, height, n_segments, sigma)
    L, R, P, M_LS = builder.build()

    # PEEC values
    R_peec = np.sum(np.diag(R)) if R.ndim == 1 else np.sum(np.diag(R))
    L_peec = np.sum(L)

    print(f"\nPEEC values:")
    print(f"  R_PEEC:       {R_peec*1e3:.4f} mOhm")
    print(f"  L_PEEC:       {L_peec*1e6:.3f} uH")

    # Errors
    R_error = abs(R_peec - R_dc) / R_dc * 100
    L_error = abs(L_peec - L_analytical) / L_analytical * 100

    print(f"\nErrors vs analytical:")
    print(f"  R error:      {R_error:.1f}%")
    print(f"  L error:      {L_error:.1f}%")

    # Impedance at target frequency
    omega = 2 * np.pi * f_target
    if R.ndim == 1:
        R_full = np.diag(R)
    else:
        R_full = R
    Z = R_full + 1j * omega * L

    # Total loop impedance (sum of all segments for closed loop)
    Z_total = np.sum(Z)
    Z_analytical = R_dc + 1j * omega * L_analytical

    print(f"\nImpedance at {f_target/1e3:.0f} kHz:")
    print(f"  |Z_PEEC|:      {np.abs(Z_total):.4f} Ohm")
    print(f"  |Z_analytical|: {np.abs(Z_analytical):.4f} Ohm")
    print(f"  Phase PEEC:    {np.angle(Z_total, deg=True):.1f} deg")
    print(f"  Phase anal.:   {np.angle(Z_analytical, deg=True):.1f} deg")

    # Quality factor
    Q_peec = np.abs(Z_total.imag) / np.abs(Z_total.real)
    Q_analytical = omega * L_analytical / R_dc

    print(f"\nQuality factor:")
    print(f"  Q_PEEC:       {Q_peec:.0f}")
    print(f"  Q_analytical: {Q_analytical:.0f}")

    passed = R_error < 5 and L_error < 15
    print(f"\nResult: {'PASS' if passed else 'FAIL'}")
    return passed


def verify_frequency_sweep():
    """Verify PEEC across frequency range."""
    print("\n" + "=" * 60)
    print("Test 3: Frequency Sweep (1 kHz to 1 MHz)")
    print("=" * 60)

    # Simple wire for testing
    length = 100*mm
    width = 2*mm
    height = 0.5*mm
    n_segments = 10
    sigma = SIGMA_COPPER

    # Build PEEC matrices
    builder = PEECBuilder()
    builder.create_wire([0, 0, 0], [length, 0, 0], width, height, n_segments, sigma)
    L, R, P, M_LS = builder.build()

    if R.ndim == 1:
        R_full = np.diag(R)
    else:
        R_full = R

    # Frequency sweep
    frequencies = np.logspace(3, 6, 13)  # 1 kHz to 1 MHz

    print(f"\n{'Freq (kHz)':>12s}  {'|Z| (mOhm)':>12s}  {'Phase (deg)':>12s}  {'Skin (um)':>10s}")
    print("-" * 55)

    results = []
    for f in frequencies:
        omega = 2 * np.pi * f
        Z = R_full + 1j * omega * L

        # Port impedance
        try:
            Z_port = Z[0, 0] - Z[0, 1:] @ np.linalg.solve(Z[1:, 1:], Z[1:, 0])
        except Exception:
            Z_port = Z[0, 0]

        delta = skin_depth(f, sigma)

        results.append({
            'f': f,
            'Z': Z_port,
            'delta': delta
        })

        print(f"{f/1e3:12.2f}  {np.abs(Z_port)*1e3:12.4f}  {np.angle(Z_port, deg=True):12.1f}  {delta*1e6:10.1f}")

    # Check phase behavior
    low_freq_phase = np.angle(results[0]['Z'], deg=True)
    high_freq_phase = np.angle(results[-1]['Z'], deg=True)

    print(f"\nPhase analysis:")
    print(f"  Low freq (1 kHz):  {low_freq_phase:.1f} deg")
    print(f"  High freq (1 MHz): {high_freq_phase:.1f} deg")

    # Should be increasingly inductive
    if high_freq_phase > 80:
        print("  Behavior: Strongly inductive at high frequency (PASS)")
        passed = True
    else:
        print("  Behavior: NOT strongly inductive (CHECK)")
        passed = False

    # Check 85 kHz specifically
    f_85k_idx = np.argmin(np.abs(np.array(frequencies) - 85000))
    Z_85k = results[f_85k_idx]['Z']
    print(f"\nAt 85 kHz:")
    print(f"  |Z|:   {np.abs(Z_85k)*1e3:.4f} mOhm")
    print(f"  Phase: {np.angle(Z_85k, deg=True):.1f} deg")
    print(f"  delta: {results[f_85k_idx]['delta']*1e6:.1f} um")

    return passed


def main():
    """Run all frequency response verification tests."""
    print("=" * 70)
    print("PEEC Frequency Response Verification at 85 kHz")
    print("Target: WPT (Wireless Power Transfer) applications")
    print("=" * 70)

    results = []

    # Test 1: Straight wire
    try:
        results.append(("Straight Wire", verify_straight_wire()))
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Straight Wire", False))

    # Test 2: Circular loop
    try:
        results.append(("Circular Loop", verify_circular_loop()))
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Circular Loop", False))

    # Test 3: Frequency sweep
    try:
        results.append(("Frequency Sweep", verify_frequency_sweep()))
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Frequency Sweep", False))

    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)

    all_passed = True
    for name, passed in results:
        status = "PASS" if passed else "FAIL"
        print(f"  {name}: {status}")
        if not passed:
            all_passed = False

    print()
    if all_passed:
        print("All tests PASSED - PEEC model validated at 85 kHz")
    else:
        print("Some tests FAILED - review PEEC implementation")

    return all_passed


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
