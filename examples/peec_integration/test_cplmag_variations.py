"""
Test CplMag solver with various configurations.

Tests:
1. Different core sizes
2. Different core positions
3. Different mu_r values
4. Different frequencies
"""

import sys
sys.path.insert(0, '../../src/radia')
import radia as rad
import numpy as np

def test_core_size_variation():
    """Test inductance vs core size."""
    print("\n" + "="*60)
    print("Test 1: Core Size Variation")
    print("="*60)

    rad.UtiDelAll()
    rad.FldUnits('m')

    loop_radius = 0.05  # 50 mm
    freq = 1000  # 1 kHz
    mu_r = 1000

    # Get baseline (air core) inductance
    coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)
    rad.CndSetFrequency(coil, freq)
    rad.CndSetVoltage(coil, 1.0, 0.0)
    rad.CndSolve(coil)
    Z_air = rad.CndGetImpedance(coil)
    L_air = Z_air.imag / (2 * np.pi * freq)

    print(f"Air core L = {L_air * 1e9:.2f} nH")
    print()

    core_sizes = [0.01, 0.02, 0.03, 0.04]  # 10mm to 40mm cube

    print(f"{'Core Size [mm]':>15} {'L [nH]':>12} {'Ratio':>10}")
    print("-"*40)

    for size in core_sizes:
        rad.UtiDelAll()
        rad.FldUnits('m')

        # Create coil
        coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)

        # Create core (cube centered at origin)
        core = rad.ObjRecMag([0, 0, 0], [size, size, size], [0, 0, 0])
        mat = rad.MatLin(mu_r)
        rad.MatApl(core, mat)

        # Create and run CplMag solver
        solver = rad.CplMagCreate(coil, core)
        rad.CplMagSetFrequency(solver, freq)
        rad.CplMagSetVoltage(solver, 1.0, 0.0)
        rad.CplMagSetMu(solver, mu_r, 0)

        result = rad.CplMagSolve(solver)
        Z = result['Z']
        L = Z.imag / (2 * np.pi * freq)
        ratio = L / L_air

        print(f"{size*1000:>15.0f} {L*1e9:>12.2f} {ratio:>10.2f}")

        rad.CplMagDelete(solver)

    print()


def test_mu_r_variation():
    """Test inductance vs permeability."""
    print("\n" + "="*60)
    print("Test 2: Permeability Variation")
    print("="*60)

    rad.UtiDelAll()
    rad.FldUnits('m')

    loop_radius = 0.05  # 50 mm
    freq = 1000  # 1 kHz
    core_size = 0.03  # 30mm cube

    # Get baseline (air core) inductance
    coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)
    rad.CndSetFrequency(coil, freq)
    rad.CndSetVoltage(coil, 1.0, 0.0)
    rad.CndSolve(coil)
    Z_air = rad.CndGetImpedance(coil)
    L_air = Z_air.imag / (2 * np.pi * freq)

    print(f"Air core L = {L_air * 1e9:.2f} nH")
    print()

    mu_r_values = [10, 100, 500, 1000, 5000]

    print(f"{'mu_r':>10} {'L [nH]':>12} {'Ratio':>10}")
    print("-"*35)

    for mu_r in mu_r_values:
        rad.UtiDelAll()
        rad.FldUnits('m')

        # Create coil
        coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)

        # Create core
        core = rad.ObjRecMag([0, 0, 0], [core_size, core_size, core_size], [0, 0, 0])
        mat = rad.MatLin(mu_r)
        rad.MatApl(core, mat)

        # Create and run CplMag solver
        solver = rad.CplMagCreate(coil, core)
        rad.CplMagSetFrequency(solver, freq)
        rad.CplMagSetVoltage(solver, 1.0, 0.0)
        rad.CplMagSetMu(solver, mu_r, 0)

        result = rad.CplMagSolve(solver)
        Z = result['Z']
        L = Z.imag / (2 * np.pi * freq)
        ratio = L / L_air

        print(f"{mu_r:>10} {L*1e9:>12.2f} {ratio:>10.2f}")

        rad.CplMagDelete(solver)

    print()


def test_frequency_sweep():
    """Test impedance vs frequency."""
    print("\n" + "="*60)
    print("Test 3: Frequency Sweep")
    print("="*60)

    rad.UtiDelAll()
    rad.FldUnits('m')

    loop_radius = 0.05  # 50 mm
    core_size = 0.03    # 30mm cube
    mu_r = 1000

    # Create coil
    coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)

    # Create core
    core = rad.ObjRecMag([0, 0, 0], [core_size, core_size, core_size], [0, 0, 0])
    mat = rad.MatLin(mu_r)
    rad.MatApl(core, mat)

    # Create CplMag solver
    solver = rad.CplMagCreate(coil, core)
    rad.CplMagSetMu(solver, mu_r, 50)  # Add some magnetic loss

    frequencies = [100, 1000, 10000, 100000]

    print(f"{'Freq [Hz]':>12} {'R [mOhm]':>12} {'L [nH]':>12} {'Q':>10}")
    print("-"*50)

    for freq in frequencies:
        rad.CplMagSetFrequency(solver, freq)
        rad.CplMagSetVoltage(solver, 1.0, 0.0)

        result = rad.CplMagSolve(solver)
        Z = result['Z']
        omega = 2 * np.pi * freq
        L = Z.imag / omega
        Q = Z.imag / Z.real if Z.real > 0 else 0

        print(f"{freq:>12.0f} {Z.real*1000:>12.4f} {L*1e9:>12.2f} {Q:>10.1f}")

    rad.CplMagDelete(solver)
    print()


def test_core_position():
    """Test inductance vs core z-position."""
    print("\n" + "="*60)
    print("Test 4: Core Position Variation (z-offset)")
    print("="*60)

    rad.UtiDelAll()
    rad.FldUnits('m')

    loop_radius = 0.05  # 50 mm
    freq = 1000  # 1 kHz
    core_size = 0.02    # 20mm cube
    mu_r = 1000

    # Get baseline (air core) inductance
    coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)
    rad.CndSetFrequency(coil, freq)
    rad.CndSetVoltage(coil, 1.0, 0.0)
    rad.CndSolve(coil)
    Z_air = rad.CndGetImpedance(coil)
    L_air = Z_air.imag / (2 * np.pi * freq)

    print(f"Air core L = {L_air * 1e9:.2f} nH")
    print()

    z_offsets = [0, 0.01, 0.02, 0.03, 0.05]  # 0 to 50mm

    print(f"{'Z offset [mm]':>15} {'L [nH]':>12} {'Ratio':>10}")
    print("-"*40)

    for z_off in z_offsets:
        rad.UtiDelAll()
        rad.FldUnits('m')

        # Create coil (at z=0)
        coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)

        # Create core at z-offset
        core = rad.ObjRecMag([0, 0, z_off], [core_size, core_size, core_size], [0, 0, 0])
        mat = rad.MatLin(mu_r)
        rad.MatApl(core, mat)

        # Create and run CplMag solver
        solver = rad.CplMagCreate(coil, core)
        rad.CplMagSetFrequency(solver, freq)
        rad.CplMagSetVoltage(solver, 1.0, 0.0)
        rad.CplMagSetMu(solver, mu_r, 0)

        result = rad.CplMagSolve(solver)
        Z = result['Z']
        L = Z.imag / (2 * np.pi * freq)
        ratio = L / L_air

        print(f"{z_off*1000:>15.0f} {L*1e9:>12.2f} {ratio:>10.2f}")

        rad.CplMagDelete(solver)

    print()


def test_complex_permeability():
    """Test with complex permeability (magnetic loss)."""
    print("\n" + "="*60)
    print("Test 5: Complex Permeability (Magnetic Loss)")
    print("="*60)

    rad.UtiDelAll()
    rad.FldUnits('m')

    loop_radius = 0.05  # 50 mm
    freq = 10000  # 10 kHz
    core_size = 0.03    # 30mm cube
    mu_r_real = 1000

    # Create coil
    coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)

    # Create core
    core = rad.ObjRecMag([0, 0, 0], [core_size, core_size, core_size], [0, 0, 0])
    mat = rad.MatLin(mu_r_real)
    rad.MatApl(core, mat)

    # Create CplMag solver
    solver = rad.CplMagCreate(coil, core)
    rad.CplMagSetFrequency(solver, freq)
    rad.CplMagSetVoltage(solver, 1.0, 0.0)

    mu_r_imag_values = [0, 10, 50, 100, 200]

    print(f"mu'_r = {mu_r_real}, f = {freq/1000:.0f} kHz")
    print()
    print(f"{'mu\"_r':>10} {'R [mOhm]':>12} {'L [nH]':>12} {'P_mag [W]':>12}")
    print("-"*50)

    for mu_r_imag in mu_r_imag_values:
        rad.CplMagSetMu(solver, mu_r_real, mu_r_imag)

        result = rad.CplMagSolve(solver)
        Z = result['Z']
        omega = 2 * np.pi * freq
        L = Z.imag / omega
        P_mag = result['P_magnet']

        print(f"{mu_r_imag:>10} {Z.real*1000:>12.4f} {L*1e9:>12.2f} {P_mag:>12.4f}")

    rad.CplMagDelete(solver)
    print()


def main():
    """Run all tests."""
    print("="*60)
    print("CplMag Solver Verification Tests")
    print("="*60)

    test_core_size_variation()
    test_mu_r_variation()
    test_frequency_sweep()
    test_core_position()
    test_complex_permeability()

    print("="*60)
    print("All tests completed!")
    print("="*60)


if __name__ == '__main__':
    main()
