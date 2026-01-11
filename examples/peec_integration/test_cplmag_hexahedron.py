"""
Test CplMag solver with ObjHexahedron (MSC method).

This verifies that the coupling matrices work correctly with hexahedral
elements using the Magnetic Surface Charge (MSC) method.
"""

import sys
sys.path.insert(0, '../../src/radia')
import radia as rad
import numpy as np

def create_hexahedron_core(center, size):
    """Create a hexahedral core using ObjHexahedron."""
    cx, cy, cz = center
    sx, sy, sz = size if isinstance(size, (list, tuple)) else (size, size, size)

    # 8 vertices of hexahedron
    vertices = [
        [cx - sx/2, cy - sy/2, cz - sz/2],
        [cx + sx/2, cy - sy/2, cz - sz/2],
        [cx + sx/2, cy + sy/2, cz - sz/2],
        [cx - sx/2, cy + sy/2, cz - sz/2],
        [cx - sx/2, cy - sy/2, cz + sz/2],
        [cx + sx/2, cy - sy/2, cz + sz/2],
        [cx + sx/2, cy + sy/2, cz + sz/2],
        [cx - sx/2, cy + sy/2, cz + sz/2],
    ]

    return rad.ObjHexahedron(vertices, [0, 0, 0])


def test_hexahedron_vs_recmag():
    """Compare ObjHexahedron vs ObjRecMag results."""
    print("\n" + "="*60)
    print("Test: ObjHexahedron vs ObjRecMag Comparison")
    print("="*60)

    loop_radius = 0.05  # 50 mm
    freq = 1000  # 1 kHz
    core_size = 0.03  # 30mm cube
    mu_r = 1000

    # Get baseline (air core) inductance
    rad.UtiDelAll()
    rad.FldUnits('m')
    coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)
    rad.CndSetFrequency(coil, freq)
    rad.CndSetVoltage(coil, 1.0, 0.0)
    rad.CndSolve(coil)
    Z_air = rad.CndGetImpedance(coil)
    L_air = Z_air.imag / (2 * np.pi * freq)
    print(f"Air core L = {L_air * 1e9:.2f} nH")

    # Test with ObjRecMag
    rad.UtiDelAll()
    rad.FldUnits('m')
    coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)
    core_rec = rad.ObjRecMag([0, 0, 0], [core_size, core_size, core_size], [0, 0, 0])
    mat = rad.MatLin(mu_r)
    rad.MatApl(core_rec, mat)

    solver_rec = rad.CplMagCreate(coil, core_rec)
    rad.CplMagSetFrequency(solver_rec, freq)
    rad.CplMagSetVoltage(solver_rec, 1.0, 0.0)
    rad.CplMagSetMu(solver_rec, mu_r, 0)

    result_rec = rad.CplMagSolve(solver_rec)
    Z_rec = result_rec['Z']
    L_rec = Z_rec.imag / (2 * np.pi * freq)
    rad.CplMagDelete(solver_rec)

    # Test with ObjHexahedron
    rad.UtiDelAll()
    rad.FldUnits('m')
    coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)
    core_hex = create_hexahedron_core([0, 0, 0], core_size)
    mat = rad.MatLin(mu_r)
    rad.MatApl(core_hex, mat)

    solver_hex = rad.CplMagCreate(coil, core_hex)
    rad.CplMagSetFrequency(solver_hex, freq)
    rad.CplMagSetVoltage(solver_hex, 1.0, 0.0)
    rad.CplMagSetMu(solver_hex, mu_r, 0)

    result_hex = rad.CplMagSolve(solver_hex)
    Z_hex = result_hex['Z']
    L_hex = Z_hex.imag / (2 * np.pi * freq)
    rad.CplMagDelete(solver_hex)

    print()
    print(f"{'':20} {'ObjRecMag':>15} {'ObjHexahedron':>15} {'Diff [%]':>10}")
    print("-"*65)
    print(f"{'L [nH]':20} {L_rec*1e9:>15.2f} {L_hex*1e9:>15.2f} {abs(L_rec-L_hex)/L_rec*100:>10.2f}")
    print(f"{'R [mOhm]':20} {Z_rec.real*1000:>15.4f} {Z_hex.real*1000:>15.4f} {abs(Z_rec.real-Z_hex.real)/Z_rec.real*100:>10.2f}")
    print(f"{'Ratio (L/L_air)':20} {L_rec/L_air:>15.2f} {L_hex/L_air:>15.2f}")
    print()


def test_different_aspect_ratios():
    """Test with different core aspect ratios using ObjHexahedron."""
    print("\n" + "="*60)
    print("Test: Different Aspect Ratios (ObjHexahedron)")
    print("="*60)

    loop_radius = 0.05  # 50 mm
    freq = 1000  # 1 kHz
    mu_r = 1000

    # Get baseline (air core) inductance
    rad.UtiDelAll()
    rad.FldUnits('m')
    coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)
    rad.CndSetFrequency(coil, freq)
    rad.CndSetVoltage(coil, 1.0, 0.0)
    rad.CndSolve(coil)
    Z_air = rad.CndGetImpedance(coil)
    L_air = Z_air.imag / (2 * np.pi * freq)
    print(f"Air core L = {L_air * 1e9:.2f} nH")
    print()

    # Different aspect ratios (same volume)
    base_volume = 0.03**3  # 27e-6 m^3

    aspect_ratios = [
        ('Cube', [0.03, 0.03, 0.03]),
        ('Flat disk', [0.06, 0.06, 0.0075]),
        ('Tall rod', [0.015, 0.015, 0.12]),
        ('Wide plate', [0.09, 0.03, 0.01]),
    ]

    print(f"{'Shape':15} {'Size [mm]':>25} {'L [nH]':>12} {'Ratio':>10}")
    print("-"*65)

    for name, size in aspect_ratios:
        rad.UtiDelAll()
        rad.FldUnits('m')

        coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)
        core = create_hexahedron_core([0, 0, 0], size)
        mat = rad.MatLin(mu_r)
        rad.MatApl(core, mat)

        solver = rad.CplMagCreate(coil, core)
        rad.CplMagSetFrequency(solver, freq)
        rad.CplMagSetVoltage(solver, 1.0, 0.0)
        rad.CplMagSetMu(solver, mu_r, 0)

        result = rad.CplMagSolve(solver)
        Z = result['Z']
        L = Z.imag / (2 * np.pi * freq)
        ratio = L / L_air

        size_str = f"{size[0]*1000:.0f}x{size[1]*1000:.0f}x{size[2]*1000:.1f}"
        print(f"{name:15} {size_str:>25} {L*1e9:>12.2f} {ratio:>10.2f}")

        rad.CplMagDelete(solver)

    print()


def test_multiple_cores():
    """Test with multiple magnetic cores (if container support exists)."""
    print("\n" + "="*60)
    print("Test: Single vs Multiple Cores")
    print("="*60)
    print("Note: Multiple cores require ObjCnt support in CplMag")
    print("      (Currently testing single core only)")

    # For now, just test that single core works
    rad.UtiDelAll()
    rad.FldUnits('m')

    loop_radius = 0.05
    freq = 1000
    mu_r = 1000

    coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)

    # Single larger core
    core = create_hexahedron_core([0, 0, 0], 0.03)
    mat = rad.MatLin(mu_r)
    rad.MatApl(core, mat)

    solver = rad.CplMagCreate(coil, core)
    rad.CplMagSetFrequency(solver, freq)
    rad.CplMagSetVoltage(solver, 1.0, 0.0)
    rad.CplMagSetMu(solver, mu_r, 0)

    result = rad.CplMagSolve(solver)
    Z = result['Z']
    L = Z.imag / (2 * np.pi * freq)

    print(f"Single 30mm core: L = {L * 1e9:.2f} nH")

    rad.CplMagDelete(solver)
    print()


def main():
    """Run all hexahedron tests."""
    print("="*60)
    print("CplMag Solver - ObjHexahedron Tests")
    print("="*60)

    test_hexahedron_vs_recmag()
    test_different_aspect_ratios()
    test_multiple_cores()

    print("="*60)
    print("All ObjHexahedron tests completed!")
    print("="*60)


if __name__ == '__main__':
    main()
