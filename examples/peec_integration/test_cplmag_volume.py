"""
Test CplMag: Verify volume dependence (single dipole approximation).

Since each magnetic object is treated as a single dipole,
the coupling should depend on volume, not shape.
"""

import sys
sys.path.insert(0, '../../src/radia')
import radia as rad
import numpy as np

def test_same_volume_different_shape():
    """Same volume, different shapes - should give SAME L."""
    print("\n" + "="*60)
    print("Test: Same Volume, Different Shapes")
    print("(Single dipole: volume matters, shape doesn't)")
    print("="*60)

    loop_radius = 0.05  # 50 mm
    freq = 1000
    mu_r = 1000

    # Fixed volume = 27e-6 m^3 (30mm cube)
    base_volume = 0.03**3

    shapes = [
        ('30mm cube', [0.03, 0.03, 0.03]),                    # V = 27e-6
        ('Flat (60x60x7.5)', [0.06, 0.06, 0.0075]),           # V = 27e-6
        ('Tall (15x15x120)', [0.015, 0.015, 0.12]),           # V = 27e-6
    ]

    print(f"\nAll shapes have volume = {base_volume*1e6:.1f} mm^3")
    print()
    print(f"{'Shape':25} {'Actual V [mm^3]':>15} {'L [nH]':>12}")
    print("-"*55)

    for name, size in shapes:
        rad.UtiDelAll()
        rad.FldUnits('m')

        actual_vol = size[0] * size[1] * size[2]

        coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)

        # Create hexahedron
        cx, cy, cz = 0, 0, 0
        sx, sy, sz = size
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
        core = rad.ObjHexahedron(vertices, [0, 0, 0])
        mat = rad.MatLin(mu_r)
        rad.MatApl(core, mat)

        solver = rad.CplMagCreate(coil, core)
        rad.CplMagSetFrequency(solver, freq)
        rad.CplMagSetVoltage(solver, 1.0, 0.0)
        rad.CplMagSetMu(solver, mu_r, 0)

        result = rad.CplMagSolve(solver)
        Z = result['Z']
        L = Z.imag / (2 * np.pi * freq)

        print(f"{name:25} {actual_vol*1e9:>15.3f} {L*1e9:>12.2f}")

        rad.CplMagDelete(solver)

    print()
    print("Expected: All L values should be EQUAL (same volume)")
    print()


def test_different_volume_same_shape():
    """Different volumes, same shape - should give DIFFERENT L."""
    print("\n" + "="*60)
    print("Test: Different Volumes, Same Shape (cube)")
    print("(L should scale with volume)")
    print("="*60)

    loop_radius = 0.05  # 50 mm
    freq = 1000
    mu_r = 1000

    cube_sizes = [0.02, 0.025, 0.03, 0.035, 0.04]  # 20mm to 40mm cube

    print()
    print(f"{'Size [mm]':>12} {'V [mm^3]':>12} {'L [nH]':>12} {'L/V ratio':>15}")
    print("-"*55)

    results = []
    for size in cube_sizes:
        rad.UtiDelAll()
        rad.FldUnits('m')

        vol = size**3

        coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)
        core = rad.ObjRecMag([0, 0, 0], [size, size, size], [0, 0, 0])
        mat = rad.MatLin(mu_r)
        rad.MatApl(core, mat)

        solver = rad.CplMagCreate(coil, core)
        rad.CplMagSetFrequency(solver, freq)
        rad.CplMagSetVoltage(solver, 1.0, 0.0)
        rad.CplMagSetMu(solver, mu_r, 0)

        result = rad.CplMagSolve(solver)
        Z = result['Z']
        L = Z.imag / (2 * np.pi * freq)

        L_over_V = L / vol
        results.append((size, vol, L, L_over_V))

        print(f"{size*1000:>12.0f} {vol*1e9:>12.3f} {L*1e9:>12.2f} {L_over_V:>15.2e}")

        rad.CplMagDelete(solver)

    print()
    # Check if L/V is approximately constant
    L_over_V_values = [r[3] for r in results]
    ratio_variation = (max(L_over_V_values) - min(L_over_V_values)) / np.mean(L_over_V_values) * 100
    print(f"L/V ratio variation: {ratio_variation:.1f}%")
    if ratio_variation < 10:
        print("-> L scales approximately linearly with V (single dipole model)")
    else:
        print("-> L does NOT scale linearly with V")
    print()


def main():
    """Run volume dependence tests."""
    print("="*60)
    print("CplMag Solver - Volume Dependence Tests")
    print("="*60)
    print()
    print("NOTE: Current implementation treats each magnetic object")
    print("      as a SINGLE DIPOLE. Therefore:")
    print("      - Same volume, different shape -> same L")
    print("      - Different volume, same shape -> different L")
    print()

    test_same_volume_different_shape()
    test_different_volume_same_shape()

    print("="*60)
    print("Volume dependence tests completed!")
    print("="*60)


if __name__ == '__main__':
    main()
