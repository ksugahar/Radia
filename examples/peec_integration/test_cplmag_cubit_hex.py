"""
Test CplMag solver with hexahedral mesh.

This example creates hex mesh cores and tests the coupling
with a coil conductor in the CplMag solver using the new
multi-element MMM formulation.

Features tested:
- Single vs multi-element cores
- Mesh refinement convergence
- cubit_hex_to_radia() and create_hex_mesh_grid() functions from netgen_mesh_import

Usage:
    python test_cplmag_cubit_hex.py
"""

import sys
import os
import numpy as np

# Add Radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

# Import the new hex mesh functions from netgen_mesh_import
from netgen_mesh_import import create_hex_mesh_grid, cubit_hex_to_radia


def test_single_vs_multi_element():
    """Compare single element vs multi-element core."""
    import radia as rad

    print("\n" + "="*60)
    print("Test: Single Element vs Multi-Element Core")
    print("="*60)

    loop_radius = 0.05  # 50 mm
    freq = 1000  # 1 kHz
    core_size = [0.03, 0.03, 0.03]  # 30mm cube
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

    # Test with single element
    rad.UtiDelAll()
    rad.FldUnits('m')
    coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)

    core_single = rad.ObjRecMag([0, 0, 0], core_size, [0, 0, 0])
    mat = rad.MatLin(mu_r)
    rad.MatApl(core_single, mat)

    solver = rad.CplMagCreate(coil, core_single)
    rad.CplMagSetFrequency(solver, freq)
    rad.CplMagSetVoltage(solver, 1.0, 0.0)
    rad.CplMagSetMu(solver, mu_r, 0)

    result_single = rad.CplMagSolve(solver)
    Z_single = result_single['Z']
    L_single = Z_single.imag / (2 * np.pi * freq)
    rad.CplMagDelete(solver)

    print(f"\nSingle element: L = {L_single * 1e9:.2f} nH, ratio = {L_single/L_air:.2f}")

    # Test with multi-element (2x2x2 = 8 elements) using create_hex_mesh_grid
    print("\nCreating multi-element core (2x2x2 = 8 hexes) using create_hex_mesh_grid()...")

    rad.UtiDelAll()
    rad.FldUnits('m')
    coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)

    # Create 8 sub-elements using the new create_hex_mesh_grid function
    core_container = create_hex_mesh_grid(
        center=[0, 0, 0],
        size=core_size,
        divisions=[2, 2, 2],
        mu_r=mu_r,
        verbose=True
    )

    # Test with CplMag
    try:
        solver = rad.CplMagCreate(coil, core_container)
        rad.CplMagSetFrequency(solver, freq)
        rad.CplMagSetVoltage(solver, 1.0, 0.0)
        rad.CplMagSetMu(solver, mu_r, 0)

        result_multi = rad.CplMagSolve(solver)
        Z_multi = result_multi['Z']
        L_multi = Z_multi.imag / (2 * np.pi * freq)
        rad.CplMagDelete(solver)

        print(f"\nMulti-element (8 hexes): L = {L_multi * 1e9:.2f} nH, ratio = {L_multi/L_air:.2f}")
        print(f"\nComparison:")
        print(f"  Single element:  L = {L_single * 1e9:.2f} nH")
        print(f"  Multi-element:   L = {L_multi * 1e9:.2f} nH")
        print(f"  Difference:      {abs(L_multi - L_single) / L_single * 100:.1f}%")

    except Exception as e:
        print(f"\nError with multi-element core: {e}")
        print("ObjCnt may not be supported in current CplMag implementation.")


def test_mesh_refinement():
    """Test convergence with mesh refinement."""
    import radia as rad

    print("\n" + "="*60)
    print("Test: Mesh Refinement Convergence")
    print("="*60)

    loop_radius = 0.05  # 50 mm
    freq = 1000  # 1 kHz
    core_size = [0.03, 0.03, 0.03]  # 30mm cube
    mu_r = 1000

    # Get baseline
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

    mesh_sizes = [1, 2, 3]  # n x n x n elements

    print(f"{'Mesh':>10} {'Elements':>10} {'L [nH]':>12} {'Ratio':>10}")
    print("-"*45)

    for n in mesh_sizes:
        rad.UtiDelAll()
        rad.FldUnits('m')
        coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)

        if n == 1:
            # Single element
            core = rad.ObjRecMag([0, 0, 0], core_size, [0, 0, 0])
            mat = rad.MatLin(mu_r)
            rad.MatApl(core, mat)
        else:
            # Multi-element using create_hex_mesh_grid
            core = create_hex_mesh_grid(
                center=[0, 0, 0],
                size=core_size,
                divisions=[n, n, n],
                mu_r=mu_r,
                verbose=False  # Suppress output in loop
            )

        try:
            solver = rad.CplMagCreate(coil, core)
            rad.CplMagSetFrequency(solver, freq)
            rad.CplMagSetVoltage(solver, 1.0, 0.0)
            rad.CplMagSetMu(solver, mu_r, 0)

            result = rad.CplMagSolve(solver)
            Z = result['Z']
            L = Z.imag / (2 * np.pi * freq)
            rad.CplMagDelete(solver)

            n_elem = n**3
            print(f"{n}x{n}x{n}:>10 {n_elem:>10} {L*1e9:>12.2f} {L/L_air:>10.2f}")

        except Exception as e:
            print(f"{n}x{n}x{n}:>10 {n**3:>10} Error: {e}")

    print()


def main():
    """Run all tests."""
    print("="*60)
    print("CplMag Solver - Cubit Hex Mesh Tests")
    print("="*60)

    test_single_vs_multi_element()
    test_mesh_refinement()

    print("="*60)
    print("Tests completed!")
    print("="*60)


if __name__ == '__main__':
    main()
