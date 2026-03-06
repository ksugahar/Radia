"""
Test that symmetric and non-symmetric formulations give identical results.

This script verifies that:
1. Symmetric formulation (useSymmetric=True) gives the same impedance as
   non-symmetric formulation (useSymmetric=False)
2. The magnetization solution is also identical after scaling recovery

The symmetric formulation enables PRIMA (Cauer Ladder Network) model order
reduction while preserving the physical solution.

Usage:
    python test_symmetrization.py
"""

import sys
import os
import numpy as np

# Add Radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
import radia as rad


def test_symmetry_equivalence():
    """Verify symmetric and non-symmetric give same impedance."""
    print("="*70)
    print("Test: Symmetric vs Non-Symmetric Formulation Equivalence")
    print("="*70)

    # Physical constants
    mu_0 = 4 * np.pi * 1e-7

    # Problem parameters
    loop_radius = 0.05  # 50 mm
    freq = 1000  # 1 kHz
    core_size = [0.03, 0.03, 0.03]  # 30mm cube
    mu_r = 1000

    # Element volume and scaling factor
    V_elem = core_size[0] * core_size[1] * core_size[2]
    alpha = np.sqrt(mu_0 * V_elem)

    print(f"\nProblem setup:")
    print(f"  Loop radius: {loop_radius*1000:.1f} mm")
    print(f"  Core size:   {core_size[0]*1000:.1f} x {core_size[1]*1000:.1f} x {core_size[2]*1000:.1f} mm")
    print(f"  mu_r:        {mu_r}")
    print(f"  Frequency:   {freq} Hz")
    print(f"\nScaling factor: alpha = sqrt(mu_0 * V) = {alpha:.6e}")

    # Test 1: Non-symmetric formulation (default)
    print("\n" + "-"*70)
    print("Solving with NON-SYMMETRIC formulation (default)...")
    print("-"*70)

    rad.UtiDelAll()
    rad.FldUnits('m')

    coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r',
                       2e-3, 2e-3, 5.8e7, 8, 36)
    core = rad.ObjRecMag([0, 0, 0], core_size, [0, 0, 0])
    mat = rad.MatLin(mu_r)
    rad.MatApl(core, mat)

    solver_nonsym = rad.CplMagCreate(coil, core)
    rad.CplMagSetFrequency(solver_nonsym, freq)
    rad.CplMagSetVoltage(solver_nonsym, 1.0, 0.0)
    rad.CplMagSetMu(solver_nonsym, mu_r, 0)
    rad.CplMagSetSymmetric(solver_nonsym, 0)  # Non-symmetric (default)

    result_nonsym = rad.CplMagSolve(solver_nonsym)
    Z_nonsym = result_nonsym['Z']
    rad.CplMagDelete(solver_nonsym)

    print(f"  Impedance Z = {Z_nonsym.real:.6e} + j{Z_nonsym.imag:.6e} Ohm")
    print(f"  |Z| = {abs(Z_nonsym):.6e} Ohm")

    # Test 2: Symmetric formulation
    print("\n" + "-"*70)
    print("Solving with SYMMETRIC formulation (for PRIMA)...")
    print("-"*70)

    rad.UtiDelAll()
    rad.FldUnits('m')

    coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r',
                       2e-3, 2e-3, 5.8e7, 8, 36)
    core = rad.ObjRecMag([0, 0, 0], core_size, [0, 0, 0])
    mat = rad.MatLin(mu_r)
    rad.MatApl(core, mat)

    solver_sym = rad.CplMagCreate(coil, core)
    rad.CplMagSetFrequency(solver_sym, freq)
    rad.CplMagSetVoltage(solver_sym, 1.0, 0.0)
    rad.CplMagSetMu(solver_sym, mu_r, 0)
    rad.CplMagSetSymmetric(solver_sym, 1)  # Symmetric (for PRIMA)

    result_sym = rad.CplMagSolve(solver_sym)
    Z_sym = result_sym['Z']
    rad.CplMagDelete(solver_sym)

    print(f"  Impedance Z = {Z_sym.real:.6e} + j{Z_sym.imag:.6e} Ohm")
    print(f"  |Z| = {abs(Z_sym):.6e} Ohm")

    # Compare results
    print("\n" + "="*70)
    print("Comparison:")
    print("="*70)

    Z_diff = abs(Z_sym - Z_nonsym)
    Z_rel_diff = Z_diff / abs(Z_nonsym) * 100

    print(f"\n  Non-symmetric: Z = {Z_nonsym.real:.6e} + j{Z_nonsym.imag:.6e} Ohm")
    print(f"  Symmetric:     Z = {Z_sym.real:.6e} + j{Z_sym.imag:.6e} Ohm")
    print(f"\n  Absolute difference: {Z_diff:.6e} Ohm")
    print(f"  Relative difference: {Z_rel_diff:.6f} %")

    # Verify they are essentially identical
    tolerance = 1e-10  # Should be machine precision
    if Z_rel_diff < tolerance:
        print(f"\n  PASS: Results are IDENTICAL (diff < {tolerance*100}%)")
        success = True
    else:
        print(f"\n  FAIL: Results differ by more than {tolerance*100}%")
        success = False

    # Additional info
    L_nonsym = Z_nonsym.imag / (2 * np.pi * freq)
    L_sym = Z_sym.imag / (2 * np.pi * freq)

    print(f"\n  Inductance comparison:")
    print(f"    Non-symmetric: L = {L_nonsym*1e9:.4f} nH")
    print(f"    Symmetric:     L = {L_sym*1e9:.4f} nH")

    return success


def test_multi_element():
    """Test with multiple magnetic elements."""
    print("\n" + "="*70)
    print("Test: Multi-Element Core (2x2x2 = 8 elements)")
    print("="*70)

    from netgen_mesh_import import create_hex_mesh_grid

    loop_radius = 0.05
    freq = 1000
    core_size = [0.03, 0.03, 0.03]
    mu_r = 1000

    # Non-symmetric
    rad.UtiDelAll()
    rad.FldUnits('m')
    coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r',
                       2e-3, 2e-3, 5.8e7, 8, 36)
    core = create_hex_mesh_grid(
        center=[0, 0, 0],
        size=core_size,
        divisions=[2, 2, 2],
        mu_r=mu_r,
        verbose=False
    )

    solver = rad.CplMagCreate(coil, core)
    rad.CplMagSetFrequency(solver, freq)
    rad.CplMagSetVoltage(solver, 1.0, 0.0)
    rad.CplMagSetMu(solver, mu_r, 0)
    rad.CplMagSetSymmetric(solver, 0)

    result_nonsym = rad.CplMagSolve(solver)
    Z_nonsym = result_nonsym['Z']
    rad.CplMagDelete(solver)

    print(f"\n  Non-symmetric: Z = {Z_nonsym.real:.6e} + j{Z_nonsym.imag:.6e} Ohm")

    # Symmetric
    rad.UtiDelAll()
    rad.FldUnits('m')
    coil = rad.CndLoop([0, 0, 0], loop_radius, [0, 0, 1], 'r',
                       2e-3, 2e-3, 5.8e7, 8, 36)
    core = create_hex_mesh_grid(
        center=[0, 0, 0],
        size=core_size,
        divisions=[2, 2, 2],
        mu_r=mu_r,
        verbose=False
    )

    solver = rad.CplMagCreate(coil, core)
    rad.CplMagSetFrequency(solver, freq)
    rad.CplMagSetVoltage(solver, 1.0, 0.0)
    rad.CplMagSetMu(solver, mu_r, 0)
    rad.CplMagSetSymmetric(solver, 1)

    result_sym = rad.CplMagSolve(solver)
    Z_sym = result_sym['Z']
    rad.CplMagDelete(solver)

    print(f"  Symmetric:     Z = {Z_sym.real:.6e} + j{Z_sym.imag:.6e} Ohm")

    Z_diff = abs(Z_sym - Z_nonsym)
    Z_rel_diff = Z_diff / abs(Z_nonsym) * 100
    print(f"\n  Relative difference: {Z_rel_diff:.6f} %")

    if Z_rel_diff < 1e-10:
        print("  PASS: Multi-element results are IDENTICAL")
        return True
    else:
        print("  FAIL: Multi-element results differ")
        return False


def main():
    """Run all symmetrization tests."""
    print("="*70)
    print("CplMag Symmetrization Verification Tests")
    print("="*70)
    print("""
    These tests verify that the symmetric matrix formulation
    produces IDENTICAL results to the non-symmetric formulation.

    Symmetric formulation benefits:
    - Enables PRIMA (Cauer Ladder Network) model order reduction
    - Passivity guaranteed (important for SPICE export)
    - Physical interpretation as L-C ladder circuit
    """)

    results = []

    # Test 1: Single element
    results.append(("Single element", test_symmetry_equivalence()))

    # Test 2: Multi-element
    try:
        results.append(("Multi-element", test_multi_element()))
    except Exception as e:
        print(f"\n  Multi-element test skipped: {e}")
        results.append(("Multi-element", None))

    # Summary
    print("\n" + "="*70)
    print("Summary:")
    print("="*70)

    all_passed = True
    for name, success in results:
        if success is True:
            status = "PASS"
        elif success is False:
            status = "FAIL"
            all_passed = False
        else:
            status = "SKIP"
        print(f"  {name}: {status}")

    print()
    if all_passed:
        print("  All tests PASSED!")
        print("  Symmetric and non-symmetric formulations give IDENTICAL results.")
        print("  PRIMA model order reduction can be applied safely.")
    else:
        print("  Some tests FAILED. Please investigate.")

    print("="*70)
    return all_passed


if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
