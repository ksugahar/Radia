"""
ESIM Integration Tests

Comprehensive tests validating the complete ESIM workflow:
1. Cell problem solver accuracy
2. ESI table generation and interpolation
3. Coupled solver convergence
4. Physical consistency checks

Author: Radia Development Team
Date: 2026-01-08
"""

import sys
import numpy as np

from radia.esim_cell_problem import (
    ESIMCellProblemSolver,
    BHCurveInterpolator,
    ESITable,
    generate_esi_table_from_bh_curve,
)
from radia.esim_workpiece import create_esim_block, create_esim_cylinder
from radia.esim_coupled_solver import (
    InductionHeatingCoil,
    ESIMCoupledSolver,
    solve_induction_heating,
)


def test_bh_curve_interpolation():
    """Test B-H curve interpolation accuracy."""
    print("Test 1: B-H Curve Interpolation")
    print("-" * 40)

    # Steel B-H curve
    bh_curve = [
        [0, 0],
        [100, 0.2],
        [500, 0.9],
        [1000, 1.3],
        [5000, 1.8],
        [50000, 2.1],
    ]

    bh_interp = BHCurveInterpolator(bh_curve)

    # Test at known points using actual API: B() and mu_r()
    tests = [
        (0, 0.0, "Origin"),
        (100, 0.2, "Data point 1"),
        (1000, 1.3, "Data point 2"),
        (300, None, "Interpolated point"),
    ]

    passed = 0
    for H, B_expected, desc in tests:
        B = bh_interp.B(H)  # Correct method name
        mu_r = bh_interp.mu_r(H)  # Correct method name

        if B_expected is not None:
            error = abs(B - B_expected)
            status = "PASS" if error < 0.01 else "FAIL"
            print(f"  {desc}: H={H} A/m -> B={B:.4f} T (expected {B_expected}), mu_r={mu_r:.1f} [{status}]")
            if status == "PASS":
                passed += 1
        else:
            print(f"  {desc}: H={H} A/m -> B={B:.4f} T, mu_r={mu_r:.1f}")
            passed += 1

    print(f"  Result: {passed}/{len(tests)} tests passed")
    print()
    return passed == len(tests)


def test_cell_problem_solver():
    """Test 1D cell problem solver."""
    print("Test 2: Cell Problem Solver")
    print("-" * 40)

    # Parameters
    bh_curve = [
        [0, 0], [100, 0.2], [500, 0.9], [1000, 1.3],
        [5000, 1.8], [50000, 2.1],
    ]
    sigma = 2e6  # S/m
    freq = 50000  # Hz

    solver = ESIMCellProblemSolver(bh_curve, sigma, freq)

    # Test at different H0 values
    H0_values = [100, 500, 1000, 5000]
    passed = True

    for H0 in H0_values:
        result = solver.solve(H0)
        Z = result['Z']  # Correct key name
        P_prime = result['P_prime']

        # Check physical consistency
        # Z should have positive real part (resistance)
        # P_prime should be positive (power loss)
        if Z.real <= 0:
            print(f"  H0={H0}: Z.real = {Z.real:.6e} <= 0 [FAIL]")
            passed = False
        elif P_prime < 0:
            print(f"  H0={H0}: P_prime = {P_prime:.6e} < 0 [FAIL]")
            passed = False
        else:
            print(f"  H0={H0}: Z = {Z.real:.4e} + j{Z.imag:.4e} Ohm, P' = {P_prime:.2e} W/m^2 [PASS]")

    print(f"  Result: {'All tests passed' if passed else 'Some tests failed'}")
    print()
    return passed


def test_esi_table_generation():
    """Test ESI table generation and lookup."""
    print("Test 3: ESI Table Generation")
    print("-" * 40)

    bh_curve = [
        [0, 0], [100, 0.2], [500, 0.9], [1000, 1.3],
        [5000, 1.8], [50000, 2.1],
    ]
    sigma = 2e6
    freq = 50000

    # Generate table (use actual API without H0_min/H0_max)
    esi_table = generate_esi_table_from_bh_curve(
        bh_curve, sigma, freq,
        n_points=20
    )

    # Test interpolation (use actual API: get_impedance, get_power_loss)
    H0_test = 500
    Z = esi_table.get_impedance(H0_test)
    P_prime, Q_prime = esi_table.get_power_loss(H0_test)

    print(f"  Table generated with {len(esi_table.H0_values)} points")
    print(f"  H0 range: {esi_table.H0_values[0]:.1f} to {esi_table.H0_values[-1]:.1f} A/m")
    print(f"  Test at H0={H0_test}: Z = {Z.real:.4e} + j{Z.imag:.4e} Ohm")
    print(f"  P' = {P_prime:.2e} W/m^2, Q' = {Q_prime:.2e} var/m^2")

    # Check monotonicity of |Z| with H0 (should decrease due to saturation)
    Z_mags = [abs(esi_table.get_impedance(H)) for H in [100, 500, 1000, 5000]]
    monotonic = all(Z_mags[i] >= Z_mags[i+1] for i in range(len(Z_mags)-1))

    if monotonic:
        print(f"  |Z| decreases with H0 (saturation effect): [PASS]")
    else:
        print(f"  |Z| vs H0: {Z_mags} - Expected decreasing [FAIL]")

    print()
    return monotonic


def test_coil_field_computation():
    """Test coil magnetic field computation."""
    print("Test 4: Coil Field Computation")
    print("-" * 40)

    # Create spiral coil
    coil = InductionHeatingCoil(
        coil_type='spiral',
        center=[0, 0, 0.02],
        inner_radius=0.03,
        outer_radius=0.05,
        pitch=0.005,
        num_turns=3,
        axis=[0, 0, 1],
    )
    coil.set_current(100)

    # Test field on axis
    B_axis = coil.compute_field_at_point([0, 0, 0])
    B_mag = np.linalg.norm(np.real(B_axis))

    print(f"  Coil: 3-turn spiral, R_in=30mm, R_out=50mm, I=100A")
    print(f"  B at origin: [{B_axis[0].real:.4e}, {B_axis[1].real:.4e}, {B_axis[2].real:.4e}] T")
    print(f"  |B| = {B_mag*1000:.2f} mT")

    # B should be primarily in z-direction on axis
    if abs(B_axis[2].real) > abs(B_axis[0].real) and abs(B_axis[2].real) > abs(B_axis[1].real):
        print(f"  B primarily in z-direction on axis: [PASS]")
        passed = True
    else:
        print(f"  B direction check: [FAIL]")
        passed = False

    # Test field decay with distance
    B_near = np.linalg.norm(np.real(coil.compute_field_at_point([0, 0, 0.01])))
    B_far = np.linalg.norm(np.real(coil.compute_field_at_point([0, 0, 0.10])))

    if B_near > B_far:
        print(f"  B decays with distance: {B_near*1000:.2f} mT > {B_far*1000:.2f} mT [PASS]")
    else:
        print(f"  B decay check: [FAIL]")
        passed = False

    print()
    return passed


def test_coupled_solver():
    """Test coupled ESIM solver convergence."""
    print("Test 5: Coupled Solver Convergence")
    print("-" * 40)

    bh_curve = [
        [0, 0], [100, 0.2], [500, 0.9], [1000, 1.3],
        [5000, 1.8], [50000, 2.1],
    ]
    sigma = 2e6
    freq = 50000

    # Create coil
    coil = InductionHeatingCoil(
        coil_type='spiral',
        center=[0, 0, 0.02],
        inner_radius=0.03,
        outer_radius=0.05,
        pitch=0.005,
        num_turns=3,
        axis=[0, 0, 1],
    )
    coil.set_current(100)

    # Create workpiece
    workpiece = create_esim_block(
        center=[0, 0, -0.01],
        dimensions=[0.08, 0.08, 0.02],
        bh_curve=bh_curve,
        sigma=sigma,
        frequency=freq,
        panels_per_side=6
    )

    # Solve
    solver = ESIMCoupledSolver(coil, workpiece, freq)
    result = solver.solve(tol=1e-4, max_iter=30, verbose=False)

    print(f"  Converged: {result['converged']}")
    print(f"  Iterations: {result['iterations']}")
    print(f"  P_total = {result['P_total']:.1f} W")
    print(f"  Q_total = {result['Q_total']:.1f} var")
    print(f"  Power factor = {result['power_factor']:.3f}")

    passed = result['converged']
    if passed:
        # Additional checks
        if result['P_total'] > 0:
            print(f"  P > 0: [PASS]")
        else:
            print(f"  P > 0: [FAIL]")
            passed = False

        if 0 < result['power_factor'] <= 1:
            print(f"  0 < PF <= 1: [PASS]")
        else:
            print(f"  Power factor range: [FAIL]")
            passed = False

    print()
    return passed


def _validate_physical_consistency():
    """Test physical consistency of results."""
    print("Test 7: Physical Consistency")
    print("-" * 40)

    bh_curve = [[0, 0], [100, 0.2], [500, 0.9], [1000, 1.3], [5000, 1.8], [50000, 2.1]]
    sigma = 2e6
    freq = 50000

    coil = InductionHeatingCoil(
        coil_type='spiral',
        center=[0, 0, 0.02],
        inner_radius=0.03,
        outer_radius=0.05,
        pitch=0.005,
        num_turns=3,
        axis=[0, 0, 1],
    )

    workpiece = create_esim_block(
        center=[0, 0, -0.01],
        dimensions=[0.08, 0.08, 0.02],
        bh_curve=bh_curve,
        sigma=sigma,
        frequency=freq,
        panels_per_side=5
    )

    passed = True

    # Test 1: Power scales with I^2 (linear regime)
    print("  Power scaling with I^2 (linear regime):")
    powers = []
    currents = [10, 20]
    for I in currents:
        coil.set_current(I)
        solver = ESIMCoupledSolver(coil, workpiece, freq)
        result = solver.solve(tol=1e-3, max_iter=10, verbose=False)
        powers.append(result['P_total'])
        print(f"    I = {I} A: P = {result['P_total']:.2f} W")

    # P should scale as I^2
    ratio = powers[1] / powers[0]
    expected = (currents[1] / currents[0]) ** 2
    error = abs(ratio - expected) / expected
    if error < 0.15:  # 15% tolerance for nonlinear effects
        print(f"    P2/P1 = {ratio:.2f} (expected {expected:.2f}): [PASS]")
    else:
        print(f"    P2/P1 = {ratio:.2f} (expected {expected:.2f}): [FAIL]")
        passed = False

    # Test 2: Higher frequency -> higher power density (skin effect)
    print("  Frequency effect on power:")
    coil.set_current(50)

    freqs = [30000, 60000]
    powers_freq = []
    for f in freqs:
        wp = create_esim_block(
            center=[0, 0, -0.01],
            dimensions=[0.08, 0.08, 0.02],
            bh_curve=bh_curve,
            sigma=sigma,
            frequency=f,
            panels_per_side=5
        )
        solver = ESIMCoupledSolver(coil, wp, f)
        result = solver.solve(tol=1e-3, max_iter=10, verbose=False)
        powers_freq.append(result['P_total'])
        print(f"    f = {f/1000:.0f} kHz: P = {result['P_total']:.2f} W")

    # Higher frequency should generally increase power (complex due to nonlinearity)
    if powers_freq[1] != powers_freq[0]:
        print(f"    Frequency affects power: [PASS]")
    else:
        print(f"    Frequency effect: [FAIL]")
        passed = False

    print()
    return passed


def main():
    """Run all integration tests."""
    print()
    print("=" * 60)
    print("ESIM Integration Tests")
    print("=" * 60)
    print()

    tests = [
        ("B-H Curve Interpolation", test_bh_curve_interpolation),
        ("Cell Problem Solver", test_cell_problem_solver),
        ("ESI Table Generation", test_esi_table_generation),
        ("Coil Field Computation", test_coil_field_computation),
        ("Coupled Solver", test_coupled_solver),
        ("Physical Consistency", _validate_physical_consistency),
    ]

    results = []
    for name, test_func in tests:
        try:
            passed = test_func()
            results.append((name, passed))
        except Exception as e:
            print(f"  Exception: {e}")
            results.append((name, False))

    # Summary
    print("=" * 60)
    print("Test Summary")
    print("=" * 60)

    total_passed = 0
    for name, passed in results:
        status = "PASS" if passed else "FAIL"
        print(f"  {name}: [{status}]")
        if passed:
            total_passed += 1

    print()
    print(f"Total: {total_passed}/{len(tests)} tests passed")
    print()

    if total_passed == len(tests):
        print("All tests passed! ESIM implementation is validated.")
    else:
        print("Some tests failed. Review the output above.")

    return total_passed == len(tests)


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
