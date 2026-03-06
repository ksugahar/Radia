"""
Test PyKAN Surface Impedance Learning

Tests the PyKANSurfaceImpedance class for learning nonlinear surface impedance
Z_s(H, f) from ESIM-like data.

Tests:
1. Basic training from synthetic ESIM data
2. Prediction accuracy validation
3. SPICE PWL export format
4. Edge cases and error handling
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import numpy as np


def test_pykan_surface_impedance():
    """Test PyKAN surface impedance learning."""

    # Check if dependencies are available
    try:
        import torch
        from kan import KAN
        PYKAN_AVAILABLE = True
    except ImportError as e:
        print(f"PyKAN or PyTorch not available: {e}")
        print("Install with: pip install pykan torch")
        return {'passed': True, 'skipped': True, 'reason': 'PyKAN not installed'}

    from lanczos_reduction import PyKANSurfaceImpedance

    results = {
        'test_1_training': False,
        'test_2_prediction': False,
        'test_3_spice_export': False,
        'test_4_edge_cases': False,
    }

    # =========================================================================
    # Test 1: Basic Training
    # =========================================================================
    print("=" * 60)
    print("Test 1: Basic KAN Training from Synthetic ESIM Data")
    print("=" * 60)

    # Generate synthetic ESIM data
    # Z_s = (1+j) * sqrt(omega * mu / (2 * sigma)) * (1 + alpha * H^2)
    # where alpha represents saturation effect

    sigma = 1e7  # S/m (copper-like)
    mu0 = 4 * np.pi * 1e-7
    mu_r = 100  # Relative permeability
    alpha = 1e-8  # Saturation coefficient

    H_values = np.array([10, 50, 100, 200, 500, 1000, 2000, 5000])  # A/m
    f_values = np.array([1e3, 5e3, 10e3, 50e3, 100e3, 500e3])  # Hz

    n_H = len(H_values)
    n_f = len(f_values)
    Zs_data = np.zeros((n_H, n_f), dtype=complex)

    for i, H in enumerate(H_values):
        for j, f in enumerate(f_values):
            omega = 2 * np.pi * f
            # Skin effect impedance with saturation
            mu_eff = mu0 * mu_r / (1 + alpha * H**2)  # Saturation
            delta = np.sqrt(2 / (omega * mu_eff * sigma))
            R_s = 1 / (sigma * delta)
            X_s = R_s  # Equal at skin depth limit
            Zs_data[i, j] = R_s + 1j * X_s

    # Create and train KAN model
    kan_model = PyKANSurfaceImpedance(width=[2, 8, 8, 2], grid=5, k=3)

    try:
        train_result = kan_model.train_from_esim_data(
            H_values, f_values, Zs_data,
            epochs=200, lr=0.02, verbose=True
        )

        if train_result['final_loss'] < 0.01:
            print(f"  Training passed: final_loss = {train_result['final_loss']:.4e}")
            results['test_1_training'] = True
        else:
            print(f"  Training loss too high: {train_result['final_loss']:.4e}")
    except Exception as e:
        print(f"  Training failed: {e}")
        import traceback
        traceback.print_exc()

    # =========================================================================
    # Test 2: Prediction Accuracy
    # =========================================================================
    print("\n" + "=" * 60)
    print("Test 2: Prediction Accuracy")
    print("=" * 60)

    if results['test_1_training']:
        # Test at interpolation points
        test_cases = [
            (100, 10e3),   # Within training range
            (500, 50e3),   # Within training range
            (1000, 100e3), # Within training range
        ]

        max_error = 0.0
        all_passed = True

        for H_test, f_test in test_cases:
            # True value (from formula)
            omega = 2 * np.pi * f_test
            mu_eff = mu0 * mu_r / (1 + alpha * H_test**2)
            delta = np.sqrt(2 / (omega * mu_eff * sigma))
            R_s_true = 1 / (sigma * delta)
            X_s_true = R_s_true
            Zs_true = R_s_true + 1j * X_s_true

            # KAN prediction
            Zs_pred = kan_model.predict(H_test, f_test)

            # Error
            error = abs(Zs_pred - Zs_true) / abs(Zs_true) * 100
            max_error = max(max_error, error)

            print(f"  H={H_test:5.0f} A/m, f={f_test/1e3:6.1f} kHz:")
            print(f"    True:  {Zs_true.real:.4e}+{Zs_true.imag:.4e}j")
            print(f"    KAN:   {Zs_pred.real:.4e}+{Zs_pred.imag:.4e}j")
            print(f"    Error: {error:.2f}%")

            if error > 10.0:  # 10% threshold
                all_passed = False

        if all_passed and max_error < 10.0:
            print(f"  Prediction test passed: max_error = {max_error:.2f}%")
            results['test_2_prediction'] = True
        else:
            print(f"  Prediction test failed: max_error = {max_error:.2f}%")
    else:
        print("  Skipped (training failed)")

    # =========================================================================
    # Test 3: SPICE PWL Export
    # =========================================================================
    print("\n" + "=" * 60)
    print("Test 3: SPICE PWL Export")
    print("=" * 60)

    if results['test_1_training']:
        try:
            # Generate SPICE PWL
            H_export = np.linspace(10, 5000, 20)
            f_ref = 100e3  # Reference frequency

            spice_pwl = kan_model.to_spice_pwl(H_export, f_ref, element_name="Zskin")

            # Check format
            lines = spice_pwl.strip().split('\n')
            has_subckt = any('.SUBCKT' in line.upper() for line in lines)
            has_ends = any('.ENDS' in line.upper() for line in lines)
            has_pwl = any('PWL' in line.upper() for line in lines)

            print("  SPICE PWL output (first 15 lines):")
            for line in lines[:15]:
                print(f"    {line}")
            if len(lines) > 15:
                print(f"    ... ({len(lines) - 15} more lines)")

            if has_subckt and has_ends:
                print(f"  SPICE export passed: {len(lines)} lines generated")
                results['test_3_spice_export'] = True
            else:
                print(f"  SPICE export incomplete: SUBCKT={has_subckt}, ENDS={has_ends}")

        except Exception as e:
            print(f"  SPICE export failed: {e}")
            import traceback
            traceback.print_exc()
    else:
        print("  Skipped (training failed)")

    # =========================================================================
    # Test 4: Edge Cases
    # =========================================================================
    print("\n" + "=" * 60)
    print("Test 4: Edge Cases and Error Handling")
    print("=" * 60)

    edge_passed = True

    # Test 4a: Prediction before training
    try:
        untrained_model = PyKANSurfaceImpedance()
        untrained_model.predict(100, 10e3)
        print("  4a: Should have raised error for untrained model")
        edge_passed = False
    except ValueError as e:
        print(f"  4a: Correctly raised ValueError: {e}")

    # Test 4b: Single H value
    try:
        if results['test_1_training']:
            single_pred = kan_model.predict(100.0, 10000.0)
            if isinstance(single_pred, complex):
                print(f"  4b: Single prediction works: {single_pred}")
            else:
                print(f"  4b: Single prediction wrong type: {type(single_pred)}")
                edge_passed = False
    except Exception as e:
        print(f"  4b: Single prediction failed: {e}")
        edge_passed = False

    # Test 4c: Invalid input dimensions
    try:
        bad_model = PyKANSurfaceImpedance()
        bad_model.train_from_esim_data(
            np.array([1, 2]),  # Too few H values
            np.array([1000, 2000]),
            np.zeros((2, 2), dtype=complex),
            epochs=1, verbose=False
        )
        print("  4c: Training with minimal data completed")
    except Exception as e:
        print(f"  4c: Training with minimal data failed (expected): {e}")

    results['test_4_edge_cases'] = edge_passed

    # =========================================================================
    # Summary
    # =========================================================================
    print("\n" + "=" * 60)
    print("Test Summary")
    print("=" * 60)

    all_passed = all(results.values())
    for test_name, passed in results.items():
        status = "PASSED" if passed else "FAILED"
        print(f"  {test_name}: {status}")

    print("-" * 60)
    overall = "ALL TESTS PASSED" if all_passed else "SOME TESTS FAILED"
    print(f"  Overall: {overall}")

    return {'passed': all_passed, 'skipped': False, 'results': results}


def test_verilog_a_export():
    """Test Verilog-A export functionality."""
    print("\n" + "=" * 60)
    print("Bonus Test: Verilog-A Export")
    print("=" * 60)

    try:
        import torch
        from kan import KAN
    except ImportError:
        print("  Skipped: PyKAN not available")
        return

    from lanczos_reduction import PyKANSurfaceImpedance

    # Quick training
    sigma = 1e7
    mu0 = 4 * np.pi * 1e-7
    mu_r = 100

    H_values = np.array([10, 100, 1000])
    f_values = np.array([1e3, 100e3])
    Zs_data = np.zeros((3, 2), dtype=complex)

    for i, H in enumerate(H_values):
        for j, f in enumerate(f_values):
            omega = 2 * np.pi * f
            delta = np.sqrt(2 / (omega * mu0 * mu_r * sigma))
            R_s = 1 / (sigma * delta)
            Zs_data[i, j] = R_s * (1 + 1j)

    model = PyKANSurfaceImpedance(width=[2, 4, 2], grid=3)
    model.train_from_esim_data(H_values, f_values, Zs_data, epochs=50, verbose=False)

    # Generate Verilog-A
    verilog_code = model.to_verilog_a(name="kan_zskin_test")

    lines = verilog_code.strip().split('\n')
    print("  Verilog-A output (first 20 lines):")
    for line in lines[:20]:
        print(f"    {line}")
    if len(lines) > 20:
        print(f"    ... ({len(lines) - 20} more lines)")

    # Basic checks
    has_module = any('module' in line.lower() for line in lines)
    has_endmodule = any('endmodule' in line.lower() for line in lines)

    if has_module and has_endmodule:
        print(f"  Verilog-A export passed: {len(lines)} lines")
    else:
        print("  Verilog-A export incomplete")


if __name__ == '__main__':
    print("PyKAN Surface Impedance Test Suite")
    print("=" * 60)

    result = test_pykan_surface_impedance()

    if not result.get('skipped', False):
        test_verilog_a_export()

    print("\n" + "=" * 60)
    if result['passed']:
        print("Test suite completed successfully!")
    else:
        print("Test suite completed with failures.")
        sys.exit(1)
