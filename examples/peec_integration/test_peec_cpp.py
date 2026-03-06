"""
C++ PEEC pybind11 module test

Tests the pybind11-based PEEC matrix construction module (peec_matrices.pyd).
"""

import numpy as np
import sys

# Add paths for the pybind11 module
sys.path.insert(0, r'S:\Radia\01_GitHub\src\radia')

# PRIMA library path (Python implementation for comparison)
sys.path.insert(0, r'W:\30_CauerLadderNetwork\2021_01_22_CauerI_to_CauerII\Python')


def test_cpp_wire():
    """Test C++ PEEC Wire Matrix Construction (pybind11)"""
    print("=" * 60)
    print("Test: C++ PEEC Wire Matrix Construction (pybind11)")
    print("=" * 60)

    try:
        # Import pybind11 module
        from peec_matrices import PEECBuilder, create_wire_peec

        # Wire parameters
        start = np.array([0, 0, 0])
        end = np.array([0.1, 0, 0])  # 10cm
        width = 1e-3  # 1mm
        height = 1e-3
        n_segments = 10
        sigma = 5.8e7

        # Method 1: Using PEECBuilder class
        builder = PEECBuilder()
        n_created = builder.create_wire(start, end, width, height, n_segments, sigma)

        print(f"\nWire: {start} -> {end}")
        print(f"Cross section: {width*1000:.1f} x {height*1000:.1f} mm")
        print(f"Segments created: {n_created}")
        print(f"n_loop: {builder.n_loop}")

        # Build matrices
        L, R, P, M_LS = builder.build(include_star=True)

        print(f"\nL matrix shape: {L.shape}")
        print(f"R vector shape: {R.shape}")
        print(f"P matrix shape: {P.shape if P is not None else 'None'}")
        print(f"M_LS matrix shape: {M_LS.shape if M_LS is not None else 'None'}")

        # L matrix statistics
        L_diag_sum = np.sum(np.diag(L))
        L_total = np.sum(L)
        print(f"\nL self-inductance sum: {L_diag_sum*1e9:.3f} nH")
        print(f"L total (with mutual): {L_total*1e9:.3f} nH")

        # R matrix statistics
        R_total = np.sum(R)
        R_expected = 0.1 / (sigma * width * height)
        print(f"\nR total: {R_total*1000:.3f} mOhm")
        print(f"R expected: {R_expected*1000:.3f} mOhm")
        print(f"R match: {abs(R_total - R_expected) / R_expected < 0.01}")

        # Port impedance
        port_vector = np.array([1.0] + [0.0] * (n_segments - 1))
        Z_1kHz = builder.compute_impedance(1e3, port_vector)
        Z_1MHz = builder.compute_impedance(1e6, port_vector)

        print(f"\nPort impedance:")
        print(f"  1 kHz: {Z_1kHz:.6e} Ohm")
        print(f"  1 MHz: {Z_1MHz:.6e} Ohm")

        # Method 2: Using convenience function
        print("\n--- Convenience function test ---")
        L2, R2, P2, M_LS2 = create_wire_peec(start, end, width, height, n_segments, sigma)
        print(f"create_wire_peec: L shape = {L2.shape}")

        print("\n*** C++ Wire Test PASSED ***")
        return True

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_cpp_loop():
    """Test C++ PEEC Loop Matrix Construction (pybind11)"""
    print("\n" + "=" * 60)
    print("Test: C++ PEEC Loop Matrix Construction (pybind11)")
    print("=" * 60)

    try:
        from peec_matrices import PEECBuilder, create_loop_peec

        # Loop parameters
        center = np.array([0, 0, 0])
        radius = 0.05  # 5cm
        normal = np.array([0, 0, 1])
        width = 2e-3  # 2mm
        height = 1e-3
        n_segments = 36
        sigma = 5.8e7

        # Method 1: Using PEECBuilder class
        builder = PEECBuilder()
        n_created = builder.create_loop(center, radius, normal,
                                        width, height, n_segments, sigma)

        print(f"\nLoop: R={radius*1000:.1f}mm, {n_segments} segments")
        print(f"Segments created: {n_created}")

        # Build matrices
        L, R, P, M_LS = builder.build(include_star=True)

        print(f"\nL matrix shape: {L.shape}")
        print(f"R vector shape: {R.shape}")
        print(f"P matrix shape: {P.shape if P is not None else 'None'}")

        # Inductance calculation
        L_total = np.sum(L)

        # Analytical: L = mu_0 * R * (ln(8R/a) - 2)
        a_equiv = np.sqrt(width * height / np.pi)
        L_analytical = 4 * np.pi * 1e-7 * radius * (np.log(8 * radius / a_equiv) - 2)

        print(f"\nL total (PEEC): {L_total*1e6:.3f} uH")
        print(f"L analytical: {L_analytical*1e6:.3f} uH")
        print(f"Ratio: {L_total/L_analytical:.3f}")

        # Port impedance
        port_vector = np.array([1.0] + [0.0] * (n_segments - 1))
        Z_100kHz = builder.compute_impedance(1e5, port_vector)

        print(f"\nPort impedance at 100 kHz: {Z_100kHz:.6e} Ohm")
        print(f"  |Z|: {abs(Z_100kHz):.3f} Ohm")
        print(f"  angle: {np.angle(Z_100kHz, deg=True):.1f} deg")

        # Method 2: Using convenience function
        print("\n--- Convenience function test ---")
        L2, R2, P2, M_LS2 = create_loop_peec(center, radius, normal, width, height, n_segments, sigma)
        print(f"create_loop_peec: L shape = {L2.shape}")

        print("\n*** C++ Loop Test PASSED ***")
        return True

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_compare_python():
    """Compare C++ and Python implementations"""
    print("\n" + "=" * 60)
    print("Test: Compare C++ (pybind11) vs Python Implementation")
    print("=" * 60)

    try:
        # Import C++ pybind11 module
        from peec_matrices import PEECBuilder as CppPEECBuilder

        # Import Python implementation
        from prima.peec_matrices import (
            create_simple_wire_segments,
            build_peec_matrices as py_build_peec_matrices
        )

        # Wire parameters
        length = 0.1
        width = 1e-3
        height = 1e-3
        n_segments = 10
        sigma = 5.8e7

        # C++ implementation
        cpp_builder = CppPEECBuilder()
        cpp_builder.create_wire(np.array([0, 0, 0]), np.array([length, 0, 0]),
                                width, height, n_segments, sigma)
        L_cpp, R_cpp, P_cpp, M_LS_cpp = cpp_builder.build(include_star=True)

        # Python implementation
        segments = create_simple_wire_segments(length, width, height, n_segments, sigma)
        matrices_py = py_build_peec_matrices(segments, include_star=True)
        L_py = matrices_py.L
        R_py = np.diag(matrices_py.R)  # Diagonal elements
        P_py = matrices_py.P
        M_LS_py = matrices_py.M_LS

        print(f"\nComparing L matrices:")
        print(f"  C++ L total: {np.sum(L_cpp)*1e9:.3f} nH")
        print(f"  Python L total: {np.sum(L_py)*1e9:.3f} nH")

        L_diff = np.abs(np.sum(L_cpp) - np.sum(L_py)) / np.abs(np.sum(L_py))
        print(f"  Relative difference: {L_diff:.2e}")

        print(f"\nComparing R matrices:")
        print(f"  C++ R total: {np.sum(R_cpp)*1000:.3f} mOhm")
        print(f"  Python R total: {np.sum(R_py)*1000:.3f} mOhm")

        R_diff = np.abs(np.sum(R_cpp) - np.sum(R_py)) / np.abs(np.sum(R_py))
        print(f"  Relative difference: {R_diff:.2e}")

        print(f"\nComparing P matrices:")
        print(f"  C++ P shape: {P_cpp.shape if P_cpp is not None else 'None'}")
        print(f"  Python P shape: {P_py.shape}")

        if L_diff < 0.01 and R_diff < 0.01:
            print("\n*** C++ vs Python Comparison PASSED ***")
            return True
        else:
            print("\n*** WARNING: Significant difference detected ***")
            return False

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_surface_impedance():
    """Test surface impedance (skin effect) functionality"""
    print("\n" + "=" * 60)
    print("Test: Surface Impedance (Skin Effect)")
    print("=" * 60)

    try:
        from peec_matrices import PEECBuilder

        # Wire parameters
        length = 0.1
        width = 1e-3
        height = 1e-3
        n_segments = 10
        sigma = 5.8e7
        mu_0 = 4 * np.pi * 1e-7

        builder = PEECBuilder()
        builder.create_wire(np.array([0, 0, 0]), np.array([length, 0, 0]),
                            width, height, n_segments, sigma)
        L, R, P, M_LS = builder.build(include_star=True)

        # Compute surface impedance at 1 MHz
        freq = 1e6
        omega = 2 * np.pi * freq
        delta = np.sqrt(2 / (omega * mu_0 * sigma))  # Skin depth

        # Per-unit-length surface impedance
        # Z_s = (1+j) / (sigma * delta) per unit length
        Z_s_per_length = (1 + 1j) / (sigma * delta)
        seg_length = length / n_segments
        Z_s_segment = Z_s_per_length * seg_length

        # Create surface impedance array
        Zs_diag = np.full(n_segments, Z_s_segment)

        print(f"\nFrequency: {freq/1e6:.1f} MHz")
        print(f"Skin depth: {delta*1e6:.1f} um")
        print(f"Z_s per segment: {Z_s_segment:.6e} Ohm")

        # Set surface impedance
        builder.set_surface_impedance(Zs_diag)

        # Compute impedance with and without surface impedance
        port_vector = np.array([1.0] + [0.0] * (n_segments - 1))
        Z_without_Zs = builder.compute_impedance(freq, port_vector)
        Z_with_Zs = builder.compute_impedance_with_Zs(freq, port_vector)

        print(f"\nImpedance at {freq/1e6:.1f} MHz:")
        print(f"  Without Z_s: {Z_without_Zs:.6e} Ohm")
        print(f"  With Z_s:    {Z_with_Zs:.6e} Ohm")
        print(f"  |Z| ratio:   {abs(Z_with_Zs)/abs(Z_without_Zs):.3f}")

        print("\n*** Surface Impedance Test PASSED ***")
        return True

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests"""
    print("C++ PEEC Matrix Construction Tests (pybind11)")
    print("=" * 60)

    results = []

    # Test 1: Wire
    results.append(("C++ Wire (pybind11)", test_cpp_wire()))

    # Test 2: Loop
    results.append(("C++ Loop (pybind11)", test_cpp_loop()))

    # Test 3: Python comparison
    results.append(("C++ vs Python", test_compare_python()))

    # Test 4: Surface impedance
    results.append(("Surface Impedance", test_surface_impedance()))

    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    for name, success in results:
        print(f"  {name}: {'PASS' if success else 'FAIL'}")

    if all(s for _, s in results):
        print("\n*** All tests PASSED ***")
    else:
        print("\n*** Some tests FAILED ***")


if __name__ == "__main__":
    main()
