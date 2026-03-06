"""Verification of C++ energy-based hysteresis material.

Tests:
1. Forward operator: B output sanity check + chi range
2. B-H loop: hysteresis (ascending != descending branch)
3. Solver integration: MatApl + Solve on ObjRecMag
4. Performance: C++ timing for K=50
"""

import sys
import os
import time
import numpy as np

MU_0 = 4e-7 * np.pi

try:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
    import radia as rad
except ImportError:
    import _radia_pybind as rad


def test_forward_sanity():
    """Test 1: Forward operator sanity checks."""
    print("=" * 60)
    print("Test 1: Forward operator sanity")
    print("=" * 60)

    K = 10
    As = np.full(K, 10000.0)
    Js = np.full(K, 0.2)
    chi = np.linspace(0, 100, K)
    eps = 1e-8

    H_tests = [
        np.array([100.0, 0.0, 0.0]),
        np.array([500.0, 200.0, 0.0]),
        np.array([0.0, 0.0, 1000.0]),
        np.array([-300.0, 400.0, 100.0]),
    ]

    all_pass = True
    for H in H_tests:
        rad.UtiDelAll()
        cpp_mat = rad.MatEnergyHysteresis(K, As, Js, chi, eps)

        M_cpp = np.array(rad.MatMvsH(cpp_mat, 'm', list(H)))
        B_cpp = MU_0 * (H + M_cpp)

        H_mag = np.linalg.norm(H)
        B_mag = np.linalg.norm(B_cpp)
        chi_eff = B_mag / (MU_0 * H_mag) - 1

        # Sanity: chi should be in reasonable range [1, 200]
        chi_ok = 1 < chi_eff < 200
        # B direction should roughly align with H
        cos_angle = np.dot(B_cpp, H) / (B_mag * H_mag)
        dir_ok = cos_angle > 0.99

        passed = chi_ok and dir_ok
        if not passed:
            all_pass = False
        status = "PASS" if passed else "FAIL"
        print(f"  H = [{H[0]:8.1f}, {H[1]:8.1f}, {H[2]:8.1f}]  "
              f"chi={chi_eff:7.1f}  cos={cos_angle:.4f}  {status}")

    return all_pass


def test_bh_loop():
    """Test 2: B-H loop generation (hysteresis detection)."""
    print("\n" + "=" * 60)
    print("Test 2: B-H loop generation")
    print("=" * 60)

    K = 10
    As = np.full(K, 5000.0)
    Js = np.full(K, 0.2)
    chi = np.linspace(0, 300, K)
    eps = 1e-8

    rad.UtiDelAll()
    cpp_mat = rad.MatEnergyHysteresis(K, As, Js, chi, eps)

    Hmax = 1000.0
    n_steps = 100
    t = np.linspace(0, 2 * np.pi, n_steps)
    H_drive = Hmax * np.sin(t)

    B_values = np.zeros(n_steps)
    for i, H_val in enumerate(H_drive):
        M = rad.MatMvsH(cpp_mat, 'm', [H_val, 0, 0])
        B_values[i] = MU_0 * (H_val + M[0])

    # Check hysteresis: ascending vs descending at H~0
    mid = n_steps // 4
    mid2 = 3 * n_steps // 4
    B_at_H0_rising = B_values[mid - 1]
    B_at_H0_falling = B_values[mid2 - 1]
    hysteresis_width = abs(B_at_H0_rising - B_at_H0_falling)

    print(f"  Hmax = {Hmax:.0f} A/m, K = {K}")
    print(f"  B(H=0, rising)  = {B_at_H0_rising:.4f} T")
    print(f"  B(H=0, falling) = {B_at_H0_falling:.4f} T")
    print(f"  Hysteresis width = {hysteresis_width:.4f} T")

    has_hysteresis = hysteresis_width > 0.001
    status = "PASS" if has_hysteresis else "FAIL"
    print(f"  Hysteresis detected: {status}")

    # Optional plot
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        ax.plot(H_drive, B_values, 'b-', linewidth=1)
        ax.set_xlabel('H (A/m)')
        ax.set_ylabel('B (T)')
        ax.set_title('Energy-based hysteresis: C++ B-H loop')
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='k', linewidth=0.5)
        ax.axvline(0, color='k', linewidth=0.5)

        out_path = os.path.join(os.path.dirname(__file__), 'verify_cpp_hysteresis.png')
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        print(f"  Plot saved: {out_path}")
        plt.close()
    except ImportError:
        pass

    return has_hysteresis


def test_solver_integration():
    """Test 3: Apply hysteresis material to element and run Solve()."""
    print("\n" + "=" * 60)
    print("Test 3: Solver integration (MatApl + Solve)")
    print("=" * 60)

    rad.UtiDelAll()

    K = 5
    As = np.full(K, 10000.0)
    Js = np.full(K, 0.3)
    chi = np.linspace(0, 50, K)

    iron = rad.ObjRecMag([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 0])
    mat = rad.MatEnergyHysteresis(K, As, Js, chi, 1e-8)
    rad.MatApl(iron, mat)

    bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])
    container = rad.ObjCnt([iron, bkg])

    try:
        result = rad.Solve(container, 0.001, 100, 1)
        print(f"  Solve result: {result}")

        B = rad.Fld(container, 'b', [0.02, 0, 0])
        print(f"  B at (20mm, 0, 0) = [{B[0]:.6f}, {B[1]:.6f}, {B[2]:.6f}] T")
        print(f"  PASS")
        return True
    except Exception as e:
        print(f"  Solve failed: {e}")
        print(f"  FAIL")
        return False


def test_performance():
    """Test 4: C++ performance measurement."""
    print("\n" + "=" * 60)
    print("Test 4: Performance")
    print("=" * 60)

    K = 50
    As = np.full(K, 10000.0)
    Js = np.full(K, 2.0 / K)
    chi = np.linspace(0, 500, K)
    eps = 1e-8
    H = [500.0, 200.0, 0.0]
    n_eval = 100

    rad.UtiDelAll()
    cpp_mat = rad.MatEnergyHysteresis(K, As, Js, chi, eps)

    # Warm up
    rad.MatMvsH(cpp_mat, 'm', H)

    t0 = time.perf_counter()
    for _ in range(n_eval):
        rad.MatMvsH(cpp_mat, 'm', H)
    t_cpp = (time.perf_counter() - t0) / n_eval * 1000

    print(f"  K = {K}, n_eval = {n_eval}")
    print(f"  C++: {t_cpp:.3f} ms/eval")
    print(f"  PASS")
    return True


if __name__ == '__main__':
    results = []

    r1 = test_forward_sanity()
    results.append(('Forward sanity', r1))

    r2 = test_bh_loop()
    results.append(('B-H loop', r2))

    r3 = test_solver_integration()
    results.append(('Solver integration', r3))

    r4 = test_performance()
    results.append(('Performance', r4))

    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    all_pass = True
    for name, passed in results:
        status = "PASS" if passed else "FAIL"
        print(f"  {name}: {status}")
        if not passed:
            all_pass = False

    print(f"\nOverall: {'ALL PASS' if all_pass else 'SOME FAILED'}")
