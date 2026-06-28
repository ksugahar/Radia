"""Verification of C++ energy-based hysteresis material with real B-input data.

Uses shape functions built from measured B-H loops (Jiles-Atherton model,
40 major loops) via the B-input Play model pipeline:
    BH.mat -> build_shape_functions() -> f_k_tables -> MatEnergyHysteresis()

Direction convention (B-input Play model):
    Forward  = B -> H  (natural direction for B-input)
    Inverse  = H -> B  (via H -> M, then B = mu_0*(H+M))

Tests:
1. Material creation + H->M evaluation with real shape functions
2. B-H loop: hysteresis detection (ascending != descending)
3. C++ energy model vs Python Play model comparison (diagnostic)
4. Solver integration: MatApl + Solve on ObjRecMag
5. Performance: C++ timing for K=80 real operators

NOTE on Test 3 (C++ vs Python Play):
    The B-input Play model computes H = sum_k f_k(|p_k|) (shape functions
    contribute to H). The C++ energy model computes M = sum_k J_k (operators
    contribute to M). Using B-input Play shape functions directly as energy
    model parameters produces unrealistically large M because:
    - Each J_k can reach r_max = 2.0 T
    - With K=20 operators, total J can reach 40 T
    - Real materials saturate at J_s ~ 2 T
    A proper conversion from B-input Play to energy model is needed.
"""

import sys
import os
import time
import numpy as np

MU_0 = 4e-7 * np.pi

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(script_dir))
sys.path.insert(0, os.path.join(repo_root, 'src', 'radia'))

import radia as rad
from hysteresis_io import play_hysteron
from scipy.interpolate import interp1d


def load_fixture(use_subset=True):
    """Load real B-input Play model shape functions from fixture.

    The fixture was generated from 40 Jiles-Atherton B-H loops via
    build_shape_functions() in hysteresis_io.py.

    Parameters
    ----------
    use_subset : bool
        If True, use K=20 subset for faster tests.
        If False, use full K=80 operators.

    Returns
    -------
    K : int
    chi : ndarray
    f_k_tables : list of (r, f) tuples
    """
    fix_path = os.path.join(script_dir, 'binput_play_fixture.npz')
    fix = np.load(fix_path)

    if use_subset:
        K = int(fix['K_sub'])
        chi = fix['chi_sub']
        f_k_tables = []
        for ki in range(K):
            r = fix[f'r_sub_{ki}']
            f = fix[f'f_sub_{ki}']
            f_k_tables.append((r, f))
    else:
        K = int(fix['K'])
        chi = fix['chi']
        f_k_tables = []
        for k in range(K):
            r = fix[f'r_{k}']
            f = fix[f'f_{k}']
            f_k_tables.append((r, f))

    return K, chi, f_k_tables


def _to_tables_arg(f_k_tables):
    """Convert f_k_tables to list of (r_list, f_list) for C++ API."""
    return [(r.tolist(), f.tolist()) for r, f in f_k_tables]


def test_material_creation():
    """Test 1: Material creation and H->M evaluation (Inverse direction)."""
    print("=" * 60)
    print("Test 1: Material creation + Inverse H->M (real data)")
    print("=" * 60)

    K, chi, f_k_tables = load_fixture(use_subset=True)
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
        tables = _to_tables_arg(f_k_tables)
        cpp_mat = rad.MatEnergyHysteresis(K, chi.tolist(), tables, eps)

        # MatMvsH: H -> M (Inverse direction in B-input convention)
        M_cpp = np.array(rad.MatMvsH(cpp_mat, 'm', list(H)))
        B_cpp = MU_0 * (H + M_cpp)

        H_mag = np.linalg.norm(H)
        M_mag = np.linalg.norm(M_cpp)

        # Minimal sanity: M should be nonzero and ferromagnetic (M parallel to H)
        dot = np.dot(M_cpp, H)
        ferromagnetic = dot > 0

        passed = M_mag > 0 and ferromagnetic
        if not passed:
            all_pass = False
        status = "PASS" if passed else "FAIL"
        J_mag = MU_0 * M_mag  # magnetic polarization in Tesla
        print(f"  H=[{H[0]:7.0f},{H[1]:7.0f},{H[2]:7.0f}]  "
              f"M={M_mag:.0f} A/m  J={J_mag:.2f} T  {status}")

    return all_pass


def test_bh_loop():
    """Test 2: B-H loop generation with real data (hysteresis detection)."""
    print("\n" + "=" * 60)
    print("Test 2: B-H loop (hysteresis detection)")
    print("=" * 60)

    K, chi, f_k_tables = load_fixture(use_subset=True)
    eps = 1e-8

    rad.UtiDelAll()
    tables = _to_tables_arg(f_k_tables)
    cpp_mat = rad.MatEnergyHysteresis(K, chi.tolist(), tables, eps)

    Hmax = 3000.0
    n_steps = 200
    t = np.linspace(0, 2 * np.pi, n_steps)
    H_drive = Hmax * np.sin(t)

    B_values = np.zeros(n_steps)
    for i, H_val in enumerate(H_drive):
        # MatMvsH: H -> M (Inverse direction)
        M = rad.MatMvsH(cpp_mat, 'm', [H_val, 0, 0])
        B_values[i] = MU_0 * (H_val + M[0])

    # Check hysteresis: ascending vs descending at H~0
    mid = n_steps // 4
    mid2 = 3 * n_steps // 4
    B_at_H0_rising = B_values[mid - 1]
    B_at_H0_falling = B_values[mid2 - 1]
    hysteresis_width = abs(B_at_H0_rising - B_at_H0_falling)

    print(f"  Hmax = {Hmax:.0f} A/m, K = {K}")
    print(f"  B(H~0, rising)  = {B_at_H0_rising:.4f} T")
    print(f"  B(H~0, falling) = {B_at_H0_falling:.4f} T")
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
        ax.set_title(f'C++ Energy model B-H loop (K={K}, real data)')
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='k', linewidth=0.5)
        ax.axvline(0, color='k', linewidth=0.5)

        out_path = os.path.join(script_dir, 'verify_cpp_hysteresis.png')
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        print(f"  Plot saved: {out_path}")
        plt.close()
    except ImportError:
        pass

    return has_hysteresis


def test_cpp_vs_python_play():
    """Test 3: Compare C++ energy model vs Python Play model (diagnostic).

    The B-input Play model: B -> p_k -> H = sum f_k(|p_k|)
    The C++ energy model: H -> J_k -> M = sum J_k

    These are different formulations. This test quantifies the discrepancy.
    """
    print("\n" + "=" * 60)
    print("Test 3: C++ energy vs Python Play (diagnostic)")
    print("=" * 60)

    K, chi, f_k_tables = load_fixture(use_subset=True)
    eta = chi.copy()  # chi_k = eta_k for B-input Play
    eps = 1e-8

    # --- Python Play model: Forward B -> H ---
    f_k_interps = []
    for r, f in f_k_tables:
        interp = interp1d(r, f, kind='linear', bounds_error=False,
                          fill_value=(f[0], f[-1]))
        f_k_interps.append(interp)

    # Ascending branch
    n_steps = 20
    Bmax = 1.5
    B_sweep = np.linspace(0.05, Bmax, n_steps)

    px = np.zeros(K)
    py = np.zeros(K)
    H_py = np.zeros(n_steps)
    for i, Bx in enumerate(B_sweep):
        Hx, Hy, px, py = play_hysteron(f_k_interps, Bx, 0.0, eta, px, py)
        H_py[i] = Hx

    # --- C++ energy model: Inverse H -> M -> B ---
    rad.UtiDelAll()
    tables = _to_tables_arg(f_k_tables)
    cpp_mat = rad.MatEnergyHysteresis(K, chi.tolist(), tables, eps)

    B_cpp = np.zeros(n_steps)
    for i, H_val in enumerate(H_py):
        M = rad.MatMvsH(cpp_mat, 'm', [H_val, 0, 0])
        B_cpp[i] = MU_0 * (H_val + M[0])

    # Report
    print(f"  {'B_play':>8s}  {'H(A/m)':>10s}  {'B_cpp(T)':>10s}  {'ratio':>8s}")
    print("  " + "-" * 45)
    for i in range(0, n_steps, 2):
        ratio = B_cpp[i] / B_sweep[i] if B_sweep[i] > 0 else 0
        print(f"  {B_sweep[i]:8.3f}  {H_py[i]:10.1f}  {B_cpp[i]:10.2f}  {ratio:8.1f}x")

    avg_ratio = np.mean(B_cpp / B_sweep)
    print(f"\n  Average B_cpp/B_play ratio: {avg_ratio:.1f}x")
    print(f"  (Expected ~{K:.0f}x due to K operators each contributing J_k)")
    print(f"  INFO: B-input Play -> energy conversion needs proper scaling")

    # Informational test, always passes
    return True


def test_solver_integration():
    """Test 4: Apply hysteresis material to element and run Solve()."""
    print("\n" + "=" * 60)
    print("Test 4: Solver integration (MatApl + Solve)")
    print("=" * 60)

    rad.UtiDelAll()

    K, chi, f_k_tables = load_fixture(use_subset=True)
    tables = _to_tables_arg(f_k_tables)

    iron = rad.magnet_box([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 0])
    mat = rad.MatEnergyHysteresis(K, chi.tolist(), tables, 1e-8)
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
    """Test 5: C++ performance measurement with full K=80."""
    print("\n" + "=" * 60)
    print("Test 5: Performance (K=80 real operators)")
    print("=" * 60)

    K, chi, f_k_tables = load_fixture(use_subset=False)
    eps = 1e-8
    H = [500.0, 200.0, 0.0]
    n_eval = 100

    rad.UtiDelAll()
    tables = _to_tables_arg(f_k_tables)
    cpp_mat = rad.MatEnergyHysteresis(K, chi.tolist(), tables, eps)

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

    r1 = test_material_creation()
    results.append(('Material creation + Inverse H->M', r1))

    r2 = test_bh_loop()
    results.append(('B-H loop (hysteresis)', r2))

    r3 = test_cpp_vs_python_play()
    results.append(('C++ vs Python Play (diagnostic)', r3))

    r4 = test_solver_integration()
    results.append(('Solver integration', r4))

    r5 = test_performance()
    results.append(('Performance K=80', r5))

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
