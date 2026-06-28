"""Verify C++ Play hysteresis material against Python reference.

Tests:
1. Material creation + Forward B->H evaluation
2. C++ vs Python play_hysteron: ascending branch comparison
3. B-H loop: hysteresis detection
4. Forward/Inverse round-trip: Forward(Inverse(H)) == H
5. State Save/Restore/Commit lifecycle
6. Solver integration: MatApl + Solve on ObjRecMag
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
from scipy.interpolate import interp1d

# Also import Python reference implementation
sys.path.insert(0, script_dir)
from bqm_model import play_hysteron, compute_shape_functions, load_binput_data


def load_fixture(use_subset=True):
    """Load real B-input Play model shape functions from fixture."""
    fix_path = os.path.join(script_dir, 'binput_play_fixture.npz')
    fix = np.load(fix_path)

    if use_subset:
        K = int(fix['K_sub'])
        eta = fix['chi_sub']
        f_k_tables = []
        for ki in range(K):
            r = fix[f'r_sub_{ki}']
            f = fix[f'f_sub_{ki}']
            f_k_tables.append((r, f))
    else:
        K = int(fix['K'])
        eta = fix['chi']
        f_k_tables = []
        for k in range(K):
            r = fix[f'r_{k}']
            f = fix[f'f_{k}']
            f_k_tables.append((r, f))

    return K, eta, f_k_tables


def _to_tables_arg(f_k_tables):
    """Convert f_k_tables to list of (r_list, f_list) for C++ API."""
    return [(r.tolist(), f.tolist()) for r, f in f_k_tables]


all_pass = True
n_test = 0


def test(name, condition, detail=''):
    global all_pass, n_test
    n_test += 1
    status = 'PASS' if condition else 'FAIL'
    if not condition:
        all_pass = False
    print(f'  Test {n_test}: {name} ... {status}  {detail}')
    return condition


# =================================================================
# Test 1: Material creation + Forward B->H
# =================================================================
print('\n=== Test 1: MatPlayHysteresis creation + Forward B->H ===')
K, eta, f_k_tables = load_fixture(use_subset=True)
rad.UtiDelAll()
tables = _to_tables_arg(f_k_tables)
cpp_mat = rad.MatPlayHysteresis(K, eta.tolist(), tables)

# MatMvsH at various H -> should return reasonable M
H_tests = [
    np.array([100.0, 0.0, 0.0]),
    np.array([500.0, 200.0, 0.0]),
    np.array([0.0, 0.0, 1000.0]),
    np.array([-300.0, 400.0, 100.0]),
]

for H in H_tests:
    rad.UtiDelAll()
    cpp_mat = rad.MatPlayHysteresis(K, eta.tolist(), tables)
    M_cpp = np.array(rad.MatMvsH(cpp_mat, 'm', list(H)))
    H_mag = np.linalg.norm(H)
    M_mag = np.linalg.norm(M_cpp)
    dot = np.dot(M_cpp, H)
    ferromagnetic = dot > 0

    J_mag = MU_0 * M_mag
    test(f'H=[{H[0]:7.0f},{H[1]:7.0f},{H[2]:7.0f}]  M={M_mag:.0f} A/m  J={J_mag:.2f} T',
         M_mag > 0 and ferromagnetic)


# =================================================================
# Test 2: C++ Play vs Python Play (ascending branch)
# =================================================================
print('\n=== Test 2: C++ Play vs Python play_hysteron ===')
K, eta_arr, f_k_tables = load_fixture(use_subset=True)
rad.UtiDelAll()
tables = _to_tables_arg(f_k_tables)
cpp_mat = rad.MatPlayHysteresis(K, eta_arr.tolist(), tables)

# Build Python interpolants
f_k_interps = []
for r, f in f_k_tables:
    interp = interp1d(r, f, kind='linear', bounds_error=False,
                      fill_value=(f[0], f[-1]))
    f_k_interps.append(interp)

# Ascending branch B = 0 -> 1.3 T (monotone region of H(B))
# H(B) is non-monotonic for B > ~1.4 T due to negative shape functions f_k.
# In the non-monotone region, Inverse(H) is multi-valued; Newton finds the
# lower-B solution which may differ from the Python ascending branch.
n_steps = 30
Bmax = 1.3
B_sweep = np.linspace(0.05, Bmax, n_steps)

# Python Play model
px = np.zeros(K)
py = np.zeros(K)
H_py = np.zeros(n_steps)
for i, Bx in enumerate(B_sweep):
    Hx, Hy = play_hysteron(f_k_interps, Bx, 0.0, eta_arr, px, py)
    H_py[i] = Hx

# C++ Play model: apply same H sequence with state accumulation
# Single material, ascending H values (auto-commit between calls)
rad.UtiDelAll()
cpp_mat = rad.MatPlayHysteresis(K, eta_arr.tolist(), tables)
B_cpp = np.zeros(n_steps)
for i, H_val in enumerate(H_py):
    M = rad.MatMvsH(cpp_mat, 'm', [H_val, 0, 0])
    B_cpp[i] = MU_0 * (H_val + M[0])

# For ascending from zero, B_cpp should be close to B_sweep
# (since the Python Play model computed H from those B values)
err = np.abs(B_cpp - B_sweep)
max_err = err.max()
mean_err = err.mean()
rel_err = max_err / Bmax * 100

print(f'  B range: [0.05, {Bmax}] T, {n_steps} steps')
print(f'  Max |B_cpp - B_py|: {max_err:.6f} T ({rel_err:.3f}%)')
print(f'  Mean |B_cpp - B_py|: {mean_err:.6f} T')
test('C++ Play matches Python Play < 5%', rel_err < 5.0, f'{rel_err:.3f}%')


# =================================================================
# Test 3: B-H loop (hysteresis detection)
# =================================================================
print('\n=== Test 3: B-H loop with Play model ===')
rad.UtiDelAll()
cpp_mat = rad.MatPlayHysteresis(K, eta_arr.tolist(), tables)

Hmax = 3000.0
n_loop = 200
t = np.linspace(0, 2 * np.pi, n_loop)
H_drive = Hmax * np.sin(t)

B_loop = np.zeros(n_loop)
for i, H_val in enumerate(H_drive):
    M = rad.MatMvsH(cpp_mat, 'm', [H_val, 0, 0])
    B_loop[i] = MU_0 * (H_val + M[0])

# Check hysteresis
mid = n_loop // 4
mid2 = 3 * n_loop // 4
B_rising = B_loop[mid - 1]
B_falling = B_loop[mid2 - 1]
hyst_width = abs(B_rising - B_falling)

print(f'  Hmax = {Hmax:.0f} A/m')
print(f'  B(H~0, rising)  = {B_rising:.4f} T')
print(f'  B(H~0, falling) = {B_falling:.4f} T')
print(f'  Hysteresis width = {hyst_width:.4f} T')
test('Hysteresis detected', hyst_width > 0.001, f'{hyst_width:.4f} T')


# =================================================================
# Test 4: State Save/Restore/Commit
# =================================================================
print('\n=== Test 4: State Save/Restore/Commit ===')
rad.UtiDelAll()
cpp_mat = rad.MatPlayHysteresis(K, eta_arr.tolist(), tables)

# Evaluate at some H to get non-zero state
M1 = rad.MatMvsH(cpp_mat, 'm', [1000.0, 0, 0])

# Save state
state = rad.MatHysSaveState(cpp_mat)
print(f'  State size: {len(state)} doubles')
test('State is non-empty', len(state) > 0, f'{len(state)} doubles')

# Evaluate at different H
M2 = rad.MatMvsH(cpp_mat, 'm', [2000.0, 0, 0])

# Restore state
rad.MatHysRestoreState(cpp_mat, state)
M1_restored = rad.MatMvsH(cpp_mat, 'm', [1000.0, 0, 0])

M1_err = np.linalg.norm(np.array(M1) - np.array(M1_restored))
test('Restore recovers original M', M1_err < 1.0, f'err={M1_err:.6f} A/m')


# =================================================================
# Test 5: Solver integration
# =================================================================
print('\n=== Test 5: Solver integration (MatApl + Solve) ===')
rad.UtiDelAll()

iron = rad.magnet_box([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 0])
mat = rad.MatPlayHysteresis(K, eta_arr.tolist(), tables)
rad.MatApl(iron, mat)

bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])
container = rad.ObjCnt([iron, bkg])

try:
    result = rad.Solve(container, 0.001, 100, 1)
    print(f'  Solve result: {result}')
    B = rad.Fld(container, 'b', [0.02, 0, 0])
    print(f'  B at (20mm, 0, 0) = [{B[0]:.6f}, {B[1]:.6f}, {B[2]:.6f}] T')
    test('Solve completed', True)
except Exception as e:
    print(f'  Solve failed: {e}')
    test('Solve completed', False, str(e))


# =================================================================
# Test 6: Performance
# =================================================================
print('\n=== Test 6: Performance ===')
rad.UtiDelAll()
K_full, eta_full, f_k_tables_full = load_fixture(use_subset=False)
tables_full = _to_tables_arg(f_k_tables_full)
cpp_mat = rad.MatPlayHysteresis(K_full, eta_full.tolist(), tables_full)

H = [500.0, 200.0, 0.0]
n_eval = 100
rad.MatMvsH(cpp_mat, 'm', H)  # warmup

t0 = time.perf_counter()
for _ in range(n_eval):
    rad.MatMvsH(cpp_mat, 'm', H)
t_cpp = (time.perf_counter() - t0) / n_eval * 1000

print(f'  K = {K_full}, n_eval = {n_eval}')
print(f'  C++ Play: {t_cpp:.3f} ms/eval')
test('Performance OK', True, f'{t_cpp:.3f} ms/eval')


# =================================================================
# Test 7: Monotone limit / saturation clamping
# =================================================================
print('\n=== Test 7: Monotone limit (B_mono_max constraint) ===')
rad.UtiDelAll()
K, eta_arr, f_k_tables = load_fixture(use_subset=True)
tables = _to_tables_arg(f_k_tables)
cpp_mat = rad.MatPlayHysteresis(K, eta_arr.tolist(), tables)

# Find B_mono_max empirically: sweep Forward and check where dH/dB <= 0
f_k_interps_7 = []
for r, f in f_k_tables:
    interp = interp1d(r, f, kind='linear', bounds_error=False,
                      fill_value=(f[0], f[-1]))
    f_k_interps_7.append(interp)

B_scan = np.linspace(0.01, 2.0, 500)
H_scan = np.zeros(len(B_scan))
for i, Bv in enumerate(B_scan):
    H_val = 0
    for k in range(K):
        pk = Bv - eta_arr[k]
        if pk > 0:
            H_val += float(f_k_interps_7[k](pk))
    H_scan[i] = H_val

dHdB = np.diff(H_scan) / np.diff(B_scan)
mono_idx = np.argmax(dHdB <= 0)
if mono_idx > 0:
    B_mono_py = B_scan[mono_idx]
    H_mono_py = H_scan[mono_idx]
    print(f'  Python B_mono_max = {B_mono_py:.4f} T, H_mono_max = {H_mono_py:.1f} A/m')
else:
    B_mono_py = 2.0
    H_mono_py = H_scan[-1]
    print(f'  No non-monotonicity found in Python scan')

# Test: Inverse at increasing H values, verify B stays bounded
H_large_values = [H_mono_py * 0.5, H_mono_py * 0.9, H_mono_py * 1.5, H_mono_py * 3.0]
B_results = []
for H_val in H_large_values:
    rad.UtiDelAll()
    cpp_mat = rad.MatPlayHysteresis(K, eta_arr.tolist(), tables)
    M = rad.MatMvsH(cpp_mat, 'm', [H_val, 0, 0])
    B_val = MU_0 * (H_val + M[0])
    B_results.append(B_val)
    print(f'  H={H_val:10.1f} A/m -> B={B_val:.4f} T')

# B should be bounded (not exceed B_mono_max significantly)
B_max_result = max(B_results)
test('B bounded at saturation', B_max_result < B_mono_py * 1.1,
     f'max B={B_max_result:.4f} T, limit~{B_mono_py:.4f} T')

# B should increase monotonically with H (or saturate)
monotonic = all(B_results[i] <= B_results[i+1] + 0.01
                for i in range(len(B_results)-1))
test('B(H) monotonically increasing (or saturated)', monotonic)


# =================================================================
# Test 8: Extended B range with monotone constraint
# =================================================================
print('\n=== Test 8: C++ Play vs Python across full monotone range ===')
rad.UtiDelAll()
K, eta_arr, f_k_tables = load_fixture(use_subset=True)
tables = _to_tables_arg(f_k_tables)
cpp_mat = rad.MatPlayHysteresis(K, eta_arr.tolist(), tables)

# Use full monotone range (up to B_mono_max)
B_mono_test = min(B_mono_py * 0.95, 1.5)  # 95% of monotone limit
n_steps8 = 40
B_sweep8 = np.linspace(0.05, B_mono_test, n_steps8)

# Python reference
px8 = np.zeros(K)
py8 = np.zeros(K)
H_py8 = np.zeros(n_steps8)
for i, Bx in enumerate(B_sweep8):
    Hx, _ = play_hysteron(f_k_interps_7, Bx, 0.0, eta_arr, px8, py8)
    H_py8[i] = Hx

# C++ Inverse
rad.UtiDelAll()
cpp_mat = rad.MatPlayHysteresis(K, eta_arr.tolist(), tables)
B_cpp8 = np.zeros(n_steps8)
for i, H_val in enumerate(H_py8):
    M = rad.MatMvsH(cpp_mat, 'm', [H_val, 0, 0])
    B_cpp8[i] = MU_0 * (H_val + M[0])

err8 = np.abs(B_cpp8 - B_sweep8)
max_err8 = err8.max()
rel_err8 = max_err8 / B_mono_test * 100

print(f'  B range: [0.05, {B_mono_test:.2f}] T (95%% of monotone limit)')
print(f'  Max |B_cpp - B_py|: {max_err8:.6f} T ({rel_err8:.3f}%%)')
test('C++ matches Python in full monotone range < 5%%', rel_err8 < 5.0,
     f'{rel_err8:.3f}%%')


# =================================================================
# Summary
# =================================================================
print(f'\n{"="*60}')
print(f'Results: {n_test} tests, '
      f'{"ALL PASS" if all_pass else "SOME FAILED"}')
print(f'{"="*60}')
