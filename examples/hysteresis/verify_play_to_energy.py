"""Verify B-input Play model from Potter-Schmulian data (Step 3).

Tests:
1. Shape function computation (JMAG identification method)
2. Shape function properties (anti-symmetry, positivity)
3. Play model sinusoidal trajectory
4. Play model matches MATLAB B_input_case_1 reference
5. Descending branch reproduction
6. Multiple amplitude loops
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from bqm_model import (
    compute_shape_functions, play_hysteron, load_binput_data,
    run_play_model,
)

DATA_DIR = (r'W:\999_菅原賢悟\19_磁気ヒステリシス'
            r'\2024_IGTE_共同研究\2024_03_08_H-input_B-input')
MAT_FILE = os.path.join(DATA_DIR, 'Potter_Schmulian', 'B_input.mat')
REF_FILE = os.path.join(DATA_DIR, 'B_input_case_1_B-input.mat')

all_pass = True
n_test = 0


def test(name, condition, detail=''):
    global all_pass, n_test
    n_test += 1
    status = 'PASS' if condition else 'FAIL'
    if not condition:
        all_pass = False
    print(f'  Test {n_test}: {name} ... {status}  {detail}')


# =================================================================
# Test 1: Load data and compute shape functions
# =================================================================
print('\n=== Test 1: Shape function computation ===')
data = load_binput_data(MAT_FILE)
print(f'  K={data["K"]} BMax levels, dB={data["dB"]:.4f} T')
print(f'  BMax range: [{data["BMax"][0]:.3f}, {data["BMax"][-1]:.3f}] T')
print(f'  HdataLUT shape: {data["HdataLUT"].shape}')
print(f'  Bplay: {len(data["Bplay"])} points')
print(f'  eta: {len(data["eta"])} operators')

f_interps = compute_shape_functions(data['HdataLUT'], data['Bplay'])
Np = len(f_interps)
eta = data['eta']
Bmax = data['Bmax']
print(f'  Number of play operators Np: {Np}')
test('Np = 2*K (operators doubled by half-step grid)',
     Np == 2 * data['K'], f'{Np} vs {2*data["K"]}')


# =================================================================
# Test 2: Anti-symmetry of shape functions
# =================================================================
print('\n=== Test 2: Shape function properties ===')

r_test = 0.5
max_asym = 0.0
for k in range(Np):
    fp = float(f_interps[k](r_test))
    fm = float(f_interps[k](-r_test))
    if abs(fp) > 1e-10:
        asym = abs(fp + fm) / abs(fp) * 100
    else:
        asym = abs(fp + fm)
    max_asym = max(max_asym, asym)
test('Anti-symmetry f(-r) = -f(r)', max_asym < 0.01,
     f'max asym={max_asym:.6f}%')

f0_at_tip = float(f_interps[0](Bmax))
test('f_0 at Bmax is positive', f0_at_tip > 0, f'{f0_at_tip:.2f} A/m')

# Count operators with negative regions
n_negative = 0
for k in range(Np):
    r_check = np.linspace(0.01, max(0.02, Bmax - eta[k] - 0.01), 50)
    f_vals = np.array([float(f_interps[k](r)) for r in r_check])
    if np.any(f_vals < -1e-6):
        n_negative += 1
print(f'  Operators with negative regions: {n_negative}/{Np}')


# =================================================================
# Test 3: Play model sinusoidal trajectory
# =================================================================
print('\n=== Test 3: Play model sinusoidal ===')

n_periods = 5
t = np.linspace(0, n_periods * np.pi, 1000)
Bm = 1.05
Bx_traj = Bm * np.sin(2 * t)
By_traj = np.zeros_like(t)

Hx, Hy = run_play_model(f_interps, eta, Bx_traj, By_traj)

H_range = Hx.max() - Hx.min()
print(f'  B amplitude: {Bm:.2f} T, H range: {H_range:.1f} A/m')
test('H range is physical (100-3000 A/m)', 100 < H_range < 3000,
     f'{H_range:.1f} A/m')

# Check hysteresis via coercivity
idx_up = idx_down = None
for i in range(500, len(Bx_traj) - 1):
    if Bx_traj[i] <= 0 < Bx_traj[i + 1] and idx_up is None:
        idx_up = i
    if Bx_traj[i] >= 0 > Bx_traj[i + 1] and idx_down is None:
        idx_down = i
    if idx_up is not None and idx_down is not None:
        break

if idx_up is not None and idx_down is not None:
    coercivity = abs(Hx[idx_up] - Hx[idx_down]) / 2
    print(f'  Coercivity: {coercivity:.1f} A/m')
    test('Coercivity > 0', coercivity > 10, f'{coercivity:.1f} A/m')


# =================================================================
# Test 4: Compare with MATLAB B_input_case_1 reference
# =================================================================
print('\n=== Test 4: MATLAB B_input_case_1 reference ===')

Hx_py = None
Hx_ref = None
if os.path.exists(REF_FILE):
    import scipy.io as sio
    ref = sio.loadmat(REF_FILE)
    Bx_ref = ref['Bx'].flatten()
    By_ref = ref['By'].flatten()
    Hx_ref = ref['Hx'].flatten()

    Hx_py, _ = run_play_model(f_interps, eta, Bx_ref, By_ref)

    H_err = np.abs(Hx_py - Hx_ref)
    H_ref_range = Hx_ref.max() - Hx_ref.min()
    max_rel_err = H_err.max() / H_ref_range * 100
    mean_rel_err = H_err.mean() / H_ref_range * 100

    print(f'  Max error: {H_err.max():.4f} A/m ({max_rel_err:.4f}%)')
    print(f'  Mean error: {H_err.mean():.4f} A/m ({mean_rel_err:.4f}%)')
    test('Match MATLAB < 0.5%', max_rel_err < 0.5,
         f'{max_rel_err:.4f}%')
else:
    print(f'  Reference not found: {REF_FILE}')
    test('MATLAB reference', False, 'File not found')


# =================================================================
# Test 5: Descending branch reproduction
# =================================================================
print('\n=== Test 5: Descending branch reproduction ===')

BH = data['BH']
H_data, B_data = BH[0]  # largest BMax curve

# Initialize: ramp up to Bmax
n_init = 200
B_init = np.linspace(0, Bmax, n_init)
px0 = np.zeros(Np)
py0 = np.zeros(Np)
for i in range(n_init):
    play_hysteron(f_interps, B_init[i], 0.0, eta, px0, py0)

# Descend along data curve
n_desc = len(B_data)
Hx_desc = np.zeros(n_desc)
Bx_desc = B_data.copy()
for i in range(n_desc):
    hx, _ = play_hysteron(f_interps, Bx_desc[i], 0.0, eta, px0, py0)
    Hx_desc[i] = hx

desc_err = np.abs(Hx_desc - H_data)
desc_range = abs(H_data.max() - H_data.min())
desc_rel = desc_err.max() / desc_range * 100 if desc_range > 0 else 0
print(f'  H_play tip: {Hx_desc[0]:.2f} vs H_data: {H_data[0]:.2f} A/m')
print(f'  Max desc error: {desc_err.max():.4f} A/m ({desc_rel:.4f}%)')
test('Descending branch < 1%', desc_rel < 1.0, f'{desc_rel:.4f}%')


# =================================================================
# Test 6: Multiple amplitudes (monotone H range)
# =================================================================
print('\n=== Test 6: Multiple amplitude loops ===')

amplitudes = [0.2, 0.5, 1.0, 1.5, 1.8]
loop_data = {}

for Bm_loop in amplitudes:
    t_loop = np.linspace(0, 6 * np.pi, 500)
    Bx_loop = Bm_loop * np.sin(2 * t_loop)
    By_loop = np.zeros_like(t_loop)
    Hx_loop, _ = run_play_model(f_interps, eta, Bx_loop, By_loop)
    H_range_loop = Hx_loop[200:].max() - Hx_loop[200:].min()
    loop_data[Bm_loop] = (Bx_loop[200:], Hx_loop[200:], H_range_loop)

ranges = [loop_data[Bm_loop][2] for Bm_loop in amplitudes]
monotone = all(ranges[i] < ranges[i + 1] for i in range(len(ranges) - 1))
for Bm_loop in amplitudes:
    print(f'  Bm={Bm_loop:.1f}T: H_range={loop_data[Bm_loop][2]:.1f} A/m')
test('H range increases with B amplitude', monotone)


# =================================================================
# Structural analysis
# =================================================================
print('\n=== Play vs Energy model structural analysis ===')
r_max_sum = sum(Bmax - eta[k] for k in range(Np))
print(f'  Sum of operator r_max = {r_max_sum:.2f} T')
print(f'  Physical J_s ~ 2 T (ratio {r_max_sum/2:.0f}x)')
print(f'  Play: H = sum g_k(|play_k(B)|) [independent operators on B]')
print(f'  Energy: B = mu_0*H + sum J_k [coupled through constitutive]')
print(f'  => Play model needs direct C++ implementation (not energy)')


# =================================================================
# Plot
# =================================================================
print('\n=== Generating plots ===')

fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# Shape functions
ax = axes[0, 0]
r_plot = np.linspace(0, Bmax / 3, 500)
for k in range(0, min(Np, 20), 2):
    f_vals = [float(f_interps[k](r)) for r in r_plot]
    ax.plot(r_plot, f_vals, linewidth=0.8)
ax.set_xlabel('|p| (T)')
ax.set_ylabel('f(|p|) (A/m)')
ax.set_title('Shape functions g_k (B-input)')
ax.grid(True, alpha=0.3)

# Sinusoidal B-H loop
ax = axes[0, 1]
ax.plot(Hx[400:] / 1000, Bx_traj[400:], 'r-', linewidth=0.8)
ax.set_xlabel('H (kA/m)')
ax.set_ylabel('B (T)')
ax.set_title(f'Play model (Bm=1.05T)')
ax.grid(True, alpha=0.3)

# MATLAB comparison
ax = axes[0, 2]
if Hx_ref is not None and Hx_py is not None:
    ax.plot(Hx_ref / 1000, Bx_ref, 'k-', linewidth=1.0, label='MATLAB')
    ax.plot(Hx_py / 1000, Bx_ref, 'r--', linewidth=0.8, label='Python')
    ax.legend(fontsize=8)
ax.set_xlabel('H (kA/m)')
ax.set_ylabel('B (T)')
ax.set_title('Python vs MATLAB reference')
ax.grid(True, alpha=0.3)

# Multiple amplitudes
ax = axes[1, 0]
colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(amplitudes)))
for i, Bm_loop in enumerate(amplitudes):
    Bx_lp, Hx_lp, _ = loop_data[Bm_loop]
    ax.plot(Hx_lp / 1000, Bx_lp, color=colors[i],
            linewidth=0.8, label=f'{Bm_loop:.1f}T')
ax.set_xlabel('H (kA/m)')
ax.set_ylabel('B (T)')
ax.set_title('Multiple amplitude loops')
ax.legend(fontsize=7, ncol=2)
ax.grid(True, alpha=0.3)

# Descending branch
ax = axes[1, 1]
ax.plot(H_data / 1000, B_data, 'ko', markersize=3, label='Data')
ax.plot(Hx_desc / 1000, Bx_desc, 'r-', linewidth=1.0, label='Play')
ax.set_xlabel('H (kA/m)')
ax.set_ylabel('B (T)')
ax.set_title('Descending branch')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Error vs MATLAB
ax = axes[1, 2]
if Hx_ref is not None:
    ax.plot(ref['t'].flatten(), H_err, 'b-', linewidth=0.5)
    ax.set_xlabel('t (s)')
    ax.set_ylabel('|H_py - H_ref| (A/m)')
    ax.set_title(f'Error vs MATLAB (max={max_rel_err:.4f}%)')
ax.grid(True, alpha=0.3)

plt.tight_layout()
out_png = os.path.join(os.path.dirname(__file__), 'verify_play_to_energy.png')
plt.savefig(out_png, dpi=150)
print(f'  Saved: {out_png}')

# =================================================================
# Summary
# =================================================================
print(f'\n{"="*60}')
print(f'Results: {n_test} tests, '
      f'{"ALL PASS" if all_pass else "SOME FAILED"}')
print(f'{"="*60}')
