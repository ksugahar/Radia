"""Verify feasibility of coenergy model for B-input play model.

Mathematical analysis:
- Coenergy model requires f_k >= 0 and monotonically increasing (convex V_k)
- Play model: H = sum f_k(|p_k|) * sign(p_k)
- On descending branch at B=0: all p_k > 0 for k >= 1
- Therefore H(B=0) = sum f_k(eta_k) >= 0 with non-negative f_k
- But physically H(B=0) < 0 (past remanence on descending branch)
- CONTRADICTION: coenergy model cannot reproduce coercivity!

This script verifies the analysis empirically and tests NNLS identification.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import nnls

sys.path.insert(0, os.path.dirname(__file__))
from bqm_model import (
    compute_shape_functions, play_hysteron, load_binput_data,
    run_play_model,
)

DATA_DIR = (r'W:\999_菅原賢悟\19_磁気ヒステリシス'
            r'\2024_IGTE_共同研究\2024_03_08_H-input_B-input')
MAT_FILE = os.path.join(DATA_DIR, 'Potter_Schmulian', 'B_input.mat')

# =====================================================================
# 1. Load data and compute JMAG shape functions (reference play model)
# =====================================================================
print('=== Loading B-input data ===')
data = load_binput_data(MAT_FILE)
f_interps = compute_shape_functions(data['HdataLUT'], data['Bplay'])
eta = data['eta']
Np = len(eta)
Bmax = data['Bmax']
print(f'  Np={Np} operators, Bmax={Bmax:.3f} T, dB={data["dB"]:.4f} T')


# =====================================================================
# 2. Fundamental impossibility proof at B=0
# =====================================================================
print('\n=== Fundamental analysis at B=0 on descending branch ===')

# Initialize at saturation
px0 = np.full(Np, Bmax)
py0 = np.zeros(Np)

# Descend to B=0
n_desc = 200
B_desc = np.linspace(Bmax, 0.0, n_desc)
for i in range(n_desc):
    hx, hy = play_hysteron(f_interps, B_desc[i], 0.0, eta, px0, py0)

# At B=0, check all p_k values
print(f'  H(B=0) = {hx:.2f} A/m (should be negative)')
print(f'  Number of p_k > 0: {np.sum(px0 > 0)}/{Np}')
print(f'  Number of p_k = 0: {np.sum(px0 == 0)}/{Np}')
print(f'  Number of p_k < 0: {np.sum(px0 < 0)}/{Np}')
print(f'  p_k values: min={px0.min():.4f}, max={px0.max():.4f}')

# With non-negative f_k, what would H be?
# H_nonneg = sum |f_k(|p_k|)| * sign(p_k) where f_k is forced non-negative
f_original = np.array([float(f_interps[k](abs(px0[k]))) for k in range(Np)])
print(f'\n  Original f_k values at p_k(B=0):')
print(f'    f_0(0) = 0 (p_0=0, does not contribute)')
n_pos = np.sum(f_original[1:] > 0)
n_neg = np.sum(f_original[1:] < 0)
print(f'    f_k > 0: {n_pos}/{Np-1} operators')
print(f'    f_k < 0: {n_neg}/{Np-1} operators')
print(f'    sum f_k * sign(p_k) = {sum(f_original[k] * np.sign(px0[k]) for k in range(Np) if abs(px0[k]) > 1e-10):.2f} A/m')

# If all f_k were non-negative:
f_abs = np.abs(f_original)
H_nonneg = sum(f_abs[k] * np.sign(px0[k]) for k in range(Np) if abs(px0[k]) > 1e-10)
print(f'    sum |f_k| * sign(p_k) = {H_nonneg:.2f} A/m (all positive!)')
print(f'\n  CONCLUSION: Non-negative f_k gives H(B=0) > 0,')
print(f'  but physics requires H(B=0) < 0 (coercive region)')
print(f'  => Coenergy model IMPOSSIBLE for B-input play')


# =====================================================================
# 3. NNLS attempt (empirical confirmation)
# =====================================================================
print('\n=== NNLS identification attempt ===')

# Use fewer operators for better conditioning
K_co = 20  # number of coenergy operators
eta_co = np.linspace(0, Bmax * 0.95, K_co)
M_grid = 15  # grid points per operator
r_grid = np.linspace(0, Bmax, M_grid)

# Build system: for each descending branch data point,
# compute p_k and assemble the equation H = sum f_k(|p_k|) * sign(p_k)
BH_curves = data['BH']
BMax_arr = data['BMax']

rows = []  # system matrix rows
rhs = []   # right-hand side (H values)

for n_branch in range(len(BH_curves)):
    H_br, B_br = BH_curves[n_branch]
    BMax_n = BMax_arr[-(n_branch + 1)]  # BMax ascending order

    # Initialize play operators at saturation
    p_k = np.full(K_co, BMax_n)

    for j in range(len(B_br)):
        B_j = B_br[j]

        # Update play operators
        for k in range(K_co):
            if eta_co[k] <= 0:
                p_k[k] = B_j
            else:
                diff = abs(B_j - p_k[k])
                if diff <= eta_co[k]:
                    pass  # stuck
                else:
                    p_k[k] = B_j + eta_co[k] * (p_k[k] - B_j) / diff

        # Build row: H_j = sum_k sign(p_k) * f_k(|p_k|)
        # f_k is piecewise-linear on r_grid, so f_k(r) is linear in grid values
        row = np.zeros(K_co * M_grid)
        for k in range(K_co):
            r_k = abs(p_k[k])
            s_k = np.sign(p_k[k]) if abs(p_k[k]) > 1e-15 else 0.0
            if s_k == 0:
                continue

            # Find grid interval
            if r_k <= r_grid[0]:
                row[k * M_grid] += s_k  # f_k(0)
            elif r_k >= r_grid[-1]:
                row[k * M_grid + M_grid - 1] += s_k  # f_k(r_max)
            else:
                idx = np.searchsorted(r_grid, r_k) - 1
                idx = max(0, min(idx, M_grid - 2))
                t = (r_k - r_grid[idx]) / (r_grid[idx + 1] - r_grid[idx])
                row[k * M_grid + idx] += s_k * (1 - t)
                row[k * M_grid + idx + 1] += s_k * t

        rows.append(row)
        rhs.append(H_br[j])

A = np.array(rows)
b = np.array(rhs)
print(f'  System size: {A.shape[0]} equations, {A.shape[1]} unknowns')

# Solve NNLS
x_nnls, residual_nnls = nnls(A, b)
H_pred = A @ x_nnls
nnls_err = np.sqrt(np.mean((H_pred - b)**2))
H_range = b.max() - b.min()
nnls_rel = nnls_err / H_range * 100

print(f'  NNLS residual: {residual_nnls:.2f}')
print(f'  RMSE: {nnls_err:.2f} A/m ({nnls_rel:.1f}% of H range)')

# Check H(B=0) prediction
# Find data point closest to B=0 on largest descending branch
H_largest, B_largest = BH_curves[0]
idx_zero = np.argmin(np.abs(B_largest))
print(f'  H_data at B~0: {H_largest[idx_zero]:.2f} A/m (B={B_largest[idx_zero]:.4f})')

# Get the prediction
if idx_zero < len(H_pred):
    # Find the right row index for this data point
    offset = 0
    for n_branch in range(len(BH_curves)):
        H_br, B_br = BH_curves[n_branch]
        idx_z = np.argmin(np.abs(B_br))
        if n_branch == 0:
            print(f'  H_nnls at B~0 (branch 0): {H_pred[offset + idx_z]:.2f} A/m')
            break
        offset += len(B_br)


# =====================================================================
# 4. Compare with unconstrained least squares
# =====================================================================
print('\n=== Unconstrained least squares (allows negative f_k) ===')
x_lstsq, res_lstsq, _, _ = np.linalg.lstsq(A, b, rcond=None)
H_lstsq = A @ x_lstsq
lstsq_err = np.sqrt(np.mean((H_lstsq - b)**2))
lstsq_rel = lstsq_err / H_range * 100
print(f'  RMSE: {lstsq_err:.2f} A/m ({lstsq_rel:.1f}% of H range)')
n_neg_lstsq = np.sum(x_lstsq < -1e-10)
print(f'  Negative coefficients: {n_neg_lstsq}/{len(x_lstsq)}')


# =====================================================================
# 5. Test the NNLS model on sinusoidal trajectory
# =====================================================================
print('\n=== NNLS model sinusoidal test ===')

# Reconstruct f_k from NNLS solution
from scipy.interpolate import interp1d
f_co = []
for k in range(K_co):
    f_vals = x_nnls[k * M_grid:(k + 1) * M_grid]
    f_co.append(interp1d(r_grid, f_vals, kind='linear',
                         bounds_error=False,
                         fill_value=(f_vals[0], f_vals[-1])))

# Run sinusoidal test with coenergy model
Bm_test = 1.5
t_test = np.linspace(0, 6 * np.pi, 1000)
Bx_test = Bm_test * np.sin(2 * t_test)

# Run original play model
Hx_play, _ = run_play_model(f_interps, eta, Bx_test, np.zeros_like(t_test))

# Run coenergy model
p_k_co = np.zeros(K_co)
Hx_co = np.zeros(len(t_test))
for nt in range(len(t_test)):
    Bx = Bx_test[nt]
    for k in range(K_co):
        if eta_co[k] <= 0:
            p_k_co[k] = Bx
        else:
            diff = abs(Bx - p_k_co[k])
            if diff > eta_co[k]:
                p_k_co[k] = Bx + eta_co[k] * (p_k_co[k] - Bx) / diff
    hx = 0.0
    for k in range(K_co):
        if abs(p_k_co[k]) < 1e-15:
            continue
        r_k = abs(p_k_co[k])
        s_k = np.sign(p_k_co[k])
        hx += float(f_co[k](r_k)) * s_k
    Hx_co[nt] = hx

# Compare
H_err_traj = np.abs(Hx_co[400:] - Hx_play[400:])
H_range_traj = Hx_play[400:].max() - Hx_play[400:].min()
traj_err = H_err_traj.max() / H_range_traj * 100
print(f'  Play model H range: {H_range_traj:.1f} A/m')
print(f'  NNLS model max error: {H_err_traj.max():.1f} A/m ({traj_err:.1f}%)')

# Check if hysteresis is opening in correct direction
# On descending branch, H should be below ascending branch at same B
idx_desc = slice(300, 550)  # descending part
idx_asc = slice(550, 800)   # ascending part
print(f'  Play: H_mean desc={np.mean(Hx_play[idx_desc]):.1f}, '
      f'asc={np.mean(Hx_play[idx_asc]):.1f} A/m')
print(f'  NNLS: H_mean desc={np.mean(Hx_co[idx_desc]):.1f}, '
      f'asc={np.mean(Hx_co[idx_asc]):.1f} A/m')


# =====================================================================
# 6. Plot comparison
# =====================================================================
fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# Shape function signs (original JMAG)
ax = axes[0, 0]
r_plot = np.linspace(0.01, Bmax * 0.9, 200)
for k in range(0, min(Np, 20), 2):
    f_vals = [float(f_interps[k](r)) for r in r_plot]
    ax.plot(r_plot, f_vals, linewidth=0.8, label=f'k={k}' if k < 6 else None)
ax.axhline(y=0, color='k', linewidth=0.5, linestyle='--')
ax.set_xlabel('|p| (T)')
ax.set_ylabel('f_k(|p|) (A/m)')
ax.set_title('JMAG shape functions (signed)')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)

# NNLS shape functions
ax = axes[0, 1]
for k in range(0, K_co, max(1, K_co // 10)):
    f_vals = x_nnls[k * M_grid:(k + 1) * M_grid]
    ax.plot(r_grid, f_vals, linewidth=0.8, label=f'k={k}')
ax.set_xlabel('|p| (T)')
ax.set_ylabel('f_k(|p|) (A/m)')
ax.set_title(f'NNLS shape functions (K={K_co}, non-negative)')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)

# Descending branch fit
ax = axes[0, 2]
H_br0, B_br0 = BH_curves[0]
ax.plot(H_br0 / 1000, B_br0, 'k-', linewidth=1.5, label='Data')
# Reconstruct prediction for branch 0
n_pts_0 = len(B_br0)
ax.plot(H_pred[:n_pts_0] / 1000, B_br0, 'r--', linewidth=1.0, label='NNLS')
ax.axhline(y=0, color='gray', linewidth=0.5, linestyle=':')
ax.axvline(x=0, color='gray', linewidth=0.5, linestyle=':')
ax.set_xlabel('H (kA/m)')
ax.set_ylabel('B (T)')
ax.set_title(f'Descending branch (NNLS err={nnls_rel:.1f}%)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Sinusoidal B-H: play model
ax = axes[1, 0]
ax.plot(Hx_play[400:] / 1000, Bx_test[400:], 'k-', linewidth=0.8)
ax.set_xlabel('H (kA/m)')
ax.set_ylabel('B (T)')
ax.set_title(f'Play model (Bm={Bm_test}T)')
ax.grid(True, alpha=0.3)

# Sinusoidal B-H: NNLS coenergy
ax = axes[1, 1]
ax.plot(Hx_co[400:] / 1000, Bx_test[400:], 'r-', linewidth=0.8)
ax.set_xlabel('H (kA/m)')
ax.set_ylabel('B (T)')
ax.set_title(f'NNLS coenergy model (err={traj_err:.1f}%)')
ax.grid(True, alpha=0.3)

# Overlay comparison
ax = axes[1, 2]
ax.plot(Hx_play[400:] / 1000, Bx_test[400:], 'k-', linewidth=1.0, label='Play')
ax.plot(Hx_co[400:] / 1000, Bx_test[400:], 'r--', linewidth=0.8, label='NNLS')
ax.set_xlabel('H (kA/m)')
ax.set_ylabel('B (T)')
ax.set_title('Play vs NNLS coenergy')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.suptitle('Coenergy Model Feasibility Analysis\n'
             '(Non-negative f_k cannot produce H < 0 when all p_k > 0)',
             fontsize=12, fontweight='bold')
plt.tight_layout()
out_png = os.path.join(os.path.dirname(__file__), 'verify_coenergy_feasibility.png')
plt.savefig(out_png, dpi=150)
print(f'\nSaved: {out_png}')

# =====================================================================
# Summary
# =====================================================================
print(f'\n{"="*70}')
print('CONCLUSION:')
print('  The coenergy model (convex V_k, non-negative f_k) CANNOT reproduce')
print('  the B-input play model behavior.')
print()
print('  Proof: On the descending branch at B=0, all play operators p_k > 0')
print('  (for k >= 1). With non-negative f_k, H(B=0) = sum f_k(eta_k) >= 0.')
print('  But physics requires H(B=0) < 0 (past remanence point).')
print()
print('  The JMAG decomposition overcomes this by allowing f_k < 0 for k > 0,')
print('  but negative f_k means non-convex V_k, which breaks the energy')
print('  minimization formulation (no unique minimum).')
print()
print('  RECOMMENDATION: Implement B-input play model directly in C++ (Option A)')
print('  - Forward (B -> H): direct play evaluation, O(K) per step')
print('  - Inverse (H -> B): Newton iteration on forward model')
print(f'{"="*70}')
